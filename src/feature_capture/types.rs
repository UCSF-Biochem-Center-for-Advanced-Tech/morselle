use anyhow::{ensure, Context, Result};
use clap::Parser;
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};

/// Nominal length for each synthetic per-feature contig. Coordinates on
/// this BAM are never real genomic positions — downstream (`molecules`,
/// `matrix`) only read tags, never position/CIGAR — so this just needs to
/// comfortably exceed any real R2 read length.
pub const SYNTHETIC_CONTIG_LEN: i64 = 1_000_000;

/// One row of a `feature_reference.csv`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub id: String,
    pub name: String,
    pub read: String,
    pub feature_type: String,
    pub(crate) anchor: Vec<u8>,
    pub(crate) barcode: Vec<u8>,
}

#[derive(Debug)]
pub struct FeatureReference {
    pub features: Vec<Feature>,
}

impl FeatureReference {
    /// Parse the real cellranger `feature_reference.csv` schema:
    /// `id,name,read,pattern,sequence,feature_type`. `pattern` must be
    /// `<anchor>(BC)` — confirmed layout against real CAT06 data. A pattern
    /// starting with `(BC)` (barcode before anchor) is a real, if rarer,
    /// 10x-supported layout this parser doesn't handle; flagged explicitly
    /// rather than silently mis-parsed as garbage.
    pub fn load(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("reading feature reference {}", path.display()))?;
        let mut lines = text.lines();
        let header = lines.next().context("empty feature reference")?;
        let cols: Vec<&str> = header.split(',').collect();
        let col = |name: &str| {
            cols.iter()
                .position(|c| *c == name)
                .with_context(|| format!("feature reference missing '{name}' column"))
        };
        let (id_i, name_i, read_i, pattern_i, seq_i, type_i) = (
            col("id")?,
            col("name")?,
            col("read")?,
            col("pattern")?,
            col("sequence")?,
            col("feature_type")?,
        );

        let mut features = Vec::new();
        for (lineno, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split(',').collect();
            let max = [id_i, name_i, read_i, pattern_i, seq_i, type_i]
                .into_iter()
                .max()
                .unwrap();
            ensure!(
                f.len() > max,
                "{}:{}: row has too few columns: '{line}'",
                path.display(),
                lineno + 2
            );
            let pattern = f[pattern_i].trim();
            ensure!(
                !pattern.starts_with("(BC)"),
                "{}:{}: pattern '{pattern}' is (BC)<anchor> (barcode-then-anchor) -- \
                 not supported by this parser, which only handles <anchor>(BC)",
                path.display(),
                lineno + 2
            );
            let anchor = pattern
                .strip_suffix("(BC)")
                .with_context(|| {
                    format!(
                        "{}:{}: unsupported pattern (expected <anchor>(BC)): '{pattern}'",
                        path.display(),
                        lineno + 2
                    )
                })?
                .as_bytes()
                .to_vec();
            ensure!(
                !anchor.is_empty(),
                "{}:{}: empty anchor in pattern '{pattern}'",
                path.display(),
                lineno + 2
            );
            features.push(Feature {
                id: f[id_i].trim().to_string(),
                name: f[name_i].trim().to_string(),
                read: f[read_i].trim().to_string(),
                feature_type: f[type_i].trim().to_string(),
                anchor,
                barcode: f[seq_i].trim().as_bytes().to_vec(),
            });
        }
        ensure!(
            !features.is_empty(),
            "feature reference {} has no rows",
            path.display()
        );
        Ok(Self { features })
    }
}

/// Assign reads to features by anchor+capture-barcode matching.
#[derive(Debug, Parser)]
#[command(
    name = "morselle-guides",
    version,
    about = "Assign 10x Feature Barcoding reads (CRISPR guides, antibody/CMO tags, ...) to features"
)]
pub struct Cli {
    /// `--minimap2`-preset R2 FASTQ from `extract` (carries `CR/CB/UR/UB`,
    /// and `PR/PB` when the run also located a Flex probe barcode -- `PB` is
    /// optional, not required, since plain non-Flex GEX + Feature Barcoding
    /// carries no probe barcode at all).
    #[arg(value_name = "R2_FASTQ")]
    pub r2_fastq: PathBuf,

    /// `feature_reference.csv`: `id,name,read,pattern,sequence,feature_type`.
    #[arg(long)]
    pub feature_reference: PathBuf,

    /// Output BAM: mapped (synthetic per-feature contig, `GX`/`GN` tags) on
    /// a match, unmapped (correctly skipped by `molecules`/`matrix`) on none.
    #[arg(short, long, value_name = "BAM")]
    pub output: PathBuf,

    /// Canonical `id\tname\tfeature_type` manifest for `matrix --features`.
    /// Created if absent; verified (not overwritten) if it already exists
    /// and doesn't match -- mirrors `annotate`'s `--genes-mtx`.
    #[arg(long)]
    pub features_out: Option<PathBuf>,

    /// Optional metrics JSON (hit rate, per-feature call counts).
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    /// Max mismatches for a confident anchor match.
    #[arg(long, default_value_t = 2)]
    pub anchor_max_mm: u32,

    /// Max mismatches for a confident capture-barcode match.
    #[arg(long, default_value_t = 1)]
    pub capture_max_mm: u32,

    /// Max bp from the read start searched for the anchor. Default is wider
    /// than the plan's original placeholder (5): the Flex design notes'
    /// Phase 0 empirically confirmed the real CAT06 CRISPR anchor sits at
    /// emitted-R2 offset ~10 with some jitter (46,938 vs. 2,301 hits at the
    /// next-best offset) -- this default covers that with real margin for
    /// other constructs, rather than re-guessing a narrower one.
    #[arg(long, default_value_t = 30)]
    pub anchor_window: usize,
}

#[derive(Debug, Serialize)]
pub struct FeatureCount {
    pub id: String,
    pub calls: u64,
}

#[derive(Debug, Serialize)]
pub struct Metrics {
    pub reads: u64,
    pub matched: u64,
    pub hit_rate_pct: f64,
    pub per_feature: Vec<FeatureCount>,
}
