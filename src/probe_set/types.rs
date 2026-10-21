use anyhow::{ensure, Context, Result};
use clap::Parser;
use serde::Serialize;
use std::fs;
use std::path::{Path, PathBuf};

/// Nominal length for each synthetic per-gene contig — see
/// `feature_capture.rs`'s identical constant for why the exact value
/// doesn't matter (downstream tools never read real coordinates off this
/// BAM, only tags).
pub const SYNTHETIC_CONTIG_LEN: i64 = 1_000_000;
/// Center of the probe-sequence window, in bp from the start of the emitted
/// R2 read.
///
/// Was empirically 69 (the Flex design notes Phase 6) when `extract()` bounded
/// R2's far end with the generic Illumina adapter search, which mostly
/// failed to find anything (13.5% hit rate on real data) and fell back to
/// the raw read's own end -- an arbitrary, per-read-variable amount of
/// trailing junk past the real probe insert, hence the wide empirical spread
/// and the offset landing well past 0.
///
/// Now 0: `extract()`'s `ExtractCfg::flex_r2_boundary` (2026-09-17) bounds R2
/// at the real structural landmark -- the constant sequence immediately
/// following `Ligated_Probe_insert` -- so the probe insert is the first
/// thing in the emitted R2, always, for any read where `extract()` located
/// that boundary at all. `--window` still absorbs the residual few-bp jitter
/// from `extract()`'s own anchor-match tolerance, just centered at 0 now
/// instead of a wide, mostly-arbitrary spread.
pub const DEFAULT_PROBE_SEQ_OFFSET: usize = 0;

/// A probe-set reference row.
#[derive(Debug, Clone)]
pub struct Probe {
    pub gene_id: String,
    pub gene_name: String,
    pub probe_id: String,
    pub region: String,
    pub seq: Vec<u8>,
}

#[derive(Debug)]
pub struct ProbeSetReference {
    pub probes: Vec<Probe>,
}

impl ProbeSetReference {
    /// Parse a real cellranger probe-set CSV. Leading `#key=value` metadata
    /// lines (panel name, reference genome, ...) are skipped verbatim.
    pub fn load(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("reading probe set {}", path.display()))?;
        let mut lines = text.lines().filter(|l| !l.starts_with('#'));
        let header = lines.next().context("empty probe set (no header)")?;
        let cols: Vec<&str> = header.split(',').collect();
        let col = |name: &str| {
            cols.iter()
                .position(|c| *c == name)
                .with_context(|| format!("probe set missing '{name}' column"))
        };
        let (gene_id_i, seq_i, probe_id_i, included_i, region_i, gene_name_i) = (
            col("gene_id")?,
            col("probe_seq")?,
            col("probe_id")?,
            col("included")?,
            col("region")?,
            col("gene_name")?,
        );

        let mut probes = Vec::new();
        let mut seq_len = None;
        for (lineno, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split(',').collect();
            let max = [
                gene_id_i,
                seq_i,
                probe_id_i,
                included_i,
                region_i,
                gene_name_i,
            ]
            .into_iter()
            .max()
            .unwrap();
            ensure!(
                f.len() > max,
                "{}:{}: row has too few columns: '{line}'",
                path.display(),
                lineno + 2
            );
            if f[included_i] != "TRUE" {
                continue;
            }
            let seq = f[seq_i].as_bytes().to_vec();
            match seq_len {
                None => seq_len = Some(seq.len()),
                Some(len) => ensure!(
                    seq.len() == len,
                    "{}:{}: probe_seq length {} != other probes' length {len} -- this parser assumes one fixed length",
                    path.display(),
                    lineno + 2,
                    seq.len()
                ),
            }
            probes.push(Probe {
                gene_id: f[gene_id_i].to_string(),
                gene_name: f[gene_name_i].to_string(),
                probe_id: f[probe_id_i].to_string(),
                region: f[region_i].to_string(),
                seq,
            });
        }
        ensure!(
            !probes.is_empty(),
            "probe set {} has no included=TRUE rows",
            path.display()
        );
        Ok(Self { probes })
    }
}

/// Assign Flex GEX reads to genes by probe-set lookup, not genome alignment.
#[derive(Debug, Parser)]
#[command(
    name = "morselle-probeset-align",
    version,
    about = "Assign 10x Flex GEX reads to genes via the real probe-set reference (not genome alignment)"
)]
pub struct Cli {
    /// `--minimap2`-preset R2 FASTQ from `extract` (carries `CR/CB/UR/UB`,
    /// and `PR/PB` when a Flex probe barcode was also located).
    #[arg(value_name = "R2_FASTQ")]
    pub r2_fastq: PathBuf,

    /// Real cellranger probe-set CSV (`gene_id,probe_seq,probe_id,included,region,gene_name`).
    #[arg(long)]
    pub probe_set: PathBuf,

    /// Output BAM: mapped (synthetic per-gene contig, `GX`/`GN` tags) on a
    /// match, unmapped (correctly skipped by `molecules`/`matrix`) on none.
    #[arg(short, long, value_name = "BAM")]
    pub output: PathBuf,

    /// Canonical `id\tname\tGene Expression` manifest for `matrix --features`.
    /// Created if absent; verified (not overwritten) if it already exists
    /// and doesn't match -- mirrors `annotate`'s `--genes-mtx`.
    #[arg(long)]
    pub genes_out: Option<PathBuf>,

    /// Optional metrics JSON (hit rate).
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    /// Max mismatches for a confident probe-sequence match. Kept small by
    /// default: `HammingIndex` precomputes every ≤max_mm neighbor of every
    /// probe, and neighbor count grows combinatorially (`C(50,k)*3^k`) --
    /// 1 is ~8M entries (fast, a few hundred MB); 2 would be ~590M (tens of
    /// GB), impractical to build.
    #[arg(long, default_value_t = 1)]
    pub max_mm: u32,

    /// Max bp each side of [`DEFAULT_PROBE_SEQ_OFFSET`] searched, to absorb
    /// per-read jitter in the measured dT-run boundary (same mechanism as
    /// `--probe-barcode-window`/`--anchor-window` elsewhere in this plan).
    #[arg(long, default_value_t = 15)]
    pub window: isize,

    /// Worker threads for parallel probe matching and BAM record construction (0 = auto/num_cpus).
    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,
}

#[derive(Debug, Serialize)]
pub struct Metrics {
    pub reads: u64,
    pub matched: u64,
    pub hit_rate_pct: f64,
    pub distinct_genes_hit: usize,
}
