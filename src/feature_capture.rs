//! 10x Feature Barcoding read -> feature assignment (FLEX_SUPPORT_PLAN.md
//! Phase 4, generalized per the plan's own amendment).
//!
//! Not CRISPR- or Flex-specific by design: 10x's `<anchor>(BC)` read
//! structure + `feature_reference.csv` schema (`id,name,read,pattern,
//! sequence,feature_type`) is the same mechanism behind CRISPR Guide
//! Capture, Antibody Capture (CITE-seq), and Cell Multiplexing (CMO
//! hashing), and runs on ordinary (non-Flex) 3'/5' GEX libraries with no
//! probe barcode at all just as often as on Flex. `FeatureMatcher` only
//! does anchor+capture-barcode matching; it carries `feature_type` through
//! as data rather than branching on it — any per-`feature_type` calling
//! logic (e.g. Phase 5's CRISPR-specific GMM protospacer calling) is a
//! separate, later stage.

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write as _};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{ensure, Context, Result};
use clap::Parser;
use needletail::parse_fastx_file;
use rust_htslib::bam::{
    self, header::HeaderRecord, record::Aux, record::Cigar, record::CigarString,
};
use serde::Serialize;
use tracing::info;

use crate::hamming_index::HammingIndex;

/// Nominal length for each synthetic per-feature contig. Coordinates on
/// this BAM are never real genomic positions — downstream (`molecules`,
/// `matrix`) only read tags, never position/CIGAR — so this just needs to
/// comfortably exceed any real R2 read length.
const SYNTHETIC_CONTIG_LEN: i64 = 1_000_000;

/// One row of a `feature_reference.csv`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Feature {
    pub id: String,
    pub name: String,
    pub read: String,
    pub feature_type: String,
    anchor: Vec<u8>,
    barcode: Vec<u8>,
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

/// All features sharing one exact anchor sequence + capture-barcode length
/// (the unit a single `HammingIndex` can serve, since Hamming needs a fixed
/// key length).
struct PatternGroup {
    capture_len: usize,
    capture_index: HammingIndex<usize>,
}

/// Matches a read against a [`FeatureReference`]: locate the anchor (a
/// windowed, Hamming-tolerant forward scan from the read start, mirroring
/// `demux::assign_windowed`/`probe_barcode::locate`), then the following
/// capture barcode against the anchor's matching feature group.
///
/// Handles a reference with more than one distinct anchor sequence or
/// capture-barcode length (real for a panel mixing feature types, e.g.
/// CRISPR Guide Capture + Antibody Capture in one experiment) by grouping
/// features into anchor-length buckets, each with its own windowed search.
pub struct FeatureMatcher {
    anchor_by_len: Vec<(usize, HammingIndex<usize>)>,
    groups: Vec<PatternGroup>,
    anchor_window: usize,
}

impl FeatureMatcher {
    pub fn build(
        features: &[Feature],
        anchor_max_mm: u32,
        capture_max_mm: u32,
        anchor_window: usize,
    ) -> Self {
        let mut group_of: HashMap<(&[u8], usize), usize> = HashMap::new();
        let mut groups: Vec<PatternGroup> = Vec::new();
        let mut anchor_by_len: HashMap<usize, HammingIndex<usize>> = HashMap::new();

        for (fi, f) in features.iter().enumerate() {
            let key = (f.anchor.as_slice(), f.barcode.len());
            let group_id = *group_of.entry(key).or_insert_with(|| {
                groups.push(PatternGroup {
                    capture_len: f.barcode.len(),
                    capture_index: HammingIndex::new(),
                });
                let id = groups.len() - 1;
                anchor_by_len.entry(f.anchor.len()).or_default().insert(
                    &f.anchor,
                    id,
                    anchor_max_mm,
                );
                id
            });
            groups[group_id]
                .capture_index
                .insert(&f.barcode, fi, capture_max_mm);
        }

        Self {
            anchor_by_len: anchor_by_len.into_iter().collect(),
            groups,
            anchor_window,
        }
    }

    /// Locate the anchor within the first `anchor_window` bp of `read`
    /// (leftmost, smallest anchor length first when several are tried at
    /// the same position — real references have one anchor length in
    /// practice, so this rarely matters), then match the following capture
    /// barcode. Returns the matched feature's index into the original
    /// `features` slice.
    pub fn locate(&self, read: &[u8]) -> Option<usize> {
        for (anchor_len, index) in &self.anchor_by_len {
            let hay_end = read.len().min(self.anchor_window + anchor_len);
            if hay_end < *anchor_len {
                continue;
            }
            for offset in 0..=(hay_end - anchor_len) {
                let Some(group_id) = index.get(&read[offset..offset + anchor_len]) else {
                    continue;
                };
                let group = &self.groups[group_id];
                let bc_start = offset + anchor_len;
                let bc_end = bc_start + group.capture_len;
                if bc_end > read.len() {
                    continue;
                }
                if let Some(feature_idx) = group.capture_index.get(&read[bc_start..bc_end]) {
                    return Some(feature_idx);
                }
            }
        }
        None
    }
}

/// Assign reads to features by anchor+capture-barcode matching (10x
/// Feature Barcoding: CRISPR Guide Capture, Antibody Capture, Cell
/// Multiplexing, ...). Named `guides` as its concrete CAT06 entry point;
/// the matching core (`FeatureMatcher`) is assay-agnostic.
#[derive(Debug, Parser)]
#[command(
    name = "sbxpress-guides",
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
    /// than the plan's original placeholder (5): FLEX_SUPPORT_PLAN.md's
    /// Phase 0 empirically confirmed the real CAT06 CRISPR anchor sits at
    /// emitted-R2 offset ~10 with some jitter (46,938 vs. 2,301 hits at the
    /// next-best offset) -- this default covers that with real margin for
    /// other constructs, rather than re-guessing a narrower one.
    #[arg(long, default_value_t = 30)]
    pub anchor_window: usize,
}

fn parse_tag<'a>(comment: &'a str, tag: &str) -> Option<&'a str> {
    let start = comment.find(tag)? + tag.len();
    let rest = &comment[start..];
    let end = rest.find('\t').unwrap_or(rest.len());
    Some(&rest[..end])
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

fn features_tsv_bytes(features: &[Feature]) -> Vec<u8> {
    let mut out = Vec::new();
    for f in features {
        out.extend_from_slice(f.id.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(f.name.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(f.feature_type.as_bytes());
        out.push(b'\n');
    }
    out
}

/// Create-or-verify `--features-out`, mirroring `annotate_bam.rs`'s
/// `establish_genes_mtx`.
fn establish_features_tsv(path: &Path, features: &[Feature]) -> Result<()> {
    let expected = features_tsv_bytes(features);
    if path.exists() {
        let observed = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        ensure!(
            observed == expected,
            "existing feature manifest {} does not match this feature reference; \
             use a matching CSV or choose a new path",
            path.display()
        );
        return Ok(());
    }
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, &expected).with_context(|| format!("writing {}", path.display()))
}

pub fn run(cli: &Cli) -> Result<()> {
    let started = Instant::now();
    let reference = FeatureReference::load(&cli.feature_reference)?;
    info!(features = reference.features.len(), path = %cli.feature_reference.display(), "loaded feature reference");
    let matcher = FeatureMatcher::build(
        &reference.features,
        cli.anchor_max_mm,
        cli.capture_max_mm,
        cli.anchor_window,
    );

    if let Some(path) = &cli.features_out {
        establish_features_tsv(path, &reference.features)?;
    }

    let mut header = bam::Header::new();
    for f in &reference.features {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", &f.id)
                .push_tag(b"LN", SYNTHETIC_CONTIG_LEN),
        );
    }
    if let Some(parent) = cli.output.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    let mut writer = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)
        .with_context(|| format!("creating {}", cli.output.display()))?;

    let mut reads = 0u64;
    let mut matched = 0u64;
    let mut per_feature = vec![0u64; reference.features.len()];

    // See `probe_set.rs`'s identical handling: a genuinely empty R2 input
    // is a valid "zero reads" case (e.g. an orchestrator's chunk that
    // happened to get no reads), not a fatal error.
    let mut fastq = match parse_fastx_file(&cli.r2_fastq) {
        Ok(reader) => Some(reader),
        Err(err) if err.kind == needletail::errors::ParseErrorKind::EmptyFile => {
            info!(path = %cli.r2_fastq.display(), "R2 input is empty; writing header-only output");
            None
        }
        Err(err) => {
            return Err(err).with_context(|| format!("opening {}", cli.r2_fastq.display()));
        }
    };
    while let Some(rec) = fastq.as_mut().and_then(|f| f.next()) {
        let rec = rec.context("reading R2 FASTQ record")?;
        reads += 1;
        let header_line = String::from_utf8_lossy(rec.id());
        let (qname, comment) = header_line.split_once('\t').unwrap_or((&header_line, ""));
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);

        let mut record = bam::Record::new();
        let hit = matcher.locate(&seq);
        if let Some(feature_idx) = hit {
            matched += 1;
            per_feature[feature_idx] += 1;
            let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
            record.set(qname.as_bytes(), Some(&cigar), &seq, qual);
            record.set_tid(feature_idx as i32);
            record.set_pos(0);
            record.set_mapq(60);
            record.unset_unmapped();
            let f = &reference.features[feature_idx];
            record.push_aux(b"GX", Aux::String(&f.id))?;
            record.push_aux(b"GN", Aux::String(&f.name))?;
        } else {
            record.set(qname.as_bytes(), None, &seq, qual);
        }
        for tag in ["CR", "CB", "UR", "UB", "PR", "PB"] {
            if let Some(value) = parse_tag(comment, &format!("{tag}:Z:")) {
                if !value.is_empty() {
                    record.push_aux(tag.as_bytes(), Aux::String(value))?;
                }
            }
        }
        writer.write(&record)?;
    }
    drop(writer);

    let metrics = Metrics {
        reads,
        matched,
        hit_rate_pct: if reads > 0 {
            100.0 * matched as f64 / reads as f64
        } else {
            0.0
        },
        per_feature: reference
            .features
            .iter()
            .zip(per_feature)
            .map(|(f, calls)| FeatureCount {
                id: f.id.clone(),
                calls,
            })
            .collect(),
    };
    if let Some(path) = &cli.metrics {
        let mut w = BufWriter::new(File::create(path)?);
        w.write_all(serde_json::to_string_pretty(&metrics)?.as_bytes())?;
    }

    info!(
        reads = metrics.reads,
        matched = metrics.matched,
        hit_rate_pct = format!("{:.1}", metrics.hit_rate_pct),
        elapsed = format!("{:.1}s", started.elapsed().as_secs_f64()),
        out = %cli.output.display(),
        "RUN SUMMARY"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bam::Read as _;

    fn write_csv(text: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        f.write_all(text.as_bytes()).unwrap();
        f
    }

    /// Real CAT06 shape: two CRISPR guides sharing one anchor.
    const REAL_SHAPE: &str = "\
id,name,read,pattern,sequence,feature_type
CD45,CD45,R2,TTCCAGCTTAGCTCTTAAAC(BC),CTGTTACCCTAAGAACAAAC,CRISPR Guide Capture
inactive,inactive,R2,TTTCTGGACTTCAGTGGACG(BC),GGTGCACTCTGGAGAAGTCC,CRISPR Guide Capture
";

    #[test]
    fn loads_real_shape() {
        let f = write_csv(REAL_SHAPE);
        let r = FeatureReference::load(f.path()).unwrap();
        assert_eq!(r.features.len(), 2);
        assert_eq!(r.features[0].id, "CD45");
        assert_eq!(r.features[0].anchor, b"TTCCAGCTTAGCTCTTAAAC");
        assert_eq!(r.features[0].barcode, b"CTGTTACCCTAAGAACAAAC");
        assert_eq!(r.features[0].feature_type, "CRISPR Guide Capture");
    }

    #[test]
    fn rejects_barcode_then_anchor_pattern_explicitly() {
        let f = write_csv(
            "id,name,read,pattern,sequence,feature_type\nA,A,R2,(BC)ACGTACGTACGTACGTACGT,ACGTACGTACGTACGTACGT,CRISPR Guide Capture\n",
        );
        let err = FeatureReference::load(f.path()).unwrap_err();
        assert!(err.to_string().contains("(BC)<anchor>"));
    }

    #[test]
    fn rejects_pattern_without_bc_suffix() {
        let f = write_csv(
            "id,name,read,pattern,sequence,feature_type\nA,A,R2,ACGTACGTACGTACGTACGT,ACGTACGTACGTACGTACGT,CRISPR Guide Capture\n",
        );
        assert!(FeatureReference::load(f.path()).is_err());
    }

    fn read_with_anchor_capture(prefix_len: usize, anchor: &[u8], capture: &[u8]) -> Vec<u8> {
        let mut seq = vec![b'N'; prefix_len];
        seq.extend_from_slice(anchor);
        seq.extend_from_slice(capture);
        seq.extend_from_slice(b"NNNNNNNNNN");
        seq
    }

    #[test]
    fn locates_exact_match_at_confirmed_offset() {
        let f = write_csv(REAL_SHAPE);
        let r = FeatureReference::load(f.path()).unwrap();
        let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
        let read = read_with_anchor_capture(10, b"TTCCAGCTTAGCTCTTAAAC", b"CTGTTACCCTAAGAACAAAC");
        assert_eq!(matcher.locate(&read), Some(0));
    }

    #[test]
    fn tolerates_mismatches_within_budget() {
        let f = write_csv(REAL_SHAPE);
        let r = FeatureReference::load(f.path()).unwrap();
        let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
        let mut anchor = b"TTCCAGCTTAGCTCTTAAAC".to_vec();
        anchor[3] = b'G'; // 1 mismatch, within anchor_max_mm=2
        let mut capture = b"CTGTTACCCTAAGAACAAAC".to_vec();
        capture[0] = b'A'; // 1 mismatch, within capture_max_mm=1
        let read = read_with_anchor_capture(10, &anchor, &capture);
        assert_eq!(matcher.locate(&read), Some(0));
    }

    #[test]
    fn no_match_on_unrelated_sequence() {
        let f = write_csv(REAL_SHAPE);
        let r = FeatureReference::load(f.path()).unwrap();
        let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
        let read = vec![b'A'; 60];
        assert_eq!(matcher.locate(&read), None);
    }

    #[test]
    fn distinguishes_two_guides_sharing_one_anchor() {
        let f = write_csv(REAL_SHAPE);
        let r = FeatureReference::load(f.path()).unwrap();
        let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
        let read = read_with_anchor_capture(10, b"TTTCTGGACTTCAGTGGACG", b"GGTGCACTCTGGAGAAGTCC");
        assert_eq!(matcher.locate(&read), Some(1)); // "inactive"
    }

    /// Proves the generalization: a mixed panel (a CRISPR-flavored feature
    /// and an "Antibody Capture"-flavored one with a *different* anchor)
    /// resolves correctly -- the matcher doesn't special-case feature_type
    /// or assume one universal anchor.
    #[test]
    fn handles_a_mixed_panel_with_different_anchors_per_feature_type() {
        let csv = "\
id,name,read,pattern,sequence,feature_type
CD45,CD45,R2,TTCCAGCTTAGCTCTTAAAC(BC),CTGTTACCCTAAGAACAAAC,CRISPR Guide Capture
CD3,CD3-APC,R2,AAAAAAAAAAAAAAAAAAAA(BC),TCGATCGATCGATCGATCGA,Antibody Capture
";
        let f = write_csv(csv);
        let r = FeatureReference::load(f.path()).unwrap();
        let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);

        let crispr_read =
            read_with_anchor_capture(10, b"TTCCAGCTTAGCTCTTAAAC", b"CTGTTACCCTAAGAACAAAC");
        assert_eq!(matcher.locate(&crispr_read), Some(0));

        let ab_read =
            read_with_anchor_capture(10, b"AAAAAAAAAAAAAAAAAAAA", b"TCGATCGATCGATCGATCGA");
        assert_eq!(matcher.locate(&ab_read), Some(1));
    }

    #[test]
    fn features_tsv_create_then_verify_round_trips() {
        let f = write_csv(REAL_SHAPE);
        let r = FeatureReference::load(f.path()).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("features.tsv");
        establish_features_tsv(&out, &r.features).unwrap();
        let text = fs::read_to_string(&out).unwrap();
        assert_eq!(
            text,
            "CD45\tCD45\tCRISPR Guide Capture\ninactive\tinactive\tCRISPR Guide Capture\n"
        );
        // idempotent re-run against the same reference
        establish_features_tsv(&out, &r.features).unwrap();
    }

    #[test]
    fn features_tsv_rejects_mismatched_existing_file() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("features.tsv");
        fs::write(&out, "stale\tSTALE\tCRISPR Guide Capture\n").unwrap();
        let f = write_csv(REAL_SHAPE);
        let r = FeatureReference::load(f.path()).unwrap();
        assert!(establish_features_tsv(&out, &r.features).is_err());
    }

    #[test]
    fn parse_tag_extracts_value_and_stops_at_tab() {
        let comment = "CR:Z:ACGT\tCB:Z:ACGT-1\tPB:Z:D-A01";
        assert_eq!(parse_tag(comment, "CB:Z:"), Some("ACGT-1"));
        assert_eq!(parse_tag(comment, "PB:Z:"), Some("D-A01"));
        assert_eq!(parse_tag(comment, "XX:Z:"), None);
    }

    #[test]
    fn run_end_to_end_writes_mapped_and_unmapped_records_with_carried_tags() {
        let feature_csv = write_csv(REAL_SHAPE);
        let dir = tempfile::tempdir().unwrap();
        let fastq_path = dir.path().join("r2.fastq");
        let matched_seq =
            read_with_anchor_capture(10, b"TTCCAGCTTAGCTCTTAAAC", b"CTGTTACCCTAAGAACAAAC");
        let unmatched_seq = vec![b'A'; 60];
        let mut fastq = File::create(&fastq_path).unwrap();
        writeln!(
            fastq,
            "@read1\tCR:Z:ACGTACGTACGTACGT\tCB:Z:ACGTACGTACGTACGT-1\tUR:Z:AAAACCCCGGGG\tUB:Z:AAAACCCCGGGG\tPB:Z:D-A01"
        ).unwrap();
        writeln!(fastq, "{}", String::from_utf8_lossy(&matched_seq)).unwrap();
        writeln!(fastq, "+").unwrap();
        writeln!(fastq, "{}", "I".repeat(matched_seq.len())).unwrap();
        writeln!(
            fastq,
            "@read2\tCR:Z:TTTTGGGGCCCC\tCB:Z:TTTTGGGGCCCC-1\tUR:Z:GGGGAAAACCCC\tUB:Z:GGGGAAAACCCC"
        )
        .unwrap();
        writeln!(fastq, "{}", String::from_utf8_lossy(&unmatched_seq)).unwrap();
        writeln!(fastq, "+").unwrap();
        writeln!(fastq, "{}", "I".repeat(unmatched_seq.len())).unwrap();
        drop(fastq);

        let out_bam = dir.path().join("guides.bam");
        let features_out = dir.path().join("features.tsv");
        let metrics_path = dir.path().join("metrics.json");
        let cli = Cli {
            r2_fastq: fastq_path,
            feature_reference: feature_csv.path().to_path_buf(),
            output: out_bam.clone(),
            features_out: Some(features_out.clone()),
            metrics: Some(metrics_path.clone()),
            anchor_max_mm: 2,
            capture_max_mm: 1,
            anchor_window: 30,
        };
        run(&cli).unwrap();

        assert!(features_out.exists());
        let metrics: serde_json::Value =
            serde_json::from_str(&fs::read_to_string(&metrics_path).unwrap()).unwrap();
        assert_eq!(metrics["reads"], 2);
        assert_eq!(metrics["matched"], 1);

        let mut reader = bam::Reader::from_path(&out_bam).unwrap();
        let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
        assert_eq!(records.len(), 2);

        let mapped = &records[0];
        assert!(!mapped.is_unmapped());
        assert_eq!(mapped.qname(), b"read1");
        assert_eq!(mapped.aux(b"GX").unwrap(), Aux::String("CD45"));
        assert_eq!(
            mapped.aux(b"CB").unwrap(),
            Aux::String("ACGTACGTACGTACGT-1")
        );
        assert_eq!(mapped.aux(b"PB").unwrap(), Aux::String("D-A01"));

        let unmapped = &records[1];
        assert!(unmapped.is_unmapped());
        assert_eq!(unmapped.qname(), b"read2");
        assert!(unmapped.aux(b"GX").is_err());
        assert_eq!(unmapped.aux(b"CB").unwrap(), Aux::String("TTTTGGGGCCCC-1"));
        assert!(unmapped.aux(b"PB").is_err()); // PB optional, absent here
    }
}
