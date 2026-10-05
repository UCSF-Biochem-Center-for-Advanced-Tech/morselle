//! 10x Flex probe-set alignment (FLEX_SUPPORT_PLAN.md Phase 6 prerequisite):
//! assigns a Flex GEX read to a gene by direct lookup of its probe sequence
//! against the real probe-set reference, **not** genome alignment.
//!
//! Flex doesn't sequence free cDNA: each read's payload is (approximately)
//! one of a fixed, closed set of ~54k known 50bp probe sequences (one of a
//! few probes per targeted gene — see `Chromium_Human_Transcriptome_Probe_Set*.csv`).
//! Confirmed empirically this session (see FLEX_SUPPORT_PLAN.md Phase 6):
//! full-genome splice-aware alignment of real CAT06 Flex reads found almost
//! nothing (59 of 1.65M records) — the payload isn't cDNA a genome mapper
//! can place. This module is the real substitute: a windowed, Hamming-tolerant
//! lookup of the probe sequence directly against the reference, mirroring
//! `probe_barcode.rs`'s own windowed-locate design (same shared
//! `HammingIndex` primitive), and emitting a BAM shaped exactly like
//! `feature_capture.rs`'s (synthetic per-gene contigs, `GX`/`GN` tags) so it
//! feeds `molecules`/`matrix`/`cells`/`probe-demux` completely unchanged.

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

/// Nominal length for each synthetic per-gene contig — see
/// `feature_capture.rs`'s identical constant for why the exact value
/// doesn't matter (downstream tools never read real coordinates off this
/// BAM, only tags).
const SYNTHETIC_CONTIG_LEN: i64 = 1_000_000;

/// Center of the probe-sequence window, in bp from the start of the emitted
/// R2 read.
///
/// Was empirically 69 (FLEX_SUPPORT_PLAN.md Phase 6) when `extract()` bounded
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

/// A probe-set reference row (real cellranger schema:
/// `gene_id,probe_seq,probe_id,included,region,gene_name`), `included=TRUE`
/// only -- `FALSE` rows are 10x's own curation of probes too
/// cross-reactive/ambiguous to count on, already excluded upstream of us.
#[derive(Debug, Clone)]
pub struct Probe {
    pub gene_id: String,
    pub gene_name: String,
    pub probe_id: String,
    pub region: String,
    seq: Vec<u8>,
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

/// All single-position deletions of `seq` (`seq.len() - 1` variants, one per
/// deleted position). Real CAT06 data (2026-09-18 diagnostic, real
/// edit-distance recovery check against 20K Hamming-miss reads): of misses
/// recoverable within edit distance <=2, deletions dominate 85.2% vs 9.2%
/// insertion-only -- a real base the read is missing relative to the probe,
/// not a substitution a mismatch budget could ever cover.
fn single_deletion_variants(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    (0..seq.len()).map(move |i| {
        let mut v = Vec::with_capacity(seq.len() - 1);
        v.extend_from_slice(&seq[..i]);
        v.extend_from_slice(&seq[i + 1..]);
        v
    })
}

/// All single-position insertions into `seq` (up to `(seq.len() + 1) * 3`
/// variants: every position, every non-redundant base). The rarer indel
/// shape in real data (9.2% of recoverable indel-shaped misses vs 85.2%
/// deletion) but real, and cheap enough to cover alongside deletion.
fn single_insertion_variants(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    (0..=seq.len()).flat_map(move |i| {
        let existing = if i < seq.len() { seq[i] } else { b'\0' };
        b"ACGT"
            .iter()
            .filter(move |&&b| b != existing)
            .map(move |&b| {
                let mut v = Vec::with_capacity(seq.len() + 1);
                v.extend_from_slice(&seq[..i]);
                v.push(b);
                v.extend_from_slice(&seq[i..]);
                v
            })
    })
}

/// Matches a read's probe sequence against the reference via a windowed,
/// Hamming-tolerant lookup (see [`DEFAULT_PROBE_SEQ_OFFSET`]) -- structurally
/// identical to `probe_barcode.rs::ProbeWhitelist::locate`, just with a
/// much larger reference and no revcomp step (probe_seq is confirmed to
/// appear forward in the emitted read, unlike the probe barcode).
///
/// Also tries single-indel-shifted windows (`del_index`/`ins_index`) at
/// each offset -- real SBX errors are dominated by single-base indels, not
/// substitutions, inside the probe span itself (not just jitter in where it
/// starts, which the offset window already absorbs): of real Hamming-miss
/// reads recoverable by edit distance <=2, 98.5% were genuinely
/// indel-shaped. Rather than a per-read full-database edit-distance scan
/// (intractable at billions-of-reads scale), this extends the same
/// precomputed-neighborhood trick [`HammingIndex`] already uses for
/// substitutions to single indels too -- still an O(1) hash lookup per
/// attempt, just at two more fixed window lengths (49bp/51bp alongside the
/// exact 50bp).
pub struct ProbeMatcher {
    index: HammingIndex<usize>,
    del_index: HammingIndex<usize>,
    ins_index: HammingIndex<usize>,
    probe_len: usize,
    window: isize,
}

impl ProbeMatcher {
    pub fn build(probes: &[Probe], max_mm: u32, window: isize) -> Self {
        let probe_len = probes.first().map_or(0, |p| p.seq.len());
        let mut index = HammingIndex::new();
        let mut del_index = HammingIndex::new();
        let mut ins_index = HammingIndex::new();
        for (i, p) in probes.iter().enumerate() {
            index.insert(&p.seq, i, max_mm);
            for v in single_deletion_variants(&p.seq) {
                del_index.insert(&v, i, 0);
            }
            for v in single_insertion_variants(&p.seq) {
                ins_index.insert(&v, i, 0);
            }
        }
        Self {
            index,
            del_index,
            ins_index,
            probe_len,
            window,
        }
    }

    /// Locate the probe sequence within a window centered on
    /// [`DEFAULT_PROBE_SEQ_OFFSET`], smallest offset first, returning the
    /// matched probe's index into the original `probes` slice. At each
    /// offset, tries the exact-length substitution-tolerant lookup first,
    /// then the single-deletion- and single-insertion-shifted windows.
    pub fn locate(&self, read: &[u8]) -> Option<usize> {
        for k in crate::hamming_index::offset_order(self.window) {
            let start = DEFAULT_PROBE_SEQ_OFFSET as isize + k;
            if start < 0 {
                continue;
            }
            let start = start as usize;
            if let Some(end) = start.checked_add(self.probe_len) {
                if end <= read.len() {
                    if let Some(probe_idx) = self.index.get(&read[start..end]) {
                        return Some(probe_idx);
                    }
                }
            }
            // Deletion-shaped: the read is missing a base the probe has, so
            // the real payload is one bp SHORTER than probe_len here.
            if let Some(end) = start.checked_add(self.probe_len.saturating_sub(1)) {
                if end <= read.len() {
                    if let Some(probe_idx) = self.del_index.get(&read[start..end]) {
                        return Some(probe_idx);
                    }
                }
            }
            // Insertion-shaped: the read carries an extra base the probe
            // doesn't, so the real payload is one bp LONGER than probe_len.
            let end = start + self.probe_len + 1;
            if end <= read.len() {
                if let Some(probe_idx) = self.ins_index.get(&read[start..end]) {
                    return Some(probe_idx);
                }
            }
        }
        None
    }
}

/// Assign Flex GEX reads to genes by probe-set lookup, not genome alignment.
#[derive(Debug, Parser)]
#[command(
    name = "sbxpress-probeset-align",
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
}

fn parse_tag<'a>(comment: &'a str, tag: &str) -> Option<&'a str> {
    let start = comment.find(tag)? + tag.len();
    let rest = &comment[start..];
    let end = rest.find('\t').unwrap_or(rest.len());
    Some(&rest[..end])
}

#[derive(Debug, Serialize)]
struct Metrics {
    reads: u64,
    matched: u64,
    hit_rate_pct: f64,
    distinct_genes_hit: usize,
}

fn genes_manifest_bytes(genes: &[(&str, &str)]) -> Vec<u8> {
    let mut out = Vec::new();
    for (id, name) in genes {
        out.extend_from_slice(id.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(b"\tGene Expression\n");
    }
    out
}

fn establish_genes_manifest(path: &Path, genes: &[(&str, &str)]) -> Result<()> {
    let expected = genes_manifest_bytes(genes);
    if path.exists() {
        let observed = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        ensure!(
            observed == expected,
            "existing gene manifest {} does not match this probe set; use a matching CSV or choose a new path",
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
    let reference = ProbeSetReference::load(&cli.probe_set)?;
    info!(probes = reference.probes.len(), path = %cli.probe_set.display(), "loaded probe set");
    let matcher = ProbeMatcher::build(&reference.probes, cli.max_mm, cli.window);

    // Distinct genes, first-seen order, for the synthetic contig header and
    // the (optional) canonical gene manifest.
    let mut gene_order: Vec<&str> = Vec::new();
    let mut gene_index: HashMap<&str, usize> = HashMap::new();
    let mut gene_name: HashMap<&str, &str> = HashMap::new();
    for p in &reference.probes {
        gene_index.entry(p.gene_id.as_str()).or_insert_with(|| {
            gene_order.push(p.gene_id.as_str());
            gene_order.len() - 1
        });
        gene_name.insert(p.gene_id.as_str(), p.gene_name.as_str());
    }

    if let Some(path) = &cli.genes_out {
        let genes: Vec<(&str, &str)> = gene_order.iter().map(|&id| (id, gene_name[id])).collect();
        establish_genes_manifest(path, &genes)?;
    }

    let mut header = bam::Header::new();
    for &gene_id in &gene_order {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", gene_id)
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
    let mut genes_hit: std::collections::HashSet<&str> = std::collections::HashSet::new();

    let mut fastq = parse_fastx_file(&cli.r2_fastq)
        .with_context(|| format!("opening {}", cli.r2_fastq.display()))?;
    while let Some(rec) = fastq.next() {
        let rec = rec.context("reading R2 FASTQ record")?;
        reads += 1;
        let header_line = String::from_utf8_lossy(rec.id());
        let (qname, comment) = header_line.split_once('\t').unwrap_or((&header_line, ""));
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);

        let mut record = bam::Record::new();
        if let Some(probe_idx) = matcher.locate(&seq) {
            let probe = &reference.probes[probe_idx];
            let gene_id = probe.gene_id.as_str();
            matched += 1;
            genes_hit.insert(gene_id);
            let tid = gene_index[gene_id] as i32;
            let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
            record.set(qname.as_bytes(), Some(&cigar), &seq, qual);
            record.set_tid(tid);
            record.set_pos(0);
            record.set_mapq(60);
            record.unset_unmapped();
            record.push_aux(b"GX", Aux::String(gene_id))?;
            record.push_aux(b"GN", Aux::String(gene_name[gene_id]))?;
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
        distinct_genes_hit: genes_hit.len(),
    };
    if let Some(path) = &cli.metrics {
        let mut w = BufWriter::new(File::create(path)?);
        w.write_all(serde_json::to_string_pretty(&metrics)?.as_bytes())?;
    }

    info!(
        reads = metrics.reads,
        matched = metrics.matched,
        hit_rate_pct = format!("{:.1}", metrics.hit_rate_pct),
        distinct_genes_hit = metrics.distinct_genes_hit,
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

    /// Real header + shape (25bp probes here, not the real 50bp, purely to
    /// keep test literals short -- ProbeSetReference doesn't hardcode a length).
    const REAL_SHAPE: &str = "\
#probe_set_file_format=3.0
#panel_name=test
gene_id,probe_seq,probe_id,included,region,gene_name
ENSG1,AAAACCCCGGGGTTTTACGTACGTA,ENSG1|GENE1|abc,TRUE,spliced,GENE1
ENSG1,CCCCGGGGTTTTAAAAACGTACGT,ENSG1|GENE1|def,FALSE,unspliced,GENE1
ENSG2,TTTTGGGGCCCCAAAATGCATGCAT,ENSG2|GENE2|ghi,TRUE,unspliced,GENE2
";

    #[test]
    fn loads_only_included_true_rows() {
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        assert_eq!(r.probes.len(), 2); // the FALSE row is excluded
        assert_eq!(r.probes[0].gene_id, "ENSG1");
        assert_eq!(r.probes[0].region, "spliced");
        assert_eq!(r.probes[1].gene_id, "ENSG2");
    }

    #[test]
    fn rejects_mixed_probe_lengths() {
        let csv = "\
gene_id,probe_seq,probe_id,included,region,gene_name
ENSG1,AAAACCCCGGGGTTTTACGTACGTA,p1,TRUE,spliced,GENE1
ENSG2,SHORT,p2,TRUE,spliced,GENE2
";
        let f = write_csv(csv);
        assert!(ProbeSetReference::load(f.path()).is_err());
    }

    fn read_with_probe_at(offset: usize, probe: &[u8]) -> Vec<u8> {
        let mut seq = vec![b'N'; offset];
        seq.extend_from_slice(probe);
        seq.extend_from_slice(b"NNNNNNNNNN");
        seq
    }

    #[test]
    fn locates_exact_match_at_confirmed_offset() {
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        let matcher = ProbeMatcher::build(&r.probes, 1, 15);
        let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, b"AAAACCCCGGGGTTTTACGTACGTA");
        assert_eq!(matcher.locate(&read), Some(0));
    }

    #[test]
    fn tolerates_a_mismatch_within_budget() {
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        let matcher = ProbeMatcher::build(&r.probes, 1, 15);
        let mut probe = b"AAAACCCCGGGGTTTTACGTACGTA".to_vec();
        probe[3] = b'G'; // 1 mismatch, within max_mm=1
        let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, &probe);
        assert_eq!(matcher.locate(&read), Some(0));
    }

    #[test]
    fn windowed_search_absorbs_offset_jitter() {
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        let matcher = ProbeMatcher::build(&r.probes, 1, 15);
        // 5bp late: still within the +-15 window. (Not "early" -- offset 0
        // has no room to go earlier; residual extract()-side anchor jitter
        // only ever pushes the probe insert a few bp *later* than 0.)
        let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET + 5, b"AAAACCCCGGGGTTTTACGTACGTA");
        assert_eq!(matcher.locate(&read), Some(0));
    }

    #[test]
    fn distinguishes_two_genes_probes() {
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        let matcher = ProbeMatcher::build(&r.probes, 1, 15);
        let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, b"TTTTGGGGCCCCAAAATGCATGCAT");
        assert_eq!(matcher.locate(&read), Some(1)); // GENE2's probe
    }

    #[test]
    fn recovers_a_single_deletion_in_the_probe_span() {
        // Real SBX error mode (2026-09-18 diagnostic): the read is missing a
        // base the probe has -- del_index should still resolve it, at the
        // *shorter* (probe_len - 1) window.
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        let matcher = ProbeMatcher::build(&r.probes, 1, 15);
        let probe = b"AAAACCCCGGGGTTTTACGTACGTA";
        let mut deleted = probe[..10].to_vec();
        deleted.extend_from_slice(&probe[11..]); // drop the base at position 10
        assert_eq!(deleted.len(), probe.len() - 1);
        let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, &deleted);
        assert_eq!(matcher.locate(&read), Some(0));
    }

    #[test]
    fn recovers_a_single_insertion_in_the_probe_span() {
        // The rarer real shape (9.2% of recoverable indel-shaped misses vs
        // 85.2% deletion) but still real -- ins_index resolves it at the
        // *longer* (probe_len + 1) window.
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        let matcher = ProbeMatcher::build(&r.probes, 1, 15);
        let probe = b"AAAACCCCGGGGTTTTACGTACGTA";
        let mut inserted = probe[..10].to_vec();
        inserted.push(b'T'); // extra base not present in the real probe
        inserted.extend_from_slice(&probe[10..]);
        assert_eq!(inserted.len(), probe.len() + 1);
        let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, &inserted);
        assert_eq!(matcher.locate(&read), Some(0));
    }

    #[test]
    fn indel_recovery_still_distinguishes_two_genes() {
        // A deletion shouldn't blur two real, distinct probes together.
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        let matcher = ProbeMatcher::build(&r.probes, 1, 15);
        let probe2 = b"TTTTGGGGCCCCAAAATGCATGCAT";
        let mut deleted = probe2[..5].to_vec();
        deleted.extend_from_slice(&probe2[6..]);
        let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, &deleted);
        assert_eq!(matcher.locate(&read), Some(1)); // GENE2's probe, not GENE1's
    }

    #[test]
    fn no_match_on_unrelated_sequence() {
        let f = write_csv(REAL_SHAPE);
        let r = ProbeSetReference::load(f.path()).unwrap();
        let matcher = ProbeMatcher::build(&r.probes, 1, 15);
        let read = vec![b'A'; 100];
        assert_eq!(matcher.locate(&read), None);
    }

    #[test]
    fn genes_manifest_create_then_verify_round_trips() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("genes.tsv");
        let genes = vec![("ENSG1", "GENE1"), ("ENSG2", "GENE2")];
        establish_genes_manifest(&out, &genes).unwrap();
        let text = fs::read_to_string(&out).unwrap();
        assert_eq!(
            text,
            "ENSG1\tGENE1\tGene Expression\nENSG2\tGENE2\tGene Expression\n"
        );
        establish_genes_manifest(&out, &genes).unwrap(); // idempotent
    }

    #[test]
    fn genes_manifest_rejects_mismatched_existing_file() {
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("genes.tsv");
        fs::write(&out, "stale\tSTALE\tGene Expression\n").unwrap();
        assert!(establish_genes_manifest(&out, &[("ENSG1", "GENE1")]).is_err());
    }

    #[test]
    fn run_end_to_end_writes_mapped_and_unmapped_records_with_carried_tags() {
        let probe_csv = write_csv(REAL_SHAPE);
        let dir = tempfile::tempdir().unwrap();
        let fastq_path = dir.path().join("r2.fastq");
        let matched_seq =
            read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, b"AAAACCCCGGGGTTTTACGTACGTA");
        let unmatched_seq = vec![b'A'; 120];
        let mut fastq = File::create(&fastq_path).unwrap();
        writeln!(
            fastq,
            "@read1\tCR:Z:ACGTACGTACGTACGT\tCB:Z:ACGTACGTACGTACGT-1\tUR:Z:AAAACCCCGGGG\tUB:Z:AAAACCCCGGGG\tPB:Z:D-A01"
        )
        .unwrap();
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

        let out_bam = dir.path().join("probeset.bam");
        let genes_out = dir.path().join("genes.tsv");
        let metrics_path = dir.path().join("metrics.json");
        let cli = Cli {
            r2_fastq: fastq_path,
            probe_set: probe_csv.path().to_path_buf(),
            output: out_bam.clone(),
            genes_out: Some(genes_out.clone()),
            metrics: Some(metrics_path.clone()),
            max_mm: 1,
            window: 15,
        };
        run(&cli).unwrap();

        assert!(genes_out.exists());
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
        assert_eq!(mapped.aux(b"GX").unwrap(), Aux::String("ENSG1"));
        assert_eq!(mapped.aux(b"GN").unwrap(), Aux::String("GENE1"));
        assert_eq!(mapped.aux(b"PB").unwrap(), Aux::String("D-A01"));

        let unmapped = &records[1];
        assert!(unmapped.is_unmapped());
        assert!(unmapped.aux(b"GX").is_err());
        assert!(unmapped.aux(b"PB").is_err()); // PB optional, absent here
    }
}
