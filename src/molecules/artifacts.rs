//! Per-gene artifact report (plan Layer D; no reads dropped).
//!
//! DZEC05's dead read 2 gave LINC00486 18-21% of all molecules before anyone
//! noticed. This tallies, for every gene-assigned primary record, whether its
//! cDNA was flagged low-complexity at conversion (`YL:i:1`) and whether it is
//! below MAPQ 60 without annotate's `mm:i:1` rescue -- the reads Cell
//! Ranger's rule would not count, measured whatever `--min-mapq` this run
//! used, so an unfiltered run still names the gene. Genes where most reads
//! are flagged are listed in `molecules.metrics.json` under `artifacts`, and
//! a large one is printed as a warning.
//!
//! In a normal library a few intronic-repeat genes are listed (DZEC04:
//! AMBRA1, NBEAL1, PLEKHA7 at 34-88% low-MAPQ), each well under 1% of reads;
//! the warning is for a listed gene holding a real share of the library.
use rust_htslib::bam;
use rustc_hash::FxHashMap as HashMap;
use serde::Serialize;

use super::helpers::MoleculeTags;

/// Below this MAPQ, and without `mm:i:1`, a read counts as low-MAPQ here
/// (minimap2's confident value; Cell Ranger's rule).
pub const CONFIDENT_MAPQ: u8 = 60;
/// A gene is listed when it has at least this many reads ...
pub const MIN_READS: u64 = 100;
/// ... and at least this fraction of them are low-complexity or low-MAPQ.
pub const MIN_FLAGGED_FRACTION: f64 = 0.5;
/// Listed genes, largest flagged count first.
pub const MAX_LISTED: usize = 25;
/// Warn when a listed gene holds at least this share of gene-assigned reads.
pub const WARN_SHARE: f64 = 0.01;

#[derive(Clone, Default)]
struct Tally {
    name: String,
    reads: u64,
    low_complexity: u64,
    low_mapq: u64,
    flagged: u64,
}

/// Per-gene counts gathered during evidence collection.
#[derive(Clone, Default)]
pub struct GeneTally(HashMap<String, Tally>);

impl GeneTally {
    pub fn observe(&mut self, gene: &str, name: &str, record: &bam::Record, tags: &MoleculeTags) {
        let low_complexity = tags.yl == Some(1);
        let low_mapq = record.mapq() < CONFIDENT_MAPQ && tags.mm != Some(1);
        let tally = match self.0.get_mut(gene) {
            Some(tally) => tally,
            None => self.0.entry(gene.to_owned()).or_insert_with(|| Tally {
                name: name.to_owned(),
                ..Tally::default()
            }),
        };
        tally.reads += 1;
        tally.low_complexity += low_complexity as u64;
        tally.low_mapq += low_mapq as u64;
        tally.flagged += (low_complexity || low_mapq) as u64;
    }

    pub fn merge_from(&mut self, other: &GeneTally) {
        for (gene, source) in &other.0 {
            let target = self.0.entry(gene.clone()).or_insert_with(|| Tally {
                name: source.name.clone(),
                ..Tally::default()
            });
            target.reads += source.reads;
            target.low_complexity += source.low_complexity;
            target.low_mapq += source.low_mapq;
            target.flagged += source.flagged;
        }
    }

    pub fn report(&self) -> ArtifactReport {
        let mut report = ArtifactReport::default();
        for tally in self.0.values() {
            report.gene_assigned_reads += tally.reads;
            report.low_complexity_reads += tally.low_complexity;
            report.low_mapq_reads += tally.low_mapq;
        }
        let total = report.gene_assigned_reads.max(1) as f64;
        let mut listed = self
            .0
            .iter()
            .filter(|(_, t)| {
                t.reads >= MIN_READS && t.flagged as f64 >= MIN_FLAGGED_FRACTION * t.reads as f64
            })
            .map(|(gene, t)| ArtifactGene {
                gene_id: gene.clone(),
                gene_name: t.name.clone(),
                reads: t.reads,
                low_complexity_reads: t.low_complexity,
                low_mapq_reads: t.low_mapq,
                flagged_fraction: t.flagged as f64 / t.reads as f64,
                share_of_reads: t.reads as f64 / total,
            })
            .collect::<Vec<_>>();
        listed.sort_by(|a, b| {
            let flagged = |g: &ArtifactGene| (g.flagged_fraction * g.reads as f64).round() as u64;
            flagged(b)
                .cmp(&flagged(a))
                .then_with(|| a.gene_id.cmp(&b.gene_id))
        });
        listed.truncate(MAX_LISTED);
        report.genes = listed;
        report
    }
}

/// `artifacts` in `molecules.metrics.json`.
#[derive(Debug, Default, Serialize)]
pub struct ArtifactReport {
    /// Gene-assigned primary records with CB and UMI (before `--min-mapq`).
    pub gene_assigned_reads: u64,
    /// Of those, cDNA flagged low-complexity at conversion (`YL:i:1`).
    pub low_complexity_reads: u64,
    /// Of those, MAPQ < 60 without `mm:i:1`.
    pub low_mapq_reads: u64,
    /// Genes with >= 100 reads, at least half flagged; largest first.
    pub genes: Vec<ArtifactGene>,
}

#[derive(Debug, Serialize)]
pub struct ArtifactGene {
    pub gene_id: String,
    pub gene_name: String,
    pub reads: u64,
    pub low_complexity_reads: u64,
    pub low_mapq_reads: u64,
    /// Share of this gene's reads that are low-complexity or low-MAPQ.
    pub flagged_fraction: f64,
    /// This gene's share of all gene-assigned reads.
    pub share_of_reads: f64,
}

impl ArtifactReport {
    /// Print a warning for each listed gene holding >= [`WARN_SHARE`] of reads.
    pub fn warn(&self) {
        for gene in self.genes.iter().filter(|g| g.share_of_reads >= WARN_SHARE) {
            eprintln!(
                "[sbxpress molecules] WARNING: {} ({}) has {:.1}% of gene-assigned reads and \
                 {:.0}% of them are low-complexity ({}) or low-MAPQ ({}); likely a library \
                 artifact (see `artifacts` in the metrics JSON)",
                gene.gene_name,
                gene.gene_id,
                100.0 * gene.share_of_reads,
                100.0 * gene.flagged_fraction,
                gene.low_complexity_reads,
                gene.low_mapq_reads
            );
        }
    }
}
