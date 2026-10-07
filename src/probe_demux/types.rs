use clap::Parser;
use serde::Serialize;
use std::path::PathBuf;

/// Assign cell barcodes to samples by probe barcode (10x Flex).
#[derive(Debug, Parser)]
#[command(
    name = "sbxpress-probe-demux",
    version,
    about = "Assign cell barcodes to samples by probe barcode (10x Flex)"
)]
pub struct Cli {
    /// GEX-library BAM after `annotate`+`molecules` (must carry `CB`, `PB`,
    /// and the `xf` representative-molecule bit).
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,

    /// Probe-barcode sample sheet: a real cellranger `multi` `config.csv`
    /// (only `[samples]` is read) or a flat `sample_id,probe_barcode_ids` CSV.
    #[arg(long)]
    pub probe_sample_sheet: PathBuf,

    /// Output directory.
    #[arg(short, long, default_value = "outs")]
    pub output: PathBuf,

    /// Minimum molecules a sample's declared probe-barcode set must capture
    /// for a GEM to be Assigned to that sample. A GEM can independently clear
    /// this floor for more than one sample — see [`assign`]'s module-level
    /// doc comment on composite `(GEM, ProbeBarcode)` cell barcodes.
    #[arg(long, default_value_t = 10)]
    pub min_molecules: u64,
}

/// Tuning for [`assign`].
#[derive(Debug, Clone, Copy)]
pub struct DemuxCfg {
    pub min_molecules: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Status {
    /// This GEM independently cleared `min_molecules` for `sample_id`'s
    /// declared set. Ordinary, expected outcome — including when another
    /// sample also has its own `Assigned` row for the same `cell_barcode`
    /// (real Flex multiplex pooling: see the module doc comment).
    Assigned,
    /// Below `min_molecules` total, or no sample's declared set captured
    /// enough of it to call — low-confidence/noisy.
    Unassigned,
    /// This GEM's number of independently-clearing samples exceeds the
    /// negative-binomial 99.9th-percentile cutoff fit to the whole run's
    /// observed per-GEM occupancy distribution — the rare chip-loading-
    /// artifact tail 10x's own "High occupancy GEM filtering" step targets,
    /// not ordinary multi-sample pooling. All of this GEM's molecules are
    /// excluded from every sample. See [`gem_overlap`] for the
    /// whole-experiment version of the underlying overlap diagnostic.
    HighOccupancyGem,
}

impl Status {
    pub fn as_str(&self) -> &'static str {
        match self {
            Status::Assigned => "Assigned",
            Status::Unassigned => "Unassigned",
            Status::HighOccupancyGem => "HighOccupancyGem",
        }
    }
}

/// One cell barcode's sample-assignment result. A single `cell_barcode` can
/// produce more than one `Assigned` row (one per sample whose declared set
/// independently clears the floor) — see the module doc comment.
#[derive(Debug, Clone)]
pub struct Assignment {
    pub cell_barcode: String,
    /// `Some` only when `status == Assigned`.
    pub sample_id: Option<String>,
    /// Distinct probe-barcode ids observed for this CB, sorted.
    pub probe_barcodes: Vec<String>,
    /// This row's own molecule count: `sample_id`'s declared-set sum when
    /// `Assigned`, the whole GEM's total otherwise.
    pub molecules: u64,
    /// This row's sample's molecules / the whole GEM's total molecules
    /// (dominant-sample fraction when not `Assigned`).
    pub purity: f64,
    /// Number of samples whose declared set independently cleared
    /// `min_molecules` for this GEM (10x's own "Probe barcodes per GEM").
    pub occupancy: usize,
    pub status: Status,
}

/// Whole-experiment probe-barcode-pair co-occurrence diagnostic, mirroring
/// cellranger's real `frp_gem_barcode_overlap.csv`: for every pair of probe
/// barcodes with at least one shared cell barcode, `overlap =
/// common_gems / min(barcode1_gems, barcode2_gems)`.
#[derive(Debug, Clone, Serialize)]
pub struct OverlapRow {
    pub barcode1: String,
    pub barcode2: String,
    pub barcode1_gems: usize,
    pub barcode2_gems: usize,
    pub common_gems: usize,
    pub overlap: f64,
}

/// Run-level metrics, written to `probe_demux_metrics.json`.
#[derive(Debug, Serialize)]
pub struct Metrics {
    pub elapsed_secs: f64,
    /// Distinct GEM barcodes seen (not row count — a GEM legitimately
    /// present in more than one sample produces more than one row).
    pub cell_barcodes: usize,
    /// `Assigned` rows (one GEM present in 2 samples counts as 2).
    pub assigned: usize,
    pub unassigned: usize,
    /// Distinct GEMs discarded by the negative-binomial high-occupancy-GEM
    /// tail filter — see the module doc comment. Expect this to be small; a
    /// large value here (or `occupancy_cutoff` far below the sheet's sample
    /// count) means the run doesn't fit 10x's own "usually >90% of cells
    /// pass" expectation and is worth a closer look.
    pub high_occupancy_gems: usize,
    /// NB(r, p=0.5) fit's `r` (the observed mean occupancy across occupied
    /// GEMs) and the resulting 99.9th-percentile cutoff actually applied.
    pub occupancy_r: f64,
    pub occupancy_cutoff: u64,
    pub per_sample: Vec<SampleCount>,
}

#[derive(Debug, Serialize)]
pub struct SampleCount {
    pub sample: String,
    pub cell_barcodes: usize,
}
