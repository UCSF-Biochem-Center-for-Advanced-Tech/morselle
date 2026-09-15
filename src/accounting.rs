//! Summary JSON metrics (plan §6, M4).
//!
//! The per-read TSV stream is deferred; M4 emits `metrics_summary.json`.

use std::path::Path;

use anyhow::{Context, Result};
use serde::Serialize;

/// Per-sample assigned-pair count.
#[derive(Serialize)]
pub struct SampleCount {
    pub sample: String,
    pub pairs: u64,
}

/// Demultiplexing breakdown (present only when a sample sheet was given).
#[derive(Serialize)]
pub struct DemuxMetrics {
    pub assigned: u64,
    pub assigned_pct: f64,
    pub undetermined: u64,
    pub via_dual: u64,
    pub via_i5_only: u64,
    pub via_i7_only: u64,
    pub via_conflict: u64,
    pub per_sample: Vec<SampleCount>,
    /// Early (pre-gate) demux: assigned the moment R1 was located, independent
    /// of CB/cDNA outcome -- the closest SBX-side analog to Illumina's own
    /// post-BCL-Convert index-only demultiplex count. Denominator is
    /// `total_reads` (matching how Illumina reports its own demux rate,
    /// against all reads sequenced), not `pairs_emitted`.
    pub assigned_early: u64,
    pub assigned_early_pct: f64,
    pub undetermined_early: u64,
    /// R1 itself was never located -- demux is structurally unreachable for
    /// these (SBX has no separate, fixed-position index read the way
    /// Illumina does; this is the real chemistry floor, not a code gap).
    pub index_unreachable: u64,
    pub per_sample_early: Vec<SampleCount>,
}

/// Top-level run metrics, serialized to `metrics_summary.json`.
#[derive(Serialize)]
pub struct Metrics {
    /// Reads written to the `--failed-output` side channel. When that option is
    /// set, `total_reads == pairs_emitted + reads_failed` for every input.
    pub reads_failed: u64,
    /// FASTQ layout selected for this run.
    pub output_mode: String,
    /// Whether SBX per-read processing provenance tags were emitted.
    pub processing_tags: bool,
    pub total_reads: u64,
    /// Wall-clock seconds for the whole run (load + process).
    pub elapsed_secs: f64,
    /// One-time setup: sample sheet, file discovery, whitelist load, target creation.
    pub load_secs: f64,
    /// Per-read parallel processing (orient/extract/demux/gate + single-member write).
    pub process_secs: f64,
    /// Throughput = total_reads / elapsed_secs.
    pub reads_per_sec: u64,
    /// Worker threads used.
    pub threads: usize,
    pub oriented: u64,
    pub oriented_pct: f64,

    pub class_full: u64,
    pub class_partial: u64,
    pub class_unparseable: u64,

    pub pairs_emitted: u64,
    pub pairs_pct: f64,
    pub r1_eq_28: u64,
    pub r1_not_28: u64,
    pub adapter_trimmed: u64,
    pub r2_len_mean: u64,
    pub r2_len_max: u64,

    /// `null` when no whitelist was provided.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub whitelist: Option<WhitelistMetrics>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub demux: Option<DemuxMetrics>,

    /// `null` unless `--probe-barcode-whitelist` was given.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub probe_barcode: Option<ProbeBarcodeMetrics>,
}

/// Cell-barcode gate metrics (present only when a whitelist was given).
#[derive(Serialize)]
pub struct WhitelistMetrics {
    pub size: usize,
    pub cb_on_whitelist: u64,
    pub cb_within1: u64,
    pub cb_fail: u64,
    /// (on_whitelist + within1) / classified.
    pub cb_hit_rate: f64,
    /// CBs whose location is dual-anchor confirmed (R1-tail + polyT flanks).
    pub cb_located_confident: u64,
    /// Fraction of classified CBs that are confidently located.
    pub cb_located_confident_pct: f64,
    /// CB-fails that are confidently located = genuine base errors (not mislocated).
    pub cb_fail_confident: u64,
    /// Fraction of CB-fails that are confidently located.
    pub cb_fail_confident_pct: f64,
    pub gated_out: u64,
    pub realigned: u64,
    pub reframed: u64,
    pub corrected_hamming2: u64,
}

/// Flex probe-barcode location/identification metrics (FLEX_SUPPORT_PLAN.md
/// Phase 1; present only when `--probe-barcode-whitelist` was given).
#[derive(Serialize)]
pub struct ProbeBarcodeMetrics {
    pub size: usize,
    /// Pairs where the probe-barcode window was located at all (in-bounds).
    pub located: u64,
    /// Pairs where the located window matched a whitelist entry confidently.
    pub hit: u64,
    /// hit / located.
    pub hit_rate: f64,
}

impl Metrics {
    pub fn write_json(&self, path: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(self)?;
        std::fs::write(path, json)
            .with_context(|| format!("cannot write metrics {}", path.display()))
    }
}

/// Percentage helper guarding divide-by-zero.
pub fn pct(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        (1000.0 * num as f64 / den as f64).round() / 10.0
    }
}
