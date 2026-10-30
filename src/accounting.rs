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
    /// Of `per_sample_early`, how many also located a real R2/cDNA 3' boundary
    /// (adapter or Flex anchor found -- see `Extracted::adapter_trimmed`),
    /// excluding reads too short/garbled to ever locate one regardless of
    /// matching-algorithm quality. The fair denominator for scoring a
    /// downstream matcher (probe-set, guide capture, ...) in isolation from
    /// upstream capture-length loss.
    pub per_sample_early_boundary_found: Vec<SampleCount>,
    /// Per-index evidence for every read whose R1 was located (the early
    /// demux population): which index regions are in the read at all and
    /// which matched a sample.
    pub index_evidence: IndexEvidenceCounts,
    /// One-indel index matching (`--index-indels`) active, per index.
    pub index_indels_i5: bool,
    pub index_indels_i7: bool,
}

/// i5/i7 presence and identification counts (early demux, R1 located).
///
/// "Seen" = the index region lies inside the read (i5: R1 located far enough
/// in; i7: 3' adapter located and the read runs past it). "Identified" = it
/// matched a sample. `full` = both seen and both identified to the same
/// sample (= `via_dual` early); every other R1-located read is `partial`.
/// `identified_both` includes conflicts (both matched, different samples).
#[derive(Serialize, Default, Clone, Debug, PartialEq, Eq)]
pub struct IndexEvidenceCounts {
    pub reads: u64,
    pub full: u64,
    pub partial: u64,
    pub seen_both: u64,
    pub seen_i5_only: u64,
    pub seen_i7_only: u64,
    pub seen_none: u64,
    pub identified_both: u64,
    pub identified_i5_only: u64,
    pub identified_i7_only: u64,
    pub identified_none: u64,
    /// Identified only through the one-indel pass (`--index-indels`).
    pub identified_i5_via_indel: u64,
    pub identified_i7_via_indel: u64,
}

impl IndexEvidenceCounts {
    /// Count one R1-located read.
    pub fn observe(&mut self, ev: &crate::demux::IndexEvidence, dual: bool) {
        self.reads += 1;
        if dual {
            self.full += 1;
        } else {
            self.partial += 1;
        }
        match (ev.i5_seen, ev.i7_seen) {
            (true, true) => self.seen_both += 1,
            (true, false) => self.seen_i5_only += 1,
            (false, true) => self.seen_i7_only += 1,
            (false, false) => self.seen_none += 1,
        }
        match (ev.i5.sample.is_some(), ev.i7.sample.is_some()) {
            (true, true) => self.identified_both += 1,
            (true, false) => self.identified_i5_only += 1,
            (false, true) => self.identified_i7_only += 1,
            (false, false) => self.identified_none += 1,
        }
        self.identified_i5_via_indel += u64::from(ev.i5.via_indel);
        self.identified_i7_via_indel += u64::from(ev.i7.via_indel);
    }

    pub fn add(&mut self, o: &Self) {
        self.reads += o.reads;
        self.full += o.full;
        self.partial += o.partial;
        self.seen_both += o.seen_both;
        self.seen_i5_only += o.seen_i5_only;
        self.seen_i7_only += o.seen_i7_only;
        self.seen_none += o.seen_none;
        self.identified_both += o.identified_both;
        self.identified_i5_only += o.identified_i5_only;
        self.identified_i7_only += o.identified_i7_only;
        self.identified_none += o.identified_none;
        self.identified_i5_via_indel += o.identified_i5_via_indel;
        self.identified_i7_via_indel += o.identified_i7_via_indel;
    }
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

    /// Pairs whose cDNA failed the low-complexity check (`YL:i:1`, or
    /// failed as `low_complexity` with `drop_low_complexity`).
    pub low_complexity: u64,
    /// `low_complexity` by the read's most frequent base.
    pub low_complexity_by_base: crate::complexity::ByBase,
    pub low_complexity_base_fraction: f64,
    pub low_complexity_run_fraction: f64,
    pub drop_low_complexity: bool,

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
    /// Rescued by `--r1-tail-relock` (XB:i:7).
    pub relocked: u64,
    pub corrected_hamming2: u64,
}

/// Flex probe-barcode location/identification metrics (the Flex design notes
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
    /// Of `hit`, how many also landed on an assigned (non-Undetermined)
    /// sample from the outer index-based demux -- the real 4-way join (CB +
    /// UMI + Sample Barcode + Probe Assigned), a tighter "pre-dedup" proxy
    /// than `hit` alone (which doesn't require the outer sample/library
    /// demux to have succeeded). `0` when no `--bclconvert-samplesheet` was
    /// given (no outer demux to join against).
    pub hit_sample_assigned: u64,
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
