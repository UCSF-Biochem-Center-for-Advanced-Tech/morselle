//! Run summary: metrics JSON assembly and the human-readable report.
use tracing::info;

use crate::accounting::{
    pct, DemuxMetrics, Metrics, ProbeBarcodeMetrics, SampleCount, WhitelistMetrics,
};
use crate::barcode::Whitelist;
use crate::demux::Demuxer;
use crate::probe_barcode::ProbeWhitelist;

use super::types::{Stats, Timing};

pub(crate) fn build_metrics(
    stats: &Stats,
    whitelist: Option<&Whitelist>,
    demuxer: Option<&Demuxer>,
    probe_wl: Option<&ProbeWhitelist>,
    t: &Timing,
    minimap2: bool,
    processing_tags: bool,
) -> Metrics {
    let r2_len_mean = if stats.pairs_emitted > 0 {
        stats.r2_len_sum / stats.pairs_emitted
    } else {
        0
    };
    let reads_per_sec = if t.elapsed > 0.0 {
        (stats.reads as f64 / t.elapsed) as u64
    } else {
        0
    };
    Metrics {
        output_mode: if minimap2 {
            "minimap2_r2_only"
        } else {
            "paired_fastq"
        }
        .into(),
        processing_tags,
        total_reads: stats.reads,
        reads_failed: stats.failed_emitted,
        elapsed_secs: t.elapsed,
        load_secs: t.load_secs,
        process_secs: t.process_secs,
        reads_per_sec,
        threads: t.threads,
        oriented: stats.oriented,
        oriented_pct: pct(stats.oriented, stats.reads),
        class_full: stats.class_full,
        class_partial: stats.class_partial,
        class_unparseable: stats.class_unparseable,
        pairs_emitted: stats.pairs_emitted,
        pairs_pct: pct(stats.pairs_emitted, stats.reads),
        r1_eq_28: stats.r1_eq_28,
        r1_not_28: stats.r1_not_28,
        adapter_trimmed: stats.adapter_trimmed,
        r2_len_mean,
        r2_len_max: stats.r2_len_max,
        // Set by the caller, which holds the CLI thresholds.
        low_complexity: 0,
        low_complexity_by_base: Default::default(),
        low_complexity_base_fraction: 0.0,
        low_complexity_run_fraction: 0.0,
        drop_low_complexity: false,
        whitelist: whitelist.map(|w| {
            let classified = stats.cb_on_whitelist + stats.cb_within1 + stats.cb_fail;
            WhitelistMetrics {
                size: w.len(),
                cb_on_whitelist: stats.cb_on_whitelist,
                cb_within1: stats.cb_within1,
                cb_fail: stats.cb_fail,
                cb_hit_rate: pct(stats.cb_on_whitelist + stats.cb_within1, classified),
                cb_located_confident: stats.cb_located_confident,
                cb_located_confident_pct: pct(stats.cb_located_confident, classified),
                cb_fail_confident: stats.cb_fail_confident,
                cb_fail_confident_pct: pct(stats.cb_fail_confident, stats.cb_fail),
                gated_out: stats.gated_out,
                realigned: stats.realigned,
                reframed: stats.reframed,
                corrected_hamming2: stats.corrected2,
            }
        }),
        demux: demuxer.map(|d| {
            let assigned = stats.pairs_emitted - stats.undetermined;
            let assigned_early: u64 = stats.per_sample_early.iter().sum();
            DemuxMetrics {
                assigned,
                assigned_pct: pct(assigned, stats.pairs_emitted),
                undetermined: stats.undetermined,
                via_dual: stats.via_dual,
                via_i5_only: stats.via_i5,
                via_i7_only: stats.via_i7,
                via_conflict: stats.via_conflict,
                per_sample: d
                    .sample_ids()
                    .iter()
                    .enumerate()
                    .map(|(i, id)| SampleCount {
                        sample: id.clone(),
                        pairs: stats.per_sample[i],
                    })
                    .collect(),
                assigned_early,
                assigned_early_pct: pct(assigned_early, stats.reads),
                undetermined_early: stats.undetermined_early,
                index_unreachable: stats.index_unreachable,
                per_sample_early: d
                    .sample_ids()
                    .iter()
                    .enumerate()
                    .map(|(i, id)| SampleCount {
                        sample: id.clone(),
                        pairs: stats.per_sample_early[i],
                    })
                    .collect(),
                per_sample_early_boundary_found: d
                    .sample_ids()
                    .iter()
                    .enumerate()
                    .map(|(i, id)| SampleCount {
                        sample: id.clone(),
                        pairs: stats.per_sample_early_boundary_found[i],
                    })
                    .collect(),
            }
        }),
        probe_barcode: probe_wl.map(|w| ProbeBarcodeMetrics {
            size: w.len(),
            located: stats.probe_located,
            hit: stats.probe_hit,
            hit_rate: pct(stats.probe_hit, stats.probe_located),
            hit_sample_assigned: stats.probe_hit_sample_assigned,
        }),
    }
}

/// Log the run summary.
pub(crate) fn report(stats: &Stats, demuxer: Option<&Demuxer>, t: &Timing) {
    let p = |x: u64| 100.0 * x as f64 / stats.reads.max(1) as f64;
    info!(
        reads = stats.reads,
        oriented_pct = format!("{:.1}", p(stats.oriented)),
        full = stats.class_full,
        partial = stats.class_partial,
        unparseable = stats.class_unparseable,
        pairs = stats.pairs_emitted,
        pairs_pct = format!("{:.1}", p(stats.pairs_emitted)),
        r2_len_mean = (stats.pairs_emitted > 0).then(|| stats.r2_len_sum / stats.pairs_emitted),
        r2_len_max = stats.r2_len_max,
        "pipeline complete"
    );
    if stats.cb_on_whitelist + stats.cb_within1 + stats.cb_fail > 0 {
        let classified = stats.cb_on_whitelist + stats.cb_within1 + stats.cb_fail;
        info!(
            on_whitelist = stats.cb_on_whitelist,
            within1 = stats.cb_within1,
            fail = stats.cb_fail,
            hit_rate = format!(
                "{:.1}",
                100.0 * (stats.cb_on_whitelist + stats.cb_within1) as f64
                    / classified.max(1) as f64
            ),
            located_confident_pct = format!(
                "{:.1}",
                100.0 * stats.cb_located_confident as f64 / classified.max(1) as f64
            ),
            fail_confident_pct = format!(
                "{:.1}",
                100.0 * stats.cb_fail_confident as f64 / stats.cb_fail.max(1) as f64
            ),
            realigned = stats.realigned,
            reframed = stats.reframed,
            corrected2 = stats.corrected2,
            gated_out = stats.gated_out,
            "CB gate"
        );
    }
    if stats.probe_located > 0 {
        info!(
            located = stats.probe_located,
            hit = stats.probe_hit,
            hit_rate = format!(
                "{:.1}",
                100.0 * stats.probe_hit as f64 / stats.probe_located.max(1) as f64
            ),
            hit_sample_assigned = stats.probe_hit_sample_assigned,
            "probe barcode"
        );
    }
    if let Some(d) = demuxer {
        let assigned = stats.pairs_emitted - stats.undetermined;
        info!(
            assigned,
            assigned_pct = format!("{:.1}", p(assigned)),
            undetermined = stats.undetermined,
            via_dual = stats.via_dual,
            via_i5_only = stats.via_i5,
            via_i7_only = stats.via_i7,
            via_conflict = stats.via_conflict,
            "demux summary"
        );
        for (i, id) in d.sample_ids().iter().enumerate() {
            info!(sample = %id, pairs = stats.per_sample[i], "  per-sample");
        }
    }
    if stats.r1_not_28 > 0 {
        tracing::warn!(count = stats.r1_not_28, "R1 records not exactly 28 bp");
    }

    // Headline run summary: time, throughput, and where reads ended up.
    let has_demux = demuxer.is_some();
    let assigned = if has_demux {
        stats.pairs_emitted.saturating_sub(stats.undetermined)
    } else {
        stats.pairs_emitted
    };
    let final_ok = if has_demux {
        assigned
    } else {
        stats.pairs_emitted
    };
    let lost = stats.reads.saturating_sub(final_ok);
    let no_cdna = stats.class_partial.saturating_sub(stats.gated_out);
    let rps = if t.elapsed > 0.0 {
        (stats.reads as f64 / t.elapsed) as u64
    } else {
        0
    };
    info!(
        elapsed = format!("{:.1}s ({:.1} min)", t.elapsed, t.elapsed / 60.0),
        load_secs = format!("{:.1}", t.load_secs),
        process_secs = format!("{:.1}", t.process_secs),
        threads = t.threads,
        reads_per_sec = rps,
        processed = stats.reads,
        emitted = stats.pairs_emitted,
        assigned = if has_demux { Some(assigned) } else { None },
        lost,
        lost_pct = format!("{:.1}", p(lost)),
        "RUN SUMMARY"
    );
    info!(
        unparseable = stats.class_unparseable,
        no_cdna,
        cb_fail_gated = stats.gated_out,
        index_undetermined = has_demux.then_some(stats.undetermined),
        "  lost breakdown"
    );
    if has_demux {
        let assigned_early: u64 = stats.per_sample_early.iter().sum();
        info!(
            assigned_early,
            assigned_early_pct = format!("{:.1}", 100.0 * assigned_early as f64 / stats.reads.max(1) as f64),
            undetermined_early = stats.undetermined_early,
            index_unreachable = stats.index_unreachable,
            "  index demux (early, pre-gate -- comparable to Illumina's own post-BCL-Convert demux count)"
        );
        info!(
            via_dual = stats.via_dual_early,
            via_i5_only = stats.via_i5_only_early,
            via_i7_only = stats.via_i7_only_early,
            via_conflict = stats.via_conflict_early,
            "    index demux (early) via breakdown"
        );
    }
}
