use std::sync::atomic::AtomicUsize;
use std::time::Instant;

use anyhow::Result;
use rayon::prelude::*;
use tracing::info;

use crate::accounting::{
    pct, DemuxMetrics, Metrics, ProbeBarcodeMetrics, SampleCount, WhitelistMetrics,
};
use crate::barcode::Whitelist;
use crate::config::Cli;
use crate::demux::Demuxer;
use crate::io::{discover_fastqs, PairShardWriter, ShardWriter};
use crate::probe_barcode::ProbeWhitelist;
use crate::sample_sheet::SampleSheet;

use super::types::{Stats, Timing};
use super::worker::{process_file, run_pipelined, Shared};

/// Entry point invoked by the binary.
pub fn run(cli: &Cli) -> Result<()> {
    if cli.params.cellranger && cli.params.minimap2 {
        anyhow::bail!("--cellranger and --minimap2 are mutually exclusive output presets");
    }
    let t_start = Instant::now();
    let sheet = match &cli.bclconvert_samplesheet {
        Some(path) => {
            let sheet = SampleSheet::load(path)?;
            info!(
                samples = sheet.samples.len(),
                i5_revcomp = sheet.i5_reverse_complemented,
                chemistry = ?sheet.chemistry,
                "loaded BCL Convert sample sheet"
            );
            if sheet.source == crate::sample_sheet::SourceFormat::BclConvert {
                let flat_path = cli.out_dir.join("bclconvert_samplesheet.csv");
                sheet.write_flat_csv(&flat_path)?;
                info!(out = %flat_path.display(), "wrote normalized flat sample sheet");
            }
            Some(sheet)
        }
        None => None,
    };

    let inputs = discover_fastqs(&cli.fastqs)?;
    info!(n_files = inputs.len(), "discovered FASTQ inputs");
    if inputs.is_empty() {
        anyhow::bail!("no FASTQ inputs found");
    }

    let orient_cfg = cli.params.orient_cfg();
    let extract_cfg = cli.params.extract_cfg();
    let cb_len = cli.params.cb_len;
    let umi_len = cli.params.umi_len;

    let whitelist = match &cli.whitelist {
        Some(path) => {
            let w = Whitelist::load(path)?;
            info!(barcodes = w.len(), path = %path.display(), "loaded CB whitelist");
            Some(w)
        }
        None => None,
    };

    let probe_wl = match &cli.params.probe_barcode_whitelist {
        Some(path) => {
            let w = ProbeWhitelist::load(
                path,
                cli.params.probe_barcode_max_mm,
                cli.params.probe_barcode_window,
            )?;
            info!(barcodes = w.len(), path = %path.display(), "loaded probe-barcode whitelist");
            Some(w)
        }
        None => None,
    };

    let demuxer = sheet
        .as_ref()
        .map(|s| Demuxer::new(s, cli.params.demux_cfg()))
        .transpose()?;
    let n_samples = demuxer.as_ref().map_or(1, |d| d.sample_ids().len());
    let gex_by_sample: Vec<bool> = demuxer
        .as_ref()
        .map(|d| {
            d.sample_ids()
                .iter()
                .map(|id| id.to_ascii_uppercase().contains("GEX"))
                .collect()
        })
        .unwrap_or_default();
    let glevel = cli.params.output_gzip_level();
    let split = cli.params.split.max(1);

    let mk = |prefix: &str, s_num: usize, out: &mut Vec<PairShardWriter>| -> Result<()> {
        for c in 0..split {
            out.push(PairShardWriter::create(
                &cli.out_dir,
                prefix,
                s_num,
                c + 1,
                glevel,
                !cli.params.minimap2,
            )?);
        }
        Ok(())
    };
    let targets: Vec<PairShardWriter> = match &demuxer {
        Some(d) => {
            let mut v = Vec::with_capacity((n_samples + 1) * split);
            for (i, id) in d.sample_ids().iter().enumerate() {
                mk(id, i + 1, &mut v)?;
            }
            mk("Undetermined", 0, &mut v)?;
            v
        }
        None => {
            let mut v = Vec::with_capacity(split);
            mk("sbxpress", 1, &mut v)?;
            v
        }
    };
    let index_len = demuxer.as_ref().map_or(10, |d| d.index_len());

    let failed_writer = cli
        .params
        .failed_output
        .as_ref()
        .map(|path| ShardWriter::create(path, glevel))
        .transpose()?;
    let chunk_counter = AtomicUsize::new(0);
    let shared = Shared {
        orient_cfg,
        orient_mode: cli.params.orientation,
        extract_cfg,
        demuxer: demuxer.as_ref(),
        gex_by_sample: &gex_by_sample,
        whitelist: whitelist.as_ref(),
        probe_wl: probe_wl.as_ref(),
        cb_len,
        umi_len,
        cb_gate: cli.params.cb_gate,
        realign: cli.params.realign,
        cb_reframe: cli.params.cb_reframe,
        reframe_window: cli.params.reframe_window,
        lenient: cli.params.lenient,
        index_len,
        index_window: cli.params.index_window,
        n_samples,
        split,
        minimap2: cli.params.minimap2,
        processing_tags: cli.params.processing_tags,
        targets: &targets,
        failed: failed_writer.as_ref(),
        chunk_counter: &chunk_counter,
    };

    let mut inputs = inputs;
    inputs.sort_by_key(|p| std::cmp::Reverse(std::fs::metadata(p).map(|m| m.len()).unwrap_or(0)));

    let load_secs = t_start.elapsed().as_secs_f64();
    let n_workers = rayon::current_num_threads();
    let forced_readers = cli.params.reader_threads;
    let use_pipeline = forced_readers > 0 || inputs.len() < n_workers;
    info!(decoder = %crate::gzin::backend_description(), "gzip input decoder");

    let t_proc = Instant::now();
    let locals = if use_pipeline {
        let n_readers = if forced_readers > 0 {
            forced_readers.min(inputs.len())
        } else {
            inputs.len()
        }
        .max(1);
        info!(
            mode = "pipelined",
            readers = n_readers,
            workers = n_workers,
            n_files = inputs.len(),
            "starting parallel processing"
        );
        run_pipelined(&inputs, &shared, n_readers, n_workers)?
    } else {
        info!(
            mode = "per-file",
            rayon_threads = n_workers,
            n_files = inputs.len(),
            "starting parallel processing"
        );
        inputs
            .par_iter()
            .map(|input| process_file(input, &shared))
            .collect::<Result<Vec<_>>>()?
    };

    let mut stats = Stats::new(n_samples);
    for s in locals {
        stats.merge(s);
    }

    for t in &targets {
        t.finish()?;
    }
    if let Some(writer) = &failed_writer {
        writer.finish()?;
    }
    let process_secs = t_proc.elapsed().as_secs_f64();

    let elapsed = t_start.elapsed().as_secs_f64();
    let threads = rayon::current_num_threads();
    let timing = Timing {
        elapsed,
        load_secs,
        process_secs,
        threads,
    };
    report(&stats, demuxer.as_ref(), &timing);
    let metrics = build_metrics(
        &stats,
        whitelist.as_ref(),
        demuxer.as_ref(),
        probe_wl.as_ref(),
        &timing,
        cli.params.minimap2,
        cli.params.processing_tags,
    );
    let metrics_path = cli.out_dir.join("metrics_summary.json");
    metrics.write_json(&metrics_path)?;
    info!(out = %metrics_path.display(), "wrote metrics summary");
    Ok(())
}

fn build_metrics(
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

fn report(stats: &Stats, demuxer: Option<&Demuxer>, t: &Timing) {
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
