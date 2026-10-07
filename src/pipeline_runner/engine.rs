use std::sync::atomic::AtomicUsize;
use std::time::Instant;

use anyhow::Result;
use rayon::prelude::*;
use tracing::info;

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

    // Load the cell-barcode whitelist (M4 gate), if provided.
    let whitelist = match &cli.whitelist {
        Some(path) => {
            let w = Whitelist::load(path)?;
            info!(barcodes = w.len(), path = %path.display(), "loaded CB whitelist");
            Some(w)
        }
        None => None,
    };

    // Load the Flex probe-barcode whitelist (docs/archive/FLEX_SUPPORT_PLAN.md Phase 1), if
    // provided. `--minimap2` is required (clap `requires`) since PR/PB ride the
    // same FASTQ-comment -> SAM-tag mechanism as CR/CB.
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

    // Output targets: per-sample then Undetermined (last) when demuxing, else one
    // combined target. Each file is a single thread-safe gzip member.
    let demuxer = sheet
        .as_ref()
        .map(|s| Demuxer::new(s, cli.params.demux_cfg()))
        .transpose()?;
    let n_samples = demuxer.as_ref().map_or(1, |d| d.sample_ids().len());
    // Naming-convention classification (matches this project's real sample
    // sheets throughout, e.g. TSeq100_L1_GEX / TSeq100_L1_CRI): used only to
    // route ExtractCfg::flex_r2_boundary's hard-gate default away from
    // libraries the sheet itself says aren't Flex GEX. A sample the sheet
    // doesn't clearly mark this way is treated as GEX by default (the
    // `--flex` "assume GEX unless told otherwise" rule) -- never a
    // correctness dependency, just an efficiency/precision routing hint.
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

    // Each logical target gets `split` chunk files (_001.._00N). Layout is
    // [target0 chunk0..N-1][target1 chunk0..N-1]...; physical = t*split + c.
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

    // Process files in parallel — each plain-gzip file inflates serially on one
    // worker, so parallelism comes from running many SBX files at once. Sort
    // largest-first (LPT) for better load balance across the rayon pool.
    let mut inputs = inputs;
    inputs.sort_by_key(|p| std::cmp::Reverse(std::fs::metadata(p).map(|m| m.len()).unwrap_or(0)));

    // Everything up to here (sample sheet, file discovery, whitelist load, target
    // creation) is one-time setup; time it separately from per-read processing.
    let load_secs = t_start.elapsed().as_secs_f64();
    // Path selection. With many files (files ≥ workers) the one-file-per-thread
    // path is optimal: each file inflates+processes on its own core, giving up to
    // `workers` concurrent decodes. With few files (e.g. a single large one) that
    // path leaves cores idle — serial gzip inflate caps one file to one core — so
    // the pipelined path dedicates a reader per file and hands record batches to the
    // full worker pool. `--reader-threads N` forces the pipelined path with N readers.
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

use super::report::{build_metrics, report};
