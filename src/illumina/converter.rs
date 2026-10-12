use std::fs::File;
use std::io::BufWriter;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;

use crate::barcode::Whitelist;
use crate::extract::ExtractCfg;
use crate::gzin::backend_description;
use crate::overlap::OverlapCfg;

use super::pair::convert_pair_full;
use super::types::{failed_paths, metrics_path, shard_allocation, shard_path, Cli, Metrics};

pub fn run(cli: Cli) -> Result<()> {
    if cli.r1.len() != cli.r2.len() {
        bail!("--r1 and --r2 must be supplied the same number of times");
    }
    if cli.gem_group == 0 {
        bail!("--gem-group must be greater than zero");
    }
    if cli.split == 0 {
        bail!("--split must be greater than zero");
    }
    let started = Instant::now();
    let overlap_cfg = cli.overlap.then_some(OverlapCfg {
        kmer: cli.overlap_kmer,
        kmer_offset: cli.overlap_kmer_offset,
        check_lengths: cli.overlap_check_lengths,
        min_overlap: cli.overlap_min_overlap,
        max_mismatch: cli.overlap_max_mismatch,
        max_mismatch_density: cli.overlap_max_mismatch_density,
        qual_cap: 40,
    });
    let extract_cfg = ExtractCfg {
        min_cdna_len: cli.min_cdna_len,
        adapter_max_mm: cli.adapter_max_mm,
        trim_tso: cli.trim_tso,
        max_r2_len: 0,
        // No Flex support in this binary (already-demultiplexed Illumina
        // R1/R2 input, no probe-barcode concept here). flex_r2_boundary is
        // false, so this floor is never consulted; 65 just matches Flex's
        // own default for consistency.
        flex_r2_boundary: false,
        flex_min_cdna_len_no_anchor: 65,
    };
    let emit_yc = !cli.no_emit_yc;
    let cli_trim_tso = cli.trim_tso;
    let r2_polya_min_run = cli.trim_r2_polya.then_some(cli.r2_polya_min_run);
    let whitelist = Whitelist::load(&cli.whitelist)?;
    if whitelist.is_empty() {
        bail!("whitelist is empty: {}", cli.whitelist.display());
    }
    let targets_per_input = 1 + usize::from(cli.failed_output.is_some());
    // Preserve the historical default: every target gets one writer even when
    // the caller leaves --split at 1. Explicit larger values are a total budget.
    let allocation = shard_allocation(
        cli.split.max(cli.r1.len() * targets_per_input),
        cli.r1.len() * targets_per_input,
    )?;
    let main_counts = (0..cli.r1.len())
        .map(|i| allocation[i * targets_per_input])
        .collect::<Vec<_>>();
    let total_main: usize = main_counts.iter().sum();
    let failed_counts = (0..cli.r1.len())
        .map(|i| {
            if targets_per_input == 2 {
                allocation[i * 2 + 1]
            } else {
                0
            }
        })
        .collect::<Vec<_>>();
    let total_failed: usize = failed_counts.iter().sum();
    let mut main_offset = 0usize;
    let per_input_outputs = (0..cli.r1.len())
        .map(|input| {
            let n = main_counts[input];
            let start = main_offset;
            main_offset += n;
            (0..n)
                .map(|i| shard_path(&cli.output, start + i, total_main))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let per_input_failed = cli
        .failed_output
        .as_ref()
        .map(|prefix| {
            let mut failed_offset = 0usize;
            (0..cli.r1.len())
                .map(|input| {
                    let n = failed_counts[input];
                    let start = failed_offset;
                    failed_offset += n;
                    (0..n)
                        .map(|i| failed_paths(prefix, start + i, total_failed))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    let output_shards = per_input_outputs
        .iter()
        .flatten()
        .cloned()
        .collect::<Vec<_>>();
    let failed_shards = per_input_failed
        .iter()
        .flatten()
        .cloned()
        .collect::<Vec<_>>();
    let mut metrics = Metrics {
        r1_files: cli.r1.clone(),
        r2_files: cli.r2.clone(),
        whitelist: cli.whitelist.clone(),
        whitelist_size: whitelist.len(),
        output: cli.output.clone(),
        output_shards: output_shards.clone(),
        split: cli.split,
        gem_group: cli.gem_group,
        overlap_enabled: cli.overlap,
        low_complexity_base_fraction: cli.complexity.low_complexity_base_fraction,
        low_complexity_run_fraction: cli.complexity.low_complexity_run_fraction,
        drop_low_complexity: cli.complexity.drop_low_complexity,
        failed_outputs: failed_shards
            .iter()
            .flat_map(|(r1, r2)| [r1.clone(), r2.clone()])
            .collect(),
        ..Metrics::default()
    };
    let mut per_barcode = crate::barcode_metrics::BarcodeTable::<12>::new(cli.barcode_metrics);
    let workers = if cli.threads == 0 {
        cli.r1.len()
    } else {
        cli.threads.min(cli.r1.len())
    }
    .max(1);
    eprintln!(
        "[sbxpress illumina] converting {} pair(s) with {} worker(s)",
        cli.r1.len(),
        workers
    );
    eprintln!(
        "[sbxpress illumina] gzip input decoder: {}",
        backend_description()
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .context("cannot create Illumina conversion thread pool")?;
    let pair_metrics = pool.install(|| {
        cli.r1
            .par_iter()
            .zip(cli.r2.par_iter())
            .enumerate()
            .map(|(index, (r1, r2))| {
                let result = convert_pair_full(
                    r1,
                    r2,
                    &per_input_outputs[index],
                    &whitelist,
                    cli.gem_group,
                    overlap_cfg.as_ref(),
                    &extract_cfg,
                    emit_yc,
                    r2_polya_min_run,
                    per_input_failed
                        .get(index)
                        .map(Vec::as_slice)
                        .unwrap_or(&[]),
                    cli.barcode_metrics,
                    &cli.complexity,
                );
                if let Ok(value) = &result {
                    eprintln!(
                        "[sbxpress illumina] shard {}/{} complete: read={} emitted={} \
                         sins={} mins={} lins={} output={}",
                        index + 1,
                        cli.r1.len(),
                        value.records_read,
                        value.records_emitted,
                        value.overlap_short,
                        value.overlap_medium,
                        value.overlap_long,
                        cli.output.display()
                    );
                }
                result
            })
            .collect::<Result<Vec<_>>>()
    })?;
    for value in pair_metrics {
        per_barcode.merge(value.per_barcode);
        metrics.records_read += value.records_read;
        metrics.records_emitted += value.records_emitted;
        metrics.barcode_exact += value.barcode_exact;
        metrics.barcode_hamming1_corrected += value.barcode_hamming1_corrected;
        metrics.barcode_rejected += value.barcode_rejected;
        metrics.short_r1 += value.short_r1;
        metrics.overlap_short_sins += value.overlap_short;
        metrics.overlap_medium_mins += value.overlap_medium;
        metrics.overlap_long_lins += value.overlap_long;
        metrics.merged_without_cdna += value.merged_without_cdna;
        metrics.merged_bases += value.merged_bases;
        metrics.discordant_bases += value.discordant_bases;
        metrics.yc_n_columns_coerced += value.yc_n_columns_coerced;
        metrics.yc_unencodable += value.yc_unencodable;
        metrics.tso_trimmed += value.tso_trimmed;
        metrics.r2_polya_trimmed += value.r2_polya_trimmed;
        metrics.r2_polya_discarded += value.r2_polya_discarded;
        metrics.r2_polya_bases_trimmed += value.r2_polya_bases_trimmed;
        metrics.low_complexity += value.low_complexity;
        metrics
            .low_complexity_by_base
            .add(&value.low_complexity_by_base);
        metrics.records_failed += value.records_failed;
        for (reason, count) in value.failed_by_reason {
            *metrics.failed_by_reason.entry(reason).or_insert(0) += count;
        }
        for (length, count) in value.insert_length_histogram {
            *metrics.insert_length_histogram.entry(length).or_insert(0) += count;
        }
    }
    metrics.elapsed_secs = started.elapsed().as_secs_f64();
    let metrics_path = cli.metrics.unwrap_or_else(|| metrics_path(&cli.output));
    let metrics_file = File::create(&metrics_path)
        .with_context(|| format!("cannot create {}", metrics_path.display()))?;
    serde_json::to_writer_pretty(BufWriter::new(metrics_file), &metrics)?;
    per_barcode.write(
        &metrics_path,
        &[
            "reads",
            "barcode_exact",
            "barcode_corrected",
            "emitted",
            "failed_merged_too_short",
            "failed_no_cdna",
            "failed_short_r1",
            "failed_r2_polya_short",
            "failed_cb_invalid",
            "failed_cb_unmatched",
            "low_complexity",
            "failed_low_complexity",
        ],
    )?;
    eprintln!(
        "[sbxpress illumina] read={} emitted={} rejected={} metrics={}",
        metrics.records_read,
        metrics.records_emitted,
        metrics.barcode_rejected,
        metrics_path.display()
    );
    // Poly-A read-through is expected unless R2 is poly-A trimmed, so with
    // --trim-r2-polya false every flagged read counts toward the warning;
    // otherwise G/N only.
    crate::complexity::warn_if_high(
        "sbxpress illumina",
        &metrics.low_complexity_by_base,
        metrics.records_emitted + metrics.failed_by_reason.get("low_complexity").unwrap_or(&0),
        cli.trim_r2_polya,
        cli.complexity.drop_low_complexity,
    );
    if !metrics.failed_outputs.is_empty() {
        let reasons = metrics
            .failed_by_reason
            .iter()
            .map(|(reason, count)| format!("{reason}={count}"))
            .collect::<Vec<_>>()
            .join(" ");
        eprintln!(
            "[sbxpress illumina] failed: {} pair(s) [{}] -> {}",
            metrics.records_failed,
            reasons,
            metrics
                .failed_outputs
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        // Every input record is now accounted for in exactly one place.
        debug_assert_eq!(
            metrics.records_read,
            metrics.records_emitted + metrics.records_failed
        );
    }
    if metrics.overlap_enabled {
        let merged = metrics.overlap_short_sins + metrics.overlap_medium_mins;
        let considered = merged + metrics.overlap_long_lins;
        let pct = |n: u64| {
            if considered > 0 {
                100.0 * n as f64 / considered as f64
            } else {
                0.0
            }
        };
        eprintln!(
            "[sbxpress illumina] overlap: sins(short)={} ({:.1}%) mins(medium)={} ({:.1}%) \
             lins(long)={} ({:.1}%)",
            metrics.overlap_short_sins,
            pct(metrics.overlap_short_sins),
            metrics.overlap_medium_mins,
            pct(metrics.overlap_medium_mins),
            metrics.overlap_long_lins,
            pct(metrics.overlap_long_lins),
        );
        let insert_mean = if merged > 0 {
            metrics
                .insert_length_histogram
                .iter()
                .map(|(len, n)| *len as u64 * n)
                .sum::<u64>() as f64
                / merged as f64
        } else {
            0.0
        };
        eprintln!(
            "[sbxpress illumina] merged: fragments={} mean_cdna={:.1}bp discordance={:.4}% \
             no_cdna={}",
            merged,
            insert_mean,
            if metrics.merged_bases > 0 {
                100.0 * metrics.discordant_bases as f64 / metrics.merged_bases as f64
            } else {
                0.0
            },
            metrics.merged_without_cdna,
        );
        if cli_trim_tso {
            eprintln!(
                "[sbxpress illumina] merged: tso_trimmed={} ({:.1}% of merged)",
                metrics.tso_trimmed,
                if merged > 0 {
                    100.0 * metrics.tso_trimmed as f64 / merged as f64
                } else {
                    0.0
                }
            );
        }
    }
    if cli.trim_r2_polya {
        eprintln!(
            "[sbxpress illumina] r2 poly-A trim: trimmed={} discarded={} bases_trimmed={}",
            metrics.r2_polya_trimmed, metrics.r2_polya_discarded, metrics.r2_polya_bases_trimmed,
        );
    }
    Ok(())
}
