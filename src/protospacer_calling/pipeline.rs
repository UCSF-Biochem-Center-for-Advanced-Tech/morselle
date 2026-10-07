use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write as _};
use std::time::Instant;

use anyhow::{ensure, Context, Result};
use tracing::info;

use super::gmm::fit_weighted_threshold;
use super::types::{
    feature_counts, resolve_filtered_h5, weighted_median, weighted_stddev, Cli, SummaryRow,
};
use crate::cell_calling::read_raw;

pub fn run(cli: &Cli) -> Result<()> {
    let started = Instant::now();
    let h5 = resolve_filtered_h5(&cli.filtered_matrix);
    let raw = read_raw(&h5)?;
    let n_cells = raw.barcodes.len();
    ensure!(n_cells > 0, "{} has no cells", h5.display());

    let mut thresholds: HashMap<usize, u64> = HashMap::new();
    let mut counts_by_feature: HashMap<usize, Vec<f64>> = HashMap::new();
    for (fi, _feature) in raw.features.iter().enumerate() {
        let counts = feature_counts(&raw, fi);
        if let Some(t) = fit_weighted_threshold(&counts, cli.min_umi) {
            thresholds.insert(fi, t);
        }
        counts_by_feature.insert(fi, counts);
    }
    let called_features: Vec<usize> = {
        let mut v: Vec<usize> = thresholds.keys().copied().collect();
        v.sort_unstable();
        v
    };

    // Per-cell calls: which called features clear their threshold.
    let mut per_cell_calls: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_cells];
    for &fi in &called_features {
        let threshold = thresholds[&fi];
        let counts = &counts_by_feature[&fi];
        for (bc, &count) in counts.iter().enumerate() {
            if count >= threshold as f64 {
                per_cell_calls[bc].push((fi, count));
            }
        }
    }

    fs::create_dir_all(&cli.output)
        .with_context(|| format!("creating {}", cli.output.display()))?;

    // protospacer_calls_per_cell.csv
    {
        let mut w = BufWriter::new(File::create(
            cli.output.join("protospacer_calls_per_cell.csv"),
        )?);
        writeln!(w, "cell_barcode,num_features,feature_call,num_umis")?;
        for (bc, calls) in per_cell_calls.iter().enumerate() {
            if calls.is_empty() {
                continue;
            }
            let ids: Vec<&str> = calls
                .iter()
                .map(|&(fi, _)| raw.features[fi].id.as_str())
                .collect();
            let umis: Vec<String> = calls.iter().map(|&(_, c)| c.to_string()).collect();
            writeln!(
                w,
                "{},{},{},{}",
                raw.barcodes[bc],
                calls.len(),
                ids.join("|"),
                umis.join("|"),
            )?;
        }
    }

    // protospacer_umi_thresholds.csv / .json
    {
        let mut csv = BufWriter::new(File::create(
            cli.output.join("protospacer_umi_thresholds.csv"),
        )?);
        writeln!(csv, "Protospacer,UMI threshold")?;
        let mut json_map = serde_json::Map::new();
        for &fi in &called_features {
            writeln!(csv, "{},{}", raw.features[fi].id, thresholds[&fi])?;
            json_map.insert(
                raw.features[fi].id.clone(),
                serde_json::json!(thresholds[&fi]),
            );
        }
        fs::write(
            cli.output.join("protospacer_umi_thresholds.json"),
            serde_json::to_string_pretty(&json_map)?,
        )?;
    }

    // protospacer_calls_summary.csv
    {
        let mut rows = Vec::new();
        let no_guide_molecules = (0..n_cells)
            .filter(|&bc| {
                called_features
                    .iter()
                    .all(|&fi| counts_by_feature[&fi][bc] == 0.0)
            })
            .count();
        let no_confident_call = per_cell_calls.iter().filter(|c| c.is_empty()).count();
        let one_assigned: Vec<usize> = per_cell_calls
            .iter()
            .enumerate()
            .filter(|(_, c)| c.len() == 1)
            .map(|(bc, _)| bc)
            .collect();
        let more_than_one = per_cell_calls.iter().filter(|c| c.len() > 1).count();

        // Fraction (0-1 scale), matching real cellranger output exactly
        // (e.g. `0.70568911653102`, not a percentage).
        let pct = |n: usize| n as f64 / n_cells as f64;
        rows.push(SummaryRow {
            feature_call: "No guide molecules".into(),
            num_cells: no_guide_molecules,
            pct_cells: pct(no_guide_molecules),
            median_umis: None,
            stddev_umis: None,
        });
        rows.push(SummaryRow {
            feature_call: "No confident call".into(),
            num_cells: no_confident_call,
            pct_cells: pct(no_confident_call),
            median_umis: None,
            stddev_umis: None,
        });
        rows.push(SummaryRow {
            feature_call: "1 protospacer assigned".into(),
            num_cells: one_assigned.len(),
            pct_cells: pct(one_assigned.len()),
            median_umis: None,
            stddev_umis: None,
        });
        rows.push(SummaryRow {
            feature_call: "More than 1 protospacer assigned".into(),
            num_cells: more_than_one,
            pct_cells: pct(more_than_one),
            median_umis: None,
            stddev_umis: None,
        });
        for &fi in &called_features {
            let sole_umis: Vec<f64> = one_assigned
                .iter()
                .filter_map(|&bc| {
                    per_cell_calls[bc]
                        .iter()
                        .find(|&&(f, _)| f == fi)
                        .map(|&(_, c)| c)
                })
                .collect();
            rows.push(SummaryRow {
                feature_call: raw.features[fi].id.clone(),
                num_cells: sole_umis.len(),
                pct_cells: pct(sole_umis.len()),
                median_umis: (!sole_umis.is_empty()).then(|| weighted_median(&sole_umis)),
                stddev_umis: (!sole_umis.is_empty()).then(|| weighted_stddev(&sole_umis)),
            });
        }

        let mut w = BufWriter::new(File::create(
            cli.output.join("protospacer_calls_summary.csv"),
        )?);
        writeln!(
            w,
            "feature_call,num_cells,pct_cells,median_umis,stddev_umis"
        )?;
        for r in &rows {
            writeln!(
                w,
                "{},{},{},{},{}",
                r.feature_call,
                r.num_cells,
                r.pct_cells,
                r.median_umis.map_or("None".to_string(), |v| v.to_string()),
                r.stddev_umis.map_or("None".to_string(), |v| v.to_string()),
            )?;
        }
    }

    info!(
        cells = n_cells,
        called_features = called_features.len(),
        called_cells = per_cell_calls.iter().filter(|c| !c.is_empty()).count(),
        elapsed = format!("{:.1}s", started.elapsed().as_secs_f64()),
        out = %cli.output.display(),
        "RUN SUMMARY"
    );
    Ok(())
}
