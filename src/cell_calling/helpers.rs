use std::collections::{BTreeSet, HashMap};
use std::fs;
use std::path::Path;

use anyhow::{bail, Context, Result};
use rand::Rng;
use rand_chacha::ChaCha8Rng;

use super::types::{MIN_RECOVERED_CELLS, ORDMAG_BOOTSTRAPS, ORDMAG_QUANTILE};

/// Resolve a --tagseq-barcodes file into raw-matrix column indices.
pub fn resolve_tagseq_barcodes(path: &Path, barcodes: &[String]) -> Result<Vec<usize>> {
    let by_prefix: HashMap<&str, usize> = barcodes
        .iter()
        .enumerate()
        .map(|(i, bc)| (bc.split('-').next().unwrap_or(bc.as_str()), i))
        .collect();
    let text = fs::read_to_string(path)
        .with_context(|| format!("reading --tagseq-barcodes {}", path.display()))?;
    let mut indices = Vec::new();
    for (lineno, line) in text.lines().enumerate() {
        let bc = line.trim();
        if bc.is_empty() {
            continue;
        }
        let bare = bc.split('-').next().unwrap_or(bc);
        let &index = by_prefix.get(bare).with_context(|| {
            format!(
                "{}:{}: barcode {bare} is not present in the raw matrix",
                path.display(),
                lineno + 1
            )
        })?;
        indices.push(index);
    }
    if indices.is_empty() {
        bail!("--tagseq-barcodes {} contained no barcodes", path.display());
    }
    indices.sort_unstable();
    indices.dedup();
    Ok(indices)
}

pub fn within_ordmag(sorted_ascending: &[u64], baseline_idx: usize) -> usize {
    let idx = baseline_idx.min(sorted_ascending.len() - 1);
    let baseline = sorted_ascending[sorted_ascending.len() - idx - 1];
    let cutoff = ((baseline as f64 * 0.1).round() as u64).max(1);
    let first = sorted_ascending.partition_point(|&x| x < cutoff);
    sorted_ascending.len() - first
}

pub fn bootstrap_counts(source: &[u64], rng: &mut ChaCha8Rng) -> Vec<u64> {
    let mut sample = Vec::with_capacity(source.len());
    for _ in 0..source.len() {
        sample.push(source[rng.gen_range(0..source.len())]);
    }
    sample.sort_unstable();
    sample
}

pub fn recovered_grid(maximum: usize) -> Vec<usize> {
    let max_log = (maximum.max(2) as f64).log2();
    let mut values = BTreeSet::new();
    for i in 0..2000 {
        let exponent = 1.0 + (max_log - 1.0) * i as f64 / 1999.0;
        values.insert(2f64.powf(exponent).round() as usize);
    }
    values.into_iter().collect()
}

pub fn estimate_recovered(nonzero: &[u64], max_expected: usize, rng: &mut ChaCha8Rng) -> usize {
    let grid = recovered_grid(max_expected);
    let mut estimates = Vec::with_capacity(ORDMAG_BOOTSTRAPS);
    for _ in 0..ORDMAG_BOOTSTRAPS {
        let sample = bootstrap_counts(nonzero, rng);
        let mut best = (f64::INFINITY, MIN_RECOVERED_CELLS);
        for &expected in &grid {
            let baseline = ((expected as f64 * (1.0 - ORDMAG_QUANTILE)).round() as usize)
                .min(sample.len() - 1);
            let observed = within_ordmag(&sample, baseline);
            let loss = (observed as f64 - expected as f64).powi(2) / expected as f64;
            if loss < best.0 {
                best = (loss, expected);
            }
        }
        estimates.push(best.1);
    }
    ((estimates.iter().sum::<usize>() as f64 / estimates.len() as f64).round() as usize)
        .max(MIN_RECOVERED_CELLS)
}
