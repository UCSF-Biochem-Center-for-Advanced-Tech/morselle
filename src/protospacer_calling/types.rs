use std::path::{Path, PathBuf};
use clap::Parser;
use serde::Serialize;

use crate::cell_calling::RawMatrix;

/// Call CRISPR guide identity per cell from a filtered feature-barcode
/// matrix, fitting a per-guide adaptive UMI threshold.
#[derive(Debug, Parser)]
#[command(
    name = "sbxpress-protospacer-call",
    version,
    about = "Call CRISPR guide identity per cell (Gaussian-mixture adaptive UMI threshold)"
)]
pub struct Cli {
    /// Per-sample filtered CRISPR feature-barcode matrix.
    #[arg(value_name = "FILTERED_MATRIX")]
    pub filtered_matrix: PathBuf,

    /// Output directory.
    #[arg(short, long, default_value = "outs")]
    pub output: PathBuf,

    /// Unconditional UMI floor.
    #[arg(long, default_value_t = 3)]
    pub min_umi: u64,
}

pub fn resolve_filtered_h5(path: &Path) -> PathBuf {
    if path.is_dir() {
        path.join("filtered_feature_bc_matrix.h5")
    } else {
        path.to_owned()
    }
}

/// Per-barcode UMI count for one feature row, in the matrix's barcode order.
pub(crate) fn feature_counts(raw: &RawMatrix, feature_row: usize) -> Vec<u64> {
    (0..raw.barcodes.len())
        .map(|bc| {
            let (start, end) = (raw.indptr[bc] as usize, raw.indptr[bc + 1] as usize);
            (start..end)
                .find(|&i| raw.indices[i] as usize == feature_row)
                .map_or(0, |i| crate::cell_calling::round_count(raw.data[i]))
        })
        .collect()
}

#[derive(Debug, Clone, Copy)]
pub struct Component {
    pub weight: f64,
    pub mean: f64,
    pub std: f64,
}

/// Floor on a fitted component's std-dev.
pub const MIN_STD: f64 = 1e-3;

/// Initial (low-percentile, high-percentile) pairs tried as EM starting points.
pub const INIT_PERCENTILES: &[(f64, f64)] = &[
    (0.1, 0.9),
    (0.25, 0.75),
    (0.05, 0.95),
    (0.4, 0.6),
    (0.5, 0.95),
];

/// Number of steps used to numerically scan for the density crossing.
pub const BOUNDARY_SCAN_STEPS: usize = 20_000;

#[derive(Debug, Serialize)]
pub struct SummaryRow {
    pub feature_call: String,
    pub num_cells: usize,
    pub pct_cells: f64,
    pub median_umis: Option<f64>,
    pub stddev_umis: Option<f64>,
}

pub fn median(values: &[u64]) -> f64 {
    let mut v: Vec<u64> = values.to_vec();
    v.sort_unstable();
    let mid = v.len() / 2;
    if v.len().is_multiple_of(2) {
        (v[mid - 1] + v[mid]) as f64 / 2.0
    } else {
        v[mid] as f64
    }
}

pub fn stddev(values: &[u64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().map(|&v| v as f64).sum::<f64>() / n;
    (values
        .iter()
        .map(|&v| (v as f64 - mean).powi(2))
        .sum::<f64>()
        / n)
        .sqrt()
}
