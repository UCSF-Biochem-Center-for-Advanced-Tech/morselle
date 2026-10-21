use crate::cr_bam_to_mtx::Feature;
use clap::Parser;
use serde::Serialize;
use std::path::PathBuf;

pub const ORDMAG_BOOTSTRAPS: usize = 100;
pub const ORDMAG_QUANTILE: f64 = 0.99;
pub const MIN_RECOVERED_CELLS: usize = 50;

#[derive(Debug, Parser)]
#[command(about = "Call cells with 10x OrdMag + EmptyDrops and write filtered matrices")]
pub struct Cli {
    /// raw_feature_bc_matrix.h5, or an outs directory containing it.
    #[arg(value_name = "RAW_MATRIX")]
    pub raw_matrix: PathBuf,

    /// 10x-style outs directory receiving filtered outputs.
    #[arg(short, long, default_value = "outs")]
    pub output: PathBuf,

    /// Expected recovered cells. If omitted, estimate with the 10x OrdMag grid search.
    #[arg(long, conflicts_with = "force_cells")]
    pub expect_cells: Option<usize>,

    /// Select exactly the top N nonzero barcodes and skip EmptyDrops.
    #[arg(long, conflicts_with = "expect_cells")]
    pub force_cells: Option<usize>,

    /// Skip OrdMag and EmptyDrops entirely and call exactly these barcodes
    /// as cells: one per line, matching the extraction whitelist format
    /// (a gem-group "-1" suffix, if present, is ignored).
    #[arg(long, conflicts_with_all = ["expect_cells", "force_cells"])]
    pub tagseq_barcodes: Option<PathBuf>,

    /// First barcode rank used for the ambient background (NextGEM default).
    #[arg(long, default_value_t = 45_000)]
    pub ambient_lower: usize,

    /// Exclusive last barcode rank used for the ambient background.
    #[arg(long, default_value_t = 90_000)]
    pub ambient_upper: usize,

    /// Minimum EmptyDrops candidate UMIs; raised above the ambient maximum.
    #[arg(long, default_value_t = 500)]
    pub emptydrops_minimum_umis: u64,

    /// Monte Carlo multinomial simulations per tested UMI total.
    #[arg(long, default_value_t = 100_000)]
    pub simulations: usize,

    /// Benjamini-Hochberg FDR for 5' / NextGEM EmptyDrops calls.
    #[arg(long, default_value_t = 0.001)]
    pub fdr: f64,

    /// Reproducible bootstrap and Monte Carlo seed.
    #[arg(long, default_value_t = 0)]
    pub seed: u64,
}

/// 10x-style CSC feature-barcode matrix (features x barcodes, `indptr`
/// column/barcode-indexed, `indices` = feature row per nonzero entry).
///
/// `data` is `f64` so that fractional matrices (`morselle molecules
/// --fractional-collisions`, written by `cr-bam-to-mtx` as an `f32`
/// dataset) are read without truncation; integer (`u32`) matrices convert
/// exactly. Algorithms that need integer UMI totals (OrdMag ranking,
/// EmptyDrops' multinomial draws) use [`round_count`] on per-barcode or
/// per-feature *sums*, never on individual fractional entries, so a 0.25
/// molecule still contributes to its barcode's total.
pub(crate) struct RawMatrix {
    pub(crate) features: Vec<Feature>,
    pub(crate) barcodes: Vec<String>,
    pub(crate) genome: String,
    pub(crate) data: Vec<f64>,
    pub(crate) indices: Vec<u32>,
    pub(crate) indptr: Vec<u32>,
}

#[derive(Serialize)]
pub struct Metrics {
    pub algorithm: &'static str,
    pub raw_matrix: String,
    pub filtered_h5: String,
    pub raw_barcodes: usize,
    pub nonzero_barcodes: usize,
    pub recovered_cells: usize,
    pub ordmag_cells: usize,
    pub ordmag_cutoff_umis: u64,
    pub ambient_rank_lower: usize,
    pub ambient_rank_upper: usize,
    pub ambient_barcodes_nonzero: usize,
    pub max_ambient_umis: u64,
    pub emptydrops_minimum_umis: u64,
    pub emptydrops_candidates: usize,
    pub emptydrops_cells: usize,
    pub emptydrops_simulations: usize,
    pub emptydrops_fdr: f64,
    pub filtered_cells: usize,
    pub filtered_molecules: u64,
    pub elapsed_secs: f64,
}

/// Round a (possibly fractional) summed count to the nearest whole UMI.
pub(crate) fn round_count(value: f64) -> u64 {
    value.round().max(0.0) as u64
}
