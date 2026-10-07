use std::collections::BTreeMap;
use std::path::PathBuf;
use clap::Parser;
use serde::Serialize;

/// Reconstruct 10x reference's raw GEX matrix from a 10x reference BAM.
#[derive(Debug, Parser)]
#[command(
    name = "cr-bam-to-mtx",
    version,
    about = "Reconstruct a 10x reference raw gene-expression matrix from its barcoded BAM"
)]
pub struct Cli {
    /// 10x reference possorted_genome_bam.bam input.
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,

    /// 10x-style outs directory.
    #[arg(short, long, default_value = "outs")]
    pub output: PathBuf,

    /// Genome label stored for every feature in raw_feature_bc_matrix.h5.
    #[arg(long, default_value = "unknown")]
    pub genome: String,

    /// Optional 10x reference features.tsv[.gz] template, including zero-count rows.
    #[arg(long)]
    pub features: Option<PathBuf>,

    /// Optional 10x reference barcodes.tsv[.gz] template, including absent columns.
    #[arg(long, conflicts_with = "cell_barcodes")]
    pub barcodes: Option<PathBuf>,

    /// 10x reference barcode whitelist (plain or .gz). Preserves whitelist order,
    /// adds the gem-group suffix, and emits zero-count columns for absent cells.
    #[arg(long, value_name = "WHITELIST", conflicts_with = "barcodes")]
    pub cell_barcodes: Option<PathBuf>,

    /// Gem-group suffix added to entries read via --cell-barcodes.
    #[arg(long, default_value_t = 1)]
    pub gem_group: u32,

    /// With --barcodes/--cell-barcodes, silently skip a counted molecule
    /// whose barcode isn't in the template instead of erroring.
    #[arg(long, default_value_t = false)]
    pub skip_unlisted_barcodes: bool,

    /// Write plain .tsv/.mtx instead of gzip-compressed files.
    #[arg(long)]
    pub no_gzip: bool,

    /// Optional per-sample probe-barcode allow-list.
    #[arg(long, value_name = "PROBE_BARCODES")]
    pub probe_barcodes: Option<PathBuf>,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub(crate) struct Feature {
    pub(crate) id: String,
    pub(crate) name: String,
}

#[derive(Debug, Default)]
pub(crate) struct Matrix {
    pub(crate) features: Vec<Feature>,
    pub(crate) barcodes: Vec<String>,
    /// (barcode index, feature index) -> molecule count. Keeping the barcode
    /// first makes Matrix Market output match 10x reference's column-major order.
    pub(crate) counts: BTreeMap<(usize, usize), u64>,
    pub(crate) records: u64,
    pub(crate) molecules: u64,
    pub(crate) barcodes_skipped: u64,
    pub(crate) probe_barcodes_skipped: u64,
}

/// Machine-readable run accounting, following the main `sbx10x` app.
#[derive(Debug, Serialize)]
pub struct Metrics {
    pub input_bam: String,
    pub output_dir: String,
    pub mex_dir: String,
    pub h5_file: String,
    pub elapsed_secs: f64,
    pub records_read: u64,
    pub records_per_sec: u64,
    pub molecules_counted: u64,
    pub features: usize,
    pub barcodes: usize,
    pub nonzero_entries: usize,
    pub used_feature_template: bool,
    pub used_barcode_template: bool,
    pub used_cell_barcode_whitelist: bool,
    pub gem_group: Option<u32>,
    pub skip_unlisted_barcodes: bool,
    pub barcodes_skipped: u64,
    pub used_probe_barcode_allowlist: bool,
    pub probe_barcodes_skipped: u64,
}
