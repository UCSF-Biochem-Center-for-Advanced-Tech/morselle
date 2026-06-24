//! CLI definition, chemistry constants, and runtime parameters.

use std::path::PathBuf;

use clap::Parser;

/// Convert raw Roche SBX reads into per-sample 10x R1/R2 FASTQs.
#[derive(Parser, Debug, Clone)]
#[command(name = "sbx10x", version, about)]
pub struct Cli {
    /// Raw SBX FASTQ(s): a directory of `PrimaryAnalysis_*.fastq.gz` or explicit files.
    #[arg(long, num_args = 1.., required = true)]
    pub fastqs: Vec<PathBuf>,

    /// Output directory for per-sample FASTQs and metrics.
    #[arg(long, short = 'o', default_value = "sbx10x_out")]
    pub out_dir: PathBuf,

    /// Sample sheet CSV with columns: sample,i7,i5. (Used from M3 onward.)
    #[arg(long)]
    pub sample_sheet: Option<PathBuf>,

    /// 10x cell-barcode whitelist (e.g. 3M-february-2018.txt[.gz]). (Used from M4.)
    #[arg(long)]
    pub whitelist: Option<PathBuf>,

    /// Worker threads (0 = rayon default = num cpus).
    #[arg(long, default_value_t = 0)]
    pub threads: usize,

    #[command(flatten)]
    pub params: Params,
}

/// Chemistry / tuning parameters with sane defaults (plan §5).
#[derive(Parser, Debug, Clone)]
pub struct Params {
    /// Cell-barcode length (bp).
    #[arg(long, default_value_t = 16)]
    pub cb_len: usize,

    /// UMI length (bp).
    #[arg(long, default_value_t = 12)]
    pub umi_len: usize,

    /// Minimum homopolymer run length to accept the dT-primer anchor (M1).
    #[arg(long, default_value_t = 12)]
    pub min_homopolymer: usize,

    /// Max mismatches for the cheap 10 bp R1-primer-tail bracket confirm (M1).
    #[arg(long, default_value_t = 2)]
    pub r1_tail_max_mm: u32,

    /// Max edits allowed when fuzzy-locating adapters/anchors (fallback path).
    #[arg(long, default_value_t = 6)]
    pub anchor_max_edits: u32,

    /// Max mismatches allowed when matching sample-sheet indices.
    #[arg(long, default_value_t = 1)]
    pub index_max_mm: u32,

    /// Enable the indel-realign branch (M5). Off by default.
    #[arg(long, default_value_t = false)]
    pub realign: bool,

    /// Keep partial reads (those missing one library end).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub keep_partials: bool,

    /// Minimum cDNA length to emit an R2.
    #[arg(long, default_value_t = 20)]
    pub min_cdna_len: usize,
}

impl Cli {
    /// R1 length = CB + UMI (28 bp for v3 GEX).
    pub fn r1_len(&self) -> usize {
        self.params.cb_len + self.params.umi_len
    }
}
