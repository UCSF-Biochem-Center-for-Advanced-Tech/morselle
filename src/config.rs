//! CLI definition, chemistry constants, and runtime parameters.

use std::path::PathBuf;

use clap::Parser;

use crate::demux::DemuxCfg;
use crate::extract::ExtractCfg;
use crate::orient::{FallbackMode, OrientCfg, OrientMode};

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

    /// Orientation fallback for reads lacking a homopolymer (M1). `hamming-tail`
    /// (default) is cheap and high-recall; `levenshtein` is thorough but ~20×
    /// slower; `off` disables the fallback.
    #[arg(long, value_enum, default_value_t = FallbackMode::HammingTail)]
    pub orient_fallback: FallbackMode,

    /// Run orientation. `auto` (default) detects + enforces the dominant
    /// orientation per file; `p5-led`/`p7-led` force it (run_5p / run1).
    #[arg(long, value_enum, default_value_t = OrientMode::Auto)]
    pub orientation: OrientMode,

    /// Max mismatches allowed when matching sample-sheet indices.
    #[arg(long, default_value_t = 1)]
    pub index_max_mm: u32,

    /// Allow demux assignment from one index when the other is truncated (M3).
    /// Safe here (min pairwise Hamming i7=3/i5=4); disable for strict dual-index.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub single_index: bool,

    /// Enable the indel-realign branch (M5). Off by default.
    #[arg(long, default_value_t = false)]
    pub realign: bool,

    /// Lenient CB recovery: for confidently-located CB-fails, accept a unique
    /// Hamming-2 whitelist hit and rewrite the CB. Off by default (slower, lower
    /// precision). Only applied when the CB location is dual-anchor confirmed.
    #[arg(long, default_value_t = false)]
    pub lenient: bool,

    /// Keep partial reads (those missing one library end).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub keep_partials: bool,

    /// Minimum cDNA length to emit an R2. SBX cDNA is long (~220 bp mean) and the
    /// no-R2 reads are cDNA-*absent* rather than short, so lowering this recovers
    /// almost nothing (measured +0.02%); kept as a knob at a safe default.
    #[arg(long, default_value_t = 20)]
    pub min_cdna_len: usize,

    /// With a whitelist, drop pairs whose CB is neither on-whitelist nor within
    /// Hamming-1 (M4). Disable to emit all and let Cell Ranger filter.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub cb_gate: bool,

    /// Max mismatches for the 13 bp 3' adapter seed used to trim cDNA (M2).
    #[arg(long, default_value_t = 2)]
    pub adapter_max_mm: u32,
}

impl Cli {
    /// R1 length = CB + UMI (28 bp for v3 GEX).
    pub fn r1_len(&self) -> usize {
        self.params.cb_len + self.params.umi_len
    }
}

impl Params {
    /// Build the M1 orientation config from these parameters.
    pub fn orient_cfg(&self) -> OrientCfg {
        OrientCfg {
            cb_len: self.cb_len,
            umi_len: self.umi_len,
            min_homopolymer: self.min_homopolymer,
            r1_tail_max_mm: self.r1_tail_max_mm,
            anchor_max_edits: self.anchor_max_edits,
            fallback: self.orient_fallback,
            expected: None,
        }
    }

    /// Build the M2 extraction config from these parameters.
    pub fn extract_cfg(&self) -> ExtractCfg {
        ExtractCfg {
            min_cdna_len: self.min_cdna_len,
            adapter_max_mm: self.adapter_max_mm,
        }
    }

    /// Build the M3 demux config from these parameters.
    pub fn demux_cfg(&self) -> DemuxCfg {
        DemuxCfg {
            index_max_mm: self.index_max_mm,
            allow_single_index: self.single_index,
        }
    }
}
