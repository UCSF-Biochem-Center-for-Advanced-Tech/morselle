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
    #[arg(long, short = 'o', default_value = "sbxpress_out")]
    pub out_dir: PathBuf,

    /// BCL Convert v2 sample sheet (or flat `sample,i7,i5` CSV) for i5/i7
    /// index-based per-library demultiplexing. Scope is deliberately narrow:
    /// this only applies to reads that passed through an Illumina library
    /// conversion step and therefore carry i5/i7 indices at all. A direct
    /// Roche library (no Illumina BCL Convert step) has no i5/i7 present --
    /// its CB/UMI live at a different position and it must not be given this
    /// flag. Demux is entirely opt-in: with this unset, no index-based
    /// splitting is attempted regardless of what the raw reads contain.
    #[arg(long)]
    pub bclconvert_samplesheet: Option<PathBuf>,

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

    /// Max mismatches allowed when matching BCL Convert sample-sheet indices.
    #[arg(long, default_value_t = 1)]
    pub index_max_mm: u32,

    /// Allow demux assignment from one index when the other is truncated (M3).
    /// Safe here (min pairwise Hamming i7=3/i5=4); disable for strict dual-index.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub single_index: bool,

    /// Max frame shift (each side) tried when matching i5/i7, to absorb indels in
    /// the R1-primer/adapter region that shift the index. 0 = rigid frame.
    #[arg(long, default_value_t = 3)]
    pub index_window: isize,

    /// Single-indel CB realignment (M5): when the measured CB+UMI span is off by
    /// exactly ±1, the CB carries a homopolymer indel that frame-shifted the 16 bp
    /// grab — recover it by trying the single-indel corrections for a unique
    /// whitelist hit. Gated on that structural span signal, so it is precision-safe
    /// (validated: recoveries land on real-cell barcodes ~4–10× above chance, and
    /// genuine single-indels are correct-or-abstain in simulation). Recovers ~+1% of
    /// reads. On by default; disable with `--realign false`.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub realign: bool,

    /// CB reframe rescue: for confidently-located CB-fails, try nearby frames
    /// (a homopolymer indel shifts the 16 bp grab) and accept a unique within-1
    /// whitelist hit, re-registering the UMI off the dT anchor. On by default —
    /// confident-location + unique-hit guards keep precision.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub cb_reframe: bool,

    /// Max CB frame shift (each side) tried by the reframe rescue. Offset +1 (the
    /// dominant SBX homopolymer-deletion mode) is tried first.
    #[arg(long, default_value_t = 2)]
    pub reframe_window: usize,

    /// Lenient CB recovery: for confidently-located CB-fails, accept a unique
    /// Hamming-2 whitelist hit and rewrite the CB. Off by default (slower, lower
    /// precision). Only applied when the CB location is dual-anchor confirmed.
    #[arg(long, default_value_t = false)]
    pub lenient: bool,

    /// Keep partial reads (those missing one library end).
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub keep_partials: bool,

    /// Minimum cDNA length to emit an R2 — the shortest reliably mappable insert.
    /// SBX cDNA is long (~220 bp mean); reads below this are treated as no-cDNA.
    #[arg(long, default_value_t = 50)]
    pub min_cdna_len: usize,

    /// Cap emitted R2 (cDNA) length; 0 = uncapped (full-length SBX cDNA).
    /// Compatibility workflows expect standard short R2 (~90 bp) and can
    /// mis-trim very long reads (seq/qual desync), so set this to your chemistry's
    /// R2 length (e.g. 90). STAR uses ~90 bp for gene mapping regardless.
    #[arg(long, default_value_t = 0)]
    pub max_r2_len: usize,

    /// With a whitelist, drop pairs whose CB is neither on-whitelist nor within
    /// Hamming-1 (M4). Disable to emit all and let 10x reference filter.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub cb_gate: bool,

    /// Max mismatches for the 13 bp 3' adapter seed used to trim cDNA (M2).
    #[arg(long, default_value_t = 2)]
    pub adapter_max_mm: u32,

    /// Emit every read that fails to produce a pair to this gzipped FASTQ, as
    /// the original (pre-orientation) read with an `XF:Z:<reason>` tag appended.
    /// With this set, every input read is accounted for in exactly one place.
    #[arg(long)]
    pub failed_output: Option<PathBuf>,

    /// Trim the template-switch oligo when it sits between the cDNA and the 3'
    /// read-through adapter. A fragment carries the TSO only when it reached the
    /// cDNA's 5' end, so it is enriched in short fragments. Off by default:
    /// enabling it changes emitted cDNA, so measure the delta before adopting.
    #[arg(long, default_value_t = false)]
    pub trim_tso: bool,

    /// gzip compression level (0–9) for output FASTQs. Each file is written as a
    /// single gzip member (10x reference mis-frames records at multi-member gzip
    /// boundaries). Lower is faster but larger; the default (4) favors throughput,
    /// since output compression is the run's wall-clock floor at scale.
    #[arg(long, default_value_t = 4)]
    pub gzip_level: u32,

    /// Write plain, uncompressed FASTQ instead of gzip (overrides
    /// --gzip-level). Measured to remove ~54% of extract's CPU time (see
    /// PROFILING.md) — gzip encoding is pure overhead when the output is
    /// about to be read once by a co-located process (e.g. piped straight
    /// into a mapper via a named pipe) rather than stored.
    #[arg(long, default_value_t = false)]
    pub no_gzip: bool,

    /// Split each sample's output across N gzip files (10x reference reads them all:
    /// `_001`..`_00N`). N>1 lets N writer threads compress one sample in parallel,
    /// removing the single-encoder bottleneck on hot samples. Keep N ≤ --threads
    /// and mind the open-file limit (`ulimit -n`): N × (samples+1) × 2 files.
    #[arg(long, default_value_t = 1)]
    pub split: usize,

    /// I/O parallelism strategy. `0` (auto): when there are fewer input files than
    /// worker threads (e.g. a single large file), decouple gzip decode from
    /// processing — one reader thread per file streams record batches to the full
    /// worker pool, so a 1-file run still uses every core (serial gzip inflate is
    /// then the only ceiling). With many files the classic one-file-per-thread path
    /// is used instead. Set `N>0` to force the pipelined path with N reader threads
    /// (capped at the file count); useful to tune decode vs. process balance.
    #[arg(long, default_value_t = 0)]
    pub reader_threads: usize,

    /// Compatibility preset for 10x-style paired FASTQ output. Applies the needed
    /// settings (currently `--max-r2-len 649`: compatible consumers have a 649 bp read cap and
    /// fails at 650 bp with "quality string length is not equal to sequence
    /// length"). Explicit flags override the preset. More presets may be added.
    #[arg(long, default_value_t = false)]
    pub cellranger: bool,

    /// Minimap2-ready output preset: emit R2 only and put CB/UMI SAM tags in the
    /// FASTQ comment. Align with minimap2 `-y` to preserve those tags in SAM/BAM.
    #[arg(long, default_value_t = false)]
    pub minimap2: bool,

    /// Add SBX processing-provenance SAM tags to each R2 FASTQ comment. These
    /// tags are copied into SAM/BAM by minimap2 `-y`. Enabled automatically by
    /// `--minimap2`.
    #[arg(long, default_value_t = false)]
    pub processing_tags: bool,

    /// 10x Flex probe-barcode whitelist (e.g. `flex-v2-384.txt`). When set,
    /// locates + identifies the probe barcode on each read and adds
    /// `PR:Z:<raw>`/`PB:Z:<corrected id>` to the FASTQ comment (requires
    /// `--minimap2`, same as CR/CB). See FLEX_SUPPORT_PLAN.md Phase 1.
    #[arg(long, requires = "minimap2")]
    pub probe_barcode_whitelist: Option<PathBuf>,

    /// Max mismatches for a confident probe-barcode whitelist match.
    #[arg(long, default_value_t = 1)]
    pub probe_barcode_max_mm: u32,

    /// Max bp each side of the empirically-confirmed offset (see
    /// FLEX_SUPPORT_PLAN.md Phase 0) searched for the probe barcode, to absorb
    /// per-read jitter in the measured dT-run boundary.
    #[arg(long, default_value_t = 10)]
    pub probe_barcode_window: isize,
}

/// R2 length cap accepted by the compatibility preset.
pub const COMPAT_MAX_R2_LEN: usize = 649;

impl Cli {
    /// R1 length = CB + UMI (28 bp for v3 GEX).
    pub fn r1_len(&self) -> usize {
        self.params.cb_len + self.params.umi_len
    }
}

impl Params {
    /// Apply convenience presets (e.g. `--cellranger`) to the parsed parameters.
    /// A preset only fills in values the user left at their default, so an
    /// explicit flag always wins. Call once, right after parsing.
    pub fn apply_presets(&mut self) {
        if self.minimap2 {
            self.processing_tags = true;
        }
        if self.cellranger {
            // The compatibility preset caps reads at 649 bp. Only set it
            // if the user hasn't already chosen a cap.
            if self.max_r2_len == 0 {
                self.max_r2_len = COMPAT_MAX_R2_LEN;
            }
        }
    }

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
            trim_tso: self.trim_tso,
            max_r2_len: self.max_r2_len,
        }
    }

    /// Build the M3 demux config from these parameters.
    pub fn demux_cfg(&self) -> DemuxCfg {
        DemuxCfg {
            index_max_mm: self.index_max_mm,
            allow_single_index: self.single_index,
        }
    }

    /// Effective output compression level: `None` (plain, uncompressed
    /// `.fastq`) when `--no-gzip` is set, else `Some(--gzip-level)`.
    pub fn output_gzip_level(&self) -> Option<u32> {
        (!self.no_gzip).then_some(self.gzip_level)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Params {
        // The Params group carries no required fields, so it parses standalone.
        let mut p = Params::parse_from(std::iter::once("sbx10x").chain(args.iter().copied()));
        p.apply_presets();
        p
    }

    #[test]
    fn no_preset_leaves_r2_uncapped() {
        assert_eq!(parse(&[]).max_r2_len, 0);
    }

    #[test]
    fn compat_preset_sets_r2_cap() {
        assert_eq!(parse(&["--cellranger"]).max_r2_len, COMPAT_MAX_R2_LEN);
    }

    #[test]
    fn explicit_max_r2_len_overrides_preset() {
        assert_eq!(
            parse(&["--cellranger", "--max-r2-len", "90"]).max_r2_len,
            90
        );
    }

    #[test]
    fn minimap2_preset_leaves_r2_uncapped() {
        let p = parse(&["--minimap2"]);
        assert_eq!(p.max_r2_len, 0);
        assert!(p.processing_tags);
    }

    #[test]
    fn processing_tags_are_opt_in_without_minimap2() {
        assert!(!parse(&[]).processing_tags);
        assert!(parse(&["--processing-tags"]).processing_tags);
    }
}
