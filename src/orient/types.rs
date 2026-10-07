use std::borrow::Cow;
use crate::pipeline::Orientation;

/// How orientation was determined, for accounting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrientMethod {
    /// Oriented from the dT-primer homopolymer (common case).
    Homopolymer,
    /// Oriented from a windowed Hamming scan for the R1-primer tail.
    TailFallback,
    /// Oriented from a full R1-primer Levenshtein search (thorough mode).
    LevenshteinFallback,
    /// Could not be oriented.
    None,
}

/// Strategy for reads with no qualifying homopolymer.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum FallbackMode {
    /// No fallback — mark such reads Unknown.
    Off,
    /// Cheap: windowed, allocation-free Hamming scan for the 10 bp R1-primer
    /// tail in the read ends. Default — no Levenshtein on the hot path.
    #[default]
    HammingTail,
    /// Thorough: full R1-primer Levenshtein search on both strands.
    Levenshtein,
}

/// Per-run orientation handling (CLI `--orientation`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, clap::ValueEnum)]
pub enum OrientMode {
    /// Detect the dominant orientation per file and enforce it.
    #[default]
    Auto,
    /// Force P5-led (no flip).
    P5Led,
    /// Force P7-led (flip).
    P7Led,
}

/// Tuning for [`orient`].
#[derive(Debug, Clone, Copy)]
pub struct OrientCfg {
    pub cb_len: usize,
    pub umi_len: usize,
    /// Minimum homopolymer run length to accept the dT-primer anchor.
    pub min_homopolymer: usize,
    /// Max mismatches for the 10 bp R1-tail bracket confirm.
    pub r1_tail_max_mm: u32,
    /// Max edits for the Levenshtein fallback R1-primer search.
    pub anchor_max_edits: u32,
    /// Strategy when no qualifying homopolymer is found.
    pub fallback: FallbackMode,
    /// Known run orientation.
    pub expected: Option<Orientation>,
}

impl OrientCfg {
    pub(crate) fn bracket(&self) -> usize {
        self.cb_len + self.umi_len
    }
}

/// Result of orienting one read to canonical P5→P7.
#[derive(Debug)]
pub struct Oriented<'a> {
    pub seq: Cow<'a, [u8]>,
    pub qual: Cow<'a, [u8]>,
    pub orientation: Orientation,
    pub method: OrientMethod,
    pub leader_present: bool,
    pub polyt_start: Option<usize>,
    pub cb_umi: Option<(usize, usize)>,
    pub r1_tail_confirmed: bool,
}
