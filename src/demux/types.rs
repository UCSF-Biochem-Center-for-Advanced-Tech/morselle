/// A single index, for real-distance measurement (`min_pairwise_hamming`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum IndexSource {
    I5,
    I7,
}

/// How to combine i5/i7 for demux assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum DemuxMode {
    /// i5 OR i7: assigned when either matches.
    Dual,
    /// i5 AND i7: assigned only when both are found and agree.
    Both,
    /// i5 alone; i7 is never even consulted.
    I5,
    /// i7 alone; i5 is never even consulted.
    I7,
}

/// Tuning for demultiplexing.
#[derive(Debug, Clone, Copy)]
pub struct DemuxCfg {
    /// Max mismatches when matching an index.
    pub index_max_mm: u32,
    /// `None` = auto mode.
    pub mode: Option<DemuxMode>,
}

/// How a read was assigned (for accounting).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DemuxVia {
    /// Both indices present and agree.
    Dual,
    /// Only i5 available (i7 truncated/absent).
    I5Only,
    /// Only i7 available (i5 truncated/absent).
    I7Only,
    /// Both present but disagree.
    Conflict,
    /// No usable index.
    None,
}

/// Demux outcome: assigned sample index + how.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Outcome {
    pub sample: Option<usize>,
    pub via: DemuxVia,
}
