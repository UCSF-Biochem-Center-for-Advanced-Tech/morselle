//! Per-read pipeline state machine and core types (plan §3, §4).
//!
//! M0: only the types and a passthrough are defined. Stages (orient, demux,
//! extract, gate) are wired in M1-M5.

/// A borrowed view of one raw FASTQ record.
pub struct RawRead<'a> {
    pub id: &'a [u8],
    pub seq: &'a [u8],
    pub qual: &'a [u8],
}

/// Read orientation relative to the canonical P5->P7 layout.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Orientation {
    /// 5' library end first (canonical after any flip).
    P5Led,
    /// 3' library end first (will be reverse-complemented to canonical).
    P7Led,
    Unknown,
}

/// Parse completeness tier (plan §4 step 6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReadClass {
    Full,
    Partial,
    Unparseable,
}

/// Cell-barcode validation status (plan §4 step 5).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CbStatus {
    OnWhitelist,
    Within1,
    Indel,
    Fail,
}

/// Identifies a sample after demultiplexing.
pub type SampleId = String;

/// The result of processing one read.
pub struct Parsed {
    pub sample: Option<SampleId>,
    pub r1: Option<(Vec<u8>, Vec<u8>)>,
    pub r2: Option<(Vec<u8>, Vec<u8>)>,
    pub class: ReadClass,
    pub cb_status: CbStatus,
    pub bracket_len: Option<u8>,
    pub flags: u16,
}

// Flag bits for `Parsed::flags`.
pub const FLAG_LEADER_PRESENT: u16 = 1 << 0;
pub const FLAG_REALIGNED: u16 = 1 << 1;
pub const FLAG_FLIPPED: u16 = 1 << 2;
