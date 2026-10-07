/// Reverse-complement of [`crate::anchors::TSO`], as the template-switch oligo
/// appears in the canonical frame: between the cDNA and the 3' read-through
/// adapter.
pub const TSO_RC: &[u8] = b"CCCATGTACTCTGCGTTGATACCACTGCTT";

/// Shortest partial TSO trimmed at the 3' end of a read that never reached an
/// adapter.
pub const MIN_PARTIAL_TSO: usize = 12;

/// 10x Flex GEX/CRISPR construct: the constant sequence immediately following
/// `Ligated_Probe_insert` (the 50bp gene/guide probe), in canonical (pre-RC)
/// coordinates.
pub const FLEX_R2_END_ANCHOR: &[u8] = b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC";
pub const FLEX_R2_END_ANCHOR_PREFIX_LEN: usize = 15;
pub const FLEX_R2_END_ANCHOR_MAX_MM: u32 = 2;

/// Max bp past `cdna_start` searched for [`FLEX_R2_END_ANCHOR`].
pub const FLEX_R2_END_ANCHOR_MAX_SEARCH: usize = 150;

/// Minimum 3' read-through adapter overlap to trim a partial adapter at the read end.
pub const MIN_PARTIAL_ADAPTER: usize = 7;

/// Tuning for [`extract`].
#[derive(Debug, Clone, Copy)]
pub struct ExtractCfg {
    /// Minimum cDNA length to emit an R2.
    pub min_cdna_len: usize,
    /// Max mismatches for the 3' adapter trim seed.
    pub adapter_max_mm: u32,
    /// Trim the template-switch oligo when it sits between the cDNA and the 3' adapter.
    pub trim_tso: bool,
    /// Cap emitted R2 length; 0 = uncapped.
    pub max_r2_len: usize,
    /// 10x Flex mode: require [`FLEX_R2_END_ANCHOR`] to bound R2's far end.
    pub flex_r2_boundary: bool,
    /// When `flex_r2_boundary` is on and [`FLEX_R2_END_ANCHOR`] can't be located.
    pub flex_min_cdna_len_no_anchor: usize,
}

/// Extraction result: ranges `[start, end)` into the canonical sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extracted {
    /// CB+UMI range (R1), 28 bp, or `None` if not located.
    pub r1: Option<(usize, usize)>,
    /// cDNA range (R2) in the canonical (polyT) frame.
    pub r2: Option<(usize, usize)>,
    /// End of the dT-primer run = cDNA start.
    pub polyt_end: Option<usize>,
    /// Length of the dT-primer run actually observed (QC).
    pub polyt_run_len: u16,
    /// Position of the 3' read-through adapter when found.
    pub adapter_pos: Option<usize>,
    /// The 3' read-through adapter was found and trimmed.
    pub adapter_trimmed: bool,
    /// Where the reverse-complemented TSO starts, when found and trimmed.
    pub tso_pos: Option<usize>,
    /// The TSO was found abutting the cDNA's 3' end and trimmed.
    pub tso_trimmed: bool,
}
