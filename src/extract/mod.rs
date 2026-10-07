//! R1/R2 bracket extraction and adapter/homopolymer trimming (plan §4 step 4, M2).
//!
//! Operates on the canonical (P5→P7) read produced by M1, reusing the coordinates
//! it already found:
//! - **R1** = the 28 bp CB+UMI block located by [`crate::orient`].
//! - **R2** = cDNA between the dT-primer run end and the 3' read-through adapter
//!   (`AGATCGGAAGAGC`), found by a cheap Hamming seed; if the read ends mid-cDNA
//!   (no adapter), R2 runs to the read end. **Standard 3' capture only** --
//!   see [`ExtractCfg::flex_r2_boundary`] for 10x Flex, which bounds R2 by a
//!   real structural landmark ([`FLEX_R2_END_ANCHOR`]) instead.
//!
//! All positions are returned as ranges into the canonical sequence — zero-copy;
//! the caller slices `seq`/`qual`.
pub mod pipeline;
pub mod types;

#[cfg(test)]
mod flex_tests;
#[cfg(test)]
mod tests;

pub use pipeline::{dt_run_end, extract, find_adapter};
pub use types::{
    ExtractCfg, Extracted, FLEX_R2_END_ANCHOR, FLEX_R2_END_ANCHOR_MAX_MM,
    FLEX_R2_END_ANCHOR_MAX_SEARCH, FLEX_R2_END_ANCHOR_PREFIX_LEN, MIN_PARTIAL_ADAPTER,
    MIN_PARTIAL_TSO, TSO_RC,
};
