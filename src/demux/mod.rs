//! i5/i7 extraction + sample-sheet matching (plan §4 step 3, M3).
//!
//! Empirically (validated against the 14 known sample indices), in a canonical
//! P5→P7 read:
//! - **i5** sits 10 bp immediately 5' of the R1 primer, as `revcomp(sheet i5)`.
//! - **i7** sits 10 bp after `rc(R2 primer)` (i.e. after the 3' adapter), as the
//!   sheet i7 forward.
//!
//! Both positions are derived from coordinates M1/M2 already found, so demux adds
//! no new search. Indices are matched via a precomputed ≤`index_max_mm` neighbor
//! map (O(1) lookup, with ambiguity detection). The sample index set here has min
//! pairwise Hamming i7=3 / i5=4, so single-index assignment (when the other index
//! is truncated — common for i7 on long reads) is unambiguous.
mod engine;
mod indel;
mod types;

#[cfg(test)]
mod indel_tests;
#[cfg(test)]
mod tests;

pub use engine::*;
pub use indel::{edit_distance, min_pairwise_edit, IndelSafety, IndexMatch};
pub use types::*;
