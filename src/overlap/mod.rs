pub mod alignment;
pub mod types;

#[cfg(test)]
mod realistic_tests;
#[cfg(test)]
mod tests;

pub use alignment::{check_overlap, overlap, overlap_in};
pub use types::{
    base_code, build_index, seed_positions, usable_kmer, InsertClass, Merged, OverlapCfg, Scratch,
    QUAL_OFFSET,
};
