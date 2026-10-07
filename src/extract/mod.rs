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
