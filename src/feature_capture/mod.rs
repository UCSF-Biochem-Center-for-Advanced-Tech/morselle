//! 10x Feature Barcoding read -> feature assignment (the Flex design notes
//! Phase 4, generalized per the plan's own amendment).
//!
//! Not CRISPR- or Flex-specific by design: 10x's `<anchor>(BC)` read
//! structure + `feature_reference.csv` schema (`id,name,read,pattern,
//! sequence,feature_type`) is the same mechanism behind CRISPR Guide
//! Capture, Antibody Capture (CITE-seq), and Cell Multiplexing (CMO
//! hashing), and runs on ordinary (non-Flex) 3'/5' GEX libraries with no
//! probe barcode at all just as often as on Flex. `FeatureMatcher` only
//! does anchor+capture-barcode matching; it carries `feature_type` through
//! as data rather than branching on it — any per-`feature_type` calling
//! logic (e.g. Phase 5's CRISPR-specific GMM protospacer calling) is a
//! separate, later stage.
mod matcher;
mod pipeline;
mod types;

#[cfg(test)]
mod tests;

pub use matcher::*;
pub use pipeline::*;
pub use types::*;
