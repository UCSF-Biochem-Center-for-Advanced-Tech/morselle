//! `morselle extract` pipeline runner: convert raw Roche SBX reads into
//! per-sample 10x R1/R2 FASTQs (orientation, demux, extraction, gating, and
//! realignment).
pub(crate) mod barcode_metrics;
pub(crate) mod engine;
pub(crate) mod formatting;
pub(crate) mod orientation_scan;
pub(crate) mod process;
pub(crate) mod relock;
pub(crate) mod report;
#[cfg(test)]
mod tests;
pub(crate) mod types;
pub(crate) mod worker;

pub use engine::run;
