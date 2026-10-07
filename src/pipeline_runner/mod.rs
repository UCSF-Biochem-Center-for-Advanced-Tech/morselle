//! sbx10x — convert raw Roche SBX reads into per-sample 10x R1/R2 FASTQs.
//!
//! See `sbx10x_implementation_plan.md` for the milestone-phased design. M0 wires
//! up the crate, CLI, logging, and a FASTQ read->write passthrough; subsequent
//! milestones fill in orientation, demux, extraction, gating, and realignment.
pub(crate) mod engine;
pub(crate) mod formatting;
pub(crate) mod process;
pub(crate) mod report;
#[cfg(test)]
mod tests;
pub(crate) mod types;
pub(crate) mod worker;

pub use engine::run;
