//! 10x Flex probe-set alignment (docs/archive/FLEX_SUPPORT_PLAN.md Phase 6 prerequisite):
//! assigns a Flex GEX read to a gene by direct lookup of its probe sequence
//! against the real probe-set reference, **not** genome alignment.
//!
//! Flex doesn't sequence free cDNA: each read's payload is (approximately)
//! one of a fixed, closed set of ~54k known 50bp probe sequences (one of a
//! few probes per targeted gene — see `Chromium_Human_Transcriptome_Probe_Set*.csv`).
//! Confirmed empirically this session (see docs/archive/FLEX_SUPPORT_PLAN.md Phase 6):
//! full-genome splice-aware alignment of real CAT06 Flex reads found almost
//! nothing (59 of 1.65M records) — the payload isn't cDNA a genome mapper
//! can place. This module is the real substitute: a windowed, Hamming-tolerant
//! lookup of the probe sequence directly against the reference, mirroring
//! `probe_barcode.rs`'s own windowed-locate design (same shared
//! `HammingIndex` primitive), and emitting a BAM shaped exactly like
//! `feature_capture.rs`'s (synthetic per-gene contigs, `GX`/`GN` tags) so it
//! feeds `molecules`/`matrix`/`cells`/`probe-demux` completely unchanged.
mod aligner;
mod pipeline;
mod types;

#[cfg(test)]
mod tests;

pub use aligner::*;
pub use pipeline::*;
pub use types::*;
