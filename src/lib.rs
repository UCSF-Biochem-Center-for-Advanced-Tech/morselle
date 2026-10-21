//! MORSELLE — MOlecular Read Structure, Extraction, Labeling & Enumeration.
//!
//! Library behind the `morselle` command: converts raw Roche SBX long reads and
//! Illumina paired-end reads into per-sample 10x-style FASTQs, annotated BAMs,
//! molecule counts, and expression matrices.

pub mod accounting;
pub mod anchors;
pub mod annotate;
/// Compatibility with the public module path used before 0.5.7.
pub mod annotate_bam {
    pub use crate::annotate::{run, Cli};
}
pub mod bam_audit;
pub mod bam_helpers;
pub mod barcode;
pub mod barcode_metrics;
pub mod cell_calling;
pub mod complexity;
pub mod config;
pub mod cr_bam_to_mtx;
pub mod demux;
pub mod extract;
pub mod feature_capture;
pub mod gzin;
pub mod hamming_index;
pub mod homopolymer;
pub mod illumina;
pub mod io;
pub mod molecules;
pub mod orient;
pub mod overlap;
pub mod pipeline;
pub mod pipeline_runner;
pub mod probe_barcode;
pub mod probe_demux;
pub mod probe_sample_sheet;
pub mod probe_set;
pub mod protospacer_calling;
pub mod sample_sheet;
pub mod yc;

pub use pipeline_runner::run;
