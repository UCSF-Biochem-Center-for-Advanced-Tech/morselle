//! sbx10x — convert raw Roche SBX reads into per-sample 10x R1/R2 FASTQs.
//!
//! See `sbx10x_implementation_plan.md` for the milestone-phased design. M0 wires
//! up the crate, CLI, logging, and a FASTQ read->write passthrough; subsequent
//! milestones fill in orientation, demux, extraction, gating, and realignment.

pub mod accounting;
pub mod anchors;
pub mod annotate;
pub mod annotate_bam;
pub mod bam_audit;
pub mod bam_helpers;
pub mod barcode;
pub mod cell_calling;
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
