//! `bam-audit` — single-pass, O(1)-memory audit of a 10x reference BAM.
//!
//! Gathers everything needed for the R2 contamination audit and the long-read
//! mapping analysis, so a full (multi-hundred-GB) BAM can be summarized on the
//! server and only a small JSON moved back for analysis.
//!
//! Per BAM it reports, globally and per read-length bin:
//!   * mapping    — unique (NH==1), MAPQ 255, unmapped
//!   * annotation — RE tag region (Exonic / iNtronic / Intergenic), GX gene call,
//!                  multi-gene (ambiguous) GX
//!   * fidelity   — nM mismatches per 100 aligned bases
//!   * trimming   — 5'/3' soft clips, and WHAT they contain (adapter / TSO / polyA /
//!                  polyT / unexplained) using indel-robust 8-mer matching. Residual
//!                  adapter or polyT cannot align, so STAR soft-clips it; that is the
//!                  contamination signal.
//!   * molecules  — counted molecules (xf bit 8 representatives), reads/molecule,
//!                  and UMIs per read.
//!
//! Geometry note: sbx10x emits R2 reverse-complemented, so residual 3' read-through
//! adapter lands at the 5' end of R2 and residual polyT lands at the 3' end as polyA.
//! Both motifs are tested at BOTH ends so the result isn't biased by that expectation.
//!
//! Cross-BAM molecule concordance (does a long read call the same gene as a short
//! one?) would normally need a hash map over hundreds of millions of molecules.
//! Instead `--molecules` dumps a CONSISTENT hash-subsample: a molecule is emitted iff
//! `hash(CB,UMI) % --molecule-sample == 0`. The same molecules are therefore selected
//! in every BAM, so the dumps can simply be joined, at 1/N the size.
//!
//! Usage:
//!   bam-audit possorted_genome_bam.bam --label "SBX p5" -o p5_audit.json \
//!       --threads 8 --molecules p5_molecules.tsv.gz --molecule-sample 64
pub mod engine;
pub mod helpers;
pub mod types;

pub use engine::run;
pub use types::Cli;
