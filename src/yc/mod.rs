//! Roche SBX `YC` duplex-consensus tag — encode and decode (plan: docs/archive/OVERLAP_PLAN.md §Phase 2).
//!
//! A `YC` tag rides along with a duplex consensus read (as a BAM auxiliary tag,
//! or as a FASTQ comment field) and encodes losslessly how the consensus differs
//! from the two strands it was built from. Given the consensus sequence, its
//! qualities, and the tag, both original strands can be reconstructed exactly.
//!
//! ```text
//! {left_tail}{delimiter}{duplex}{delimiter}{right_tail}
//! ```
//!
//! * `+` means that tail's simplex bases came from R1, `-` from R2. A tag always
//!   carries exactly two delimiters; `+` is the default when a tail is empty.
//! * Tails are bare integers — counts of simplex bases. Never letters.
//! * The duplex segment alternates integers (runs of concordant bases) with
//!   single-letter discordance codes.
//!
//! We use this format for Illumina paired-end overlap merges as well as SBX-D
//! reads: an overlap merge has exactly the same structure (two opposite strands
//! of one duplex fragment, a concordant middle, simplex tails).
//!
//! Written from the public specification at
//! <https://roche-axelios.gitbook.io/xoos/latest/tutorials/interpreting-sbx-reads/interpreting-sbx-d-reads>.
//! Deliberately NOT derived from the XOOS reference implementation, whose licence
//! restricts derivatives to data originating from Roche instruments.
pub mod codec;
pub mod tag;
pub mod types;

#[cfg(test)]
mod tests;

pub use codec::{decode, encode_overlap};
pub use types::{
    code_bases, complement_code, pair_code, parse_u64, BaseType, DecodedPair, DuplexOp,
    EncodeStats, YcError, YcTag, CONSENSUS_GAP_QUAL, GAP,
};
