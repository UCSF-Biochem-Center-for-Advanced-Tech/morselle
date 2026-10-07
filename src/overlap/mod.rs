//! Paired-end overlap merging for Illumina reads (plan: docs/archive/OVERLAP_PLAN.md §Phase 1).
//!
//! A port of HTStream's `hts_Overlapper` seed-and-extend algorithm, adapted for
//! the long-R1 10x GEX geometry.
//!
//! **Merge the raw mates, before any trimming.** The RT primer's dT block is the
//! polyA complement, so the molecule R1 reads and the molecule `rc(R2)` reads are
//! the same sequence:
//!
//! ```text
//! molecule = CB(16) | UMI(12) | T(30) | rc(cDNA)
//! R1       = molecule                        | adapter        (truncated)
//! R2       = rc(molecule)                    | rc(R1 primer)  (truncated)
//! rc(R2)   = R1 primer | CB | UMI | T(30) | rc(cDNA)
//! ```
//!
//! Merging before trimming is what makes the alignment gapless — strip CB+UMI+dT
//! from R1 first and the mates no longer describe the same span. It buys two more
//! things:
//!
//! * **Adapter removal by geometry, not by sequence matching.** The R1 primer
//!   sits *before* the CB, so `rc(R2)`'s 5' read-through is primer and `R1`'s 3'
//!   read-through is adapter — and both fall outside the overlap, so the merge
//!   drops them without ever consulting an adapter sequence. Each mate is
//!   independent evidence of where the fragment ends, which a Hamming-seeded
//!   trimmer like [`crate::extract::find_adapter`] can only approximate; it has
//!   a false-trim rate, and this has none. Running a matching trimmer over a
//!   merged read can therefore only add false positives.
//! * **Duplex confirmation of the barcode.** On short fragments R2 reads back
//!   through the UMI and CB, so CB+UMI lands inside the overlap and comes out
//!   with summed quality — and 28 bp of complex, per-read-unique seed material.
//!
//! This module takes the **raw** R1 and R2 and returns the merged molecule in
//! **R1 orientation** with a [`crate::yc`] tag describing the merge losslessly.
//! The caller then locates CB/UMI/cDNA in the merged read and trims the tag to
//! match.
//!
//! Deliberate divergences from `hts_Overlapper`, all covered by tests:
//!
//! * **Seeding is decoupled from the merge frame.** HTStream swaps the mates so
//!   R1 is the longer one, which also flips the merge frame. We cannot flip the
//!   frame — the merge must stay in R1's orientation or the `YC` tag's notion of
//!   which strand is R1 silently inverts. So we choose which mate to *index*
//!   independently: index the longer, probe the shorter. The shorter read's ends
//!   are always inside the overlap; `rc(R2)`'s need not be. In 10x geometry R1
//!   loses 28 bp of CB+UMI and the whole dT run, so its insert is routinely
//!   *contained* within `rc(R2)`, with overhangs on both sides — probing
//!   `rc(R2)`'s ends would find nothing at all.
//! * **Correct 3' extension offset.** HTStream appends `seq2.substr(maxLoop, ..)`,
//!   which is only right when the `rc(R2)` offset is zero — its swap normally
//!   guarantees that. Without the swap we must use `loc2 + maxLoop`.
//! * **Discordant quality** is `max(1, |q1 - q2|)`. HTStream writes
//!   `min(q1 - q2 + offset, 1 + offset)`, collapsing every disagreement to Q1.
//! * **Mismatch budget** uses [`crate::anchors::hamming_leq`] (`<= max`), which
//!   early-exits correctly. HTStream's `threshold_mismatches` decrements before
//!   testing and underflows `size_t` when the budget is zero.
//! * **Seed masking.** k-mers containing a non-ACGT base, or lying inside a
//!   homopolymer, are kept out of the index — the dT run abutting the R1 insert
//!   would otherwise seed against any poly-A in R2.
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
