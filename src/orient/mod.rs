//! Per-read orientation → canonical P5→P7 (plan §4 step 2, M1).
//!
//! Homopolymer-first: the 3' GEX dT-primer run (poly-T in canonical, poly-A when
//! the read is P7-led/reversed) is found by a cheap O(n) scan and used to both
//! orient the read and bracket the CB+UMI block. A fixed-offset 10 bp Hamming
//! check of the R1-primer tail confirms the bracket for free. Only when no
//! qualifying homopolymer is found do we fall back to a `triple_accel` Levenshtein
//! search for the R1 primer (the expensive path, kept off the common case).
//!
//! Scope: single 3' assay. The 5' (TSO-spacer) path is stubbed — see
//! [`orient_5p_stub`].
//!
//! Cost note: `levenshtein_search` is deliberately *not* on the default path.
//! Benchmarks (300k real reads/run) show the windowed Hamming R1-tail fallback
//! recovers as many or more reads at ~1/20th the cost — see [`FallbackMode`].
mod matcher;
mod types;

#[cfg(test)]
mod tests;

pub use matcher::*;
pub use types::*;
