//! Assign cell barcodes to samples by their declared probe-barcode set
//! (FLEX_SUPPORT_PLAN.md Phase 2).
//!
//! Published-algorithm parity notes (see FLEX_SUPPORT_PLAN.md "Where this
//! plan corrects the initial design"): assignment sums molecule counts across
//! a sample's *whole declared set* (which can be >1 probe barcode — real
//! `config.csv` example: `Pos,D-A01|D-A02`) before comparing samples, not a
//! naive top-single-barcode vote.
//!
//! Composite cell barcode (verified against 10x's own "Cell Ranger's Flex
//! Algorithm" support doc, "Multiplex workflow" + "High occupancy GEM
//! filtering"): in real Flex multiplex pooling, the true cell barcode is the
//! *(GEM barcode, Probe Barcode)* pair, not the GEM barcode alone — multiple
//! fixed, permeabilized samples are pooled before loading, so it is normal
//! and expected for one physical GEM droplet to legitimately hold cells from
//! more than one sample. `assign` therefore emits one `Assigned` row *per
//! sample* whose declared set clears `min_molecules` for a GEM (not a single
//! winner-take-all call), and only discards a GEM outright (`HighOccupancyGem`)
//! when its number of clearing samples exceeds a negative-binomial-fit
//! 99.9th-percentile cutoff on the observed per-GEM occupancy distribution —
//! mirroring 10x's real "High occupancy GEM filtering" step, which estimates
//! NB(r, p=0.5) from the data (r from the observed mean, p fixed at 0.5) and
//! discards only that rare chip-loading-artifact tail (their own worked
//! examples: >7 probe barcodes/GEM for a 16-plex run, >19 for 96-plex).
//! Earlier revision of this module treated *any* second sample clearing the
//! floor as an all-or-nothing reject — real CAT06 data showed that fires on
//! 60%+ of GEMs with enough signal to assign, which cannot be real
//! multiplet/doublet biology at that rate and was discarding the majority of
//! otherwise-real, correctly multiplexed cells.
mod engine;
mod types;

#[cfg(test)]
mod tests;

pub use engine::*;
pub use types::*;
