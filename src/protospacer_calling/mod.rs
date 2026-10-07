//! CRISPR protospacer calling (docs/archive/FLEX_SUPPORT_PLAN.md Phase 5) — the
//! published-algorithm gap a raw UMI-count passthrough would miss.
//!
//! Cellranger's real published algorithm (`cr-crispr-algorithm`): "models
//! two cell populations for each guide RNA... fits a Gaussian Mixture Model
//! to the log-transformed distribution of Molecules/Cell," with "a minimum
//! of three guide RNA UMIs to be considered positive." Confirmed genuinely
//! adaptive against real CAT06 data (`protospacer_umi_thresholds.csv` shows
//! `CD45,3` for sample `Pos` but `CD45,6` for `RZ211`), which is why this
//! module fits a real 2-component Gaussian mixture rather than always
//! applying the floor.
//!
//! **Validated against real ground truth with one honestly-documented gap**:
//! this implementation reproduces `Pos`'s real threshold (3) exactly. It
//! does *not* reproduce `RZ211`'s (6) — investigated thoroughly (see
//! `fit_threshold`'s doc comment) rather than shrugged off: this module's
//! fit on `RZ211`'s real data scores a *higher* log-likelihood than a
//! reference implementation (`sklearn.mixture.GaussianMixture`) achieves on
//! the identical data, so the gap is a real difference between "the
//! maximum-likelihood 2-Gaussian-mixture fit" and whatever cellranger's
//! real, undocumented algorithm does -- not a bug or an unescaped local
//! optimum in this code. Pinned as a regression test, not silently
//! papered over.
//!
//! CRISPR-specific by design (unlike `feature_capture.rs`): Antibody/CMO
//! Capture use different downstream calling (normalization, sample
//! assignment), not this GMM.
mod gmm;
mod pipeline;
#[cfg(test)]
mod tests;
mod types;

pub use gmm::*;
pub use pipeline::*;
pub use types::*;
