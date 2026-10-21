//! Cross-gene UMI collision resolution: what happens when one
//! (cell, sample, corrected UMI) is claimed by more than one gene.
//!
//! | Policy                         | Outcome per claimant          | Counted as          |
//! |--------------------------------|-------------------------------|---------------------|
//! | default (single-cell)          | `Won` / `Lost`, or `Tie`      | 1 for the winner    |
//! | `--tagseq --no-fractional-...` | `Shared`                      | 1 for every claimant|
//! | `--fractional-collisions` (and default under `--tagseq`) | `Fractional` | `support / total` |
//!
//! Fractional weights for one collision group always sum to 1.0, so the
//! group contributes exactly one molecule to the matrix, split across genes
//! in proportion to read support.
// FxHash: these maps hold billions of fixed-size integer keys at real
// scale; SipHash's DoS resistance buys nothing here and cost ~6% of runtime.
use rustc_hash::FxHashMap as HashMap;

use super::dedup::Molecule;
use rayon::prelude::*;

use super::types::{Metrics, Outcome, RawKey};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Resolved {
    pub(crate) outcome: Outcome,
    pub(crate) degree: u32,
    pub(crate) weight: f32,
}

impl Resolved {
    /// Placeholder before `resolve_collisions` runs; every molecule is
    /// overwritten there (each root is in exactly one collision group).
    pub(crate) const UNRESOLVED: Self = Self {
        outcome: Outcome::Lost,
        degree: 0,
        weight: 0.0,
    };

    fn new(outcome: Outcome, degree: u32) -> Self {
        let weight = if outcome.is_counted() { 1.0 } else { 0.0 };
        Self {
            outcome,
            degree,
            weight,
        }
    }
}

/// Resolve every collision group and store each molecule's outcome in its
/// own `Molecule::resolved` (no separate root -> outcome table).
pub(crate) fn resolve_collisions(
    molecules: &mut HashMap<RawKey, Molecule>,
    metrics: &mut Metrics,
    tagseq: bool,
    fractional: bool,
) {
    eprintln!("[morselle molecules] phase=gene_collision_index start");
    // (cell, sample, umi): without `sample`, two different samples'
    // molecules sharing a cell+UMI would look like one real molecule
    // colliding across two genes. Sorted runs instead of a
    // HashMap<_, Vec<RawKey>>: no per-molecule Vec allocation (nearly every
    // group has one claimant). Claimant order never affects outcomes.
    let mut roots: Vec<RawKey> = molecules.keys().copied().collect();
    roots.par_sort_unstable_by_key(|r| (r.cell, r.sample, r.umi, r.gene));
    let same_umi = |a: &RawKey, b: &RawKey| (a.cell, a.sample, a.umi) == (b.cell, b.sample, b.umi);
    eprintln!(
        "[morselle molecules] phase=gene_collision_index complete collision_keys={}",
        roots.chunk_by(same_umi).count()
    );
    eprintln!("[morselle molecules] phase=gene_collision_resolve start");
    // Outcomes for one group are computed from shared reads of `molecules`,
    // then written back; `pending` is reused across groups.
    let mut pending: Vec<(RawKey, Resolved)> = Vec::new();
    for claimants in roots.chunk_by(same_umi) {
        let degree = claimants.len() as u32;
        if degree == 1 {
            set(molecules, claimants[0], Resolved::new(Outcome::Unique, 1));
            metrics.molecules_counted += 1;
            super::barcode_metrics::bump(metrics, claimants[0].cell, claimants[0].sample, 2);
            continue;
        }
        metrics.collision_groups += 1;
        super::barcode_metrics::bump(metrics, claimants[0].cell, claimants[0].sample, 3);
        pending.clear();
        if fractional {
            resolve_fractional(claimants, molecules, metrics, &mut pending);
        } else if tagseq {
            // Pooled-sample mode: a shared UMI is usually two genuinely
            // distinct genes' molecules colliding by chance in a 12bp UMI
            // space, not a mapping artifact -- count every claimant.
            for &root in claimants {
                pending.push((root, Resolved::new(Outcome::Shared, degree)));
                metrics.collision_molecules_shared += 1;
                metrics.molecules_counted += 1;
                super::barcode_metrics::bump(metrics, root.cell, root.sample, 2);
            }
        } else {
            resolve_winner_take_all(claimants, molecules, metrics, &mut pending);
        }
        for &(root, resolved) in &pending {
            set(molecules, root, resolved);
        }
    }
    eprintln!(
        "[morselle molecules] phase=gene_collision_resolve complete molecules={} collision_groups={}",
        metrics.molecules_counted, metrics.collision_groups
    );
}

fn set(molecules: &mut HashMap<RawKey, Molecule>, root: RawKey, resolved: Resolved) {
    molecules
        .get_mut(&root)
        .expect("collision claimant is a molecule")
        .resolved = resolved;
}

fn resolve_winner_take_all(
    claimants: &[RawKey],
    molecules: &HashMap<RawKey, Molecule>,
    metrics: &mut Metrics,
    out: &mut Vec<(RawKey, Resolved)>,
) {
    let degree = claimants.len() as u32;
    let best = claimants
        .iter()
        .map(|root| molecules[root].support)
        .max()
        .unwrap();
    let winners: Vec<RawKey> = claimants
        .iter()
        .copied()
        .filter(|root| molecules[root].support == best)
        .collect();
    if winners.len() == 1 {
        for &root in claimants {
            if root == winners[0] {
                out.push((root, Resolved::new(Outcome::Won, degree)));
                metrics.collision_molecules_won += 1;
                metrics.molecules_counted += 1;
                super::barcode_metrics::bump(metrics, root.cell, root.sample, 2);
            } else {
                out.push((root, Resolved::new(Outcome::Lost, degree)));
                metrics.collision_molecules_lost += 1;
            }
        }
    } else {
        metrics.collision_ties += 1;
        for &root in claimants {
            out.push((root, Resolved::new(Outcome::Tie, degree)));
        }
    }
}

fn resolve_fractional(
    claimants: &[RawKey],
    molecules: &HashMap<RawKey, Molecule>,
    metrics: &mut Metrics,
    out: &mut Vec<(RawKey, Resolved)>,
) {
    let degree = claimants.len() as u32;
    let total: u64 = claimants.iter().map(|root| molecules[root].support).sum();
    for &root in claimants {
        let weight = if total > 0 {
            (molecules[&root].support as f64 / total as f64) as f32
        } else {
            1.0 / degree as f32
        };
        out.push((
            root,
            Resolved {
                outcome: Outcome::Fractional,
                degree,
                weight,
            },
        ));
        metrics.collision_molecules_fractional += 1;
    }
    // The group's weights sum to 1.0: it is one molecule split across genes.
    metrics.molecules_counted += 1;
    super::barcode_metrics::bump(metrics, claimants[0].cell, claimants[0].sample, 2);
}
