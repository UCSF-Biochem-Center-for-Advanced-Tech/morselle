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
use std::collections::HashMap;

use super::dedup::Molecule;
use super::types::{Dna, Metrics, Outcome, RawKey};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct Resolved {
    pub(crate) outcome: Outcome,
    pub(crate) degree: u32,
    pub(crate) weight: f32,
}

impl Resolved {
    fn new(outcome: Outcome, degree: u32) -> Self {
        let weight = if outcome.is_counted() { 1.0 } else { 0.0 };
        Self { outcome, degree, weight }
    }
}

pub(crate) fn resolve_collisions(
    molecules: &HashMap<RawKey, Molecule>,
    metrics: &mut Metrics,
    tagseq: bool,
    fractional: bool,
) -> HashMap<RawKey, Resolved> {
    eprintln!("[sbxpress molecules] phase=gene_collision_index start");
    // (cell, sample, umi): without `sample`, two different samples'
    // molecules sharing a cell+UMI would look like one real molecule
    // colliding across two genes.
    let mut collisions = HashMap::<(Dna, u32, Dna), Vec<RawKey>>::new();
    for &root in molecules.keys() {
        collisions
            .entry((root.cell, root.sample, root.umi))
            .or_default()
            .push(root);
    }
    eprintln!(
        "[sbxpress molecules] phase=gene_collision_index complete collision_keys={}",
        collisions.len()
    );
    eprintln!("[sbxpress molecules] phase=gene_collision_resolve start");
    let mut out = HashMap::<RawKey, Resolved>::with_capacity(molecules.len());
    for claimants in collisions.values() {
        let degree = claimants.len() as u32;
        if degree == 1 {
            out.insert(claimants[0], Resolved::new(Outcome::Unique, 1));
            metrics.molecules_counted += 1;
            continue;
        }
        metrics.collision_groups += 1;
        if fractional {
            resolve_fractional(claimants, molecules, metrics, &mut out);
        } else if tagseq {
            // Pooled-sample mode: a shared UMI is usually two genuinely
            // distinct genes' molecules colliding by chance in a 12bp UMI
            // space, not a mapping artifact -- count every claimant.
            for &root in claimants {
                out.insert(root, Resolved::new(Outcome::Shared, degree));
                metrics.collision_molecules_shared += 1;
                metrics.molecules_counted += 1;
            }
        } else {
            resolve_winner_take_all(claimants, molecules, metrics, &mut out);
        }
    }
    eprintln!(
        "[sbxpress molecules] phase=gene_collision_resolve complete molecules={} collision_groups={}",
        metrics.molecules_counted, metrics.collision_groups
    );
    out
}

fn resolve_winner_take_all(
    claimants: &[RawKey],
    molecules: &HashMap<RawKey, Molecule>,
    metrics: &mut Metrics,
    out: &mut HashMap<RawKey, Resolved>,
) {
    let degree = claimants.len() as u32;
    let best = claimants.iter().map(|root| molecules[root].support).max().unwrap();
    let winners: Vec<RawKey> = claimants
        .iter()
        .copied()
        .filter(|root| molecules[root].support == best)
        .collect();
    if winners.len() == 1 {
        for &root in claimants {
            if root == winners[0] {
                out.insert(root, Resolved::new(Outcome::Won, degree));
                metrics.collision_molecules_won += 1;
                metrics.molecules_counted += 1;
            } else {
                out.insert(root, Resolved::new(Outcome::Lost, degree));
                metrics.collision_molecules_lost += 1;
            }
        }
    } else {
        metrics.collision_ties += 1;
        for &root in claimants {
            out.insert(root, Resolved::new(Outcome::Tie, degree));
        }
    }
}

fn resolve_fractional(
    claimants: &[RawKey],
    molecules: &HashMap<RawKey, Molecule>,
    metrics: &mut Metrics,
    out: &mut HashMap<RawKey, Resolved>,
) {
    let degree = claimants.len() as u32;
    let total: u64 = claimants.iter().map(|root| molecules[root].support).sum();
    for &root in claimants {
        let weight = if total > 0 {
            (molecules[&root].support as f64 / total as f64) as f32
        } else {
            1.0 / degree as f32
        };
        out.insert(root, Resolved { outcome: Outcome::Fractional, degree, weight });
        metrics.collision_molecules_fractional += 1;
    }
    // The group's weights sum to 1.0: it is one molecule split across genes.
    metrics.molecules_counted += 1;
}
