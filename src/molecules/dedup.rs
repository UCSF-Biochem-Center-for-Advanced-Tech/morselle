//! UMI correction and cross-gene collision resolution.
//!
//! Restored from the pre-0.5.7 `sbx10x-molecules` implementation (commit
//! `4d42a18`), with fractional collision resolution layered on top as an
//! opt-in fourth policy. See `resolve_collisions` for the policy table.
use std::collections::HashMap;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use rayon::prelude::*;
use tracing::info;

use super::saturation::Saturation;
use super::types::{Decision, Dna, Evidence, Metrics, RawKey};

pub fn correct_group(
    mut keys: Vec<RawKey>,
    evidence: &HashMap<RawKey, Evidence>,
) -> HashMap<RawKey, RawKey> {
    keys.sort_by_key(|key| (std::cmp::Reverse(evidence[key].reads), key.umi));
    // Every raw UMI in a (cell, gene) group is already distinct (RawKey
    // includes the UMI, and `evidence` is keyed on the full RawKey), so
    // indexing by UMI once and generating each child's 3*len single-mismatch
    // neighbors directly replaces an O(n^2) all-pairs Hamming scan with O(n)
    // hash lookups. This matters because raw UMI diversity per (cell, gene)
    // is heavily skewed toward highly expressed genes, where n^2 previously
    // dominated the whole molecules stage.
    let by_umi: HashMap<Dna, RawKey> = keys.iter().map(|&k| (k.umi, k)).collect();
    let mut roots = HashMap::with_capacity(keys.len());
    for &child in &keys {
        let mut parent = child;
        let threshold = evidence[&child].reads.saturating_mul(2).saturating_sub(1);
        for neighbor in child.umi.one_mismatch_neighbors() {
            let Some(&candidate) = by_umi.get(&neighbor) else {
                continue;
            };
            if evidence[&candidate].reads < threshold {
                continue;
            }
            if (evidence[&candidate].reads, std::cmp::Reverse(candidate.umi))
                > (evidence[&parent].reads, std::cmp::Reverse(parent.umi))
            {
                parent = candidate;
            }
        }
        let root = roots.get(&parent).copied().unwrap_or(parent);
        roots.insert(child, root);
    }
    roots
}

pub fn merge_roots(
    mut left: HashMap<RawKey, RawKey>,
    mut right: HashMap<RawKey, RawKey>,
) -> HashMap<RawKey, RawKey> {
    if left.len() < right.len() {
        std::mem::swap(&mut left, &mut right);
    }
    left.extend(right.drain());
    left
}

/// One corrected UMI node (a molecule candidate) after aggregation.
#[derive(Clone, Copy)]
pub(crate) struct Molecule {
    pub(crate) support: u64,
    pub(crate) representative: u64,
    /// MAPQ (high 8 bits) followed by aligned bases (low 24 bits). Preserves
    /// the representative ordering without growing this large, temporary
    /// hash-map value.
    pub(crate) best_rank: u32,
    pub(crate) min_read_hash: u32,
}

fn rank(item: &Evidence) -> u32 {
    (u32::from(item.best_mapq) << 24) | item.best_aligned.min(0x00ff_ffff)
}

/// Directional UMI correction within each (cell, sample, gene) group.
///
/// Grouped by (cell, sample, gene), not just (cell, gene): the same CB
/// whitelist is shared across every sample in a pooled Flex run, so omitting
/// `sample` would let `correct_group()`'s by-UMI map silently overwrite one
/// sample's evidence with another's whenever two samples share a UMI value.
///
/// Under `--tagseq` correction is skipped entirely: at pooled-sample UMI
/// diversity two genuinely distinct molecules land 1 Hamming-distance apart
/// in a 12bp space by chance often enough that "low-count neighbour =
/// sequencing error" no longer holds. Trade-off: true 1bp UMI errors are not
/// corrected either, so counts trend slightly high rather than low.
fn correct_umis(
    evidence: &HashMap<RawKey, Evidence>,
    groups: HashMap<(Dna, u32, u32), Vec<RawKey>>,
    workers: usize,
    tagseq: bool,
) -> HashMap<RawKey, RawKey> {
    if tagseq {
        // Pooled-sample mode: directional UMI correction assumes a
        // low-count 1-mismatch neighbor is a sequencing error of a nearby
        // high-count node. At pooled-sample UMI diversity that assumption
        // breaks down the same way cross-gene collision does -- two
        // genuinely distinct molecules can land 1 Hamming-distance apart
        // in a 12bp space purely by chance and get merged into one. Skip
        // correction entirely: every raw (already exact-duplicate-collapsed)
        // node is its own molecule. Trade-off: genuine 1bp sequencing
        // errors in the UMI itself are no longer corrected either, so
        // molecule counts trend slightly high rather than slightly low.
        eprintln!(
            "[sbxpress molecules] phase=umi_correct skipped (--tagseq) raw_nodes={}",
            evidence.len()
        );
        return evidence.keys().map(|&key| (key, key)).collect();
    }
    if workers <= 1 {
        eprintln!(
            "[sbxpress molecules] phase=umi_correct mode=serial start groups={}",
            groups.len()
        );
        return groups
            .into_values()
            .map(|keys| correct_group(keys, evidence))
            .fold(HashMap::new(), merge_roots);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .expect("creating UMI-correction worker pool");
    info!(workers, groups = groups.len(), "correcting UMI groups in parallel");
    eprintln!(
        "[sbxpress molecules] phase=umi_correct start workers={} groups={}",
        workers,
        groups.len()
    );
    let total_groups = groups.len();
    let completed = Arc::new(AtomicUsize::new(0));
    let checkpoint = (total_groups / 100).max(1_000);
    pool.install(|| {
        groups
            .into_par_iter()
            .map(|(_, keys)| {
                let roots = correct_group(keys, evidence);
                let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                if done % checkpoint == 0 || done == total_groups {
                    eprintln!("[sbxpress molecules] phase=umi_correct groups={done}/{total_groups}");
                }
                roots
            })
            .reduce(HashMap::new, merge_roots)
    })
}

/// Sum each corrected root's read support and pick its representative read.
fn aggregate_molecules(
    evidence: &HashMap<RawKey, Evidence>,
    roots: &HashMap<RawKey, RawKey>,
) -> HashMap<RawKey, Molecule> {
    let mut molecules = HashMap::<RawKey, Molecule>::new();
    for (&key, item) in evidence {
        let root = roots[&key];
        molecules
            .entry(root)
            .and_modify(|m| {
                m.support += item.reads as u64;
                m.min_read_hash = m.min_read_hash.min(item.min_read_hash);
                let r = rank(item);
                // On an exact rank tie (common for true PCR/optical
                // duplicates), fall back to `best_ordinal` as an explicit,
                // deterministic tiebreak instead of HashMap iteration order.
                if r > m.best_rank || (r == m.best_rank && item.best_ordinal < m.representative) {
                    m.representative = item.best_ordinal;
                    m.best_rank = r;
                }
            })
            .or_insert(Molecule {
                support: item.reads as u64,
                representative: item.best_ordinal,
                best_rank: rank(item),
                min_read_hash: item.min_read_hash,
            });
    }
    molecules
}

pub fn build_decisions(
    evidence: &HashMap<RawKey, Evidence>,
    metrics: &mut Metrics,
    saturation_interval: Option<u64>,
    workers: usize,
    tagseq: bool,
    fractional: bool,
) -> (HashMap<RawKey, Decision>, Option<Saturation>) {
    // Grouped by (cell, sample, gene), not just (cell, gene): the same CB
    // whitelist is shared across every sample in a pooled Flex run, so
    // omitting `sample` here would let correct_group()'s internal by-UMI
    // map (keyed on UMI alone within a group) silently overwrite one
    // sample's evidence with another's whenever two samples' reads for the
    // same cell+gene happen to share a UMI value.
    let mut groups = HashMap::<(Dna, u32, u32), Vec<RawKey>>::new();
    for &key in evidence.keys() {
        groups
            .entry((key.cell, key.sample, key.gene))
            .or_default()
            .push(key);
    }
    eprintln!(
        "[sbxpress molecules] phase=umi_group_index complete groups={} raw_nodes={}",
        groups.len(),
        evidence.len()
    );
    let mut saturation =
        saturation_interval.map(|interval| Saturation::new(metrics.eligible_records, interval));
    if let Some(curve) = saturation.as_mut() {
        for (&(cell, _sample, gene), keys) in &groups {
            let min_hash = keys
                .iter()
                .map(|key| evidence[key].min_read_hash)
                .min()
                .expect("non-empty UMI group");
            curve.add_cell_gene(cell, gene, min_hash);
        }
    }

    let roots = correct_umis(evidence, groups, workers, tagseq);
    eprintln!("[sbxpress molecules] phase=umi_correct complete roots={}", roots.len());
    // `+=`, not `=`: the sharded path calls build_decisions() once per shard
    // against the same accumulating Metrics.
    metrics.raw_cell_gene_umis += evidence.len() as u64;
    metrics.umi_nodes_collapsed += roots.iter().filter(|(key, root)| key.umi != root.umi).count() as u64;

    eprintln!("[sbxpress molecules] phase=molecule_aggregate start");
    let molecules = aggregate_molecules(evidence, &roots);
    metrics.corrected_umi_nodes += molecules.len() as u64;
    metrics.molecules_before_gene_resolution += molecules.len() as u64;
    eprintln!(
        "[sbxpress molecules] phase=molecule_aggregate complete corrected_nodes={}",
        molecules.len()
    );

    let root_outcomes = super::collisions::resolve_collisions(&molecules, metrics, tagseq, fractional);

    if let Some(curve) = saturation.as_mut() {
        for (&root, molecule) in &molecules {
            if root_outcomes[&root].outcome.is_counted() {
                curve.add_molecule(root.cell, molecule.min_read_hash);
            }
        }
    }
    eprintln!("[sbxpress molecules] phase=decision_materialize start");
    let decisions: HashMap<RawKey, Decision> = evidence
        .keys()
        .map(|&key| {
            let root = roots[&key];
            let resolved = root_outcomes[&root];
            (
                key,
                Decision {
                    corrected: root.umi,
                    representative: molecules[&root].representative,
                    collision_degree: resolved.degree,
                    outcome: resolved.outcome,
                    weight: resolved.weight,
                },
            )
        })
        .collect();
    eprintln!(
        "[sbxpress molecules] phase=decision_materialize complete decisions={}",
        decisions.len()
    );
    (decisions, saturation)
}
