//! UMI correction and cross-gene collision resolution.
//!
//! Restored from the pre-0.5.7 `morselle-molecules` implementation (commit
//! `4d42a18`), with fractional collision resolution layered on top as an
//! opt-in fourth policy. See `resolve_collisions` for the policy table.
// FxHash: these maps hold billions of fixed-size integer keys at real
// scale; SipHash's DoS resistance buys nothing here and cost ~6% of runtime.
use rustc_hash::FxHashMap as HashMap;
use std::collections::HashMap as StdHashMap;
use std::hash::BuildHasher;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use rayon::prelude::*;
use tracing::info;

use super::collisions::Resolved;
use super::saturation::Saturation;
use super::types::{Decision, Dna, Evidence, Metrics, RawKey};
use super::umi_rules::{fold_deep_singletons, one_indel_neighbors, UmiRules};

pub fn correct_group<S: BuildHasher>(
    keys: Vec<RawKey>,
    evidence: &StdHashMap<RawKey, Evidence, S>,
) -> HashMap<RawKey, RawKey> {
    correct_group_with(keys, evidence, UmiRules::default()).0
}

/// `correct_group` plus the opt-in SBX rules in `umi_rules`: one-indel
/// neighbours as extra correction candidates, then folding 1-read molecules
/// of deep families. Returns the roots and how many singletons were folded.
pub fn correct_group_with<S: BuildHasher>(
    mut keys: Vec<RawKey>,
    evidence: &StdHashMap<RawKey, Evidence, S>,
    rules: UmiRules,
) -> (HashMap<RawKey, RawKey>, u64) {
    keys.sort_by_key(|key| (std::cmp::Reverse(evidence[key].reads), key.umi));
    // Every raw UMI in a (cell, gene) group is already distinct (RawKey
    // includes the UMI, and `evidence` is keyed on the full RawKey), so
    // indexing by UMI once and generating each child's 3*len single-mismatch
    // neighbors directly replaces an O(n^2) all-pairs Hamming scan with O(n)
    // hash lookups. This matters because raw UMI diversity per (cell, gene)
    // is heavily skewed toward highly expressed genes, where n^2 previously
    // dominated the whole molecules stage.
    let by_umi: HashMap<Dna, RawKey> = keys.iter().map(|&k| (k.umi, k)).collect();
    let mut roots = HashMap::with_capacity_and_hasher(keys.len(), Default::default());
    for &child in &keys {
        let mut parent = child;
        let threshold = evidence[&child].reads.saturating_mul(2).saturating_sub(1);
        let indels = rules
            .indels
            .then(|| one_indel_neighbors(child.umi))
            .into_iter()
            .flatten();
        for neighbor in child.umi.one_mismatch_neighbors().chain(indels) {
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
    let folded = fold_deep_singletons(&keys, evidence, &mut roots, rules.fold_singletons_min_top);
    (roots, folded)
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
    /// Collision outcome, filled in by `resolve_collisions` (kept here
    /// rather than in a second root-keyed table to cut peak memory).
    pub(crate) resolved: Resolved,
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
///
/// Exception, `--position-aware-umis`: the UMI carries the read's 5'
/// position tag, so 1-mismatch neighbours are only ever looked up at the
/// same position, where a 1bp-different UMI is almost always an error of
/// the same molecule. Correction then runs even under `--tagseq`.
/// (cell, sample, gene): the UMI-correction group of a raw key.
fn same_umi_group(a: &RawKey, b: &RawKey) -> bool {
    (a.cell, a.sample, a.gene) == (b.cell, b.sample, b.gene)
}

/// Roots for one (cell, sample, gene) group. Most groups hold a single raw
/// UMI, whose root is itself -- skip `correct_group`'s allocations for them.
/// Returns how many singletons `rules` folded.
fn correct_run<S: BuildHasher>(
    run: &[RawKey],
    evidence: &StdHashMap<RawKey, Evidence, S>,
    roots: &mut HashMap<RawKey, RawKey>,
    rules: UmiRules,
) -> u64 {
    if let [only] = run {
        roots.insert(*only, *only);
        0
    } else {
        let (group, folded) = correct_group_with(run.to_vec(), evidence, rules);
        roots.extend(group);
        folded
    }
}

fn correct_umis<S: BuildHasher + Sync>(
    evidence: &StdHashMap<RawKey, Evidence, S>,
    sorted_keys: &[RawKey],
    total_groups: usize,
    workers: usize,
    skip: bool,
    rules: UmiRules,
) -> (HashMap<RawKey, RawKey>, u64) {
    if skip {
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
            "[morselle molecules] phase=umi_correct skipped (--tagseq without --position-aware-umis) raw_nodes={}",
            evidence.len()
        );
        return (evidence.keys().map(|&key| (key, key)).collect(), 0);
    }
    if workers <= 1 {
        eprintln!("[morselle molecules] phase=umi_correct mode=serial start groups={total_groups}");
        let mut roots = HashMap::with_capacity_and_hasher(sorted_keys.len(), Default::default());
        let mut folded = 0;
        for run in sorted_keys.chunk_by(same_umi_group) {
            folded += correct_run(run, evidence, &mut roots, rules);
        }
        return (roots, folded);
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .expect("creating UMI-correction worker pool");
    info!(
        workers,
        groups = total_groups,
        "correcting UMI groups in parallel"
    );
    eprintln!(
        "[morselle molecules] phase=umi_correct start workers={workers} groups={total_groups}"
    );
    let completed = Arc::new(AtomicUsize::new(0));
    let checkpoint = (total_groups / 100).max(1_000);
    pool.install(|| {
        sorted_keys
            .par_chunk_by(same_umi_group)
            .fold(
                || (HashMap::default(), 0_u64),
                |(mut roots, mut folded), run| {
                    folded += correct_run(run, evidence, &mut roots, rules);
                    let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                    if done % checkpoint == 0 || done == total_groups {
                        eprintln!(
                            "[morselle molecules] phase=umi_correct groups={done}/{total_groups}"
                        );
                    }
                    (roots, folded)
                },
            )
            .reduce(
                || (HashMap::default(), 0),
                |(a, fa), (b, fb)| (merge_roots(a, b), fa + fb),
            )
    })
}

/// Sum each corrected root's read support and pick its representative read.
fn aggregate_molecules<S: BuildHasher>(
    evidence: &StdHashMap<RawKey, Evidence, S>,
    roots: &HashMap<RawKey, RawKey>,
) -> HashMap<RawKey, Molecule> {
    let mut molecules = HashMap::<RawKey, Molecule>::default();
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
                resolved: Resolved::UNRESOLVED,
            });
    }
    molecules
}

/// `evidence` may use any hasher (callers and tests pass std `HashMap`s);
/// every internal table is FxHash.
pub fn build_decisions<S: BuildHasher + Sync>(
    evidence: &StdHashMap<RawKey, Evidence, S>,
    metrics: &mut Metrics,
    saturation_interval: Option<u64>,
    workers: usize,
    tagseq: bool,
    fractional: bool,
) -> (HashMap<RawKey, Decision>, Option<Saturation>) {
    let policy = DedupPolicy {
        tagseq,
        fractional,
        correct_umis: !tagseq,
        umi_rules: UmiRules::default(),
    };
    build_decisions_with(evidence, metrics, saturation_interval, workers, policy)
}

/// How molecules are formed and collisions resolved.
#[derive(Clone, Copy, Debug)]
pub struct DedupPolicy {
    /// Count every gene in a collision (`Shared`) unless `fractional`.
    pub tagseq: bool,
    /// Split collisions by read support (`Fractional`).
    pub fractional: bool,
    /// Directional 1-mismatch UMI correction within (cell, sample, gene).
    pub correct_umis: bool,
    /// Opt-in SBX rules applied with the correction (no effect when
    /// `correct_umis` is false). See `umi_rules`.
    pub umi_rules: UmiRules,
}

impl DedupPolicy {
    /// The policy `morselle-molecules` uses for its flags.
    pub fn from_flags(tagseq: bool, fractional: bool, position_aware: bool) -> Self {
        Self {
            tagseq,
            fractional,
            correct_umis: !tagseq || position_aware,
            umi_rules: UmiRules::default(),
        }
    }

    /// The same policy with the opt-in SBX UMI rules.
    pub fn with_umi_rules(self, umi_rules: UmiRules) -> Self {
        Self { umi_rules, ..self }
    }
}

pub fn build_decisions_with<S: BuildHasher + Sync>(
    evidence: &StdHashMap<RawKey, Evidence, S>,
    metrics: &mut Metrics,
    saturation_interval: Option<u64>,
    workers: usize,
    policy: DedupPolicy,
) -> (HashMap<RawKey, Decision>, Option<Saturation>) {
    let DedupPolicy {
        tagseq,
        fractional,
        correct_umis: correct,
        umi_rules,
    } = policy;
    // Grouped by (cell, sample, gene), not just (cell, gene): the same CB
    // whitelist is shared across every sample in a pooled Flex run, so
    // omitting `sample` here would let correct_group()'s internal by-UMI
    // map (keyed on UMI alone within a group) silently overwrite one
    // sample's evidence with another's whenever two samples' reads for the
    // same cell+gene happen to share a UMI value.
    // One sorted key vector instead of a HashMap<(cell, sample, gene),
    // Vec<RawKey>>: no per-group heap allocation (most groups hold a single
    // UMI), and groups are contiguous runs. Group order never matters --
    // correct_group sorts its own input and all tie-breaks are explicit.
    let mut sorted_keys: Vec<RawKey> = evidence.keys().copied().collect();
    sorted_keys.par_sort_unstable_by_key(|k| (k.cell, k.sample, k.gene, k.umi));
    let total_groups = sorted_keys.chunk_by(same_umi_group).count();
    eprintln!(
        "[morselle molecules] phase=umi_group_index complete groups={} raw_nodes={}",
        total_groups,
        evidence.len()
    );
    let mut saturation =
        saturation_interval.map(|interval| Saturation::new(metrics.eligible_records, interval));
    if let Some(curve) = saturation.as_mut() {
        for run in sorted_keys.chunk_by(same_umi_group) {
            let min_hash = run
                .iter()
                .map(|key| evidence[key].min_read_hash)
                .min()
                .expect("non-empty UMI group");
            curve.add_cell_gene(run[0].cell, run[0].gene, min_hash);
        }
    }

    let (roots, folded) = correct_umis(
        evidence,
        &sorted_keys,
        total_groups,
        workers,
        !correct,
        umi_rules,
    );
    metrics.umi_singletons_folded += folded;
    drop(sorted_keys);
    eprintln!(
        "[morselle molecules] phase=umi_correct complete roots={}",
        roots.len()
    );
    // `+=`, not `=`: the sharded path calls build_decisions() once per shard
    // against the same accumulating Metrics.
    metrics.raw_cell_gene_umis += evidence.len() as u64;
    for key in evidence.keys() {
        super::barcode_metrics::bump(metrics, key.cell, key.sample, 1);
    }
    metrics.umi_nodes_collapsed += roots
        .iter()
        .filter(|(key, root)| key.umi != root.umi)
        .count() as u64;

    eprintln!("[morselle molecules] phase=molecule_aggregate start");
    let mut molecules = aggregate_molecules(evidence, &roots);
    metrics.corrected_umi_nodes += molecules.len() as u64;
    metrics.molecules_before_gene_resolution += molecules.len() as u64;
    eprintln!(
        "[morselle molecules] phase=molecule_aggregate complete corrected_nodes={}",
        molecules.len()
    );

    super::collisions::resolve_collisions(&mut molecules, metrics, tagseq, fractional);

    if let Some(curve) = saturation.as_mut() {
        // Fractional claimants form ONE resolved molecule. Its first
        // supporting read is the earliest event across all its genes.
        let mut fractional_events = HashMap::<_, u32>::default();
        for (&root, molecule) in &molecules {
            if molecule.resolved.outcome == super::types::Outcome::Fractional {
                fractional_events
                    .entry((root.cell, root.sample, root.umi))
                    .and_modify(|hash| *hash = (*hash).min(molecule.min_read_hash))
                    .or_insert(molecule.min_read_hash);
            } else if molecule.resolved.outcome.is_counted() {
                curve.add_molecule(root.cell, molecule.min_read_hash);
            }
        }
        for ((cell, _, _), hash) in fractional_events {
            curve.add_molecule(cell, hash);
        }
    }
    eprintln!("[morselle molecules] phase=decision_materialize start");
    let decisions: HashMap<RawKey, Decision> = evidence
        .keys()
        .map(|&key| {
            let root = roots[&key];
            let molecule = &molecules[&root];
            let resolved = molecule.resolved;
            (
                key,
                Decision {
                    corrected: root.umi,
                    representative: molecule.representative,
                    collision_degree: resolved.degree,
                    outcome: resolved.outcome,
                    weight: resolved.weight,
                },
            )
        })
        .collect();
    eprintln!(
        "[morselle molecules] phase=decision_materialize complete decisions={}",
        decisions.len()
    );
    (decisions, saturation)
}
