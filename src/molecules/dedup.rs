use std::collections::HashMap;
use std::path::Path;
use anyhow::Result;
use rayon::prelude::*;

use super::saturation::Saturation;
use super::types::{
    Decision, Evidence, Metrics, Outcome, RawKey, Dna,
};


pub fn correct_group(
    mut keys: Vec<RawKey>,
    evidence: &HashMap<RawKey, Evidence>,
) -> HashMap<RawKey, RawKey> {
    keys.sort_by_key(|key| (std::cmp::Reverse(evidence[key].reads), key.umi));
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

pub fn build_decisions(
    evidence: &HashMap<RawKey, Evidence>,
    metrics: &mut Metrics,
    saturation_interval: Option<u64>,
    workers: usize,
    tagseq: bool,
    fractional_collisions: bool,
) -> (HashMap<RawKey, Decision>, Option<Saturation>) {
    let mut groups = HashMap::<(Dna, u32, u32), Vec<RawKey>>::new();
    for &key in evidence.keys() {
        groups
            .entry((key.cell, key.sample, key.gene))
            .or_default()
            .push(key);
    }
    metrics.raw_cell_gene_umis += evidence.len() as u64;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers.max(1))
        .build()
        .expect("building rayon threadpool");
    let group_vec: Vec<Vec<RawKey>> = groups.into_values().collect();
    let roots = pool.install(|| {
        group_vec
            .into_par_iter()
            .map(|group| correct_group(group, evidence))
            .reduce_with(merge_roots)
            .unwrap_or_default()
    });

    let mut corrected_nodes = HashMap::<RawKey, (RawKey, Evidence)>::new();
    for (&child, &root) in &roots {
        let child_ev = evidence[&child];
        if child.umi != root.umi {
            metrics.umi_nodes_collapsed += 1;
        }
        corrected_nodes
            .entry(root)
            .and_modify(|(winner_key, winner_ev)| {
                winner_ev.reads += child_ev.reads;
                if (child_ev.best_mapq, child_ev.best_aligned, std::cmp::Reverse(child_ev.best_ordinal))
                    > (winner_ev.best_mapq, winner_ev.best_aligned, std::cmp::Reverse(winner_ev.best_ordinal))
                {
                    *winner_key = child;
                    winner_ev.best_ordinal = child_ev.best_ordinal;
                    winner_ev.best_mapq = child_ev.best_mapq;
                    winner_ev.best_aligned = child_ev.best_aligned;
                }
            })
            .or_insert((child, child_ev));
    }
    metrics.corrected_umi_nodes += corrected_nodes.len() as u64;

    let mut collision_groups = HashMap::<(Dna, u32, Dna), Vec<RawKey>>::new();
    for &root in corrected_nodes.keys() {
        collision_groups
            .entry((root.cell, root.sample, root.umi))
            .or_default()
            .push(root);
    }
    metrics.molecules_before_gene_resolution += collision_groups.len() as u64;

    let mut saturation =
        saturation_interval.map(|interval| Saturation::new(metrics.eligible_records, interval));

    let use_fractional = fractional_collisions;
    let mut decisions = HashMap::<RawKey, Decision>::new();
    for group in collision_groups.values() {
        let degree = group.len() as u32;
        if degree > 1 {
            metrics.collision_groups += 1;
        }
        if use_fractional {
            let total_reads: u64 = group
                .iter()
                .map(|&root| corrected_nodes[&root].1.reads as u64)
                .sum();
            for &root in group {
                let (_, winner_ev) = corrected_nodes[&root];
                let root_reads = winner_ev.reads as u64;
                let weight = if total_reads > 0 {
                    root_reads as f32 / total_reads as f32
                } else {
                    1.0 / degree as f32
                };
                let outcome = if degree > 1 {
                    metrics.collision_molecules_fractional += 1;
                    Outcome::Fractional
                } else {
                    Outcome::Unique
                };
                for (&child, &child_root) in &roots {
                    if child_root == root {
                        decisions.insert(
                            child,
                            Decision {
                                corrected: root.umi,
                                representative: winner_ev.best_ordinal,
                                collision_degree: degree,
                                outcome,
                                weight,
                            },
                        );
                    }
                }
                metrics.molecules_counted += 1;
            }
        } else {
            let mut gene_support = HashMap::<u32, (u64, RawKey, Evidence)>::new();
            for &root in group {
                let (winner_key, winner_ev) = corrected_nodes[&root];
                gene_support
                    .entry(root.gene)
                    .and_modify(|(sum, win_k, win_ev)| {
                        *sum += winner_ev.reads as u64;
                        if (winner_ev.best_mapq, winner_ev.best_aligned, std::cmp::Reverse(winner_ev.best_ordinal))
                            > (win_ev.best_mapq, win_ev.best_aligned, std::cmp::Reverse(win_ev.best_ordinal))
                        {
                            *win_k = winner_key;
                            *win_ev = winner_ev;
                        }
                    })
                    .or_insert((winner_ev.reads as u64, winner_key, winner_ev));
            }

            let max_support = gene_support.values().map(|(sum, ..)| *sum).max().unwrap_or(0);
            let winners_count = gene_support.values().filter(|(sum, ..)| *sum == max_support).count();

            for &root in group {
                let (support, _, winner_ev) = gene_support[&root.gene];
                let outcome = if winners_count > 1 && support == max_support {
                    metrics.collision_ties += 1;
                    Outcome::Tie
                } else if support == max_support {
                    if degree > 1 {
                        metrics.collision_molecules_won += 1;
                        Outcome::Won
                    } else {
                        Outcome::Unique
                    }
                } else {
                    if degree > 1 {
                        metrics.collision_molecules_lost += 1;
                    }
                    Outcome::Lost
                };

                let weight = if matches!(outcome, Outcome::Won | Outcome::Unique) {
                    1.0
                } else {
                    0.0
                };

                if matches!(outcome, Outcome::Won | Outcome::Unique) {
                    metrics.molecules_counted += 1;
                }

                for (&child, &child_root) in &roots {
                    if child_root == root {
                        decisions.insert(
                            child,
                            Decision {
                                corrected: root.umi,
                                representative: winner_ev.best_ordinal,
                                collision_degree: degree,
                                outcome,
                                weight,
                            },
                        );
                    }
                }
            }
        }
    }

    if let Some(sat) = &mut saturation {
        for (&key, ev) in evidence {
            sat.add_cell_gene(key.cell, key.gene, ev.min_read_hash);
            sat.add_molecule(key.cell, ev.min_read_hash);
        }
    }

    (decisions, saturation)
}

pub fn deduplicate_evidence(
    evidence: &HashMap<RawKey, Evidence>,
    metrics: &mut Metrics,
    saturation_interval: Option<u64>,
    workers: usize,
    tagseq: bool,
    fractional_collisions: bool,
) -> (HashMap<RawKey, Decision>, Option<Saturation>) {
    build_decisions(
        evidence,
        metrics,
        saturation_interval,
        workers,
        tagseq,
        fractional_collisions,
    )
}

pub fn process_shard(
    path: &Path,
    threads: usize,
    require_pb: bool,
    saturation_interval: Option<u64>,
    workers: usize,
    tagseq: bool,
    fractional_collisions: bool,
) -> Result<(
    Metrics,
    HashMap<RawKey, Decision>,
    Vec<String>,
    Vec<String>,
    Option<Saturation>,
)> {
    let mut metrics = Metrics::default();
    let (evidence, genes, _gene_ids, samples, _sample_ids) =
        super::evidence::collect_evidence(path, threads, require_pb, &mut metrics)?;
    let (decisions, saturation) = deduplicate_evidence(
        &evidence,
        &mut metrics,
        saturation_interval,
        workers,
        tagseq,
        fractional_collisions,
    );
    Ok((metrics, decisions, genes, samples, saturation))
}

pub fn process_and_write_shard(
    input_path: &Path,
    output_path: &Path,
    threads: usize,
    require_pb: bool,
    saturation_interval: Option<u64>,
    workers: usize,
    tagseq: bool,
    fractional_collisions: bool,
) -> Result<(Metrics, Option<Saturation>)> {
    let mut metrics = Metrics::default();
    let (evidence, _genes, gene_ids, _samples, sample_ids) =
        super::evidence::collect_evidence(input_path, threads, require_pb, &mut metrics)?;
    let (decisions, saturation) = deduplicate_evidence(
        &evidence,
        &mut metrics,
        saturation_interval,
        workers,
        tagseq,
        fractional_collisions,
    );

    // Drop evidence map to free memory before writing output
    drop(evidence);

    super::pipeline::write_shard_output(
        input_path,
        output_path,
        threads,
        require_pb,
        &decisions,
        &gene_ids,
        &sample_ids,
        &mut metrics,
    )?;

    // `decisions` dropped at end of scope
    Ok((metrics, saturation))
}
