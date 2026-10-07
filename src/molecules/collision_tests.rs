//! Collision-policy tests: winner-take-all, --tagseq shared, and fractional.
//! The first two were dropped by the 0.5.7 refactor and are restored from
//! commit 4d42a18; the fractional ones pin the policy added on top.
use std::collections::HashMap;

use super::dedup::build_decisions;
use super::types::{Dna, Evidence, Metrics, Outcome, RawKey, NO_SAMPLE};

fn evidence(reads: u32, ordinal: u64) -> Evidence {
    Evidence {
        reads,
        best_ordinal: ordinal,
        best_mapq: 60,
        best_aligned: 100,
        min_read_hash: ordinal as u32,
    }
}

fn key(gene: u32, umi: &str) -> RawKey {
    RawKey {
        cell: Dna::pack("AACCGGTTAACCGGTT").unwrap(),
        sample: NO_SAMPLE,
        gene,
        umi: Dna::pack(umi).unwrap(),
    }
}

#[test]
fn tagseq_without_fractional_counts_every_gene_in_a_collision() {
    // Same raw evidence as corrects_within_gene_then_resolves_cross_gene_collision,
    // but under --tagseq: UMI correction is also skipped (pooled-sample UMI
    // diversity makes the 1-mismatch "sequencing error" assumption unsafe
    // the same way winner-take-all collision resolution is unsafe), so
    // gene0's "error" UMI no longer merges into "major" -- it stands as
    // its own uncollided molecule. gene0_major and gene1_major still
    // collide on the shared "major" UMI, and --tagseq keeps both instead
    // of picking a winner.
    let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
    let major = Dna::pack("AAAAAAAAAAAA").unwrap();
    let error = Dna::pack("AAAAAAAAAAAC").unwrap();
    let mut input = HashMap::new();
    let gene0_major = RawKey {
        cell,
        sample: NO_SAMPLE,
        gene: 0,
        umi: major,
    };
    let gene0_error = RawKey {
        cell,
        sample: NO_SAMPLE,
        gene: 0,
        umi: error,
    };
    let gene1_major = RawKey {
        cell,
        sample: NO_SAMPLE,
        gene: 1,
        umi: major,
    };
    input.insert(gene0_major, evidence(10, 0));
    input.insert(gene0_error, evidence(2, 1));
    input.insert(gene1_major, evidence(3, 2));

    let mut metrics = Metrics::default();
    let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, true, false);
    // No UMI correction under --tagseq: "error" stays uncorrected and
    // uncollided (no other gene/UMI shares it).
    assert_eq!(decisions[&gene0_error].corrected, error);
    assert_eq!(decisions[&gene0_error].outcome, Outcome::Unique);
    // "major" is shared by two genes -- both kept, neither discarded.
    assert_eq!(decisions[&gene0_major].outcome, Outcome::Shared);
    assert_eq!(decisions[&gene1_major].outcome, Outcome::Shared);
    assert_eq!(decisions[&gene0_major].collision_degree, 2);
    assert_eq!(decisions[&gene1_major].collision_degree, 2);
    assert_eq!(metrics.umi_nodes_collapsed, 0);
    // Diagnostics still reflect that a collision happened...
    assert_eq!(metrics.collision_groups, 1);
    // ...but nothing is discarded: all three raw nodes are counted.
    assert_eq!(metrics.collision_molecules_shared, 2);
    assert_eq!(metrics.collision_molecules_won, 0);
    assert_eq!(metrics.collision_molecules_lost, 0);
    assert_eq!(metrics.molecules_counted, 3);

    let mut parallel_metrics = Metrics::default();
    let (parallel, _) = build_decisions(&input, &mut parallel_metrics, None, 2, true, false);
    assert_eq!(parallel[&gene0_error].corrected, error);
    assert_eq!(parallel[&gene0_major].outcome, Outcome::Shared);
    assert_eq!(parallel[&gene1_major].outcome, Outcome::Shared);
    assert_eq!(
        parallel_metrics.molecules_counted,
        metrics.molecules_counted
    );
}

#[test]
fn equal_cross_gene_support_is_not_counted() {
    let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
    let umi = Dna::pack("ACGTACGTACGT").unwrap();
    let a = RawKey {
        cell,
        sample: NO_SAMPLE,
        gene: 0,
        umi,
    };
    let b = RawKey {
        cell,
        sample: NO_SAMPLE,
        gene: 1,
        umi,
    };
    let input = HashMap::from([(a, evidence(2, 0)), (b, evidence(2, 1))]);
    let mut metrics = Metrics::default();
    let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, false, false);
    assert_eq!(decisions[&a].outcome, Outcome::Tie);
    assert_eq!(decisions[&b].outcome, Outcome::Tie);
    assert_eq!(metrics.molecules_counted, 0);
    assert_eq!(metrics.collision_ties, 1);
}

/// --tagseq defaults to fractional: still no UMI correction, and a
/// collision is one molecule split by read support, not two full molecules.
#[test]
fn tagseq_fractional_skips_correction_and_splits_one_molecule() {
    let g0_major = key(0, "AAAAAAAAAAAA");
    let g0_error = key(0, "AAAAAAAAAAAC");
    let g1_major = key(1, "AAAAAAAAAAAA");
    let input = HashMap::from([
        (g0_major, evidence(9, 0)),
        (g0_error, evidence(2, 1)),
        (g1_major, evidence(3, 2)),
    ]);
    let mut metrics = Metrics::default();
    let (d, _) = build_decisions(&input, &mut metrics, None, 1, true, true);
    assert_eq!(d[&g0_error].corrected, g0_error.umi);
    assert_eq!(d[&g0_error].outcome, Outcome::Unique);
    assert_eq!(d[&g0_error].weight, 1.0);
    assert_eq!(d[&g0_major].outcome, Outcome::Fractional);
    assert_eq!(d[&g1_major].outcome, Outcome::Fractional);
    assert!((d[&g0_major].weight - 0.75).abs() < 1e-6);
    assert!((d[&g1_major].weight - 0.25).abs() < 1e-6);
    assert_eq!(metrics.umi_nodes_collapsed, 0);
    assert_eq!(metrics.collision_groups, 1);
    assert_eq!(metrics.collision_molecules_fractional, 2);
    // One unique molecule + one collision group (weights sum to 1).
    assert_eq!(metrics.molecules_counted, 2);
}

/// Without --tagseq, fractional mode still corrects UMIs first and only
/// changes how the surviving collision is resolved; an exact tie splits
/// 50/50 instead of being discarded.
#[test]
fn fractional_without_tagseq_corrects_first_and_splits_ties() {
    let g0_major = key(0, "AAAAAAAAAAAA");
    let g0_error = key(0, "AAAAAAAAAAAC");
    let g1_major = key(1, "AAAAAAAAAAAA");
    let input = HashMap::from([
        (g0_major, evidence(4, 0)),
        (g0_error, evidence(1, 1)),
        (g1_major, evidence(5, 2)),
    ]);
    let mut metrics = Metrics::default();
    let (d, _) = build_decisions(&input, &mut metrics, None, 1, false, true);
    assert_eq!(d[&g0_error].corrected, g0_major.umi);
    assert_eq!(metrics.umi_nodes_collapsed, 1);
    // gene0 support 4+1 = 5 ties gene1's 5.
    for k in [g0_major, g0_error, g1_major] {
        assert_eq!(d[&k].outcome, Outcome::Fractional);
        assert!((d[&k].weight - 0.5).abs() < 1e-6);
    }
    assert_eq!(metrics.collision_ties, 0);
    assert_eq!(metrics.molecules_counted, 1);
}

#[test]
fn winner_take_all_weights_are_one_or_zero() {
    let a = key(0, "ACGTACGTACGT");
    let b = key(1, "ACGTACGTACGT");
    let input = HashMap::from([(a, evidence(3, 0)), (b, evidence(1, 1))]);
    let mut metrics = Metrics::default();
    let (d, _) = build_decisions(&input, &mut metrics, None, 1, false, false);
    assert_eq!((d[&a].outcome, d[&a].weight), (Outcome::Won, 1.0));
    assert_eq!((d[&b].outcome, d[&b].weight), (Outcome::Lost, 0.0));
}
