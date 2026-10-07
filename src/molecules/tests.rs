use std::collections::HashMap;
use std::fs::File;
use rust_htslib::bam::{self, record::Aux};


use super::dedup::build_decisions;
use super::evidence::collect_evidence;
use super::helpers::intern_sample;
use super::saturation::Saturation;
use super::types::{
    Dna, Evidence, Metrics, Outcome, RawKey, NO_SAMPLE,
};

fn evidence(reads: u32, ordinal: u64) -> Evidence {
    Evidence {
        reads,
        best_ordinal: ordinal,
        best_mapq: 60,
        best_aligned: 100,
        min_read_hash: ordinal as u32,
    }
}

#[test]
fn dna_round_trip_and_hamming() {
    let a = Dna::pack("AACCGGTT").unwrap();
    let b = Dna::pack("AACCGGTA").unwrap();
    assert_eq!(a.unpack(), "AACCGGTT");
    assert_eq!(a.hamming(b), Some(1));
    assert_eq!(Dna::pack("AACCGGTT-1"), Some(a));
    assert_eq!(std::mem::size_of::<Evidence>(), 24);
}

#[test]
fn one_mismatch_neighbors_match_brute_force_hamming() {
    use rand::{Rng, SeedableRng};
    let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(42);
    let bases = [b'A', b'C', b'G', b'T'];
    for _ in 0..200 {
        let seq: String = (0..12)
            .map(|_| bases[rng.gen_range(0..4)] as char)
            .collect();
        let target = Dna::pack(&seq).unwrap();

        let neighbors: std::collections::HashSet<Dna> =
            target.one_mismatch_neighbors().collect();
        assert_eq!(neighbors.len(), 36, "3*len distinct neighbors for {seq}");
        assert!(!neighbors.contains(&target));
        for &n in &neighbors {
            assert_eq!(target.hamming(n), Some(1));
        }

        for _ in 0..50 {
            let candidate_seq: String = (0..12)
                .map(|_| bases[rng.gen_range(0..4)] as char)
                .collect();
            let candidate = Dna::pack(&candidate_seq).unwrap();
            assert_eq!(
                neighbors.contains(&candidate),
                target.hamming(candidate) == Some(1)
            );
        }
    }
}

#[test]
fn corrects_within_gene_then_resolves_cross_gene_collision() {
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
    let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, false, false);
    assert_eq!(decisions[&gene0_error].corrected, major);
    assert_eq!(decisions[&gene0_major].outcome, Outcome::Won);
    assert_eq!(decisions[&gene0_error].outcome, Outcome::Won);
    assert_eq!(decisions[&gene1_major].outcome, Outcome::Lost);
    assert_eq!(metrics.umi_nodes_collapsed, 1);
    assert_eq!(metrics.collision_groups, 1);
    assert_eq!(metrics.molecules_counted, 1);

    let mut parallel_metrics = Metrics::default();
    let (parallel, _) = build_decisions(&input, &mut parallel_metrics, None, 2, false, false);
    assert_eq!(parallel[&gene0_error].corrected, major);
    assert_eq!(parallel[&gene0_major].outcome, Outcome::Won);
    assert_eq!(parallel[&gene1_major].outcome, Outcome::Lost);
    assert_eq!(
        parallel_metrics.molecules_counted,
        metrics.molecules_counted
    );
}

#[test]
fn different_samples_sharing_cell_gene_umi_are_never_merged() {
    let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
    let umi = Dna::pack("AAAAAAAAAAAA").unwrap();
    let mut input = HashMap::new();
    let sample_a = RawKey {
        cell,
        sample: 0,
        gene: 0,
        umi,
    };
    let sample_b = RawKey {
        cell,
        sample: 1,
        gene: 0,
        umi,
    };
    input.insert(sample_a, evidence(10, 0));
    input.insert(sample_b, evidence(7, 1));

    let mut metrics = Metrics::default();
    let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, false, false);
    assert_eq!(decisions[&sample_a].outcome, Outcome::Unique);
    assert_eq!(decisions[&sample_b].outcome, Outcome::Unique);
    assert_eq!(metrics.molecules_counted, 2);

    let mut parallel_metrics = Metrics::default();
    let (parallel, _) = build_decisions(&input, &mut parallel_metrics, None, 2, false, false);
    assert_eq!(parallel[&sample_a].outcome, Outcome::Unique);
    assert_eq!(parallel[&sample_b].outcome, Outcome::Unique);
    assert_eq!(parallel_metrics.molecules_counted, 2);
}

#[test]
fn intern_sample_assigns_stable_ids_per_distinct_pb() {
    let mut samples = Vec::new();
    let mut sample_ids = HashMap::new();
    let a = intern_sample("D-A01", &mut samples, &mut sample_ids);
    let b = intern_sample("D-A02", &mut samples, &mut sample_ids);
    let a_again = intern_sample("D-A01", &mut samples, &mut sample_ids);
    assert_ne!(a, b);
    assert_eq!(a, a_again);
    assert_eq!(samples.len(), 2);
}

#[test]
fn records_without_a_confident_sample_barcode_are_dropped_not_bucketed() {
    use rust_htslib::bam::{header::HeaderRecord, record::Cigar, record::CigarString};

    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("input.bam");
    let mut header = bam::Header::new();
    header.push_record(
        HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "ENSG1")
            .push_tag(b"LN", 1_000_000),
    );
    let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    let cigar = CigarString(vec![Cigar::Match(20)]);
    for (name, cb, umi, pb) in [
        ("no_pb_tag", "AAAACCCCGGGGTTTT", "AAAAAAAAAAAA", None),
        ("empty_pb", "AAAACCCCGGGGTTTT", "CCCCCCCCCCCC", Some("")),
        ("real_pb", "AAAACCCCGGGGTTTT", "GGGGGGGGGGGG", Some("D-A01")),
    ] {
        let mut record = bam::Record::new();
        record.set(
            name.as_bytes(),
            Some(&cigar),
            b"ACGTACGTACGTACGTACGT",
            &[30; 20],
        );
        record.set_tid(0);
        record.set_pos(1);
        record.set_flags(0);
        record.push_aux(b"CB", Aux::String(cb)).unwrap();
        record.push_aux(b"UR", Aux::String(umi)).unwrap();
        record.push_aux(b"GX", Aux::String("ENSG1")).unwrap();
        if let Some(pb) = pb {
            record.push_aux(b"PB", Aux::String(pb)).unwrap();
        }
        writer.write(&record).unwrap();
    }
    drop(writer);

    let mut metrics = Metrics::default();
    let (evidence, ..) = collect_evidence(&path, 1, true, &mut metrics).unwrap();
    assert_eq!(
        metrics.eligible_records, 1,
        "only the real-PB record is eligible"
    );
    assert_eq!(
        metrics.no_sample_records, 2,
        "both no-PB and empty-PB records are counted as dropped"
    );
    assert_eq!(evidence.len(), 1);

    let mut tagseq_metrics = Metrics::default();
    let (tagseq_evidence, ..) = collect_evidence(&path, 1, false, &mut tagseq_metrics).unwrap();
    assert_eq!(
        tagseq_metrics.eligible_records, 3,
        "tag-seq must not require PB"
    );
    assert_eq!(tagseq_metrics.no_sample_records, 0);
    assert_eq!(tagseq_evidence.len(), 3);
}

#[test]
fn saturation_is_nested_and_reports_run_and_cells() {
    let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
    let mut curve = Saturation::new(10, 5);
    curve.add_cell_gene(cell, 7, 0);
    curve.add_molecule(cell, 0);
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("curve.tsv.gz");
    let summary = curve.write(&path).unwrap();
    assert_eq!(summary.points, 2);
    assert_eq!(summary.observed_cells, 1);
    let contents = {
        use flate2::read::MultiGzDecoder;
        use std::io::Read;
        let mut decoder = MultiGzDecoder::new(File::open(path).unwrap());
        let mut text = String::new();
        decoder.read_to_string(&mut text).unwrap();
        text
    };
    assert!(contents.contains("run\t.\t5\t0.5\t1\t1\t1\t1"));
    assert!(contents.contains("run\t.\t10\t1\t1\t1\t0\t0"));
    assert_eq!(contents.lines().count(), 5);
}

#[test]
fn fractional_umi_collision_resolution_distributes_weights_proportionally() {
    let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
    let umi = Dna::pack("AAAAAAAAAAAA").unwrap();
    let mut input = HashMap::new();
    let g0 = RawKey {
        cell,
        sample: 0,
        gene: 0,
        umi,
    };
    let g1 = RawKey {
        cell,
        sample: 0,
        gene: 1,
        umi,
    };
    let g2 = RawKey {
        cell,
        sample: 0,
        gene: 2,
        umi,
    };
    input.insert(g0, evidence(5, 0));
    input.insert(g1, evidence(3, 1));
    input.insert(g2, evidence(2, 2));

    let mut metrics = Metrics::default();
    let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, false, true);
    assert_eq!(decisions.len(), 3);
    assert_eq!(decisions[&g0].outcome, Outcome::Fractional);
    assert_eq!(decisions[&g1].outcome, Outcome::Fractional);
    assert_eq!(decisions[&g2].outcome, Outcome::Fractional);
    assert!((decisions[&g0].weight - 0.5).abs() < 1e-5);
    assert!((decisions[&g1].weight - 0.3).abs() < 1e-5);
    assert!((decisions[&g2].weight - 0.2).abs() < 1e-5);
    assert_eq!(metrics.collision_groups, 1);
    assert_eq!(metrics.collision_molecules_fractional, 3);
}

