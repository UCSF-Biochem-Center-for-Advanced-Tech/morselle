#[cfg(test)]
mod tests {
    use rust_htslib::bam::{self, record::Aux, Read};
    use std::collections::HashMap;
    use std::path::{Path, PathBuf};

    use crate::molecules::dedup::build_decisions;
    use crate::molecules::evidence::{
        collate, collect_evidence, collect_evidence_indexed, collect_region, indexed_regions,
    };
    use crate::molecules::helpers::shard_hash;
    use crate::molecules::types::{Dna, Evidence, Metrics, RawKey};

    fn write_test_bam(path: &Path, records: &[(&str, i64, &str, &str, &str)]) {
        use rust_htslib::bam::{header::HeaderRecord, record::Cigar, record::CigarString};
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "ENSG1")
                .push_tag(b"LN", 1_000_000),
        );
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
        let cigar = CigarString(vec![Cigar::Match(20)]);
        for (name, pos, cb, umi, gene) in records {
            let mut record = bam::Record::new();
            record.set(
                name.as_bytes(),
                Some(&cigar),
                b"ACGTACGTACGTACGTACGT",
                &[30; 20],
            );
            record.set_tid(0);
            record.set_pos(*pos);
            record.set_flags(0);
            record.push_aux(b"CB", Aux::String(cb)).unwrap();
            record.push_aux(b"UR", Aux::String(umi)).unwrap();
            record.push_aux(b"GX", Aux::String(gene)).unwrap();
            writer.write(&record).unwrap();
        }
    }

    fn read_tag_set(path: &Path) -> std::collections::BTreeSet<(String, String, String)> {
        let mut reader = bam::Reader::from_path(path).unwrap();
        reader
            .records()
            .map(|r| {
                let record = r.unwrap();
                (
                    String::from_utf8(record.qname().to_vec()).unwrap(),
                    crate::molecules::helpers::aux_string(&record, b"CB").unwrap(),
                    crate::molecules::helpers::aux_string(&record, b"UR").unwrap(),
                )
            })
            .collect()
    }

    #[test]
    fn indexed_regions_collect_each_alignment_once() {
        use rust_htslib::bam::{header::HeaderRecord, record::Cigar, record::CigarString};

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("input.bam");
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 100),
        );
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        let cigar = CigarString(vec![Cigar::Match(20)]);
        for (name, pos, umi) in [
            ("a", 1, "AAAAAAAAAAAA"),
            ("b", 9, "AAAAAAAAAAAA"),
            ("c", 10, "CCCCCCCCCCCC"),
        ] {
            let mut record = bam::Record::new();
            record.set(
                name.as_bytes(),
                Some(&cigar),
                b"ACGTACGTACGTACGTACGT",
                &[30; 20],
            );
            record.set_tid(0);
            record.set_pos(pos);
            record.set_flags(0);
            record
                .push_aux(b"CB", Aux::String("AACCGGTTAACCGGTT"))
                .unwrap();
            record.push_aux(b"UR", Aux::String(umi)).unwrap();
            record.push_aux(b"GX", Aux::String("ENSG1")).unwrap();
            record.push_aux(b"PB", Aux::String("D-A01")).unwrap();
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        assert_eq!(indexed_regions(&path, 10).unwrap().len(), 10);
        assert_eq!(
            collect_region(&path, 0, 0, 10, false)
                .unwrap()
                .eligible_records,
            2
        );

        let mut serial_metrics = Metrics::default();
        let (serial, _, _, _, _) = collect_evidence(&path, 1, false, &mut serial_metrics).unwrap();
        let mut parallel_metrics = Metrics::default();
        let (parallel, _, _, _, _) =
            collect_evidence_indexed(&path, 2, 10, false, &mut parallel_metrics).unwrap();
        assert_eq!(serial_metrics.eligible_records, 3);
        assert_eq!(parallel_metrics.eligible_records, 3);
        assert_eq!(serial.len(), 2);
        assert_eq!(parallel.len(), 2);
        assert_eq!(
            serial.values().map(|item| item.reads as u64).sum::<u64>(),
            3
        );
        assert_eq!(
            parallel.values().map(|item| item.reads as u64).sum::<u64>(),
            3
        );
        let mut serial_decision_metrics = Metrics::default();
        let mut parallel_decision_metrics = Metrics::default();
        let _ = build_decisions(&serial, &mut serial_decision_metrics, None, 2, false, false);
        let _ = build_decisions(
            &parallel,
            &mut parallel_decision_metrics,
            None,
            2,
            false,
            false,
        );
        assert_eq!(
            serial_decision_metrics.molecules_counted,
            parallel_decision_metrics.molecules_counted
        );
    }

    /// Core correctness guarantee behind `--shards`: every real grouping
    /// (`(cell, sample, gene)`) and collision-resolution (`(cell, sample,
    /// umi)`) key in `build_decisions` is scoped to `(cell, sample)`, so
    /// partitioning evidence by `(cell, sample)` into independent shards and
    /// merging each shard's own `build_decisions()` output must produce
    /// exactly the decisions (and exactly the accumulated metrics) that a
    /// single global `build_decisions()` call over the same evidence would
    /// -- this is what makes it safe to process shards one at a time,
    /// discarding memory in between, instead of materializing one global
    /// evidence/roots/decisions table.
    #[test]
    fn sharding_by_cell_sample_matches_global_build_decisions() {
        let cells: Vec<Dna> = [
            "AAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCC",
            "GGGGGGGGGGGGGGGG",
            "ACGTACGTACGTACGT",
            "TGCATGCATGCATGCA",
            "AACCGGTTAACCGGTT",
        ]
        .iter()
        .map(|s| Dna::pack(s).unwrap())
        .collect();
        let umis: Vec<Dna> = [
            "AAAAAAAAAAAA",
            "AAAAAAAAAAAC",
            "GGGGGGGGGGGG",
            "TTTTTTTTTTTT",
        ]
        .iter()
        .map(|s| Dna::pack(s).unwrap())
        .collect();

        // A synthetic but real-shaped evidence set: several cells, some with
        // more than one sample (matching the pooled-Flex architecture this
        // whole fix targets), several genes, a within-group UMI-correction
        // pair (umis[0]/umis[1] differ by 1 base), and a cross-gene
        // collision (umis[2] shared by two genes in the same cell+sample).
        let mut evidence = HashMap::new();
        let mut next_ordinal = 0_u64;
        let mut push = |evidence: &mut HashMap<RawKey, Evidence>,
                        cell: Dna,
                        sample: u32,
                        gene: u32,
                        umi: Dna,
                        reads: u32| {
            let ordinal = next_ordinal;
            next_ordinal += 1;
            evidence.insert(
                RawKey {
                    cell,
                    sample,
                    gene,
                    umi,
                },
                Evidence {
                    reads,
                    best_ordinal: ordinal,
                    best_mapq: 60,
                    // Distinct per entry (real reads essentially never tie
                    // exactly): an exact rank tie makes `molecule_aggregate`'s
                    // representative choice depend on HashMap iteration
                    // order, which a shard's own map instance need not share
                    // with the global one -- a real but orthogonal,
                    // pre-existing property of tie-breaking, not something
                    // this test is meant to exercise.
                    best_aligned: 100 + ordinal as u32,
                    min_read_hash: ordinal as u32,
                },
            );
        };
        for (i, &cell) in cells.iter().enumerate() {
            let sample = (i % 3) as u32; // several cells reuse the same sample id, and vice versa
            push(&mut evidence, cell, sample, 0, umis[0], 10); // major
            push(&mut evidence, cell, sample, 0, umis[1], 2); // 1bp error of major -> corrects into it
            push(&mut evidence, cell, sample, 1, umis[2], 4); // gene 1, umi X
            push(&mut evidence, cell, sample, 2, umis[2], 3); // gene 2, same umi X, same (cell,sample) -> real cross-gene collision
                                                              // Same cell, a DIFFERENT sample, reusing the same UMI+a new gene:
                                                              // must NOT be pulled into the collision group above.
            push(&mut evidence, cell, (sample + 1) % 3, 3, umis[2], 5);
            push(&mut evidence, cell, sample, 4, umis[3], 1); // plain unique molecule
        }

        let mut global_metrics = Metrics::default();
        let (global_decisions, _) =
            build_decisions(&evidence, &mut global_metrics, None, 1, false, false);

        const SHARDS: u64 = 4;
        let mut shard_maps: Vec<HashMap<RawKey, Evidence>> =
            (0..SHARDS).map(|_| HashMap::new()).collect();
        for (&key, &item) in &evidence {
            // Mirrors shard_hash's real invariant (route by (cell, sample)
            // alone) without needing raw bytes -- any valid partition
            // function proves the same property.
            let idx = (key.cell.0 ^ (key.sample as u64).wrapping_mul(0x9E3779B97F4A7C15)) % SHARDS;
            shard_maps[idx as usize].insert(key, item);
        }

        let mut sharded_decisions = HashMap::new();
        let mut sharded_metrics = Metrics::default();
        for shard in &shard_maps {
            if shard.is_empty() {
                continue;
            }
            let (decisions, _) =
                build_decisions(shard, &mut sharded_metrics, None, 1, false, false);
            sharded_decisions.extend(decisions);
        }

        assert_eq!(sharded_decisions.len(), global_decisions.len());
        for (key, global_decision) in &global_decisions {
            assert_eq!(
                sharded_decisions.get(key),
                Some(global_decision),
                "decision for {key:?} differs between sharded and global runs"
            );
        }
        assert_eq!(
            sharded_metrics.molecules_counted,
            global_metrics.molecules_counted
        );
        assert_eq!(
            sharded_metrics.collision_groups,
            global_metrics.collision_groups
        );
        assert_eq!(
            sharded_metrics.umi_nodes_collapsed,
            global_metrics.umi_nodes_collapsed
        );
        assert_eq!(
            sharded_metrics.raw_cell_gene_umis,
            global_metrics.raw_cell_gene_umis
        );
        // Sanity: this synthetic input actually exercises correction and
        // collision resolution, not just a trivial all-unique case.
        assert!(global_metrics.umi_nodes_collapsed > 0);
        assert!(global_metrics.collision_groups > 0);
    }

    #[test]
    fn shard_hash_is_stable_and_distinguishes_cb_pb_boundary() {
        let a = shard_hash(b"AB", Some(b"C"));
        let b = shard_hash(b"A", Some(b"BC"));
        assert_ne!(
            a, b,
            "CB/PB boundary must affect the hash, not just concatenation"
        );
        assert_eq!(
            a,
            shard_hash(b"AB", Some(b"C")),
            "hash must be deterministic"
        );
        assert_eq!(
            shard_hash(b"AB", None),
            shard_hash(b"AB", Some(b"")),
            "no PB tag and an empty PB tag hash the same"
        );
    }

    /// Core correctness guarantee behind reading multiple un-merged input
    /// BAMs directly (the real motivation: skip samtools cat/sort/index
    /// entirely, since collate() only ever does a plain sequential scan --
    /// see the `bam` CLI field doc comment). Splitting the exact same
    /// records across two "chunk" files must route every record into the
    /// same shard, in the same multiset, as collating from one file holding
    /// all of them -- collate() doesn't care how many files or what order
    /// the records arrive in, only their CB+PB content.
    #[test]
    fn collate_from_multiple_inputs_matches_collate_from_one() {
        let dir = tempfile::tempdir().unwrap();
        let all_records = [
            ("r1", 1, "AAAACCCCGGGGTTTTAAAA", "AAAAAAAAAAAA", "ENSG1"),
            ("r2", 2, "AAAACCCCGGGGTTTTAAAA", "CCCCCCCCCCCC", "ENSG1"),
            ("r3", 3, "CCCCGGGGTTTTAAAACCCC", "GGGGGGGGGGGG", "ENSG1"),
            ("r4", 4, "GGGGTTTTAAAACCCCGGGG", "TTTTTTTTTTTT", "ENSG1"),
            ("r5", 5, "TTTTAAAACCCCGGGGTTTT", "AAAACCCCGGGG", "ENSG1"),
        ];

        // One pre-merged file, all records together.
        let merged_path = dir.path().join("merged.bam");
        write_test_bam(&merged_path, &all_records);
        let merged_shards: Vec<PathBuf> = (0..3)
            .map(|i| dir.path().join(format!("merged_shard_{i}.bam")))
            .collect();
        collate(&[merged_path], 1, &merged_shards).unwrap();

        // The same records, split across three un-merged "chunk" files
        // (simulating raw per-chunk probeset-align output, no cat/sort/index).
        let chunk_paths: Vec<PathBuf> = (0..3)
            .map(|i| dir.path().join(format!("chunk_{i}.bam")))
            .collect();
        write_test_bam(&chunk_paths[0], &all_records[0..2]);
        write_test_bam(&chunk_paths[1], &all_records[2..3]);
        write_test_bam(&chunk_paths[2], &all_records[3..5]);
        let split_shards: Vec<PathBuf> = (0..3)
            .map(|i| dir.path().join(format!("split_shard_{i}.bam")))
            .collect();
        collate(&chunk_paths, 1, &split_shards).unwrap();

        for i in 0..3 {
            assert_eq!(
                read_tag_set(&merged_shards[i]),
                read_tag_set(&split_shards[i]),
                "shard {i} differs between single-file and multi-file collate"
            );
        }
    }
}
