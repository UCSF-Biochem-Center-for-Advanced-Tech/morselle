//! Full-chain contract tests: annotated BAM -> `molecules` -> `matrix` ->
//! `cells` (and `probe-demux` for Flex), each stage run as the real binary
//! on the previous stage's real output.
//!
//! Why this exists: before 0.5.9 no test ran `molecules`' BAM writer at all,
//! and every `matrix`/`probe-demux` test built its own `xf:25` records by
//! hand. Each stage was tested against what it *assumed* its producer
//! wrote, so 0.5.7 could stop writing `xf` (and crash on `UB`) with every
//! test still green. These tests check the stage-to-stage invariants on
//! actual outputs; tests/golden_outputs.rs pins the exact values.
mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use common::{
    barcode_table, bin, molecule_counters, read_bam, read_mex, run, write_fixture, FixtureSpec,
};
use rust_htslib::bam::{self, record::Aux};

fn molecules() -> std::path::PathBuf {
    bin("sbx10x-molecules", env!("CARGO_BIN_EXE_sbx10x-molecules"))
}
fn matrix() -> std::path::PathBuf {
    bin("cr-bam-to-mtx", env!("CARGO_BIN_EXE_cr-bam-to-mtx"))
}
fn cells() -> std::path::PathBuf {
    bin("sbxpress-cells", env!("CARGO_BIN_EXE_sbxpress-cells"))
}
fn probe_demux() -> std::path::PathBuf {
    bin(
        "sbxpress-probe-demux",
        env!("CARGO_BIN_EXE_sbxpress-probe-demux"),
    )
}

fn s(record: &bam::Record, tag: &[u8]) -> Option<String> {
    match record.aux(tag) {
        Ok(Aux::String(v)) => Some(v.to_owned()),
        _ => None,
    }
}
fn xf(record: &bam::Record) -> Option<i32> {
    match record.aux(b"xf") {
        Ok(Aux::I32(v)) => Some(v),
        _ => None,
    }
}
fn xw(record: &bam::Record) -> Option<u8> {
    match record.aux(b"XW") {
        Ok(Aux::Char(v)) => Some(v),
        _ => None,
    }
}

/// Invariants every `molecules` output must satisfy, whatever the policy.
fn check_molecule_bam(input: &Path, output: &Path, metrics: &BTreeMap<String, u64>) {
    let before = read_bam(input);
    let after = read_bam(output);
    assert_eq!(before.len(), after.len(), "every input record is written");
    for (a, b) in before.iter().zip(&after) {
        assert_eq!(a.qname(), b.qname(), "input order is preserved");
        assert_eq!(
            a.seq().as_bytes(),
            b.seq().as_bytes(),
            "sequence is preserved"
        );
        assert_eq!(a.qual(), b.qual(), "qualities are preserved");
        assert_eq!(
            a.cigar().to_string(),
            b.cigar().to_string(),
            "CIGAR is preserved"
        );
        assert_eq!(
            (a.mapq(), a.mtid(), a.mpos(), a.insert_size()),
            (b.mapq(), b.mtid(), b.mpos(), b.insert_size()),
            "alignment and mate fields are preserved"
        );
        // annotate's tags survive untouched.
        assert_eq!(s(a, b"ZG"), s(b, b"ZG"), "ZG belongs to annotate");
        assert_eq!(s(a, b"UR"), s(b, b"UR"), "UR is the raw UMI");
        if let Some(zg) = s(a, b"ZG") {
            if xf(b).is_some() {
                assert_eq!(
                    s(b, b"GX"),
                    Some(zg),
                    "multi-gene GX resolved to the ZG gene"
                );
            }
        }
    }
    let reps: Vec<&bam::Record> = after.iter().filter(|r| xf(r) == Some(25)).collect();
    assert!(!reps.is_empty(), "molecules produced representatives");
    assert_eq!(reps.len() as u64, metrics["representatives_tagged"]);
    for r in &reps {
        assert!(s(r, b"CB").is_some() && s(r, b"UB").is_some());
        assert!(
            !s(r, b"GX").unwrap().contains(';'),
            "counted record has one gene"
        );
        assert!(matches!(xw(r), Some(b'U' | b'W' | b'S' | b'F')));
    }
    for r in &after {
        assert_eq!(
            r.aux(b"ZW").is_ok(),
            xw(r) == Some(b'F'),
            "ZW only on fractional molecules"
        );
    }
}

/// Matrix totals must equal what molecules says it counted.
fn check_matrix_total(counts: &BTreeMap<(String, String), f64>, metrics: &BTreeMap<String, u64>) {
    let total: f64 = counts.values().sum();
    assert!(
        (total - metrics["molecules_counted"] as f64).abs() < 1e-3,
        "matrix total {total} != molecules_counted {}",
        metrics["molecules_counted"]
    );
}

/// Filtered matrix = raw matrix restricted to the called cells, unchanged.
fn check_filtered(raw: &BTreeMap<(String, String), f64>, filtered_dir: &Path) {
    let (filtered, _) = read_mex(filtered_dir);
    let called: BTreeSet<&String> = filtered.keys().map(|(bc, _)| bc).collect();
    assert!(!called.is_empty(), "cells called something");
    for ((bc, feature), value) in &filtered {
        let raw_value = raw[&(bc.clone(), feature.clone())];
        assert!(
            (value - raw_value).abs() < 1e-4,
            "{bc}/{feature}: {value} vs raw {raw_value}"
        );
    }
    for ((bc, feature), value) in raw {
        if called.contains(bc) {
            assert!(
                filtered.contains_key(&(bc.clone(), feature.clone())),
                "lost {bc}/{feature}={value}"
            );
        }
    }
}

fn molecules_matrix(
    dir: &Path,
    input: &Path,
    name: &str,
    args: &[&str],
) -> (
    BTreeMap<String, u64>,
    BTreeMap<(String, String), f64>,
    bool,
    std::path::PathBuf,
) {
    let out = dir.join(format!("{name}.molecules.bam"));
    let mut a = vec![
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--barcode-metrics",
    ];
    a.extend_from_slice(args);
    run(&molecules(), &a);
    let metrics = molecule_counters(&dir.join(format!("{name}.molecules.bam.metrics.json")));
    check_molecule_bam(input, &out, &metrics);
    let outs = dir.join(format!("{name}.outs"));
    run(
        &matrix(),
        &[
            out.to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--no-gzip",
        ],
    );
    let (counts, is_real) = read_mex(&outs.join("raw_feature_bc_matrix"));
    let table_path = dir.join(format!("{name}.molecules.bam.metrics.barcodes.tsv.gz"));
    let rows = barcode_table(&table_path);
    for (column, metric) in [
        ("eligible_reads", "eligible_records"),
        ("raw_cell_gene_umis", "raw_cell_gene_umis"),
        ("molecules_counted", "molecules_counted"),
        ("collision_groups", "collision_groups"),
        ("representatives_tagged", "representatives_tagged"),
    ] {
        assert_eq!(
            rows.values().map(|r| r[column]).sum::<u64>(),
            metrics[metric],
            "{column}"
        );
    }
    let mut per_cb = BTreeMap::<String, u64>::new();
    for ((cb, _sample), row) in &rows {
        *per_cb.entry(cb.clone()).or_default() += row["molecules_counted"];
    }
    let mut matrix_by_cb = BTreeMap::<String, f64>::new();
    for ((cb, _gene), value) in &counts {
        *matrix_by_cb.entry(cb.clone()).or_default() += value;
    }
    for (cb, molecules) in per_cb {
        assert!(
            (matrix_by_cb.get(&cb).copied().unwrap_or(0.0) - molecules as f64).abs() < 1e-3,
            "{cb}"
        );
    }
    (metrics, counts, is_real, outs)
}

#[test]
fn gex_chain_counts_corrected_molecules_and_calls_cells() {
    let dir = tempfile::tempdir().unwrap();
    let fx = write_fixture(dir.path(), FixtureSpec::default());
    let (metrics, counts, is_real, outs) = molecules_matrix(dir.path(), &fx.bam, "gex", &[]);
    assert!(!is_real);
    check_matrix_total(&counts, &metrics);
    assert!(
        metrics["umi_nodes_collapsed"] > 0 && metrics["collision_groups"] > 0,
        "fixture exercises correction and collisions"
    );

    // UMI correction is visible as UB != UR on corrected reads.
    let out = read_bam(&dir.path().join("gex.molecules.bam"));
    let corrected: BTreeSet<(String, String, String)> = out
        .iter()
        .filter(|r| xf(r).is_some() && s(r, b"UB") != s(r, b"UR"))
        .map(|r| {
            (
                s(r, b"CB").unwrap(),
                s(r, b"GX").unwrap(),
                s(r, b"UR").unwrap(),
            )
        })
        .collect();
    assert_eq!(corrected.len() as u64, metrics["umi_nodes_collapsed"]);

    let raw_h5 = outs.join("raw_feature_bc_matrix.h5");
    run(
        &cells(),
        &[
            raw_h5.to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--force-cells",
            "12",
        ],
    );
    let filtered = outs.join("filtered_feature_bc_matrix");
    check_filtered(&counts, &filtered);
}

#[test]
fn tagseq_fractional_chain_keeps_fractions_through_cells() {
    let dir = tempfile::tempdir().unwrap();
    let fx = write_fixture(
        dir.path(),
        FixtureSpec {
            umi_pool: 150,
            ..FixtureSpec::default()
        },
    );
    let (metrics, counts, is_real, outs) = molecules_matrix(
        dir.path(),
        &fx.bam,
        "tagseq",
        &["--tagseq", "--fractional-collisions"],
    );
    assert!(is_real, "fractional collisions give a real-valued matrix");
    assert!(metrics["collision_molecules_fractional"] > 0);
    check_matrix_total(&counts, &metrics);

    let out = read_bam(&dir.path().join("tagseq.molecules.bam"));
    // --tagseq never corrects UMIs.
    assert!(out
        .iter()
        .filter(|r| xf(r).is_some())
        .all(|r| s(r, b"UB") == s(r, b"UR")));
    // Fractional weights of one (cell, UMI) collision sum to 1.
    let mut groups: BTreeMap<(String, String), f64> = BTreeMap::new();
    for r in out
        .iter()
        .filter(|r| xf(r) == Some(25) && xw(r) == Some(b'F'))
    {
        let Ok(Aux::Float(w)) = r.aux(b"ZW") else {
            panic!("F without ZW")
        };
        *groups
            .entry((s(r, b"CB").unwrap(), s(r, b"UB").unwrap()))
            .or_default() += w as f64;
    }
    assert!(!groups.is_empty());
    for (key, sum) in &groups {
        assert!((sum - 1.0).abs() < 1e-5, "{key:?} weights sum to {sum}");
    }

    let known = dir.path().join("known_cells.txt");
    std::fs::write(&known, fx.cells[..12].join("\n") + "\n").unwrap();
    let raw_h5 = outs.join("raw_feature_bc_matrix.h5");
    run(
        &cells(),
        &[
            raw_h5.to_str().unwrap(),
            "-o",
            outs.to_str().unwrap(),
            "--tagseq-barcodes",
            known.to_str().unwrap(),
        ],
    );
    let filtered = outs.join("filtered_feature_bc_matrix");
    check_filtered(&counts, &filtered);
    let (f, _) = read_mex(&filtered);
    assert!(
        f.values().any(|v| v.fract() != 0.0),
        "fractions survive cell calling"
    );
}

#[test]
fn tagseq_shared_chain_counts_every_claimant() {
    let dir = tempfile::tempdir().unwrap();
    let fx = write_fixture(
        dir.path(),
        FixtureSpec {
            umi_pool: 150,
            ..FixtureSpec::default()
        },
    );
    let (metrics, counts, is_real, _) = molecules_matrix(
        dir.path(),
        &fx.bam,
        "shared",
        &["--tagseq", "--no-fractional-collisions"],
    );
    assert!(!is_real);
    assert!(metrics["collision_molecules_shared"] > 0);
    check_matrix_total(&counts, &metrics);
}

#[test]
fn flex_chain_demultiplexes_samples_by_probe_barcode() {
    let dir = tempfile::tempdir().unwrap();
    let fx = write_fixture(
        dir.path(),
        FixtureSpec {
            probe_barcodes: Some(3),
            ..FixtureSpec::default()
        },
    );
    let (metrics, _, _, _) = molecules_matrix(dir.path(), &fx.bam, "flex", &["--flex"]);
    let molecules_bam = dir.path().join("flex.molecules.bam");
    let reps: Vec<bam::Record> = read_bam(&molecules_bam)
        .into_iter()
        .filter(|r| xf(r) == Some(25))
        .collect();
    assert!(
        reps.iter().all(|r| s(r, b"PB").is_some()),
        "--flex counts only PB-tagged molecules"
    );
    assert_eq!(reps.len() as u64, metrics["molecules_counted"]);

    let sheet = dir.path().join("samples.csv");
    std::fs::write(
        &sheet,
        "sample_id,probe_barcode_ids\nS1,BC001\nS2,BC002|BC003\n",
    )
    .unwrap();
    let demux = dir.path().join("demux");
    run(
        &probe_demux(),
        &[
            molecules_bam.to_str().unwrap(),
            "--probe-sample-sheet",
            sheet.to_str().unwrap(),
            "-o",
            demux.to_str().unwrap(),
            "--min-molecules",
            "1",
        ],
    );
    for sample in ["S1", "S2"] {
        let barcodes = demux.join(format!("{sample}.barcodes.txt"));
        let probes = demux.join(format!("{sample}.probe_barcodes.txt"));
        let outs = dir.path().join(format!("{sample}.outs"));
        run(
            &matrix(),
            &[
                molecules_bam.to_str().unwrap(),
                "-o",
                outs.to_str().unwrap(),
                "--no-gzip",
                "--cell-barcodes",
                barcodes.to_str().unwrap(),
                "--probe-barcodes",
                probes.to_str().unwrap(),
                "--skip-unlisted-barcodes",
            ],
        );
        let (counts, _) = read_mex(&outs.join("raw_feature_bc_matrix"));
        // Independent recount straight from the molecules BAM.
        let listed: BTreeSet<String> = std::fs::read_to_string(&barcodes)
            .unwrap()
            .lines()
            .map(|l| l.split('-').next().unwrap().to_owned())
            .collect();
        let pbs: BTreeSet<String> = std::fs::read_to_string(&probes)
            .unwrap()
            .lines()
            .map(str::to_owned)
            .collect();
        let expected = reps
            .iter()
            .filter(|r| {
                pbs.contains(&s(r, b"PB").unwrap())
                    && listed.contains(s(r, b"CB").unwrap().split('-').next().unwrap())
            })
            .count();
        assert!(expected > 0, "{sample} has molecules");
        assert_eq!(
            counts.values().sum::<f64>() as usize,
            expected,
            "{sample} matrix total"
        );
    }
}

#[test]
fn compression_level_changes_file_size_not_records() {
    let dir = tempfile::tempdir().unwrap();
    let fx = write_fixture(dir.path(), FixtureSpec::default());
    let mut digests = Vec::new();
    let mut sizes = Vec::new();
    for (name, extra) in [
        ("default", vec![]),
        ("fast", vec!["--compression-level", "1"]),
        ("none", vec!["--compression-level", "0"]),
        // Sharded output is in shard order, so compare it with itself.
        ("sharded", vec!["--shards", "2"]),
        (
            "sharded_fast",
            vec!["--compression-level", "1", "--shards", "2"],
        ),
    ] {
        let out = dir.path().join(format!("{name}.bam"));
        let mut args = vec![fx.bam.to_str().unwrap(), "-o", out.to_str().unwrap()];
        args.extend(extra);
        run(&molecules(), &args);
        digests.push(common::bam_digest(&out));
        sizes.push(std::fs::metadata(&out).unwrap().len());
    }
    assert!(
        digests[..3].iter().all(|d| d == &digests[0]) && digests[3] == digests[4],
        "records identical at every level"
    );
    assert!(
        sizes[2] > sizes[1] && sizes[1] >= sizes[0],
        "sizes {sizes:?}"
    );
    let bad = std::process::Command::new(molecules())
        .args([
            fx.bam.to_str().unwrap(),
            "-o",
            "/dev/null",
            "--compression-level",
            "10",
        ])
        .output()
        .unwrap();
    assert!(!bad.status.success(), "levels above 9 are rejected");
}

#[path = "regressions/pipeline.rs"]
mod regressions;
