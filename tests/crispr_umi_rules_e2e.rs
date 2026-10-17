//! Flex CRISPR chain with the SBX UMI rules: guide BAM -> `molecules --flex`
//! -> `matrix`, each stage the real binary on the previous stage's output.
//!
//! The fixture (`common::guides`) is shaped like real `sbxpress guides`
//! output and carries the errors that made CAT06 SBX CRISPR overcount 6.1x:
//! 1-base UMI indels and garbled 1-read UMIs in saturated families. Each
//! rule must remove exactly its own error class, shallow families' genuine
//! 1-read molecules must survive, and the matrix must agree with molecules.
//! (A separate file only because pipeline_chain_e2e.rs is at the size limit.)
mod common;

use std::collections::BTreeMap;
use std::path::Path;

use common::guides::{write_guide_fixture, GARBLED_PER_DEEP_FAMILY};
use common::{bin, molecule_counters, read_bam, read_mex, run};
use rust_htslib::bam::record::Aux;

fn molecules() -> std::path::PathBuf {
    bin("sbx10x-molecules", env!("CARGO_BIN_EXE_sbx10x-molecules"))
}
fn matrix() -> std::path::PathBuf {
    bin("cr-bam-to-mtx", env!("CARGO_BIN_EXE_cr-bam-to-mtx"))
}

/// Run molecules (with `args`) then matrix; return counters, matrix total,
/// and the folded-singleton count.
fn chain(dir: &Path, input: &Path, name: &str, args: &[&str]) -> (BTreeMap<String, u64>, f64, u64) {
    let out = dir.join(format!("{name}.molecules.bam"));
    let mut a = vec![
        input.to_str().unwrap(),
        "-o",
        out.to_str().unwrap(),
        "--flex",
    ];
    a.extend_from_slice(args);
    run(&molecules(), &a);
    let metrics_path = dir.join(format!("{name}.molecules.bam.metrics.json"));
    let counters = molecule_counters(&metrics_path);
    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&metrics_path).unwrap()).unwrap();
    let folded = json["umi_singletons_folded"].as_u64().unwrap();

    let records = read_bam(&out);
    assert_eq!(
        records.len(),
        read_bam(input).len(),
        "every record is written"
    );
    let reps = records
        .iter()
        .filter(|r| matches!(r.aux(b"xf"), Ok(Aux::I32(25))))
        .count() as u64;
    assert_eq!(
        reps, counters["molecules_counted"],
        "{name}: one representative per molecule"
    );

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
    let (counts, _) = read_mex(&outs.join("raw_feature_bc_matrix"));
    let total: f64 = counts.values().sum();
    assert_eq!(
        total as u64, counters["molecules_counted"],
        "{name}: matrix agrees with molecules"
    );
    (counters, total, folded)
}

#[test]
fn crispr_umi_rules_remove_sbx_error_molecules_only() {
    let dir = tempfile::tempdir().unwrap();
    let fx = write_guide_fixture(dir.path(), 0x5EED_C415);
    assert!(
        fx.deep_families > 0 && fx.true_singletons > 0,
        "fixture has both regimes"
    );
    let garbled = (fx.deep_families * GARBLED_PER_DEEP_FAMILY) as f64;
    let truth = fx.molecules as f64;
    let deep = fx.deep_molecules as f64;

    // Default: 1-mismatch copies merge, but each deep molecule keeps its two
    // indel variants and every garbled read is its own molecule.
    let (_, base, f0) = chain(dir.path(), &fx.bam, "default", &[]);
    assert_eq!(f0, 0);
    assert_eq!(
        base,
        truth + 2.0 * deep + garbled,
        "default overcounts as SBX did"
    );

    // Indel merge alone removes exactly the indel variants.
    let (_, indel, f1) = chain(dir.path(), &fx.bam, "indel", &["--umi-indel-correction"]);
    assert_eq!(f1, 0);
    assert_eq!(indel, truth + garbled);

    // Fold alone removes 1-read molecules of deep families only: the 1-read
    // indel variants and garbled reads go, the 2-read indel variants stay,
    // and shallow families' genuine singletons are kept.
    let (_, fold, f2) = chain(
        dir.path(),
        &fx.bam,
        "fold",
        &["--fold-singletons-min-top", "10"],
    );
    assert_eq!(f2 as f64, deep + garbled);
    assert_eq!(fold, truth + deep);

    // Both: the true molecule count, singletons included.
    let (_, both, _) = chain(
        dir.path(),
        &fx.bam,
        "both",
        &["--umi-indel-correction", "--fold-singletons-min-top", "10"],
    );
    assert_eq!(both, truth, "both rules recover the true molecules");

    // A threshold above every family's top leaves the counts untouched.
    let (_, high, f3) = chain(
        dir.path(),
        &fx.bam,
        "high",
        &["--fold-singletons-min-top", "1000"],
    );
    assert_eq!((high, f3), (base, 0));
}
