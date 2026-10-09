//! Ground-truth comparison of the molecule-counting policies on a
//! tag-seq-shaped fixture: a heavily oversubscribed UMI space (so distinct
//! molecules often share a UMI, within and across genes), PCR duplicates
//! that share their 5' start, and 1-mismatch UMI errors on duplicates.
//!
//! The fixture knows the truth: every read's name carries its source
//! molecule, so the true count per (cell, gene) is the number of distinct
//! molecules with at least one countable read. Each policy's matrix is
//! scored against it. This is a simulation -- it encodes the assumption
//! that UMI collisions are chance events between distinct molecules -- so it
//! proves each mode does what it claims, not which mode real data needs
//! (scripts/compare_tagseq_modes.sh does that against an Illumina matrix).
mod common;

use std::collections::{BTreeMap, BTreeSet};

use common::{bin, molecule_id, read_bam, read_mex, run, write_fixture, FixtureSpec};
use rust_htslib::bam::record::Aux;

fn tag(record: &rust_htslib::bam::Record, name: &[u8]) -> Option<String> {
    match record.aux(name) {
        Ok(Aux::String(v)) => Some(v.to_owned()),
        _ => None,
    }
}

/// (barcode, gene) -> number of distinct source molecules with a countable
/// read (mapped, MAPQ >= 60 as `molecules --min-mapq` requires -- the fixture
/// writes no `mm` rescues -- with CB; gene = the single GX or the ZG-resolved
/// one).
fn truth(bam: &std::path::Path) -> BTreeMap<(String, String), f64> {
    let mut molecules: BTreeMap<(String, String), BTreeSet<u64>> = BTreeMap::new();
    for r in read_bam(bam) {
        if r.is_unmapped() || r.mapq() < 60 {
            continue;
        }
        let Some(cb) = tag(&r, b"CB") else { continue };
        let gx = tag(&r, b"GX").unwrap();
        let gene = if gx.contains(';') {
            tag(&r, b"ZG").unwrap()
        } else {
            gx
        };
        molecules
            .entry((cb, gene))
            .or_default()
            .insert(molecule_id(&r));
    }
    molecules
        .into_iter()
        .map(|(k, v)| (k, v.len() as f64))
        .collect()
}

/// Sum of |estimate - truth| over all (cell, gene), relative to total truth.
fn relative_error(
    est: &BTreeMap<(String, String), f64>,
    truth: &BTreeMap<(String, String), f64>,
) -> f64 {
    let keys: BTreeSet<_> = est.keys().chain(truth.keys()).collect();
    let err: f64 = keys
        .iter()
        .map(|k| {
            (est.get(*k).copied().unwrap_or(0.0) - truth.get(*k).copied().unwrap_or(0.0)).abs()
        })
        .sum();
    err / truth.values().sum::<f64>()
}

#[test]
fn position_aware_mode_recovers_true_molecule_counts() {
    let dir = tempfile::tempdir().unwrap();
    let fx = write_fixture(
        dir.path(),
        FixtureSpec {
            seed: 0x7A65_0001,
            molecules: 6_000,
            umi_pool: 40, // ~150 molecules per UMI: deep oversubscription
            shared_start: true,
            ..FixtureSpec::default()
        },
    );
    let truth = truth(&fx.bam);
    let total_truth: f64 = truth.values().sum();
    let molecules = bin("sbx10x-molecules", env!("CARGO_BIN_EXE_sbx10x-molecules"));
    let matrix = bin("cr-bam-to-mtx", env!("CARGO_BIN_EXE_cr-bam-to-mtx"));

    let mut errors = BTreeMap::new();
    let mut totals = BTreeMap::new();
    for (name, args) in [
        ("single_cell_default", vec![]),
        ("tagseq_shared", vec!["--tagseq"]),
        (
            "tagseq_fractional",
            vec!["--tagseq", "--fractional-collisions"],
        ),
        (
            "tagseq_position_aware",
            vec!["--tagseq", "--position-aware-umis"],
        ),
        (
            "tagseq_position_aware_sharded",
            vec!["--tagseq", "--position-aware-umis", "--shards", "3"],
        ),
    ] {
        let out = dir.path().join(format!("{name}.bam"));
        let outs = dir.path().join(format!("{name}.outs"));
        let mut a = vec![fx.bam.to_str().unwrap(), "-o", out.to_str().unwrap()];
        a.extend(args);
        run(&molecules, &a);
        run(
            &matrix,
            &[
                out.to_str().unwrap(),
                "-o",
                outs.to_str().unwrap(),
                "--no-gzip",
            ],
        );
        let (est, _) = read_mex(&outs.join("raw_feature_bc_matrix"));
        errors.insert(name, relative_error(&est, &truth));
        totals.insert(name, est.values().sum::<f64>());
    }
    eprintln!("true molecules: {total_truth}");
    for (name, err) in &errors {
        eprintln!(
            "{name:32} total {:8.1}  relative L1 error {:6.2}%",
            totals[name],
            err * 100.0
        );
    }
    let pos = errors["tagseq_position_aware"];
    assert!(
        pos < 0.01,
        "position-aware error {:.2}% (expected <1%)",
        pos * 100.0
    );
    assert_eq!(
        pos, errors["tagseq_position_aware_sharded"],
        "sharding changes nothing"
    );
    for other in ["single_cell_default", "tagseq_shared", "tagseq_fractional"] {
        assert!(
            pos < errors[other],
            "position-aware should beat {other}: {errors:?}"
        );
    }
    // Under chance collisions, fractional undercounts and shared does not
    // undercount collisions (it may overcount from uncorrected UMI errors).
    assert!(totals["tagseq_fractional"] < total_truth);
    assert!(totals["tagseq_fractional"] < totals["tagseq_shared"]);
}

#[test]
fn position_aware_umis_keep_ub_clean_and_reject_overlong_umis() {
    let dir = tempfile::tempdir().unwrap();
    let fx = write_fixture(
        dir.path(),
        FixtureSpec {
            shared_start: true,
            ..FixtureSpec::default()
        },
    );
    let out = dir.path().join("pos.bam");
    run(
        &bin("sbx10x-molecules", env!("CARGO_BIN_EXE_sbx10x-molecules")),
        &[
            fx.bam.to_str().unwrap(),
            "-o",
            out.to_str().unwrap(),
            "--position-aware-umis",
        ],
    );
    // The position tag never leaks into the written UMI.
    for r in read_bam(&out) {
        if let Some(ub) = tag(&r, b"UB") {
            assert_eq!(ub.len(), 12);
            assert!(ub.bytes().all(|b| b"ACGT".contains(&b)));
        }
    }
    let json: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(dir.path().join("pos.bam.metrics.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(json["position_aware_umis"], true);

    // A 20bp UMI cannot carry the position tag: fail loudly, never drop reads.
    let long = dir.path().join("long.bam");
    {
        use rust_htslib::bam::{self, Read as _};
        let mut reader = bam::Reader::from_path(&fx.bam).unwrap();
        let header = bam::Header::from_template(reader.header());
        let mut writer = bam::Writer::from_path(&long, &header, bam::Format::Bam).unwrap();
        for r in reader.records() {
            let mut r = r.unwrap();
            if r.remove_aux(b"UR").is_ok() {
                r.push_aux(b"UR", Aux::String("ACGTACGTACGTACGTACGT"))
                    .unwrap();
            }
            writer.write(&r).unwrap();
        }
    }
    let failed = std::process::Command::new(bin(
        "sbx10x-molecules",
        env!("CARGO_BIN_EXE_sbx10x-molecules"),
    ))
    .args([
        long.to_str().unwrap(),
        "-o",
        dir.path().join("x.bam").to_str().unwrap(),
        "--position-aware-umis",
    ])
    .output()
    .unwrap();
    assert!(!failed.status.success());
    assert!(String::from_utf8_lossy(&failed.stderr).contains("at most 16 bases"));
}

/// scripts/compare_tagseq_modes.sh end to end, with the fixture's true
/// counts as the "reference" matrix: the tool must rank position-aware
/// closest to truth and show fractional's expression-dependent undercount.
#[test]
fn compare_tagseq_modes_tool_ranks_modes_against_a_reference() {
    let numpy = std::process::Command::new("python3")
        .args(["-c", "import numpy"])
        .output();
    if !numpy.is_ok_and(|o| o.status.success()) {
        eprintln!("python3 with numpy not available; skipping");
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let fx = write_fixture(
        dir.path(),
        FixtureSpec {
            seed: 0x7A65_0002,
            molecules: 6_000,
            umi_pool: 40,
            shared_start: true,
            ..FixtureSpec::default()
        },
    );
    // Truth as a MEX reference.
    let reference = dir.path().join("truth");
    std::fs::create_dir_all(&reference).unwrap();
    let truth = truth(&fx.bam);
    let barcodes: Vec<&String> = truth
        .keys()
        .map(|(b, _)| b)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let genes: Vec<&String> = truth
        .keys()
        .map(|(_, g)| g)
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    std::fs::write(
        reference.join("barcodes.tsv"),
        barcodes
            .iter()
            .map(|b| format!("{b}\n"))
            .collect::<String>(),
    )
    .unwrap();
    std::fs::write(
        reference.join("features.tsv"),
        genes
            .iter()
            .map(|g| format!("{g}\t{g}\tGene Expression\n"))
            .collect::<String>(),
    )
    .unwrap();
    let mut mtx = format!(
        "%%MatrixMarket matrix coordinate integer general\n{} {} {}\n",
        genes.len(),
        barcodes.len(),
        truth.len()
    );
    for ((b, g), v) in &truth {
        let gi = genes.iter().position(|x| *x == g).unwrap() + 1;
        let bi = barcodes.iter().position(|x| *x == b).unwrap() + 1;
        mtx.push_str(&format!("{gi} {bi} {v}\n"));
    }
    std::fs::write(reference.join("matrix.mtx"), mtx).unwrap();

    let out = dir.path().join("compare");
    let result = std::process::Command::new("bash")
        .arg(
            std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("scripts/compare_tagseq_modes.sh"),
        )
        .args([
            "--bam",
            fx.bam.to_str().unwrap(),
            "--reference",
            reference.to_str().unwrap(),
        ])
        .args(["--out-dir", out.to_str().unwrap(), "--threads", "1"])
        .env("SBXPRESS", env!("CARGO_BIN_EXE_sbxpress"))
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
    let summary = std::fs::read_to_string(out.join("summary.tsv")).unwrap();
    eprintln!("{summary}");
    let rows: BTreeMap<String, BTreeMap<String, f64>> = {
        let mut lines = summary.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        lines
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                let values = header
                    .iter()
                    .zip(&f)
                    .skip(1)
                    .map(|(h, v)| (h.to_string(), v.parse().unwrap()))
                    .collect();
                (f[0].to_owned(), values)
            })
            .collect()
    };
    let total = |m: &str| rows[m]["total_vs_reference"];
    let err = |m: &str| rows[m]["median_abs_log2_cpm_ratio"];
    assert!(
        (total("position") - 1.0).abs() < 0.01,
        "position total {}",
        total("position")
    );
    assert!(err("position") <= err("shared") && err("position") < err("fractional"));
    assert!(total("fractional") < total("shared"));
    // Collision loss grows with expression: fractional's ratio falls with level.
    assert!(
        rows["fractional"]["log2_ratio_vs_expression_slope"]
            < rows["position"]["log2_ratio_vs_expression_slope"]
    );
}
