use std::fs;
use crate::cr_bam_to_mtx::{write_h5, Feature, Matrix};
use super::algo::{bh_adjust, run};
use super::helpers::within_ordmag;

use super::types::Cli;

#[test]
fn bh_is_monotone_in_rank() {
    let q = bh_adjust(&[0.01, 0.04, 0.03, 0.002]);
    assert!((q[3] - 0.008).abs() < 1e-12);
    assert!(q[0] <= q[2] && q[2] <= q[1]);
}

#[test]
fn ordmag_threshold_uses_top_one_percent() {
    let mut x = vec![1; 1000];
    x.extend([1000; 10]);
    x.sort_unstable();
    assert_eq!(within_ordmag(&x, 1), 10);
}

#[test]
fn force_cells_writes_filtered_outputs() {
    let temp = tempfile::tempdir().unwrap();
    let raw_path = temp.path().join("raw.h5");
    let matrix = Matrix {
        features: vec![Feature {
            id: "g1".into(),
            name: "G1".into(),
        }],
        barcodes: vec!["bc1-1".into(), "bc2-1".into()],
        counts: [((0, 0), 10.0), ((1, 0), 2.0)].into_iter().collect(),
        records: 0,
        molecules: 12,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    };
    write_h5(&matrix, &raw_path, "GRCh38").unwrap();
    let output = temp.path().join("outs");
    run(&Cli {
        raw_matrix: raw_path,
        output: output.clone(),
        expect_cells: None,
        force_cells: Some(1),
        tagseq_barcodes: None,
        ambient_lower: 1,
        ambient_upper: 2,
        emptydrops_minimum_umis: 1,
        simulations: 10,
        fdr: 0.001,
        seed: 0,
    })
    .unwrap();
    assert!(output.join("filtered_feature_bc_matrix.h5").is_file());
    assert_eq!(
        fs::read_to_string(output.join("filtered_barcodes.csv")).unwrap(),
        "genome,barcode\nGRCh38,bc1-1\n"
    );
}

#[test]
fn tagseq_barcodes_selects_exact_list_ignoring_umi_rank() {
    let temp = tempfile::tempdir().unwrap();
    let raw_path = temp.path().join("raw.h5");
    let matrix = Matrix {
        features: vec![Feature {
            id: "g1".into(),
            name: "G1".into(),
        }],
        barcodes: vec!["bc1-1".into(), "bc2-1".into()],
        counts: [((0, 0), 10.0), ((1, 0), 2.0)].into_iter().collect(),
        records: 0,
        molecules: 12,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    };
    write_h5(&matrix, &raw_path, "GRCh38").unwrap();
    let barcodes_path = temp.path().join("samples.txt");
    fs::write(&barcodes_path, "bc2\n").unwrap();
    let output = temp.path().join("outs");
    run(&Cli {
        raw_matrix: raw_path,
        output: output.clone(),
        expect_cells: None,
        force_cells: None,
        tagseq_barcodes: Some(barcodes_path),
        ambient_lower: 1,
        ambient_upper: 2,
        emptydrops_minimum_umis: 1,
        simulations: 10,
        fdr: 0.001,
        seed: 0,
    })
    .unwrap();
    assert_eq!(
        fs::read_to_string(output.join("filtered_barcodes.csv")).unwrap(),
        "genome,barcode\nGRCh38,bc2-1\n"
    );
    let calls = fs::read_to_string(output.join("cell_calls.csv")).unwrap();
    assert!(calls.contains("bc2-1,2,tagseq"));
    let metrics_json = fs::read_to_string(output.join("cell_calling_metrics.json")).unwrap();
    assert!(metrics_json.contains("explicit tag-seq sample list"));
}

#[test]
fn tagseq_barcodes_rejects_unknown_barcode() {
    let temp = tempfile::tempdir().unwrap();
    let raw_path = temp.path().join("raw.h5");
    let matrix = Matrix {
        features: vec![Feature {
            id: "g1".into(),
            name: "G1".into(),
        }],
        barcodes: vec!["bc1-1".into()],
        counts: [((0, 0), 10.0)].into_iter().collect(),
        records: 0,
        molecules: 10,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    };
    write_h5(&matrix, &raw_path, "GRCh38").unwrap();
    let barcodes_path = temp.path().join("samples.txt");
    fs::write(&barcodes_path, "not_in_the_matrix\n").unwrap();
    let output = temp.path().join("outs");
    let err = run(&Cli {
        raw_matrix: raw_path,
        output,
        expect_cells: None,
        force_cells: None,
        tagseq_barcodes: Some(barcodes_path),
        ambient_lower: 1,
        ambient_upper: 2,
        emptydrops_minimum_umis: 1,
        simulations: 10,
        fdr: 0.001,
        seed: 0,
    })
    .unwrap_err();
    assert!(err.to_string().contains("not_in_the_matrix"));
}

/// Fractional raw matrices (`molecules --fractional-collisions`) must be
/// read and filtered without truncating 2.5 -> 2 or dropping 0.25 -> 0.
#[test]
fn fractional_raw_matrix_is_read_and_filtered_without_truncation() {
    let temp = tempfile::tempdir().unwrap();
    let raw_path = temp.path().join("raw.h5");
    let matrix = Matrix {
        features: vec![
            Feature { id: "g1".into(), name: "G1".into() },
            Feature { id: "g2".into(), name: "G2".into() },
        ],
        barcodes: vec!["bc1-1".into(), "bc2-1".into()],
        counts: [((0, 0), 2.5), ((0, 1), 0.25), ((1, 0), 1.0)].into_iter().collect(),
        records: 0,
        molecules: 4,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    };
    write_h5(&matrix, &raw_path, "GRCh38").unwrap();
    let raw = super::h5::read_raw(&raw_path).unwrap();
    assert_eq!(raw.data, vec![2.5, 0.25, 1.0]);
    assert_eq!(super::h5::totals(&raw), vec![3, 1]); // round(2.75), round(1.0)

    let output = temp.path().join("outs");
    run(&Cli {
        raw_matrix: raw_path,
        output: output.clone(),
        expect_cells: None,
        force_cells: Some(1),
        tagseq_barcodes: None,
        ambient_lower: 1,
        ambient_upper: 2,
        emptydrops_minimum_umis: 1,
        simulations: 10,
        fdr: 0.001,
        seed: 0,
    })
    .unwrap();
    let filtered = super::h5::read_raw(&output.join("filtered_feature_bc_matrix.h5")).unwrap();
    assert_eq!(filtered.barcodes, vec!["bc1-1".to_string()]);
    assert_eq!(filtered.data, vec![2.5, 0.25]);
}
