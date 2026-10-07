//! Precision and invalid-count regression cases for weighted matrices.
use super::{write_h5, write_matrix, Feature, Matrix};
use std::fs;

fn matrix(counts: &[f64]) -> Matrix {
    Matrix {
        features: counts
            .iter()
            .enumerate()
            .map(|(i, _)| Feature {
                id: format!("G{i}"),
                name: format!("g{i}"),
            })
            .collect(),
        barcodes: vec!["AAAC-1".into()],
        counts: counts
            .iter()
            .enumerate()
            .map(|(i, &v)| ((0, i), v))
            .collect(),
        records: 0,
        molecules: 0,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    }
}

#[test]
fn tiny_fractional_counts_and_precision_survive_mex_and_hdf5() {
    let values = [
        0.000001_f32 as f64,
        0.1234567_f32 as f64,
        1.000001_f32 as f64,
    ];
    let matrix = matrix(&values);
    let dir = tempfile::tempdir().unwrap();
    write_matrix(&matrix, dir.path(), false).unwrap();
    let text = fs::read_to_string(dir.path().join("matrix.mtx")).unwrap();
    assert!(text.lines().next().unwrap().contains(" real "));
    let mex: Vec<f64> = text
        .lines()
        .filter(|l| !l.starts_with('%'))
        .skip(1)
        .map(|l| l.split_whitespace().nth(2).unwrap().parse().unwrap())
        .collect();
    assert_eq!(mex, values);
    let path = dir.path().join("matrix.h5");
    write_h5(&matrix, &path, "GRCh38").unwrap();
    let h5 = hdf5::File::open(path)
        .unwrap()
        .dataset("matrix/data")
        .unwrap()
        .read_raw::<f64>()
        .unwrap();
    assert_eq!(h5, mex);
}

#[test]
fn negative_and_nonfinite_matrix_counts_fail_before_creating_outputs() {
    let dir = tempfile::tempdir().unwrap();
    for value in [-0.25, f64::NAN, f64::INFINITY] {
        let matrix = matrix(&[value]);
        assert!(write_matrix(&matrix, &dir.path().join("mex"), false).is_err());
        assert!(write_h5(&matrix, &dir.path().join("out.h5"), "G").is_err());
        assert!(!dir.path().join("mex").exists());
        assert!(!dir.path().join("out.h5").exists());
    }
}

#[test]
fn malformed_csc_pointers_and_indices_return_errors_instead_of_panicking() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("raw.h5");
    for (name, data) in [
        ("indptr", vec![0_u32, 2]),
        ("indptr", vec![1, 1]),
        ("indices", vec![7]),
    ] {
        write_h5(&matrix(&[0.25]), &path, "G").unwrap();
        {
            let file = hdf5::File::open_rw(&path).unwrap();
            file.dataset(&format!("matrix/{name}"))
                .unwrap()
                .write_raw(&data)
                .unwrap();
        }
        assert!(
            crate::cell_calling::read_raw(&path).is_err(),
            "{name}: {data:?}"
        );
    }
}
