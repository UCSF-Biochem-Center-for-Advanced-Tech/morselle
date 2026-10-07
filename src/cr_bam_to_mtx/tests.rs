use std::fs;
use std::path::{Path, PathBuf};
use hdf5::types::VarLenAscii;
use rust_htslib::bam::{self, record::Aux};

use super::io_helpers::read_cell_barcodes;
use super::matrix_ops::{read_matrix, write_h5, write_matrix};
use super::types::{Feature, Matrix};

#[test]
fn tenx_whitelist_gets_gem_group_suffix() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("barcodes.txt");
    fs::write(&path, "AACCGGTT\nTTGGCCAA\n").unwrap();
    assert_eq!(
        read_cell_barcodes(&path, 2).unwrap(),
        ["AACCGGTT-2", "TTGGCCAA-2"]
    );
}

#[test]
fn an_existing_gem_group_suffix_is_accepted_but_gem_group_flag_wins() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("barcodes.txt");
    fs::write(&path, "AACCGGTT-1\nTTGGCCAA-2\n").unwrap();
    assert_eq!(
        read_cell_barcodes(&path, 3).unwrap(),
        ["AACCGGTT-3", "TTGGCCAA-3"]
    );
}

#[test]
fn writes_one_based_sparse_matrix() {
    let matrix = Matrix {
        features: vec![Feature {
            id: "g1".into(),
            name: "G1".into(),
        }],
        barcodes: vec!["bc-1".into()],
        counts: [((0, 0), 2.0)].into_iter().collect(),
        records: 3,
        molecules: 2,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    };
    let dir = tempfile::tempdir().unwrap();
    write_matrix(&matrix, dir.path(), false).unwrap();
    let text = fs::read_to_string(dir.path().join("matrix.mtx")).unwrap();
    assert!(text.ends_with("1 1 1\n1 1 2\n"));
}

#[test]
fn writes_compatibility_barcode_major_order() {
    let matrix = Matrix {
        features: vec![
            Feature {
                id: "g1".into(),
                name: "G1".into(),
            },
            Feature {
                id: "g2".into(),
                name: "G2".into(),
            },
        ],
        barcodes: vec!["bc1-1".into(), "bc2-1".into()],
        counts: [((1, 0), 3.0), ((0, 1), 2.0), ((0, 0), 1.0)]
            .into_iter()
            .collect(),
        records: 6,
        molecules: 6,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    };
    let dir = tempfile::tempdir().unwrap();
    write_matrix(&matrix, dir.path(), false).unwrap();
    let text = fs::read_to_string(dir.path().join("matrix.mtx")).unwrap();
    let data: Vec<_> = text.lines().filter(|line| !line.starts_with('%')).collect();
    assert_eq!(data, ["2 2 3", "1 1 1", "2 1 2", "1 2 3"]);
}

#[test]
fn writes_compatibility_csc_h5() {
    let matrix = Matrix {
        features: vec![
            Feature {
                id: "g1".into(),
                name: "G1".into(),
            },
            Feature {
                id: "g2".into(),
                name: "G2".into(),
            },
        ],
        barcodes: vec!["bc1-1".into(), "bc2-1".into(), "bc3-1".into()],
        counts: [((0, 1), 2.0), ((0, 0), 1.0), ((2, 0), 3.0)]
            .into_iter()
            .collect(),
        records: 6,
        molecules: 6,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    };
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("raw_feature_bc_matrix.h5");
    write_h5(&matrix, &path, "GRCh38").unwrap();

    let file = hdf5::File::open(path).unwrap();
    let group = file.group("matrix").unwrap();
    assert_eq!(
        group.dataset("shape").unwrap().read_raw::<u64>().unwrap(),
        [2, 3]
    );
    assert_eq!(
        group.dataset("data").unwrap().read_raw::<u32>().unwrap(),
        [1, 2, 3]
    );
    assert_eq!(
        group.dataset("indices").unwrap().read_raw::<u32>().unwrap(),
        [0, 1, 0]
    );
    assert_eq!(
        group.dataset("indptr").unwrap().read_raw::<u32>().unwrap(),
        [0, 2, 2, 3]
    );
    let features = group.group("features").unwrap();
    let genomes = features
        .dataset("genome")
        .unwrap()
        .read_raw::<VarLenAscii>()
        .unwrap();
    assert_eq!(
        genomes.iter().map(|x| x.as_str()).collect::<Vec<_>>(),
        ["GRCh38", "GRCh38"]
    );
}

fn two_barcode_bam(dir: &Path) -> PathBuf {
    let path = dir.join("in.bam");
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 1000),
    );
    let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(10)]);
    for (i, cb) in ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1"]
        .into_iter()
        .enumerate()
    {
        let mut r = bam::Record::new();
        r.set(b"read", Some(&cigar), b"ACGTACGTAC", &[30; 10]);
        r.set_tid(0);
        r.set_pos(i as i64);
        r.push_aux(b"CB", Aux::String(cb)).unwrap();
        r.push_aux(b"UB", Aux::String("UMI0000AAAA")).unwrap();
        r.push_aux(b"GX", Aux::String("GENE1")).unwrap();
        r.push_aux(b"xf", Aux::I32(25)).unwrap();
        writer.write(&r).unwrap();
    }
    drop(writer);
    path
}

#[test]
fn unlisted_barcode_errors_by_default() {
    let dir = tempfile::tempdir().unwrap();
    let bam = two_barcode_bam(dir.path());
    let barcodes_file = dir.path().join("only_one.txt");
    fs::write(&barcodes_file, "AAAAAAAAAAAAAAAA\n").unwrap();

    let err = read_matrix(&bam, None, None, Some(&barcodes_file), 1, false, None).unwrap_err();
    assert!(err.to_string().contains("absent from barcode template"));
}

#[test]
fn skip_unlisted_barcodes_drops_them_instead_of_erroring() {
    let dir = tempfile::tempdir().unwrap();
    let bam = two_barcode_bam(dir.path());
    let barcodes_file = dir.path().join("only_one.txt");
    fs::write(&barcodes_file, "AAAAAAAAAAAAAAAA\n").unwrap();

    let matrix = read_matrix(&bam, None, None, Some(&barcodes_file), 1, true, None).unwrap();
    assert_eq!(matrix.barcodes, ["AAAAAAAAAAAAAAAA-1"]);
    assert_eq!(matrix.barcodes_skipped, 1);
    assert_eq!(matrix.molecules, 1);
    assert_eq!(matrix.counts.len(), 1);
}

fn shared_gem_two_probe_barcodes_bam(dir: &Path) -> PathBuf {
    let path = dir.join("in.bam");
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 1000),
    );
    let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(10)]);
    for (i, (pb, gene)) in [("D-A01", "GENE1"), ("D-A02", "GENE2")]
        .into_iter()
        .enumerate()
    {
        let mut r = bam::Record::new();
        r.set(b"read", Some(&cigar), b"ACGTACGTAC", &[30; 10]);
        r.set_tid(0);
        r.set_pos(i as i64);
        r.push_aux(b"CB", Aux::String("AAAAAAAAAAAAAAAA-1"))
            .unwrap();
        r.push_aux(b"UB", Aux::String("UMI0000AAAA")).unwrap();
        r.push_aux(b"GX", Aux::String(gene)).unwrap();
        r.push_aux(b"PB", Aux::String(pb)).unwrap();
        r.push_aux(b"xf", Aux::I32(25)).unwrap();
        writer.write(&r).unwrap();
    }
    drop(writer);
    path
}

#[test]
fn without_a_probe_barcode_allowlist_a_shared_gem_counts_both_samples_molecules() {
    let dir = tempfile::tempdir().unwrap();
    let bam = shared_gem_two_probe_barcodes_bam(dir.path());
    let matrix = read_matrix(&bam, None, None, None, 1, false, None).unwrap();
    assert_eq!(matrix.molecules, 2);
    assert_eq!(matrix.counts.len(), 2);
}

#[test]
fn probe_barcode_allowlist_keeps_only_that_samples_own_molecules() {
    let dir = tempfile::tempdir().unwrap();
    let bam = shared_gem_two_probe_barcodes_bam(dir.path());
    let allowlist: std::collections::HashSet<String> = ["D-A01".to_string()].into();

    let matrix = read_matrix(&bam, None, None, None, 1, false, Some(&allowlist)).unwrap();
    assert_eq!(matrix.molecules, 1);
    assert_eq!(matrix.probe_barcodes_skipped, 1);
    assert_eq!(matrix.features.len(), 1);
    assert_eq!(matrix.features[0].id, "GENE1");
}

#[test]
fn reads_fractional_zw_tags_and_emits_real_matrix_market_and_f32_h5() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("fractional.bam");
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "chr1")
            .push_tag(b"LN", 1000),
    );
    let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
    let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(10)]);

    for (i, (gene, weight)) in [("GENE1", 0.6_f32), ("GENE2", 0.4_f32)].into_iter().enumerate() {
        let mut r = bam::Record::new();
        r.set(b"read", Some(&cigar), b"ACGTACGTAC", &[30; 10]);
        r.set_tid(0);
        r.set_pos(i as i64);
        r.push_aux(b"CB", Aux::String("AAAAAAAAAAAAAAAA-1")).unwrap();
        r.push_aux(b"UB", Aux::String("UMI0000AAAA")).unwrap();
        r.push_aux(b"GX", Aux::String(gene)).unwrap();
        r.push_aux(b"ZG", Aux::String("fractional")).unwrap();
        r.push_aux(b"ZW", Aux::Float(weight)).unwrap();
        r.push_aux(b"xf", Aux::I32(25)).unwrap();
        writer.write(&r).unwrap();
    }
    drop(writer);

    let matrix = read_matrix(&path, None, None, None, 1, false, None).unwrap();
    assert_eq!(matrix.molecules, 2);
    assert_eq!(matrix.counts.len(), 2);
    let c0 = matrix.counts[&(0, 0)];
    let c1 = matrix.counts[&(0, 1)];
    assert!((c0 - 0.6).abs() < 1e-4);
    assert!((c1 - 0.4).abs() < 1e-4);

    let mex_dir = dir.path().join("mex");
    write_matrix(&matrix, &mex_dir, false).unwrap();
    let mtx_text = fs::read_to_string(mex_dir.join("matrix.mtx")).unwrap();
    assert!(mtx_text.contains("%%MatrixMarket matrix coordinate real general"));

    let h5_path = dir.path().join("raw_feature_bc_matrix.h5");
    write_h5(&matrix, &h5_path, "GRCh38").unwrap();
    let file = hdf5::File::open(h5_path).unwrap();
    let group = file.group("matrix").unwrap();
    let data_vals = group.dataset("data").unwrap().read_raw::<f32>().unwrap();
    assert_eq!(data_vals.len(), 2);
    assert!((data_vals[0] - 0.6).abs() < 1e-4);
    assert!((data_vals[1] - 0.4).abs() < 1e-4);
}

