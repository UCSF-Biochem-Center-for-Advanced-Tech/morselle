use std::collections::BTreeMap;
use rust_htslib::bam::{self, record::Aux};

use crate::probe_sample_sheet::{ProbeSample, ProbeSampleSheet};
use super::engine::*;
use super::types::*;

fn sheet() -> ProbeSampleSheet {
    ProbeSampleSheet {
        samples: vec![
            ProbeSample {
                id: "Pos".into(),
                probe_barcode_ids: vec!["D-A01".into(), "D-A02".into()],
            },
            ProbeSample {
                id: "RZ211".into(),
                probe_barcode_ids: vec!["D-A03".into(), "D-A04".into()],
            },
            ProbeSample {
                id: "RZ213NR".into(),
                probe_barcode_ids: vec!["D-A07".into()],
            },
        ],
    }
}

fn tally(entries: &[(&str, &str, u64)]) -> BTreeMap<(String, String), u64> {
    entries
        .iter()
        .map(|&(cb, pb, n)| ((cb.to_string(), pb.to_string()), n))
        .collect()
}

#[test]
fn sums_across_a_samples_whole_declared_set() {
    let t = tally(&[("CB1", "D-A01", 30), ("CB1", "D-A02", 25)]);
    let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
    assert_eq!(out.len(), 1);
    assert_eq!(out[0].status, Status::Assigned);
    assert_eq!(out[0].sample_id.as_deref(), Some("Pos"));
    assert_eq!(out[0].molecules, 55);
    assert!((out[0].purity - 1.0).abs() < 1e-9);
}

#[test]
fn below_floor_is_unassigned() {
    let t = tally(&[("CB1", "D-A01", 3)]);
    let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
    assert_eq!(out[0].status, Status::Unassigned);
    assert_eq!(out[0].sample_id, None);
}

#[test]
fn two_samples_each_clearing_the_floor_both_get_assigned() {
    let t = tally(&[("CB1", "D-A01", 50), ("CB1", "D-A03", 40)]);
    let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
    assert_eq!(out.len(), 2);
    let pos = out
        .iter()
        .find(|a| a.sample_id.as_deref() == Some("Pos"))
        .unwrap();
    assert_eq!(pos.status, Status::Assigned);
    assert_eq!(pos.molecules, 50);
    assert_eq!(pos.occupancy, 2);
    assert!((pos.purity - 50.0 / 90.0).abs() < 1e-9);
    let rz = out
        .iter()
        .find(|a| a.sample_id.as_deref() == Some("RZ211"))
        .unwrap();
    assert_eq!(rz.status, Status::Assigned);
    assert_eq!(rz.molecules, 40);
    assert!((rz.purity - 40.0 / 90.0).abs() < 1e-9);
}

#[test]
fn a_little_bleed_from_a_second_sample_does_not_trigger_high_occupancy() {
    let t = tally(&[("CB1", "D-A01", 50), ("CB1", "D-A03", 2)]);
    let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
    assert_eq!(out[0].status, Status::Assigned);
    assert_eq!(out[0].sample_id.as_deref(), Some("Pos"));
    assert!((out[0].purity - 50.0 / 52.0).abs() < 1e-9);
}

#[test]
fn molecules_on_an_undeclared_probe_barcode_count_against_purity_but_not_toward_any_sample() {
    let t = tally(&[("CB1", "D-A01", 20), ("CB1", "A-F04", 15)]);
    let out = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
    assert_eq!(out[0].status, Status::Assigned);
    assert_eq!(out[0].molecules, 20);
    assert!((out[0].purity - 20.0 / 35.0).abs() < 1e-9);
}

#[test]
fn negbinom_quantile_matches_hand_computed_values() {
    assert_eq!(negbinom_quantile(2.0, 0.5, 0.999), 12);
    assert_eq!(negbinom_quantile(0.0, 0.5, 0.999), u64::MAX);
}

#[test]
fn high_occupancy_gem_tail_is_filtered_but_ordinary_multi_sample_gems_are_not() {
    let samples: Vec<ProbeSample> = (1..=40)
        .map(|i| ProbeSample {
            id: format!("S{i}"),
            probe_barcode_ids: vec![format!("P{i}")],
        })
        .collect();
    let big_sheet = ProbeSampleSheet { samples };

    let mut entries: Vec<(String, String, u64)> = Vec::new();
    for i in 1..=39 {
        entries.push((format!("NORMAL{i}"), format!("P{i}"), 20));
    }
    for i in 1..=40 {
        entries.push(("ARTIFACT".to_string(), format!("P{i}"), 20));
    }
    let t: BTreeMap<(String, String), u64> = entries
        .into_iter()
        .map(|(cb, pb, n)| ((cb, pb), n))
        .collect();

    let out = assign(&big_sheet, &t, &DemuxCfg { min_molecules: 10 });

    let normal_rows: Vec<&Assignment> = out
        .iter()
        .filter(|a| a.cell_barcode.starts_with("NORMAL"))
        .collect();
    assert_eq!(normal_rows.len(), 39);
    assert!(normal_rows.iter().all(|a| a.status == Status::Assigned));

    let artifact_rows: Vec<&Assignment> = out
        .iter()
        .filter(|a| a.cell_barcode == "ARTIFACT")
        .collect();
    assert_eq!(
        artifact_rows.len(),
        1,
        "the artifact GEM must be filtered as one row, not 40 Assigned rows"
    );
    assert_eq!(artifact_rows[0].status, Status::HighOccupancyGem);
    assert_eq!(artifact_rows[0].sample_id, None);
    assert_eq!(artifact_rows[0].occupancy, 40);
}

#[test]
fn gem_overlap_reports_only_pairs_with_shared_cells() {
    let t = tally(&[
        ("CB1", "D-A01", 10),
        ("CB1", "D-A02", 10),
        ("CB2", "D-A01", 10),
        ("CB3", "D-A03", 10),
    ]);
    let declared: std::collections::HashSet<&str> = ["D-A01", "D-A02", "D-A03"].into();
    let rows = gem_overlap(&t, &declared);
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!(
        (r.barcode1.as_str(), r.barcode2.as_str()),
        ("D-A01", "D-A02")
    );
    assert_eq!(r.barcode1_gems, 2);
    assert_eq!(r.barcode2_gems, 1);
    assert_eq!(r.common_gems, 1);
    assert!((r.overlap - 1.0).abs() < 1e-9);
}

#[test]
fn gem_overlap_ignores_undeclared_noise_barcodes() {
    let t = tally(&[("CB1", "D-A01", 10), ("CB1", "A-F04", 1)]);
    let declared: std::collections::HashSet<&str> = ["D-A01", "D-A02"].into();
    let rows = gem_overlap(&t, &declared);
    assert!(rows.is_empty());
}

#[test]
fn bam_tally_gates_on_xf_bit_8_and_skips_missing_or_empty_pb() {
    let temp = tempfile::tempdir().unwrap();
    let bam_path = temp.path().join("in.bam");
    let mut header = bam::Header::new();
    header.push_record(
        bam::header::HeaderRecord::new(b"SQ")
            .push_tag(b"SN", "guide_contig")
            .push_tag(b"LN", 100),
    );
    let mut writer = bam::Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(10)]);

    let mut push = |cb: Option<&str>, pb: Option<&str>, xf: i32| {
        let mut r = bam::Record::new();
        r.set(b"read", Some(&cigar), b"ACGTACGTAC", &[30; 10]);
        r.set_tid(0);
        r.set_pos(0);
        if let Some(cb) = cb {
            r.push_aux(b"CB", Aux::String(cb)).unwrap();
        }
        if let Some(pb) = pb {
            r.push_aux(b"PB", Aux::String(pb)).unwrap();
        }
        r.push_aux(b"xf", Aux::I32(xf)).unwrap();
        writer.write(&r).unwrap();
    };
    push(Some("CB1-1"), Some("D-A01"), 25);
    push(Some("CB1-1"), Some("D-A01"), 17);
    push(Some("CB1-1"), Some(""), 25);
    push(Some("CB2-1"), None, 25);
    drop(writer);

    let tally = read_bam_tally(&bam_path).unwrap();
    assert_eq!(tally.len(), 1);
    assert_eq!(tally[&("CB1-1".to_string(), "D-A01".to_string())], 1);
}

#[test]
fn strip_gem_suffix_removes_trailing_digits_only() {
    assert_eq!(strip_gem_suffix("ACGTACGTACGTACGT-1"), "ACGTACGTACGTACGT");
    assert_eq!(strip_gem_suffix("ACGTACGTACGTACGT-12"), "ACGTACGTACGTACGT");
    assert_eq!(strip_gem_suffix("ACGTACGTACGTACGT"), "ACGTACGTACGTACGT");
    assert_eq!(strip_gem_suffix("D-A01"), "D-A01");
}

#[test]
fn per_sample_barcode_files_are_bare_acgt_matrix_cell_barcodes_ready() {
    let t = tally(&[
        ("AAACCCAAGAAACCCA-1", "D-A01", 50),
        ("TTTGGGTTCTTTGGGT-1", "D-A02", 50),
    ]);
    let assignments = assign(&sheet(), &t, &DemuxCfg { min_molecules: 10 });
    let dir = tempfile::tempdir().unwrap();
    write_per_sample_barcodes(&sheet(), &assignments, dir.path()).unwrap();

    let pos_barcodes = std::fs::read_to_string(dir.path().join("Pos.barcodes.txt")).unwrap();
    let lines: Vec<&str> = pos_barcodes.lines().collect();
    assert_eq!(lines, vec!["AAACCCAAGAAACCCA", "TTTGGGTTCTTTGGGT"]);
    for line in &lines {
        assert!(
            line.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')),
            "barcode {line} is not bare ACGT"
        );
    }
}
