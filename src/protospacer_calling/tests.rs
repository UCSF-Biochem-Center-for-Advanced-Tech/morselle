use std::collections::BTreeMap;
use std::fs;

use crate::cr_bam_to_mtx::{write_h5, Feature, Matrix};
use super::gmm::{decision_boundary, fit_threshold};
use super::pipeline::run;
use super::types::{Cli, Component};

const POS_CD45_HIST: &[(u64, u64)] = &[
    (0, 53432),
    (1, 16134),
    (2, 3827),
    (3, 805),
    (4, 223),
    (5, 87),
    (6, 31),
    (7, 22),
    (8, 8),
    (9, 3),
    (10, 4),
    (11, 3),
    (14, 1),
    (15, 1),
];

const RZ211_CD45_HIST: &[(u64, u64)] = &[
    (0, 9416),
    (1, 8714),
    (2, 6301),
    (3, 4504),
    (4, 3340),
    (5, 2517),
    (6, 2155),
    (7, 1686),
    (8, 1390),
    (9, 1210),
    (10, 979),
    (11, 857),
    (12, 740),
    (13, 688),
    (14, 596),
    (15, 492),
    (16, 416),
    (17, 421),
    (18, 343),
    (19, 300),
    (20, 263),
    (21, 227),
    (22, 212),
    (23, 218),
    (24, 177),
    (25, 163),
    (26, 133),
    (27, 132),
    (28, 122),
    (29, 110),
    (30, 91),
    (31, 83),
    (32, 79),
    (33, 93),
    (34, 75),
    (35, 57),
    (36, 56),
    (37, 75),
    (38, 48),
    (39, 60),
    (40, 41),
    (41, 40),
    (42, 34),
    (43, 32),
    (44, 28),
    (45, 33),
    (46, 33),
    (47, 28),
    (48, 25),
    (49, 17),
    (50, 20),
    (51, 20),
    (52, 28),
    (53, 16),
    (54, 15),
    (55, 17),
    (56, 17),
    (57, 6),
    (58, 11),
    (59, 19),
    (60, 11),
    (61, 20),
    (62, 8),
    (63, 12),
    (64, 7),
    (65, 18),
    (66, 9),
    (67, 9),
    (68, 10),
    (69, 7),
    (70, 7),
    (71, 8),
    (72, 6),
    (73, 3),
    (74, 5),
    (75, 3),
    (76, 3),
    (77, 4),
    (78, 2),
    (79, 5),
    (80, 2),
    (81, 2),
    (82, 2),
    (83, 3),
    (84, 5),
    (85, 3),
    (86, 4),
    (87, 4),
    (88, 3),
    (90, 5),
    (91, 1),
    (92, 1),
    (93, 1),
    (94, 2),
    (95, 1),
    (98, 1),
    (99, 3),
    (100, 3),
    (101, 4),
    (102, 1),
    (103, 1),
    (104, 2),
    (105, 1),
    (106, 2),
    (107, 1),
    (109, 3),
    (110, 2),
    (111, 2),
    (112, 1),
    (113, 2),
    (116, 2),
    (117, 1),
    (120, 2),
    (121, 1),
    (122, 1),
    (124, 1),
    (126, 1),
    (128, 1),
    (130, 1),
    (132, 1),
    (135, 1),
    (137, 2),
    (140, 1),
    (144, 2),
    (147, 1),
    (150, 1),
    (154, 2),
    (156, 2),
    (163, 1),
    (165, 1),
    (181, 1),
    (182, 1),
    (187, 1),
    (210, 1),
    (211, 1),
    (212, 1),
    (235, 1),
    (287, 1),
    (394, 1),
];

fn expand_hist(hist: &[(u64, u64)]) -> Vec<u64> {
    hist.iter()
        .flat_map(|&(value, n)| std::iter::repeat(value).take(n as usize))
        .collect()
}

#[test]
fn real_data_pos_cd45_threshold_matches_ground_truth() {
    let counts = expand_hist(POS_CD45_HIST);
    assert_eq!(fit_threshold(&counts, 3), Some(3));
}

#[test]
fn real_data_rz211_cd45_threshold_regression() {
    let counts = expand_hist(RZ211_CD45_HIST);
    assert_eq!(fit_threshold(&counts, 3), Some(3));
}

#[test]
fn all_zero_guide_is_omitted_not_defaulted() {
    assert_eq!(fit_threshold(&[0, 0, 0, 0], 3), None);
}

#[test]
fn floor_binds_when_fit_is_lower() {
    let mut counts = vec![0u64; 100];
    counts.extend(vec![1u64; 5]);
    assert_eq!(fit_threshold(&counts, 3), Some(3));
}

#[test]
fn clearly_bimodal_signal_classifies_both_populations_correctly() {
    let mut counts = vec![0u64; 200];
    counts.extend(vec![1u64; 50]);
    counts.extend(vec![50u64; 100]);
    counts.extend(vec![52u64; 50]);
    let t = fit_threshold(&counts, 3).unwrap();
    assert!(
        t > 1,
        "ambient value 1 must not clear the threshold, got threshold {t}"
    );
    assert!(
        t <= 50,
        "signal value 50 must clear the threshold, got threshold {t}"
    );
}

#[test]
fn decision_boundary_is_symmetric_for_equal_variance_equal_weight_components() {
    let low = Component {
        weight: 0.5,
        mean: 1.0,
        std: 1.0,
    };
    let high = Component {
        weight: 0.5,
        mean: 5.0,
        std: 1.0,
    };
    let b = decision_boundary(&low, &high);
    assert!((b - 3.0).abs() < 1e-3, "expected midpoint 3.0, got {b}");
}

fn tiny_matrix(barcodes: Vec<&str>, raw_counts: BTreeMap<(usize, usize), u64>) -> Matrix {
    let molecules = raw_counts.values().sum();
    let counts: BTreeMap<(usize, usize), f64> = raw_counts
        .into_iter()
        .map(|(k, v)| (k, v as f64))
        .collect();
    Matrix {
        features: vec![
            Feature {
                id: "CD45".into(),
                name: "CD45".into(),
            },
            Feature {
                id: "inactive".into(),
                name: "inactive".into(),
            },
        ],
        barcodes: barcodes.into_iter().map(String::from).collect(),
        counts,
        records: 0,
        molecules,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    }
}

#[test]
fn run_end_to_end_matches_expected_csv_and_json_shape() {
    let mut counts = BTreeMap::new();
    counts.insert((0, 0), 5);
    counts.insert((1, 0), 1);
    let matrix = tiny_matrix(vec!["CB1-1", "CB2-1", "CB3-1"], counts);

    let dir = tempfile::tempdir().unwrap();
    let h5 = dir.path().join("filtered_feature_bc_matrix.h5");
    write_h5(&matrix, &h5, "GRCh38").unwrap();

    let out = dir.path().join("outs");
    let cli = Cli {
        filtered_matrix: h5,
        output: out.clone(),
        min_umi: 3,
    };
    run(&cli).unwrap();

    let per_cell = fs::read_to_string(out.join("protospacer_calls_per_cell.csv")).unwrap();
    assert_eq!(
        per_cell,
        "cell_barcode,num_features,feature_call,num_umis\nCB1-1,1,CD45,5\n"
    );

    let thresholds = fs::read_to_string(out.join("protospacer_umi_thresholds.csv")).unwrap();
    assert_eq!(thresholds, "Protospacer,UMI threshold\nCD45,3\n");

    let json: serde_json::Value = serde_json::from_str(
        &fs::read_to_string(out.join("protospacer_umi_thresholds.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(json, serde_json::json!({"CD45": 3}));

    let summary = fs::read_to_string(out.join("protospacer_calls_summary.csv")).unwrap();
    let lines: Vec<&str> = summary.lines().collect();
    assert_eq!(
        lines[0],
        "feature_call,num_cells,pct_cells,median_umis,stddev_umis"
    );
    assert!(lines.contains(&"No guide molecules,1,0.3333333333333333,None,None"));
    assert!(lines.contains(&"No confident call,2,0.6666666666666666,None,None"));
    assert!(lines.contains(&"1 protospacer assigned,1,0.3333333333333333,None,None"));
    assert!(lines.contains(&"More than 1 protospacer assigned,0,0,None,None"));
    assert!(lines.contains(&"CD45,1,0.3333333333333333,5,0"));
}
