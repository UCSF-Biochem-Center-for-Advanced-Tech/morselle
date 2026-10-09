//! `--failed-output` on the SBX arm: every input read is accounted for exactly
//! once, and the reads written to the side channel are the untouched originals.

use std::io::Read;
use std::path::PathBuf;

use clap::Parser;
use flate2::read::MultiGzDecoder;
use sbx10x::config::Cli;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn read_gz(path: &PathBuf) -> String {
    let mut text = String::new();
    MultiGzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

#[test]
fn extract_per_barcode_values_include_barcode_less_reads() {
    let temp = tempfile::tempdir().unwrap();
    let input = temp.path().join("three.fastq");
    let wanted = [
        "@instrument:date_group_0:chip:0:3:1z1G2Y",
        "@instrument:date_group_0:chip:0:3:2uAEXJ",
        "@instrument:date_group_0:chip:0:3:1zBptJ",
    ];
    let source = read_gz(&fixture("smoke.fastq.gz"));
    let lines: Vec<&str> = source.lines().collect();
    let selected = lines
        .chunks_exact(4)
        .filter(|record| wanted.contains(&record[0]))
        .flat_map(|record| record.iter().copied())
        .collect::<Vec<_>>()
        .join("\n")
        + "\n";
    assert_eq!(selected.lines().count(), 12);
    std::fs::write(&input, selected).unwrap();
    let out_dir = temp.path().join("out");
    let cli = Cli::parse_from([
        "sbxpress",
        "--fastqs",
        input.to_str().unwrap(),
        "--out-dir",
        out_dir.to_str().unwrap(),
        "--threads",
        "1",
        "--minimap2",
        "--barcode-metrics",
    ]);
    sbx10x::run(&cli).unwrap();
    let table = read_gz(&out_dir.join("metrics_summary.barcodes.tsv.gz"));
    assert!(table.contains("TGGCCATCACCTGCCT-1\t.\t1\t0\t0\t1\t1\t0\t0"));
    assert!(table.contains("GGTAACTTCTTGATGA-1\t.\t1\t0\t0\t1\t1\t0\t0"));
    assert!(table.contains("__no_barcode__\t.\t1\t0\t0\t0\t0\t0\t1"));
}

#[test]
fn failed_output_accounts_for_every_read() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let failed = temp.path().join("failed.fastq.gz");

    let cli = Cli::parse_from([
        "sbxpress",
        "--fastqs",
        fixture("smoke.fastq.gz").to_str().unwrap(),
        "--out-dir",
        out_dir.to_str().unwrap(),
        "--threads",
        "1",
        "--failed-output",
        failed.to_str().unwrap(),
        "--barcode-metrics",
    ]);
    sbx10x::run(&cli).unwrap();

    assert!(failed.is_file(), "side channel must be created");

    let metrics: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(out_dir.join("metrics_summary.json")).unwrap(),
    )
    .unwrap();
    let total = metrics["total_reads"].as_u64().unwrap();
    let emitted = metrics["pairs_emitted"].as_u64().unwrap();
    let failed_n = metrics["reads_failed"].as_u64().unwrap();
    let table = read_gz(&out_dir.join("metrics_summary.barcodes.tsv.gz"));
    let mut lines = table.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let rows: Vec<Vec<&str>> = lines.map(|line| line.split('\t').collect()).collect();
    let sum = |name: &str| {
        let index = header.iter().position(|&column| column == name).unwrap();
        rows.iter()
            .map(|row| row[index].parse::<u64>().unwrap())
            .sum::<u64>()
    };
    assert_eq!(sum("reads"), total);
    assert_eq!(sum("emitted"), emitted);
    for (column, field) in [
        ("class_full", "class_full"),
        ("class_partial", "class_partial"),
        ("class_unparseable", "class_unparseable"),
    ] {
        assert_eq!(sum(column), metrics[field].as_u64().unwrap());
    }
    assert!(total > 0, "fixture must contain reads");
    assert_eq!(
        total,
        emitted + failed_n,
        "every input read must be emitted or written to the side channel"
    );

    let text = read_gz(&failed);
    assert_eq!(
        text.lines().count() as u64,
        failed_n * 4,
        "one four-line record per failed read"
    );
    if failed_n > 0 {
        assert!(
            text.lines().next().unwrap().contains("XF:Z:"),
            "each failed read carries its reason"
        );
    }
}

#[test]
fn without_the_flag_no_side_channel_is_written() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let cli = Cli::parse_from([
        "sbxpress",
        "--fastqs",
        fixture("smoke.fastq.gz").to_str().unwrap(),
        "--out-dir",
        out_dir.to_str().unwrap(),
        "--threads",
        "1",
    ]);
    sbx10x::run(&cli).unwrap();

    let metrics: serde_json::Value = serde_json::from_str(
        &std::fs::read_to_string(out_dir.join("metrics_summary.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(metrics["reads_failed"], 0);
}

#[test]
fn barcode_table_is_written_only_with_barcode_metrics() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let cli = Cli::parse_from([
        "sbxpress",
        "--fastqs",
        fixture("smoke.fastq.gz").to_str().unwrap(),
        "--out-dir",
        out_dir.to_str().unwrap(),
        "--threads",
        "1",
    ]);
    sbx10x::run(&cli).unwrap();
    assert!(out_dir.join("metrics_summary.json").is_file());
    assert!(!out_dir.join("metrics_summary.barcodes.tsv.gz").exists());
}
