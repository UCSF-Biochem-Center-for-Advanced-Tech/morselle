//! Low-complexity cDNA at Illumina conversion (crate::complexity): flagged
//! and counted by default, failed with `--drop-low-complexity`.
use super::converter::run;
use super::tests::{base_cli, read_gz, table_sums};
use super::types::metrics_path;
use std::path::{Path, PathBuf};

/// One pair with a dead (poly-G) read 2 and one ordinary pair, same cell.
fn write_pairs(dir: &Path) -> (PathBuf, PathBuf) {
    let r1 = "AAACCCAAGAAACACTACGTACGTACGT";
    let normal = "GATTACACCATGGCTTGAACGTTAGCCATGCAAGTCCGATACGGATTTACAGG";
    let dead = "G".repeat(60);
    let fastq = |reads: &[(&str, &str)]| {
        reads
            .iter()
            .map(|(name, seq)| format!("@{name}\n{seq}\n+\n{}\n", "I".repeat(seq.len())))
            .collect::<String>()
    };
    let r1_path = dir.join("lc_R1.fastq");
    let r2_path = dir.join("lc_R2.fastq");
    std::fs::write(&r1_path, fastq(&[("dead/1", r1), ("normal/1", r1)])).unwrap();
    std::fs::write(&r2_path, fastq(&[("dead/2", &dead), ("normal/2", normal)])).unwrap();
    (r1_path, r2_path)
}

fn metrics(output: &Path) -> serde_json::Value {
    serde_json::from_str(&std::fs::read_to_string(metrics_path(output)).unwrap()).unwrap()
}

#[test]
fn low_complexity_read_is_tagged_and_counted_by_default() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let output = temp.path().join("tagged.fastq.gz");
    let (r1, r2) = write_pairs(temp.path());
    let mut cli = base_cli(&root, output.clone());
    cli.r1 = vec![r1];
    cli.r2 = vec![r2];
    run(cli).unwrap();

    let text = read_gz(&output);
    let header = |name: &str| {
        text.lines()
            .find(|l| l.starts_with(name))
            .unwrap()
            .to_owned()
    };
    assert!(header("@dead\t").ends_with("\tYL:i:1"));
    assert!(!header("@normal\t").contains("YL:"));
    let metrics = metrics(&output);
    assert_eq!(metrics["records_emitted"], 2);
    assert_eq!(metrics["low_complexity"], 1);
    assert_eq!(metrics["low_complexity_by_base"]["G"], 1);
    assert_eq!(metrics["low_complexity_by_base"]["A"], 0);
    assert_eq!(metrics["low_complexity_base_fraction"], 0.8);
    assert_eq!(metrics["drop_low_complexity"], false);
    let table = read_gz(&crate::barcode_metrics::barcode_path(&metrics_path(
        &output,
    )));
    let sums = table_sums(&table);
    assert_eq!(sums["low_complexity"], 1);
    assert_eq!(sums["failed_low_complexity"], 0);
    assert_eq!(sums["emitted"], 2);
}

#[test]
fn drop_low_complexity_fails_the_pair_with_its_reason() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let output = temp.path().join("tagged.fastq.gz");
    let failed = temp.path().join("failed");
    let (r1, r2) = write_pairs(temp.path());
    let mut cli = base_cli(&root, output.clone());
    cli.r1 = vec![r1];
    cli.r2 = vec![r2];
    cli.failed_output = Some(failed.clone());
    cli.complexity.drop_low_complexity = true;
    run(cli).unwrap();

    let text = read_gz(&output);
    assert!(!text.contains("@dead"));
    assert!(text.contains("@normal\t"));
    let metrics = metrics(&output);
    assert_eq!(metrics["records_read"], 2);
    assert_eq!(metrics["records_emitted"], 1);
    assert_eq!(metrics["low_complexity"], 1);
    assert_eq!(metrics["failed_by_reason"]["low_complexity"], 1);
    let failed_r2 = read_gz(&temp.path().join("failed_R2.fastq.gz"));
    assert!(failed_r2.starts_with("@dead/2\tXF:Z:low_complexity\n"));
    let table = read_gz(&crate::barcode_metrics::barcode_path(&metrics_path(
        &output,
    )));
    let sums = table_sums(&table);
    assert_eq!(sums["reads"], 2);
    assert_eq!(sums["emitted"], 1);
    assert_eq!(sums["low_complexity"], 1);
    assert_eq!(sums["failed_low_complexity"], 1);
}
