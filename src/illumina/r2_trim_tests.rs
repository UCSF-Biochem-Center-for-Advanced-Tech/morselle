//! R2 trimming in the non-merged path: poly-A read-through (`--trim-r2-polya`)
//! and a leading TSO (`--trim-tso`). Moved out of tests.rs to keep it under
//! 500 lines.
use super::converter::run;
use super::tests::{base_cli, read_gz, table_sums};
use super::types::metrics_path;
use std::path::{Path, PathBuf};

/// A short-R1 (28 bp CB+UMI only, no overlap possible) pair whose R2 reads
/// past a `cdna_len`-bp insert into a 15 bp poly-A tail and 10 bp of junk
/// beyond it — the shape a short fragment takes on a standard 10x-style
/// read layout, where nothing merges because R1 never reaches the cDNA.
fn write_r2_readthrough_fixture(dir: &Path, cdna_len: usize) -> (PathBuf, PathBuf, Vec<u8>) {
    const CB: &[u8] = b"AAACCCAAGAAACACT";
    const UMI: &[u8] = b"ACGTACGTACGT";
    let mut r1 = CB.to_vec();
    r1.extend_from_slice(UMI);
    assert_eq!(r1.len(), 28);

    let mut cdna: Vec<u8> = b"CTGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGC"
        .iter()
        .cycle()
        .take(cdna_len)
        .copied()
        .collect();
    // End on C: the pattern can land on A at an arbitrary cdna_len, which
    // would abut the poly-A block appended below and shift the run's true
    // start one base earlier than intended.
    if let Some(last) = cdna.last_mut() {
        *last = b'C';
    }
    let mut r2 = cdna.clone();
    r2.extend(std::iter::repeat_n(b'A', 15));
    r2.extend_from_slice(b"GGCGTAGCTA");

    let write = |name: &str, seq: &[u8]| -> PathBuf {
        let path = dir.join(name);
        let text = format!(
            "@short_frag 1:N:0:1\n{}\n+\n{}\n",
            String::from_utf8_lossy(seq),
            "I".repeat(seq.len())
        );
        std::fs::write(&path, text).unwrap();
        path
    };
    (write("rt_R1.fastq", &r1), write("rt_R2.fastq", &r2), cdna)
}

#[test]
fn trim_r2_polya_trims_readthrough_when_cdna_survives_min_len() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let (r1, r2, cdna) = write_r2_readthrough_fixture(temp.path(), 60);

    // Without the flag: R2 is emitted verbatim, poly-A and junk included.
    let plain_out = temp.path().join("plain.fastq.gz");
    let mut cli = base_cli(&root, plain_out.clone());
    cli.r1 = vec![r1.clone()];
    cli.r2 = vec![r2.clone()];
    run(cli).unwrap();
    let plain_seq = read_gz(&plain_out).lines().nth(1).unwrap().to_string();
    assert!(
        plain_seq.len() > cdna.len(),
        "untrimmed R2 carries the tail"
    );

    // With the flag: R2 is cut back to just the cDNA before the poly-A run.
    let trimmed_out = temp.path().join("trimmed.fastq.gz");
    let mut cli = base_cli(&root, trimmed_out.clone());
    cli.r1 = vec![r1];
    cli.r2 = vec![r2];
    cli.trim_r2_polya = true;
    run(cli).unwrap();
    let trimmed_seq = read_gz(&trimmed_out).lines().nth(1).unwrap().to_string();
    assert_eq!(trimmed_seq.as_bytes(), &cdna[..]);

    let metrics = std::fs::read_to_string(metrics_path(&trimmed_out)).unwrap();
    assert!(metrics.contains("\"r2_polya_trimmed\": 1"), "{metrics}");
    assert!(metrics.contains("\"r2_polya_discarded\": 0"), "{metrics}");
}

#[test]
fn trim_r2_polya_discards_when_cdna_left_is_too_short() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    // 20 bp of cDNA before the poly-A run, below the default 50 bp min.
    let (r1, r2, _cdna) = write_r2_readthrough_fixture(temp.path(), 20);

    let output = temp.path().join("tagged.fastq.gz");
    let failed = temp.path().join("failed");
    let mut cli = base_cli(&root, output.clone());
    cli.r1 = vec![r1];
    cli.r2 = vec![r2];
    cli.trim_r2_polya = true;
    cli.failed_output = Some(failed);
    run(cli).unwrap();

    let metrics = std::fs::read_to_string(metrics_path(&output)).unwrap();
    assert!(metrics.contains("\"records_emitted\": 0"), "{metrics}");
    assert!(metrics.contains("\"r2_polya_discarded\": 1"), "{metrics}");
    assert!(metrics.contains("\"r2_polya_short\": 1"), "{metrics}");
}

/// A non-merged pair whose R2 starts with the whole TSO, then `cdna_len` bp
/// of cDNA: a 3' GEX fragment that reached the transcript's 5' end.
fn write_r2_tso_fixture(dir: &Path, cdna_len: usize) -> (PathBuf, PathBuf, Vec<u8>) {
    let r1 = b"AAACCCAAGAAACACTACGTACGTACGT".to_vec();
    let cdna: Vec<u8> = b"CTGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGC"
        .iter()
        .cycle()
        .take(cdna_len)
        .copied()
        .collect();
    let mut r2 = crate::anchors::TSO.to_vec();
    r2.extend_from_slice(&cdna);
    let write = |name: &str, seq: &[u8]| -> PathBuf {
        let path = dir.join(name);
        let text = format!(
            "@tso_frag 1:N:0:1\n{}\n+\n{}\n",
            String::from_utf8_lossy(seq),
            "I".repeat(seq.len())
        );
        std::fs::write(&path, text).unwrap();
        path
    };
    (write("tso_R1.fastq", &r1), write("tso_R2.fastq", &r2), cdna)
}

#[test]
fn trim_tso_trims_the_start_of_a_non_merged_r2() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let (r1, r2, cdna) = write_r2_tso_fixture(temp.path(), 60);

    // Off by default: R2 is emitted with the TSO.
    let plain_out = temp.path().join("plain.fastq.gz");
    let mut cli = base_cli(&root, plain_out.clone());
    cli.r1 = vec![r1.clone()];
    cli.r2 = vec![r2.clone()];
    run(cli).unwrap();
    let plain_seq = read_gz(&plain_out).lines().nth(1).unwrap().to_string();
    assert!(plain_seq.as_bytes().starts_with(crate::anchors::TSO));

    let trimmed_out = temp.path().join("trimmed.fastq.gz");
    let mut cli = base_cli(&root, trimmed_out.clone());
    cli.r1 = vec![r1];
    cli.r2 = vec![r2];
    cli.trim_tso = true;
    run(cli).unwrap();
    let text = read_gz(&trimmed_out);
    let mut lines = text.lines();
    lines.next();
    assert_eq!(lines.next().unwrap().as_bytes(), &cdna[..]);
    lines.next();
    assert_eq!(
        lines.next().unwrap().len(),
        cdna.len(),
        "quality trimmed with the sequence"
    );

    let metrics = std::fs::read_to_string(metrics_path(&trimmed_out)).unwrap();
    assert!(metrics.contains("\"r2_tso_trimmed\": 1"), "{metrics}");
    assert!(
        metrics.contains("\"r2_tso_bases_trimmed\": 30"),
        "{metrics}"
    );
    assert!(metrics.contains("\"r2_tso_discarded\": 0"), "{metrics}");
}

#[test]
fn trim_tso_discards_when_cdna_left_is_too_short() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    // 20 bp of cDNA after the TSO, below the default 50 bp min.
    let (r1, r2, _cdna) = write_r2_tso_fixture(temp.path(), 20);

    let output = temp.path().join("tagged.fastq.gz");
    let mut cli = base_cli(&root, output.clone());
    cli.r1 = vec![r1];
    cli.r2 = vec![r2];
    cli.trim_tso = true;
    cli.failed_output = Some(temp.path().join("failed"));
    run(cli).unwrap();

    let metrics = std::fs::read_to_string(metrics_path(&output)).unwrap();
    assert!(metrics.contains("\"records_emitted\": 0"), "{metrics}");
    assert!(metrics.contains("\"r2_tso_discarded\": 1"), "{metrics}");
    assert!(metrics.contains("\"r2_tso_short\": 1"), "{metrics}");
    let table = read_gz(&crate::barcode_metrics::barcode_path(&metrics_path(
        &output,
    )));
    assert_eq!(table_sums(&table)["failed_r2_tso_short"], 1, "{table}");
}
