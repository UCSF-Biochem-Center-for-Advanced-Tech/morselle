use super::converter::run;
use super::types::{canonical_id, metrics_path, Cli};
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

#[test]
fn normalizes_common_mate_suffixes() {
    assert_eq!(canonical_id(b"read/1"), b"read");
    assert_eq!(canonical_id(b"read/2 extra"), b"read");
    assert_eq!(canonical_id(b"read"), b"read");
}

fn base_cli(root: &Path, output: PathBuf) -> Cli {
    Cli {
        r1: vec![root.join("tests/fixtures/illumina_R1.fastq")],
        r2: vec![root.join("tests/fixtures/illumina_R2.fastq")],
        whitelist: root.join("tests/fixtures/illumina_whitelist.txt"),
        output,
        metrics: None,
        gem_group: 1,
        threads: 1,
        split: 1,
        overlap: false,
        overlap_min_overlap: 20,
        overlap_kmer: 8,
        overlap_kmer_offset: 1,
        overlap_check_lengths: 20,
        overlap_max_mismatch: 100,
        overlap_max_mismatch_density: 0.25,
        no_emit_yc: false,
        min_cdna_len: 50,
        adapter_max_mm: 2,
        trim_tso: false,
        trim_r2_polya: false,
        r2_polya_min_run: 12,
        failed_output: None,
    }
}

fn read_gz(path: &Path) -> String {
    let mut text = String::new();
    MultiGzDecoder::new(File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

#[test]
fn converts_exact_and_hamming1_and_gates_failure() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let output = temp.path().join("tagged.fastq.gz");
    run(base_cli(&root, output.clone())).unwrap();

    let text = read_gz(&output);
    assert!(text.starts_with("@read_exact\tCR:Z:"));
    assert!(!text.contains("read_exact/2"));
    assert_eq!(text.matches("\n+").count(), 2);
    assert!(text.contains("CR:Z:AAACCCAAGAAACACA"));
    assert!(text.contains("CB:Z:AAACCCAAGAAACACT-1"));
    assert!(!text.contains("read_rejected"));
    assert!(!text.contains("YC:Z:"));
}

#[test]
fn converts_multiple_pairs_to_numbered_shards() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let output = temp.path().join("tagged.fastq.gz");
    let mut cli = base_cli(&root, output.clone());
    cli.r1 = vec![cli.r1[0].clone(), cli.r1[0].clone()];
    cli.r2 = vec![cli.r2[0].clone(), cli.r2[0].clone()];
    cli.threads = 2;
    run(cli).unwrap();

    assert!(!output.exists());
    assert!(temp.path().join("tagged.shard-001.fastq.gz").is_file());
    assert!(temp.path().join("tagged.shard-002.fastq.gz").is_file());
    let metrics = std::fs::read_to_string(metrics_path(&output)).unwrap();
    assert!(metrics.contains("tagged.shard-001.fastq.gz"));
    assert!(metrics.contains("\"records_read\": 6"));
}

fn write_overlapping_fixture(dir: &Path) -> (PathBuf, PathBuf, Vec<u8>) {
    const CB: &[u8] = b"AAACCCAAGAAACACT";
    const UMI: &[u8] = b"ACGTACGTACGT";
    let mut core: Vec<u8> =
        b"CTGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATCC"
            .to_vec();
    core.truncate(82);
    *core.last_mut().unwrap() = b'C';
    assert_eq!(core.len(), 82);

    let mut molecule = Vec::new();
    molecule.extend_from_slice(CB);
    molecule.extend_from_slice(UMI);
    molecule.extend(std::iter::repeat_n(b'T', 30));
    molecule.extend_from_slice(&crate::anchors::revcomp(&core));
    assert_eq!(molecule.len(), 140);

    let mut r1 = molecule.clone();
    r1.truncate(100);
    let mut r2 = crate::anchors::revcomp(&molecule);
    r2.extend_from_slice(&crate::anchors::revcomp(crate::anchors::R1_PRIMER));
    r2.truncate(150);

    let write = |name: &str, seq: &[u8]| -> PathBuf {
        let path = dir.join(name);
        let text = format!(
            "@frag 1:N:0:1\n{}\n+\n{}\n",
            String::from_utf8_lossy(seq),
            "I".repeat(seq.len())
        );
        std::fs::write(&path, text).unwrap();
        path
    };
    (write("ov_R1.fastq", &r1), write("ov_R2.fastq", &r2), core)
}

#[test]
fn overlap_merges_and_emits_the_full_cdna_with_a_yc_tag() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let (r1, r2, core) = write_overlapping_fixture(temp.path());
    let output = temp.path().join("tagged.fastq.gz");

    let mut cli = base_cli(&root, output.clone());
    cli.r1 = vec![r1];
    cli.r2 = vec![r2];
    cli.overlap = true;
    run(cli).unwrap();

    let text = read_gz(&output);
    let mut lines = text.lines();
    let header = lines.next().unwrap();
    let seq = lines.next().unwrap();

    assert!(header.contains("CR:Z:AAACCCAAGAAACACT"));
    assert!(header.contains("CB:Z:AAACCCAAGAAACACT-1"));
    assert_eq!(seq.as_bytes(), &core[..]);
    assert!(header.ends_with("\tYC:Z:40-42+"), "header was {header}");

    let metrics = std::fs::read_to_string(metrics_path(&output)).unwrap();
    assert!(metrics.contains("\"overlap_enabled\": true"));
    assert!(metrics.contains("\"overlap_short_sins\": 1"));
    assert!(metrics.contains("\"overlap_medium_mins\": 0"));
    assert!(metrics.contains("\"overlap_long_lins\": 0"));
    assert!(metrics.contains("\"82\": 1"), "insert histogram: {metrics}");
}

#[test]
fn failed_output_accounts_for_every_input_read() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let output = temp.path().join("tagged.fastq.gz");
    let failed = temp.path().join("failed");
    let mut cli = base_cli(&root, output.clone());
    cli.failed_output = Some(failed);
    run(cli).unwrap();

    let failed_r1 = temp.path().join("failed_R1.fastq.gz");
    let failed_r2 = temp.path().join("failed_R2.fastq.gz");
    assert!(failed_r1.is_file() && failed_r2.is_file());

    let r1_text = read_gz(&failed_r1);
    let r2_text = read_gz(&failed_r2);
    assert!(r1_text.contains("@read_rejected"), "R1 was: {r1_text}");
    assert!(r2_text.contains("@read_rejected"), "R2 was: {r2_text}");
    assert!(r1_text.contains("XF:Z:cb_unmatched"));
    assert_eq!(r1_text.lines().count(), 4);
    assert_eq!(r2_text.lines().count(), 4);

    let metrics: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(metrics_path(&output)).unwrap()).unwrap();
    let read = metrics["records_read"].as_u64().unwrap();
    let emitted = metrics["records_emitted"].as_u64().unwrap();
    let failed_n = metrics["records_failed"].as_u64().unwrap();
    assert_eq!(read, emitted + failed_n, "every read accounted for");
    assert_eq!(failed_n, 1);
    assert_eq!(metrics["failed_by_reason"]["cb_unmatched"], 1);
}

#[test]
fn failed_reads_are_emitted_verbatim() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let failed = temp.path().join("failed");
    let mut cli = base_cli(&root, temp.path().join("tagged.fastq.gz"));
    cli.failed_output = Some(failed);
    run(cli).unwrap();

    let source =
        std::fs::read_to_string(root.join("tests/fixtures/illumina_R1.fastq")).unwrap();
    let rejected: Vec<&str> = source
        .lines()
        .skip_while(|l| !l.starts_with("@read_rejected"))
        .take(4)
        .collect();
    let emitted: Vec<String> = read_gz(&temp.path().join("failed_R1.fastq.gz"))
        .lines()
        .map(str::to_string)
        .collect();
    assert_eq!(emitted[1], rejected[1], "sequence unchanged");
    assert_eq!(emitted[3], rejected[3], "qualities unchanged");
    assert!(
        emitted[0].starts_with(rejected[0]),
        "header preserved, tag appended: {} vs {}",
        emitted[0],
        rejected[0]
    );
}

#[test]
fn both_paths_emit_the_same_orientation() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let (r1, r2, core) = write_overlapping_fixture(temp.path());

    let merged_out = temp.path().join("merged.fastq.gz");
    let mut cli = base_cli(&root, merged_out.clone());
    cli.r1 = vec![r1.clone()];
    cli.r2 = vec![r2.clone()];
    cli.overlap = true;
    run(cli).unwrap();

    let lins_out = temp.path().join("lins.fastq.gz");
    let mut cli = base_cli(&root, lins_out.clone());
    cli.r1 = vec![r1];
    cli.r2 = vec![r2];
    cli.overlap = true;
    cli.overlap_min_overlap = 200;
    run(cli).unwrap();

    let merged_seq = read_gz(&merged_out).lines().nth(1).unwrap().to_string();
    let lins_seq = read_gz(&lins_out).lines().nth(1).unwrap().to_string();

    assert_eq!(
        merged_seq.as_bytes(),
        &core[..],
        "merged is sense-strand cDNA"
    );
    assert!(
        lins_seq.starts_with(&merged_seq),
        "the lins path emits R2 verbatim, which must begin with the same bases"
    );
    assert!(lins_seq.len() > merged_seq.len());
}

#[test]
fn a_non_overlapping_pair_falls_back_to_the_standard_path() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let plain = temp.path().join("plain.fastq.gz");
    run(base_cli(&root, plain.clone())).unwrap();

    let overlapped = temp.path().join("overlapped.fastq.gz");
    let mut cli = base_cli(&root, overlapped.clone());
    cli.overlap = true;
    run(cli).unwrap();

    assert_eq!(read_gz(&plain), read_gz(&overlapped));
    let metrics = std::fs::read_to_string(metrics_path(&overlapped)).unwrap();
    assert!(metrics.contains("\"overlap_short_sins\": 0"));
    assert!(metrics.contains("\"overlap_medium_mins\": 0"));
    assert!(metrics.contains("\"overlap_long_lins\": 3"));
}

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

    let plain_out = temp.path().join("plain.fastq.gz");
    let mut cli = base_cli(&root, plain_out.clone());
    cli.r1 = vec![r1.clone()];
    cli.r2 = vec![r2.clone()];
    run(cli).unwrap();
    let plain_seq = read_gz(&plain_out).lines().nth(1).unwrap().to_string();
    assert!(plain_seq.len() > cdna.len(), "untrimmed R2 carries the tail");

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
