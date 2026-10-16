//! Convert standard paired Illumina 10x GEX FASTQs into the tagged R2-only
//! FASTQ contract consumed by `minimap2 -y` in the sbxpress workflow.
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

/// Clap's `default_value_t` does not feed `Default`, so mirror the defaults
/// here; tests flip only what they are exercising.
pub(super) fn base_cli(root: &Path, output: PathBuf) -> Cli {
    Cli {
        r1: vec![root.join("tests/fixtures/illumina_R1.fastq")],
        r2: vec![root.join("tests/fixtures/illumina_R2.fastq")],
        whitelist: root.join("tests/fixtures/illumina_whitelist.txt"),
        output,
        metrics: None,
        barcode_metrics: true,
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
        complexity: Default::default(),
    }
}

pub(super) fn read_gz(path: &Path) -> String {
    let mut text = String::new();
    MultiGzDecoder::new(File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    text
}

pub(super) fn table_sums(table: &str) -> std::collections::BTreeMap<String, u64> {
    let mut lines = table.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let mut sums = std::collections::BTreeMap::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        for (&column, value) in header[2..].iter().zip(&fields[2..]) {
            *sums.entry(column.to_owned()).or_insert(0) += value.parse::<u64>().unwrap();
        }
    }
    sums
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
    // Overlap is off by default, so nothing about the output may change.
    assert!(!text.contains("YC:Z:"));
    let table_path = crate::barcode_metrics::barcode_path(&metrics_path(&output));
    let table = read_gz(&table_path);
    assert!(table.contains("AAACCCAAGAAACACT-1\t.\t2\t1\t1\t2\t0\t0\t0\t0\t0\t0"));
    assert!(table.contains("__unresolved_barcode__\t.\t1\t0\t0\t0\t0\t0\t0\t0\t0\t1"));
    let sums = table_sums(&table);
    let metrics: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(metrics_path(&output)).unwrap()).unwrap();
    for (column, field) in [
        ("reads", "records_read"),
        ("barcode_exact", "barcode_exact"),
        ("barcode_corrected", "barcode_hamming1_corrected"),
        ("emitted", "records_emitted"),
    ] {
        assert_eq!(sums[column], metrics[field].as_u64().unwrap());
    }
    for (reason, count) in metrics["failed_by_reason"].as_object().unwrap() {
        assert_eq!(sums[&format!("failed_{reason}")], count.as_u64().unwrap());
    }
}

#[test]
fn converts_multiple_pairs_to_numbered_shards() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let output = temp.path().join("tagged.fastq.gz");
    let mut cli = base_cli(&root, output.clone());
    cli.r1 = vec![cli.r1[0].clone(), cli.r1[0].clone()];
    cli.r2 = vec![cli.r2[0].clone(), cli.r2[0].clone()];
    cli.threads = 8;
    cli.gem_group = 2;
    run(cli).unwrap();

    assert!(!output.exists());
    assert!(temp.path().join("tagged.shard-001.fastq.gz").is_file());
    assert!(temp.path().join("tagged.shard-002.fastq.gz").is_file());
    let metrics = std::fs::read_to_string(metrics_path(&output)).unwrap();
    assert!(metrics.contains("tagged.shard-001.fastq.gz"));
    assert!(metrics.contains("\"records_read\": 6"));
    let table = read_gz(&crate::barcode_metrics::barcode_path(&metrics_path(
        &output,
    )));
    assert!(table.contains("AAACCCAAGAAACACT-2\t.\t4\t2\t2\t4"));
    assert!(table.contains("__unresolved_barcode__\t.\t2"));
    let serial_output = temp.path().join("serial.fastq.gz");
    let mut serial = base_cli(&root, serial_output.clone());
    serial.r1 = vec![serial.r1[0].clone(), serial.r1[0].clone()];
    serial.r2 = vec![serial.r2[0].clone(), serial.r2[0].clone()];
    serial.gem_group = 2;
    run(serial).unwrap();
    let serial_table = read_gz(&crate::barcode_metrics::barcode_path(&metrics_path(
        &serial_output,
    )));
    assert_eq!(table, serial_table);
}

/// A long-R1 pair that merges: R1 100 bp, R2 150 bp reading past the CB into
/// the R1 primer, over an 82 bp cDNA.
fn write_overlapping_fixture(dir: &Path) -> (PathBuf, PathBuf, Vec<u8>) {
    const CB: &[u8] = b"AAACCCAAGAAACACT"; // on the whitelist fixture
    const UMI: &[u8] = b"ACGTACGTACGT";
    let mut core: Vec<u8> =
        b"CTGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATCC"
            .to_vec();
    core.truncate(82);
    // End on C so rc(core) does not start with T and blur the dT boundary.
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

    // The barcode is still tagged exactly as on the unmerged path...
    assert!(header.contains("CR:Z:AAACCCAAGAAACACT"));
    assert!(header.contains("CB:Z:AAACCCAAGAAACACT-1"));
    // ...and the emitted read is now the whole 82 bp cDNA, sense strand,
    // rather than the 150 bp R2 with its barcode read-through.
    assert_eq!(seq.as_bytes(), &core[..]);
    // 5'-most 40 bp came from R2 alone; the remaining 42 bp are duplex.
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
    // The stock fixture has 3 pairs, one of which carries an unmatchable
    // barcode. With --failed-output nothing may be silently dropped.
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
    // Both mates of the failed pair, with the reason and the original header.
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
    // The side channel must be re-runnable, so the mates come back byte for
    // byte apart from the appended reason tag.
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let temp = tempfile::tempdir().unwrap();
    let failed = temp.path().join("failed");
    let mut cli = base_cli(&root, temp.path().join("tagged.fastq.gz"));
    cli.failed_output = Some(failed);
    run(cli).unwrap();

    let source = std::fs::read_to_string(root.join("tests/fixtures/illumina_R1.fastq")).unwrap();
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
    // The same fragment down both routes: once merged, once forced onto the
    // lins path by an unreachable --overlap-min-overlap. The merged read must
    // be the sense-strand cDNA, and R2 (emitted verbatim by the lins path)
    // begins with exactly those bases — so the two paths agree.
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
    cli.overlap_min_overlap = 200; // longer than either read: cannot merge
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
    // The merge recovers the cDNA and drops the polyA/barcode read-through
    // that the lins path leaves on the end of R2.
    assert!(lins_seq.len() > merged_seq.len());
}

#[test]
fn a_non_overlapping_pair_falls_back_to_the_standard_path() {
    // The stock fixtures are far too short to overlap, so with --overlap on
    // they must all be counted as lins and emitted exactly as before.
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
