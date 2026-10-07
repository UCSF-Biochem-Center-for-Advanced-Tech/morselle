use std::fs::{self, File};
use std::io::Write as _;
use rust_htslib::bam::{self, record::Aux, Read as _};

use super::aligner::ProbeMatcher;
use super::pipeline::{establish_genes_manifest, run};
use super::types::{Cli, ProbeSetReference, DEFAULT_PROBE_SEQ_OFFSET};

fn write_csv(text: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(text.as_bytes()).unwrap();
    f
}

/// Real header + shape (25bp probes here, not the real 50bp, purely to
/// keep test literals short -- ProbeSetReference doesn't hardcode a length).
const REAL_SHAPE: &str = "\
#probe_set_file_format=3.0
#panel_name=test
gene_id,probe_seq,probe_id,included,region,gene_name
ENSG1,AAAACCCCGGGGTTTTACGTACGTA,ENSG1|GENE1|abc,TRUE,spliced,GENE1
ENSG1,CCCCGGGGTTTTAAAAACGTACGT,ENSG1|GENE1|def,FALSE,unspliced,GENE1
ENSG2,TTTTGGGGCCCCAAAATGCATGCAT,ENSG2|GENE2|ghi,TRUE,unspliced,GENE2
";

#[test]
fn loads_only_included_true_rows() {
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    assert_eq!(r.probes.len(), 2);  // the FALSE row is excluded
    assert_eq!(r.probes[0].gene_id, "ENSG1");
    assert_eq!(r.probes[0].region, "spliced");
    assert_eq!(r.probes[1].gene_id, "ENSG2");
}

#[test]
fn rejects_mixed_probe_lengths() {
    let csv = "\
gene_id,probe_seq,probe_id,included,region,gene_name
ENSG1,AAAACCCCGGGGTTTTACGTACGTA,p1,TRUE,spliced,GENE1
ENSG2,SHORT,p2,TRUE,spliced,GENE2
";
    let f = write_csv(csv);
    assert!(ProbeSetReference::load(f.path()).is_err());
}

fn read_with_probe_at(offset: usize, probe: &[u8]) -> Vec<u8> {
    let mut seq = vec![b'N'; offset];
    seq.extend_from_slice(probe);
    seq.extend_from_slice(b"NNNNNNNNNN");
    seq
}

#[test]
fn locates_exact_match_at_confirmed_offset() {
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    let matcher = ProbeMatcher::build(&r.probes, 1, 15);
    let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, b"AAAACCCCGGGGTTTTACGTACGTA");
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn tolerates_a_mismatch_within_budget() {
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    let matcher = ProbeMatcher::build(&r.probes, 1, 15);
    let mut probe = b"AAAACCCCGGGGTTTTACGTACGTA".to_vec();
    probe[3] = b'G';  // 1 mismatch, within max_mm=1
    let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, &probe);
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn windowed_search_absorbs_offset_jitter() {
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    let matcher = ProbeMatcher::build(&r.probes, 1, 15);
    // 5bp late: still within the +-15 window. (Not "early" -- offset 0
    // has no room to go earlier; residual extract()-side anchor jitter
    // only ever pushes the probe insert a few bp *later* than 0.)
    let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET + 5, b"AAAACCCCGGGGTTTTACGTACGTA");
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn distinguishes_two_genes_probes() {
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    let matcher = ProbeMatcher::build(&r.probes, 1, 15);
    let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, b"TTTTGGGGCCCCAAAATGCATGCAT");
    assert_eq!(matcher.locate(&read), Some(1));  // GENE2's probe
}

#[test]
fn recovers_a_single_deletion_in_the_probe_span() {
    // Real SBX error mode (2026-09-18 diagnostic): the read is missing a
    // base the probe has -- del_index should still resolve it, at the
    // *shorter* (probe_len - 1) window.
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    let matcher = ProbeMatcher::build(&r.probes, 1, 15);
    let probe = b"AAAACCCCGGGGTTTTACGTACGTA";
    let mut deleted = probe[..10].to_vec();
    deleted.extend_from_slice(&probe[11..]);  // drop the base at position 10
    assert_eq!(deleted.len(), probe.len() - 1);
    let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, &deleted);
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn recovers_a_single_insertion_in_the_probe_span() {
    // The rarer real shape (9.2% of recoverable indel-shaped misses vs
    // 85.2% deletion) but still real -- ins_index resolves it at the
    // *longer* (probe_len + 1) window.
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    let matcher = ProbeMatcher::build(&r.probes, 1, 15);
    let probe = b"AAAACCCCGGGGTTTTACGTACGTA";
    let mut inserted = probe[..10].to_vec();
    inserted.push(b'T');  // extra base not present in the real probe
    inserted.extend_from_slice(&probe[10..]);
    assert_eq!(inserted.len(), probe.len() + 1);
    let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, &inserted);
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn indel_recovery_still_distinguishes_two_genes() {
    // A deletion shouldn't blur two real, distinct probes together.
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    let matcher = ProbeMatcher::build(&r.probes, 1, 15);
    let probe2 = b"TTTTGGGGCCCCAAAATGCATGCAT";
    let mut deleted = probe2[..5].to_vec();
    deleted.extend_from_slice(&probe2[6..]);
    let read = read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, &deleted);
    assert_eq!(matcher.locate(&read), Some(1));  // GENE2's probe, not GENE1's
}

#[test]
fn no_match_on_unrelated_sequence() {
    let f = write_csv(REAL_SHAPE);
    let r = ProbeSetReference::load(f.path()).unwrap();
    let matcher = ProbeMatcher::build(&r.probes, 1, 15);
    let read = vec![b'A'; 100];
    assert_eq!(matcher.locate(&read), None);
}

#[test]
fn genes_manifest_create_then_verify_round_trips() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("genes.tsv");
    let genes = vec![("ENSG1", "GENE1"), ("ENSG2", "GENE2")];
    establish_genes_manifest(&out, &genes).unwrap();
    let text = fs::read_to_string(&out).unwrap();
    assert_eq!(
        text,
        "ENSG1\tGENE1\tGene Expression\nENSG2\tGENE2\tGene Expression\n"
    );
    establish_genes_manifest(&out, &genes).unwrap();  // idempotent
}

#[test]
fn genes_manifest_rejects_mismatched_existing_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("genes.tsv");
    fs::write(&out, "stale\tSTALE\tGene Expression\n").unwrap();
    assert!(establish_genes_manifest(&out, &[("ENSG1", "GENE1")]).is_err());
}

#[test]
fn run_end_to_end_writes_mapped_and_unmapped_records_with_carried_tags() {
    let probe_csv = write_csv(REAL_SHAPE);
    let dir = tempfile::tempdir().unwrap();
    let fastq_path = dir.path().join("r2.fastq");
    let matched_seq =
        read_with_probe_at(DEFAULT_PROBE_SEQ_OFFSET, b"AAAACCCCGGGGTTTTACGTACGTA");
    let unmatched_seq = vec![b'A'; 120];
    let mut fastq = File::create(&fastq_path).unwrap();
    writeln!(
        fastq,
        "@read1\tCR:Z:ACGTACGTACGTACGT\tCB:Z:ACGTACGTACGTACGT-1\tUR:Z:AAAACCCCGGGG\tUB:Z:AAAACCCCGGGG\tPB:Z:D-A01"
    )
    .unwrap();
    writeln!(fastq, "{}", String::from_utf8_lossy(&matched_seq)).unwrap();
    writeln!(fastq, "+").unwrap();
    writeln!(fastq, "{}", "I".repeat(matched_seq.len())).unwrap();
    writeln!(
        fastq,
        "@read2\tCR:Z:TTTTGGGGCCCC\tCB:Z:TTTTGGGGCCCC-1\tUR:Z:GGGGAAAACCCC\tUB:Z:GGGGAAAACCCC"
    )
    .unwrap();
    writeln!(fastq, "{}", String::from_utf8_lossy(&unmatched_seq)).unwrap();
    writeln!(fastq, "+").unwrap();
    writeln!(fastq, "{}", "I".repeat(unmatched_seq.len())).unwrap();
    drop(fastq);

    let out_bam = dir.path().join("probeset.bam");
    let genes_out = dir.path().join("genes.tsv");
    let metrics_path = dir.path().join("metrics.json");
    let cli = Cli {
        r2_fastq: fastq_path,
        probe_set: probe_csv.path().to_path_buf(),
        output: out_bam.clone(),
        genes_out: Some(genes_out.clone()),
        metrics: Some(metrics_path.clone()),
        max_mm: 1,
        window: 15,
        threads: 1,
    };
    run(&cli).unwrap();

    assert!(genes_out.exists());
    let metrics: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&metrics_path).unwrap()).unwrap();
    assert_eq!(metrics["reads"], 2);
    assert_eq!(metrics["matched"], 1);

    let mut reader = bam::Reader::from_path(&out_bam).unwrap();
    let records: Vec<_> = reader.records().map(|r| r.unwrap()).collect();
    assert_eq!(records.len(), 2);

    let mapped = &records[0];
    assert!(!mapped.is_unmapped());
    assert_eq!(mapped.qname(), b"read1");
    assert_eq!(mapped.aux(b"GX").unwrap(), Aux::String("ENSG1"));
    assert_eq!(mapped.aux(b"GN").unwrap(), Aux::String("GENE1"));
    assert_eq!(mapped.aux(b"PB").unwrap(), Aux::String("D-A01"));

    let unmapped = &records[1];
    assert!(unmapped.is_unmapped());
    assert!(unmapped.aux(b"GX").is_err());
    assert!(unmapped.aux(b"PB").is_err());  // PB optional, absent here
}
