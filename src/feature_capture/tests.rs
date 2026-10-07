use rust_htslib::bam::{self, record::Aux, Read as _};
use std::fs::{self, File};
use std::io::Write as _;

use super::matcher::FeatureMatcher;
use super::pipeline::{establish_features_tsv, parse_tag, run};
use super::types::{Cli, FeatureReference};

fn write_csv(text: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::NamedTempFile::new().unwrap();
    f.write_all(text.as_bytes()).unwrap();
    f
}

/// Real CAT06 shape: two CRISPR guides sharing one anchor.
const REAL_SHAPE: &str = "\
id,name,read,pattern,sequence,feature_type
CD45,CD45,R2,TTCCAGCTTAGCTCTTAAAC(BC),CTGTTACCCTAAGAACAAAC,CRISPR Guide Capture
inactive,inactive,R2,TTTCTGGACTTCAGTGGACG(BC),GGTGCACTCTGGAGAAGTCC,CRISPR Guide Capture
";

#[test]
fn loads_real_shape() {
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    assert_eq!(r.features.len(), 2);
    assert_eq!(r.features[0].id, "CD45");
    assert_eq!(r.features[0].anchor, b"TTCCAGCTTAGCTCTTAAAC");
    assert_eq!(r.features[0].barcode, b"CTGTTACCCTAAGAACAAAC");
    assert_eq!(r.features[0].feature_type, "CRISPR Guide Capture");
}

#[test]
fn rejects_barcode_then_anchor_pattern_explicitly() {
    let f = write_csv(
        "id,name,read,pattern,sequence,feature_type\nA,A,R2,(BC)ACGTACGTACGTACGTACGT,ACGTACGTACGTACGTACGT,CRISPR Guide Capture\n",
    );
    let err = FeatureReference::load(f.path()).unwrap_err();
    assert!(err.to_string().contains("(BC)<anchor>"));
}

#[test]
fn rejects_pattern_without_bc_suffix() {
    let f = write_csv(
        "id,name,read,pattern,sequence,feature_type\nA,A,R2,ACGTACGTACGTACGTACGT,ACGTACGTACGTACGTACGT,CRISPR Guide Capture\n",
    );
    assert!(FeatureReference::load(f.path()).is_err());
}

fn read_with_anchor_capture(prefix_len: usize, anchor: &[u8], capture: &[u8]) -> Vec<u8> {
    let mut seq = vec![b'N'; prefix_len];
    seq.extend_from_slice(anchor);
    seq.extend_from_slice(capture);
    seq.extend_from_slice(b"NNNNNNNNNN");
    seq
}

#[test]
fn locates_exact_match_at_confirmed_offset() {
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
    let read = read_with_anchor_capture(10, b"TTCCAGCTTAGCTCTTAAAC", b"CTGTTACCCTAAGAACAAAC");
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn tolerates_mismatches_within_budget() {
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
    let mut anchor = b"TTCCAGCTTAGCTCTTAAAC".to_vec();
    anchor[3] = b'G'; // 1 mismatch, within anchor_max_mm=2
    let mut capture = b"CTGTTACCCTAAGAACAAAC".to_vec();
    capture[0] = b'A'; // 1 mismatch, within capture_max_mm=1
    let read = read_with_anchor_capture(10, &anchor, &capture);
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn no_match_on_unrelated_sequence() {
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
    let read = vec![b'A'; 60];
    assert_eq!(matcher.locate(&read), None);
}

#[test]
fn distinguishes_two_guides_sharing_one_anchor() {
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
    let read = read_with_anchor_capture(10, b"TTTCTGGACTTCAGTGGACG", b"GGTGCACTCTGGAGAAGTCC");
    assert_eq!(matcher.locate(&read), Some(1)); // "inactive"
}

#[test]
fn recovers_a_single_deletion_in_the_anchor_span() {
    // Mirrors probe_set::tests::recovers_a_single_deletion_in_the_probe_span:
    // the read is missing a base the anchor has -- del index should still
    // resolve it, at the *shorter* (anchor_len - 1) window, and the
    // capture barcode should still be found starting right after it.
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
    let anchor = b"TTCCAGCTTAGCTCTTAAAC";
    let mut deleted = anchor[..10].to_vec();
    deleted.extend_from_slice(&anchor[11..]); // drop the base at position 10
    assert_eq!(deleted.len(), anchor.len() - 1);
    let read = read_with_anchor_capture(10, &deleted, b"CTGTTACCCTAAGAACAAAC");
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn recovers_a_single_insertion_in_the_anchor_span() {
    // Mirrors probe_set::tests::recovers_a_single_insertion_in_the_probe_span.
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
    let anchor = b"TTCCAGCTTAGCTCTTAAAC";
    let mut inserted = anchor[..10].to_vec();
    inserted.push(b'T'); // extra base not present in the real anchor
    inserted.extend_from_slice(&anchor[10..]);
    assert_eq!(inserted.len(), anchor.len() + 1);
    let read = read_with_anchor_capture(10, &inserted, b"CTGTTACCCTAAGAACAAAC");
    assert_eq!(matcher.locate(&read), Some(0));
}

#[test]
fn indel_recovery_still_distinguishes_two_guides() {
    // A deletion in one guide's anchor shouldn't cause it to be confused
    // with the other guide sharing the same anchor sequence.
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);
    let anchor = b"TTTCTGGACTTCAGTGGACG";
    let mut deleted = anchor[..8].to_vec();
    deleted.extend_from_slice(&anchor[9..]);
    let read = read_with_anchor_capture(10, &deleted, b"GGTGCACTCTGGAGAAGTCC");
    assert_eq!(matcher.locate(&read), Some(1)); // "inactive", not "CD45"
}

/// Proves the generalization: a mixed panel (a CRISPR-flavored feature
/// and an "Antibody Capture"-flavored one with a *different* anchor)
/// resolves correctly -- the matcher doesn't special-case feature_type
/// or assume one universal anchor.
#[test]
fn handles_a_mixed_panel_with_different_anchors_per_feature_type() {
    let csv = "\
id,name,read,pattern,sequence,feature_type
CD45,CD45,R2,TTCCAGCTTAGCTCTTAAAC(BC),CTGTTACCCTAAGAACAAAC,CRISPR Guide Capture
CD3,CD3-APC,R2,AAAAAAAAAAAAAAAAAAAA(BC),TCGATCGATCGATCGATCGA,Antibody Capture
";
    let f = write_csv(csv);
    let r = FeatureReference::load(f.path()).unwrap();
    let matcher = FeatureMatcher::build(&r.features, 2, 1, 30);

    let crispr_read =
        read_with_anchor_capture(10, b"TTCCAGCTTAGCTCTTAAAC", b"CTGTTACCCTAAGAACAAAC");
    assert_eq!(matcher.locate(&crispr_read), Some(0));

    let ab_read = read_with_anchor_capture(10, b"AAAAAAAAAAAAAAAAAAAA", b"TCGATCGATCGATCGATCGA");
    assert_eq!(matcher.locate(&ab_read), Some(1));
}

#[test]
fn features_tsv_create_then_verify_round_trips() {
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("features.tsv");
    establish_features_tsv(&out, &r.features).unwrap();
    let text = fs::read_to_string(&out).unwrap();
    assert_eq!(
        text,
        "CD45\tCD45\tCRISPR Guide Capture\ninactive\tinactive\tCRISPR Guide Capture\n"
    );
    // idempotent re-run against the same reference
    establish_features_tsv(&out, &r.features).unwrap();
}

#[test]
fn features_tsv_rejects_mismatched_existing_file() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("features.tsv");
    fs::write(&out, "stale\tSTALE\tCRISPR Guide Capture\n").unwrap();
    let f = write_csv(REAL_SHAPE);
    let r = FeatureReference::load(f.path()).unwrap();
    assert!(establish_features_tsv(&out, &r.features).is_err());
}

#[test]
fn parse_tag_extracts_value_and_stops_at_tab() {
    let comment = "CR:Z:ACGT\tCB:Z:ACGT-1\tPB:Z:D-A01";
    assert_eq!(parse_tag(comment, "CB:Z:"), Some("ACGT-1"));
    assert_eq!(parse_tag(comment, "PB:Z:"), Some("D-A01"));
    assert_eq!(parse_tag(comment, "XX:Z:"), None);
}

#[test]
fn run_end_to_end_writes_mapped_and_unmapped_records_with_carried_tags() {
    let feature_csv = write_csv(REAL_SHAPE);
    let dir = tempfile::tempdir().unwrap();
    let fastq_path = dir.path().join("r2.fastq");
    let matched_seq =
        read_with_anchor_capture(10, b"TTCCAGCTTAGCTCTTAAAC", b"CTGTTACCCTAAGAACAAAC");
    let unmatched_seq = vec![b'A'; 60];
    let mut fastq = File::create(&fastq_path).unwrap();
    writeln!(
        fastq,
        "@read1\tCR:Z:ACGTACGTACGTACGT\tCB:Z:ACGTACGTACGTACGT-1\tUR:Z:AAAACCCCGGGG\tUB:Z:AAAACCCCGGGG\tPB:Z:D-A01"
    ).unwrap();
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

    let out_bam = dir.path().join("guides.bam");
    let features_out = dir.path().join("features.tsv");
    let metrics_path = dir.path().join("metrics.json");
    let cli = Cli {
        r2_fastq: fastq_path,
        feature_reference: feature_csv.path().to_path_buf(),
        output: out_bam.clone(),
        features_out: Some(features_out.clone()),
        metrics: Some(metrics_path.clone()),
        anchor_max_mm: 2,
        capture_max_mm: 1,
        anchor_window: 30,
    };
    run(&cli).unwrap();

    assert!(features_out.exists());
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
    assert_eq!(mapped.aux(b"GX").unwrap(), Aux::String("CD45"));
    assert_eq!(
        mapped.aux(b"CB").unwrap(),
        Aux::String("ACGTACGTACGTACGT-1")
    );
    assert_eq!(mapped.aux(b"PB").unwrap(), Aux::String("D-A01"));

    let unmapped = &records[1];
    assert!(unmapped.is_unmapped());
    assert_eq!(unmapped.qname(), b"read2");
    assert!(unmapped.aux(b"GX").is_err());
    assert_eq!(unmapped.aux(b"CB").unwrap(), Aux::String("TTTTGGGGCCCC-1"));
    assert!(unmapped.aux(b"PB").is_err()); // PB optional, absent here
}
