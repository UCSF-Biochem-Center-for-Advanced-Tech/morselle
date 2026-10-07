use super::alignment::{overlap, overlap_in};
use super::types::{build_index, InsertClass, OverlapCfg, Scratch};
use crate::anchors::revcomp;
use crate::yc::decode;

/// A 60 bp fragment with no long homopolymers and unique 8-mers.
const FRAGMENT: &[u8] = b"ACGTTGCAAGCTTAGCCATGGTACCAGTCTGAACGTACGGTCATGGCGTTACAGCTAGCA";
/// The 5' end of the universal Illumina 3' read-through adapter.
const ADAPTER: &[u8] = b"AGATCGGAAG";

/// R1 quality Q40, R2 quality Q20, so R1 wins every disagreement.
const Q40: u8 = b'I';
const Q20: u8 = b'5';

fn quals(base: u8, len: usize) -> Vec<u8> {
    vec![base; len]
}

/// Build a mate pair from the R1 insert and the desired `rc(R2)`.
fn mates(r1: &[u8], rc2: &[u8]) -> (Vec<u8>, Vec<u8>) {
    (r1.to_vec(), revcomp(rc2))
}

#[test]
fn fixture_is_the_length_the_tests_assume() {
    assert_eq!(FRAGMENT.len(), 60);
    assert_eq!(ADAPTER.len(), 10);
}

#[test]
fn merges_a_medium_insert_into_the_whole_fragment() {
    let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q40, r2.len()),
        &OverlapCfg::default(),
    )
    .expect("mates overlap by 20 bp");

    assert_eq!(merged.seq, FRAGMENT);
    assert_eq!(merged.class, InsertClass::Medium);
    assert_eq!(merged.overlap_len, 20);
    assert_eq!(merged.left_overhang, 20);
    assert_eq!(merged.right_overhang, 20);
    assert_eq!(merged.discordant, 0);
    assert_eq!(merged.yc.unwrap().to_string(), "20+20-20");
}

#[test]
fn merges_a_short_insert_and_drops_read_through_adapter() {
    let insert = &FRAGMENT[..30];
    let mut r1 = insert.to_vec();
    r1.extend_from_slice(ADAPTER);
    let mut rc2 = revcomp(ADAPTER);
    rc2.extend_from_slice(insert);
    let r2 = revcomp(&rc2);

    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q40, r2.len()),
        &OverlapCfg::default(),
    )
    .expect("mates overlap across the whole insert");

    assert_eq!(merged.seq, insert);
    assert_eq!(merged.class, InsertClass::Short);
    assert_eq!(merged.left_overhang, 0);
    assert_eq!(merged.right_overhang, 0);
    assert_eq!(merged.yc.unwrap().to_string(), "+30+");
}

#[test]
fn handles_a_mate_longer_than_r1() {
    let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[10..60]);
    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q40, r2.len()),
        &OverlapCfg::default(),
    )
    .expect("30 bp overlap");

    assert_eq!(merged.seq, FRAGMENT);
    assert_eq!(merged.overlap_len, 30);
    assert_eq!(merged.left_overhang, 10);
    assert_eq!(merged.right_overhang, 20);
    assert_eq!(merged.yc.unwrap().to_string(), "10+30-20");
}

#[test]
fn extends_from_the_correct_offset_when_rc_r2_starts_first() {
    let mut rc2 = revcomp(ADAPTER);
    rc2.extend_from_slice(&FRAGMENT[..50]);
    let r1 = FRAGMENT[..30].to_vec();
    let r2 = revcomp(&rc2);

    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q40, r2.len()),
        &OverlapCfg::default(),
    )
    .expect("30 bp overlap");

    assert_eq!(merged.seq, &FRAGMENT[..50]);
    assert_eq!(merged.seq.len(), 50);
    assert_eq!(merged.left_overhang, 0);
    assert_eq!(merged.right_overhang, 20);
    assert_eq!(merged.yc.unwrap().to_string(), "+30-20");
}

#[test]
fn agreeing_bases_sum_quality_up_to_the_cap() {
    let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q20, r2.len()),
        &OverlapCfg::default(),
    )
    .unwrap();

    assert!(merged.qual[..20].iter().all(|&q| q == Q40));
    assert!(merged.qual[20..40].iter().all(|&q| q == Q40));
    assert!(merged.qual[40..].iter().all(|&q| q == Q20));
}

#[test]
fn disagreeing_bases_take_the_better_call_and_the_quality_delta() {
    let mut rc2 = FRAGMENT[20..60].to_vec();
    rc2[10] = if rc2[10] == b'A' { b'C' } else { b'A' };
    let r1 = FRAGMENT[..40].to_vec();
    let r2 = revcomp(&rc2);

    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q20, r2.len()),
        &OverlapCfg::default(),
    )
    .unwrap();

    assert_eq!(merged.seq, FRAGMENT);
    assert_eq!(merged.discordant, 1);
    assert_eq!(merged.qual[30], Q20);
    assert_eq!(merged.qual[29], Q40);
    assert_eq!(merged.yc.unwrap().count_discordant(), 1);
}

#[test]
fn unrelated_reads_do_not_merge() {
    let r1 = b"ACGTTGCAAGCTTAGCCATGGTACCAGTCTGA".to_vec();
    let r2 = b"TTTGGCCATATCAGCATCAGGATTCCGGATAC".to_vec();
    assert!(overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q40, r2.len()),
        &OverlapCfg::default()
    )
    .is_none());
}

#[test]
fn overlap_below_min_overlap_is_refused() {
    let (r1, r2) = mates(&FRAGMENT[..30], &FRAGMENT[20..50]);
    let r1q = quals(Q40, r1.len());
    let r2q = quals(Q40, r2.len());

    assert!(overlap(&r1, &r1q, &r2, &r2q, &OverlapCfg::default()).is_none());

    let permissive = OverlapCfg {
        min_overlap: 10,
        ..OverlapCfg::default()
    };
    let merged = overlap(&r1, &r1q, &r2, &r2q, &permissive).unwrap();
    assert_eq!(merged.seq, &FRAGMENT[..50]);
    assert_eq!(merged.overlap_len, 10);
}

#[test]
fn overlap_over_the_mismatch_budget_is_refused() {
    let mut rc2 = FRAGMENT[20..60].to_vec();
    for base in &mut rc2[8..14] {
        *base = if *base == b'A' { b'C' } else { b'A' };
    }
    let r1 = FRAGMENT[..40].to_vec();
    let r2 = revcomp(&rc2);

    assert!(overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q40, r2.len()),
        &OverlapCfg::default()
    )
    .is_none());
}

#[test]
fn mismatched_sequence_and_quality_lengths_are_refused() {
    let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
    assert!(overlap(
        &r1,
        &quals(Q40, 5),
        &r2,
        &quals(Q40, r2.len()),
        &OverlapCfg::default()
    )
    .is_none());
}

#[test]
fn short_mates_do_not_underflow() {
    let cfg = OverlapCfg {
        min_overlap: 5,
        ..OverlapCfg::default()
    };
    for len in 0..30usize {
        let r2: Vec<u8> = FRAGMENT[20..20 + len].to_vec();
        let r1 = FRAGMENT[..40].to_vec();
        let _ = overlap(&r1, &quals(Q40, r1.len()), &r2, &quals(Q40, len), &cfg);
    }
    for len in 0..12usize {
        let r1: Vec<u8> = FRAGMENT[..len].to_vec();
        let r2 = revcomp(&FRAGMENT[20..60]);
        let _ = overlap(&r1, &quals(Q40, len), &r2, &quals(Q40, r2.len()), &cfg);
    }
}

#[test]
fn homopolymer_seeds_are_masked() {
    let mut index = Vec::new();
    build_index(b"AAAAAAAAAAAAAAAAAAAA", 8, 1, &mut index);
    assert!(index.is_empty(), "a pure homopolymer must index nothing");

    let mut index = Vec::new();
    build_index(b"TTTTTTTTTTTTACGTTGCAAGCTTAGC", 8, 1, &mut index);
    assert!(index.iter().all(|&(_, pos)| pos >= 5));
}

#[test]
fn seeds_containing_non_acgt_are_masked() {
    let mut index = Vec::new();
    build_index(b"ACGTTGCAAGCTTAGC", 8, 1, &mut index);
    let clean = index.len();

    let mut index = Vec::new();
    build_index(b"ACGTNGCAAGCTTAGC", 8, 1, &mut index);
    assert_eq!(index.len(), clean - 5);
    assert!(index.iter().all(|&(_, pos)| pos > 4));
}

#[test]
fn kmer_offset_strides_the_index() {
    let mut every = Vec::new();
    build_index(&FRAGMENT[..40], 8, 1, &mut every);
    let mut strided = Vec::new();
    build_index(&FRAGMENT[..40], 8, 4, &mut strided);
    assert_eq!(every.len(), 33);
    assert!(strided.iter().all(|&(_, pos)| pos % 4 == 0));
    assert!(strided.len() < every.len());
}

#[test]
fn the_tag_reconstructs_both_mates_from_the_merge() {
    let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q20, r2.len()),
        &OverlapCfg::default(),
    )
    .unwrap();

    let tag = merged.yc.as_ref().unwrap();
    assert_eq!(tag.sequence_length(), merged.seq.len() as u64);

    let decoded = decode(tag, &merged.seq, &merged.qual).unwrap();
    assert_eq!(decoded.r1_seq, r1, "R1 recovered from the merge");
    assert_eq!(
        decoded.r2_seq,
        revcomp(&r2),
        "rc(R2) recovered from the merge"
    );
}

#[test]
fn the_tag_records_disagreements_losslessly() {
    let mut rc2 = FRAGMENT[20..60].to_vec();
    rc2[10] = if rc2[10] == b'A' { b'C' } else { b'A' };
    let r1 = FRAGMENT[..40].to_vec();
    let r2 = revcomp(&rc2);

    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &quals(Q20, r2.len()),
        &OverlapCfg::default(),
    )
    .unwrap();

    let decoded = decode(merged.yc.as_ref().unwrap(), &merged.seq, &merged.qual).unwrap();
    assert_eq!(decoded.r1_seq, r1);
    assert_eq!(decoded.r2_seq, rc2);
}

#[test]
fn n_bases_are_coerced_and_counted() {
    let mut rc2 = FRAGMENT[20..60].to_vec();
    rc2[5] = b'N';
    let r1 = FRAGMENT[..40].to_vec();
    let r2 = revcomp(&rc2);
    let mut r2q = quals(Q40, r2.len());
    r2q[r2.len() - 1 - 5] = b'#';

    let merged = overlap(
        &r1,
        &quals(Q40, r1.len()),
        &r2,
        &r2q,
        &OverlapCfg::default(),
    )
    .expect("one N does not break the seed or the budget");

    assert_eq!(merged.seq, FRAGMENT);
    assert_eq!(merged.yc_stats.n_columns_coerced, 1);
    assert_eq!(merged.yc.unwrap().to_string(), "20+20-20");
}

#[test]
fn scratch_buffers_are_reusable_across_pairs() {
    let mut scratch = Scratch::default();
    let cfg = OverlapCfg::default();
    for _ in 0..3 {
        let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
        let merged = overlap_in(
            &mut scratch,
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q40, r2.len()),
            &cfg,
        )
        .unwrap();
        assert_eq!(merged.seq, FRAGMENT);
    }
}
