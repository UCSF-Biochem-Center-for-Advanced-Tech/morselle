//! `--index-indels`: one-indel i5/i7 matching and its sheet safety check.
use super::engine::*;
use super::indel::*;
use super::types::*;
use crate::anchors::{revcomp, P5, P7, R1_PRIMER, R2_PRIMER};
use crate::sample_sheet::SampleSheet;
use std::path::PathBuf;

fn sheet() -> SampleSheet {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/CE15204_bclconvert_input.csv");
    SampleSheet::load(&p).unwrap()
}

fn flat_sheet(rows: &[(&str, &str, &str)]) -> SampleSheet {
    let dir = tempfile::tempdir().unwrap();
    let p = dir.path().join("sheet.csv");
    let mut text = String::from("sample,i7,i5\n");
    for (id, i7, i5) in rows {
        text.push_str(&format!("{id},{i7},{i5}\n"));
    }
    std::fs::write(&p, text).unwrap();
    SampleSheet::load(&p).unwrap()
}

fn cfg() -> DemuxCfg {
    DemuxCfg {
        index_max_mm: 1,
        mode: None,
    }
}

fn idx(sh: &SampleSheet, id: &str) -> usize {
    sh.samples.iter().position(|s| s.id == id).unwrap()
}

/// Canonical read region: `P5 + CTACAC + i5_read + R1_PRIMER`, i5 nominal start.
fn i5_context(i5_read: &[u8]) -> (Vec<u8>, usize) {
    let mut seq = [P5, b"CTACAC"].concat();
    let start = seq.len();
    seq.extend_from_slice(i5_read);
    seq.extend_from_slice(R1_PRIMER);
    (seq, start)
}

/// Canonical read region: `rc(R2_PRIMER) + i7 + rc(P7)`, i7 nominal start.
fn i7_context(i7_read: &[u8]) -> (Vec<u8>, usize) {
    let mut seq = revcomp(R2_PRIMER);
    let start = seq.len();
    seq.extend_from_slice(i7_read);
    seq.extend_from_slice(&revcomp(P7));
    (seq, start)
}

fn with_deletion(key: &[u8], at: usize) -> Vec<u8> {
    let mut v = key.to_vec();
    v.remove(at);
    v
}

fn with_insertion(key: &[u8], at: usize, base: u8) -> Vec<u8> {
    let mut v = key.to_vec();
    v.insert(at, base);
    v
}

#[test]
fn edit_distance_counts_substitutions_and_indels() {
    assert_eq!(edit_distance(b"ACGTACGTAC", b"ACGTACGTAC"), 0);
    assert_eq!(edit_distance(b"ACGTACGTAC", b"ACGAACGTAC"), 1);
    assert_eq!(edit_distance(b"ACGTACGTAC", b"ACGACGTAC"), 1);
    assert_eq!(edit_distance(b"ACGTACGTAC", b"ACGTTACGTAC"), 1);
    // Hamming 10 (a one-base shift), edit distance 2.
    assert_eq!(edit_distance(b"ACGTACGTAC", b"CGTACGTACG"), 2);
}

/// Every one-indel variant of `key`.
fn indel_variants(key: &[u8]) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = (0..key.len()).map(|i| with_deletion(key, i)).collect();
    for i in 0..=key.len() {
        for &b in b"ACGT" {
            out.push(with_insertion(key, i, b));
        }
    }
    out
}

/// Pins the real CE15204 sheet's safety numbers -- if the fixture changes,
/// this must be re-verified, not silently pass with other numbers. Counting
/// shifts, the i5 set's minimum edit distance is 3 (Hamming says 4): V4_GEX
/// and M2_GEX read-oriented i5s are one-base shifts of each other, so one
/// V4 insertion lands on M2 as a one-mismatch hit in a shifted frame. That
/// hit used to be taken (wrong 1, with or without indels); the float now
/// abstains when such a hit has a competitor (here: V4 one insertion away),
/// so wrong is 0 and 2 of the 1,176 i5 variants abstain (one of them was
/// recovered before). The indel pass adds no wrong call.
#[test]
fn real_sheet_is_safe_for_indel_matching_on_both_indexes() {
    let sh = sheet();
    let mut d = Demuxer::new(&sh, cfg()).unwrap();
    let (s5, s7) = d.enable_indels(3);
    assert_eq!((s5.min_edit, s7.min_edit), (3, 3));
    assert_eq!((s5.wrong, s5.wrong_mismatch_only), (0, 0));
    assert_eq!((s7.wrong, s7.wrong_mismatch_only), (0, 0));
    assert_eq!(s5.variants, 14 * (30 + 10 + 44));
    assert_eq!((s5.recovered, s7.recovered), (1170, 1175));
    assert!(s5.safe() && s7.safe());
    // With a rigid frame nothing collides and every variant is recovered.
    let mut rigid = Demuxer::new(&sh, cfg()).unwrap();
    let (r5, r7) = rigid.enable_indels(0);
    assert_eq!((r5.wrong, r7.wrong), (0, 0));
    assert_eq!((r5.recovered, r7.recovered), (r5.variants, r7.variants));
}

#[test]
fn the_indel_pass_never_makes_a_wrong_call_of_its_own() {
    let sh = sheet();
    let mut on = Demuxer::new(&sh, cfg()).unwrap();
    on.enable_indels(3);
    let off = Demuxer::new(&sh, cfg()).unwrap();
    for (n, s) in sh.samples.iter().enumerate() {
        for v in indel_variants(&revcomp(&s.i5)) {
            let (seq, start) = i5_context(&v);
            let (o, ev) = on.assign_windowed_detail(&seq, Some(start), None, 3);
            assert!(ev.i5_seen);
            if o.sample.is_some_and(|x| x != n) {
                // Only the pre-existing shifted-frame mismatch hit, unchanged.
                assert!(
                    !ev.i5.via_indel,
                    "{} i5 {:?}",
                    s.id,
                    String::from_utf8_lossy(&v)
                );
                assert_eq!(
                    off.assign_windowed(&seq, Some(start), None, 3).sample,
                    o.sample
                );
            }
        }
        for v in indel_variants(&s.i7) {
            let (seq, start) = i7_context(&v);
            let (o, _) = on.assign_windowed_detail(&seq, None, Some(start), 3);
            assert!(o.sample.is_none() || o.sample == Some(n), "{} i7", s.id);
        }
    }
}

#[test]
fn indel_in_i5_is_unassigned_by_default_and_assigned_when_enabled() {
    let sh = sheet();
    let m2 = idx(&sh, "M2_GEX");
    let off = Demuxer::new(&sh, cfg()).unwrap();
    let mut on = Demuxer::new(&sh, cfg()).unwrap();
    on.enable_indels(3);
    let mut rescued = 0;
    for v in indel_variants(&revcomp(b"TACGATCAAG")) {
        let (seq, start) = i5_context(&v);
        let before = off.assign_windowed(&seq, Some(start), None, 3);
        let (after, ev) = on.assign_windowed_detail(&seq, Some(start), None, 3);
        // Every M2 indel variant is assigned to M2 with indels on ...
        assert_eq!(after.sample, Some(m2), "{:?}", String::from_utf8_lossy(&v));
        assert_eq!(after.via, DemuxVia::I5Only);
        // ... the ones mismatch-only already caught are unchanged ...
        if before.sample.is_some() {
            assert_eq!(before.sample, Some(m2));
            assert!(!ev.i5.via_indel);
        } else {
            // ... and the rest come from the indel pass.
            assert!(ev.i5.via_indel);
            rescued += 1;
        }
    }
    assert!(rescued >= 20, "rescued {rescued} of 54");
}

/// First indel variant of `key` that only the indel pass assigns.
fn indel_only_variant(
    on: &Demuxer,
    off: &Demuxer,
    key: &[u8],
    context: fn(&[u8]) -> (Vec<u8>, usize),
    i5: bool,
) -> Vec<u8> {
    indel_variants(key)
        .into_iter()
        .find(|v| {
            let (seq, st) = context(v);
            let (a, b) = if i5 {
                (Some(st), None)
            } else {
                (None, Some(st))
            };
            let (o, ev) = on.assign_windowed_detail(&seq, a, b, 3);
            o.sample.is_some()
                && (ev.i5.via_indel || ev.i7.via_indel)
                && off.assign_windowed(&seq, a, b, 3).sample.is_none()
        })
        .expect("an indel-only variant")
}

#[test]
fn indel_in_i7_upgrades_a_single_index_read_to_dual() {
    let sh = sheet();
    let m2 = idx(&sh, "M2_GEX");
    let mut d = Demuxer::new(&sh, cfg()).unwrap();
    d.enable_indels(3);
    let off = Demuxer::new(&sh, cfg()).unwrap();
    let i7v = indel_only_variant(&d, &off, b"TGTAGTCATT", i7_context, false);
    // One canonical stretch holding both regions: i5 context, then i7 context.
    let (mut seq, i5s) = i5_context(&revcomp(b"TACGATCAAG"));
    let (tail, i7_off) = i7_context(&i7v);
    let i7s = seq.len() + i7_off;
    seq.extend_from_slice(&tail);
    assert_eq!(
        off.assign_windowed(&seq, Some(i5s), Some(i7s), 3).via,
        DemuxVia::I5Only
    );
    let (o, ev) = d.assign_windowed_detail(&seq, Some(i5s), Some(i7s), 3);
    assert_eq!(o.sample, Some(m2));
    assert_eq!(o.via, DemuxVia::Dual);
    assert!(!ev.i5.via_indel && ev.i7.via_indel);
}

#[test]
fn a_mismatch_hit_is_never_overridden_by_the_indel_pass() {
    let sh = sheet();
    let m2 = idx(&sh, "M2_GEX");
    let mut d = Demuxer::new(&sh, cfg()).unwrap();
    d.enable_indels(3);
    let (seq, start) = i5_context(&revcomp(b"TACGATCAAG"));
    let (o, ev) = d.assign_windowed_detail(&seq, Some(start), None, 3);
    assert_eq!(o.sample, Some(m2));
    assert!(!ev.i5.via_indel);
}

#[test]
fn indel_hits_on_different_samples_stay_a_conflict() {
    let sh = sheet();
    let mut d = Demuxer::new(&sh, cfg()).unwrap();
    d.enable_indels(3);
    let off = Demuxer::new(&sh, cfg()).unwrap();
    // i5 of M2_GEX and i7 of M7_GEX, each with an indel only the indel pass reads.
    let i5v = indel_only_variant(&d, &off, &revcomp(b"TACGATCAAG"), i5_context, true);
    let i7v = indel_only_variant(&d, &off, b"TAGCATAGTG", i7_context, false);
    let (mut seq, i5s) = i5_context(&i5v);
    let (tail, i7_off) = i7_context(&i7v);
    let i7s = seq.len() + i7_off;
    seq.extend_from_slice(&tail);
    let o = d.assign_windowed(&seq, Some(i5s), Some(i7s), 3);
    assert_eq!(o.sample, None);
    assert_eq!(o.via, DemuxVia::Conflict);
}

#[test]
fn indexes_one_indel_apart_keep_that_index_mismatch_only() {
    // i7s: B is A with one base deleted and a base appended (edit distance 2);
    // i5s are far apart. i7 must be refused, i5 allowed.
    let sh = flat_sheet(&[
        ("A", "ACGTTGCATG", "AAAACCCCGG"),
        ("B", "ACGTGCATGA", "TTTTGGGGAA"),
        ("C", "GGCCAATTCC", "CCAAGGTTCA"),
    ]);
    let mut d = Demuxer::new(
        &sh,
        DemuxCfg {
            index_max_mm: 1,
            mode: Some(DemuxMode::Both),
        },
    )
    .unwrap();
    let (s5, s7) = d.enable_indels(3);
    assert!(s5.safe());
    assert_eq!(s7.min_edit, 2);
    assert!(!s7.safe());
    // No i7 call ever comes from the (refused) indel table.
    for v in indel_variants(b"ACGTTGCATG") {
        let (seq, start) = i7_context(&v);
        let (_, ev) = d.assign_windowed_detail(&seq, None, Some(start), 3);
        assert!(!ev.i7.via_indel);
    }
}

#[test]
fn evidence_reports_regions_outside_the_read_as_not_seen() {
    let sh = sheet();
    let d = Demuxer::new(&sh, cfg()).unwrap();
    let (seq, start) = i5_context(&revcomp(b"TACGATCAAG"));
    // i7 start past the end of the read: not seen, not identified.
    let (o, ev) = d.assign_windowed_detail(&seq, Some(start), Some(seq.len() - 4), 3);
    assert!(ev.i5_seen && !ev.i7_seen);
    assert_eq!(ev.i7.sample, None);
    assert_eq!(o.via, DemuxVia::I5Only);
    let (_, ev) = d.assign_windowed_detail(&seq, None, None, 3);
    assert!(!ev.i5_seen && !ev.i7_seen);
}
