use super::engine::*;
use super::types::*;
use crate::anchors::revcomp;
use crate::sample_sheet::SampleSheet;
use std::path::PathBuf;

fn sheet() -> SampleSheet {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/CE15204_bclconvert_input.csv");
    SampleSheet::load(&p).unwrap()
}

/// Auto-resolves to `Dual` on the real fixture (both i5 and i7 are
/// standalone-safe there -- see `real_fixture_min_pairwise_hamming_matches_doc_comment`).
fn cfg() -> DemuxCfg {
    DemuxCfg {
        index_max_mm: 1,
        mode: None,
    }
}

#[test]
fn dual_index_exact() {
    let sh = sheet();
    let d = Demuxer::new(&sh, cfg()).unwrap();
    // M2_GEX: i7=TGTAGTCATT, i5=TACGATCAAG. Read carries rc(i5) and i7-fwd.
    let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
    let i5_read = revcomp(b"TACGATCAAG");
    let i7_read = b"TGTAGTCATT";
    let o = d.assign(Some(&i5_read), Some(i7_read));
    assert_eq!(o.sample, Some(m2));
    assert_eq!(o.via, DemuxVia::Dual);
    assert_eq!(d.sample_ids()[m2], "M2_GEX");
}

/// Pins the real fixture's own min pairwise Hamming distances, matching
/// this module's own doc comment ("min pairwise Hamming i7=3 / i5=4") --
/// if the fixture ever changes, this must be re-verified, not just
/// silently pass with different numbers.
#[test]
fn real_fixture_min_pairwise_hamming_matches_doc_comment() {
    let sh = sheet();
    assert_eq!(min_pairwise_hamming(&sh, IndexSource::I5), 4);
    assert_eq!(min_pairwise_hamming(&sh, IndexSource::I7), 3);
}

#[test]
fn single_index_i5_ignores_a_conflicting_i7() {
    let sh = sheet();
    let d = Demuxer::new(
        &sh,
        DemuxCfg {
            index_max_mm: 1,
            mode: Some(DemuxMode::I5),
        },
    )
    .unwrap();
    let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
    let i5_read = revcomp(b"TACGATCAAG"); // M2_GEX's real i5
    let wrong_i7 = b"AAAAAAAAAA"; // not any real sample's i7 -- must be ignored entirely
    let o = d.assign(Some(&i5_read), Some(wrong_i7));
    assert_eq!(o.sample, Some(m2));
    assert_eq!(o.via, DemuxVia::I5Only);
}

#[test]
fn single_index_refused_when_sheet_is_not_safe_at_this_mm() {
    let sh = sheet();
    // Real i5 min pairwise distance is 4 (see the pinned test above), so
    // unambiguous decoding needs index_max_mm <= 1 (> 2*mm required).
    // mm=2 needs min_dist > 4, which the real sheet does not have --
    // must be refused, not silently allowed.
    let result = Demuxer::new(
        &sh,
        DemuxCfg {
            index_max_mm: 2,
            mode: Some(DemuxMode::I5),
        },
    );
    assert!(result.is_err());
}

#[test]
fn one_mismatch_tolerated() {
    let sh = sheet();
    let d = Demuxer::new(&sh, cfg()).unwrap();
    let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
    let mut i5_read = revcomp(b"TACGATCAAG");
    i5_read[3] = if i5_read[3] == b'A' { b'C' } else { b'A' }; // 1 mm
    let o = d.assign(Some(&i5_read), Some(b"TGTAGTCATT"));
    assert_eq!(o.sample, Some(m2));
}

#[test]
fn single_index_when_i7_truncated() {
    let sh = sheet();
    let d = Demuxer::new(&sh, cfg()).unwrap();
    let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
    let i5_read = revcomp(b"TACGATCAAG");
    let o = d.assign(Some(&i5_read), None);
    assert_eq!(o.sample, Some(m2));
    assert_eq!(o.via, DemuxVia::I5Only);
}

#[test]
fn conflict_goes_undetermined() {
    let sh = sheet();
    let d = Demuxer::new(&sh, cfg()).unwrap();
    // i5 of M2_GEX but i7 of M7_GEX -> conflict
    let o = d.assign(Some(&revcomp(b"TACGATCAAG")), Some(b"TAGCATAGTG"));
    assert_eq!(o.sample, None);
    assert_eq!(o.via, DemuxVia::Conflict);
}

#[test]
fn unknown_index_undetermined() {
    let sh = sheet();
    let d = Demuxer::new(&sh, cfg()).unwrap();
    let o = d.assign(Some(b"GGGGGGGGGG"), Some(b"GGGGGGGGGG"));
    assert_eq!(o.sample, None);
    assert_eq!(o.via, DemuxVia::None);
}

#[test]
fn both_mode_requires_agreement_rejects_lone_i5() {
    let sh = sheet();
    let d = Demuxer::new(
        &sh,
        DemuxCfg {
            index_max_mm: 1,
            mode: Some(DemuxMode::Both),
        },
    )
    .unwrap();
    let o = d.assign(Some(&revcomp(b"TACGATCAAG")), None);
    assert_eq!(o.sample, None);
}

#[test]
fn both_mode_assigns_on_real_agreement() {
    let sh = sheet();
    let d = Demuxer::new(
        &sh,
        DemuxCfg {
            index_max_mm: 1,
            mode: Some(DemuxMode::Both),
        },
    )
    .unwrap();
    let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
    let o = d.assign(Some(&revcomp(b"TACGATCAAG")), Some(b"TGTAGTCATT"));
    assert_eq!(o.sample, Some(m2));
    assert_eq!(o.via, DemuxVia::Dual);
}

#[test]
fn auto_picks_dual_when_both_indices_are_standalone_safe() {
    // Real fixture: i5 min=4, i7 min=3. At index_max_mm=1 (required=3),
    // both are individually safe, so auto must pick Dual.
    let sh = sheet();
    let d = Demuxer::new(
        &sh,
        DemuxCfg {
            index_max_mm: 1,
            mode: None,
        },
    )
    .unwrap();
    assert_eq!(d.resolved_mode(), DemuxMode::Dual);
}

#[test]
fn auto_picks_both_when_neither_index_is_standalone_safe() {
    // Same real fixture, but index_max_mm=2 (required=5): neither i5
    // (4) nor i7 (3) clears that margin alone, so auto must fall back
    // to Both -- the only mode still safe by construction.
    let sh = sheet();
    let d = Demuxer::new(
        &sh,
        DemuxCfg {
            index_max_mm: 2,
            mode: None,
        },
    )
    .unwrap();
    assert_eq!(d.resolved_mode(), DemuxMode::Both);
}

#[test]
fn windowed_recovers_shifted_index() {
    let sh = sheet();
    let d = Demuxer::new(&sh, cfg()).unwrap();
    let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
    // M2 read-oriented i5 = revcomp(TACGATCAAG); place it at position 30.
    let i5_read = revcomp(b"TACGATCAAG");
    let mut seq = vec![b'G'; 80];
    seq[30..40].copy_from_slice(&i5_read);
    // Anchor is off by one (true index starts at 30, anchor says 29).
    let rigid = d.assign_windowed(&seq, Some(29), None, 0);
    assert_eq!(rigid.sample, None); // rigid frame misses it
    let windowed = d.assign_windowed(&seq, Some(29), None, 3);
    assert_eq!(windowed.sample, Some(m2)); // offset search recovers it
    assert_eq!(windowed.via, DemuxVia::I5Only);
}

#[test]
fn index_regions_slicing() {
    // R1_PRIMER (33) + i5 (10) before cb_start; rc(R2_PRIMER) (34) before i7.
    let cb_start = 43;
    let mut seq = vec![b'N'; 200];
    seq[0..10].copy_from_slice(b"AAAACCCCGG"); // i5 region [0,10)
    let adapter_pos = 100;
    seq[134..144].copy_from_slice(b"TTTTGGGGCC"); // i7 region [134,144)
    let (i5, i7) = index_regions(&seq, cb_start, Some(adapter_pos), 10);
    assert_eq!(i5, Some(&b"AAAACCCCGG"[..]));
    assert_eq!(i7, Some(&b"TTTTGGGGCC"[..]));
}
