use std::path::PathBuf;
use crate::anchors::revcomp;
use crate::sample_sheet::SampleSheet;
use super::engine::*;
use super::types::*;

fn sheet() -> SampleSheet {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/CE15204_bclconvert_input.csv");
    SampleSheet::load(&p).unwrap()
}

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
    let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
    let i5_read = revcomp(b"TACGATCAAG");
    let i7_read = b"TGTAGTCATT";
    let o = d.assign(Some(&i5_read), Some(i7_read));
    assert_eq!(o.sample, Some(m2));
    assert_eq!(o.via, DemuxVia::Dual);
    assert_eq!(d.sample_ids()[m2], "M2_GEX");
}

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
    let i5_read = revcomp(b"TACGATCAAG");
    let wrong_i7 = b"AAAAAAAAAA";
    let o = d.assign(Some(&i5_read), Some(wrong_i7));
    assert_eq!(o.sample, Some(m2));
    assert_eq!(o.via, DemuxVia::I5Only);
}

#[test]
fn single_index_refused_when_sheet_is_not_safe_at_this_mm() {
    let sh = sheet();
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
    i5_read[3] = if i5_read[3] == b'A' { b'C' } else { b'A' };
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
    let i5_read = revcomp(b"TACGATCAAG");
    let mut seq = vec![b'G'; 80];
    seq[30..40].copy_from_slice(&i5_read);
    let rigid = d.assign_windowed(&seq, Some(29), None, 0);
    assert_eq!(rigid.sample, None);
    let windowed = d.assign_windowed(&seq, Some(29), None, 3);
    assert_eq!(windowed.sample, Some(m2));
    assert_eq!(windowed.via, DemuxVia::I5Only);
}

#[test]
fn index_regions_slicing() {
    let cb_start = 43;
    let mut seq = vec![b'N'; 200];
    seq[0..10].copy_from_slice(b"AAAACCCCGG");
    let adapter_pos = 100;
    seq[134..144].copy_from_slice(b"TTTTGGGGCC");
    let (i5, i7) = index_regions(&seq, cb_start, Some(adapter_pos), 10);
    assert_eq!(i5, Some(&b"AAAACCCCGG"[..]));
    assert_eq!(i7, Some(&b"TTTTGGGGCC"[..]));
}
