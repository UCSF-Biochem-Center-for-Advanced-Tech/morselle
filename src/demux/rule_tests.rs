//! The per-read rule: look up both indexes; assign on one index when it
//! names exactly one sample, otherwise on the single sample both agree on.
use super::engine::*;
use super::types::*;
use crate::anchors::{revcomp, P5, P7, R1_PRIMER, R2_PRIMER};
use crate::sample_sheet::SampleSheet;

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

fn demuxer(sh: &SampleSheet, mode: Option<DemuxMode>) -> Demuxer {
    Demuxer::new(
        sh,
        DemuxCfg {
            index_max_mm: 1,
            mode,
        },
    )
    .unwrap()
}

fn idx(sh: &SampleSheet, id: &str) -> usize {
    sh.samples.iter().position(|s| s.id == id).unwrap()
}

const X: &str = "TTATTCGAGG";
const Y: &str = "TGCGCGGTTT";
const P: &str = "AGCAGGACAG";
const Q: &str = "TTTATCCTTG";
const R: &str = "TTCGACAAGC";

/// A canonical read: P5 + CTACAC + rc(i5) + R1 primer + insert +
/// rc(R2 primer) + i7 + rc(P7), with the i5/i7 anchor starts.
fn read(i5_read: &[u8], i7_read: &[u8]) -> (Vec<u8>, usize, usize) {
    let mut seq = [P5, b"CTACAC"].concat();
    let i5 = seq.len();
    seq.extend_from_slice(i5_read);
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(&b"ACGTTGCA".repeat(8));
    seq.extend_from_slice(&revcomp(R2_PRIMER));
    let i7 = seq.len();
    seq.extend_from_slice(i7_read);
    seq.extend_from_slice(&revcomp(P7));
    (seq, i5, i7)
}

fn i5r(sheet_i5: &str) -> Vec<u8> {
    revcomp(sheet_i5.as_bytes())
}

/// Combinatorial sheet: one i7 shared by three samples, told apart by i5.
fn shared_i7_sheet() -> SampleSheet {
    flat_sheet(&[
        ("A", X, P),
        ("B", X, Q),
        ("C", X, R),
        ("D", Y, "CAGGGTTGGC"),
    ])
}

#[test]
fn a_shared_i7_alone_never_assigns_but_i5_picks_among_its_samples() {
    let sh = shared_i7_sheet();
    let d = demuxer(&sh, None);
    // Before the rule auto fell back to i5 alone here (i7 min Hamming 0).
    assert_eq!(d.resolved_mode(), DemuxMode::Dual);
    let b = idx(&sh, "B");
    let o = d.assign(Some(&i5r(Q)), Some(X.as_bytes()));
    assert_eq!((o.sample, o.via), (Some(b), DemuxVia::Dual));
    let o = d.assign(None, Some(X.as_bytes()));
    assert_eq!((o.sample, o.via), (None, DemuxVia::None));
    let o = d.assign(Some(&i5r(Q)), None);
    assert_eq!((o.sample, o.via), (Some(b), DemuxVia::I5Only));
    // D's i7 is its own: alone it decides.
    let o = d.assign(None, Some(Y.as_bytes()));
    assert_eq!((o.sample, o.via), (Some(idx(&sh, "D")), DemuxVia::I7Only));
    // i5 of A with D's i7: no sample has both.
    let o = d.assign(Some(&i5r(P)), Some(Y.as_bytes()));
    assert_eq!((o.sample, o.via), (None, DemuxVia::Conflict));
    // Unknown i5 with the shared i7: still undetermined.
    let o = d.assign(Some(b"GGGGGGGGGG"), Some(X.as_bytes()));
    assert_eq!(o.sample, None);
}

#[test]
fn a_full_grid_of_shared_indexes_needs_both_for_every_sample() {
    let sh = flat_sheet(&[("A", X, P), ("B", Y, P), ("C", X, Q), ("D", Y, Q)]);
    let d = demuxer(&sh, None);
    for (id, i7, i5) in [("A", X, P), ("B", Y, P), ("C", X, Q), ("D", Y, Q)] {
        let (seq, s5, s7) = read(&i5r(i5), i7.as_bytes());
        let (o, ev) = d.assign_windowed_detail(&seq, Some(s5), Some(s7), 3);
        assert_eq!(
            (o.sample, o.via),
            (Some(idx(&sh, id)), DemuxVia::Dual),
            "{id}"
        );
        assert!(ev.i5.shared.is_some() && ev.i7.shared.is_some());
        assert!(ev.i5.matched() && ev.i7.matched());
        let o = d.assign(Some(&i5r(i5)), None);
        assert_eq!(o.sample, None, "{id} on i5 alone");
        let o = d.assign(None, Some(i7.as_bytes()));
        assert_eq!(o.sample, None, "{id} on i7 alone");
    }
}

#[test]
fn overlapping_mismatch_neighbourhoods_are_resolved_by_the_other_index() {
    // i5s 2 mismatches apart: at --index-max-mm 1 the midpoint is one
    // mismatch from both. It used to match nothing; now it is {A, B}.
    let (pa, pb, mid) = ("AGCAGGACAG", "AGCAGGACTT", "AGCAGGACTG");
    let sh = flat_sheet(&[("A", X, pa), ("B", Y, pb)]);
    let d = demuxer(&sh, None);
    let o = d.assign(Some(&i5r(mid)), Some(Y.as_bytes()));
    assert_eq!((o.sample, o.via), (Some(idx(&sh, "B")), DemuxVia::Dual));
    let o = d.assign(Some(&i5r(mid)), None);
    assert_eq!(o.sample, None);
    // The exact indexes still decide alone.
    let o = d.assign(Some(&i5r(pa)), None);
    assert_eq!((o.sample, o.via), (Some(idx(&sh, "A")), DemuxVia::I5Only));
}

#[test]
fn both_mode_and_single_index_modes_keep_their_meaning() {
    let sh = shared_i7_sheet();
    let both = demuxer(&sh, Some(DemuxMode::Both));
    let o = both.assign(Some(&i5r(Q)), Some(X.as_bytes()));
    assert_eq!(o.sample, Some(idx(&sh, "B")));
    let o = both.assign(Some(&i5r(Q)), None);
    assert_eq!(o.sample, None);
    // i7 alone cannot separate A/B/C: refused, as before.
    assert!(Demuxer::new(
        &sh,
        DemuxCfg {
            index_max_mm: 1,
            mode: Some(DemuxMode::I7)
        }
    )
    .is_err());
}

#[test]
fn an_indel_in_a_shared_i7_still_gives_its_candidate_set() {
    let sh = shared_i7_sheet();
    let mut d = demuxer(&sh, None);
    let (s5, s7) = d.enable_indels(3);
    assert!(s5.safe() && s7.safe(), "{s5:?} {s7:?}");
    let mut i7 = X.as_bytes().to_vec();
    i7.remove(4);
    let (seq, a5, a7) = read(&i5r(R), &i7);
    let (o, ev) = d.assign_windowed_detail(&seq, Some(a5), Some(a7), 3);
    assert_eq!((o.sample, o.via), (Some(idx(&sh, "C")), DemuxVia::Dual));
    assert!(ev.i7.via_indel && ev.i7.shared.is_some());
}
