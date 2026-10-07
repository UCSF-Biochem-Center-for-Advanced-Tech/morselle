use super::pipeline::extract;
use super::tests::{ecfg, ocfg};
use super::types::{ExtractCfg, FLEX_R2_END_ANCHOR, FLEX_R2_END_ANCHOR_MAX_SEARCH};
use crate::anchors::{P5, R1_PRIMER};
use crate::orient::orient;

fn flex_ecfg() -> ExtractCfg {
    ExtractCfg {
        flex_r2_boundary: true,
        ..ecfg()
    }
}

fn synthetic_flex(cb: &[u8], umi: &[u8], pre_anchor: &[u8], dt: usize) -> Vec<u8> {
    let mut v = Vec::new();
    v.extend_from_slice(b"AACAAT");
    v.extend_from_slice(P5);
    v.extend_from_slice(b"ACGTACGTAC");
    v.extend_from_slice(R1_PRIMER);
    v.extend_from_slice(cb);
    v.extend_from_slice(umi);
    v.extend(std::iter::repeat(b'T').take(dt));
    v.extend_from_slice(pre_anchor);
    v.extend_from_slice(FLEX_R2_END_ANCHOR);
    v.extend_from_slice(b"ACGTACGTACACGTACGTAC");
    v
}

#[test]
fn flex_r2_boundary_stops_exactly_at_the_anchor() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
    let seq = synthetic_flex(cb, umi, pre_anchor, 30);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);

    let (r2s, r2e) = e.r2.expect("flex R2 should be located");
    assert_eq!(&o.seq[r2s..r2e], pre_anchor.as_slice());
    assert!(e.adapter_trimmed);
    assert_eq!(e.adapter_pos, Some(r2e));
}

#[test]
fn flex_r2_boundary_tolerates_mismatches_within_budget() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
    let mut anchor = FLEX_R2_END_ANCHOR.to_vec();
    anchor[3] = if anchor[3] == b'A' { b'C' } else { b'A' };
    anchor[9] = if anchor[9] == b'A' { b'C' } else { b'A' };
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AACAAT");
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(pre_anchor);
    seq.extend_from_slice(&anchor);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
    let (r2s, r2e) =
        e.r2.expect("anchor within mismatch budget should still be found");
    assert_eq!(&o.seq[r2s..r2e], pre_anchor.as_slice());
}

#[test]
fn flex_r2_boundary_recovers_a_single_deletion_in_the_anchor() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
    let mut anchor = FLEX_R2_END_ANCHOR.to_vec();
    anchor.remove(7);
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AACAAT");
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(pre_anchor);
    seq.extend_from_slice(&anchor);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
    let (r2s, _r2e) =
        e.r2.expect("single-deletion-shifted anchor should still be found");
    assert_eq!(&o.seq[r2s..r2s + pre_anchor.len()], pre_anchor.as_slice());
}

#[test]
fn flex_r2_boundary_recovers_a_single_insertion_in_the_anchor() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
    let mut anchor = FLEX_R2_END_ANCHOR.to_vec();
    anchor.insert(7, b'G');
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AACAAT");
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(pre_anchor);
    seq.extend_from_slice(&anchor);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
    let (r2s, _r2e) =
        e.r2.expect("single-insertion-shifted anchor should still be found");
    assert_eq!(&o.seq[r2s..r2s + pre_anchor.len()], pre_anchor.as_slice());
}

#[test]
fn flex_r2_boundary_rejects_when_anchor_absent_and_payload_too_short_for_a_probe() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let payload = b"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AACAAT");
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(payload);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
    assert!(e.r2.is_none());
    assert!(!e.adapter_trimmed);
}

#[test]
fn flex_r2_boundary_falls_back_to_a_length_floor_when_anchor_absent_but_payload_long_enough() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let payload = b"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
    assert!(payload.len() >= 65);
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AACAAT");
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(payload);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
    let (r2s, r2e) =
        e.r2.expect("long-enough payload should fall back to a length floor");
    assert_eq!(&o.seq[r2s..r2e], payload.as_slice());
    assert!(!e.adapter_trimmed);
}

#[test]
fn flex_min_cdna_len_no_anchor_is_configurable() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let payload = b"GATTACAGATTACAGATTACA";
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AACAAT");
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(payload);
    let qual = vec![b'I'; seq.len()];
    let o = orient(&seq, &qual, &ocfg());

    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
    assert!(e.r2.is_none());

    let lenient = ExtractCfg {
        flex_min_cdna_len_no_anchor: payload.len(),
        ..flex_ecfg()
    };
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &lenient, false);
    assert!(e.r2.is_some());
}

#[test]
fn flex_r2_boundary_skips_search_for_known_non_gex_reads() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let payload = b"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AACAAT");
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(payload);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), true);
    let (r2s, r2e) = e.r2.expect("known-non-GEX reads must not be dropped");
    assert_eq!(r2e, o.seq.len());
    assert_eq!(&o.seq[r2s..r2e], payload.as_slice());
    assert!(!e.adapter_trimmed);
}

#[test]
fn flex_r2_boundary_search_is_bounded_past_the_real_structural_range() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let padding = vec![b'G'; FLEX_R2_END_ANCHOR_MAX_SEARCH + 20];
    let mut seq = Vec::new();
    seq.extend_from_slice(b"AACAAT");
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(&padding);
    seq.extend_from_slice(FLEX_R2_END_ANCHOR);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
    let (r2s, r2e) =
        e.r2.expect("length floor should still fire past the search bound");
    assert_eq!(
        r2e - r2s,
        FLEX_R2_END_ANCHOR_MAX_SEARCH,
        "bounded at the search limit, not the real anchor"
    );
    assert!(!e.adapter_trimmed, "not a verified landmark match");
}

#[test]
fn flex_r2_boundary_is_off_by_default_non_flex_path_unaffected() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
    let seq = synthetic_flex(cb, umi, pre_anchor, 30);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);
    let (r2s, r2e) =
        e.r2.expect("non-flex path should still emit R2 (long enough)");
    assert_eq!(r2e, o.seq.len());
    assert!(&o.seq[r2s..r2e] != pre_anchor.as_slice());
}
