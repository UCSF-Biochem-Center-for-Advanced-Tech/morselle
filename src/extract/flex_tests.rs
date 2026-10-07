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

/// Build a canonical Flex read: [leader][P5][i5][R1 primer][CB][UMI]
/// [dT run][junk before the anchor, standing in for BC10+constant29]
/// [FLEX_R2_END_ANCHOR][trailing junk, standing in for i7/P7]. Mirrors
/// `synthetic()`, but with the Flex-specific tail instead of the
/// standard-3'-capture `ILLUMINA_3P_ADAPTER` tail.
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
    v.extend_from_slice(b"ACGTACGTACACGTACGTAC");  // trailing i7/P7-ish junk
    v
}

#[test]
fn flex_r2_boundary_stops_exactly_at_the_anchor() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    // Stand-in for BC10(10)+constant29(29): any 39bp, content doesn't matter.
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
    // 2 mismatches inside the searched prefix (positions 3 and 9 of the
    // first FLEX_R2_END_ANCHOR_PREFIX_LEN bases) -- within max_mm budget.
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
    // Real SBX error mode (85.2% of recoverable Hamming-miss shapes):
    // a genuine anchor with one base deleted from inside the searched
    // prefix. A pure substitution-only scan can never find this at any
    // offset, since every base after the deletion point is shifted by 1.
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
    let mut anchor = FLEX_R2_END_ANCHOR.to_vec();
    anchor.remove(7);  // delete one base from inside the 15bp prefix
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
    // Rarer real shape (9.2% of recoverable Hamming-miss shapes) but
    // real: one extra base inserted inside the searched prefix.
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
    let mut anchor = FLEX_R2_END_ANCHOR.to_vec();
    anchor.insert(7, b'G');  // insert one extra base inside the 15bp prefix
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
    // likely_non_gex=false (the default -- caller has no reason to think
    // otherwise): no FLEX_R2_END_ANCHOR anywhere, AND the payload (63bp)
    // is below flex_min_cdna_len_no_anchor's real 65bp floor (50bp probe
    // + 15bp real search jitter) -- too short to structurally contain a
    // real probe match either way, so correctly rejected even under the
    // length-floor fallback.
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
    // Same shape as the too-short case, but 70bp of payload -- clears
    // flex_min_cdna_len_no_anchor's real 65bp floor even with no
    // FLEX_R2_END_ANCHOR anywhere. Real rationale: this read's
    // hit-rate-among-survivors on real CAT06 data is statistically the
    // same as the verified-anchor population (96.4-96.5% either way),
    // so it deserves a real shot at probeset-align/probe_barcode
    // matching instead of being dropped outright.
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
    // Not a verified landmark -- honest accounting, unlike the real
    // anchor-found path.
    assert!(!e.adapter_trimmed);
}

#[test]
fn flex_min_cdna_len_no_anchor_is_configurable() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let payload = b"GATTACAGATTACAGATTACA";  // 22bp -- below the real 65bp default
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

    // Default (65) rejects a 22bp payload.
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
    assert!(e.r2.is_none());

    // Lowering the floor to <= payload.len() accepts it.
    let lenient = ExtractCfg {
        flex_min_cdna_len_no_anchor: payload.len(),
        ..flex_ecfg()
    };
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &lenient, false);
    assert!(e.r2.is_some());
}

#[test]
fn flex_r2_boundary_skips_search_for_known_non_gex_reads() {
    // likely_non_gex=true: the caller's early i5-only sample-sheet
    // lookup (lib.rs::process_one) already identified this as e.g. a
    // CRISPR/Feature-Barcoding read sharing the same flex_r2_boundary=
    // true invocation as real Flex GEX reads. Must go straight to the
    // standard generic-adapter/read-end/min_cdna_len path -- same
    // outcome as non-Flex mode -- without ever attempting the (futile)
    // Flex anchor search.
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
    // The anchor sits real bp past cdna_start (FLEX_R2_END_ANCHOR_MAX_SEARCH,
    // ~150bp) -- a structurally-impossible position for the real
    // Ligated_Probe_insert/constant junction. Under the hard-gate default
    // (likely_non_gex=false) this must not be matched as the verified
    // landmark: confirms the search is actually bounded (an unbounded
    // scan would wrongly find it and mis-locate the boundary far too
    // late). The padding alone clears flex_min_cdna_len_no_anchor's real
    // 65bp floor though, so the length-floor fallback correctly still
    // emits an R2 -- bounded at the search limit, never at the real
    // (too-far-away) anchor position.
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
    // The same read that gets correctly bounded under flex_ecfg() should,
    // under the plain (non-Flex) ecfg(), fall through to the standard
    // generic-adapter/read-end behavior -- it has neither
    // ILLUMINA_3P_ADAPTER nor a length below min_cdna_len, so R2 just
    // runs to the read's own end, ignoring FLEX_R2_END_ANCHOR entirely.
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
