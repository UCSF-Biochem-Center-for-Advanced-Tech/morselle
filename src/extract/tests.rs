use super::pipeline::{extract, find_adapter};
use super::types::{ExtractCfg, TSO_RC};
use crate::anchors::{revcomp, ILLUMINA_3P_ADAPTER, P5, R1_PRIMER};
use crate::orient::{orient, FallbackMode, OrientCfg};
use crate::pipeline::Orientation;

#[test]
fn adapter_exact_and_substitution() {
    let cdna = b"ACGTACGTACGTACGTACGT";  // 20bp, no adapter
    let mut s = cdna.to_vec();
    let astart = s.len();
    s.extend_from_slice(ILLUMINA_3P_ADAPTER);  // exact adapter
    s.extend_from_slice(b"CACGTCTGAACTCC");
    assert_eq!(find_adapter(&s, 0, 2), Some(astart));

    // one substitution in the adapter still found
    let mut s2 = cdna.to_vec();
    s2.extend_from_slice(b"AGATCGGAAGTGC");  // A->T at pos 10
    assert_eq!(find_adapter(&s2, 0, 2), Some(cdna.len()));
}

#[test]
fn adapter_homopolymer_insertion() {
    // The real compatibility-breaking case: AAGAGC -> AAAGAGC (extra A).
    let cdna = b"ACGTACGTACGTACGTACGT";
    let mut s = cdna.to_vec();
    let astart = s.len();
    s.extend_from_slice(b"AGATCGGAAAGAGCACACGTCTGAACTCC");  // indel adapter
                                                           // pure-Hamming (no indel handling) would miss this; find_adapter catches it
    assert_eq!(find_adapter(&s, 0, 2), Some(astart));
}

#[test]
fn adapter_partial_at_end() {
    let cdna = b"ACGTACGTACGTACGTACGT";
    let mut s = cdna.to_vec();
    let astart = s.len();
    s.extend_from_slice(&ILLUMINA_3P_ADAPTER[..11]);  // 11bp partial (>= MIN_PARTIAL)
    assert_eq!(find_adapter(&s, 0, 2), Some(astart));

    // a 7bp exact partial (the dominant real residual) IS called
    let mut s7 = cdna.to_vec();
    s7.extend_from_slice(&ILLUMINA_3P_ADAPTER[..7]);  // "AGATCGG"
    assert_eq!(find_adapter(&s7, 0, 2), Some(cdna.len()));

    // a too-short tail (< MIN_PARTIAL) is not called
    let mut s2 = cdna.to_vec();
    s2.extend_from_slice(&ILLUMINA_3P_ADAPTER[..5]);
    assert_eq!(find_adapter(&s2, 0, 2), None);
}

#[test]
fn adapter_absent_not_falsely_trimmed() {
    let s = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    assert_eq!(find_adapter(s, 0, 2), None);
}

/// Build a canonical read: CB+UMI | dT | rc(cDNA) | [rc(TSO)] | [adapter].
fn canonical(cdna: &[u8], with_tso: bool, with_adapter: bool) -> Vec<u8> {
    // A cDNA ending in A makes revcomp start with T, which the dT-run walker
    // legitimately absorbs — a fixture trap, not a bug. Refuse to build one.
    assert_ne!(*cdna.last().unwrap(), b'A', "cDNA must not end in A");
    let mut v = vec![b'A'; 28];
    v.extend(std::iter::repeat_n(b'T', 30));
    v.extend_from_slice(&revcomp(cdna));
    if with_tso {
        v.extend_from_slice(TSO_RC);
    }
    if with_adapter {
        v.extend_from_slice(ILLUMINA_3P_ADAPTER);
    }
    v
}

pub(crate) fn ecfg() -> ExtractCfg {
    ExtractCfg {
        min_cdna_len: 20,
        adapter_max_mm: 2,
        trim_tso: false,
        max_r2_len: 0,
        flex_r2_boundary: false,
        flex_min_cdna_len_no_anchor: 65,
    }
}

fn tso_cfg() -> ExtractCfg {
    ExtractCfg {
        trim_tso: true,
        ..ecfg()
    }
}

pub(crate) fn ocfg() -> OrientCfg {
    OrientCfg {
        cb_len: 16,
        umi_len: 12,
        min_homopolymer: 12,
        r1_tail_max_mm: 2,
        anchor_max_edits: 6,
        fallback: FallbackMode::HammingTail,
        expected: None,
    }
}

/// Canonical 3' read with a real 3' read-through adapter after the cDNA.
pub(crate) fn synthetic(cb: &[u8], umi: &[u8], cdna: &[u8], dt: usize) -> Vec<u8> {
    let mut v = Vec::new();
    v.extend_from_slice(b"AACAAT");  // leader
    v.extend_from_slice(P5);
    v.extend_from_slice(b"ACGTACGTAC");  // i5
    v.extend_from_slice(R1_PRIMER);
    v.extend_from_slice(cb);
    v.extend_from_slice(umi);
    v.extend(std::iter::repeat(b'T').take(dt));  // dT primer
    v.extend_from_slice(cdna);
    v.extend_from_slice(ILLUMINA_3P_ADAPTER);
    v.extend_from_slice(b"CACGTCTGAACTCCAGTCACTGCATGCATG");  // rest of rc(R2)+i7
    v
}

#[test]
fn tso_rc_matches_the_declared_oligo() {
    assert_eq!(TSO_RC, revcomp(crate::anchors::TSO).as_slice());
}

#[test]
fn trims_a_full_tso_between_cdna_and_adapter() {
    let cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCG";
    let seq = canonical(cdna, true, true);
    let e = extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false);
    assert!(e.tso_trimmed);
    let (s, t) = e.r2.unwrap();
    assert_eq!(&seq[s..t], revcomp(cdna).as_slice());
    assert_eq!(e.tso_pos, Some(t));
}

#[test]
fn leaves_the_tso_alone_by_default() {
    let cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCG";
    let seq = canonical(cdna, true, true);
    let e = extract(&seq, Some((0, 28)), Some(28), &ecfg(), false);
    assert!(!e.tso_trimmed);
    let (s, t) = e.r2.unwrap();
    // The TSO is still sitting on the end of the "cDNA".
    assert_eq!(t - s, cdna.len() + TSO_RC.len());
}

#[test]
fn trims_a_partial_tso_running_off_the_read_end() {
    let cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCG";
    let mut seq = canonical(cdna, false, false);
    seq.extend_from_slice(&TSO_RC[..18]);  // read stopped mid-TSO
    let e = extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false);
    assert!(e.tso_trimmed);
    let (s, t) = e.r2.unwrap();
    assert_eq!(&seq[s..t], revcomp(cdna).as_slice());
}

#[test]
fn tolerates_mismatches_in_the_tso() {
    let cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCG";
    let mut seq = canonical(cdna, true, true);
    let tso_at = seq.len() - ILLUMINA_3P_ADAPTER.len() - TSO_RC.len();
    // 3 mismatches in 30 bp is inside the 10% budget.
    for offset in [4usize, 11, 19] {
        seq[tso_at + offset] = if seq[tso_at + offset] == b'A' {
            b'C'
        } else {
            b'A'
        };
    }
    assert!(extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false).tso_trimmed);

    // 6 is not.
    let mut seq = canonical(cdna, true, true);
    for offset in [2usize, 6, 10, 14, 18, 22] {
        seq[tso_at + offset] = if seq[tso_at + offset] == b'A' {
            b'C'
        } else {
            b'A'
        };
    }
    assert!(!extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false).tso_trimmed);
}

#[test]
fn does_not_trim_a_tso_buried_inside_the_cdna() {
    // The oligo appears mid-cDNA, not abutting the adapter. Anchoring is what
    // stops a chance match from truncating real sequence.
    let mut cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCG".to_vec();
    cdna.extend_from_slice(&revcomp(TSO_RC));
    cdna.extend_from_slice(b"ATCGATCGATGCTAGCTAGCATCGATCGTC");
    let seq = canonical(&cdna, false, true);
    let e = extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false);
    assert!(!e.tso_trimmed);
    let (s, t) = e.r2.unwrap();
    assert_eq!(&seq[s..t], revcomp(&cdna).as_slice());
}

#[test]
fn a_fragment_that_is_only_tso_yields_no_cdna() {
    // Trimming happens before the length gate, so an adapter-dimer-like
    // fragment is rejected rather than emitted as 30 bp of "cDNA".
    let mut seq = vec![b'A'; 28];
    seq.extend(std::iter::repeat_n(b'T', 30));
    seq.extend_from_slice(TSO_RC);
    seq.extend_from_slice(ILLUMINA_3P_ADAPTER);
    let e = extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false);
    assert!(e.tso_trimmed);
    assert!(e.r2.is_none());
}

#[test]
fn extracts_r1_and_trimmed_r2() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let cdna = b"GGCATTACTGGCATTACTGGCATTACTGGCATTACTG";  // 37 bp, no long T-run
    let seq = synthetic(cb, umi, cdna, 30);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    assert_eq!(o.orientation, Orientation::P5Led);
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);

    // R1 == 28 bp and equals CB+UMI.
    let (r1s, r1e) = e.r1.unwrap();
    assert_eq!(r1e - r1s, 28);
    assert_eq!(&o.seq[r1s..r1s + 16], cb);
    assert_eq!(&o.seq[r1s + 16..r1e], umi);

    // R2 == cDNA, adapter trimmed.
    assert!(e.adapter_trimmed);
    let (r2s, r2e) = e.r2.unwrap();
    assert_eq!(&o.seq[r2s..r2e], cdna);
}

#[test]
fn r2_runs_to_end_without_adapter() {
    let cb = b"AAACCCATCAAACCTG";
    let umi = b"TTTTGGGGCCCC";
    let cdna = b"GGCATTACTGGCATTACTGGCATTACTGGCATTACTG";
    // build a read that simply ends in cDNA (truncated, no adapter)
    let mut seq = Vec::new();
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(cb);
    seq.extend_from_slice(umi);
    seq.extend(std::iter::repeat(b'T').take(30));
    seq.extend_from_slice(cdna);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);
    assert!(!e.adapter_trimmed);
    let (r2s, r2e) = e.r2.unwrap();
    assert_eq!(r2e, o.seq.len());
    assert_eq!(&o.seq[r2s..r2e], cdna);
}

#[test]
fn reversed_read_extracts_same() {
    let cb = b"AAACCCAAGAAACCTG";
    let umi = b"ACGTACGTACGT";
    let cdna = b"GGCATTACTGGCATTACTGGCATTACTGGCATTACTG";
    let canonical = synthetic(cb, umi, cdna, 30);
    let reversed = revcomp(&canonical);
    let qual = vec![b'I'; reversed.len()];

    let o = orient(&reversed, &qual, &ocfg());
    assert_eq!(o.orientation, Orientation::P7Led);
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);
    let (r1s, r1e) = e.r1.unwrap();
    assert_eq!(r1e - r1s, 28);
    assert_eq!(&o.seq[r1s..r1s + 16], cb);
    let (r2s, r2e) = e.r2.unwrap();
    assert_eq!(&o.seq[r2s..r2e], cdna);
}

#[test]
fn short_cdna_dropped() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let cdna = b"GGCATTAC";  // 8 bp < min_cdna_len
    let seq = synthetic(cb, umi, cdna, 30);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &ocfg());
    let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);
    assert!(e.r2.is_none());
    assert!(e.r1.is_some());
}
