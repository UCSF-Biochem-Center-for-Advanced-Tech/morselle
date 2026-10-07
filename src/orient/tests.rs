use crate::anchors::{revcomp, P5, P7, R1_PRIMER, R2_PRIMER, SBX_LEADER};
use crate::pipeline::Orientation;
use super::matcher::orient;
use super::types::{FallbackMode, OrientCfg, OrientMethod};

fn cfg() -> OrientCfg {
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

fn synthetic_canonical(cb: &[u8], umi: &[u8]) -> Vec<u8> {
    let mut v = Vec::new();
    v.extend_from_slice(SBX_LEADER);
    v.extend_from_slice(P5);
    v.extend_from_slice(b"ACGTACGTAC");
    v.extend_from_slice(R1_PRIMER);
    v.extend_from_slice(cb);
    v.extend_from_slice(umi);
    v.extend_from_slice(&[b'T'; 30]);
    v.extend_from_slice(b"GACTGACTGACTGACTGACTGACT");
    v.extend_from_slice(R2_PRIMER);
    v.extend_from_slice(b"TGCATGCATG");
    v.extend_from_slice(P7);
    v
}

#[test]
fn orients_canonical_no_flip() {
    let cb = b"AAACCCAAGAAACCCA";
    let umi = b"ACGTACGTACGT";
    let seq = synthetic_canonical(cb, umi);
    let qual = vec![b'I'; seq.len()];

    let o = orient(&seq, &qual, &cfg());
    assert_eq!(o.orientation, Orientation::P5Led);
    assert_eq!(o.method, OrientMethod::Homopolymer);
    assert!(o.leader_present);
    assert!(o.r1_tail_confirmed);
    let (s, e) = o.cb_umi.unwrap();
    assert_eq!(&o.seq[s..s + 16], cb);
    assert_eq!(&o.seq[s + 16..e], umi);
    assert_eq!(o.seq.as_ref(), seq.as_slice());
}

#[test]
fn orients_reversed_and_flips_back() {
    let cb = b"AAACCCATCAAACCTG";
    let umi = b"TTTTGGGGCCCC";
    let canonical = synthetic_canonical(cb, umi);
    let reversed = revcomp(&canonical);
    let qual: Vec<u8> = (0..reversed.len() as u8).map(|i| b'!' + (i % 40)).collect();

    let o = orient(&reversed, &qual, &cfg());
    assert_eq!(o.orientation, Orientation::P7Led);
    assert_eq!(o.method, OrientMethod::Homopolymer);
    assert!(o.r1_tail_confirmed);
    assert_eq!(o.seq.as_ref(), canonical.as_slice());
    let (s, e) = o.cb_umi.unwrap();
    assert_eq!(&o.seq[s..s + 16], cb);
    assert_eq!(&o.seq[s + 16..e], umi);
    let mut q_rev = qual.clone();
    q_rev.reverse();
    assert_eq!(o.qual.as_ref(), q_rev.as_slice());
}

#[test]
fn unknown_when_no_anchor() {
    let seq = b"GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC".to_vec();
    let qual = vec![b'I'; seq.len()];
    let o = orient(&seq, &qual, &cfg());
    assert_eq!(o.orientation, Orientation::Unknown);
    assert_eq!(o.method, OrientMethod::None);
    assert!(o.cb_umi.is_none());
}

fn no_homopolymer_read() -> Vec<u8> {
    let mut seq = Vec::new();
    seq.extend_from_slice(P5);
    seq.extend_from_slice(b"ACGTACGTAC");
    seq.extend_from_slice(R1_PRIMER);
    seq.extend_from_slice(b"AAACCCAAGAAACCCA");
    seq.extend_from_slice(b"ACGTACGTACGT");
    seq.extend_from_slice(b"GGGTT");
    seq.extend_from_slice(b"GACGACGACGACGACGACGA");
    seq
}

#[test]
fn hamming_tail_fallback_orients() {
    let seq = no_homopolymer_read();
    let qual = vec![b'I'; seq.len()];
    let o = orient(&seq, &qual, &cfg());
    assert_eq!(o.method, OrientMethod::TailFallback);
    assert_eq!(o.orientation, Orientation::P5Led);
    assert!(o.r1_tail_confirmed);
    let (s, _e) = o.cb_umi.unwrap();
    assert_eq!(&o.seq[s..s + 16], b"AAACCCAAGAAACCCA");
}

#[test]
fn hamming_tail_fallback_flips_reversed() {
    let seq = revcomp(&no_homopolymer_read());
    let qual = vec![b'I'; seq.len()];
    let o = orient(&seq, &qual, &cfg());
    assert_eq!(o.method, OrientMethod::TailFallback);
    assert_eq!(o.orientation, Orientation::P7Led);
    let (s, _e) = o.cb_umi.unwrap();
    assert_eq!(&o.seq[s..s + 16], b"AAACCCAAGAAACCCA");
}

#[test]
fn levenshtein_fallback_orients() {
    let seq = no_homopolymer_read();
    let qual = vec![b'I'; seq.len()];
    let mut c = cfg();
    c.fallback = FallbackMode::Levenshtein;
    let o = orient(&seq, &qual, &c);
    assert_eq!(o.method, OrientMethod::LevenshteinFallback);
    assert_eq!(o.orientation, Orientation::P5Led);
    let (s, _e) = o.cb_umi.unwrap();
    assert_eq!(&o.seq[s..s + 16], b"AAACCCAAGAAACCCA");
}

#[test]
fn fixed_orientation_reanchors_past_false_homopolymer() {
    let cb = b"AAACCCAAGAAACACT";
    let umi = b"ACGTACGTACGT";
    let mut s = Vec::new();
    s.extend_from_slice(P5);
    s.extend_from_slice(b"ACGTACGTAC");
    s.extend_from_slice(R1_PRIMER);
    s.extend_from_slice(cb);
    s.extend_from_slice(umi);
    s.extend(std::iter::repeat(b'T').take(20));
    s.extend_from_slice(b"GACGACGACGACGACGACGACGACGACGAC");
    s.extend(std::iter::repeat(b'T').take(40));
    s.extend_from_slice(b"GACGACGACGAC");
    let qual = vec![b'I'; s.len()];

    let mut c = cfg();
    c.expected = Some(Orientation::P5Led);
    let o = orient(&s, &qual, &c);
    assert_eq!(o.orientation, Orientation::P5Led);
    let (cs, _) = o.cb_umi.unwrap();
    assert_eq!(&o.seq[cs..cs + 16], cb);
}

#[test]
fn fixed_orientation_p7led_flips_and_locates() {
    let cb = b"AAACCCATCAAACCTG";
    let umi = b"TTTTGGGGCCCC";
    let canonical = synthetic_canonical(cb, umi);
    let reversed = revcomp(&canonical);
    let qual = vec![b'I'; reversed.len()];

    let mut c = cfg();
    c.expected = Some(Orientation::P7Led);
    let o = orient(&reversed, &qual, &c);
    assert_eq!(o.orientation, Orientation::P7Led);
    let (cs, _) = o.cb_umi.unwrap();
    assert_eq!(&o.seq[cs..cs + 16], cb);
}

#[test]
fn fallback_off_marks_unknown() {
    let seq = no_homopolymer_read();
    let qual = vec![b'I'; seq.len()];
    let mut c = cfg();
    c.fallback = FallbackMode::Off;
    let o = orient(&seq, &qual, &c);
    assert_eq!(o.orientation, Orientation::Unknown);
}
