//! Per-read orientation → canonical P5→P7 (plan §4 step 2, M1).
//!
//! Homopolymer-first: the 3' GEX dT-primer run (poly-T in canonical, poly-A when
//! the read is P7-led/reversed) is found by a cheap O(n) scan and used to both
//! orient the read and bracket the CB+UMI block. A fixed-offset 10 bp Hamming
//! check of the R1-primer tail confirms the bracket for free. Only when no
//! qualifying homopolymer is found do we fall back to a `triple_accel` Levenshtein
//! search for the R1 primer (the expensive path, kept off the common case).
//!
//! Scope: single 3' assay. The 5' (TSO-spacer) path is stubbed — see
//! [`orient_5p_stub`].

use std::borrow::Cow;

use crate::anchors::{hamming, revcomp, R1_PRIMER, R1_PRIMER_TAIL, SBX_LEADER};
use crate::homopolymer::longest_run;
use crate::pipeline::Orientation;

/// How orientation was determined, for accounting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrientMethod {
    /// Oriented from the dT-primer homopolymer (common case).
    Homopolymer,
    /// Oriented from a fallback R1-primer Levenshtein search.
    R1Fallback,
    /// Could not be oriented.
    None,
}

/// Tuning for [`orient`].
#[derive(Debug, Clone, Copy)]
pub struct OrientCfg {
    pub cb_len: usize,
    pub umi_len: usize,
    /// Minimum homopolymer run length to accept the dT-primer anchor.
    pub min_homopolymer: usize,
    /// Max mismatches for the 10 bp R1-tail bracket confirm.
    pub r1_tail_max_mm: u32,
    /// Max edits for the fallback R1-primer search.
    pub anchor_max_edits: u32,
    /// Optional orientation hint (e.g. from the input dir), used ONLY to break
    /// genuine ties. `None` = anchors only.
    pub expected: Option<Orientation>,
}

impl OrientCfg {
    /// CB+UMI bracket width.
    fn bracket(&self) -> usize {
        self.cb_len + self.umi_len
    }
}

/// Result of orienting one read to canonical P5→P7.
#[derive(Debug)]
pub struct Oriented<'a> {
    /// Canonical sequence (borrowed when no flip, owned revcomp when flipped).
    pub seq: Cow<'a, [u8]>,
    /// Canonical quality, aligned to `seq`.
    pub qual: Cow<'a, [u8]>,
    pub orientation: Orientation,
    pub method: OrientMethod,
    /// SBX leader present at the 5' end of the *original* read.
    pub leader_present: bool,
    /// Start of the dT-primer run in canonical coords (= end of CB+UMI).
    pub polyt_start: Option<usize>,
    /// CB+UMI block `[start, end)` in canonical coords.
    pub cb_umi: Option<(usize, usize)>,
    /// The 10 bp R1-primer-tail bracket confirm passed.
    pub r1_tail_confirmed: bool,
}

/// Detect the SBX leader in the first ≤8 bp of the original read (≤1 mm).
fn detect_leader(seq: &[u8]) -> bool {
    let n = SBX_LEADER.len().min(seq.len());
    n == SBX_LEADER.len() && hamming(&seq[..n], SBX_LEADER) <= 1
}

/// Hamming distance of the R1-primer tail ending just before `cb_start`.
fn r1_tail_mm(canon: &[u8], cb_start: usize) -> Option<u32> {
    let tail = R1_PRIMER_TAIL.len();
    (cb_start >= tail).then(|| hamming(&canon[cb_start - tail..cb_start], R1_PRIMER_TAIL))
}

/// Window (± bases) searched around the homopolymer-implied bracket start when
/// locking it to the R1-primer tail — covers dT-boundary fuzz (e.g. a UMI ending
/// in T bleeding into the run) without a full alignment.
const BRACKET_LOCK_WINDOW: isize = 6;

/// Lock the CB+UMI bracket precisely: slide the 10 bp R1-primer tail over a small
/// window around `approx_cb_start` and return the offset-corrected `cb_start` with
/// the best (≤ `max_mm`) match. `None` if the tail can't be confirmed nearby.
fn lock_bracket(canon: &[u8], approx_cb_start: usize, max_mm: u32) -> Option<usize> {
    let mut best: Option<(u32, usize)> = None;
    for d in -BRACKET_LOCK_WINDOW..=BRACKET_LOCK_WINDOW {
        let cand = approx_cb_start as isize + d;
        if cand < 0 {
            continue;
        }
        let cand = cand as usize;
        if let Some(mm) = r1_tail_mm(canon, cand) {
            if mm <= max_mm && best.map_or(true, |(bm, _)| mm < bm) {
                best = Some((mm, cand));
            }
        }
    }
    best.map(|(_, s)| s)
}

/// Orient a single read to canonical P5→P7 (3' assay).
pub fn orient<'a>(seq: &'a [u8], qual: &'a [u8], cfg: &OrientCfg) -> Oriented<'a> {
    let leader_present = detect_leader(seq);
    let bracket = cfg.bracket();
    let n = seq.len();

    // --- Stage A: homopolymer scan (cheap) ---
    let t_run = longest_run(seq, b'T').filter(|r| r.len >= cfg.min_homopolymer);
    let a_run = longest_run(seq, b'A').filter(|r| r.len >= cfg.min_homopolymer);

    // A T-run is a canonical dT-primer candidate only if there's room for the
    // CB+UMI block before it. An A-run is the reversed dT-primer only if there's
    // room for the (revcomp) CB+UMI block after it.
    let t_ok = t_run.filter(|r| r.start >= bracket);
    let a_ok = a_run.filter(|r| n.saturating_sub(r.end()) >= bracket);

    let flip = match (t_ok, a_ok) {
        (Some(t), Some(a)) => {
            // Both plausible: longer run wins; tie → expected hint; else canonical.
            if t.len > a.len {
                false
            } else if a.len > t.len {
                true
            } else {
                matches!(cfg.expected, Some(Orientation::P7Led))
            }
        }
        (Some(_), None) => false,
        (None, Some(_)) => true,
        (None, None) => return fallback(seq, qual, cfg, leader_present),
    };

    // --- Build canonical seq/qual + dT-primer start in canonical coords ---
    let (canon_seq, canon_qual, polyt_start): (Cow<[u8]>, Cow<[u8]>, usize) = if flip {
        let a = a_ok.expect("a_ok set when flip");
        // A-run [start, end) maps to a T-run at [n-end, n-start) after revcomp.
        let polyt_start = n - a.end();
        let mut q = qual.to_vec();
        q.reverse();
        (Cow::Owned(revcomp(seq)), Cow::Owned(q), polyt_start)
    } else {
        let t = t_ok.expect("t_ok set when !flip");
        (Cow::Borrowed(seq), Cow::Borrowed(qual), t.start)
    };

    let orientation = if flip {
        Orientation::P7Led
    } else {
        Orientation::P5Led
    };

    // The homopolymer locates the bracket approximately; lock it precisely to the
    // R1-primer tail over a small window (handles dT-boundary fuzz). Fall back to
    // the raw homopolymer offset (unconfirmed) if the tail isn't found nearby.
    let approx_cb_start = polyt_start.checked_sub(bracket);
    let locked = approx_cb_start.and_then(|s| lock_bracket(&canon_seq, s, cfg.r1_tail_max_mm));
    let r1_tail_confirmed = locked.is_some();
    let cb_start = locked.or(approx_cb_start);
    let cb_umi = cb_start
        .filter(|&s| s + bracket <= canon_seq.len())
        .map(|s| (s, s + bracket));
    let polyt_refined = cb_umi.map(|(_, e)| e).unwrap_or(polyt_start);

    Oriented {
        seq: canon_seq,
        qual: canon_qual,
        orientation,
        method: OrientMethod::Homopolymer,
        leader_present,
        polyt_start: Some(polyt_refined),
        cb_umi,
        r1_tail_confirmed,
    }
}

/// Fallback: locate the R1 primer on both strands with a Levenshtein search and
/// orient from the better hit. Only reached when no qualifying homopolymer exists.
fn fallback<'a>(
    seq: &'a [u8],
    qual: &'a [u8],
    cfg: &OrientCfg,
    leader_present: bool,
) -> Oriented<'a> {
    let rc_primer = revcomp(R1_PRIMER);
    let fwd = best_search(R1_PRIMER, seq);
    let rev = best_search(&rc_primer, seq);

    let thr = cfg.anchor_max_edits;
    let bracket = cfg.bracket();

    let (orientation, flip, r1_end) = match (fwd, rev) {
        (Some((fe, fd)), Some((re, rd))) if fd.min(rd) <= thr => {
            if fd <= rd {
                (Orientation::P5Led, false, fe)
            } else {
                (Orientation::P7Led, true, re)
            }
        }
        (Some((fe, fd)), _) if fd <= thr => (Orientation::P5Led, false, fe),
        (_, Some((re, rd))) if rd <= thr => (Orientation::P7Led, true, re),
        _ => {
            return Oriented {
                seq: Cow::Borrowed(seq),
                qual: Cow::Borrowed(qual),
                orientation: Orientation::Unknown,
                method: OrientMethod::None,
                leader_present,
                polyt_start: None,
                cb_umi: None,
                r1_tail_confirmed: false,
            }
        }
    };

    let n = seq.len();
    let (canon_seq, canon_qual): (Cow<[u8]>, Cow<[u8]>) = if flip {
        let mut q = qual.to_vec();
        q.reverse();
        (Cow::Owned(revcomp(seq)), Cow::Owned(q))
    } else {
        (Cow::Borrowed(seq), Cow::Borrowed(qual))
    };

    // r1_end is on the as-given strand; map to canonical coords if we flipped.
    let cb_start = if flip { n - r1_end } else { r1_end };
    let cb_umi = (cb_start + bracket <= canon_seq.len()).then(|| (cb_start, cb_start + bracket));
    let polyt_start = cb_umi.map(|(_, e)| e);

    Oriented {
        seq: canon_seq,
        qual: canon_qual,
        orientation,
        method: OrientMethod::R1Fallback,
        leader_present,
        polyt_start,
        cb_umi,
        r1_tail_confirmed: false,
    }
}

/// Best Levenshtein match of `needle` in `haystack`: `(end_index, distance)`.
fn best_search(needle: &[u8], haystack: &[u8]) -> Option<(usize, u32)> {
    triple_accel::levenshtein_search(needle, haystack)
        .min_by_key(|m| m.k)
        .map(|m| (m.end, m.k))
}

/// 5' (TSO-spacer) assay orientation — not yet implemented.
///
/// The current data is all 3' GEX (poly-dT spacer). When 5' libraries are added,
/// orientation there keys off the TSO and the distal poly-A tail instead of the
/// dT primer. Tracked as a future extension.
pub fn orient_5p_stub<'a>(_seq: &'a [u8], _qual: &'a [u8], _cfg: &OrientCfg) -> Oriented<'a> {
    unimplemented!("5' (TSO) assay orientation not yet implemented")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchors::{P5, P7, R2_PRIMER};

    fn cfg() -> OrientCfg {
        OrientCfg {
            cb_len: 16,
            umi_len: 12,
            min_homopolymer: 12,
            r1_tail_max_mm: 2,
            anchor_max_edits: 6,
            expected: None,
        }
    }

    /// Build a synthetic canonical 3' read from the anchor consts.
    fn synthetic_canonical(cb: &[u8], umi: &[u8]) -> Vec<u8> {
        let mut v = Vec::new();
        v.extend_from_slice(SBX_LEADER);
        v.extend_from_slice(P5);
        v.extend_from_slice(b"ACGTACGTAC"); // i5 (10)
        v.extend_from_slice(R1_PRIMER);
        v.extend_from_slice(cb); // 16
        v.extend_from_slice(umi); // 12
        v.extend_from_slice(&[b'T'; 30]); // poly(dT)
        v.extend_from_slice(b"GACTGACTGACTGACTGACTGACT"); // cDNA
        v.extend_from_slice(R2_PRIMER);
        v.extend_from_slice(b"TGCATGCATG"); // i7 (10)
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
        assert_eq!(o.seq.as_ref(), seq.as_slice()); // borrowed, unchanged
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
        // Flipped back to the exact canonical sequence.
        assert_eq!(o.seq.as_ref(), canonical.as_slice());
        let (s, e) = o.cb_umi.unwrap();
        assert_eq!(&o.seq[s..s + 16], cb);
        assert_eq!(&o.seq[s + 16..e], umi);
        // Qual was reversed to track the flip.
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

    #[test]
    fn fallback_orients_without_homopolymer() {
        // Canonical structure but with the dT run too short for the cheap path,
        // forcing the R1-primer Levenshtein fallback.
        let mut seq = Vec::new();
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(b"AAACCCAAGAAACCCA"); // CB
        seq.extend_from_slice(b"ACGTACGTACGT"); // UMI
        seq.extend_from_slice(b"TTTTT"); // short dT (< min_homopolymer)
        seq.extend_from_slice(b"GACTGACTGACTGACTGACT");
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &cfg());
        assert_eq!(o.method, OrientMethod::R1Fallback);
        assert_eq!(o.orientation, Orientation::P5Led);
        let (s, _e) = o.cb_umi.unwrap();
        assert_eq!(&o.seq[s..s + 16], b"AAACCCAAGAAACCCA");
    }
}
