//! R1-tail re-lock rescue (`--r1-tail-relock`, `XB:i:7`).
//!
//! Orientation locks the CB start with a mismatch-only match of the 10 bp
//! R1-primer tail (`CTTCCGATCT`, ≤ `--r1-tail-max-mm`). One insertion or
//! deletion inside that tail can still leave a frame within budget one or two
//! bp off, so the 16 bp CB grab is shifted and fails the whitelist. For such
//! CB-fails, re-find the last 15 bp of the R1 primer allowing up to two edits
//! (insertions and deletions included, bit-parallel Myers search) and read
//! the CB and UMI right after it. Accepted only when that CB is on the
//! whitelist or a unique Hamming-1 neighbour.
//!
//! Measured on CAT02 (private recovery-pools note, 2026-10-29): 3.2%
//! (barcode-first) / 0.9% (barcode-last) of CB-fails have a whitelist CB at
//! the edit-located tail end and not at today's lock, about +0.6 / +0.2% of
//! reads; those CBs land on cell-like barcodes at the exact-route rate.

use crate::anchors::R1_PRIMER;
use crate::barcode::Whitelist;
use crate::demux::{index_anchor_starts, Demuxer, Outcome};

/// The R1-primer tail searched with indels allowed (its last 15 bp).
const TAIL_LEN: usize = 15;
/// Edits allowed on the 15 bp tail.
const MAX_EDITS: u32 = 2;
/// Canonical bases searched from the read start (as orientation's prefix
/// search), extended to cover the current lock when it lies further in.
const SEARCH_WINDOW: usize = 120;

/// End positions (exclusive) and edit scores of the semi-global best matches
/// of `pat` (≤ 64 bp) in `text`, Myers' bit-parallel algorithm.
fn myers_scores(pat: &[u8], text: &[u8]) -> Vec<u32> {
    let m = pat.len();
    debug_assert!((1..=64).contains(&m));
    let mask = if m == 64 { u64::MAX } else { (1u64 << m) - 1 };
    let high = 1u64 << (m - 1);
    let mut peq = [0u64; 256];
    for (i, &c) in pat.iter().enumerate() {
        peq[c as usize] |= 1 << i;
    }
    let (mut pv, mut mv, mut score) = (mask, 0u64, m as u32);
    let mut out = Vec::with_capacity(text.len());
    for &c in text {
        let eq = peq[c as usize];
        let xv = eq | mv;
        let xh = ((((eq & pv).wrapping_add(pv)) & mask) ^ pv) | eq;
        let mut ph = mv | (!(xh | pv) & mask);
        let mut mh = pv & xh;
        if ph & high != 0 {
            score += 1;
        } else if mh & high != 0 {
            score -= 1;
        }
        ph = (ph << 1) & mask;
        mh = (mh << 1) & mask;
        pv = mh | (!(xv | ph) & mask);
        mv = ph & xv;
        out.push(score);
    }
    out
}

/// Candidate CB starts: every end of a best-scoring (≤2-edit) match of the
/// R1-primer tail in the read's first bases. An indel near the tail's 3' end
/// can tie several adjacent ends; the caller lets the whitelist pick.
pub(crate) fn relocked_cb_starts(seq: &[u8], cb_start: usize) -> Vec<usize> {
    let tail = &R1_PRIMER[R1_PRIMER.len() - TAIL_LEN..];
    let hay = &seq[..seq.len().min(SEARCH_WINDOW.max(cb_start + TAIL_LEN))];
    let scores = myers_scores(tail, hay);
    let Some(&best) = scores.iter().min() else {
        return Vec::new();
    };
    if best > MAX_EDITS {
        return Vec::new();
    }
    scores
        .iter()
        .enumerate()
        .filter(|&(_, &s)| s == best)
        .map(|(i, _)| i + 1)
        .collect()
}

/// The rescue: a corrected 28 bp R1 (CB + UMI), its qualities, and the
/// re-located CB start, read right after the re-located tail. Of the candidate starts other than `cb_start`,
/// exactly one may give a CB on (or uniquely within 1 of) the whitelist;
/// none or several → no rescue.
pub(crate) fn relock_r1(
    w: &Whitelist,
    seq: &[u8],
    qual: &[u8],
    cb_start: usize,
    cb_len: usize,
    umi_len: usize,
) -> Option<(Vec<u8>, Vec<u8>, usize)> {
    let mut hit: Option<(usize, [u8; 16])> = None;
    for s in relocked_cb_starts(seq, cb_start) {
        let end = s + cb_len + umi_len;
        if s == cb_start || end > seq.len() || end > qual.len() {
            continue;
        }
        if let Some(cb) = w.correct_within1(&seq[s..s + cb_len]) {
            if hit.is_some() {
                return None;
            }
            hit = Some((s, cb));
        }
    }
    let (s, cb) = hit?;
    let end = s + cb_len + umi_len;
    let mut r1 = cb.to_vec();
    r1.extend_from_slice(&seq[s + cb_len..end]);
    Some((r1, qual[s..end].to_vec(), s))
}

/// Demux again for a re-locked read that the early demux left unassigned:
/// that ran on the original lock, so its i5 anchor (counted back from the CB
/// start) is off by the same base or two as the CB was. CAT02 barcode-first
/// slice before this: 15.8% of re-locked reads undetermined vs 3.4% of exact
/// ones. A read the early demux already assigned keeps its sample.
pub(crate) fn reanchored_outcome(
    d: &Demuxer,
    seq: &[u8],
    cb_start: usize,
    adapter_pos: Option<usize>,
    index_len: usize,
    window: isize,
) -> Option<Outcome> {
    let (i5s, i7s) = index_anchor_starts(cb_start, adapter_pos, index_len);
    let o = d.assign_windowed(seq, i5s, i7s, window);
    o.sample.is_some().then_some(o)
}

#[cfg(test)]
mod tests {
    use super::*;

    const CB: &[u8] = b"AAACCCAAGAAACCCA";
    const UMI: &[u8] = b"ACGTACGTACGT";

    /// Leader + P5-ish spacer + `tail` (in place of the R1 primer) + CB + UMI + dT.
    fn read_with_primer(primer: &[u8]) -> (Vec<u8>, usize) {
        let mut s = b"AACAATAATGATACGGCGACCACCGAGATCTACACCTTGATCGTA".to_vec();
        s.extend_from_slice(primer);
        let cb_start = s.len();
        s.extend_from_slice(CB);
        s.extend_from_slice(UMI);
        s.extend(std::iter::repeat_n(b'T', 30));
        s.extend_from_slice(b"GACTGACTGACTGACTGACTGACTGACTGACT");
        (s, cb_start)
    }

    fn edit(primer: &[u8], f: impl Fn(&mut Vec<u8>)) -> Vec<u8> {
        let mut p = primer.to_vec();
        f(&mut p);
        p
    }

    #[test]
    fn myers_matches_plain_edit_distance_on_small_cases() {
        // pattern "ACGT" in text: best scores at each end.
        let sc = myers_scores(b"ACGT", b"TTACGTTT");
        assert_eq!(sc[5], 0); // ACGT ends at index 5 (exclusive end 6)
        assert_eq!(*sc.iter().min().unwrap(), 0);
        let sc = myers_scores(b"ACGT", b"TTACTTT"); // one deletion of G
        assert_eq!(*sc.iter().min().unwrap(), 1);
    }

    #[test]
    fn exact_tail_relocks_to_the_true_cb_start() {
        let (seq, cb) = read_with_primer(R1_PRIMER);
        assert_eq!(relocked_cb_starts(&seq, cb), vec![cb]);
        assert_eq!(relocked_cb_starts(&seq, cb + 2), vec![cb]);
    }

    #[test]
    fn single_and_double_indels_in_the_tail_still_lock_the_cb_start() {
        let n = R1_PRIMER.len();
        for i in (n - TAIL_LEN)..n {
            let del = edit(R1_PRIMER, |p| {
                p.remove(i);
            });
            let (seq, cb) = read_with_primer(&del);
            assert!(
                relocked_cb_starts(&seq, cb).contains(&cb),
                "deletion at {i}"
            );
            for &b in b"ACGT" {
                let ins = edit(R1_PRIMER, |p| p.insert(i, b));
                let (seq, cb) = read_with_primer(&ins);
                // The true start is a candidate (ties possible near the
                // tail's 3' end; the whitelist then decides) -- except an
                // extra base equal to the tail's last base right at its end,
                // which no aligner can tell from a CB starting with it: then
                // the candidate is one base early and the whitelist rejects it.
                let got = relocked_cb_starts(&seq, cb);
                let boundary_dup = i >= n - 1 && b == R1_PRIMER[n - 1];
                if boundary_dup {
                    assert_eq!(got, vec![cb - 1], "insertion {} at {i}", b as char);
                } else {
                    assert!(got.contains(&cb), "insertion {} at {i}: {got:?}", b as char);
                }
                assert!(got.len() <= 3, "insertion {} at {i}: {got:?}", b as char);
            }
        }
        let two = edit(R1_PRIMER, |p| {
            p.remove(n - 3);
            p.remove(n - 9);
        });
        let (seq, cb) = read_with_primer(&two);
        assert!(relocked_cb_starts(&seq, cb).contains(&cb));
    }

    #[test]
    fn a_read_without_the_tail_never_relocks() {
        let mut s = 0x2545_f491_u32;
        let mut hits = 0;
        for _ in 0..2000 {
            let seq: Vec<u8> = (0..400)
                .map(|_| {
                    s = s.wrapping_mul(1_103_515_245).wrapping_add(12345);
                    b"ACGT"[((s >> 16) & 3) as usize]
                })
                .collect();
            hits += usize::from(!relocked_cb_starts(&seq, 60).is_empty());
        }
        // Random sequence: a ≤2-edit 15-mer in 120 bp is rare (< 1%).
        assert!(hits < 20, "{hits} false locks in 2000 random reads");
    }
}

#[cfg(test)]
mod rescue_tests {
    use super::*;

    fn whitelist(cbs: &[&[u8]]) -> Whitelist {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("wl.txt");
        let text: String = cbs
            .iter()
            .map(|c| format!("{}\n", String::from_utf8_lossy(c)))
            .collect();
        std::fs::write(&p, text).unwrap();
        Whitelist::load(&p).unwrap()
    }

    #[test]
    fn rescue_reads_cb_and_umi_after_an_indel_shifted_tail() {
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let mut primer = R1_PRIMER.to_vec();
        primer.remove(R1_PRIMER.len() - 4); // deletion inside CTTCCGATCT
        let mut seq = b"AACAATAATGATACGGCGACCACCGAGATCTACACCTTGATCGTA".to_vec();
        seq.extend_from_slice(&primer);
        let true_start = seq.len();
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat_n(b'T', 30));
        let qual = vec![b'I'; seq.len()];
        let w = whitelist(&[cb]);
        // Today's lock one base late: the grab starts inside the CB.
        let (r1, q, s) = relock_r1(&w, &seq, &qual, true_start + 1, 16, 12).unwrap();
        assert_eq!(s, true_start);
        assert_eq!(&r1[..16], cb);
        assert_eq!(&r1[16..], umi);
        assert_eq!(q.len(), 28);
        // Already at the true start: nothing to re-lock.
        assert!(relock_r1(&w, &seq, &qual, true_start, 16, 12).is_none());
        // CB not on the whitelist: no rescue.
        let other = whitelist(&[b"TTTTGGGGCCCCAAAA"]);
        assert!(relock_r1(&other, &seq, &qual, true_start + 1, 16, 12).is_none());
    }
}
