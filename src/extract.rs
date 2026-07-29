//! R1/R2 bracket extraction and adapter/homopolymer trimming (plan §4 step 4, M2).
//!
//! Operates on the canonical (P5→P7) read produced by M1, reusing the coordinates
//! it already found:
//! - **R1** = the 28 bp CB+UMI block located by [`crate::orient`].
//! - **R2** = cDNA between the dT-primer run end and the 3' read-through adapter
//!   (`AGATCGGAAGAGC`), found by a cheap Hamming seed; if the read ends mid-cDNA
//!   (no adapter), R2 runs to the read end.
//!
//! All positions are returned as ranges into the canonical sequence — zero-copy;
//! the caller slices `seq`/`qual`.

use crate::anchors::{hamming_leq, ILLUMINA_3P_ADAPTER};

/// Tuning for [`extract`].
#[derive(Debug, Clone, Copy)]
pub struct ExtractCfg {
    /// Minimum cDNA length to emit an R2.
    pub min_cdna_len: usize,
    /// Max mismatches for the 3' adapter trim seed.
    pub adapter_max_mm: u32,
    /// Cap emitted R2 length; 0 = uncapped. SBX cDNA can be >1 kb, but 10x reference
    /// hard-caps reads at 649 bp. NOT applied here — `extract` retains the whole
    /// fragment; the cap is applied by the emitter after reverse-complementing into
    /// R2 orientation, so it keeps the bases the R2 primer would read.
    pub max_r2_len: usize,
}

/// Extraction result: ranges `[start, end)` into the canonical sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extracted {
    /// CB+UMI range (R1), 28 bp, or `None` if not located.
    pub r1: Option<(usize, usize)>,
    /// cDNA range (R2) in the canonical (polyT) frame, or `None` if absent / below
    /// `min_cdna_len`. NOTE: this slice is on the same strand as R1; the emitter
    /// reverse-complements it (and reverses its quality) so R2 follows the Illumina
    /// paired-end convention and maps sense to the gene.
    pub r2: Option<(usize, usize)>,
    /// End of the dT-primer run = cDNA start.
    pub polyt_end: Option<usize>,
    /// Length of the dT-primer run actually observed (QC).
    pub polyt_run_len: u16,
    /// Position of the 3' read-through adapter (= start of `rc(R2 primer)`), when
    /// found. Marks the cDNA end and lets demux locate the downstream i7.
    pub adapter_pos: Option<usize>,
    /// The 3' read-through adapter was found and trimmed.
    pub adapter_trimmed: bool,
}

/// Extend the dT-primer run rightward from `approx_start` (tolerating isolated
/// mismatches) and return the cDNA start = one past the run's last T. If no T is
/// near `approx_start`, returns `approx_start` (run length 0 — polyT absent).
fn dt_run_end(seq: &[u8], approx_start: usize) -> usize {
    let n = seq.len();
    // Probe for the run start at approx_start, then ±1 (boundary fuzz).
    let start = if approx_start < n && seq[approx_start] == b'T' {
        approx_start
    } else if approx_start > 0 && approx_start - 1 < n && seq[approx_start - 1] == b'T' {
        approx_start - 1
    } else if approx_start + 1 < n && seq[approx_start + 1] == b'T' {
        approx_start + 1
    } else {
        return approx_start.min(n);
    };
    let mut last = start;
    let mut i = start + 1;
    while i < n {
        if seq[i] == b'T' {
            last = i;
        } else if seq[i - 1] != b'T' {
            break; // two consecutive non-T -> run ends
        }
        i += 1;
    }
    last + 1
}

/// Minimum 3' read-through adapter overlap to trim a partial adapter at the read
/// end. Short partials must match near-exactly (mismatch budget scales with
/// length, below), so 7 bp is safe — chance of a 7 bp exact adapter prefix at a
/// read's 3' end is ~6e-5.
const MIN_PARTIAL_ADAPTER: usize = 7;

/// Leftmost position ≥ `from` where the 3' read-through adapter matches. Catches:
///   * full-length match within `max_mm` substitutions (the common case);
///   * full-length match with a single indel + ≤`max_mm` substitutions — SBX makes
///     indel errors in homopolymers (e.g. AAGAGC→AAAGAGC), which a pure-Hamming
///     seed misses because the insertion frame-shifts it to many mismatches;
///   * a partial adapter (≥ `MIN_PARTIAL_ADAPTER` bp) running off the 3' end.
/// Trimming the adapter ourselves keeps cDNA out of downstream tools' adapter
/// trimmers (10x reference mis-trims long reads that still carry the adapter).
/// Allocation-free. The indel check runs only at a cheap 3 bp prefix anchor.
fn find_adapter(seq: &[u8], from: usize, max_mm: u32) -> Option<usize> {
    let a = ILLUMINA_3P_ADAPTER;
    let m = a.len();
    let end = seq.len();
    let mut i = from;
    while i < end {
        let remain = end - i;
        if remain >= m {
            if hamming_leq(&seq[i..i + m], a, max_mm) {
                return Some(i);
            }
            // Indel path: anchor on the indel-safe 3 bp prefix (no homopolymer),
            // then verify a single-indel alignment.
            if seq[i] == a[0]
                && seq[i + 1] == a[1]
                && seq[i + 2] == a[2]
                && adapter_single_indel(&seq[i..], a, max_mm)
            {
                return Some(i);
            }
        } else if remain >= MIN_PARTIAL_ADAPTER {
            // Partial adapter running off the 3' end: match the adapter prefix,
            // mismatch budget scaling with overlap (7 bp exact; 8+ bp allow 1) —
            // false trims stay ~0.04% and only cost a few 3' bp, while leaving an
            // adapter risks crashing downstream trimmers.
            let allowed = (remain / 8) as u32;
            if hamming_leq(&seq[i..end], &a[..remain], allowed) {
                return Some(i);
            }
        } else {
            break;
        }
        i += 1;
    }
    None
}

/// True if adapter `a` matches `s` (which starts at the candidate position) with
/// exactly one indel and ≤`max_mm` substitutions. Tries an insertion in the read
/// (extra base, e.g. homopolymer expansion) and a deletion in the read.
fn adapter_single_indel(s: &[u8], a: &[u8], max_mm: u32) -> bool {
    let m = a.len();
    // Insertion in read: skip one read base. s[0..k]≈a[0..k], s[k+1..m+1]≈a[k..m].
    if s.len() >= m + 1 {
        for k in 0..m {
            let mut mm = 0u32;
            for j in 0..k {
                if s[j] != a[j] {
                    mm += 1;
                }
            }
            for j in k..m {
                if s[j + 1] != a[j] {
                    mm += 1;
                }
            }
            if mm <= max_mm {
                return true;
            }
        }
    }
    // Deletion in read: skip one adapter base. s[0..k]≈a[0..k], s[k..m-1]≈a[k+1..m].
    if s.len() >= m - 1 {
        for k in 0..m {
            let mut mm = 0u32;
            for j in 0..k {
                if s[j] != a[j] {
                    mm += 1;
                }
            }
            for j in (k + 1)..m {
                if s[j - 1] != a[j] {
                    mm += 1;
                }
            }
            if mm <= max_mm {
                return true;
            }
        }
    }
    false
}

/// Bracket R1 (CB+UMI) and R2 (cDNA) on a canonical read.
///
/// `cb_umi` and `polyt_start` come from [`crate::orient::Oriented`].
pub fn extract(
    seq: &[u8],
    cb_umi: Option<(usize, usize)>,
    polyt_start: Option<usize>,
    cfg: &ExtractCfg,
) -> Extracted {
    let mut polyt_end = None;
    let mut polyt_run_len = 0u16;
    let mut r2 = None;
    let mut adapter_pos = None;
    let mut adapter_trimmed = false;

    if let Some(pts) = polyt_start {
        let cdna_start = dt_run_end(seq, pts);
        polyt_run_len = cdna_start.saturating_sub(pts).min(u16::MAX as usize) as u16;
        polyt_end = Some(cdna_start);

        adapter_pos = find_adapter(seq, cdna_start, cfg.adapter_max_mm);
        let cdna_end = adapter_pos.unwrap_or(seq.len());
        adapter_trimmed = adapter_pos.is_some();
        // Retain the WHOLE adapter-trimmed cDNA fragment here; any `--max-r2-len`
        // cap is applied at write time, after the fragment is reverse-complemented
        // into R2 orientation, so the cap keeps the bases the R2 primer would
        // actually read (see the emitter in lib.rs). See `ExtractCfg::max_r2_len`.
        if cdna_end > cdna_start && cdna_end - cdna_start >= cfg.min_cdna_len {
            r2 = Some((cdna_start, cdna_end));
        }
    }

    Extracted {
        r1: cb_umi,
        r2,
        polyt_end,
        polyt_run_len,
        adapter_pos,
        adapter_trimmed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchors::{revcomp, P5, R1_PRIMER};
    use crate::orient::{orient, FallbackMode, OrientCfg};
    use crate::pipeline::Orientation;

    #[test]
    fn adapter_exact_and_substitution() {
        let cdna = b"ACGTACGTACGTACGTACGT"; // 20bp, no adapter
        let mut s = cdna.to_vec();
        let astart = s.len();
        s.extend_from_slice(ILLUMINA_3P_ADAPTER); // exact adapter
        s.extend_from_slice(b"CACGTCTGAACTCC");
        assert_eq!(find_adapter(&s, 0, 2), Some(astart));
        // one substitution in the adapter still found
        let mut s2 = cdna.to_vec();
        s2.extend_from_slice(b"AGATCGGAAGTGC"); // A->T at pos 10
        assert_eq!(find_adapter(&s2, 0, 2), Some(cdna.len()));
    }

    #[test]
    fn adapter_homopolymer_insertion() {
        // The real compatibility-breaking case: AAGAGC -> AAAGAGC (extra A).
        let cdna = b"ACGTACGTACGTACGTACGT";
        let mut s = cdna.to_vec();
        let astart = s.len();
        s.extend_from_slice(b"AGATCGGAAAGAGCACACGTCTGAACTCC"); // indel adapter
                                                               // pure-Hamming (no indel handling) would miss this; find_adapter catches it
        assert_eq!(find_adapter(&s, 0, 2), Some(astart));
    }

    #[test]
    fn adapter_partial_at_end() {
        let cdna = b"ACGTACGTACGTACGTACGT";
        let mut s = cdna.to_vec();
        let astart = s.len();
        s.extend_from_slice(&ILLUMINA_3P_ADAPTER[..11]); // 11bp partial (>= MIN_PARTIAL)
        assert_eq!(find_adapter(&s, 0, 2), Some(astart));
        // a 7bp exact partial (the dominant real residual) IS called
        let mut s7 = cdna.to_vec();
        s7.extend_from_slice(&ILLUMINA_3P_ADAPTER[..7]); // "AGATCGG"
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

    fn ecfg() -> ExtractCfg {
        ExtractCfg {
            min_cdna_len: 20,
            adapter_max_mm: 2,
            max_r2_len: 0,
        }
    }

    fn ocfg() -> OrientCfg {
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
    fn synthetic(cb: &[u8], umi: &[u8], cdna: &[u8], dt: usize) -> Vec<u8> {
        let mut v = Vec::new();
        v.extend_from_slice(b"AACAAT"); // leader
        v.extend_from_slice(P5);
        v.extend_from_slice(b"ACGTACGTAC"); // i5
        v.extend_from_slice(R1_PRIMER);
        v.extend_from_slice(cb);
        v.extend_from_slice(umi);
        v.extend(std::iter::repeat(b'T').take(dt)); // dT primer
        v.extend_from_slice(cdna);
        v.extend_from_slice(ILLUMINA_3P_ADAPTER);
        v.extend_from_slice(b"CACGTCTGAACTCCAGTCACTGCATGCATG"); // rest of rc(R2)+i7
        v
    }

    #[test]
    fn extracts_r1_and_trimmed_r2() {
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let cdna = b"GGCATTACTGGCATTACTGGCATTACTGGCATTACTG"; // 37 bp, no long T-run
        let seq = synthetic(cb, umi, cdna, 30);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        assert_eq!(o.orientation, Orientation::P5Led);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg());

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
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg());
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
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg());
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
        let cdna = b"GGCATTAC"; // 8 bp < min_cdna_len
        let seq = synthetic(cb, umi, cdna, 30);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg());
        assert!(e.r2.is_none());
        assert!(e.r1.is_some());
    }
}
