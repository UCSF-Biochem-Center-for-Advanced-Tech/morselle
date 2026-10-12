//! poly(dT)/poly(A) run detection (plan §4, M1/M2).
//!
//! The 3' GEX dT-primer (~30 nt) sits immediately 3' of the CB+UMI block, so its
//! homopolymer run is the cheap O(n) structural anchor used to orient reads and
//! bracket the barcode — no adapter alignment on the hot path.

/// A detected homopolymer run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Run {
    /// Start index of the run (inclusive).
    pub start: usize,
    /// Span length: `last_matching_base - start + 1`.
    pub len: usize,
}

impl Run {
    /// One past the last base of the run.
    pub fn end(&self) -> usize {
        self.start + self.len
    }
}

/// Longest near-homopolymer run of `base` in `seq`, tolerating *isolated* single
/// mismatches (the run breaks only on two consecutive non-`base` bases). Returns
/// `None` if no base matches. O(n), single pass.
pub fn longest_run(seq: &[u8], base: u8) -> Option<Run> {
    let n = seq.len();
    let mut best: Option<Run> = None;
    let mut i = 0;
    while i < n {
        if seq[i] != base {
            i += 1;
            continue;
        }
        // Extend from i, allowing isolated mismatches.
        let start = i;
        let mut last_match = i;
        let mut k = i + 1;
        while k < n {
            if seq[k] == base {
                last_match = k;
            } else if seq[k - 1] != base {
                // two consecutive non-base bases -> end the run
                break;
            }
            k += 1;
        }
        let len = last_match - start + 1;
        if best.map_or(true, |b| len > b.len) {
            best = Some(Run { start, len });
        }
        i = last_match + 1;
    }
    best
}

/// Longest A-run and longest T-run in one convenience call.
pub fn poly_a_t(seq: &[u8]) -> (Option<Run>, Option<Run>) {
    (longest_run(seq, b'A'), longest_run(seq, b'T'))
}

/// Leftmost run of `base` at least `min_len` long (same isolated-mismatch
/// tolerance as [`longest_run`]), or `None` if no run reaches the threshold.
///
/// Unlike `longest_run`, this stops at the first qualifying run instead of
/// scanning for the overall longest. For locating where real sequence ends
/// and a homopolymer tail (e.g. read-through into a poly-A tail) begins, the
/// leftmost qualifying run is the one that marks the boundary — a longer run
/// further along could be coincidental downstream sequence.
pub fn first_run_at_least(seq: &[u8], base: u8, min_len: usize) -> Option<Run> {
    let n = seq.len();
    let mut i = 0;
    while i < n {
        if seq[i] != base {
            i += 1;
            continue;
        }
        let start = i;
        let mut last_match = i;
        let mut k = i + 1;
        while k < n {
            if seq[k] == base {
                last_match = k;
            } else if seq[k - 1] != base {
                break;
            }
            k += 1;
        }
        let len = last_match - start + 1;
        if len >= min_len {
            return Some(Run { start, len });
        }
        i = last_match + 1;
    }
    None
}

/// Leftmost poly-A tail (read-through) of `base`: at least `min_len` long,
/// at least 90% `base` (one mismatch per 10 bases), mismatches isolated.
/// Starts and ends on `base`; extends as far as it stays at 90%.
///
/// Stricter than [`first_run_at_least`], whose every-other-base tolerance
/// fires on A-rich transcript sequence. On DZEC04, `--trim-r2-polya` with
/// that rule discarded 124 K counted reads (5.8%) that aligned across all
/// 150 bases at MAPQ 60, cut at genomic runs such as ALB's
/// `AGAATAAGAGAAAGAAAA` (67% A; 37% of ALB's reads). A real read-through
/// tail is the ~30-base dT primer and passes easily.
pub fn first_tail_run(seq: &[u8], base: u8, min_len: usize) -> Option<Run> {
    let n = seq.len();
    for start in 0..n {
        if seq[start] != base {
            continue;
        }
        let (mut mismatches, mut best_end) = (0usize, None);
        for k in start..n {
            if seq[k] != base {
                if seq[k - 1] != base {
                    break;
                }
                mismatches += 1;
                continue;
            }
            let len = k - start + 1;
            if len >= min_len && mismatches * 10 <= len {
                best_end = Some(k);
            }
        }
        if let Some(end) = best_end {
            return Some(Run {
                start,
                len: end - start + 1,
            });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_clean_run() {
        let s = b"ACGTTTTTTTTTTGCA";
        let r = longest_run(s, b'T').unwrap();
        assert_eq!(r.start, 3);
        assert_eq!(r.len, 10);
        assert_eq!(r.end(), 13);
    }

    #[test]
    fn tolerates_isolated_mismatch() {
        // single C inside the T-run should not break it
        let s = b"GGTTTTTCTTTTTGG";
        let r = longest_run(s, b'T').unwrap();
        assert_eq!(r.start, 2);
        assert_eq!(r.len, 11); // spans across the single C
    }

    #[test]
    fn breaks_on_two_consecutive() {
        let s = b"TTTTTGGTTT";
        let r = longest_run(s, b'T').unwrap();
        assert_eq!(r.start, 0);
        assert_eq!(r.len, 5);
    }

    #[test]
    fn none_when_absent() {
        assert_eq!(longest_run(b"GGGCCCGGG", b'T'), None);
    }

    #[test]
    fn picks_a_and_t() {
        let s = b"AAAAAAGCTTTTTTTTT";
        let (a, t) = poly_a_t(s);
        assert_eq!(a.unwrap().len, 6);
        assert_eq!(t.unwrap().len, 9);
    }

    #[test]
    fn first_run_at_least_finds_leftmost_qualifying_run_not_the_longest() {
        // Two qualifying A-runs (len 8 and len 16); must return the first.
        let s = b"ACGTAAAAAAAAGGGGAAAAAAAAAAAAAAAA";
        let r = first_run_at_least(s, b'A', 8).unwrap();
        assert_eq!(r.start, 4);
        assert_eq!(r.len, 8);
    }

    #[test]
    fn first_run_at_least_returns_none_below_threshold() {
        let s = b"ACGTAAAGGGG";
        assert_eq!(first_run_at_least(s, b'A', 8), None);
    }

    #[test]
    fn first_run_at_least_tolerates_isolated_mismatch() {
        let s = b"GGAAAAACAAAAAGG";
        let r = first_run_at_least(s, b'A', 10).unwrap();
        assert_eq!(r.start, 2);
        assert_eq!(r.len, 11);
    }

    #[test]
    fn first_tail_run_needs_ninety_percent_base() {
        // DZEC04 transcript sequence that first_run_at_least cut at.
        for genomic in [
            &b"GCCTACCATGAGAATAAGAGAAAGAAAATGAAGATCAA"[..], // ALB
            b"CCATAGTCTAATAGAAAACAACCGAAACCAAATAAT",        // MT-CO3
            b"ACCATTGCTGAAGAGAGAGACAAAAGACTGGCGGCCA",       // RPS4X
        ] {
            assert!(first_run_at_least(genomic, b'A', 12).is_some());
            assert_eq!(first_tail_run(genomic, b'A', 12), None);
        }
        // Read-through: cDNA, then the dT tail with one error, then UMI.
        let tail = b"GATTACACCATGGCTTAAAAAAAAAAAAGAAAAAAAAAAAAAAACGTTAGC";
        assert_eq!(
            first_tail_run(tail, b'A', 12),
            Some(Run { start: 16, len: 28 })
        );
        // An A-rich lead-in does not move the cut while it keeps 90%.
        let lead = b"CCCCAGAAAAAAAAAAAAAAAAAAAAGG";
        assert_eq!(first_tail_run(lead, b'A', 12).unwrap().start, 4);
        // Exactly one mismatch in 10 is allowed; 2 in 12 is not.
        assert!(first_tail_run(b"AAAAAGAAAA", b'A', 10).is_some());
        assert!(first_tail_run(b"AAGAAAAGAAAA", b'A', 12).is_none());
    }
}
