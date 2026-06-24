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
}
