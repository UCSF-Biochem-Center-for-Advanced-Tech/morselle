//! 3' flank of the CB + UMI on 10x Flex v2 reads.
//!
//! 3' GEX reads confirm the CB location with the dT run after the UMI. Flex
//! v2 R1 has no dT run: the UMI is followed by a 0-3 bp linker and the fixed
//! sequence `TTGCTAGGACCG`, then the 10 bp sample barcode (CAT06: 20 bp after
//! the CB start + 28 in about 95% of reads, one base late in about 5%). That
//! fixed sequence, found with up to one edit where it is expected, is the
//! Flex counterpart of the dT flank for `located_confident`.

/// Fixed sequence after the UMI (+ 0-3 bp linker) on Flex v2 R1.
pub(crate) const FLEX_V2_FIXED: &[u8] = b"TTGCTAGGACCG";
/// Linker bases allowed between the UMI and the fixed sequence.
const MAX_LINKER: usize = 3;
/// Extra slack each side for indels in the CB or UMI.
const SLACK: usize = 2;

/// Whether `FLEX_V2_FIXED` occurs (at most one edit) starting between
/// `cb_umi_end - SLACK` and `cb_umi_end + MAX_LINKER + SLACK`.
pub(crate) fn flex_v2_flank(seq: &[u8], cb_umi_end: usize) -> bool {
    let from = cb_umi_end.saturating_sub(SLACK);
    let to = (cb_umi_end + MAX_LINKER + SLACK + FLEX_V2_FIXED.len() + 1).min(seq.len());
    if to <= from {
        return false;
    }
    min_semiglobal(FLEX_V2_FIXED, &seq[from..to]) <= 1
}

/// Smallest edit distance of all of `pat` against any substring of `text`.
fn min_semiglobal(pat: &[u8], text: &[u8]) -> usize {
    // Column-wise DP over `pat`; row 0 is free (a match may start anywhere).
    let mut col: Vec<usize> = (0..=pat.len()).collect();
    let mut best = col[pat.len()];
    for &t in text {
        let mut diag = col[0];
        col[0] = 0;
        for i in 1..=pat.len() {
            let left = col[i];
            col[i] = (diag + usize::from(pat[i - 1] != t))
                .min(left + 1)
                .min(col[i - 1] + 1);
            diag = left;
        }
        best = best.min(col[pat.len()]);
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read(linker: &str, fixed: &[u8]) -> (Vec<u8>, usize) {
        let mut s = b"ACGTACGGTCAGCATCGACCATGAGCGA".to_vec(); // CB + UMI
        let end = s.len();
        s.extend_from_slice(linker.as_bytes());
        s.extend_from_slice(fixed);
        s.extend_from_slice(b"CTTATTCGAT");
        (s, end)
    }

    #[test]
    fn found_after_every_linker_length() {
        for linker in ["", "A", "AC", "ACG"] {
            let (s, end) = read(linker, FLEX_V2_FIXED);
            assert!(flex_v2_flank(&s, end), "{linker:?}");
        }
    }

    #[test]
    fn found_with_one_edit_or_a_shifted_umi() {
        let mut sub = FLEX_V2_FIXED.to_vec();
        sub[5] = b'A';
        let mut del = FLEX_V2_FIXED.to_vec();
        del.remove(4);
        let mut ins = FLEX_V2_FIXED.to_vec();
        ins.insert(6, b'C');
        for fixed in [&sub, &del, &ins] {
            let (s, end) = read("A", fixed);
            assert!(flex_v2_flank(&s, end));
        }
        // A one-base indel upstream moves it by one either way.
        let (s, end) = read("", FLEX_V2_FIXED);
        assert!(flex_v2_flank(&s, end + 1) && flex_v2_flank(&s, end - 1));
    }

    #[test]
    fn absent_or_two_edits_or_far_away_is_not_found() {
        let mut two = FLEX_V2_FIXED.to_vec();
        two[2] = b'A';
        two[8] = b'T';
        let (s, end) = read("", &two);
        assert!(!flex_v2_flank(&s, end));
        let (s, end) = read("ACGTACGT", FLEX_V2_FIXED);
        assert!(!flex_v2_flank(&s, end));
        let (s, end) = read("", b"GATTACAGATTA");
        assert!(!flex_v2_flank(&s, end));
        assert!(!flex_v2_flank(&s[..end + 4], end));
    }
}
