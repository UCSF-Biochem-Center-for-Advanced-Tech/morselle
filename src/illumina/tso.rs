//! Template-switch oligo at the start of a non-merged R2.
//!
//! In 3' GEX the TSO is at the far end of the molecule from the CB/UMI and
//! dT. On a fragment that reaches the transcript's 5' end, read 2 starts in
//! the TSO, so finding it needs only R2 (no mate merging). On 1 M raw CAT02
//! M2_GEX R2 reads, 8.5% start with it: 6.4% with all 30 bases, the rest
//! with a suffix (the fragment end fell inside the TSO).
use crate::anchors::{hamming_leq, TSO};

/// Shortest TSO suffix that counts at the start of R2. Matches the extract
/// path's `MIN_PARTIAL_TSO`: shorter suffixes are too common by chance.
pub const MIN_R2_TSO: usize = 12;

/// Number of leading R2 bases that are TSO (the whole oligo or a suffix of
/// it, at most one mismatch per 10 bases), or `None`. Longest match wins.
pub fn r2_tso_prefix(r2: &[u8]) -> Option<usize> {
    let m = TSO.len();
    let longest = r2.len().min(m);
    (MIN_R2_TSO..=longest)
        .rev()
        .find(|&len| hamming_leq(&r2[..len], &TSO[m - len..], (len / 10) as u32))
}

#[cfg(test)]
mod tests {
    use super::*;

    const CDNA: &[u8] = b"GATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGC";

    fn with_prefix(prefix: &[u8]) -> Vec<u8> {
        let mut r2 = prefix.to_vec();
        r2.extend_from_slice(CDNA);
        r2
    }

    #[test]
    fn finds_the_whole_tso() {
        assert_eq!(r2_tso_prefix(&with_prefix(TSO)), Some(TSO.len()));
    }

    #[test]
    fn finds_a_tso_suffix() {
        let suffix = &TSO[TSO.len() - 17..];
        assert_eq!(r2_tso_prefix(&with_prefix(suffix)), Some(17));
    }

    #[test]
    fn tolerates_one_mismatch_per_ten_bases() {
        let mut tso = TSO.to_vec();
        tso[5] = b'T';
        tso[20] = b'C';
        assert_eq!(r2_tso_prefix(&with_prefix(&tso)), Some(TSO.len()));
    }

    #[test]
    fn ignores_short_suffixes_and_plain_cdna() {
        assert_eq!(r2_tso_prefix(&with_prefix(&TSO[TSO.len() - 11..])), None);
        assert_eq!(r2_tso_prefix(CDNA), None);
        assert_eq!(r2_tso_prefix(b"ACATGGG"), None);
    }

    #[test]
    fn tso_later_in_the_read_is_not_a_prefix() {
        let mut r2 = b"ACGTAC".to_vec();
        r2.extend_from_slice(TSO);
        r2.extend_from_slice(CDNA);
        assert_eq!(r2_tso_prefix(&r2), None);
    }
}
