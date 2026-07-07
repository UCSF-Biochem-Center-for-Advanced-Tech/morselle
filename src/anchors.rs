//! Adapter / leader sequences — the single source of truth (plan §4).
//!
//! NOTE: exact/full SBX leader & adapter set must be confirmed with Roche
//! (deck slide 14 Q1). These are the working values from the data analysis.

/// SBX leader, 5' of P5/P7 on intact reads (~40-45% of reads carry it).
pub const SBX_LEADER: &[u8] = b"AACAAT";

/// P5 library end.
pub const P5: &[u8] = b"AATGATACGGCGACCACCGAGAT";

/// R1 primer — immediately precedes the CB+UMI block.
pub const R1_PRIMER: &[u8] = b"ACACTCTTTCCCTACACGACGCTCTTCCGATCT";

/// 3' tail of the R1 primer (its last 10 bp). Sits immediately 5' of the CB+UMI
/// block, so a fixed-offset Hamming check of this short word cheaply confirms the
/// homopolymer-anchored bracket without a full adapter alignment.
pub const R1_PRIMER_TAIL: &[u8] = b"CTTCCGATCT";

/// R2 primer — immediately precedes cDNA on the 3' side.
pub const R2_PRIMER: &[u8] = b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";

/// P7 library end.
pub const P7: &[u8] = b"CAAGCAGAAGACGGCATACGAGAT";

/// Template-switch oligo (5' assay spacer).
pub const TSO: &[u8] = b"AAGCAGTGGTATCAACGCAGAGTACATGGG";

/// Universal Illumina 3' read-through adapter. In a canonical P5→P7 read this is
/// the start of `revcomp(R2 primer)` and marks the END of the cDNA (everything
/// from here — rc(R2 primer)/rc(i7)/rc(P7) — is trimmed). A cheap Hamming seed.
pub const ILLUMINA_3P_ADAPTER: &[u8] = b"AGATCGGAAGAGC";

/// Hamming distance between two equal-length slices; counts length mismatch as
/// the extra positions so callers can pass a window of the same length.
pub fn hamming(a: &[u8], b: &[u8]) -> u32 {
    let common = a.len().min(b.len());
    let mut d = (a.len().max(b.len()) - common) as u32;
    for i in 0..common {
        if a[i] != b[i] {
            d += 1;
        }
    }
    d
}

/// True iff the Hamming distance between `a` and `b` is ≤ `max`, early-exiting as
/// soon as the budget is exceeded. For short seeds scanned across a long read (the
/// 3' adapter trim), most positions blow the budget in the first few bases, so this
/// is much cheaper than summing the full distance and comparing. Same result as
/// `hamming(a, b) <= max`.
pub fn hamming_leq(a: &[u8], b: &[u8], max: u32) -> bool {
    let common = a.len().min(b.len());
    let mut d = (a.len().max(b.len()) - common) as u32;
    if d > max {
        return false;
    }
    for i in 0..common {
        if a[i] != b[i] {
            d += 1;
            if d > max {
                return false;
            }
        }
    }
    true
}

/// Complement one DNA base (IUPAC ACGTN, case-insensitive → upper); non-ACGT → `N`.
#[inline]
pub fn comp(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        b'T' | b't' => b'A',
        _ => b'N',
    }
}

/// Reverse-complement a DNA sequence (IUPAC ACGTN, case-insensitive -> upper).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| comp(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_roundtrip() {
        assert_eq!(revcomp(b"AACAAT"), b"ATTGTT");
        assert_eq!(revcomp(&revcomp(P5)), P5);
    }

    #[test]
    fn revcomp_handles_n() {
        assert_eq!(revcomp(b"ACGTN"), b"NACGT");
    }
}
