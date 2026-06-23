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

/// R2 primer — immediately precedes cDNA on the 3' side.
pub const R2_PRIMER: &[u8] = b"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT";

/// P7 library end.
pub const P7: &[u8] = b"CAAGCAGAAGACGGCATACGAGAT";

/// Template-switch oligo (5' assay spacer).
pub const TSO: &[u8] = b"AAGCAGTGGTATCAACGCAGAGTACATGGG";

/// Reverse-complement a DNA sequence (IUPAC ACGTN, case-insensitive -> upper).
pub fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
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
