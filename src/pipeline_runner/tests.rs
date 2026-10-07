#[cfg(test)]
mod tests {
    use crate::anchors::revcomp;
    use crate::pipeline_runner::formatting::{
        append_processing_tags, minimap2_header, push_fastq, push_fastq_r2_rc,
    };
    use crate::pipeline_runner::types::ProcessingTags;

    /// The optimized R2 emitter must produce byte-identical output to the original
    /// `revcomp(seq[start..end])` + reversed-qual + truncate approach it replaced.
    fn reference_r2(
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
        start: usize,
        end: usize,
        cap: usize,
    ) -> Vec<u8> {
        let mut r2_seq = revcomp(&seq[start..end]);
        let mut r2_qual = qual[start..end].to_vec();
        r2_qual.reverse();
        if cap > 0 && r2_seq.len() > cap {
            r2_seq.truncate(cap);
            r2_qual.truncate(cap);
        }
        let mut buf = Vec::new();
        push_fastq(&mut buf, id, &r2_seq, &r2_qual);
        buf
    }

    #[test]
    fn r2_rc_matches_reference() {
        let seq = b"ACGTACGTNNGGCCTTAAACGTACGTTTGGCCAA";
        let qual: Vec<u8> = (0..seq.len() as u8).map(|i| b'!' + (i % 40)).collect();
        let id = b"read/1";
        // Sweep sub-ranges and caps (including 0 = uncapped and caps > length).
        for start in 0..6 {
            for end in (start + 1)..=seq.len() {
                for cap in [0usize, 1, 5, 10, end - start, end - start + 5] {
                    let mut got = Vec::new();
                    let out = push_fastq_r2_rc(&mut got, id, seq, &qual, start, end, cap);
                    let want = reference_r2(id, seq, &qual, start, end, cap);
                    assert_eq!(got, want, "start={start} end={end} cap={cap}");
                    let expected_len = if cap > 0 {
                        (end - start).min(cap)
                    } else {
                        end - start
                    };
                    assert_eq!(out, expected_len);
                }
            }
        }
    }

    #[test]
    fn minimap2_header_uses_10x_style_tags() {
        let got = minimap2_header(
            b"read_id_001",
            b"AACCGGTT",
            b"AACCGGTT-1",
            b"ACTG",
            b"ACTG",
            b"IIIIIIII",
            b"JJJJ",
        );
        assert_eq!(
            got,
            b"read_id_001\tCR:Z:AACCGGTT\tCB:Z:AACCGGTT-1\tUR:Z:ACTG\tUB:Z:ACTG\tCY:Z:IIIIIIII\tUY:Z:JJJJ"
        );
    }

    #[test]
    fn processing_tags_are_valid_sam_style_fastq_comments() {
        let mut got = b"read_id_001".to_vec();
        append_processing_tags(
            &mut got,
            ProcessingTags {
                orientation: b'5',
                barcode_route: 3,
                index_route: 1,
                anchor_bits: 15,
                extracted_cdna_len: 427,
            },
        );
        assert_eq!(
            got,
            b"read_id_001\tXO:A:5\tXB:i:3\tXI:i:1\tXA:i:15\tXL:i:427"
        );
    }
}
