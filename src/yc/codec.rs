use super::types::{
    code_bases, pair_code, DecodedPair, DuplexOp, EncodeStats, YcError, YcTag, CONSENSUS_GAP_QUAL,
    GAP,
};

pub fn encode_overlap(
    left_overhang: u64,
    left_is_r1: bool,
    r1_duplex: &[u8],
    r2_duplex: &[u8],
    right_overhang: u64,
    right_is_r1: bool,
) -> Result<(YcTag, EncodeStats), YcError> {
    if r1_duplex.len() != r2_duplex.len() {
        return Err(YcError::RaggedAlignment {
            r1: r1_duplex.len(),
            r2: r2_duplex.len(),
        });
    }
    let mut stats = EncodeStats::default();
    let mut duplex = Vec::new();
    let mut run = 0u64;

    for index in 0..r1_duplex.len() {
        let r1 = r1_duplex[index].to_ascii_uppercase();
        let r2 = r2_duplex[index].to_ascii_uppercase();
        if r1 == GAP && r2 == GAP {
            return Err(YcError::DoubleGap(index));
        }
        if r1 == b'N' || r2 == b'N' {
            stats.n_columns_coerced += 1;
            run += 1;
            continue;
        }
        if r1 == r2 {
            run += 1;
            continue;
        }
        if run > 0 {
            duplex.push(DuplexOp::Concordant(run));
            run = 0;
        }
        let code = pair_code(r1, r2).ok_or(YcError::UnsupportedPair {
            index,
            r1: r1 as char,
            r2: r2 as char,
        })?;
        duplex.push(DuplexOp::Discordant(code));
    }
    if run > 0 {
        duplex.push(DuplexOp::Concordant(run));
    }

    Ok((
        YcTag {
            left_overhang,
            left_is_r1: if left_overhang == 0 { true } else { left_is_r1 },
            duplex,
            right_overhang,
            right_is_r1: if right_overhang == 0 {
                true
            } else {
                right_is_r1
            },
        },
        stats,
    ))
}

pub fn decode(tag: &YcTag, seq: &[u8], qual: &[u8]) -> Result<DecodedPair, YcError> {
    if seq.len() != qual.len() {
        return Err(YcError::SeqQualMismatch {
            seq: seq.len(),
            qual: qual.len(),
        });
    }
    let span = tag.sequence_length();
    if span != seq.len() as u64 {
        return Err(YcError::LengthMismatch {
            tag: span,
            seq: seq.len(),
        });
    }

    let mut out = DecodedPair {
        r1_seq: Vec::new(),
        r1_qual: Vec::new(),
        r2_seq: Vec::new(),
        r2_qual: Vec::new(),
    };
    let mut pos = 0usize;

    let push_tail = |out: &mut DecodedPair, pos: &mut usize, len: u64, is_r1: bool| {
        for _ in 0..len {
            if is_r1 {
                out.r1_seq.push(seq[*pos]);
                out.r1_qual.push(qual[*pos]);
            } else {
                out.r2_seq.push(seq[*pos]);
                out.r2_qual.push(qual[*pos]);
            }
            *pos += 1;
        }
    };

    push_tail(&mut out, &mut pos, tag.left_overhang, tag.left_is_r1);

    for op in &tag.duplex {
        match op {
            DuplexOp::Concordant(n) => {
                for _ in 0..*n {
                    out.r1_seq.push(seq[pos]);
                    out.r1_qual.push(qual[pos]);
                    out.r2_seq.push(seq[pos]);
                    out.r2_qual.push(qual[pos]);
                    pos += 1;
                }
            }
            DuplexOp::Discordant(code) => {
                let (r1, r2) = code_bases(*code).ok_or(YcError::UnsupportedCode(*code as char))?;
                let (base_qual, consumed) = if code.is_ascii_uppercase() {
                    (qual[pos], 1)
                } else {
                    (CONSENSUS_GAP_QUAL, 0)
                };
                if r1 != GAP {
                    out.r1_seq.push(r1);
                    out.r1_qual.push(base_qual);
                }
                if r2 != GAP {
                    out.r2_seq.push(r2);
                    out.r2_qual.push(base_qual);
                }
                pos += consumed;
            }
        }
    }

    push_tail(&mut out, &mut pos, tag.right_overhang, tag.right_is_r1);
    Ok(out)
}
