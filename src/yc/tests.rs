use super::codec::{decode, encode_overlap};
use super::types::{code_bases, is_reserved, pair_code, BaseType, YcError, YcTag, DECODE_R1, DECODE_R2, COMPLEMENT, GAP};
use super::types::complement_code;
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

/// Base complement, with gaps passing through.
fn comp(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    }
}

fn parse(text: &str) -> YcTag {
    YcTag::parse(text.as_bytes()).unwrap()
}

/// Every code that actually stands for a strand disagreement.
fn usable_codes() -> Vec<u8> {
    (b'A'..=b'Z')
        .filter(|&c| !is_reserved(c))
        .filter(|&c| {
            let (r1, r2) = code_bases(c).unwrap();
            r1 != r2
        })
        .collect()
}

#[test]
fn tables_are_well_formed() {
    assert_eq!(DECODE_R1.len(), 26);
    assert_eq!(DECODE_R2.len(), 26);
    assert_eq!(COMPLEMENT.len(), 26);
    // The four identity slots really are identities, so they can never be
    // produced by `pair_code` and never encode a disagreement.
    for base in [b'A', b'C', b'G', b'T'] {
        assert_eq!(
            code_bases(base),
            Some((base, base)),
            "slot {}",
            base as char
        );
    }
    // Reserved slots are gap/gap.
    for reserved in [b'N', b'U'] {
        let index = (reserved - b'A') as usize;
        assert_eq!(DECODE_R1[index], GAP);
        assert_eq!(DECODE_R2[index], GAP);
        assert!(code_bases(reserved).is_none());
    }
}

#[test]
fn published_code_table_matches_documented_matrix() {
    // The matrix published in the XOOS docs, rows R1, columns R2.
    let expected: &[(u8, u8, u8)] = &[
        (b'A', b'C', b'M'),
        (b'A', b'G', b'R'),
        (b'A', b'T', b'W'),
        (b'A', GAP, b'I'),
        (b'C', b'A', b'B'),
        (b'C', b'G', b'S'),
        (b'C', b'T', b'Y'),
        (b'C', GAP, b'P'),
        (b'G', b'A', b'D'),
        (b'G', b'C', b'V'),
        (b'G', b'T', b'K'),
        (b'G', GAP, b'J'),
        (b'T', b'A', b'H'),
        (b'T', b'C', b'E'),
        (b'T', b'G', b'F'),
        (b'T', GAP, b'X'),
        (GAP, b'A', b'L'),
        (GAP, b'C', b'Q'),
        (GAP, b'G', b'O'),
        (GAP, b'T', b'Z'),
    ];
    for &(r1, r2, code) in expected {
        assert_eq!(
            pair_code(r1, r2),
            Some(code),
            "{}->{} should encode as {}",
            r1 as char,
            r2 as char,
            code as char
        );
        assert_eq!(code_bases(code), Some((r1, r2)), "code {}", code as char);
    }
    // The matrix is exhaustive: 20 usable codes, one per off-diagonal cell.
    assert_eq!(usable_codes().len(), expected.len());
}

#[test]
fn complement_table_complements_the_bases() {
    for code in usable_codes() {
        let (r1, r2) = code_bases(code).unwrap();
        let complemented = complement_code(code);
        assert_eq!(
            code_bases(complemented),
            Some((comp(r1), comp(r2))),
            "complement of {}",
            code as char
        );
        // ...and it is an involution, in both cases.
        assert_eq!(complement_code(complemented), code);
        let lower = code.to_ascii_lowercase();
        assert!(complement_code(lower).is_ascii_lowercase());
        assert_eq!(complement_code(complement_code(lower)), lower);
    }
}

#[test]
fn encodes_published_example_1() {
    let (tag, stats) =
        encode_overlap(4, true, b"ACGTACGGTCATG", b"ACGTACGGTCATG", 3, false).unwrap();
    assert_eq!(tag.to_string(), "4+13-3");
    assert_eq!(stats.n_columns_coerced, 0);
}

#[test]
fn encodes_published_example_2() {
    let (tag, _) = encode_overlap(0, true, b"CAGAAG-", b"CTTA-GA", 0, true).unwrap();
    assert_eq!(tag.to_string(), "+1WK1I1L+");
}

#[test]
fn encodes_published_example_3() {
    let (tag, _) = encode_overlap(0, true, b"CGACGC-", b"CTCCGCA", 2, true).unwrap();
    assert_eq!(tag.to_string(), "+1KM3L+2");
}

#[test]
fn encodes_published_example_4() {
    let (tag, _) = encode_overlap(9, true, b"", b"", 0, true).unwrap();
    assert_eq!(tag.to_string(), "9++");
    assert!(tag.is_simplex());
}

#[test]
fn decodes_published_example_1() {
    let tag = parse("4+13-3");
    let seq = b"CATGACGTACGGTCATGGCG";
    let decoded = decode(&tag, seq, &vec![b'I'; seq.len()]).unwrap();
    assert_eq!(decoded.r1_seq, b"CATGACGTACGGTCATG");
    assert_eq!(decoded.r2_seq, b"ACGTACGGTCATGGCG");
}

#[test]
fn decodes_published_example_2() {
    let tag = parse("+1WK1I1L+");
    // Concordant positions (0, 3, 5) must match; the rest are the caller's
    // choice and do not affect what the strands decode to.
    let seq = b"CAGAAGA";
    let decoded = decode(&tag, seq, &vec![b'I'; seq.len()]).unwrap();
    assert_eq!(decoded.r1_seq, b"CAGAAG");
    assert_eq!(decoded.r2_seq, b"CTTAGA");
    assert_eq!(tag.count_discordant(), 4);
}

#[test]
fn decodes_published_example_3() {
    let tag = parse("+1KM3L+2");
    let seq = b"CGACGCAGT";
    let decoded = decode(&tag, seq, &vec![b'I'; seq.len()]).unwrap();
    assert_eq!(decoded.r1_seq, b"CGACGCGT");
    assert_eq!(decoded.r2_seq, b"CTCCGCA");
}

#[test]
fn decodes_published_example_4() {
    let tag = parse("9++");
    let seq = b"CGACGCACG";
    let decoded = decode(&tag, seq, &vec![b'I'; seq.len()]).unwrap();
    assert_eq!(decoded.r1_seq, seq);
    assert!(decoded.r2_seq.is_empty());
}

#[test]
fn published_examples_round_trip_as_strings() {
    for text in ["4+13-3", "+1WK1I1L+", "+1KM3L+2", "9++"] {
        assert_eq!(parse(text).to_string(), text);
    }
}

#[test]
fn reports_span_and_provenance() {
    let tag = parse("4+13-3");
    assert_eq!(tag.sequence_length(), 20);
    assert_eq!(tag.count_duplex(), 13);
    assert_eq!(tag.count_discordant(), 0);
    assert!(tag.is_duplex());

    let types = tag.base_types();
    assert_eq!(types.len(), 20);
    assert!(types[..4].iter().all(|t| *t == BaseType::Simplex));
    assert!(types[4..17].iter().all(|t| *t == BaseType::Concordant));
    assert!(types[17..].iter().all(|t| *t == BaseType::Simplex));

    assert_eq!(
        parse("+1WK1I1L+").base_types(),
        vec![
            BaseType::Concordant,
            BaseType::Discordant,
            BaseType::Discordant,
            BaseType::Concordant,
            BaseType::Discordant,
            BaseType::Concordant,
            BaseType::Discordant,
        ]
    );
}

#[test]
fn lowercase_codes_span_no_consensus_position() {
    // Consensus has a gap: R1 keeps its A, R2 its C, neither costs a
    // consensus position, and both decode at Q0.
    let tag = parse("+1m1+");
    assert_eq!(tag.sequence_length(), 2);
    assert_eq!(tag.count_discordant(), 1);
    assert_eq!(
        tag.base_types(),
        vec![BaseType::Concordant, BaseType::Concordant]
    );
    let decoded = decode(&tag, b"GT", b"II").unwrap();
    assert_eq!(decoded.r1_seq, b"GAT");
    assert_eq!(decoded.r1_qual, b"I!I");
    assert_eq!(decoded.r2_seq, b"GCT");
    assert_eq!(decoded.r2_qual, b"I!I");
}

#[test]
fn reverse_complement_is_an_involution() {
    for text in ["4+13-3", "+1WK1I1L+", "+1KM3L+2", "9++", "+1m1+"] {
        let original = parse(text);
        let mut flipped = original.clone();
        flipped.reverse_complement();
        flipped.reverse_complement();
        assert_eq!(flipped, original, "{text}");
    }
}

#[test]
fn reverse_complement_swaps_tails_and_complements_codes() {
    let mut tag = parse("4+13-3");
    tag.reverse_complement();
    // The R1 tail was on the left and is now on the right, and vice versa.
    assert_eq!(tag.to_string(), "3-13+4");

    let mut tag = parse("+1WK1I1L+");
    tag.reverse_complement();
    // Reversed order, each code complemented: L->Z, I->X, K->B, W->H.
    assert_eq!(tag.to_string(), "+Z1X1BH1+");
}

#[test]
fn reverse_complement_agrees_with_decoding_the_flipped_consensus() {
    let tag = parse("4+13-3");
    let seq: &[u8] = b"CATGACGTACGGTCATGGCG";
    let qual = vec![b'I'; seq.len()];
    let forward = decode(&tag, seq, &qual).unwrap();

    let rc_seq: Vec<u8> = seq.iter().rev().map(|&b| comp(b)).collect();
    let mut rc_tag = tag.clone();
    rc_tag.reverse_complement();
    let reverse = decode(&rc_tag, &rc_seq, &qual).unwrap();

    let expect_r1: Vec<u8> = forward.r1_seq.iter().rev().map(|&b| comp(b)).collect();
    let expect_r2: Vec<u8> = forward.r2_seq.iter().rev().map(|&b| comp(b)).collect();
    assert_eq!(reverse.r1_seq, expect_r1);
    assert_eq!(reverse.r2_seq, expect_r2);
}

#[test]
fn trims_within_the_tails() {
    let mut tag = parse("4+13-3");
    tag.trim(2, 1);
    assert_eq!(tag.to_string(), "2+13-2");
    assert_eq!(tag.sequence_length(), 17);
}

#[test]
fn trims_into_the_duplex_segment() {
    let mut tag = parse("4+13-3");
    tag.trim(6, 0);
    assert_eq!(tag.to_string(), "+11-3");

    let mut tag = parse("4+13-3");
    tag.trim(0, 5);
    assert_eq!(tag.to_string(), "4+11+");
}

#[test]
fn trimming_everything_empties_the_tag() {
    let mut tag = parse("4+13-3");
    tag.trim(100, 100);
    assert_eq!(tag.sequence_length(), 0);
    assert!(tag.is_simplex());
}

#[test]
fn trimming_drops_zero_length_ops_at_the_boundary() {
    // The leading `m` describes a base in the region being trimmed away.
    let mut tag = parse("+m1W1+");
    assert_eq!(tag.sequence_length(), 3);
    tag.trim(1, 0);
    assert_eq!(tag.to_string(), "+W1+");
}

#[test]
fn rejects_malformed_tags() {
    let cases: &[(&str, YcError)] = &[
        ("13", YcError::DelimiterCount(0)),
        ("1+2+3+", YcError::DelimiterCount(3)),
        ("+0+", YcError::ZeroRun),
        ("+01+", YcError::LeadingZero("01".into())),
        ("04+13-3", YcError::LeadingZero("04".into())),
        ("4+13-3x", YcError::NonNumericTail("3x".into())),
        ("+N+", YcError::UnsupportedCode('N')),
        ("+1U1+", YcError::UnsupportedCode('U')),
    ];
    for (text, expected) in cases {
        assert_eq!(
            YcTag::parse(text.as_bytes()).unwrap_err(),
            *expected,
            "input {text:?}"
        );
    }
}

#[test]
fn rejects_decoding_against_the_wrong_consensus() {
    let tag = parse("4+13-3");
    assert_eq!(
        decode(&tag, b"ACGT", b"IIII").unwrap_err(),
        YcError::LengthMismatch { tag: 20, seq: 4 }
    );
    assert_eq!(
        decode(&tag, b"ACGT", b"II").unwrap_err(),
        YcError::SeqQualMismatch { seq: 4, qual: 2 }
    );
}

#[test]
fn rejects_ragged_and_doubly_gapped_alignments() {
    assert_eq!(
        encode_overlap(0, true, b"ACG", b"AC", 0, true).unwrap_err(),
        YcError::RaggedAlignment { r1: 3, r2: 2 }
    );
    assert_eq!(
        encode_overlap(0, true, b"A-G", b"A-G", 0, true).unwrap_err(),
        YcError::DoubleGap(1)
    );
}

#[test]
fn coerces_n_columns_to_concordant() {
    let (tag, stats) = encode_overlap(0, true, b"ACNT", b"ACGT", 0, true).unwrap();
    assert_eq!(tag.to_string(), "+4+");
    assert_eq!(stats.n_columns_coerced, 1);
}

#[test]
fn gapless_alignments_round_trip_losslessly() {
    let bases = b"ACGT";
    let mut rng = ChaCha8Rng::seed_from_u64(0x5BC0DE);

    for _ in 0..2000 {
        let left = rng.gen_range(0..12u64);
        let right = rng.gen_range(0..12u64);
        let width = rng.gen_range(0..40usize);
        let left_is_r1 = rng.gen_bool(0.5);
        let right_is_r1 = rng.gen_bool(0.5);

        let random_bases = |rng: &mut ChaCha8Rng, n: usize| -> Vec<u8> {
            (0..n).map(|_| bases[rng.gen_range(0..4)]).collect()
        };

        let r1_duplex = random_bases(&mut rng, width);
        // Mostly agree, so runs and codes both get exercised.
        let r2_duplex: Vec<u8> = r1_duplex
            .iter()
            .map(|&b| {
                if rng.gen_bool(0.2) {
                    bases[rng.gen_range(0..4)]
                } else {
                    b
                }
            })
            .collect();
        let left_tail = random_bases(&mut rng, left as usize);
        let right_tail = random_bases(&mut rng, right as usize);

        let (tag, stats) =
            encode_overlap(left, left_is_r1, &r1_duplex, &r2_duplex, right, right_is_r1)
                .unwrap();
        assert_eq!(stats.n_columns_coerced, 0);

        // Consensus takes R1 at disagreements; the tag is lossless either way.
        let mut consensus = Vec::new();
        consensus.extend_from_slice(&left_tail);
        consensus.extend_from_slice(&r1_duplex);
        consensus.extend_from_slice(&right_tail);
        let qual = vec![b'I'; consensus.len()];

        assert_eq!(tag.sequence_length(), consensus.len() as u64);
        assert_eq!(YcTag::parse(tag.to_string().as_bytes()).unwrap(), tag);

        let decoded = decode(&tag, &consensus, &qual).unwrap();

        let mut expect_r1 = Vec::new();
        if left_is_r1 || left == 0 {
            expect_r1.extend_from_slice(&left_tail);
        }
        expect_r1.extend_from_slice(&r1_duplex);
        if right_is_r1 || right == 0 {
            expect_r1.extend_from_slice(&right_tail);
        }

        let mut expect_r2 = Vec::new();
        if !(left_is_r1 || left == 0) {
            expect_r2.extend_from_slice(&left_tail);
        }
        expect_r2.extend_from_slice(&r2_duplex);
        if !(right_is_r1 || right == 0) {
            expect_r2.extend_from_slice(&right_tail);
        }

        assert_eq!(decoded.r1_seq, expect_r1, "tag {tag}");
        assert_eq!(decoded.r2_seq, expect_r2, "tag {tag}");
    }
}
