//! Roche SBX `YC` duplex-consensus tag — encode and decode (plan: OVERLAP_PLAN.md §Phase 2).
//!
//! A `YC` tag rides along with a duplex consensus read (as a BAM auxiliary tag,
//! or as a FASTQ comment field) and encodes losslessly how the consensus differs
//! from the two strands it was built from. Given the consensus sequence, its
//! qualities, and the tag, both original strands can be reconstructed exactly.
//!
//! ```text
//! {left_tail}{delimiter}{duplex}{delimiter}{right_tail}
//! ```
//!
//! * `+` means that tail's simplex bases came from R1, `-` from R2. A tag always
//!   carries exactly two delimiters; `+` is the default when a tail is empty.
//! * Tails are bare integers — counts of simplex bases. Never letters.
//! * The duplex segment alternates integers (runs of concordant bases) with
//!   single-letter discordance codes.
//!
//! We use this format for Illumina paired-end overlap merges as well as SBX-D
//! reads: an overlap merge has exactly the same structure (two opposite strands
//! of one duplex fragment, a concordant middle, simplex tails).
//!
//! Written from the public specification at
//! <https://roche-axelios.gitbook.io/xoos/latest/tutorials/interpreting-sbx-reads/interpreting-sbx-d-reads>.
//! Deliberately NOT derived from the XOOS reference implementation, whose licence
//! restricts derivatives to data originating from Roche instruments.

use std::fmt;

use thiserror::Error;

/// Gap placeholder in an aligned sequence column (not a tag delimiter).
pub const GAP: u8 = b'-';

/// Delimiter marking a tail contributed by R1.
const DELIM_R1: u8 = b'+';
/// Delimiter marking a tail contributed by R2.
const DELIM_R2: u8 = b'-';

/// Quality assigned to a base that exists in R1/R2 but not in the consensus:
/// the consensus carries no quality for it, so it decodes as Q0.
const CONSENSUS_GAP_QUAL: u8 = b'!';

/// R1 base for each code, indexed by `code - 'A'`.
const DECODE_R1: &[u8; 26] = b"ACCGTTGTAGG-A--C-ACT-GATC-";
/// R2 base for each code, indexed by `code - 'A'`.
const DECODE_R2: &[u8; 26] = b"AACACGGA--TAC-G-CGGT-CT-TT";
/// Code whose `(R1, R2)` pair is the base-complement of this code's, same index.
const COMPLEMENT: &[u8; 26] = b"TKGYRMCWXPBZF-QJOEVA-SHIDL";

/// Codes that are structurally unusable: `N`/`U` occupy gap-gap slots.
fn is_reserved(upper: u8) -> bool {
    upper == b'N' || upper == b'U'
}

/// Per-position provenance of a consensus base.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BaseType {
    /// Both strands covered this position and agreed.
    Concordant,
    /// Both strands covered this position and disagreed.
    Discordant,
    /// Only one strand covered this position.
    Simplex,
}

/// One operator in the duplex segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DuplexOp {
    /// A run of `n` consecutive concordant positions. Never zero.
    Concordant(u64),
    /// A single discordant position, carrying its letter code.
    ///
    /// Uppercase: the consensus has a base here (consumes one consensus
    /// position). Lowercase: the consensus has a *gap* here (consumes none), and
    /// the reconstructed base takes quality [`CONSENSUS_GAP_QUAL`].
    Discordant(u8),
}

impl DuplexOp {
    /// How many consensus positions this operator spans.
    pub fn consensus_len(self) -> u64 {
        match self {
            DuplexOp::Concordant(n) => n,
            DuplexOp::Discordant(code) => {
                if code.is_ascii_uppercase() {
                    1
                } else {
                    0
                }
            }
        }
    }
}

/// A parsed `YC` tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct YcTag {
    /// Simplex bases at the 5' end.
    pub left_overhang: u64,
    /// True if the left tail came from R1 (`+`), false for R2 (`-`).
    pub left_is_r1: bool,
    /// Operators covering the region both strands reached.
    pub duplex: Vec<DuplexOp>,
    /// Simplex bases at the 3' end.
    pub right_overhang: u64,
    /// True if the right tail came from R1 (`+`), false for R2 (`-`).
    pub right_is_r1: bool,
}

impl Default for YcTag {
    fn default() -> Self {
        // `+` is the documented default delimiter for an empty tail.
        Self {
            left_overhang: 0,
            left_is_r1: true,
            duplex: Vec::new(),
            right_overhang: 0,
            right_is_r1: true,
        }
    }
}

/// Everything that can go wrong parsing, encoding, or decoding a tag.
#[derive(Debug, Clone, PartialEq, Eq, Error)]
pub enum YcError {
    #[error("YC tag must contain exactly two delimiters (`+` or `-`), found {0}")]
    DelimiterCount(usize),
    #[error("YC tag is not ASCII")]
    NotAscii,
    #[error("YC tail segment must contain digits only, found {0:?}")]
    NonNumericTail(String),
    #[error("YC integer {0:?} has a leading zero")]
    LeadingZero(String),
    #[error("YC concordant run must not be zero")]
    ZeroRun,
    #[error("YC integer {0:?} does not fit in u64")]
    IntegerOverflow(String),
    #[error("unsupported YC code {0:?}")]
    UnsupportedCode(char),
    #[error("cannot encode column {index}: R1 {r1:?} / R2 {r2:?} has no YC code")]
    UnsupportedPair { index: usize, r1: char, r2: char },
    #[error("cannot encode column {0}: both strands have a gap")]
    DoubleGap(usize),
    #[error("aligned R1 and R2 must be the same length ({r1} vs {r2})")]
    RaggedAlignment { r1: usize, r2: usize },
    #[error("YC tag spans {tag} positions but the consensus is {seq} long")]
    LengthMismatch { tag: u64, seq: usize },
    #[error("consensus sequence and quality lengths differ ({seq} vs {qual})")]
    SeqQualMismatch { seq: usize, qual: usize },
}

/// The `(R1, R2)` base pair a code stands for, or `None` if the code is not a
/// usable YC code. Case-insensitive: case encodes consensus presence, not bases.
fn code_bases(code: u8) -> Option<(u8, u8)> {
    let upper = code.to_ascii_uppercase();
    if !upper.is_ascii_uppercase() || is_reserved(upper) {
        return None;
    }
    let index = (upper - b'A') as usize;
    Some((DECODE_R1[index], DECODE_R2[index]))
}

/// The code for an `(R1, R2)` pair, or `None` if the pair has no encoding.
///
/// Concordant pairs deliberately return `None` — they are runs, not codes. The
/// four identity slots (`A`/`C`/`G`/`T`) are therefore never emitted.
fn pair_code(r1: u8, r2: u8) -> Option<u8> {
    if r1 == r2 {
        return None;
    }
    (0..26u8)
        .filter(|i| !is_reserved(b'A' + i))
        .find(|&i| DECODE_R1[i as usize] == r1 && DECODE_R2[i as usize] == r2)
        .map(|i| b'A' + i)
}

/// The code standing for the base-complement of `code`'s pair, preserving case.
fn complement_code(code: u8) -> u8 {
    if code.is_ascii_uppercase() {
        COMPLEMENT[(code - b'A') as usize]
    } else {
        COMPLEMENT[(code - b'a') as usize].to_ascii_lowercase()
    }
}

/// Parse a non-negative integer, rejecting leading zeros.
fn parse_u64(text: &[u8]) -> Result<u64, YcError> {
    let rendered = || String::from_utf8_lossy(text).into_owned();
    if text.len() > 1 && text[0] == b'0' {
        return Err(YcError::LeadingZero(rendered()));
    }
    let mut value: u64 = 0;
    for &byte in text {
        let digit = (byte - b'0') as u64;
        value = value
            .checked_mul(10)
            .and_then(|v| v.checked_add(digit))
            .ok_or_else(|| YcError::IntegerOverflow(rendered()))?;
    }
    Ok(value)
}

impl YcTag {
    /// Parse a tag from its wire form.
    ///
    /// Deliberately lenient about the four identity codes (`A`/`C`/`G`/`T`):
    /// we never emit them, but we accept them so that tags written by other
    /// implementations round-trip rather than erroring out.
    pub fn parse(text: &[u8]) -> Result<Self, YcError> {
        if !text.is_ascii() {
            return Err(YcError::NotAscii);
        }
        let delimiters: Vec<usize> = text
            .iter()
            .enumerate()
            .filter(|(_, &b)| b == DELIM_R1 || b == DELIM_R2)
            .map(|(i, _)| i)
            .collect();
        if delimiters.len() != 2 {
            return Err(YcError::DelimiterCount(delimiters.len()));
        }
        let (first, second) = (delimiters[0], delimiters[1]);

        let left_text = &text[..first];
        let duplex_text = &text[first + 1..second];
        let right_text = &text[second + 1..];

        let tail = |segment: &[u8]| -> Result<u64, YcError> {
            if segment.iter().any(|b| !b.is_ascii_digit()) {
                return Err(YcError::NonNumericTail(
                    String::from_utf8_lossy(segment).into_owned(),
                ));
            }
            parse_u64(segment)
        };

        let mut duplex = Vec::new();
        let mut digits_from: Option<usize> = None;
        for (index, &byte) in duplex_text.iter().enumerate() {
            if byte.is_ascii_digit() {
                digits_from.get_or_insert(index);
                continue;
            }
            if let Some(start) = digits_from.take() {
                let run = parse_u64(&duplex_text[start..index])?;
                if run == 0 {
                    return Err(YcError::ZeroRun);
                }
                duplex.push(DuplexOp::Concordant(run));
            }
            if code_bases(byte).is_none() {
                return Err(YcError::UnsupportedCode(byte as char));
            }
            duplex.push(DuplexOp::Discordant(byte));
        }
        if let Some(start) = digits_from {
            let run = parse_u64(&duplex_text[start..])?;
            if run == 0 {
                return Err(YcError::ZeroRun);
            }
            duplex.push(DuplexOp::Concordant(run));
        }

        Ok(Self {
            left_overhang: tail(left_text)?,
            left_is_r1: text[first] == DELIM_R1,
            duplex,
            right_overhang: tail(right_text)?,
            right_is_r1: text[second] == DELIM_R1,
        })
    }

    /// Total consensus positions this tag accounts for.
    pub fn sequence_length(&self) -> u64 {
        self.left_overhang
            + self.duplex.iter().map(|op| op.consensus_len()).sum::<u64>()
            + self.right_overhang
    }

    /// True if any part of the read was covered by both strands.
    pub fn is_duplex(&self) -> bool {
        !self.duplex.is_empty()
    }

    /// True if no part of the read was covered by both strands. An empty tag
    /// counts as simplex.
    pub fn is_simplex(&self) -> bool {
        self.duplex.is_empty()
    }

    /// Positions where the two strands disagreed. Counts zero-length (lowercase)
    /// operators too — the strands still disagreed there, the consensus just
    /// dropped the base.
    pub fn count_discordant(&self) -> u64 {
        self.duplex
            .iter()
            .filter(|op| matches!(op, DuplexOp::Discordant(_)))
            .count() as u64
    }

    /// Consensus positions covered by both strands, concordant and discordant.
    pub fn count_duplex(&self) -> u64 {
        self.duplex.iter().map(|op| op.consensus_len()).sum()
    }

    /// Provenance of every consensus position, in order.
    pub fn base_types(&self) -> Vec<BaseType> {
        let mut types = Vec::with_capacity(self.sequence_length() as usize);
        types.extend(std::iter::repeat_n(
            BaseType::Simplex,
            self.left_overhang as usize,
        ));
        for op in &self.duplex {
            match op {
                DuplexOp::Concordant(n) => {
                    types.extend(std::iter::repeat_n(BaseType::Concordant, *n as usize));
                }
                DuplexOp::Discordant(code) => {
                    if code.is_ascii_uppercase() {
                        types.push(BaseType::Discordant);
                    }
                }
            }
        }
        types.extend(std::iter::repeat_n(
            BaseType::Simplex,
            self.right_overhang as usize,
        ));
        types
    }

    /// Reverse-complement the tag, to accompany a reverse-complemented consensus.
    ///
    /// A tag cannot be reversed by reversing its string: the tails swap ends
    /// *and* every code must go through the complement table.
    pub fn reverse_complement(&mut self) {
        std::mem::swap(&mut self.left_overhang, &mut self.right_overhang);
        std::mem::swap(&mut self.left_is_r1, &mut self.right_is_r1);
        self.duplex.reverse();
        for op in &mut self.duplex {
            if let DuplexOp::Discordant(code) = op {
                *code = complement_code(*code);
            }
        }
    }

    /// Trim `left` consensus positions off the 5' end and `right` off the 3' end,
    /// as when an aligner soft-clips the read. Trimming more than the tag spans
    /// empties it.
    ///
    /// Zero-length (lowercase) operators sitting at a trimmed boundary are
    /// dropped: they describe bases in the trimmed-away region.
    pub fn trim(&mut self, left: u64, right: u64) {
        self.trim_left(left);
        self.trim_right(right);
        self.normalize();
    }

    /// Put the tag in canonical form: a tail with no bases carries no strand
    /// information, so it takes the documented default delimiter `+`.
    pub fn normalize(&mut self) {
        if self.left_overhang == 0 {
            self.left_is_r1 = true;
        }
        if self.right_overhang == 0 {
            self.right_is_r1 = true;
        }
    }

    fn trim_left(&mut self, mut remaining: u64) {
        let taken = remaining.min(self.left_overhang);
        self.left_overhang -= taken;
        remaining -= taken;

        let mut consumed = 0usize;
        while consumed < self.duplex.len() && remaining > 0 {
            let span = self.duplex[consumed].consensus_len();
            if span <= remaining {
                remaining -= span;
                consumed += 1;
            } else if let DuplexOp::Concordant(run) = &mut self.duplex[consumed] {
                *run = span - remaining;
                remaining = 0;
            }
        }
        self.duplex.drain(..consumed);
        while matches!(self.duplex.first(), Some(DuplexOp::Discordant(c)) if c.is_ascii_lowercase())
        {
            self.duplex.remove(0);
        }
        if remaining > 0 {
            // Trimmed past the duplex segment entirely.
            let taken = remaining.min(self.right_overhang);
            self.right_overhang -= taken;
        }
    }

    fn trim_right(&mut self, mut remaining: u64) {
        let taken = remaining.min(self.right_overhang);
        self.right_overhang -= taken;
        remaining -= taken;

        while !self.duplex.is_empty() && remaining > 0 {
            let last = self.duplex.len() - 1;
            let span = self.duplex[last].consensus_len();
            if span <= remaining {
                remaining -= span;
                self.duplex.pop();
            } else if let DuplexOp::Concordant(run) = &mut self.duplex[last] {
                *run = span - remaining;
                remaining = 0;
            }
        }
        while matches!(self.duplex.last(), Some(DuplexOp::Discordant(c)) if c.is_ascii_lowercase())
        {
            self.duplex.pop();
        }
        if remaining > 0 {
            let taken = remaining.min(self.left_overhang);
            self.left_overhang -= taken;
        }
    }
}

impl fmt::Display for YcTag {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.left_overhang > 0 {
            write!(f, "{}", self.left_overhang)?;
        }
        f.write_str(if self.left_is_r1 { "+" } else { "-" })?;
        for op in &self.duplex {
            match op {
                DuplexOp::Concordant(n) => write!(f, "{n}")?,
                DuplexOp::Discordant(code) => write!(f, "{}", *code as char)?,
            }
        }
        f.write_str(if self.right_is_r1 { "+" } else { "-" })?;
        if self.right_overhang > 0 {
            write!(f, "{}", self.right_overhang)?;
        }
        Ok(())
    }
}

/// Counts describing choices the encoder had to make.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EncodeStats {
    /// Columns where at least one strand read `N` and which were therefore
    /// recorded as concordant. See [`encode_overlap`].
    pub n_columns_coerced: u64,
}

/// Build a tag from an overlap alignment.
///
/// `r1_duplex` and `r2_duplex` are the aligned columns of the region both strands
/// reached, equal length, with [`GAP`] marking a gap. The tail lengths and their
/// strands are supplied by the caller rather than inferred, because the boundary
/// between a trailing tail and a trailing indel is genuinely ambiguous from the
/// alignment alone — the published `+1WK1I1L+` example encodes a trailing
/// one-sided column as an in-duplex code, not as a tail.
///
/// **`N` handling.** There is no code for a pair involving `N`, so a column where
/// either strand read `N` is recorded as concordant and counted in
/// [`EncodeStats::n_columns_coerced`]. This is the one place the tag is not
/// strictly lossless: decoding reconstructs the consensus base on both strands
/// rather than the original `N`. `N` carries no base information, so the
/// alternative — dropping the tag for any read with an `N` in its overlap — loses
/// strictly more.
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
            // An empty tail takes the documented default delimiter.
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

/// The two strands reconstructed from a consensus read and its tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedPair {
    pub r1_seq: Vec<u8>,
    pub r1_qual: Vec<u8>,
    pub r2_seq: Vec<u8>,
    pub r2_qual: Vec<u8>,
}

/// Reconstruct both strands from a consensus read and its tag.
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
                // Uppercase: the consensus has a base here and lends its quality.
                // Lowercase: the consensus has a gap, so the bases decode at Q0.
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

#[cfg(test)]
mod tests {
    use super::*;
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

    // ── The four worked examples published in the XOOS documentation ─────────

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

    // ── Structure ───────────────────────────────────────────────────────────

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

    // ── Orientation ─────────────────────────────────────────────────────────

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

    // ── Trimming ────────────────────────────────────────────────────────────

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

    // ── Rejections ──────────────────────────────────────────────────────────

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

    // ── The test that actually proves the module ────────────────────────────

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
}
