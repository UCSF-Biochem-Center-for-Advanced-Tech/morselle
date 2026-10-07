use thiserror::Error;

/// Gap placeholder in an aligned sequence column (not a tag delimiter).
pub const GAP: u8 = b'-';

/// Delimiter marking a tail contributed by R1.
pub const DELIM_R1: u8 = b'+';
/// Delimiter marking a tail contributed by R2.
pub const DELIM_R2: u8 = b'-';

/// Quality assigned to a base that exists in R1/R2 but not in the consensus.
pub const CONSENSUS_GAP_QUAL: u8 = b'!';

/// R1 base for each code, indexed by `code - 'A'`.
pub const DECODE_R1: &[u8; 26] = b"ACCGTTGTAGG-A--C-ACT-GATC-";
/// R2 base for each code, indexed by `code - 'A'`.
pub const DECODE_R2: &[u8; 26] = b"AACACGGA--TAC-G-CGGT-CT-TT";
/// Code whose `(R1, R2)` pair is the base-complement of this code's, same index.
pub const COMPLEMENT: &[u8; 26] = b"TKGYRMCWXPBZF-QJOEVA-SHIDL";

/// Codes that are structurally unusable: `N`/`U` occupy gap-gap slots.
pub fn is_reserved(upper: u8) -> bool {
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

pub fn code_bases(code: u8) -> Option<(u8, u8)> {
    let upper = code.to_ascii_uppercase();
    if !upper.is_ascii_uppercase() || is_reserved(upper) {
        return None;
    }
    let index = (upper - b'A') as usize;
    Some((DECODE_R1[index], DECODE_R2[index]))
}

pub fn pair_code(r1: u8, r2: u8) -> Option<u8> {
    if r1 == r2 {
        return None;
    }
    (0..26u8)
        .filter(|i| !is_reserved(b'A' + i))
        .find(|&i| DECODE_R1[i as usize] == r1 && DECODE_R2[i as usize] == r2)
        .map(|i| b'A' + i)
}

pub fn complement_code(code: u8) -> u8 {
    if code.is_ascii_uppercase() {
        COMPLEMENT[(code - b'A') as usize]
    } else {
        COMPLEMENT[(code - b'a') as usize].to_ascii_lowercase()
    }
}

pub fn parse_u64(text: &[u8]) -> Result<u64, YcError> {
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

/// Counts describing choices the encoder had to make.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct EncodeStats {
    pub n_columns_coerced: u64,
}

/// The two strands reconstructed from a consensus read and its tag.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DecodedPair {
    pub r1_seq: Vec<u8>,
    pub r1_qual: Vec<u8>,
    pub r2_seq: Vec<u8>,
    pub r2_qual: Vec<u8>,
}
