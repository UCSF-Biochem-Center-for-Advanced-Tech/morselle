use std::fmt;
use super::types::{
    code_bases, complement_code, parse_u64, BaseType, DuplexOp, YcError, YcTag, DELIM_R1, DELIM_R2,
};

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
