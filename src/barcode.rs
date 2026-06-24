//! Cell-barcode whitelist load + validate-gate (plan §4 step 5, M4).
//!
//! Chemistry note: this is 10x **Single Cell Multiome (ATAC + GEX)**, GEX side,
//! whose gel-bead barcodes come from `737K-arc-v1.txt.gz` (NOT any 3' GEX list).
//! Validated: Roche's known-good CBs hit it 73% as-is (0% reverse-complemented).
//!
//! Barcodes are 2-bit packed into a `u32` (16 bp = 32 bits) and held in a sorted
//! vec for cache-friendly binary search. CB status is exact / Hamming-1 neighbor
//! (49 probes) / fail. We do NOT rewrite barcodes — Cell Ranger ARC corrects.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;

use crate::pipeline::CbStatus;

/// 2-bit pack a 16 bp barcode into a `u32` (A=00,C=01,G=10,T=11). `None` on any
/// non-ACGT base or wrong length.
pub fn encode16(cb: &[u8]) -> Option<u32> {
    if cb.len() != 16 {
        return None;
    }
    let mut code = 0u32;
    for &b in cb {
        let two = match b {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => return None,
        };
        code = (code << 2) | two;
    }
    Some(code)
}

/// A sorted set of 2-bit-encoded 16 bp barcodes.
pub struct Whitelist {
    codes: Vec<u32>,
}

impl Whitelist {
    /// Load a single-column barcode list (`.txt` or `.txt.gz`).
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("cannot open whitelist {}", path.display()))?;
        let is_gz = path.extension().is_some_and(|e| e == "gz");
        let reader: Box<dyn BufRead> = if is_gz {
            Box::new(BufReader::new(MultiGzDecoder::new(file)))
        } else {
            Box::new(BufReader::new(file))
        };

        let mut codes = Vec::with_capacity(1 << 20);
        for line in reader.lines() {
            let line = line?;
            let bc = line.trim().as_bytes();
            if bc.is_empty() {
                continue;
            }
            // whitelist files may be single-column or `barcode<TAB>...`
            let bc = bc.split(|&c| c == b'\t' || c == b' ').next().unwrap_or(bc);
            if let Some(code) = encode16(bc) {
                codes.push(code);
            } else {
                anyhow::bail!(
                    "whitelist {} has a non-16bp/non-ACGT entry: {}",
                    path.display(),
                    line.trim()
                );
            }
        }
        codes.sort_unstable();
        codes.dedup();
        Ok(Whitelist { codes })
    }

    pub fn len(&self) -> usize {
        self.codes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.codes.is_empty()
    }

    fn contains_code(&self, code: u32) -> bool {
        self.codes.binary_search(&code).is_ok()
    }

    /// Classify a 16 bp barcode: exact whitelist hit, within Hamming-1 of a
    /// whitelist member, or fail. (`Indel` is assigned later by M5.)
    pub fn classify(&self, cb: &[u8]) -> CbStatus {
        let Some(code) = encode16(cb) else {
            return CbStatus::Fail;
        };
        if self.contains_code(code) {
            return CbStatus::OnWhitelist;
        }
        // 16 positions × 3 alternative 2-bit values = 48 neighbor probes.
        for pos in 0..16 {
            let shift = pos * 2;
            let cur = (code >> shift) & 0b11;
            for alt in 0u32..4 {
                if alt != cur {
                    let neighbor = (code & !(0b11 << shift)) | (alt << shift);
                    if self.contains_code(neighbor) {
                        return CbStatus::Within1;
                    }
                }
            }
        }
        CbStatus::Fail
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_roundtrip_and_reject() {
        assert_eq!(encode16(b"AAAAAAAAAAAAAAAA"), Some(0));
        assert_eq!(encode16(b"TTTTTTTTTTTTTTTT"), Some(u32::MAX));
        // distinct barcodes give distinct codes
        assert_ne!(encode16(b"ACGTACGTACGTACGT"), encode16(b"TGCATGCATGCATGCA"));
        assert_eq!(encode16(b"AAAAAAAAAAAAAAAN"), None); // N
        assert_eq!(encode16(b"AAAA"), None); // wrong length
    }

    fn wl(barcodes: &[&[u8]]) -> Whitelist {
        let mut codes: Vec<u32> = barcodes.iter().map(|b| encode16(b).unwrap()).collect();
        codes.sort_unstable();
        Whitelist { codes }
    }

    #[test]
    fn classify_exact_within1_fail() {
        let w = wl(&[b"AAACCCAAGAAACACT", b"TTTTTTTTTTTTTTTT"]);
        assert_eq!(w.classify(b"AAACCCAAGAAACACT"), CbStatus::OnWhitelist);
        // one substitution from the first member
        assert_eq!(w.classify(b"AAACCCAAGAAACACA"), CbStatus::Within1);
        // two substitutions -> fail
        assert_eq!(w.classify(b"AAACCCAAGAAACAGG"), CbStatus::Fail);
        // N -> fail
        assert_eq!(w.classify(b"NAACCCAAGAAACACT"), CbStatus::Fail);
    }
}
