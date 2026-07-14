//! Cell-barcode whitelist load + validate-gate (plan §4 step 5, M4).
//!
//! Chemistry note: this is 10x **Single Cell Multiome (ATAC + GEX)**, GEX side,
//! whose gel-bead barcodes come from `737K-arc-v1.txt.gz` (NOT any 3' GEX list).
//! Validated: Roche's known-good CBs hit it 73% as-is (0% reverse-complemented).
//!
//! Barcodes are 2-bit packed into a `u32` (16 bp = 32 bits) and held in an
//! `FxHashSet` for O(1) membership. CB status is exact / Hamming-1 neighbor
//! (49 probes) / fail; the fail path probes all 48 single-substitution neighbors,
//! which scatter across the code space, so O(1) lookups beat a binary search whose
//! ~log₂(N) comparisons each miss cache. FxHash (multiply+rotate) is used rather
//! than the slow default SipHash — with SipHash the set can lose to binary search
//! on `u32` keys. We do NOT rewrite barcodes — Cell Ranger ARC corrects.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use anyhow::{Context, Result};
use flate2::read::MultiGzDecoder;
use rustc_hash::FxHashSet;

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

/// Helper for `correct_within2`: record a whitelist hit, flagging ambiguity if a
/// second distinct barcode is seen.
fn consider_hit(wl: &Whitelist, code: u32, found: &mut Option<u32>, ambiguous: &mut bool) {
    if wl.contains_code(code) {
        match *found {
            None => *found = Some(code),
            Some(c) if c != code => *ambiguous = true,
            _ => {}
        }
    }
}

/// Decode a 2-bit-packed 16 bp barcode back to bytes.
pub fn decode16(mut code: u32) -> [u8; 16] {
    let mut out = [0u8; 16];
    for i in (0..16).rev() {
        out[i] = b"ACGT"[(code & 0b11) as usize];
        code >>= 2;
    }
    out
}

/// A set of 2-bit-encoded 16 bp barcodes with O(1) membership.
pub struct Whitelist {
    codes: FxHashSet<u32>,
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

        let mut codes = FxHashSet::with_capacity_and_hasher(1 << 20, rustc_hash::FxBuildHasher);
        for line in reader.lines() {
            let line = line?;
            let bc = line.trim().as_bytes();
            if bc.is_empty() {
                continue;
            }
            // whitelist files may be single-column or `barcode<TAB>...`
            let bc = bc.split(|&c| c == b'\t' || c == b' ').next().unwrap_or(bc);
            if let Some(code) = encode16(bc) {
                codes.insert(code); // the set dedups
            } else {
                anyhow::bail!(
                    "whitelist {} has a non-16bp/non-ACGT entry: {}",
                    path.display(),
                    line.trim()
                );
            }
        }
        Ok(Whitelist { codes })
    }

    pub fn len(&self) -> usize {
        self.codes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.codes.is_empty()
    }

    fn contains_code(&self, code: u32) -> bool {
        self.codes.contains(&code)
    }

    /// Whether a 16 bp barcode is exactly on the whitelist (diagnostics).
    pub fn contains(&self, cb: &[u8]) -> bool {
        encode16(cb).is_some_and(|c| self.contains_code(c))
    }

    /// All 2-bit-encoded codes, sorted (diagnostics: brute-force nearest-neighbor
    /// scans and iteration). Materializes a `Vec` — not for the hot path.
    pub fn codes_sorted(&self) -> Vec<u32> {
        let mut v: Vec<u32> = self.codes.iter().copied().collect();
        v.sort_unstable();
        v
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
        if self.correct_within1(cb).is_some() {
            CbStatus::Within1
        } else {
            CbStatus::Fail
        }
    }

    /// Lenient correction: return the unique whitelist barcode within Hamming-2
    /// of `cb` (exact/within-1 included), or `None` if there are zero or multiple.
    /// ~1080 probes — opt-in (`--lenient`) and only worth applying to confidently
    /// located CBs.
    pub fn correct_within2(&self, cb: &[u8]) -> Option<[u8; 16]> {
        let Some(code) = encode16(cb) else {
            return None;
        };
        if self.contains_code(code) {
            return Some(decode16(code));
        }
        let mut found: Option<u32> = None;
        let mut ambiguous = false;
        for i in 0..16 {
            let shift_i = i * 2;
            let cur_i = (code >> shift_i) & 0b11;
            for ai in 0u32..4 {
                if ai == cur_i {
                    continue;
                }
                let c1 = (code & !(0b11 << shift_i)) | (ai << shift_i);
                // distance-1
                consider_hit(self, c1, &mut found, &mut ambiguous);
                // distance-2 (only j > i to avoid double-visiting)
                for j in (i + 1)..16 {
                    let shift_j = j * 2;
                    let cur_j = (c1 >> shift_j) & 0b11;
                    for aj in 0u32..4 {
                        if aj == cur_j {
                            continue;
                        }
                        let c2 = (c1 & !(0b11 << shift_j)) | (aj << shift_j);
                        consider_hit(self, c2, &mut found, &mut ambiguous);
                    }
                }
                if ambiguous {
                    return None;
                }
            }
        }
        match (found, ambiguous) {
            (Some(c), false) => Some(decode16(c)),
            _ => None,
        }
    }

    /// Return the unique whitelist barcode within Hamming-1 of `cb` (exact
    /// included), or `None` if zero or multiple match. Cheaper sibling of
    /// `correct_within2` (49 probes); used by the CB-reframe rescue.
    pub fn correct_within1(&self, cb: &[u8]) -> Option<[u8; 16]> {
        let code = encode16(cb)?;
        if self.contains_code(code) {
            return Some(decode16(code));
        }
        let mut found: Option<u32> = None;
        for pos in 0..16 {
            let shift = pos * 2;
            let cur = (code >> shift) & 0b11;
            for alt in 0u32..4 {
                if alt != cur {
                    let nb = (code & !(0b11 << shift)) | (alt << shift);
                    if self.contains_code(nb) {
                        match found {
                            None => found = Some(nb),
                            Some(c) if c != nb => return None,
                            _ => {}
                        }
                    }
                }
            }
        }
        found.map(decode16)
    }

    /// CB reframe rescue: a single homopolymer indel in the CB/UMI block (the
    /// dominant SBX error) frame-shifts the naive 16 bp grab anchored at the rigid
    /// `cb_start`. Try frame offsets ordered +1 first (the deletion mode that
    /// dominates, CB sits one base 3′ of the rigid frame), then ±2, returning the
    /// corrected CB and the winning offset from the first frame that yields a
    /// UNIQUE within-1 whitelist hit. Caller should gate on confident location so
    /// the shift reflects a real indel, not a mis-anchor.
    pub fn reframe_cb(
        &self,
        seq: &[u8],
        cb_start: usize,
        window: usize,
    ) -> Option<([u8; 16], i64)> {
        for mag in 1..=window as i64 {
            for k in [mag, -mag] {
                let s = cb_start as i64 + k;
                if s < 0 || s as usize + 16 > seq.len() {
                    continue;
                }
                if let Some(cb) = self.correct_within1(&seq[s as usize..s as usize + 16]) {
                    return Some((cb, k));
                }
            }
        }
        None
    }

    /// M5 indel rescue: when the measured CB+UMI span is one off 28, the CB likely
    /// carries a single indel that frame-shifted the naive 16 bp grab. `region`
    /// starts at `cb_start`; `delta` = measured_span − expected_span (+1 insertion,
    /// −1 deletion). Returns the corrected 16 bp CB iff exactly one whitelist
    /// barcode is consistent with the single indel (ambiguous/none → `None`).
    pub fn realign_cb(&self, region: &[u8], delta: i32) -> Option<[u8; 16]> {
        let mut found: Option<u32> = None;
        let mut ambiguous = false;
        let mut consider = |code: u32| {
            if self.contains_code(code) {
                match found {
                    None => found = Some(code),
                    Some(c) if c != code => ambiguous = true,
                    _ => {}
                }
            }
        };

        match delta {
            // Insertion: drop one base from the first 17 to recover a 16-mer.
            1 => {
                if region.len() < 17 {
                    return None;
                }
                let r = &region[..17];
                let mut cb = [0u8; 16];
                for skip in 0..17 {
                    let mut j = 0;
                    for (i, &b) in r.iter().enumerate() {
                        if i != skip {
                            cb[j] = b;
                            j += 1;
                        }
                    }
                    if let Some(code) = encode16(&cb) {
                        consider(code);
                    }
                }
            }
            // Deletion: insert one base into the first 15 to recover a 16-mer.
            -1 => {
                if region.len() < 15 {
                    return None;
                }
                let r = &region[..15];
                let mut cb = [0u8; 16];
                for pos in 0..=15 {
                    cb[..pos].copy_from_slice(&r[..pos]);
                    cb[pos + 1..].copy_from_slice(&r[pos..]);
                    for &b in b"ACGT" {
                        cb[pos] = b;
                        if let Some(code) = encode16(&cb) {
                            consider(code);
                        }
                    }
                }
            }
            _ => return None,
        }

        match (found, ambiguous) {
            (Some(code), false) => Some(decode16(code)),
            _ => None,
        }
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
        let codes: FxHashSet<u32> = barcodes.iter().map(|b| encode16(b).unwrap()).collect();
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

    #[test]
    fn decode_roundtrip() {
        for bc in [
            b"AAACCCAAGAAACACT",
            b"ACGTACGTACGTACGT",
            b"TTTTTTTTTTTTTTTT",
        ] {
            assert_eq!(&decode16(encode16(bc).unwrap()), bc);
        }
    }

    #[test]
    fn realign_insertion() {
        let truth = b"AAACCCAAGAAACACT";
        let w = wl(&[truth]);
        // region = CB with one inserted base ('G' after position 4) + UMI tail
        // truth: AAACC|CAAGAAACACT -> insert G -> AAACCGCAAGAAACACT (17) + UMI...
        let region = b"AAACCGCAAGAAACACTNNNNNNNNNNNN";
        assert_eq!(w.realign_cb(region, 1), Some(*truth));
    }

    #[test]
    fn realign_deletion() {
        let truth = b"AAACCCAAGAAACACT";
        let w = wl(&[truth]);
        // region = CB missing its 6th base (the 2nd 'C' group) -> 15 bases + UMI
        // truth AAACC C AAGAAACACT -> delete one C -> AAACCAAGAAACACT (15)
        let region = b"AAACCAAGAAACACTNNNNNNNNNNNN";
        assert_eq!(w.realign_cb(region, -1), Some(*truth));
    }

    #[test]
    fn within2_unique_and_ambiguous() {
        let truth = b"AAACCCAAGAAACACT";
        let w = wl(&[truth, b"TTTTTTTTTTTTTTTT"]);
        // exact
        assert_eq!(w.correct_within2(truth), Some(*truth));
        // two substitutions from truth -> unique
        assert_eq!(w.correct_within2(b"AAACCCAAGAAACAGG"), Some(*truth));
        // far from everything
        assert_eq!(w.correct_within2(b"GGGGGGGGGGGGGGGG"), None);
        // ambiguous: equidistant (1 mm) from two whitelist members
        let a = b"AAAAAAAAAAAAAAAA";
        let b = b"AAAAAAAAAAAAAAAC";
        let w2 = wl(&[a, b]);
        // "AAAAAAAAAAAAAAAG" is 1mm from both a and b -> ambiguous
        assert_eq!(w2.correct_within2(b"AAAAAAAAAAAAAAAG"), None);
    }

    #[test]
    fn within1_unique_and_ambiguous() {
        let truth = b"AAACCCAAGAAACACT";
        let w = wl(&[truth, b"TTTTTTTTTTTTTTTT"]);
        assert_eq!(w.correct_within1(truth), Some(*truth)); // exact
        assert_eq!(w.correct_within1(b"AAACCCAAGAAACACA"), Some(*truth)); // 1 mm
        assert_eq!(w.correct_within1(b"AAACCCAAGAAACAGG"), None); // 2 mm -> none
                                                                  // ambiguous: 1 mm from two whitelist members
        let w2 = wl(&[b"AAAAAAAAAAAAAAAA", b"AAAAAAAAAAAAAAAC"]);
        assert_eq!(w2.correct_within1(b"AAAAAAAAAAAAAAAG"), None);
        assert_eq!(w2.classify(b"AAAAAAAAAAAAAAAG"), CbStatus::Fail);
    }

    #[test]
    fn reframe_recovers_shifted_cb() {
        let truth = b"AAACCCAAGAAACACT";
        let w = wl(&[truth, b"TTTTTTTTTTTTTTTT"]);
        // Rigid cb_start grabs one base too early (a deletion put the CB one base
        // 3'): prepend a junk base so truth sits at offset +1. Within-1 noise on
        // the CB ('T' at the last pos) must still resolve uniquely.
        let mut seq = vec![b'G']; // [0] junk, rigid frame would grab seq[0..16]
        seq.extend_from_slice(b"AAACCCAAGAAACACT"); // truth at [1..17]
        seq.extend_from_slice(b"NNNNNNNNNNNN");
        let (cb, off) = w.reframe_cb(&seq, 0, 2).expect("unique +1 hit");
        assert_eq!(&cb, truth);
        assert_eq!(off, 1);
        // no nearby unique hit -> None
        let junk: Vec<u8> = b"GGGGGGGGGGGGGGGGGGGG".to_vec();
        assert_eq!(w.reframe_cb(&junk, 1, 2), None);
    }

    #[test]
    fn realign_ambiguous_or_none() {
        let w = wl(&[b"AAACCCAAGAAACACT", b"TTTTTTTTTTTTTTTT"]);
        // garbage region -> no whitelist hit
        assert_eq!(w.realign_cb(b"GGGGGGGGGGGGGGGGGGGGGGGGGGGG", 1), None);
        // delta 0 not handled
        assert_eq!(w.realign_cb(b"AAACCCAAGAAACACTNNNNNNNNNNNN", 0), None);
    }
}
