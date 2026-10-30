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
//! on `u32` keys. We do NOT rewrite barcodes — downstream compatibility tools correct them.

use std::io::BufRead;
use std::path::Path;

use anyhow::Result;
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
        let reader = crate::io::open_text_reader(path)?;

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

    /// Cost-1 insertion rescue (`--cb-insertion`, `XB:i:8`): the read gained
    /// one base inside the CB. `region` starts at the CB lock (≥ 17 bp).
    /// Returns the corrected CB and the dropped position (1..=15) iff exactly
    /// one whitelist barcode drops one inner base of the first 17, and no
    /// one-deletion reading (a base restored at 1..=15 of the first 15) names a
    /// different whitelist barcode -- the cost-1 uniqueness that was measured.
    ///
    /// Insertions only: on the CAT02 slices (private recovery-pools note,
    /// 2026-10-29) unique cost-1 insertion hits among CB-fails land on
    /// cell-like barcodes at 50.0% / 20.8% (barcode-first / barcode-last;
    /// Hamming-1 47.4% / 24.3%, reframe 38.1% / 19.5%), deletion hits at only
    /// 25.8% / 7.0%: a deletion reading has 4x the candidates (60 vs 15), and
    /// SBX makes about twice as many insertions as deletions.
    pub fn insertion_cb(&self, region: &[u8]) -> Option<([u8; 16], usize)> {
        if region.len() < 17 {
            return None;
        }
        let mut found: Option<(u32, usize)> = None;
        let mut cb = [0u8; 16];
        for skip in 1..16 {
            cb[..skip].copy_from_slice(&region[..skip]);
            cb[skip..].copy_from_slice(&region[skip + 1..17]);
            let Some(code) = encode16(&cb) else { continue };
            if self.contains_code(code) {
                match found {
                    None => found = Some((code, skip)),
                    Some((c, _)) if c != code => return None,
                    _ => {}
                }
            }
        }
        let (code, skip) = found?;
        for pos in 1..16 {
            cb[..pos].copy_from_slice(&region[..pos]);
            cb[pos + 1..].copy_from_slice(&region[pos..15]);
            for &b in b"ACGT" {
                cb[pos] = b;
                if encode16(&cb).is_some_and(|c| c != code && self.contains_code(c)) {
                    return None;
                }
            }
        }
        Some((decode16(code), skip))
    }
}

#[cfg(test)]
#[path = "barcode_tests.rs"]
mod tests;
