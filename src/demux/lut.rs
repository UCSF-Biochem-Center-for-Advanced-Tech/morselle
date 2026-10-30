//! All-variants i5/i7 lookup tables, precomputed once per run (XOOS-style).
//!
//! A sample sheet has tens to a few hundred indexes, so every read-side
//! sequence that should match one -- each index's ≤`index_max_mm` Hamming
//! neighbours, and (with `--index-indels`) its one-deletion and one-insertion
//! variants -- is enumerated up front. Keys are 2-bit packed and, for indexes
//! up to 10 bp, looked up in a directly addressed array (4^k slots, ≤ 8 MB for
//! the 11 bp insertion table) instead of hashing a byte slice per frame. The
//! nominal frame is tried first (most reads stop there); otherwise the anchor
//! region is packed once and every other frame is a shift, mask and array
//! read. Longer keys fall back to an `FxHashMap` on the packed code.
//! CE15204 bench vs the byte-slice hash maps: exact dual 36 → ~20 ns, no
//! match with `--index-indels` 395 → ~80 ns per read.
//!
//! Semantics match `HammingIndex` exactly: a variant claimed by two different
//! samples becomes ambiguous and never matches, and a frame holding any base
//! other than `A`/`C`/`G`/`T` matches nothing.

use rustc_hash::FxHashMap;

use super::indel::IndexMatch;
use crate::hamming_index::offset_order;

const EMPTY: u16 = 0;
const AMBIGUOUS: u16 = u16::MAX;
/// Largest sample count a slot can hold (`sample + 1`, below `AMBIGUOUS`).
pub const MAX_SAMPLES: usize = u16::MAX as usize - 1;
/// Keys up to this length use a dense array (4^11 slots × 2 B = 8 MB).
const DENSE_MAX_LEN: usize = 11;

/// 2-bit base codes; 4 = not ACGT.
const fn base_codes() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
}
static CODE: [u8; 256] = base_codes();

/// Pack `seq` (≤ 16 bp) 2 bits per base; `None` if any base is not ACGT.
#[inline]
fn pack(seq: &[u8]) -> Option<u32> {
    let mut code = 0u32;
    for &b in seq {
        let c = CODE[b as usize];
        if c > 3 {
            return None;
        }
        code = (code << 2) | c as u32;
    }
    Some(code)
}

enum Slots {
    Dense(Vec<u16>),
    Sparse(FxHashMap<u32, u16>),
}

/// Fixed-length packed-key table: sample per key, or ambiguous.
struct Table {
    len: usize,
    slots: Slots,
}

impl Table {
    fn new(len: usize) -> Self {
        assert!((1..=16).contains(&len), "index length {len} not in 1..=16");
        let slots = if len <= DENSE_MAX_LEN {
            Slots::Dense(vec![EMPTY; 1 << (2 * len)])
        } else {
            Slots::Sparse(FxHashMap::default())
        };
        Table { len, slots }
    }

    /// Claim `key` for `sample`; a key already held by another sample
    /// becomes ambiguous (and stays so).
    fn claim(&mut self, key: &[u8], sample: usize) {
        let Some(code) = pack(key) else { return };
        let v = sample as u16 + 1;
        let slot = match &mut self.slots {
            Slots::Dense(d) => &mut d[code as usize],
            Slots::Sparse(m) => m.entry(code).or_insert(EMPTY),
        };
        *slot = match *slot {
            EMPTY => v,
            s if s == v => v,
            _ => AMBIGUOUS,
        };
    }

    #[inline]
    fn get(&self, key: &[u8]) -> Option<usize> {
        self.get_code(pack(key)?)
    }

    #[inline]
    fn get_code(&self, code: u32) -> Option<usize> {
        let slot = match &self.slots {
            Slots::Dense(d) => d[code as usize],
            Slots::Sparse(m) => *m.get(&code)?,
        };
        (slot != EMPTY && slot != AMBIGUOUS).then(|| slot as usize - 1)
    }

    /// Look up the `len`-bp frame at `pos` in `seq`, if it fits.
    #[inline]
    fn at(&self, seq: &[u8], pos: isize) -> Option<usize> {
        if pos < 0 {
            return None;
        }
        self.get(seq.get(pos as usize..pos as usize + self.len)?)
    }
}

/// The bases around an index anchor, packed once so every frame's key is a
/// shift and mask instead of a fresh pack (≤ 32 bp: `2 * window + len + 1`).
struct Region {
    lo: isize,
    hi: isize,
    codes: u64,
    /// Bit `hi-1-i` set when base `i` is not ACGT.
    bad: u32,
}

impl Region {
    fn new(seq: &[u8], lo: isize, hi: isize) -> Option<Self> {
        // An anchor can sit past the read end (i7 after a late adapter):
        // clamp so the region is empty rather than out of bounds.
        let hi = hi.min(seq.len() as isize).max(0);
        let lo = lo.max(0).min(hi);
        if hi - lo > 32 {
            return None;
        }
        let (mut codes, mut bad) = (0u64, 0u32);
        for &b in &seq[lo as usize..hi as usize] {
            let c = CODE[b as usize];
            codes = (codes << 2) | u64::from(c & 3);
            bad = (bad << 1) | u32::from(c > 3);
        }
        Some(Region { lo, hi, codes, bad })
    }

    /// Packed key of the `len`-bp frame at `pos`, if inside and all ACGT.
    #[inline]
    fn key(&self, pos: isize, len: usize) -> Option<u32> {
        let end = pos + len as isize;
        if pos < self.lo || end > self.hi {
            return None;
        }
        let shift = (self.hi - end) as u32;
        if (self.bad >> shift) & ((1u32 << len) - 1) != 0 {
            return None;
        }
        Some(((self.codes >> (2 * shift)) & ((1u64 << (2 * len)) - 1)) as u32)
    }
}

/// Every matchable variant of one index type (i5 or i7) for a sample sheet.
pub struct IndexLut {
    len: usize,
    /// ≤`index_max_mm` Hamming neighbours, `len` bp.
    ham: Table,
    /// One-deletion (`len-1` bp) and one-insertion (`len+1` bp) variants,
    /// present only while `--index-indels` is active for this index.
    indel: Option<(Table, Table)>,
}

impl IndexLut {
    /// Hamming table for read-oriented `keys` (sample `i` = `keys[i]`).
    pub fn new(keys: &[Vec<u8>], max_mm: u32) -> Self {
        let len = keys.first().map_or(10, Vec::len);
        let mut ham = Table::new(len);
        for (sample, key) in keys.iter().enumerate() {
            let mut buf = key.clone();
            neighbours(&mut ham, &mut buf, 0, max_mm, sample);
        }
        IndexLut {
            len,
            ham,
            indel: None,
        }
    }

    /// Add the one-indel tables for the same `keys`.
    pub fn add_indels(&mut self, keys: &[Vec<u8>]) {
        let mut del = Table::new(self.len - 1);
        let mut ins = Table::new(self.len + 1);
        for (sample, key) in keys.iter().enumerate() {
            for i in 0..key.len() {
                let mut v = key.clone();
                v.remove(i);
                del.claim(&v, sample);
            }
            for i in 0..=key.len() {
                for &b in b"ACGT" {
                    let mut v = key.clone();
                    v.insert(i, b);
                    ins.claim(&v, sample);
                }
            }
        }
        self.indel = Some((del, ins));
    }

    pub fn drop_indels(&mut self) {
        self.indel = None;
    }

    pub fn has_indels(&self) -> bool {
        self.indel.is_some()
    }

    /// Exact-frame Hamming lookup (`key` already the index region).
    pub fn get(&self, key: &[u8]) -> Option<usize> {
        self.ham.get(key)
    }

    /// The full per-index search: mismatch-only over frames (0, ±1, …,
    /// ±`window`), first unique hit wins (unchanged behaviour); only if
    /// that finds nothing, the one-indel pass over the same frames, which
    /// must agree on a single sample across every frame and both lengths.
    /// `use_indels = false` skips the indel pass even when it is built.
    pub fn find(
        &self,
        seq: &[u8],
        start: Option<usize>,
        window: isize,
        use_indels: bool,
    ) -> IndexMatch {
        let Some(start) = start else {
            return IndexMatch::default();
        };
        let start = start as isize;
        // Most reads hit at the nominal frame: one pack, one read, done.
        if let Some(s) = self.ham.at(seq, start) {
            return IndexMatch {
                sample: Some(s),
                via_indel: false,
            };
        }
        let w = window.max(0);
        let Some(region) = Region::new(seq, start - w, start + w + self.len as isize + 1) else {
            return self.find_slow(seq, start, window, use_indels);
        };
        for k in offset_order(window) {
            let hit = region
                .key(start + k, self.len)
                .and_then(|c| self.ham.get_code(c));
            if let Some(s) = hit {
                return IndexMatch {
                    sample: Some(s),
                    via_indel: false,
                };
            }
        }
        let Some((del, ins)) = self.indel.as_ref().filter(|_| use_indels) else {
            return IndexMatch::default();
        };
        let mut found = None;
        for k in offset_order(window) {
            for table in [del, ins] {
                let hit = region
                    .key(start + k, table.len)
                    .and_then(|c| table.get_code(c));
                if let Some(v) = hit {
                    match found {
                        None => found = Some(v),
                        Some(f) if f == v => {}
                        Some(_) => return IndexMatch::default(),
                    }
                }
            }
        }
        IndexMatch {
            sample: found,
            via_indel: found.is_some(),
        }
    }

    /// [`find`](Self::find) packing each frame separately, for windows too
    /// wide to pack the whole region into one word.
    fn find_slow(&self, seq: &[u8], start: isize, window: isize, use_indels: bool) -> IndexMatch {
        for k in offset_order(window) {
            if let Some(s) = self.ham.at(seq, start + k) {
                return IndexMatch {
                    sample: Some(s),
                    via_indel: false,
                };
            }
        }
        let Some((del, ins)) = self.indel.as_ref().filter(|_| use_indels) else {
            return IndexMatch::default();
        };
        let mut found = None;
        for k in offset_order(window) {
            for table in [del, ins] {
                if let Some(v) = table.at(seq, start + k) {
                    match found {
                        None => found = Some(v),
                        Some(f) if f == v => {}
                        Some(_) => return IndexMatch::default(),
                    }
                }
            }
        }
        IndexMatch {
            sample: found,
            via_indel: found.is_some(),
        }
    }
}

/// Claim every ≤`budget`-substitution neighbour of `buf` (from `pos` on),
/// in the same order as `HammingIndex::insert`.
fn neighbours(table: &mut Table, buf: &mut [u8], pos: usize, budget: u32, sample: usize) {
    table.claim(buf, sample);
    if budget == 0 {
        return;
    }
    for i in pos..buf.len() {
        let orig = buf[i];
        for &b in b"ACGT" {
            if b != orig {
                buf[i] = b;
                neighbours(table, buf, i + 1, budget - 1, sample);
            }
        }
        buf[i] = orig;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hamming_index::HammingIndex;

    fn keys(n: usize, len: usize, seed: u32) -> Vec<Vec<u8>> {
        let mut s = seed;
        (0..n)
            .map(|_| {
                (0..len)
                    .map(|_| {
                        s = s.wrapping_mul(1_103_515_245).wrapping_add(12345);
                        b"ACGT"[((s >> 16) & 3) as usize]
                    })
                    .collect()
            })
            .collect()
    }

    /// Same answers as `HammingIndex::find_windowed` for every frame of a
    /// mixed haystack, dense and sparse, including ambiguous neighbours.
    #[test]
    fn matches_hamming_index_exactly() {
        for (n, len, mm) in [(14, 10, 1), (384, 10, 1), (96, 8, 2), (40, 12, 1)] {
            let ks = keys(n, len, 7 + n as u32);
            let mut old = HammingIndex::new();
            for (i, k) in ks.iter().enumerate() {
                old.insert(k, i, mm);
            }
            let lut = IndexLut::new(&ks, mm);
            let mut hay = Vec::new();
            for k in &ks {
                hay.extend_from_slice(b"GATTNACA");
                let mut v = k.clone();
                v[len / 2] = b'A';
                hay.extend_from_slice(&v);
            }
            for start in 0..hay.len() {
                for window in [0isize, 3, 12] {
                    assert_eq!(
                        lut.find(&hay, Some(start), window, true).sample,
                        old.find_windowed(&hay, Some(start), len, window),
                        "n={n} len={len} start={start} window={window}"
                    );
                }
            }
        }
    }

    #[test]
    fn non_acgt_frames_never_match() {
        let ks = vec![b"ACGTACGTAC".to_vec()];
        let lut = IndexLut::new(&ks, 1);
        assert_eq!(lut.get(b"ACGTACGTAC"), Some(0));
        assert_eq!(lut.get(b"ACGTNCGTAC"), None);
        assert_eq!(lut.get(b"acgtacgtac"), None);
    }
}

#[cfg(test)]
mod region_tests {
    use super::*;

    /// The packed-region path and the per-frame path agree on every start,
    /// including starts near both read ends, N bases, and the indel pass.
    #[test]
    fn region_path_matches_per_frame_path() {
        let ks: Vec<Vec<u8>> = vec![
            b"ACGTTGCATG".to_vec(),
            b"TTGACCAGGA".to_vec(),
            b"CAGTCAGTCC".to_vec(),
        ];
        let mut lut = IndexLut::new(&ks, 1);
        lut.add_indels(&ks);
        let mut hay = Vec::new();
        for k in &ks {
            hay.extend_from_slice(b"GANTACA");
            let mut d = k.clone();
            d.remove(3);
            hay.extend_from_slice(&d);
            hay.extend_from_slice(b"CC");
            hay.extend_from_slice(k);
        }
        for start in 0..hay.len() + 8 {
            for window in [0isize, 1, 3] {
                for use_indels in [false, true] {
                    assert_eq!(
                        lut.find(&hay, Some(start), window, use_indels),
                        lut.find_slow(&hay, start as isize, window, use_indels),
                        "start={start} window={window} indels={use_indels}"
                    );
                }
            }
        }
    }
}
