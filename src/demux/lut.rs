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
//! Semantics match `HammingIndex` for the single-sample answer: a frame
//! holding any base other than `A`/`C`/`G`/`T` matches nothing, and a key
//! claimed by several samples never names one of them on its own. Such a key
//! now keeps its candidate set (`super::groups`) so the other index can pick
//! among them (`Demuxer`'s per-read rule).

use rustc_hash::FxHashMap;

pub use super::groups::MAX_SAMPLES;
use super::groups::{Groups, EMPTY, OVERFLOW};
use super::indel::IndexMatch;
use crate::hamming_index::offset_order;

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

/// Fixed-length packed-key table: candidate group per key ([`Groups`]).
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

    /// Add `sample` to `key`'s candidate group.
    fn claim(&mut self, groups: &mut Groups, key: &[u8], sample: usize) {
        let Some(code) = pack(key) else { return };
        let slot = match &mut self.slots {
            Slots::Dense(d) => &mut d[code as usize],
            Slots::Sparse(m) => m.entry(code).or_insert(EMPTY),
        };
        *slot = groups.claim(*slot, sample);
    }

    /// Raw slot of a packed key (`EMPTY` when unclaimed).
    #[inline]
    fn slot(&self, code: u32) -> u16 {
        match &self.slots {
            Slots::Dense(d) => d[code as usize],
            Slots::Sparse(m) => m.get(&code).copied().unwrap_or(EMPTY),
        }
    }

    /// Slot of the `len`-bp frame at `pos` in `seq` (`EMPTY` if it does not
    /// fit or holds a non-ACGT base): the nominal-frame fast path.
    #[inline]
    fn slot_at(&self, seq: &[u8], pos: isize) -> u16 {
        if pos < 0 {
            return EMPTY;
        }
        seq.get(pos as usize..pos as usize + self.len)
            .and_then(pack)
            .map_or(EMPTY, |c| self.slot(c))
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
    /// One-deletion (`len-1` bp) and one-insertion (`len+1` bp) variants.
    /// Always built: they assign reads only while `--index-indels` is active
    /// for this index (`indel_on`), but every shifted-frame mismatch hit is
    /// checked against them (see [`search`](Self::search)).
    indel: Option<(Table, Table)>,
    indel_on: bool,
    /// Packed read-oriented keys, for mismatch counts (`None`: not ACGT).
    keys: Vec<Option<u32>>,
    groups: Groups,
}

impl IndexLut {
    /// Tables for read-oriented `keys` (sample `i` = `keys[i]`).
    pub fn new(keys: &[Vec<u8>], max_mm: u32) -> Self {
        let len = keys.first().map_or(10, Vec::len);
        let mut ham = Table::new(len);
        let mut groups = Groups::new(keys.len());
        for (sample, key) in keys.iter().enumerate() {
            let mut buf = key.clone();
            neighbours(&mut ham, &mut groups, &mut buf, 0, max_mm, sample);
        }
        let indel = (2..=15)
            .contains(&len)
            .then(|| indel_tables(keys, len, &mut groups));
        IndexLut {
            len,
            ham,
            indel,
            indel_on: false,
            keys: keys.iter().map(|k| pack(k)).collect(),
            groups,
        }
    }

    /// Let the one-indel tables assign reads (`--index-indels`).
    pub fn add_indels(&mut self) {
        self.indel_on = self.indel.is_some();
    }

    pub fn drop_indels(&mut self) {
        self.indel_on = false;
    }

    pub fn has_indels(&self) -> bool {
        self.indel_on
    }

    /// Exact-frame lookup (`key` already the index region).
    pub fn get(&self, key: &[u8]) -> IndexMatch {
        let slot = pack(key).map_or(EMPTY, |c| self.ham.slot(c));
        self.to_match(slot, false)
    }

    /// Sorted candidate samples of a shared match ([`IndexMatch::shared`]).
    pub fn members(&self, group: u16) -> &[u16] {
        self.groups.members(group)
    }

    fn to_match(&self, slot: u16, via_indel: bool) -> IndexMatch {
        let sample = self.groups.single(slot);
        let shared = self.groups.is_shared(slot).then_some(slot);
        IndexMatch {
            sample,
            via_indel: via_indel && (sample.is_some() || shared.is_some()),
            shared,
            float_tie: false,
        }
    }

    /// The full per-index search: mismatch-only over frames (0, ±1, …,
    /// ±`window`), first single-sample hit wins -- unless it is a mismatch
    /// hit in a shifted frame and another sample competes (below); with
    /// none, the first frame whose key is shared by several samples gives
    /// that candidate set. Only if neither, the one-indel pass over the same
    /// frames, which must agree on one candidate set (usually one sample)
    /// across every frame and both lengths. `use_indels = false` skips the
    /// indel pass even when it is on.
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
        let nominal = self.ham.slot_at(seq, start);
        if self.groups.single(nominal).is_some() {
            return self.to_match(nominal, false);
        }
        let w = window.max(0);
        match Region::new(seq, start - w, start + w + self.len as isize + 1) {
            Some(region) => self.search(start, window, use_indels, |pos, len| region.key(pos, len)),
            None => self.find_slow(seq, start, window, use_indels),
        }
    }

    /// [`find`](Self::find) packing each frame separately, for windows too
    /// wide to pack the whole region into one word.
    fn find_slow(&self, seq: &[u8], start: isize, window: isize, use_indels: bool) -> IndexMatch {
        self.search(start, window, use_indels, |pos, len| {
            let pos = usize::try_from(pos).ok()?;
            pack(seq.get(pos..pos + len)?)
        })
    }

    /// The frame scan behind [`find`](Self::find); `key(pos, len)` packs one
    /// frame (`None` if it does not fit or is not all ACGT).
    ///
    /// A one-mismatch hit in a shifted frame is two edits from the nominal
    /// index (the shift and the mismatch), so it is weaker evidence than a
    /// competing sample: when any other frame names a different sample, or
    /// the read is one insertion or deletion from another sample's index,
    /// the index abstains (`float_tie`). CE15204: V4_GEX's read-oriented i5
    /// with one inserted base is M2_GEX's i5 one mismatch off at frame +1.
    fn search(
        &self,
        start: isize,
        window: isize,
        use_indels: bool,
        key: impl Fn(isize, usize) -> Option<u32>,
    ) -> IndexMatch {
        let slot = |pos: isize, t: &Table| key(pos, t.len).map_or(EMPTY, |c| t.slot(c));
        let mut shared = None;
        for k in offset_order(window) {
            let s = slot(start + k, &self.ham);
            if let Some(sample) = self.groups.single(s) {
                let exact = key(start + k, self.len)
                    .zip(self.keys[sample])
                    .is_some_and(|(c, kc)| c == kc);
                if k != 0 && !exact && self.competes(start, window, sample, &slot) {
                    return IndexMatch {
                        float_tie: true,
                        ..IndexMatch::default()
                    };
                }
                return self.to_match(s, false);
            }
            if shared.is_none() && self.groups.is_shared(s) {
                shared = Some(s);
            }
        }
        if let Some(s) = shared {
            return self.to_match(s, false);
        }
        let Some((del, ins)) = self.indel.as_ref().filter(|_| use_indels && self.indel_on) else {
            return IndexMatch::default();
        };
        let mut found = None;
        for k in offset_order(window) {
            for table in [del, ins] {
                let v = slot(start + k, table);
                if v == EMPTY || v == OVERFLOW {
                    continue;
                }
                match found {
                    None => found = Some(v),
                    Some(f) if f == v => {}
                    Some(_) => return IndexMatch::default(),
                }
            }
        }
        found.map_or_else(IndexMatch::default, |v| self.to_match(v, true))
    }

    /// Whether any frame names candidates without `sample`: another Hamming
    /// hit, or a one-insertion/deletion variant of another sample's index.
    fn competes(
        &self,
        start: isize,
        window: isize,
        sample: usize,
        slot: &impl Fn(isize, &Table) -> u16,
    ) -> bool {
        let other = |s: u16| {
            s != EMPTY && s != OVERFLOW && !self.groups.members(s).contains(&(sample as u16))
        };
        let tables = std::iter::once(&self.ham).chain(self.indel.iter().flat_map(|(d, i)| [d, i]));
        tables
            .flat_map(|t| offset_order(window).map(move |k| (t, k)))
            .any(|(t, k)| other(slot(start + k, t)))
    }
}

/// One-deletion (`len-1` bp) and one-insertion (`len+1` bp) variant tables.
fn indel_tables(keys: &[Vec<u8>], len: usize, groups: &mut Groups) -> (Table, Table) {
    let mut del = Table::new(len - 1);
    let mut ins = Table::new(len + 1);
    for (sample, key) in keys.iter().enumerate() {
        for i in 0..key.len() {
            let mut v = key.clone();
            v.remove(i);
            del.claim(groups, &v, sample);
        }
        for i in 0..=key.len() {
            for &b in b"ACGT" {
                let mut v = key.clone();
                v.insert(i, b);
                ins.claim(groups, &v, sample);
            }
        }
    }
    (del, ins)
}

/// Claim every ≤`budget`-substitution neighbour of `buf` (from `pos` on),
/// in the same order as `HammingIndex::insert`.
fn neighbours(
    table: &mut Table,
    groups: &mut Groups,
    buf: &mut [u8],
    pos: usize,
    budget: u32,
    sample: usize,
) {
    table.claim(groups, buf, sample);
    if budget == 0 {
        return;
    }
    for i in pos..buf.len() {
        let orig = buf[i];
        for &b in b"ACGT" {
            if b != orig {
                buf[i] = b;
                neighbours(table, groups, buf, i + 1, budget - 1, sample);
            }
        }
        buf[i] = orig;
    }
}

#[cfg(test)]
#[path = "lut_region_tests.rs"]
mod region_tests;
#[cfg(test)]
#[path = "lut_tests.rs"]
mod tests;
