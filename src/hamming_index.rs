//! Generic fixed-length-key, Hamming-tolerant O(1) lookup index.
//!
//! Extracted from `demux.rs`'s original i5/i7 index (the Flex design notes
//! Phase 1) so the same primitive backs both sample-index demux and
//! probe-barcode/feature-barcode matching, rather than duplicating it.
//!
//! Every inserted key's ≤`max_mm` Hamming neighborhood is precomputed, so a
//! lookup is a single hash lookup regardless of `max_mm`. A neighbor claimed
//! by more than one distinct value becomes ambiguous and is refused rather
//! than guessed — this only matters for keys close enough together that their
//! error clouds overlap; see `demux.rs`'s doc comment on the real index sets'
//! min pairwise Hamming distances for why that's rare in practice.

use std::collections::HashMap;

#[derive(Debug, Clone, Copy)]
enum Slot<T> {
    Unique(T),
    Ambiguous,
}

/// Hamming-tolerant lookup from fixed-length byte keys to a small `Copy`
/// value (a sample index, a probe-barcode id, ...).
pub struct HammingIndex<T> {
    map: HashMap<Box<[u8]>, Slot<T>>,
}

impl<T: Copy + Eq> Default for HammingIndex<T> {
    fn default() -> Self {
        Self {
            map: HashMap::new(),
        }
    }
}

impl<T: Copy + Eq> HammingIndex<T> {
    pub fn new() -> Self {
        Self::default()
    }

    /// Insert all ≤`max_mm` Hamming neighbors of `key`, mapped to `value`. A
    /// neighbor already claimed by a *different* value becomes ambiguous (and
    /// stays so even if this same key is later re-inserted with that value).
    pub fn insert(&mut self, key: &[u8], value: T, max_mm: u32) {
        let mut buf = key.to_vec();
        self.enumerate(&mut buf, 0, max_mm, value);
    }

    fn enumerate(&mut self, buf: &mut [u8], pos: usize, budget: u32, value: T) {
        self.map
            .entry(buf.to_vec().into_boxed_slice())
            .and_modify(|slot| {
                if !matches!(slot, Slot::Unique(v) if *v == value) {
                    *slot = Slot::Ambiguous;
                }
            })
            .or_insert(Slot::Unique(value));
        if budget == 0 {
            return;
        }
        for i in pos..buf.len() {
            let orig = buf[i];
            for &b in b"ACGT" {
                if b != orig {
                    buf[i] = b;
                    self.enumerate(buf, i + 1, budget - 1, value);
                }
            }
            buf[i] = orig;
        }
    }

    /// Exact-frame lookup: `key` must already be the right length/position.
    pub fn get(&self, key: &[u8]) -> Option<T> {
        match self.map.get(key) {
            Some(Slot::Unique(v)) => Some(*v),
            _ => None,
        }
    }

    /// Best unique match for a `key_len`-bp key anchored at `start` in `seq`,
    /// trying offsets (0, ±1, …, ±`window`) smallest magnitude first. Absorbs
    /// a few bp of frame shift (e.g. an upstream indel) without a full
    /// alignment. Returns the first offset giving a unique hit.
    pub fn find_windowed(
        &self,
        seq: &[u8],
        start: Option<usize>,
        key_len: usize,
        window: isize,
    ) -> Option<T> {
        let start = start? as isize;
        for k in offset_order(window) {
            let pos = start + k;
            if pos < 0 {
                continue;
            }
            let pos = pos as usize;
            if pos + key_len > seq.len() {
                continue;
            }
            if let Some(v) = self.get(&seq[pos..pos + key_len]) {
                return Some(v);
            }
        }
        None
    }
}

/// Offsets to try, smallest magnitude first: 0, 1, -1, 2, -2, …
pub fn offset_order(window: isize) -> impl Iterator<Item = isize> {
    std::iter::once(0).chain((1..=window).flat_map(|k| [k, -k]))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_and_neighbor_lookup() {
        let mut idx = HammingIndex::new();
        idx.insert(b"ACGTACGTAC", 7usize, 1);
        assert_eq!(idx.get(b"ACGTACGTAC"), Some(7));
        let mut mm = b"ACGTACGTAC".to_vec();
        mm[2] = b'T';
        assert_eq!(idx.get(&mm), Some(7));
    }

    #[test]
    fn beyond_budget_misses() {
        let mut idx = HammingIndex::new();
        idx.insert(b"ACGTACGTAC", 7usize, 1);
        let mut mm2 = b"ACGTACGTAC".to_vec();
        mm2[2] = b'T';
        mm2[5] = b'A';
        assert_eq!(idx.get(&mm2), None);
    }

    #[test]
    fn overlapping_neighbors_become_ambiguous() {
        let mut idx = HammingIndex::new();
        idx.insert(b"ACGTACGTAC", 1usize, 1);
        idx.insert(b"ACGAACGTAC", 2usize, 1); // 1 mm from the first
        assert_eq!(idx.get(b"ACGTACGTAC"), None);
        assert_eq!(idx.get(b"ACGAACGTAC"), None);
    }

    #[test]
    fn windowed_recovers_shifted_key() {
        let mut idx = HammingIndex::new();
        idx.insert(b"ACGTACGTAC", 5usize, 0);
        let mut seq = vec![b'G'; 40];
        seq[10..20].copy_from_slice(b"ACGTACGTAC");
        assert_eq!(idx.find_windowed(&seq, Some(9), 10, 0), None);
        assert_eq!(idx.find_windowed(&seq, Some(9), 10, 3), Some(5));
    }

    #[test]
    fn windowed_none_when_start_missing() {
        let idx: HammingIndex<usize> = HammingIndex::new();
        assert_eq!(idx.find_windowed(b"ACGTACGTAC", None, 10, 3), None);
    }
}
