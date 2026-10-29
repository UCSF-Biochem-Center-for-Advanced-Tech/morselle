//! Single-indel i5/i7 matching (`--index-indels`).
//!
//! The mismatch-only matcher ([`HammingIndex::find_windowed`]) floats the
//! index a few bp to absorb an indel *upstream* of it, but an insertion or
//! deletion *inside* the 10 bp index leaves no frame with ≤1 mismatch. SBX
//! reads carry ~3 insertions and ~1.5 deletions per kb, and on CAT02 an
//! edit-1 search recovers ~35% of undetermined reads and turns ~50% of
//! i7-only (barcode-first) assignments into dual ones, with 0 wrong samples
//! on 2,999 dual-assigned reads (private recovery-pools measurement,
//! 2026-10-29).
//!
//! Like XOOS's edited-sequence lookup tables, every one-deletion and
//! one-insertion variant of each index is precomputed, so a lookup stays a
//! hash probe per frame: deletion variants are `L-1` bp keys, insertion
//! variants `L+1` bp keys. Exactly one indel and no further mismatch -- one
//! edit in total, the same budget as the default `--index-max-mm 1`.
//!
//! Safety is checked per index before it is enabled ([`indel_safety`]): the
//! minimum pairwise edit distance must leave the one-edit neighbourhoods
//! disjoint (≥ 3), and every one-edit variant of every index, placed in its
//! real flanking adapter sequence, is run through the full search (mismatch
//! pass, then indel pass, over the whole window): the indel pass must add no
//! wrong-sample call beyond those mismatch-only matching already makes.

use crate::hamming_index::{offset_order, HammingIndex};

/// Precomputed one-insertion / one-deletion variants of a set of keys.
pub struct IndelIndex {
    key_len: usize,
    /// Read lost one index base: `key_len - 1` bp keys.
    del: HammingIndex<usize>,
    /// Read gained one base: `key_len + 1` bp keys.
    ins: HammingIndex<usize>,
}

impl IndelIndex {
    pub fn new(key_len: usize) -> Self {
        IndelIndex {
            key_len,
            del: HammingIndex::new(),
            ins: HammingIndex::new(),
        }
    }

    /// Insert every one-indel variant of `key` (read orientation) for `value`.
    /// A variant shared with a different value becomes ambiguous and is refused.
    pub fn insert(&mut self, key: &[u8], value: usize) {
        for i in 0..key.len() {
            let mut v = key.to_vec();
            v.remove(i);
            self.del.insert(&v, value, 0);
        }
        for i in 0..=key.len() {
            for &b in b"ACGT" {
                let mut v = key.to_vec();
                v.insert(i, b);
                self.ins.insert(&v, value, 0);
            }
        }
    }

    /// Unique sample over every frame (0, ±1, …, ±`window`) and both variant
    /// lengths; `None` when nothing matches or two samples do.
    pub fn find_windowed(&self, seq: &[u8], start: usize, window: isize) -> Option<usize> {
        let mut found = None;
        for k in offset_order(window) {
            let pos = start as isize + k;
            if pos < 0 {
                continue;
            }
            let pos = pos as usize;
            for (table, len) in [(&self.del, self.key_len - 1), (&self.ins, self.key_len + 1)] {
                let Some(key) = seq.get(pos..pos + len) else {
                    continue;
                };
                if let Some(v) = table.get(key) {
                    match found {
                        None => found = Some(v),
                        Some(f) if f == v => {}
                        Some(_) => return None,
                    }
                }
            }
        }
        found
    }
}

/// How one index was matched.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct IndexMatch {
    pub sample: Option<usize>,
    /// Matched only through the one-indel pass.
    pub via_indel: bool,
}

/// The full per-index search: mismatch-only first (unchanged behaviour), then
/// the indel pass only if that found nothing anywhere in the window.
pub fn find_index(
    ham: &HammingIndex<usize>,
    indel: Option<&IndelIndex>,
    seq: &[u8],
    start: Option<usize>,
    key_len: usize,
    window: isize,
) -> IndexMatch {
    if let Some(s) = ham.find_windowed(seq, start, key_len, window) {
        return IndexMatch {
            sample: Some(s),
            via_indel: false,
        };
    }
    match (indel, start) {
        (Some(ix), Some(st)) => {
            let sample = ix.find_windowed(seq, st, window);
            IndexMatch {
                sample,
                via_indel: sample.is_some(),
            }
        }
        _ => IndexMatch::default(),
    }
}

/// Plain Levenshtein distance (short keys only).
pub fn edit_distance(a: &[u8], b: &[u8]) -> usize {
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, &x) in a.iter().enumerate() {
        let mut cur = vec![i + 1; b.len() + 1];
        for (j, &y) in b.iter().enumerate() {
            cur[j + 1] = (prev[j] + usize::from(x != y))
                .min(prev[j + 1] + 1)
                .min(cur[j] + 1);
        }
        prev = cur;
    }
    prev[b.len()]
}

/// Minimum pairwise edit distance among `keys` (`usize::MAX` for < 2 keys).
pub fn min_pairwise_edit(keys: &[Vec<u8>]) -> usize {
    let mut best = usize::MAX;
    for i in 0..keys.len() {
        for j in (i + 1)..keys.len() {
            best = best.min(edit_distance(&keys[i], &keys[j]));
        }
    }
    best
}

/// Result of [`indel_safety`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IndelSafety {
    pub min_edit: usize,
    /// One-edit variants (in flank context) that resolved to another sample.
    pub wrong: usize,
    /// The same, with the indel pass off (mismatch-only matching as it was):
    /// a shifted frame can already turn one sample's indel into another
    /// sample's one-mismatch hit when two indexes are near-shifts of each
    /// other (CE15204: V4_GEX vs M2_GEX i5, 1 of 1,176 variants at window ≥1).
    pub wrong_mismatch_only: usize,
    /// One-edit variants that resolved to their own sample.
    pub recovered: usize,
    pub variants: usize,
}

impl IndelSafety {
    /// The indel pass is safe to enable when the indexes' one-edit
    /// neighbourhoods are disjoint and it adds no wrong-sample call beyond
    /// those mismatch-only matching already makes.
    pub fn safe(&self) -> bool {
        self.min_edit >= 3 && self.wrong <= self.wrong_mismatch_only
    }
}

/// Every one-edit (substitution, insertion, deletion) variant of each key,
/// placed between its real `left`/`right` flanks at the nominal start, run
/// through [`find_index`]. `keys[i]` must be sample `i`'s read-oriented index.
pub fn indel_safety(
    keys: &[Vec<u8>],
    ham: &HammingIndex<usize>,
    indel: &IndelIndex,
    left: &[u8],
    right: &[u8],
    window: isize,
) -> IndelSafety {
    let mut out = IndelSafety {
        min_edit: min_pairwise_edit(keys),
        wrong: 0,
        wrong_mismatch_only: 0,
        recovered: 0,
        variants: 0,
    };
    for (sample, key) in keys.iter().enumerate() {
        for v in one_edit_variants(key) {
            let mut seq = left.to_vec();
            seq.extend_from_slice(&v);
            seq.extend_from_slice(right);
            let m = find_index(ham, Some(indel), &seq, Some(left.len()), key.len(), window);
            out.variants += 1;
            let base = find_index(ham, None, &seq, Some(left.len()), key.len(), window);
            if base.sample.is_some_and(|s| s != sample) {
                out.wrong_mismatch_only += 1;
            }
            match m.sample {
                Some(s) if s == sample => out.recovered += 1,
                Some(_) => out.wrong += 1,
                None => {}
            }
        }
    }
    out
}

fn one_edit_variants(key: &[u8]) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    for i in 0..key.len() {
        for &b in b"ACGT" {
            if b != key[i] {
                let mut v = key.to_vec();
                v[i] = b;
                out.push(v);
            }
        }
        let mut v = key.to_vec();
        v.remove(i);
        out.push(v);
    }
    for i in 0..=key.len() {
        for &b in b"ACGT" {
            let mut v = key.to_vec();
            v.insert(i, b);
            out.push(v);
        }
    }
    out
}
