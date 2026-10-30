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
//! one-insertion variant of each index is precomputed (`super::lut`), so a
//! lookup stays one table read per frame: deletion variants are `L-1` bp
//! keys, insertion variants `L+1` bp keys. Exactly one indel and no further
//! mismatch -- one edit in total, the same budget as the default
//! `--index-max-mm 1`.
//!
//! Safety is checked per index before it is enabled ([`indel_safety`]): the
//! minimum pairwise edit distance must leave the one-edit neighbourhoods
//! disjoint (≥ 3), and every one-edit variant of every index, placed in its
//! real flanking adapter sequence, is run through the full search (mismatch
//! pass, then indel pass, over the whole window): the indel pass must add no
//! wrong-sample call beyond those mismatch-only matching already makes.

use super::lut::IndexLut;

/// How one index was matched.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct IndexMatch {
    pub sample: Option<usize>,
    /// Matched only through the one-indel pass.
    pub via_indel: bool,
    /// The matched key is shared by several samples (`sample` is `None`):
    /// its candidate-group id in that index's table (`IndexLut::members`),
    /// for the other index to choose among.
    pub shared: Option<u16>,
    /// No match because a shifted-frame mismatch hit had a competing sample
    /// (`IndexLut::search`).
    pub float_tie: bool,
}

impl IndexMatch {
    /// Matched one sample or a shared candidate set.
    pub fn matched(&self) -> bool {
        self.sample.is_some() || self.shared.is_some()
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
    /// One-edit variants (in flank context) that resolved to another sample
    /// (or to a shared candidate set without their own).
    pub wrong: usize,
    /// The same, with the indel pass off (mismatch-only matching): a shifted
    /// frame can turn one sample's indel into another sample's one-mismatch
    /// hit when two indexes are near-shifts of each other (CE15204: V4_GEX vs
    /// M2_GEX i5, 1 of 1,176 variants at window ≥1, before the float learned
    /// to abstain on such ties -- now 0).
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
/// through `lut` (which must already hold the indel tables) with and without
/// the indel pass. `keys[i]` must be sample `i`'s read-oriented index.
pub fn indel_safety(
    keys: &[Vec<u8>],
    lut: &IndexLut,
    left: &[u8],
    right: &[u8],
    window: isize,
) -> IndelSafety {
    // Samples sharing one index (combinatorial sheets) are told apart by
    // the other index; only distinct indexes must be ≥ 3 edits apart.
    let mut distinct = keys.to_vec();
    distinct.sort();
    distinct.dedup();
    let mut out = IndelSafety {
        min_edit: min_pairwise_edit(&distinct),
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
            let m = lut.find(&seq, Some(left.len()), window, true);
            out.variants += 1;
            let base = lut.find(&seq, Some(left.len()), window, false);
            if base.sample.is_some_and(|s| s != sample) {
                out.wrong_mismatch_only += 1;
            }
            match m.sample {
                Some(s) if s == sample => out.recovered += 1,
                Some(_) => out.wrong += 1,
                // A shared candidate set must at least contain the sample.
                None if m
                    .shared
                    .is_some_and(|g| !lut.members(g).contains(&(sample as u16))) =>
                {
                    out.wrong += 1
                }
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
