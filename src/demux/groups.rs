//! Interned candidate-sample sets for the index lookup tables.
//!
//! A table slot used to hold one sample or "ambiguous" (matches nothing).
//! For the per-read rule -- one index decides when it names exactly one
//! sample, otherwise both must agree -- a slot holds a *group*: the sorted
//! set of samples whose index (or one of its tolerated variants) is that
//! key. Group `i < n_samples` is the singleton `{i}`, so a slot value
//! `sample + 1` still means "exactly this sample"; shared keys (combinatorial
//! sheets reusing an i7 across samples, or two indexes whose neighbourhoods
//! overlap) get an interned multi-sample group after them. One interner is
//! shared by an index's Hamming and indel tables, so equal ids mean equal
//! sets across both.

use rustc_hash::FxHashMap;

/// Slot value for "no sample".
pub(super) const EMPTY: u16 = 0;
/// Slot value when the group store is full (matches nothing, as before).
pub(super) const OVERFLOW: u16 = u16::MAX;
/// Largest sample count a slot can hold (`sample + 1`, below `OVERFLOW`).
pub const MAX_SAMPLES: usize = u16::MAX as usize - 1;

pub(super) struct Groups {
    n_samples: usize,
    /// Members of group `slot - 1`: `{i}` for `i < n_samples`, then the
    /// interned multi-sample sets.
    sets: Vec<Vec<u16>>,
    ids: FxHashMap<Vec<u16>, u16>,
}

impl Groups {
    pub(super) fn new(n_samples: usize) -> Self {
        assert!(n_samples <= MAX_SAMPLES);
        Groups {
            n_samples,
            sets: (0..n_samples as u16).map(|i| vec![i]).collect(),
            ids: FxHashMap::default(),
        }
    }

    /// The slot value after `sample` also claims a slot holding `slot`.
    pub(super) fn claim(&mut self, slot: u16, sample: usize) -> u16 {
        if slot == EMPTY {
            return sample as u16 + 1;
        }
        if slot == OVERFLOW {
            return OVERFLOW;
        }
        let mut set = self.members(slot).to_vec();
        match set.binary_search(&(sample as u16)) {
            Ok(_) => return slot,
            Err(i) => set.insert(i, sample as u16),
        }
        if let Some(&id) = self.ids.get(&set) {
            return id;
        }
        let id = self.sets.len() + 1;
        if id >= OVERFLOW as usize {
            return OVERFLOW;
        }
        self.sets.push(set.clone());
        self.ids.insert(set, id as u16);
        id as u16
    }

    /// The single sample a slot names, if it names exactly one.
    #[inline]
    pub(super) fn single(&self, slot: u16) -> Option<usize> {
        (slot != EMPTY && (slot as usize) <= self.n_samples).then(|| slot as usize - 1)
    }

    /// Whether a slot names two or more samples.
    #[inline]
    pub(super) fn is_shared(&self, slot: u16) -> bool {
        slot != OVERFLOW && (slot as usize) > self.n_samples
    }

    /// Sorted members of a non-empty, non-overflow slot.
    pub(super) fn members(&self, slot: u16) -> &[u16] {
        &self.sets[slot as usize - 1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn claims_build_interned_sorted_sets() {
        let mut g = Groups::new(4);
        let a = g.claim(EMPTY, 2);
        assert_eq!(g.single(a), Some(2));
        assert_eq!(g.claim(a, 2), a);
        let ab = g.claim(a, 0);
        assert!(g.is_shared(ab));
        assert_eq!(g.single(ab), None);
        assert_eq!(g.members(ab), &[0, 2]);
        // Same set reached in another order: same id.
        let b = g.claim(EMPTY, 0);
        assert_eq!(g.claim(b, 2), ab);
        let abc = g.claim(ab, 3);
        assert_eq!(g.members(abc), &[0, 2, 3]);
        assert_eq!(g.claim(abc, 0), abc);
    }
}
