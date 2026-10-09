//! Evidence merge routines shared by indexed region workers.
use super::types::{Evidence, RawKey};
use rustc_hash::FxHashMap as HashMap;

pub fn merge_item(target: &mut Evidence, source: Evidence) {
    target.reads += source.reads;
    target.min_read_hash = target.min_read_hash.min(source.min_read_hash);
    if source.best_mapq > target.best_mapq
        || (source.best_mapq == target.best_mapq
            && (source.best_aligned > target.best_aligned
                || (source.best_aligned == target.best_aligned
                    && source.best_ordinal < target.best_ordinal)))
    {
        target.best_ordinal = source.best_ordinal;
        target.best_mapq = source.best_mapq;
        target.best_aligned = source.best_aligned;
    }
}

pub fn merge_evidence_maps(
    mut left: HashMap<RawKey, Evidence>,
    mut right: HashMap<RawKey, Evidence>,
) -> HashMap<RawKey, Evidence> {
    // Keeping the larger allocation avoids repeated hash-table growth during
    // the parallel reduction tree.
    if left.len() < right.len() {
        std::mem::swap(&mut left, &mut right);
    }
    for (key, item) in right {
        left.entry(key)
            .and_modify(|current| merge_item(current, item))
            .or_insert(item);
    }
    left
}
