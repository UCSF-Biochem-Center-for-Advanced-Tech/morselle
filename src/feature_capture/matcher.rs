use std::collections::HashMap;
use crate::hamming_index::HammingIndex;
use super::types::Feature;

struct PatternGroup {
    capture_len: usize,
    capture_index: HammingIndex<usize>,
}

struct AnchorIndexSet {
    exact: HammingIndex<usize>,
    del: HammingIndex<usize>,
    ins: HammingIndex<usize>,
}

pub struct FeatureMatcher {
    anchor_by_len: Vec<(usize, AnchorIndexSet)>,
    groups: Vec<PatternGroup>,
    anchor_window: usize,
}

impl FeatureMatcher {
    pub fn build(
        features: &[Feature],
        anchor_max_mm: u32,
        capture_max_mm: u32,
        anchor_window: usize,
    ) -> Self {
        let mut group_of: HashMap<(&[u8], usize), usize> = HashMap::new();
        let mut groups: Vec<PatternGroup> = Vec::new();
        let mut anchor_by_len: HashMap<usize, AnchorIndexSet> = HashMap::new();

        for (fi, f) in features.iter().enumerate() {
            let key = (f.anchor.as_slice(), f.barcode.len());
            let group_id = *group_of.entry(key).or_insert_with(|| {
                groups.push(PatternGroup {
                    capture_len: f.barcode.len(),
                    capture_index: HammingIndex::new(),
                });
                let id = groups.len() - 1;
                let entry = anchor_by_len.entry(f.anchor.len()).or_insert_with(|| AnchorIndexSet {
                    exact: HammingIndex::new(),
                    del: HammingIndex::new(),
                    ins: HammingIndex::new(),
                });
                entry.exact.insert(&f.anchor, id, anchor_max_mm);
                for v in crate::probe_set::single_deletion_variants(&f.anchor) {
                    entry.del.insert(&v, id, 0);
                }
                for v in crate::probe_set::single_insertion_variants(&f.anchor) {
                    entry.ins.insert(&v, id, 0);
                }
                id
            });
            groups[group_id]
                .capture_index
                .insert(&f.barcode, fi, capture_max_mm);
        }

        Self {
            anchor_by_len: anchor_by_len.into_iter().collect(),
            groups,
            anchor_window,
        }
    }

    fn match_capture(&self, group_id: usize, read: &[u8], bc_start: usize) -> Option<usize> {
        let group = &self.groups[group_id];
        let bc_end = bc_start.checked_add(group.capture_len)?;
        if bc_end > read.len() {
            return None;
        }
        group.capture_index.get(&read[bc_start..bc_end])
    }

    pub fn locate(&self, read: &[u8]) -> Option<usize> {
        for (anchor_len, idx_set) in &self.anchor_by_len {
            let anchor_len = *anchor_len;
            let hay_end = read.len().min(self.anchor_window + anchor_len);
            if hay_end < anchor_len {
                continue;
            }
            for offset in 0..=(hay_end - anchor_len) {
                if let Some(group_id) = idx_set.exact.get(&read[offset..offset + anchor_len]) {
                    if let Some(feature_idx) = self.match_capture(group_id, read, offset + anchor_len) {
                        return Some(feature_idx);
                    }
                }
                if anchor_len >= 1 {
                    let del_end = offset + anchor_len - 1;
                    if del_end <= read.len() {
                        if let Some(group_id) = idx_set.del.get(&read[offset..del_end]) {
                            if let Some(feature_idx) = self.match_capture(group_id, read, del_end) {
                                return Some(feature_idx);
                            }
                        }
                    }
                }
                let ins_end = offset + anchor_len + 1;
                if ins_end <= read.len() {
                    if let Some(group_id) = idx_set.ins.get(&read[offset..ins_end]) {
                        if let Some(feature_idx) = self.match_capture(group_id, read, ins_end) {
                            return Some(feature_idx);
                        }
                    }
                }
            }
        }
        None
    }
}
