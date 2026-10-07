use crate::hamming_index::HammingIndex;
use super::types::{Probe, DEFAULT_PROBE_SEQ_OFFSET};

pub fn single_deletion_variants(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    (0..seq.len()).map(move |i| {
        let mut v = Vec::with_capacity(seq.len() - 1);
        v.extend_from_slice(&seq[..i]);
        v.extend_from_slice(&seq[i + 1..]);
        v
    })
}

pub fn single_insertion_variants(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    (0..=seq.len()).flat_map(move |i| {
        let existing = if i < seq.len() { seq[i] } else { b'\0' };
        b"ACGT"
            .iter()
            .filter(move |&&b| b != existing)
            .map(move |&b| {
                let mut v = Vec::with_capacity(seq.len() + 1);
                v.extend_from_slice(&seq[..i]);
                v.push(b);
                v.extend_from_slice(&seq[i..]);
                v
            })
    })
}

pub struct ProbeMatcher {
    index: HammingIndex<usize>,
    del_index: HammingIndex<usize>,
    ins_index: HammingIndex<usize>,
    probe_len: usize,
    window: isize,
}

impl ProbeMatcher {
    pub fn build(probes: &[Probe], max_mm: u32, window: isize) -> Self {
        let probe_len = probes.first().map_or(0, |p| p.seq.len());
        let mut index = HammingIndex::new();
        let mut del_index = HammingIndex::new();
        let mut ins_index = HammingIndex::new();
        for (i, p) in probes.iter().enumerate() {
            index.insert(&p.seq, i, max_mm);
            for v in single_deletion_variants(&p.seq) {
                del_index.insert(&v, i, 0);
            }
            for v in single_insertion_variants(&p.seq) {
                ins_index.insert(&v, i, 0);
            }
        }
        Self {
            index,
            del_index,
            ins_index,
            probe_len,
            window,
        }
    }

    pub fn locate(&self, read: &[u8]) -> Option<usize> {
        for k in crate::hamming_index::offset_order(self.window) {
            let start = DEFAULT_PROBE_SEQ_OFFSET as isize + k;
            if start < 0 {
                continue;
            }
            let start = start as usize;
            if let Some(end) = start.checked_add(self.probe_len) {
                if end <= read.len() {
                    if let Some(probe_idx) = self.index.get(&read[start..end]) {
                        return Some(probe_idx);
                    }
                }
            }
            if let Some(end) = start.checked_add(self.probe_len.saturating_sub(1)) {
                if end <= read.len() {
                    if let Some(probe_idx) = self.del_index.get(&read[start..end]) {
                        return Some(probe_idx);
                    }
                }
            }
            let end = start + self.probe_len + 1;
            if end <= read.len() {
                if let Some(probe_idx) = self.ins_index.get(&read[start..end]) {
                    return Some(probe_idx);
                }
            }
        }
        None
    }
}
