use super::types::{Probe, DEFAULT_PROBE_SEQ_OFFSET};
use crate::hamming_index::HammingIndex;

/// All single-position deletions of `seq` (`seq.len() - 1` variants, one per
/// deleted position). Real CAT06 data (2026-09-18 diagnostic, real
/// edit-distance recovery check against 20K Hamming-miss reads): of misses
/// recoverable within edit distance <=2, deletions dominate 85.2% vs 9.2%
/// insertion-only -- a real base the read is missing relative to the probe,
/// not a substitution a mismatch budget could ever cover.
pub fn single_deletion_variants(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    (0..seq.len()).map(move |i| {
        let mut v = Vec::with_capacity(seq.len() - 1);
        v.extend_from_slice(&seq[..i]);
        v.extend_from_slice(&seq[i + 1..]);
        v
    })
}

/// All single-position insertions into `seq` (up to `(seq.len() + 1) * 3`
/// variants: every position, every non-redundant base). The rarer indel
/// shape in real data (9.2% of recoverable indel-shaped misses vs 85.2%
/// deletion) but real, and cheap enough to cover alongside deletion.
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

/// Matches a read's probe sequence against the reference via a windowed,
/// Hamming-tolerant lookup (see [`DEFAULT_PROBE_SEQ_OFFSET`]) -- structurally
/// identical to `probe_barcode.rs::ProbeWhitelist::locate`, just with a
/// much larger reference and no revcomp step (probe_seq is confirmed to
/// appear forward in the emitted read, unlike the probe barcode).
///
/// Also tries single-indel-shifted windows (`del_index`/`ins_index`) at
/// each offset -- real SBX errors are dominated by single-base indels, not
/// substitutions, inside the probe span itself (not just jitter in where it
/// starts, which the offset window already absorbs): of real Hamming-miss
/// reads recoverable by edit distance <=2, 98.5% were genuinely
/// indel-shaped. Rather than a per-read full-database edit-distance scan
/// (intractable at billions-of-reads scale), this extends the same
/// precomputed-neighborhood trick [`HammingIndex`] already uses for
/// substitutions to single indels too -- still an O(1) hash lookup per
/// attempt, just at two more fixed window lengths (49bp/51bp alongside the
/// exact 50bp).
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

    /// Locate the probe sequence within a window centered on
    /// [`DEFAULT_PROBE_SEQ_OFFSET`], smallest offset first, returning the
    /// matched probe's index into the original `probes` slice. At each
    /// offset, tries the exact-length substitution-tolerant lookup first,
    /// then the single-deletion- and single-insertion-shifted windows.
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
            // Deletion-shaped: the read is missing a base the probe has, so
            // the real payload is one bp SHORTER than probe_len here.
            if let Some(end) = start.checked_add(self.probe_len.saturating_sub(1)) {
                if end <= read.len() {
                    if let Some(probe_idx) = self.del_index.get(&read[start..end]) {
                        return Some(probe_idx);
                    }
                }
            }
            // Insertion-shaped: the read carries an extra base the probe
            // doesn't, so the real payload is one bp LONGER than probe_len.
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
