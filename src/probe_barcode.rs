//! 10x Flex probe-barcode location + identification (FLEX_SUPPORT_PLAN.md
//! Phase 1).
//!
//! Empirically confirmed with the real `orient()`/`extract()` code path (see
//! FLEX_SUPPORT_PLAN.md, "Phase 0 — results"): the 10bp probe barcode sits
//! ~10bp after `Extracted::polyt_end` in canonical (P5→P7) coordinates, as
//! `revcomp(whitelist_seq)`. That offset has real per-read jitter — `polyt_end`
//! is `dt_run_end`'s *measured* T-run boundary, which isn't always exactly at
//! the physical junction — so location is a windowed search centered on the
//! confirmed offset, not a fixed exact position. Built on the same shared
//! [`crate::hamming_index::HammingIndex`] primitive `demux.rs` uses for i5/i7.

use std::fs;
use std::path::Path;

use anyhow::{ensure, Context, Result};

use crate::anchors::revcomp;
use crate::hamming_index::{offset_order, HammingIndex};

/// Probe-barcode length (10x Flex).
pub const PROBE_BARCODE_LEN: usize = 10;

/// Empirically confirmed center of the probe-barcode window, in bp past
/// `Extracted::polyt_end` (canonical coordinates). See FLEX_SUPPORT_PLAN.md
/// Phase 0 results: this offset hit ~62% of confidently-oriented reads,
/// dwarfing the next-highest offset — the windowed search below is what
/// absorbs the remaining per-read jitter around it.
pub const DEFAULT_PROBE_OFFSET: isize = 10;

/// A located probe barcode: the raw observed sequence (in whitelist-forward
/// orientation, i.e. already revcomp'd out of canonical coordinates) and the
/// corrected probe-barcode id, when a confident whitelist match was found.
pub struct ProbeHit {
    /// Raw 10bp window as observed, revcomp'd into whitelist orientation —
    /// directly comparable to `id`'s whitelist sequence, no mental revcomp
    /// needed when eyeballing a BAM's `PR`/`PB` tags side by side.
    pub raw: Vec<u8>,
    /// Corrected probe-barcode id (e.g. `"D-A01"`), if `raw` matched a
    /// whitelist entry within the configured max mismatches.
    pub id: Option<String>,
}

/// 10x Flex probe-barcode whitelist (e.g. `flex-v2-384.txt`): Hamming-tolerant
/// lookup from observed sequence to probe-barcode id.
pub struct ProbeWhitelist {
    index: HammingIndex<u16>,
    ids: Vec<String>,
    window: isize,
}

impl ProbeWhitelist {
    /// Load a whitelist file (`seq\tseq\tid` per line, as 10x ships them —
    /// only the 1st and 3rd columns are used). Builds the index over
    /// `revcomp(seq)`, the orientation the barcode appears in on a canonical
    /// (P5→P7) read.
    pub fn load(path: &Path, max_mm: u32, window: isize) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("reading probe-barcode whitelist {}", path.display()))?;
        let mut index = HammingIndex::new();
        let mut ids = Vec::new();
        for (lineno, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut cols = line.split_whitespace();
            let seq = cols.next().unwrap_or_default();
            let id = cols.nth(1).unwrap_or(seq); // 3rd column when present, else the seq itself
            ensure!(
                seq.len() == PROBE_BARCODE_LEN,
                "{}:{}: expected {PROBE_BARCODE_LEN}bp probe barcode, got {}bp ({seq})",
                path.display(),
                lineno + 1,
                seq.len()
            );
            let idx: u16 = ids.len().try_into().with_context(|| {
                format!(
                    "probe-barcode whitelist {} has >65535 entries",
                    path.display()
                )
            })?;
            index.insert(&revcomp(seq.as_bytes()), idx, max_mm);
            ids.push(id.to_string());
        }
        ensure!(
            !ids.is_empty(),
            "empty probe-barcode whitelist: {}",
            path.display()
        );
        Ok(Self { index, ids, window })
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    /// Locate + identify the probe barcode on a canonical read. Searches a
    /// window of `self.window` bp either side of
    /// `polyt_end + `[`DEFAULT_PROBE_OFFSET`], smallest offset first, and
    /// returns the first confident (≤`self.max_mm`) whitelist hit. When no
    /// offset in the window produces a confident hit, still reports the raw
    /// sequence observed at the center offset (id = `None`), so a miss is
    /// still visible in `PR` rather than silently absent — mirrors how a
    /// CB whitelist-fail still carries its raw `CR`.
    pub fn locate(&self, canon_seq: &[u8], polyt_end: usize) -> Option<ProbeHit> {
        let center = polyt_end as isize + DEFAULT_PROBE_OFFSET;
        for k in offset_order(self.window) {
            let start = center + k;
            if start < 0 {
                continue;
            }
            let start = start as usize;
            let end = start + PROBE_BARCODE_LEN;
            if end > canon_seq.len() {
                continue;
            }
            // Index keys are `revcomp(whitelist_seq)` — i.e. already in the
            // orientation the barcode physically appears in a canonical read
            // — so the lookup uses the canonical window as-is; `raw` (for
            // reporting) is revcomp'd back into whitelist-forward orientation.
            let window = &canon_seq[start..end];
            if let Some(idx) = self.index.get(window) {
                return Some(ProbeHit {
                    raw: revcomp(window),
                    id: Some(self.ids[idx as usize].clone()),
                });
            }
        }
        let start = center.max(0) as usize;
        let end = start + PROBE_BARCODE_LEN;
        (end <= canon_seq.len()).then(|| ProbeHit {
            raw: revcomp(&canon_seq[start..end]),
            id: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn whitelist_file(entries: &[(&str, &str)]) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().unwrap();
        for (seq, id) in entries {
            writeln!(f, "{seq}\t{seq}\t{id}").unwrap();
        }
        f
    }

    /// Build a canonical read whose probe barcode sits at
    /// `polyt_end + DEFAULT_PROBE_OFFSET`, embedded as `revcomp(barcode)`.
    fn read_with_probe(polyt_end: usize, barcode: &[u8]) -> Vec<u8> {
        let mut seq = vec![b'N'; polyt_end];
        let start = polyt_end + DEFAULT_PROBE_OFFSET as usize;
        seq.resize(start, b'N');
        seq.extend_from_slice(&revcomp(barcode));
        seq.extend_from_slice(b"NNNNNNNNNN"); // trailing padding
        seq
    }

    #[test]
    fn loads_and_reports_len() {
        let f = whitelist_file(&[("ACGTACGGGA", "A-A01"), ("GAATGGGAAG", "A-A02")]);
        let wl = ProbeWhitelist::load(f.path(), 1, 5).unwrap();
        assert_eq!(wl.len(), 2);
    }

    #[test]
    fn rejects_wrong_length_barcode() {
        let f = whitelist_file(&[("ACGTACGGG", "A-A01")]); // 9bp, not 10
        assert!(ProbeWhitelist::load(f.path(), 1, 5).is_err());
    }

    #[test]
    fn locates_exact_match_at_confirmed_offset() {
        let f = whitelist_file(&[("ACGTACGGGA", "A-A01"), ("GAATGGGAAG", "A-A02")]);
        let wl = ProbeWhitelist::load(f.path(), 1, 5).unwrap();
        let seq = read_with_probe(20, b"ACGTACGGGA");
        let hit = wl.locate(&seq, 20).unwrap();
        assert_eq!(hit.id.as_deref(), Some("A-A01"));
        assert_eq!(hit.raw, b"ACGTACGGGA");
    }

    #[test]
    fn windowed_search_absorbs_offset_jitter() {
        let f = whitelist_file(&[("ACGTACGGGA", "A-A01")]);
        let wl = ProbeWhitelist::load(f.path(), 1, 5).unwrap();
        // polyt_end measured 2bp early (jitter): true site is polyt_end+2+DEFAULT_PROBE_OFFSET.
        let seq = read_with_probe(22, b"ACGTACGGGA");
        let hit = wl.locate(&seq, 20).unwrap();
        assert_eq!(hit.id.as_deref(), Some("A-A01"));
    }

    #[test]
    fn no_confident_match_still_reports_raw() {
        let f = whitelist_file(&[("ACGTACGGGA", "A-A01")]);
        let wl = ProbeWhitelist::load(f.path(), 1, 5).unwrap();
        let seq = read_with_probe(20, b"TTTTTTTTTT"); // not in whitelist, far from A-A01
        let hit = wl.locate(&seq, 20).unwrap();
        assert_eq!(hit.id, None);
        assert_eq!(hit.raw.len(), PROBE_BARCODE_LEN);
    }

    #[test]
    fn out_of_bounds_returns_none() {
        let f = whitelist_file(&[("ACGTACGGGA", "A-A01")]);
        let wl = ProbeWhitelist::load(f.path(), 1, 5).unwrap();
        let short_seq = vec![b'N'; 15];
        assert!(wl.locate(&short_seq, 20).is_none());
    }
}
