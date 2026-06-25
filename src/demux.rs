//! i5/i7 extraction + sample-sheet matching (plan §4 step 3, M3).
//!
//! Empirically (validated against the 14 known sample indices), in a canonical
//! P5→P7 read:
//! - **i5** sits 10 bp immediately 5' of the R1 primer, as `revcomp(sheet i5)`.
//! - **i7** sits 10 bp after `rc(R2 primer)` (i.e. after the 3' adapter), as the
//!   sheet i7 forward.
//!
//! Both positions are derived from coordinates M1/M2 already found, so demux adds
//! no new search. Indices are matched via a precomputed ≤`index_max_mm` neighbor
//! map (O(1) lookup, with ambiguity detection). The sample index set here has min
//! pairwise Hamming i7=3 / i5=4, so single-index assignment (when the other index
//! is truncated — common for i7 on long reads) is unambiguous.

use std::collections::HashMap;

use crate::anchors::{revcomp, R1_PRIMER, R2_PRIMER};
use crate::sample_sheet::SampleSheet;

/// Tuning for demultiplexing.
#[derive(Debug, Clone, Copy)]
pub struct DemuxCfg {
    /// Max mismatches when matching an index.
    pub index_max_mm: u32,
    /// Allow assignment from a single index when the other is unavailable.
    pub allow_single_index: bool,
}

/// How a read was assigned (for accounting).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DemuxVia {
    /// Both indices present and agree.
    Dual,
    /// Only i5 available (i7 truncated/absent).
    I5Only,
    /// Only i7 available (i5 truncated/absent).
    I7Only,
    /// Both present but disagree.
    Conflict,
    /// No usable index.
    None,
}

/// Demux outcome: assigned sample index (into [`Demuxer::sample_ids`]) + how.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Outcome {
    pub sample: Option<usize>,
    pub via: DemuxVia,
}

#[derive(Debug, Clone, Copy)]
enum Slot {
    Unique(usize),
    Ambiguous,
}

/// Matches extracted index sequences to samples.
pub struct Demuxer {
    sample_ids: Vec<String>,
    index_len: usize,
    i5_map: HashMap<Box<[u8]>, Slot>,
    i7_map: HashMap<Box<[u8]>, Slot>,
    cfg: DemuxCfg,
}

impl Demuxer {
    /// Build from a sample sheet. The read-oriented index sequences are
    /// `revcomp(i5)` and `i7` (forward), each expanded to ≤`index_max_mm`
    /// neighbors for O(1), error-tolerant lookup.
    pub fn new(sheet: &SampleSheet, cfg: DemuxCfg) -> Self {
        let mut i5_map = HashMap::new();
        let mut i7_map = HashMap::new();
        let index_len = sheet.samples.first().map_or(10, |s| s.i7.len());
        for (idx, s) in sheet.samples.iter().enumerate() {
            insert_neighbors(&mut i5_map, &revcomp(&s.i5), idx, cfg.index_max_mm);
            insert_neighbors(&mut i7_map, &s.i7, idx, cfg.index_max_mm);
        }
        Demuxer {
            sample_ids: sheet.samples.iter().map(|s| s.id.clone()).collect(),
            index_len,
            i5_map,
            i7_map,
            cfg,
        }
    }

    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    pub fn index_len(&self) -> usize {
        self.index_len
    }

    fn lookup(map: &HashMap<Box<[u8]>, Slot>, region: &[u8]) -> Option<usize> {
        match map.get(region) {
            Some(Slot::Unique(i)) => Some(*i),
            _ => None,
        }
    }

    /// Combine per-index sample matches into an outcome.
    fn combine(&self, a5: Option<usize>, a7: Option<usize>) -> Outcome {
        match (a5, a7) {
            (Some(x), Some(y)) if x == y => Outcome { sample: Some(x), via: DemuxVia::Dual },
            (Some(_), Some(_)) => Outcome { sample: None, via: DemuxVia::Conflict },
            (Some(x), None) if self.cfg.allow_single_index => {
                Outcome { sample: Some(x), via: DemuxVia::I5Only }
            }
            (None, Some(y)) if self.cfg.allow_single_index => {
                Outcome { sample: Some(y), via: DemuxVia::I7Only }
            }
            _ => Outcome { sample: None, via: DemuxVia::None },
        }
    }

    /// Assign a read from its (optional) i5/i7 read-region slices (fixed frame).
    pub fn assign(&self, i5_region: Option<&[u8]>, i7_region: Option<&[u8]>) -> Outcome {
        let a5 = i5_region.and_then(|r| Self::lookup(&self.i5_map, r));
        let a7 = i7_region.and_then(|r| Self::lookup(&self.i7_map, r));
        self.combine(a5, a7)
    }

    /// Best unique sample match for an index anchored at `start`, trying offsets
    /// (0, ±1, …, ±`window`) to absorb R1-primer / adapter indels that shift the
    /// rigid frame. Returns the first offset (smallest |k|) giving a unique hit.
    fn match_windowed(
        map: &HashMap<Box<[u8]>, Slot>,
        seq: &[u8],
        start: Option<usize>,
        index_len: usize,
        window: isize,
    ) -> Option<usize> {
        let start = start? as isize;
        for k in offset_order(window) {
            let pos = start + k;
            if pos < 0 {
                continue;
            }
            let pos = pos as usize;
            if pos + index_len > seq.len() {
                continue;
            }
            if let Some(s) = Self::lookup(map, &seq[pos..pos + index_len]) {
                return Some(s);
            }
        }
        None
    }

    /// Assign a read by offset-tolerant matching around the i5/i7 anchor starts.
    /// `window` = max frame shift to try each side (0 = rigid, like [`assign`]).
    pub fn assign_windowed(
        &self,
        seq: &[u8],
        i5_start: Option<usize>,
        i7_start: Option<usize>,
        window: isize,
    ) -> Outcome {
        let a5 = Self::match_windowed(&self.i5_map, seq, i5_start, self.index_len, window);
        let a7 = Self::match_windowed(&self.i7_map, seq, i7_start, self.index_len, window);
        self.combine(a5, a7)
    }
}

/// Offsets to try, smallest magnitude first: 0, 1, -1, 2, -2, …
fn offset_order(window: isize) -> impl Iterator<Item = isize> {
    std::iter::once(0).chain((1..=window).flat_map(|k| [k, -k]))
}

/// Anchor START positions for i5 and i7 (before any offset search).
/// i5 starts `R1_PRIMER.len() + index_len` before the CB; i7 starts
/// `R2_PRIMER.len()` after the 3' adapter.
pub fn index_anchor_starts(
    cb_start: usize,
    adapter_pos: Option<usize>,
    index_len: usize,
) -> (Option<usize>, Option<usize>) {
    let i5 = cb_start.checked_sub(R1_PRIMER.len() + index_len);
    let i7 = adapter_pos.map(|p| p + R2_PRIMER.len());
    (i5, i7)
}

/// Slice the i5 and i7 read-regions from a canonical read using M1/M2 coords.
///
/// i5 = the `index_len` bp ending where the R1 primer starts (`cb_start -
/// R1_PRIMER.len()`). i7 = the `index_len` bp starting after `rc(R2 primer)`
/// (`adapter_pos + R2_PRIMER.len()`). Either is `None` if the read doesn't span it.
pub fn index_regions<'a>(
    seq: &'a [u8],
    cb_start: usize,
    adapter_pos: Option<usize>,
    index_len: usize,
) -> (Option<&'a [u8]>, Option<&'a [u8]>) {
    let i5 = cb_start
        .checked_sub(R1_PRIMER.len() + index_len)
        .map(|s| &seq[s..s + index_len]);
    let i7 = adapter_pos.and_then(|p| {
        let s = p + R2_PRIMER.len();
        (s + index_len <= seq.len()).then(|| &seq[s..s + index_len])
    });
    (i5, i7)
}

/// Insert all ≤`max_mm` Hamming neighbors of `seq` into `map`, marking entries
/// claimed by more than one sample as ambiguous.
fn insert_neighbors(map: &mut HashMap<Box<[u8]>, Slot>, seq: &[u8], sample: usize, max_mm: u32) {
    let mut buf = seq.to_vec();
    enumerate(map, &mut buf, 0, max_mm, sample);
}

fn enumerate(
    map: &mut HashMap<Box<[u8]>, Slot>,
    buf: &mut [u8],
    pos: usize,
    budget: u32,
    sample: usize,
) {
    // Record the current variant.
    map.entry(buf.to_vec().into_boxed_slice())
        .and_modify(|slot| {
            if !matches!(slot, Slot::Unique(s) if *s == sample) {
                *slot = Slot::Ambiguous;
            }
        })
        .or_insert(Slot::Unique(sample));
    if budget == 0 {
        return;
    }
    for i in pos..buf.len() {
        let orig = buf[i];
        for &b in b"ACGT" {
            if b != orig {
                buf[i] = b;
                enumerate(map, buf, i + 1, budget - 1, sample);
            }
        }
        buf[i] = orig;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn sheet() -> SampleSheet {
        let p = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/CE15204_bclconvert_input.csv");
        SampleSheet::load(&p).unwrap()
    }

    fn cfg() -> DemuxCfg {
        DemuxCfg { index_max_mm: 1, allow_single_index: true }
    }

    #[test]
    fn dual_index_exact() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg());
        // M2_GEX: i7=TGTAGTCATT, i5=TACGATCAAG. Read carries rc(i5) and i7-fwd.
        let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
        let i5_read = revcomp(b"TACGATCAAG");
        let i7_read = b"TGTAGTCATT";
        let o = d.assign(Some(&i5_read), Some(i7_read));
        assert_eq!(o.sample, Some(m2));
        assert_eq!(o.via, DemuxVia::Dual);
        assert_eq!(d.sample_ids()[m2], "M2_GEX");
    }

    #[test]
    fn one_mismatch_tolerated() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg());
        let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
        let mut i5_read = revcomp(b"TACGATCAAG");
        i5_read[3] = if i5_read[3] == b'A' { b'C' } else { b'A' }; // 1 mm
        let o = d.assign(Some(&i5_read), Some(b"TGTAGTCATT"));
        assert_eq!(o.sample, Some(m2));
    }

    #[test]
    fn single_index_when_i7_truncated() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg());
        let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
        let i5_read = revcomp(b"TACGATCAAG");
        let o = d.assign(Some(&i5_read), None);
        assert_eq!(o.sample, Some(m2));
        assert_eq!(o.via, DemuxVia::I5Only);
    }

    #[test]
    fn conflict_goes_undetermined() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg());
        // i5 of M2_GEX but i7 of M7_GEX -> conflict
        let o = d.assign(Some(&revcomp(b"TACGATCAAG")), Some(b"TAGCATAGTG"));
        assert_eq!(o.sample, None);
        assert_eq!(o.via, DemuxVia::Conflict);
    }

    #[test]
    fn unknown_index_undetermined() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg());
        let o = d.assign(Some(b"GGGGGGGGGG"), Some(b"GGGGGGGGGG"));
        assert_eq!(o.sample, None);
        assert_eq!(o.via, DemuxVia::None);
    }

    #[test]
    fn single_index_disabled() {
        let sh = sheet();
        let mut c = cfg();
        c.allow_single_index = false;
        let d = Demuxer::new(&sh, c);
        let o = d.assign(Some(&revcomp(b"TACGATCAAG")), None);
        assert_eq!(o.sample, None);
    }

    #[test]
    fn windowed_recovers_shifted_index() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg());
        let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
        // M2 read-oriented i5 = revcomp(TACGATCAAG); place it at position 30.
        let i5_read = revcomp(b"TACGATCAAG");
        let mut seq = vec![b'G'; 80];
        seq[30..40].copy_from_slice(&i5_read);
        // Anchor is off by one (true index starts at 30, anchor says 29).
        let rigid = d.assign_windowed(&seq, Some(29), None, 0);
        assert_eq!(rigid.sample, None); // rigid frame misses it
        let windowed = d.assign_windowed(&seq, Some(29), None, 3);
        assert_eq!(windowed.sample, Some(m2)); // offset search recovers it
        assert_eq!(windowed.via, DemuxVia::I5Only);
    }

    #[test]
    fn index_regions_slicing() {
        // R1_PRIMER (33) + i5 (10) before cb_start; rc(R2_PRIMER) (34) before i7.
        let cb_start = 43;
        let mut seq = vec![b'N'; 200];
        seq[0..10].copy_from_slice(b"AAAACCCCGG"); // i5 region [0,10)
        let adapter_pos = 100;
        seq[134..144].copy_from_slice(b"TTTTGGGGCC"); // i7 region [134,144)
        let (i5, i7) = index_regions(&seq, cb_start, Some(adapter_pos), 10);
        assert_eq!(i5, Some(&b"AAAACCCCGG"[..]));
        assert_eq!(i7, Some(&b"TTTTGGGGCC"[..]));
    }
}
