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

use anyhow::{bail, Result};

use crate::anchors::{revcomp, R1_PRIMER, R2_PRIMER};
use crate::hamming_index::HammingIndex;
use crate::sample_sheet::SampleSheet;

/// Which index to use for demux. `I5` needs only `r1s` (no dependency on the
/// 3' adapter ever being found), so restricting to it reaches the widest
/// possible population of reads -- but is only safe when the sample sheet's
/// i5 values alone separate every sample with enough Hamming margin.
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum IndexSource {
    I5,
    I7,
}

/// Tuning for demultiplexing.
#[derive(Debug, Clone, Copy)]
pub struct DemuxCfg {
    /// Max mismatches when matching an index.
    pub index_max_mm: u32,
    /// Allow assignment from a single index when the other is unavailable.
    pub allow_single_index: bool,
    /// Restrict demux to one index, skipping the other entirely. `None`
    /// (default) keeps today's dual-index-with-fallback behavior.
    pub single_index_source: Option<IndexSource>,
}

/// Real minimum pairwise Hamming distance among the sheet's samples, for the
/// given index alone. Determines whether single-index demux is safe:
/// unambiguous nearest-neighbor decoding at `max_mm` mismatches needs
/// `min_dist > 2 * max_mm`.
pub fn min_pairwise_hamming(sheet: &SampleSheet, source: IndexSource) -> u32 {
    let seqs: Vec<&[u8]> = sheet
        .samples
        .iter()
        .map(|s| match source {
            IndexSource::I5 => s.i5.as_slice(),
            IndexSource::I7 => s.i7.as_slice(),
        })
        .collect();
    let mut min_dist = u32::MAX;
    for i in 0..seqs.len() {
        for j in (i + 1)..seqs.len() {
            let d = seqs[i]
                .iter()
                .zip(seqs[j])
                .filter(|(a, b)| a != b)
                .count() as u32;
            min_dist = min_dist.min(d);
        }
    }
    if seqs.len() < 2 {
        0
    } else {
        min_dist
    }
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

/// Matches extracted index sequences to samples.
pub struct Demuxer {
    sample_ids: Vec<String>,
    index_len: usize,
    i5_index: HammingIndex<usize>,
    i7_index: HammingIndex<usize>,
    cfg: DemuxCfg,
}

impl Demuxer {
    /// Build from a sample sheet. The read-oriented index sequences are
    /// `revcomp(i5)` and `i7` (forward), each expanded to ≤`index_max_mm`
    /// neighbors for O(1), error-tolerant lookup.
    ///
    /// When `cfg.single_index_source` is set, refuses (rather than silently
    /// risking ambiguous assignment) unless that index alone separates every
    /// sample with enough real Hamming margin for `cfg.index_max_mm`.
    pub fn new(sheet: &SampleSheet, cfg: DemuxCfg) -> Result<Self> {
        if let Some(source) = cfg.single_index_source {
            let min_dist = min_pairwise_hamming(sheet, source);
            let required = 2 * cfg.index_max_mm + 1;
            if min_dist < required {
                bail!(
                    "--demux-index {source:?} is unsafe for this sample sheet: real min pairwise \
                     Hamming distance is {min_dist}, but unambiguous decoding at --index-max-mm \
                     {mm} needs > {mm2}. Use dual-index demux (omit --demux-index) or lower \
                     --index-max-mm.",
                    mm = cfg.index_max_mm,
                    mm2 = required - 1,
                );
            }
        }
        let mut i5_index = HammingIndex::new();
        let mut i7_index = HammingIndex::new();
        let index_len = sheet.samples.first().map_or(10, |s| s.i7.len());
        for (idx, s) in sheet.samples.iter().enumerate() {
            i5_index.insert(&revcomp(&s.i5), idx, cfg.index_max_mm);
            i7_index.insert(&s.i7, idx, cfg.index_max_mm);
        }
        Ok(Demuxer {
            sample_ids: sheet.samples.iter().map(|s| s.id.clone()).collect(),
            index_len,
            i5_index,
            i7_index,
            cfg,
        })
    }

    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    pub fn index_len(&self) -> usize {
        self.index_len
    }

    /// Combine per-index sample matches into an outcome.
    fn combine(&self, a5: Option<usize>, a7: Option<usize>) -> Outcome {
        match (a5, a7) {
            (Some(x), Some(y)) if x == y => Outcome {
                sample: Some(x),
                via: DemuxVia::Dual,
            },
            (Some(_), Some(_)) => Outcome {
                sample: None,
                via: DemuxVia::Conflict,
            },
            (Some(x), None) if self.cfg.allow_single_index => Outcome {
                sample: Some(x),
                via: DemuxVia::I5Only,
            },
            (None, Some(y)) if self.cfg.allow_single_index => Outcome {
                sample: Some(y),
                via: DemuxVia::I7Only,
            },
            _ => Outcome {
                sample: None,
                via: DemuxVia::None,
            },
        }
    }

    /// One index's outcome, standing alone (single-index mode): a hit is
    /// reported the same way dual-index reports a fallback single-index hit
    /// (`I5Only`/`I7Only`), since semantically it's the same situation --
    /// only one index was ever consulted.
    fn solo(sample: Option<usize>, via: DemuxVia) -> Outcome {
        match sample {
            Some(x) => Outcome { sample: Some(x), via },
            None => Outcome { sample: None, via: DemuxVia::None },
        }
    }

    /// Assign a read from its (optional) i5/i7 read-region slices (fixed frame).
    pub fn assign(&self, i5_region: Option<&[u8]>, i7_region: Option<&[u8]>) -> Outcome {
        match self.cfg.single_index_source {
            Some(IndexSource::I5) => {
                Self::solo(i5_region.and_then(|r| self.i5_index.get(r)), DemuxVia::I5Only)
            }
            Some(IndexSource::I7) => {
                Self::solo(i7_region.and_then(|r| self.i7_index.get(r)), DemuxVia::I7Only)
            }
            None => {
                let a5 = i5_region.and_then(|r| self.i5_index.get(r));
                let a7 = i7_region.and_then(|r| self.i7_index.get(r));
                self.combine(a5, a7)
            }
        }
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
        match self.cfg.single_index_source {
            Some(IndexSource::I5) => Self::solo(
                self.i5_index.find_windowed(seq, i5_start, self.index_len, window),
                DemuxVia::I5Only,
            ),
            Some(IndexSource::I7) => Self::solo(
                self.i7_index.find_windowed(seq, i7_start, self.index_len, window),
                DemuxVia::I7Only,
            ),
            None => {
                let a5 = self.i5_index.find_windowed(seq, i5_start, self.index_len, window);
                let a7 = self.i7_index.find_windowed(seq, i7_start, self.index_len, window);
                self.combine(a5, a7)
            }
        }
    }
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
        DemuxCfg {
            index_max_mm: 1,
            allow_single_index: true,
            single_index_source: None,
        }
    }

    #[test]
    fn dual_index_exact() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg()).unwrap();
        // M2_GEX: i7=TGTAGTCATT, i5=TACGATCAAG. Read carries rc(i5) and i7-fwd.
        let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
        let i5_read = revcomp(b"TACGATCAAG");
        let i7_read = b"TGTAGTCATT";
        let o = d.assign(Some(&i5_read), Some(i7_read));
        assert_eq!(o.sample, Some(m2));
        assert_eq!(o.via, DemuxVia::Dual);
        assert_eq!(d.sample_ids()[m2], "M2_GEX");
    }

    /// Pins the real fixture's own min pairwise Hamming distances, matching
    /// this module's own doc comment ("min pairwise Hamming i7=3 / i5=4") --
    /// if the fixture ever changes, this must be re-verified, not just
    /// silently pass with different numbers.
    #[test]
    fn real_fixture_min_pairwise_hamming_matches_doc_comment() {
        let sh = sheet();
        assert_eq!(min_pairwise_hamming(&sh, IndexSource::I5), 4);
        assert_eq!(min_pairwise_hamming(&sh, IndexSource::I7), 3);
    }

    #[test]
    fn single_index_i5_ignores_a_conflicting_i7() {
        let sh = sheet();
        let d = Demuxer::new(
            &sh,
            DemuxCfg {
                index_max_mm: 1,
                allow_single_index: true,
                single_index_source: Some(IndexSource::I5),
            },
        )
        .unwrap();
        let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
        let i5_read = revcomp(b"TACGATCAAG"); // M2_GEX's real i5
        let wrong_i7 = b"AAAAAAAAAA"; // not any real sample's i7 -- must be ignored entirely
        let o = d.assign(Some(&i5_read), Some(wrong_i7));
        assert_eq!(o.sample, Some(m2));
        assert_eq!(o.via, DemuxVia::I5Only);
    }

    #[test]
    fn single_index_refused_when_sheet_is_not_safe_at_this_mm() {
        let sh = sheet();
        // Real i5 min pairwise distance is 4 (see the pinned test above), so
        // unambiguous decoding needs index_max_mm <= 1 (> 2*mm required).
        // mm=2 needs min_dist > 4, which the real sheet does not have --
        // must be refused, not silently allowed.
        let result = Demuxer::new(
            &sh,
            DemuxCfg {
                index_max_mm: 2,
                allow_single_index: true,
                single_index_source: Some(IndexSource::I5),
            },
        );
        assert!(result.is_err());
    }

    #[test]
    fn one_mismatch_tolerated() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg()).unwrap();
        let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
        let mut i5_read = revcomp(b"TACGATCAAG");
        i5_read[3] = if i5_read[3] == b'A' { b'C' } else { b'A' }; // 1 mm
        let o = d.assign(Some(&i5_read), Some(b"TGTAGTCATT"));
        assert_eq!(o.sample, Some(m2));
    }

    #[test]
    fn single_index_when_i7_truncated() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg()).unwrap();
        let m2 = sh.samples.iter().position(|s| s.id == "M2_GEX").unwrap();
        let i5_read = revcomp(b"TACGATCAAG");
        let o = d.assign(Some(&i5_read), None);
        assert_eq!(o.sample, Some(m2));
        assert_eq!(o.via, DemuxVia::I5Only);
    }

    #[test]
    fn conflict_goes_undetermined() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg()).unwrap();
        // i5 of M2_GEX but i7 of M7_GEX -> conflict
        let o = d.assign(Some(&revcomp(b"TACGATCAAG")), Some(b"TAGCATAGTG"));
        assert_eq!(o.sample, None);
        assert_eq!(o.via, DemuxVia::Conflict);
    }

    #[test]
    fn unknown_index_undetermined() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg()).unwrap();
        let o = d.assign(Some(b"GGGGGGGGGG"), Some(b"GGGGGGGGGG"));
        assert_eq!(o.sample, None);
        assert_eq!(o.via, DemuxVia::None);
    }

    #[test]
    fn single_index_disabled() {
        let sh = sheet();
        let mut c = cfg();
        c.allow_single_index = false;
        let d = Demuxer::new(&sh, c).unwrap();
        let o = d.assign(Some(&revcomp(b"TACGATCAAG")), None);
        assert_eq!(o.sample, None);
    }

    #[test]
    fn windowed_recovers_shifted_index() {
        let sh = sheet();
        let d = Demuxer::new(&sh, cfg()).unwrap();
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
