use anyhow::{bail, Result};

use super::indel::{find_index, indel_safety, IndelIndex, IndelSafety, IndexMatch};
use super::types::{DemuxCfg, DemuxMode, DemuxVia, IndexEvidence, IndexSource, Outcome};
use crate::anchors::{revcomp, P5, P7, R1_PRIMER, R2_PRIMER};
use crate::hamming_index::HammingIndex;
use crate::sample_sheet::SampleSheet;

/// Real minimum pairwise Hamming distance among the sheet's samples, for the
/// given index alone.
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
            let d = seqs[i].iter().zip(seqs[j]).filter(|(a, b)| a != b).count() as u32;
            min_dist = min_dist.min(d);
        }
    }
    if seqs.len() < 2 {
        0
    } else {
        min_dist
    }
}

/// Matches extracted index sequences to samples.
pub struct Demuxer {
    sample_ids: Vec<String>,
    index_len: usize,
    i5_index: HammingIndex<usize>,
    i7_index: HammingIndex<usize>,
    /// Read-oriented keys, indexed like `sample_ids` (`revcomp(i5)`, `i7`).
    i5_keys: Vec<Vec<u8>>,
    i7_keys: Vec<Vec<u8>>,
    /// One-indel tables, set by [`Demuxer::enable_indels`] per safe index.
    i5_indel: Option<IndelIndex>,
    i7_indel: Option<IndelIndex>,
    resolved_mode: DemuxMode,
}

/// The P5 adapter's 3' end as read before i5 in a canonical read
/// (`...CGAGAT` + `CTACAC` + i5), seen in every CAT02 read checked.
const P5_INDEX_TAIL: &[u8] = b"CTACAC";

impl Demuxer {
    /// Build from a sample sheet. The read-oriented index sequences are
    /// `revcomp(i5)` and `i7` (forward), each expanded to ≤`index_max_mm`
    /// neighbors for O(1), error-tolerant lookup.
    ///
    /// `cfg.mode`: `Some(I5)`/`Some(I7)` refuse (rather than silently risk
    /// ambiguous assignment) unless that index alone separates every sample
    /// with enough real Hamming margin for `cfg.index_max_mm`; `Some(Dual)`/
    /// `Some(Both)` are taken as given (`Both`'s safety only needs the
    /// declared (i5,i7) *pairs* to be distinct, true of any valid sheet, not
    /// either index alone). `None` (auto) resolves to `Dual` if both indices
    /// are standalone-safe, else whichever single index is safe, else `Both`
    /// -- always logged, so the choice is never silent.
    pub fn new(sheet: &SampleSheet, cfg: DemuxCfg) -> Result<Self> {
        let i5_dist = min_pairwise_hamming(sheet, IndexSource::I5);
        let i7_dist = min_pairwise_hamming(sheet, IndexSource::I7);
        let required = 2 * cfg.index_max_mm + 1;
        let i5_safe = i5_dist >= required;
        let i7_safe = i7_dist >= required;
        let auto = cfg.mode.is_none();

        let resolved_mode = match cfg.mode {
            Some(DemuxMode::I5) if !i5_safe => bail!(unsafe_single_index_msg(
                DemuxMode::I5,
                i5_dist,
                cfg.index_max_mm,
                required
            )),
            Some(DemuxMode::I7) if !i7_safe => bail!(unsafe_single_index_msg(
                DemuxMode::I7,
                i7_dist,
                cfg.index_max_mm,
                required
            )),
            Some(explicit) => explicit,
            None if i5_safe && i7_safe => DemuxMode::Dual,
            None if i5_safe => DemuxMode::I5,
            None if i7_safe => DemuxMode::I7,
            None => DemuxMode::Both,
        };
        tracing::info!(
            mode = ?resolved_mode,
            auto,
            i5_min_hamming = i5_dist,
            i7_min_hamming = i7_dist,
            index_max_mm = cfg.index_max_mm,
            "demux mode resolved"
        );

        let mut i5_index = HammingIndex::new();
        let mut i7_index = HammingIndex::new();
        let index_len = sheet.samples.first().map_or(10, |s| s.i7.len());
        let i5_keys: Vec<Vec<u8>> = sheet.samples.iter().map(|s| revcomp(&s.i5)).collect();
        let i7_keys: Vec<Vec<u8>> = sheet.samples.iter().map(|s| s.i7.clone()).collect();
        for idx in 0..sheet.samples.len() {
            i5_index.insert(&i5_keys[idx], idx, cfg.index_max_mm);
            i7_index.insert(&i7_keys[idx], idx, cfg.index_max_mm);
        }
        Ok(Demuxer {
            sample_ids: sheet.samples.iter().map(|s| s.id.clone()).collect(),
            index_len,
            i5_index,
            i7_index,
            i5_keys,
            i7_keys,
            i5_indel: None,
            i7_indel: None,
            resolved_mode,
        })
    }

    /// Turn on one-indel matching (`--index-indels`) for each index the sheet
    /// can support safely at this `window` (see [`indel_safety`]); an unsafe
    /// index stays mismatch-only, with a warning. Returns `(i5, i7)` safety.
    pub fn enable_indels(&mut self, window: isize) -> (IndelSafety, IndelSafety) {
        let build = |keys: &[Vec<u8>]| {
            let mut ix = IndelIndex::new(self.index_len);
            for (i, k) in keys.iter().enumerate() {
                ix.insert(k, i);
            }
            ix
        };
        let i5_left = [P5, P5_INDEX_TAIL].concat();
        let i7_left = revcomp(R2_PRIMER);
        let i7_right = revcomp(P7);
        let i5 = build(&self.i5_keys);
        let i7 = build(&self.i7_keys);
        let s5 = indel_safety(
            &self.i5_keys,
            &self.i5_index,
            &i5,
            &i5_left,
            R1_PRIMER,
            window,
        );
        let s7 = indel_safety(
            &self.i7_keys,
            &self.i7_index,
            &i7,
            &i7_left,
            &i7_right,
            window,
        );
        for (name, safety) in [("i5", s5), ("i7", s7)] {
            if safety.safe() {
                tracing::info!(
                    index = name,
                    min_edit_distance = safety.min_edit,
                    one_edit_variants_recovered =
                        format!("{}/{}", safety.recovered, safety.variants),
                    one_edit_variants_to_wrong_sample = safety.wrong,
                    without_indels = safety.wrong_mismatch_only,
                    "index indel matching enabled (wrong-sample variants all come from the \
                     mismatch-only frame float, unchanged)"
                );
            } else {
                tracing::warn!(
                    index = name,
                    min_edit_distance = safety.min_edit,
                    one_edit_variants_to_wrong_sample = safety.wrong,
                    without_indels = safety.wrong_mismatch_only,
                    "index indel matching NOT enabled for this index: two sample indexes are \
                     too close (needs min edit distance >= 3, and the indel pass must add no \
                     wrong-sample call); it stays mismatch-only"
                );
            }
        }
        self.i5_indel = s5.safe().then_some(i5);
        self.i7_indel = s7.safe().then_some(i7);
        (s5, s7)
    }

    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    pub fn index_len(&self) -> usize {
        self.index_len
    }

    pub fn resolved_mode(&self) -> DemuxMode {
        self.resolved_mode
    }

    /// Whether one-indel matching is active for (i5, i7).
    pub fn indels_enabled(&self) -> (bool, bool) {
        (self.i5_indel.is_some(), self.i7_indel.is_some())
    }

    /// Combine per-index sample matches into an outcome (`Dual`/`Both` modes only).
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
            (Some(x), None) if self.resolved_mode == DemuxMode::Dual => Outcome {
                sample: Some(x),
                via: DemuxVia::I5Only,
            },
            (None, Some(y)) if self.resolved_mode == DemuxMode::Dual => Outcome {
                sample: Some(y),
                via: DemuxVia::I7Only,
            },
            _ => Outcome {
                sample: None,
                via: DemuxVia::None,
            },
        }
    }

    /// One index's outcome, standing alone (`I5`/`I7` mode): a hit is
    /// reported the same way dual-index reports a fallback single-index hit
    /// (`I5Only`/`I7Only`), since semantically it's the same situation --
    /// only one index was ever consulted.
    fn solo(sample: Option<usize>, via: DemuxVia) -> Outcome {
        match sample {
            Some(x) => Outcome {
                sample: Some(x),
                via,
            },
            None => Outcome {
                sample: None,
                via: DemuxVia::None,
            },
        }
    }

    /// Assign a read from its (optional) i5/i7 read-region slices (fixed frame).
    pub fn assign(&self, i5_region: Option<&[u8]>, i7_region: Option<&[u8]>) -> Outcome {
        match self.resolved_mode {
            DemuxMode::I5 => Self::solo(
                i5_region.and_then(|r| self.i5_index.get(r)),
                DemuxVia::I5Only,
            ),
            DemuxMode::I7 => Self::solo(
                i7_region.and_then(|r| self.i7_index.get(r)),
                DemuxVia::I7Only,
            ),
            DemuxMode::Dual | DemuxMode::Both => {
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
        self.assign_windowed_detail(seq, i5_start, i7_start, window)
            .0
    }

    /// [`assign_windowed`] plus the per-index evidence behind it: whether each
    /// index region is in the read at all ("seen") and what it matched
    /// ("identified"). An index the resolved mode never consults is reported
    /// unmatched.
    pub fn assign_windowed_detail(
        &self,
        seq: &[u8],
        i5_start: Option<usize>,
        i7_start: Option<usize>,
        window: isize,
    ) -> (Outcome, IndexEvidence) {
        let len = self.index_len;
        let seen = |s: Option<usize>| s.is_some_and(|s| s + len <= seq.len());
        let use_i5 = self.resolved_mode != DemuxMode::I7;
        let use_i7 = self.resolved_mode != DemuxMode::I5;
        let i5 = if use_i5 {
            find_index(
                &self.i5_index,
                self.i5_indel.as_ref(),
                seq,
                i5_start,
                len,
                window,
            )
        } else {
            IndexMatch::default()
        };
        let i7 = if use_i7 {
            find_index(
                &self.i7_index,
                self.i7_indel.as_ref(),
                seq,
                i7_start,
                len,
                window,
            )
        } else {
            IndexMatch::default()
        };
        let outcome = match self.resolved_mode {
            DemuxMode::I5 => Self::solo(i5.sample, DemuxVia::I5Only),
            DemuxMode::I7 => Self::solo(i7.sample, DemuxVia::I7Only),
            DemuxMode::Dual | DemuxMode::Both => self.combine(i5.sample, i7.sample),
        };
        // A match found in a floated frame counts as seen even when the
        // nominal frame runs off the read end, so identified ⊆ seen.
        let evidence = IndexEvidence {
            i5_seen: seen(i5_start) || i5.sample.is_some(),
            i7_seen: seen(i7_start) || i7.sample.is_some(),
            i5,
            i7,
        };
        (outcome, evidence)
    }
}

fn unsafe_single_index_msg(mode: DemuxMode, dist: u32, mm: u32, required: u32) -> String {
    format!(
        "--demux-index {mode:?} is unsafe for this sample sheet: real min pairwise Hamming \
         distance is {dist}, but unambiguous decoding at --index-max-mm {mm} needs > {req}. \
         Use --demux-index dual, both, or lower --index-max-mm.",
        req = required - 1,
    )
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
