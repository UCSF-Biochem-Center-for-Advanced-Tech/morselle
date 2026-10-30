use anyhow::{bail, Result};

use super::indel::{indel_safety, IndelSafety, IndexMatch};
use super::lut::{IndexLut, MAX_SAMPLES};
use super::types::{DemuxCfg, DemuxMode, DemuxVia, IndexEvidence, IndexSource, Outcome};
use crate::anchors::{revcomp, P5, P7, R1_PRIMER, R2_PRIMER};
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
    /// Every matchable variant per index, precomputed ([`IndexLut`]); the
    /// one-indel tables are added by [`Demuxer::enable_indels`] per safe index.
    i5_lut: IndexLut,
    i7_lut: IndexLut,
    /// Read-oriented keys, indexed like `sample_ids` (`revcomp(i5)`, `i7`).
    i5_keys: Vec<Vec<u8>>,
    i7_keys: Vec<Vec<u8>>,
    resolved_mode: DemuxMode,
}

/// The P5 adapter's 3' end as read before i5 in a canonical read
/// (`...CGAGAT` + `CTACAC` + i5), seen in every CAT02 read checked.
const P5_INDEX_TAIL: &[u8] = b"CTACAC";

impl Demuxer {
    /// Build from a sample sheet. The read-oriented index sequences are
    /// `revcomp(i5)` and `i7` (forward), each expanded to ≤`index_max_mm`
    /// neighbors in a precomputed table for O(1), error-tolerant lookup.
    ///
    /// `cfg.mode`: `Some(I5)`/`Some(I7)` refuse (rather than silently risk
    /// ambiguous assignment) unless that index alone separates every sample
    /// with enough real Hamming margin for `cfg.index_max_mm`; `Some(Dual)`/
    /// `Some(Both)` are taken as given (`Both`'s safety only needs the
    /// declared (i5,i7) *pairs* to be distinct, true of any valid sheet, not
    /// either index alone). `None` (auto) resolves to `Dual`, which decides
    /// per read: both indexes are looked up, one alone assigns only when it
    /// names exactly one sample, and otherwise the sample must be the single
    /// one both indexes agree on (see `combine`). Before this per-read rule,
    /// auto fell back to a single index or `Both` when an index was not
    /// standalone-safe across the whole sheet; a key shared by several
    /// samples now never assigns on its own, so `Dual` is safe on any sheet.
    /// The resolution is always logged, so the choice is never silent.
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
            None => DemuxMode::Dual,
        };
        tracing::info!(
            mode = ?resolved_mode,
            auto,
            i5_min_hamming = i5_dist,
            i7_min_hamming = i7_dist,
            i5_standalone_safe = i5_safe,
            i7_standalone_safe = i7_safe,
            index_max_mm = cfg.index_max_mm,
            "demux mode resolved"
        );

        if sheet.samples.len() > MAX_SAMPLES {
            bail!(
                "{} samples in the sheet; demultiplexing supports at most {MAX_SAMPLES}",
                sheet.samples.len()
            );
        }
        let index_len = sheet.samples.first().map_or(10, |s| s.i7.len());
        let i5_keys: Vec<Vec<u8>> = sheet.samples.iter().map(|s| revcomp(&s.i5)).collect();
        let i7_keys: Vec<Vec<u8>> = sheet.samples.iter().map(|s| s.i7.clone()).collect();
        Ok(Demuxer {
            sample_ids: sheet.samples.iter().map(|s| s.id.clone()).collect(),
            index_len,
            i5_lut: IndexLut::new(&i5_keys, cfg.index_max_mm),
            i7_lut: IndexLut::new(&i7_keys, cfg.index_max_mm),
            i5_keys,
            i7_keys,
            resolved_mode,
        })
    }

    /// Turn on one-indel matching (`--index-indels`) for each index the sheet
    /// can support safely at this `window` (see [`indel_safety`]); an unsafe
    /// index stays mismatch-only, with a warning. Returns `(i5, i7)` safety.
    pub fn enable_indels(&mut self, window: isize) -> (IndelSafety, IndelSafety) {
        self.i5_lut.add_indels();
        self.i7_lut.add_indels();
        let i5_left = [P5, P5_INDEX_TAIL].concat();
        let i7_left = revcomp(R2_PRIMER);
        let i7_right = revcomp(P7);
        let s5 = indel_safety(&self.i5_keys, &self.i5_lut, &i5_left, R1_PRIMER, window);
        let s7 = indel_safety(&self.i7_keys, &self.i7_lut, &i7_left, &i7_right, window);
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
        if !s5.safe() {
            self.i5_lut.drop_indels();
        }
        if !s7.safe() {
            self.i7_lut.drop_indels();
        }
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
        (self.i5_lut.has_indels(), self.i7_lut.has_indels())
    }

    /// Combine the two per-index matches (`Dual`/`Both` modes): the per-read
    /// rule. Each index gives no match, one sample, or a shared candidate
    /// set (a key several samples claim). With both matched, the read goes
    /// to the one sample in both sets (`Dual`); no common sample is a
    /// `Conflict`, several stay undetermined. With one matched (`Dual` mode
    /// only), it assigns only when it names exactly one sample.
    fn combine(&self, i5: &IndexMatch, i7: &IndexMatch) -> Outcome {
        let none = |via| Outcome { sample: None, via };
        match (i5.matched(), i7.matched()) {
            (true, true) => {
                let a = self.candidates(&self.i5_lut, i5);
                let b = self.candidates(&self.i7_lut, i7);
                match common(&a, &b) {
                    Common::None => none(DemuxVia::Conflict),
                    Common::One(x) => Outcome {
                        sample: Some(x),
                        via: DemuxVia::Dual,
                    },
                    Common::Several => none(DemuxVia::None),
                }
            }
            (true, false) if self.resolved_mode == DemuxMode::Dual => match i5.sample {
                Some(x) => Outcome {
                    sample: Some(x),
                    via: DemuxVia::I5Only,
                },
                None => none(DemuxVia::None),
            },
            (false, true) if self.resolved_mode == DemuxMode::Dual => match i7.sample {
                Some(y) => Outcome {
                    sample: Some(y),
                    via: DemuxVia::I7Only,
                },
                None => none(DemuxVia::None),
            },
            _ => none(DemuxVia::None),
        }
    }

    /// One matched index's candidate samples, sorted.
    fn candidates<'a>(&self, lut: &'a IndexLut, m: &IndexMatch) -> CandSet<'a> {
        match (m.sample, m.shared) {
            (Some(x), _) => CandSet::One(x as u16),
            (None, Some(g)) => CandSet::Many(lut.members(g)),
            (None, None) => CandSet::Many(&[]),
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
                i5_region.and_then(|r| self.i5_lut.get(r).sample),
                DemuxVia::I5Only,
            ),
            DemuxMode::I7 => Self::solo(
                i7_region.and_then(|r| self.i7_lut.get(r).sample),
                DemuxVia::I7Only,
            ),
            DemuxMode::Dual | DemuxMode::Both => {
                let a5 = i5_region.map(|r| self.i5_lut.get(r)).unwrap_or_default();
                let a7 = i7_region.map(|r| self.i7_lut.get(r)).unwrap_or_default();
                self.combine(&a5, &a7)
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
            self.i5_lut.find(seq, i5_start, window, true)
        } else {
            IndexMatch::default()
        };
        let i7 = if use_i7 {
            self.i7_lut.find(seq, i7_start, window, true)
        } else {
            IndexMatch::default()
        };
        let outcome = match self.resolved_mode {
            DemuxMode::I5 => Self::solo(i5.sample, DemuxVia::I5Only),
            DemuxMode::I7 => Self::solo(i7.sample, DemuxVia::I7Only),
            DemuxMode::Dual | DemuxMode::Both => self.combine(&i5, &i7),
        };
        // A match found in a floated frame counts as seen even when the
        // nominal frame runs off the read end, so identified ⊆ seen.
        let evidence = IndexEvidence {
            i5_seen: seen(i5_start) || i5.matched(),
            i7_seen: seen(i7_start) || i7.matched(),
            i5,
            i7,
        };
        (outcome, evidence)
    }
}

/// One index's candidate samples (sorted).
enum CandSet<'a> {
    One(u16),
    Many(&'a [u16]),
}

impl CandSet<'_> {
    fn as_slice(&self) -> &[u16] {
        match self {
            CandSet::One(x) => std::slice::from_ref(x),
            CandSet::Many(m) => m,
        }
    }
}

enum Common {
    None,
    One(usize),
    Several,
}

/// How many samples two sorted candidate sets share (`One` names it).
fn common(a: &CandSet, b: &CandSet) -> Common {
    let (a, b) = (a.as_slice(), b.as_slice());
    let (mut i, mut j, mut hit) = (0, 0, None);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                if hit.is_some() {
                    return Common::Several;
                }
                hit = Some(a[i] as usize);
                i += 1;
                j += 1;
            }
        }
    }
    hit.map_or(Common::None, Common::One)
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
