use anyhow::{bail, Result};

use crate::anchors::{revcomp, R1_PRIMER, R2_PRIMER};
use crate::hamming_index::HammingIndex;
use crate::sample_sheet::SampleSheet;
use super::types::{DemuxCfg, DemuxMode, DemuxVia, IndexSource, Outcome};

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
    resolved_mode: DemuxMode,
}

impl Demuxer {
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
        for (idx, s) in sheet.samples.iter().enumerate() {
            i5_index.insert(&revcomp(&s.i5), idx, cfg.index_max_mm);
            i7_index.insert(&s.i7, idx, cfg.index_max_mm);
        }
        Ok(Demuxer {
            sample_ids: sheet.samples.iter().map(|s| s.id.clone()).collect(),
            index_len,
            i5_index,
            i7_index,
            resolved_mode,
        })
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

    pub fn assign_windowed(
        &self,
        seq: &[u8],
        i5_start: Option<usize>,
        i7_start: Option<usize>,
        window: isize,
    ) -> Outcome {
        match self.resolved_mode {
            DemuxMode::I5 => Self::solo(
                self.i5_index
                    .find_windowed(seq, i5_start, self.index_len, window),
                DemuxVia::I5Only,
            ),
            DemuxMode::I7 => Self::solo(
                self.i7_index
                    .find_windowed(seq, i7_start, self.index_len, window),
                DemuxVia::I7Only,
            ),
            DemuxMode::Dual | DemuxMode::Both => {
                let a5 = self
                    .i5_index
                    .find_windowed(seq, i5_start, self.index_len, window);
                let a7 = self
                    .i7_index
                    .find_windowed(seq, i7_start, self.index_len, window);
                self.combine(a5, a7)
            }
        }
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

pub fn index_anchor_starts(
    cb_start: usize,
    adapter_pos: Option<usize>,
    index_len: usize,
) -> (Option<usize>, Option<usize>) {
    let i5 = cb_start.checked_sub(R1_PRIMER.len() + index_len);
    let i7 = adapter_pos.map(|p| p + R2_PRIMER.len());
    (i5, i7)
}

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
