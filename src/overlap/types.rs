use crate::yc::{EncodeStats, YcTag};

/// Phred offset for Illumina FASTQ qualities.
pub const QUAL_OFFSET: u8 = 33;

/// Tuning for `overlap`. Defaults follow `hts_Overlapper` except `min_overlap`.
#[derive(Debug, Clone, Copy)]
pub struct OverlapCfg {
    /// Seed length for the R1 lookup table.
    pub kmer: usize,
    /// Stride between indexed seeds. 1 indexes every position.
    pub kmer_offset: usize,
    /// How far into each end of `rc(R2)` to probe for a seed.
    pub check_lengths: usize,
    /// Shortest overlap that may be merged.
    pub min_overlap: usize,
    /// Absolute cap on mismatches within an overlap.
    pub max_mismatch: usize,
    /// Cap on mismatches as a fraction of the overlap length.
    pub max_mismatch_density: f64,
    /// Ceiling for summed quality at concordant positions, in Phred.
    pub qual_cap: u8,
}

impl Default for OverlapCfg {
    fn default() -> Self {
        Self {
            kmer: 8,
            kmer_offset: 1,
            check_lengths: 20,
            // NOT HTStream's 8. A false merge yields a chimeric read that maps
            // to the wrong place, which is strictly worse than not merging.
            min_overlap: 20,
            max_mismatch: 100,
            max_mismatch_density: 0.25,
            qual_cap: 40,
        }
    }
}

/// How the fragment length compares to the reads that produced it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InsertClass {
    /// At least one mate ran past the end of the molecule — R1 into adapter, or
    /// R2 back into the R1 primer. The merge has trimmed the read-through off.
    Short,
    /// The mates met without either running off the end of the molecule.
    Medium,
    /// Reads did not overlap.
    Long,
}

/// A merged fragment, in R1 orientation.
#[derive(Debug, Clone)]
pub struct Merged {
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    pub class: InsertClass,
    /// Positions covered by both mates.
    pub overlap_len: usize,
    /// Simplex bases at the 5' end, contributed by R1.
    pub left_overhang: usize,
    /// Simplex bases at the 3' end, contributed by R2.
    pub right_overhang: usize,
    /// Positions inside the overlap where the mates disagreed.
    pub discordant: usize,
    /// Lossless description of the merge, in the same orientation as `seq`.
    pub yc: Option<YcTag>,
    /// Choices the tag encoder had to make.
    pub yc_stats: EncodeStats,
}

/// Reusable buffers, so a per-shard worker can merge without allocating per read.
#[derive(Debug, Default)]
pub struct Scratch {
    pub(crate) index: Vec<(u64, u32)>,
    pub(crate) probes: Vec<(u32, u64)>,
    pub(crate) rc2: Vec<u8>,
    pub(crate) rc2q: Vec<u8>,
}

/// Two-bit code for a nucleotide, or `None` for anything else (`N` included).
pub fn base_code(base: u8) -> Option<u64> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Index every `kmer_offset`-th seed of `seq`, skipping seeds that contain a
/// non-ACGT base or that lie wholly inside a homopolymer run.
pub fn build_index(seq: &[u8], k: usize, offset: usize, out: &mut Vec<(u64, u32)>) {
    out.clear();
    if k == 0 || k > 32 || seq.len() < k {
        return;
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let mut value = 0u64;
    let mut valid = 0usize;
    let mut same = 1usize;
    for (i, &base) in seq.iter().enumerate() {
        match base_code(base) {
            Some(code) => {
                value = ((value << 2) | code) & mask;
                valid += 1;
            }
            None => {
                value = 0;
                valid = 0;
            }
        }
        same = if i > 0 && base.eq_ignore_ascii_case(&seq[i - 1]) {
            same + 1
        } else {
            1
        };
        if valid >= k {
            let start = i + 1 - k;
            if start.is_multiple_of(offset) && same < k {
                out.push((value, start as u32));
            }
        }
    }
    out.sort_unstable();
}

/// Encode the seed at `start`, or `None` if it is unusable as a probe.
pub fn usable_kmer(seq: &[u8], start: usize, k: usize) -> Option<u64> {
    let window = &seq[start..start + k];
    if window.iter().all(|&b| b.eq_ignore_ascii_case(&window[0])) {
        return None;
    }
    let mut value = 0u64;
    for &base in window {
        value = (value << 2) | base_code(base)?;
    }
    Some(value)
}

/// All indexed R1 positions carrying `key`.
pub fn seed_positions(index: &[(u64, u32)], key: u64) -> &[(u64, u32)] {
    let lo = index.partition_point(|&(k, _)| k < key);
    let hi = index.partition_point(|&(k, _)| k <= key);
    &index[lo..hi]
}
