//! Paired-end overlap merging for Illumina reads (plan: OVERLAP_PLAN.md §Phase 1).
//!
//! A port of HTStream's `hts_Overlapper` seed-and-extend algorithm, adapted for
//! the long-R1 10x GEX geometry.
//!
//! **Merge the raw mates, before any trimming.** The RT primer's dT block is the
//! polyA complement, so the molecule R1 reads and the molecule `rc(R2)` reads are
//! the same sequence:
//!
//! ```text
//! molecule = CB(16) | UMI(12) | T(30) | rc(cDNA)
//! R1       = molecule                        | adapter        (truncated)
//! R2       = rc(molecule)                    | rc(R1 primer)  (truncated)
//! rc(R2)   = R1 primer | CB | UMI | T(30) | rc(cDNA)
//! ```
//!
//! Merging before trimming is what makes the alignment gapless — strip CB+UMI+dT
//! from R1 first and the mates no longer describe the same span. It buys two more
//! things:
//!
//! * **Adapter removal by geometry, not by sequence matching.** The R1 primer
//!   sits *before* the CB, so `rc(R2)`'s 5' read-through is primer and `R1`'s 3'
//!   read-through is adapter — and both fall outside the overlap, so the merge
//!   drops them without ever consulting an adapter sequence. Each mate is
//!   independent evidence of where the fragment ends, which a Hamming-seeded
//!   trimmer like [`crate::extract::find_adapter`] can only approximate; it has
//!   a false-trim rate, and this has none. Running a matching trimmer over a
//!   merged read can therefore only add false positives.
//! * **Duplex confirmation of the barcode.** On short fragments R2 reads back
//!   through the UMI and CB, so CB+UMI lands inside the overlap and comes out
//!   with summed quality — and 28 bp of complex, per-read-unique seed material.
//!
//! This module takes the **raw** R1 and R2 and returns the merged molecule in
//! **R1 orientation** with a [`crate::yc`] tag describing the merge losslessly.
//! The caller then locates CB/UMI/cDNA in the merged read and trims the tag to
//! match.
//!
//! Deliberate divergences from `hts_Overlapper`, all covered by tests:
//!
//! * **Seeding is decoupled from the merge frame.** HTStream swaps the mates so
//!   R1 is the longer one, which also flips the merge frame. We cannot flip the
//!   frame — the merge must stay in R1's orientation or the `YC` tag's notion of
//!   which strand is R1 silently inverts. So we choose which mate to *index*
//!   independently: index the longer, probe the shorter. The shorter read's ends
//!   are always inside the overlap; `rc(R2)`'s need not be. In 10x geometry R1
//!   loses 28 bp of CB+UMI and the whole dT run, so its insert is routinely
//!   *contained* within `rc(R2)`, with overhangs on both sides — probing
//!   `rc(R2)`'s ends would find nothing at all.
//! * **Correct 3' extension offset.** HTStream appends `seq2.substr(maxLoop, ..)`,
//!   which is only right when the `rc(R2)` offset is zero — its swap normally
//!   guarantees that. Without the swap we must use `loc2 + maxLoop`.
//! * **Discordant quality** is `max(1, |q1 - q2|)`. HTStream writes
//!   `min(q1 - q2 + offset, 1 + offset)`, collapsing every disagreement to Q1.
//! * **Mismatch budget** uses [`crate::anchors::hamming_leq`] (`<= max`), which
//!   early-exits correctly. HTStream's `threshold_mismatches` decrements before
//!   testing and underflows `size_t` when the budget is zero.
//! * **Seed masking.** k-mers containing a non-ACGT base, or lying inside a
//!   homopolymer, are kept out of the index — the dT run abutting the R1 insert
//!   would otherwise seed against any poly-A in R2.

use crate::anchors::{comp, hamming_leq};
use crate::yc::{self, EncodeStats, YcTag};

/// Phred offset for Illumina FASTQ qualities.
const QUAL_OFFSET: u8 = 33;

/// Tuning for [`overlap`]. Defaults follow `hts_Overlapper` except `min_overlap`.
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
    /// Reads did not overlap. Never carried by a [`Merged`] — [`overlap`] returns
    /// `None` in this case; the variant exists so callers can use one enum for
    /// their counters.
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
    /// `None` only if the overlap contained a base no `YC` code covers (an IUPAC
    /// ambiguity code other than `N`).
    pub yc: Option<YcTag>,
    /// Choices the tag encoder had to make.
    pub yc_stats: EncodeStats,
}

/// Reusable buffers, so a per-shard worker can merge without allocating per read.
#[derive(Debug, Default)]
pub struct Scratch {
    index: Vec<(u64, u32)>,
    probes: Vec<(u32, u64)>,
    rc2: Vec<u8>,
    rc2q: Vec<u8>,
}

/// Two-bit code for a nucleotide, or `None` for anything else (`N` included).
fn base_code(base: u8) -> Option<u64> {
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
fn build_index(seq: &[u8], k: usize, offset: usize, out: &mut Vec<(u64, u32)>) {
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

/// Encode the seed at `start`, or `None` if it is unusable as a probe: it holds
/// a non-ACGT base, or it is a homopolymer. Homopolymer seeds are never indexed,
/// so probing one can only waste a lookup — and the dT block guarantees a long
/// run of them right where the mates meet.
fn usable_kmer(seq: &[u8], start: usize, k: usize) -> Option<u64> {
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
fn seed_positions(index: &[(u64, u32)], key: u64) -> &[(u64, u32)] {
    let lo = index.partition_point(|&(k, _)| k < key);
    let hi = index.partition_point(|&(k, _)| k <= key);
    &index[lo..hi]
}

/// Merge `r1` with `rc(r2)` anchored so that `r1[loc1]` aligns to `rc2[loc2]`,
/// or `None` if the implied overlap is too short or too divergent.
#[allow(clippy::too_many_arguments)]
fn check_overlap(
    r1: &[u8],
    r1q: &[u8],
    rc2: &[u8],
    rc2q: &[u8],
    loc1: usize,
    loc2: usize,
    cfg: &OverlapCfg,
) -> Option<Merged> {
    // Slide both reads back so the earlier one starts at zero.
    let min_loc = loc1.min(loc2);
    let start1 = loc1 - min_loc;
    let start2 = loc2 - min_loc;

    let reach1 = r1.len() - start1;
    let reach2 = rc2.len() - start2;
    let span = reach1.min(reach2);
    if span < cfg.min_overlap {
        return None;
    }

    let budget = cfg
        .max_mismatch
        .min((span as f64 * cfg.max_mismatch_density) as usize);
    if !hamming_leq(
        &r1[start1..start1 + span],
        &rc2[start2..start2 + span],
        budget as u32,
    ) {
        return None;
    }

    // 5' overhang: kept only when it is R1 running past the start of rc(R2).
    // The mirror case is rc(R2)'s adapter read-through, which we drop.
    let left_overhang = if start1 > start2 { start1 } else { 0 };
    // 3' overhang: kept only when rc(R2) runs past the end of R1, likewise.
    let right_overhang = if reach1 < reach2 { reach2 - span } else { 0 };

    let mut seq = Vec::with_capacity(left_overhang + span + right_overhang);
    let mut qual = Vec::with_capacity(seq.capacity());

    seq.extend_from_slice(&r1[..left_overhang]);
    qual.extend_from_slice(&r1q[..left_overhang]);

    let mut discordant = 0usize;
    for i in 0..span {
        let (a, b) = (r1[start1 + i], rc2[start2 + i]);
        let qa = r1q[start1 + i].saturating_sub(QUAL_OFFSET) as u16;
        let qb = rc2q[start2 + i].saturating_sub(QUAL_OFFSET) as u16;
        if a.eq_ignore_ascii_case(&b) {
            // Independent observations agree: their evidence adds, up to the cap.
            seq.push(a);
            qual.push((qa + qb).min(cfg.qual_cap as u16) as u8 + QUAL_OFFSET);
        } else {
            discordant += 1;
            let (base, delta) = if qa >= qb { (a, qa - qb) } else { (b, qb - qa) };
            seq.push(base);
            qual.push(delta.max(1) as u8 + QUAL_OFFSET);
        }
    }

    if right_overhang > 0 {
        // `start2 + span`, not `span`: without HTStream's mate swap the rc(R2)
        // offset need not be zero.
        seq.extend_from_slice(&rc2[start2 + span..]);
        qual.extend_from_slice(&rc2q[start2 + span..]);
    }

    // HTStream classifies by comparing the merge to the longest input read. That
    // is unreliable here because the mates are not symmetric: R1 may be
    // sequenced shorter than R2, and R2 additionally covers the polyA and the
    // barcode. Use the geometric signal the classification actually means —
    // read-through happened iff *either* mate ran past the end of the molecule:
    //
    //   * `reach1 > reach2` — R1 ran past rc(R2)'s 3' end, into adapter.
    //   * `start2 > 0`      — rc(R2) began before R1's first base. R1's
    //     sequencing primer starts exactly at the molecule's 5' end, so
    //     anything rc(R2) contributes before it is R1-primer read-through.
    //
    // Checking only the first misses the asymmetric case: a 100 bp R1 that never
    // reaches the adapter, paired with a 150 bp R2 that reads clean through the
    // CB. The fragment is plainly short; only `start2` says so.
    let class = if start2 > 0 || reach1 > reach2 {
        InsertClass::Short
    } else {
        InsertClass::Medium
    };

    let (tag, yc_stats) = match yc::encode_overlap(
        left_overhang as u64,
        true,
        &r1[start1..start1 + span],
        &rc2[start2..start2 + span],
        right_overhang as u64,
        false,
    ) {
        Ok((tag, stats)) => (Some(tag), stats),
        // Only reachable for IUPAC ambiguity codes beyond `N`. Keep the merge,
        // drop the tag rather than losing the read.
        Err(_) => (None, EncodeStats::default()),
    };

    Some(Merged {
        seq,
        qual,
        class,
        overlap_len: span,
        left_overhang,
        right_overhang,
        discordant,
        yc: tag,
        yc_stats,
    })
}

/// Merge a mate pair, allocating fresh buffers. See [`overlap_in`].
pub fn overlap(r1: &[u8], r1q: &[u8], r2: &[u8], r2q: &[u8], cfg: &OverlapCfg) -> Option<Merged> {
    overlap_in(&mut Scratch::default(), r1, r1q, r2, r2q, cfg)
}

/// Merge a mate pair, reusing `scratch`.
///
/// `r1` is the R1 **insert** — the caller has already stripped CB/UMI and the dT
/// run. Returns the fragment in R1 orientation, or `None` when the mates do not
/// overlap (a "long" insert), which the caller should emit as R2 alone.
pub fn overlap_in(
    scratch: &mut Scratch,
    r1: &[u8],
    r1q: &[u8],
    r2: &[u8],
    r2q: &[u8],
    cfg: &OverlapCfg,
) -> Option<Merged> {
    if r1.len() != r1q.len() || r2.len() != r2q.len() {
        return None;
    }
    if r1.len() < cfg.min_overlap || r2.len() < cfg.min_overlap {
        return None;
    }

    // Disjoint field borrows, so the index can be filled while rc(R2) is read.
    let Scratch {
        index,
        probes,
        rc2,
        rc2q,
    } = scratch;

    // rc(R2) shares R1's strand, so the two can be compared directly.
    rc2.clear();
    rc2.extend(r2.iter().rev().map(|&base| comp(base)));
    rc2q.clear();
    rc2q.extend(r2q.iter().rev().copied());
    let rc2: &[u8] = rc2;
    let rc2q: &[u8] = rc2q;

    let k = cfg.kmer.min(r1.len()).min(rc2.len());
    if k == 0 {
        return None;
    }
    let offset = cfg.kmer_offset.max(1);

    // Index the longer mate and probe the shorter one's ends. Only the shorter
    // read is guaranteed to have an end inside the overlap.
    let index_is_r1 = r1.len() >= rc2.len();
    let (indexed, probed) = if index_is_r1 { (r1, rc2) } else { (rc2, r1) };
    build_index(indexed, k, offset, index);
    if index.is_empty() {
        return None;
    }

    // Collect up to `check_lengths` *usable* seeds inward from each end.
    // HTStream probes a fixed number of positions, which fails here: the dT
    // block puts a long homopolymer exactly where the mates meet, so a
    // positional window can be entirely unindexable. Counting usable seeds
    // instead walks past the run without unbounding the work.
    let budget = cfg.check_lengths.max(1);
    let last_start = probed.len() - k;
    probes.clear();
    let mut found = 0usize;
    for probe in 0..=last_start {
        if found >= budget {
            break;
        }
        if let Some(key) = usable_kmer(probed, probe, k) {
            probes.push((probe as u32, key));
            found += 1;
        }
    }
    found = 0;
    for probe in (0..=last_start).rev() {
        if found >= budget {
            break;
        }
        if let Some(key) = usable_kmer(probed, probe, k) {
            probes.push((probe as u32, key));
            found += 1;
        }
    }

    for &(probe, key) in probes.iter() {
        let probe = probe as usize;
        for &(_, hit) in seed_positions(index, key) {
            // Map the seed back into (R1, rc(R2)) coordinates. The merge always
            // happens in R1's frame regardless of which mate was indexed.
            let (loc1, loc2) = if index_is_r1 {
                (hit as usize, probe)
            } else {
                (probe, hit as usize)
            };
            if let Some(merged) = check_overlap(r1, r1q, rc2, rc2q, loc1, loc2, cfg) {
                return Some(merged);
            }
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchors::revcomp;
    use crate::yc::decode;

    /// A 60 bp fragment with no long homopolymers and unique 8-mers.
    const FRAGMENT: &[u8] = b"ACGTTGCAAGCTTAGCCATGGTACCAGTCTGAACGTACGGTCATGGCGTTACAGCTAGCA";
    /// The 5' end of the universal Illumina 3' read-through adapter.
    const ADAPTER: &[u8] = b"AGATCGGAAG";

    /// R1 quality Q40, R2 quality Q20, so R1 wins every disagreement.
    const Q40: u8 = b'I';
    const Q20: u8 = b'5';

    fn quals(base: u8, len: usize) -> Vec<u8> {
        vec![base; len]
    }

    /// Build a mate pair from the R1 insert and the desired `rc(R2)`.
    fn mates(r1: &[u8], rc2: &[u8]) -> (Vec<u8>, Vec<u8>) {
        (r1.to_vec(), revcomp(rc2))
    }

    #[test]
    fn fixture_is_the_length_the_tests_assume() {
        assert_eq!(FRAGMENT.len(), 60);
        assert_eq!(ADAPTER.len(), 10);
    }

    #[test]
    fn merges_a_medium_insert_into_the_whole_fragment() {
        let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q40, r2.len()),
            &OverlapCfg::default(),
        )
        .expect("mates overlap by 20 bp");

        assert_eq!(merged.seq, FRAGMENT);
        assert_eq!(merged.class, InsertClass::Medium);
        assert_eq!(merged.overlap_len, 20);
        assert_eq!(merged.left_overhang, 20);
        assert_eq!(merged.right_overhang, 20);
        assert_eq!(merged.discordant, 0);
        assert_eq!(merged.yc.unwrap().to_string(), "20+20-20");
    }

    #[test]
    fn merges_a_short_insert_and_drops_read_through_adapter() {
        // Fragment shorter than the reads: both mates run into adapter.
        let insert = &FRAGMENT[..30];
        let mut r1 = insert.to_vec();
        r1.extend_from_slice(ADAPTER);
        let mut rc2 = revcomp(ADAPTER);
        rc2.extend_from_slice(insert);
        let r2 = revcomp(&rc2);

        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q40, r2.len()),
            &OverlapCfg::default(),
        )
        .expect("mates overlap across the whole insert");

        // Adapter is gone from both ends, purely from the overlap geometry.
        assert_eq!(merged.seq, insert);
        assert_eq!(merged.class, InsertClass::Short);
        assert_eq!(merged.left_overhang, 0);
        assert_eq!(merged.right_overhang, 0);
        assert_eq!(merged.yc.unwrap().to_string(), "+30+");
    }

    #[test]
    fn handles_a_mate_longer_than_r1() {
        // R2 longer than R1. HTStream swaps here; we must not, and still align.
        let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[10..60]);
        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q40, r2.len()),
            &OverlapCfg::default(),
        )
        .expect("30 bp overlap");

        assert_eq!(merged.seq, FRAGMENT);
        assert_eq!(merged.overlap_len, 30);
        assert_eq!(merged.left_overhang, 10);
        assert_eq!(merged.right_overhang, 20);
        assert_eq!(merged.yc.unwrap().to_string(), "10+30-20");
    }

    #[test]
    fn extends_from_the_correct_offset_when_rc_r2_starts_first() {
        // rc(R2) begins with adapter *and* runs past R1's 3' end. This is the
        // case HTStream mis-indexes: it would append rc2[span..] instead of
        // rc2[start2 + span..], duplicating 10 bp of the fragment.
        let mut rc2 = revcomp(ADAPTER);
        rc2.extend_from_slice(&FRAGMENT[..50]);
        let r1 = FRAGMENT[..30].to_vec();
        let r2 = revcomp(&rc2);

        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q40, r2.len()),
            &OverlapCfg::default(),
        )
        .expect("30 bp overlap");

        assert_eq!(merged.seq, &FRAGMENT[..50]);
        assert_eq!(merged.seq.len(), 50);
        assert_eq!(merged.left_overhang, 0);
        assert_eq!(merged.right_overhang, 20);
        assert_eq!(merged.yc.unwrap().to_string(), "+30-20");
    }

    // ── Quality ─────────────────────────────────────────────────────────────

    #[test]
    fn agreeing_bases_sum_quality_up_to_the_cap() {
        let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q20, r2.len()),
            &OverlapCfg::default(),
        )
        .unwrap();

        // R1-only overhang keeps R1's quality, the overlap is capped at Q40,
        // and the R2-only overhang keeps R2's.
        assert!(merged.qual[..20].iter().all(|&q| q == Q40));
        assert!(merged.qual[20..40].iter().all(|&q| q == Q40));
        assert!(merged.qual[40..].iter().all(|&q| q == Q20));
    }

    #[test]
    fn disagreeing_bases_take_the_better_call_and_the_quality_delta() {
        let mut rc2 = FRAGMENT[20..60].to_vec();
        // Corrupt one overlap position, away from the seed window.
        rc2[10] = if rc2[10] == b'A' { b'C' } else { b'A' };
        let r1 = FRAGMENT[..40].to_vec();
        let r2 = revcomp(&rc2);

        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q20, r2.len()),
            &OverlapCfg::default(),
        )
        .unwrap();

        // R1 was the higher-quality call, so the fragment is still recovered...
        assert_eq!(merged.seq, FRAGMENT);
        assert_eq!(merged.discordant, 1);
        // ...and that base carries Q40 - Q20 = Q20, not HTStream's flat Q1.
        assert_eq!(merged.qual[30], Q20);
        assert_eq!(merged.qual[29], Q40);
        // The tag records what each strand actually read.
        assert_eq!(merged.yc.unwrap().count_discordant(), 1);
    }

    // ── Rejection ───────────────────────────────────────────────────────────

    #[test]
    fn unrelated_reads_do_not_merge() {
        let r1 = b"ACGTTGCAAGCTTAGCCATGGTACCAGTCTGA".to_vec();
        let r2 = b"TTTGGCCATATCAGCATCAGGATTCCGGATAC".to_vec();
        assert!(overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q40, r2.len()),
            &OverlapCfg::default()
        )
        .is_none());
    }

    #[test]
    fn overlap_below_min_overlap_is_refused() {
        // Mates share only 10 bp.
        let (r1, r2) = mates(&FRAGMENT[..30], &FRAGMENT[20..50]);
        let r1q = quals(Q40, r1.len());
        let r2q = quals(Q40, r2.len());

        assert!(overlap(&r1, &r1q, &r2, &r2q, &OverlapCfg::default()).is_none());

        let permissive = OverlapCfg {
            min_overlap: 10,
            ..OverlapCfg::default()
        };
        let merged = overlap(&r1, &r1q, &r2, &r2q, &permissive).unwrap();
        assert_eq!(merged.seq, &FRAGMENT[..50]);
        assert_eq!(merged.overlap_len, 10);
    }

    #[test]
    fn overlap_over_the_mismatch_budget_is_refused() {
        let mut rc2 = FRAGMENT[20..60].to_vec();
        // 20 bp overlap allows floor(20 * 0.25) = 5 mismatches; make 6. They sit
        // outside the leading seed window, so a seed is still found and the
        // rejection genuinely comes from the mismatch check.
        for base in &mut rc2[8..14] {
            *base = if *base == b'A' { b'C' } else { b'A' };
        }
        let r1 = FRAGMENT[..40].to_vec();
        let r2 = revcomp(&rc2);

        assert!(overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q40, r2.len()),
            &OverlapCfg::default()
        )
        .is_none());
    }

    #[test]
    fn mismatched_sequence_and_quality_lengths_are_refused() {
        let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
        assert!(overlap(
            &r1,
            &quals(Q40, 5),
            &r2,
            &quals(Q40, r2.len()),
            &OverlapCfg::default()
        )
        .is_none());
    }

    #[test]
    fn short_mates_do_not_underflow() {
        // rc(R2) shorter than check_lengths + kmer: the probe window arithmetic
        // that underflows in HTStream.
        let cfg = OverlapCfg {
            min_overlap: 5,
            ..OverlapCfg::default()
        };
        for len in 0..30usize {
            let r2: Vec<u8> = FRAGMENT[20..20 + len].to_vec();
            let r1 = FRAGMENT[..40].to_vec();
            let _ = overlap(&r1, &quals(Q40, r1.len()), &r2, &quals(Q40, len), &cfg);
        }
        // Also exercise a degenerate R1.
        for len in 0..12usize {
            let r1: Vec<u8> = FRAGMENT[..len].to_vec();
            let r2 = revcomp(&FRAGMENT[20..60]);
            let _ = overlap(&r1, &quals(Q40, len), &r2, &quals(Q40, r2.len()), &cfg);
        }
    }

    // ── Seed masking ────────────────────────────────────────────────────────

    #[test]
    fn homopolymer_seeds_are_masked() {
        let mut index = Vec::new();
        build_index(b"AAAAAAAAAAAAAAAAAAAA", 8, 1, &mut index);
        assert!(index.is_empty(), "a pure homopolymer must index nothing");

        // The dT run abutting the R1 insert must not seed; the cDNA after it must.
        let mut index = Vec::new();
        build_index(b"TTTTTTTTTTTTACGTTGCAAGCTTAGC", 8, 1, &mut index);
        assert!(index.iter().all(|&(_, pos)| pos >= 5));
    }

    #[test]
    fn seeds_containing_non_acgt_are_masked() {
        let mut index = Vec::new();
        build_index(b"ACGTTGCAAGCTTAGC", 8, 1, &mut index);
        let clean = index.len();

        let mut index = Vec::new();
        build_index(b"ACGTNGCAAGCTTAGC", 8, 1, &mut index);
        // Every window covering the N is gone.
        assert_eq!(index.len(), clean - 5);
        assert!(index.iter().all(|&(_, pos)| pos > 4));
    }

    #[test]
    fn kmer_offset_strides_the_index() {
        let mut every = Vec::new();
        build_index(&FRAGMENT[..40], 8, 1, &mut every);
        let mut strided = Vec::new();
        build_index(&FRAGMENT[..40], 8, 4, &mut strided);
        assert_eq!(every.len(), 33);
        assert!(strided.iter().all(|&(_, pos)| pos % 4 == 0));
        assert!(strided.len() < every.len());
    }

    // ── The tag ─────────────────────────────────────────────────────────────

    #[test]
    fn the_tag_reconstructs_both_mates_from_the_merge() {
        let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q20, r2.len()),
            &OverlapCfg::default(),
        )
        .unwrap();

        let tag = merged.yc.as_ref().unwrap();
        assert_eq!(tag.sequence_length(), merged.seq.len() as u64);

        let decoded = decode(tag, &merged.seq, &merged.qual).unwrap();
        assert_eq!(decoded.r1_seq, r1, "R1 recovered from the merge");
        assert_eq!(
            decoded.r2_seq,
            revcomp(&r2),
            "rc(R2) recovered from the merge"
        );
    }

    #[test]
    fn the_tag_records_disagreements_losslessly() {
        let mut rc2 = FRAGMENT[20..60].to_vec();
        rc2[10] = if rc2[10] == b'A' { b'C' } else { b'A' };
        let r1 = FRAGMENT[..40].to_vec();
        let r2 = revcomp(&rc2);

        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &quals(Q20, r2.len()),
            &OverlapCfg::default(),
        )
        .unwrap();

        let decoded = decode(merged.yc.as_ref().unwrap(), &merged.seq, &merged.qual).unwrap();
        // Both original strands come back, disagreement and all.
        assert_eq!(decoded.r1_seq, r1);
        assert_eq!(decoded.r2_seq, rc2);
    }

    #[test]
    fn n_bases_are_coerced_and_counted() {
        let mut rc2 = FRAGMENT[20..60].to_vec();
        rc2[5] = b'N';
        let r1 = FRAGMENT[..40].to_vec();
        let r2 = revcomp(&rc2);
        let mut r2q = quals(Q40, r2.len());
        // needletail/Illumina give an N a low quality; R1 therefore wins.
        r2q[r2.len() - 1 - 5] = b'#';

        let merged = overlap(
            &r1,
            &quals(Q40, r1.len()),
            &r2,
            &r2q,
            &OverlapCfg::default(),
        )
        .expect("one N does not break the seed or the budget");

        assert_eq!(merged.seq, FRAGMENT);
        assert_eq!(merged.yc_stats.n_columns_coerced, 1);
        // Coerced to concordant, so the duplex segment stays one clean run.
        assert_eq!(merged.yc.unwrap().to_string(), "20+20-20");
    }

    #[test]
    fn scratch_buffers_are_reusable_across_pairs() {
        let mut scratch = Scratch::default();
        let cfg = OverlapCfg::default();
        for _ in 0..3 {
            let (r1, r2) = mates(&FRAGMENT[..40], &FRAGMENT[20..60]);
            let merged = overlap_in(
                &mut scratch,
                &r1,
                &quals(Q40, r1.len()),
                &r2,
                &quals(Q40, r2.len()),
                &cfg,
            )
            .unwrap();
            assert_eq!(merged.seq, FRAGMENT);
        }
    }
}

/// Synthetic 3' GEX mate pairs built to the real molecule structure, so the
/// merge is exercised against what actual data contains: the dT block sitting
/// exactly where the mates meet, internal homopolymers, and read-through past
/// the fragment into the opposite primer.
///
/// The RT primer's dT block *is* the polyA complement, so R1 and `rc(R2)`
/// describe one identical span — which is why the merge must run on the raw
/// mates, before CB/UMI/dT are stripped.
#[cfg(test)]
mod dataset {
    use crate::anchors::{revcomp, R1_PRIMER};
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    pub const CB: &[u8] = b"AAACCCAAGAAACACT";
    pub const UMI: &[u8] = b"ACGTACGTACGT";
    /// The dT block: RT-primer dT annealed to the polyA tail.
    pub const DT: usize = 30;
    /// CB(16) + UMI(12): where the dT run starts in R1 and in the merged read.
    pub const BARCODE_LEN: usize = 28;
    pub const ADAPTER: &[u8] = b"AGATCGGAAGAGC";

    /// Deterministic complex sequence — no run longer than 3.
    pub fn complex(len: usize, seed: u64) -> Vec<u8> {
        let bases = b"ACGT";
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let mut out: Vec<u8> = Vec::with_capacity(len);
        while out.len() < len {
            let base = bases[rng.gen_range(0..4)];
            let tail = out.len().saturating_sub(3);
            if out.len() >= 3 && out[tail..].iter().all(|&b| b == base) {
                continue;
            }
            out.push(base);
        }
        out
    }

    /// A cDNA core carrying four homopolymer runs of different bases. Ends on `C`
    /// so `revcomp` does not begin with `T` and blur the dT boundary.
    pub fn core_with_homopolymers() -> Vec<u8> {
        let mut core = Vec::new();
        core.extend_from_slice(&complex(24, 1));
        core.extend_from_slice(b"AAAAAAAAAA");
        core.extend_from_slice(&complex(22, 2));
        core.extend_from_slice(b"TTTTTTTTT");
        core.extend_from_slice(&complex(20, 3));
        core.extend_from_slice(b"GGGGGGGG");
        core.extend_from_slice(&complex(18, 4));
        core.extend_from_slice(b"CCCCCCC");
        core.extend_from_slice(&complex(18, 5));
        *core.last_mut().unwrap() = b'C';
        core
    }

    pub struct Pair {
        pub r1: Vec<u8>,
        pub r1q: Vec<u8>,
        pub r2: Vec<u8>,
        pub r2q: Vec<u8>,
        /// `CB | UMI | T(DT) | rc(core)` — what a complete merge yields.
        pub molecule: Vec<u8>,
        /// True when the mates reach each other at these read lengths.
        pub overlaps: bool,
    }

    /// Build a mate pair for a fragment with cDNA `core`, sequenced `r1_len` from
    /// the barcode end and `r2_len` from the cDNA end.
    pub fn build(core: &[u8], r1_len: usize, r2_len: usize) -> Pair {
        let mut molecule = Vec::new();
        molecule.extend_from_slice(CB);
        molecule.extend_from_slice(UMI);
        molecule.extend(std::iter::repeat_n(b'T', DT));
        molecule.extend_from_slice(&revcomp(core));

        let mut r1 = molecule.clone();
        r1.extend_from_slice(ADAPTER); // read-through past a short fragment
        r1.truncate(r1_len);

        let mut r2 = revcomp(&molecule); // core | A(DT) | rc(UMI) | rc(CB)
        r2.extend_from_slice(&revcomp(R1_PRIMER)); // read-through into the primer
        r2.truncate(r2_len);

        // In molecule coordinates R1 covers [0, a) and rc(R2) covers [M - b, M).
        let m = molecule.len();
        let a = r1_len.min(m);
        let b = r2_len.min(m);
        let overlaps = m.saturating_sub(b) < a;

        Pair {
            r1q: vec![b'I'; r1.len()],
            r2q: vec![b'I'; r2.len()],
            r1,
            r2,
            molecule,
            overlaps,
        }
    }
}

#[cfg(test)]
mod realistic {
    use super::dataset::*;
    use super::*;
    use crate::extract::dt_run_end;
    use crate::yc::decode;

    fn merge(pair: &Pair, cfg: &OverlapCfg) -> Option<Merged> {
        // The raw mates, untrimmed. This is the whole point.
        overlap(&pair.r1, &pair.r1q, &pair.r2, &pair.r2q, cfg)
    }

    #[test]
    fn merges_a_short_fragment_with_read_through_on_both_mates() {
        // 60 bp cDNA, 150 bp reads: R1 runs into adapter, R2 back into the primer.
        let core = complex(60, 11);
        let pair = build(&core, 150, 150);
        let merged = merge(&pair, &OverlapCfg::default()).expect("short fragment merges");

        assert_eq!(merged.seq, pair.molecule);
        assert_eq!(merged.class, InsertClass::Short);
        assert_eq!(merged.discordant, 0);
        // Both read-throughs are dropped purely by the overlap geometry.
        assert_eq!(merged.left_overhang, 0);
        assert_eq!(merged.right_overhang, 0);
    }

    #[test]
    fn merges_a_medium_fragment_into_the_whole_molecule() {
        // 130 bp cDNA: the mates meet across the tail of the dT block.
        let core = complex(130, 12);
        let pair = build(&core, 150, 150);
        let merged = merge(&pair, &OverlapCfg::default()).expect("medium fragment merges");

        assert_eq!(merged.seq, pair.molecule);
        assert_eq!(merged.class, InsertClass::Medium);
        // Merging the raw mates gives an overhang on each side, unlike the
        // pre-trimmed geometry where R1 could never start first.
        assert!(merged.left_overhang > 0);
        assert!(merged.right_overhang > 0);
    }

    #[test]
    fn a_long_fragment_does_not_merge() {
        let core = complex(400, 13);
        let pair = build(&core, 150, 150);
        assert!(!pair.overlaps, "fixture must not overlap");
        assert!(merge(&pair, &OverlapCfg::default()).is_none());
    }

    #[test]
    fn merges_when_r1_is_sequenced_shorter_than_r2() {
        // R2 spans the whole molecule while R1 stops inside it, so R1 is
        // strictly interior to rc(R2) — overhangs on both sides. Probing
        // rc(R2)'s ends would find nothing; we index it and probe R1 instead.
        let core = complex(60, 23);
        let pair = build(&core, 100, 150);
        assert!(pair.r1.len() < pair.r2.len());

        let merged = merge(&pair, &OverlapCfg::default()).expect("an interior R1 must still merge");
        assert_eq!(merged.seq, pair.molecule);
        assert_eq!(merged.left_overhang, 0);
        assert!(merged.right_overhang > 0);
    }

    /// Does `haystack` contain `needle` as a contiguous run?
    fn contains(haystack: &[u8], needle: &[u8]) -> bool {
        needle.len() <= haystack.len() && haystack.windows(needle.len()).any(|w| w == needle)
    }

    #[test]
    fn read_through_is_removed_by_geometry_not_by_matching() {
        // The R1 primer precedes the CB, so rc(R2)'s 5' read-through is primer
        // and R1's 3' read-through is adapter. Both sit outside the overlap.
        let core = complex(60, 26);
        let pair = build(&core, 150, 150);

        // The raw mates really do carry read-through, or this proves nothing.
        assert!(contains(&pair.r1, ADAPTER), "R1 must read into adapter");
        let rc2 = crate::anchors::revcomp(&pair.r2);
        assert!(
            contains(&rc2, &crate::anchors::R1_PRIMER[8..]),
            "rc(R2) must read back into the R1 primer"
        );

        let merged = merge(&pair, &OverlapCfg::default()).unwrap();
        assert_eq!(merged.seq, pair.molecule);
        assert!(!contains(&merged.seq, ADAPTER), "adapter is gone");
        assert!(
            !contains(&merged.seq, &crate::anchors::R1_PRIMER[8..]),
            "primer read-through is gone"
        );
        // Nothing was matched to achieve that — both overhangs were simply
        // outside the span the two mates agree on.
        assert_eq!(merged.left_overhang, 0);
        assert_eq!(merged.right_overhang, 0);
    }

    #[test]
    fn read_through_is_removed_even_when_the_adapter_is_misread() {
        // Corrupt the adapter badly enough to defeat a Hamming-seeded trimmer.
        // The geometry never looks at it, so the merge is unaffected.
        let core = complex(60, 27);
        let mut pair = build(&core, 150, 150);
        let end = pair.molecule.len();
        for offset in [0usize, 2, 4, 6, 8] {
            if end + offset < pair.r1.len() {
                pair.r1[end + offset] = b'G';
            }
        }
        assert!(
            !contains(&pair.r1, ADAPTER),
            "adapter is no longer matchable"
        );

        let merged = merge(&pair, &OverlapCfg::default()).unwrap();
        assert_eq!(
            merged.seq, pair.molecule,
            "an unmatchable adapter is still trimmed, because geometry decided it"
        );
    }

    #[test]
    fn short_r1_with_r2_reading_past_the_barcode() {
        // R1 100 bp, R2 150 bp reading 10 bp past the CB into the R1 primer.
        // Solving R2 = core | A(30) | rc(UMI)(12) | rc(CB)(16) | rc(primer)(10)
        // gives an 82 bp cDNA and a 140 bp molecule. R1 covers molecule [0,100)
        // and rc(R2) covers [-10,140): R1 is strictly interior, overhang on both
        // sides, and R1 never reaches the adapter.
        let core = complex(82, 28);
        let pair = build(&core, 100, 150);
        assert_eq!(pair.molecule.len(), 140);

        let merged = merge(&pair, &OverlapCfg::default()).expect("interior R1 merges");

        assert_eq!(merged.seq, pair.molecule);
        // Every base of R1 is inside the overlap: CB, UMI, dT and 42 bp of cDNA
        // all come out duplex-confirmed.
        assert_eq!(merged.overlap_len, 100);
        // The 10 bp of primer is dropped; R2 alone extends the molecule 40 bp.
        assert_eq!(merged.left_overhang, 0);
        assert_eq!(merged.right_overhang, 40);
        // R1 never saw adapter, but R2 read through — so this is a short insert.
        // A rule that only checks R1 would wrongly call it Medium.
        assert_eq!(merged.class, InsertClass::Short);
        assert_eq!(merged.yc.as_ref().unwrap().to_string(), "+100-40");

        // What Phase 3 will emit: strip CB+UMI+dT, keep the 82 bp cDNA, and
        // carry the tag across the trim and the flip into R2 orientation.
        let cdna_start = dt_run_end(&merged.seq, BARCODE_LEN);
        assert_eq!(cdna_start, 58);
        let mut tag = merged.yc.clone().unwrap();
        tag.trim(cdna_start as u64, 0);
        assert_eq!(tag.to_string(), "+42-40");
        tag.reverse_complement();
        // In sense orientation the 5'-most 40 bp are R2-only, the rest duplex.
        assert_eq!(tag.to_string(), "40-42+");
        assert_eq!(tag.sequence_length(), 82);

        let emitted = crate::anchors::revcomp(&merged.seq[cdna_start..]);
        assert_eq!(emitted, core, "the emitted read is the cDNA, sense strand");
    }

    #[test]
    fn merges_across_four_homopolymer_runs() {
        let core = core_with_homopolymers();
        let pair = build(&core, 150, 150);
        let merged = merge(&pair, &OverlapCfg::default()).expect("homopolymers must not block");

        assert_eq!(merged.seq, pair.molecule);
        assert_eq!(merged.discordant, 0);

        // The tag is lossless with respect to the merged construct: it
        // reconstructs each mate's contribution, not the read-through the merge
        // discarded. R1 contributes a prefix, rc(R2) a suffix.
        let tag = merged.yc.as_ref().unwrap();
        let decoded = decode(tag, &merged.seq, &merged.qual).unwrap();
        assert!(
            pair.r1.starts_with(&decoded.r1_seq),
            "R1 contributes a prefix"
        );
        assert!(
            crate::anchors::revcomp(&pair.r2).ends_with(&decoded.r2_seq),
            "rc(R2) contributes a suffix"
        );
    }

    #[test]
    fn the_dt_block_does_not_starve_the_seed_probe() {
        // The mates meet across the dT run, so the first usable seed on rc(R2)
        // is ~20 bp in. A fixed positional probe window would see only T.
        let core = complex(130, 24);
        let pair = build(&core, 150, 150);
        let rc2 = crate::anchors::revcomp(&pair.r2);
        assert!(
            rc2[..20].iter().all(|&b| b == b'T'),
            "fixture must start rc(R2) inside the dT block"
        );
        assert!(merge(&pair, &OverlapCfg::default()).is_some());
    }

    #[test]
    fn the_merged_read_still_yields_cb_umi_and_cdna() {
        // What Phase 3 does with the merge: barcode at a fixed offset, cDNA
        // after the dT run — now backed by both strands.
        let core = complex(60, 25);
        let mut pair = build(&core, 150, 150);
        // Q20 inputs, so the summed duplex quality is visible below the Q40 cap.
        pair.r1q = vec![b'5'; pair.r1.len()];
        pair.r2q = vec![b'5'; pair.r2.len()];
        let merged = merge(&pair, &OverlapCfg::default()).unwrap();

        assert_eq!(&merged.seq[..16], CB);
        assert_eq!(&merged.seq[16..28], UMI);
        let cdna_start = dt_run_end(&merged.seq, BARCODE_LEN);
        assert_eq!(cdna_start, BARCODE_LEN + DT);
        assert_eq!(
            &merged.seq[cdna_start..],
            &crate::anchors::revcomp(&core)[..]
        );

        // The barcode sat inside the overlap, so both strands confirm it: Q20 +
        // Q20 saturates at the Q40 cap. On a short fragment, merging before
        // trimming buys duplex confirmation of the CB and UMI themselves.
        assert!(merged.qual[..28].iter().all(|&q| q == b'I'));
    }

    #[test]
    fn tolerates_scattered_sequencing_mismatches() {
        let core = complex(130, 14);
        for every in [40usize, 20, 12] {
            let mut pair = build(&core, 150, 150);
            let mut corrupted = 0;
            for i in (30..pair.r2.len()).step_by(every) {
                pair.r2[i] = if pair.r2[i] == b'A' { b'C' } else { b'A' };
                corrupted += 1;
            }
            let merged = merge(&pair, &OverlapCfg::default())
                .unwrap_or_else(|| panic!("every {every} bp ({corrupted} errors) should merge"));
            assert_eq!(merged.seq.len(), pair.molecule.len());
            // Inside the overlap R1 is the higher-quality call and wins every
            // disagreement, so the molecule is recovered exactly. Beyond it R2 is
            // the only evidence there is and its errors are copied through —
            // which is precisely why merging is worth doing.
            let duplex_end = merged.left_overhang + merged.overlap_len;
            assert_eq!(
                merged.seq[..duplex_end],
                pair.molecule[..duplex_end],
                "every {every} bp: the overlap must be error-corrected"
            );
            assert!(merged.discordant > 0);
        }
    }

    #[test]
    fn polya_tail_does_not_seed_a_false_merge() {
        // Far too long to overlap, but both mates carry long A/T runs.
        // Homopolymer seed masking is what keeps this from merging.
        let mut core = complex(200, 15);
        core.extend_from_slice(b"AAAAAAAAAAAAAAAAAAAA");
        core.extend_from_slice(&complex(200, 16));
        *core.last_mut().unwrap() = b'C';
        let pair = build(&core, 150, 150);

        assert!(!pair.overlaps);
        assert!(
            merge(&pair, &OverlapCfg::default()).is_none(),
            "poly-A/T must not be enough to anchor a merge"
        );
    }

    #[test]
    fn an_overlap_that_is_only_homopolymer_does_not_merge() {
        // The mates meet, but only across a run of A. No complex seed exists in
        // the shared region, so we decline rather than merge on evidence that
        // means nothing. Built directly: the fixed 58 bp CB+UMI+dT offset cannot
        // place an overlap wholly inside an internal homopolymer.
        let mut r1 = complex(40, 21);
        r1.extend(std::iter::repeat_n(b'A', 40));
        let mut rc2 = vec![b'A'; 40];
        rc2.extend_from_slice(&complex(40, 22));
        let r2 = crate::anchors::revcomp(&rc2);

        assert!(overlap(
            &r1,
            &vec![b'I'; r1.len()],
            &r2,
            &vec![b'I'; r2.len()],
            &OverlapCfg::default()
        )
        .is_none());
    }

    #[test]
    fn indel_between_mates_is_characterised() {
        // The merge is gapless by construction. An indel inside the overlap
        // frame-shifts everything after it, so the mismatch budget rejects the
        // alignment: we decline rather than emit a shifted chimera. The dT block
        // makes this worth watching — homopolymers are where Illumina indels live.
        let core = complex(130, 19);
        let mut pair = build(&core, 150, 150);
        pair.r2.insert(60, b'A');
        pair.r2.truncate(150);
        pair.r2q = vec![b'I'; pair.r2.len()];

        match merge(&pair, &OverlapCfg::default()) {
            None => {} // preferred outcome
            Some(m) => assert!(
                m.seq.starts_with(&pair.molecule[..40]),
                "a partial merge must not corrupt the R1 side"
            ),
        }
    }

    #[test]
    fn indel_outside_the_overlap_is_harmless() {
        let core = complex(130, 20);
        let mut pair = build(&core, 150, 150);
        pair.r2.insert(148, b'A');
        pair.r2.truncate(150);
        pair.r2q = vec![b'I'; pair.r2.len()];

        let merged = merge(&pair, &OverlapCfg::default()).expect("overlap is untouched");
        assert_eq!(merged.seq, pair.molecule);
    }
}
