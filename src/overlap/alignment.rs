use crate::anchors::{comp, hamming_leq};
use crate::yc::{self, EncodeStats};
use super::types::{
    build_index, seed_positions, usable_kmer, InsertClass, Merged, OverlapCfg, Scratch, QUAL_OFFSET,
};

/// Merge `r1` with `rc(r2)` anchored so that `r1[loc1]` aligns to `rc2[loc2]`,
/// or `None` if the implied overlap is too short or too divergent.
#[allow(clippy::too_many_arguments)]
pub fn check_overlap(
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
