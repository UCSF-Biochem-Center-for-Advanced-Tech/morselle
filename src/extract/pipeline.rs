use super::types::{
    ExtractCfg, Extracted, FLEX_R2_END_ANCHOR, FLEX_R2_END_ANCHOR_MAX_MM,
    FLEX_R2_END_ANCHOR_MAX_SEARCH, FLEX_R2_END_ANCHOR_PREFIX_LEN, MIN_PARTIAL_ADAPTER,
    MIN_PARTIAL_TSO, TSO_RC,
};
use crate::anchors::{hamming_leq, ILLUMINA_3P_ADAPTER};

/// Where the reverse-complemented TSO starts, if this read carries one at the cDNA's 3' end.
fn find_tso(seq: &[u8], cdna_start: usize, cdna_end: usize) -> Option<usize> {
    let m = TSO_RC.len();
    // Full-length TSO immediately before the adapter.
    if cdna_end >= cdna_start + m {
        let start = cdna_end - m;
        if hamming_leq(&seq[start..cdna_end], TSO_RC, (m / 10) as u32) {
            return Some(start);
        }
    }
    // Partial TSO running off the 3' end: the read holds a prefix of it. Longest
    // match wins, so a read that stopped mid-TSO still trims everything it saw.
    if cdna_end == seq.len() {
        let longest = (cdna_end - cdna_start).min(m - 1);
        for len in (MIN_PARTIAL_TSO..=longest).rev() {
            let start = cdna_end - len;
            if hamming_leq(&seq[start..cdna_end], &TSO_RC[..len], (len / 10) as u32) {
                return Some(start);
            }
        }
    }
    None
}

/// Extend the dT-primer run rightward from `approx_start` (tolerating isolated mismatches).
pub fn dt_run_end(seq: &[u8], approx_start: usize) -> usize {
    let n = seq.len();
    // Probe for the run start at approx_start, then ±1 (boundary fuzz).
    let start = if approx_start < n && seq[approx_start] == b'T' {
        approx_start
    } else if approx_start > 0 && approx_start - 1 < n && seq[approx_start - 1] == b'T' {
        approx_start - 1
    } else if approx_start + 1 < n && seq[approx_start + 1] == b'T' {
        approx_start + 1
    } else {
        return approx_start.min(n);
    };
    let mut last = start;
    let mut i = start + 1;
    while i < n {
        if seq[i] == b'T' {
            last = i;
        } else if seq[i - 1] != b'T' {
            break; // two consecutive non-T -> run ends
        }
        i += 1;
    }
    last + 1
}

/// Leftmost position ≥ `from` where the 3' read-through adapter matches. Catches:
///   * full-length match within `max_mm` substitutions (the common case);
///   * full-length match with a single indel + ≤`max_mm` substitutions — SBX makes
///     indel errors in homopolymers (e.g. AAGAGC→AAAGAGC), which a pure-Hamming
///     seed misses because the insertion frame-shifts it to many mismatches;
///   * a partial adapter (≥ `MIN_PARTIAL_ADAPTER` bp) running off the 3' end.
/// Trimming the adapter ourselves keeps cDNA out of downstream tools' adapter
/// trimmers (10x reference mis-trims long reads that still carry the adapter).
/// Allocation-free. The indel check runs only at a cheap 3 bp prefix anchor.
pub fn find_adapter(seq: &[u8], from: usize, max_mm: u32) -> Option<usize> {
    let a = ILLUMINA_3P_ADAPTER;
    let m = a.len();
    let end = seq.len();
    let mut i = from;
    while i < end {
        let remain = end - i;
        if remain >= m {
            if hamming_leq(&seq[i..i + m], a, max_mm) {
                return Some(i);
            }
            // Indel path: anchor on the indel-safe 3 bp prefix (no homopolymer),
            // then verify a single-indel alignment.
            if seq[i] == a[0]
                && seq[i + 1] == a[1]
                && seq[i + 2] == a[2]
                && adapter_single_indel(&seq[i..], a, max_mm)
            {
                return Some(i);
            }
        } else if remain >= MIN_PARTIAL_ADAPTER {
            // Partial adapter running off the 3' end: match the adapter prefix,
            // mismatch budget scaling with overlap (7 bp exact; 8+ bp allow 1) —
            // false trims stay ~0.04% and only cost a few 3' bp, while leaving an
            // adapter risks crashing downstream trimmers.
            let allowed = (remain / 8) as u32;
            if hamming_leq(&seq[i..end], &a[..remain], allowed) {
                return Some(i);
            }
        } else {
            break;
        }
        i += 1;
    }
    None
}

/// True if adapter `a` matches `s` (which starts at the candidate position) with
/// exactly one indel and ≤`max_mm` substitutions. Tries an insertion in the read
/// (extra base, e.g. homopolymer expansion) and a deletion in the read.
fn adapter_single_indel(s: &[u8], a: &[u8], max_mm: u32) -> bool {
    let m = a.len();
    // Insertion in read: skip one read base. s[0..k]≈a[0..k], s[k+1..m+1]≈a[k..m].
    if s.len() >= m + 1 {
        for k in 0..m {
            let mut mm = 0u32;
            for j in 0..k {
                if s[j] != a[j] {
                    mm += 1;
                }
            }
            for j in k..m {
                if s[j + 1] != a[j] {
                    mm += 1;
                }
            }
            if mm <= max_mm {
                return true;
            }
        }
    }
    // Deletion in read: skip one adapter base. s[0..k]≈a[0..k], s[k..m-1]≈a[k+1..m].
    if s.len() >= m - 1 {
        for k in 0..m {
            let mut mm = 0u32;
            for j in 0..k {
                if s[j] != a[j] {
                    mm += 1;
                }
            }
            for j in (k + 1)..m {
                if s[j - 1] != a[j] {
                    mm += 1;
                }
            }
            if mm <= max_mm {
                return true;
            }
        }
    }
    false
}

/// All single-position deletions of `seq` (one per position).
fn single_deletion_variants(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    (0..seq.len()).map(move |i| {
        let mut v = Vec::with_capacity(seq.len() - 1);
        v.extend_from_slice(&seq[..i]);
        v.extend_from_slice(&seq[i + 1..]);
        v
    })
}

/// All single-position insertions into `seq` (every position, every
/// non-redundant base).
fn single_insertion_variants(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
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

/// Leftmost position ≥ `from` where a confident (possibly single-indel-
/// shifted) prefix of [`FLEX_R2_END_ANCHOR`] matches -- the real far
/// boundary of R2 on a Flex read (end of `Ligated_Probe_insert`).
///
/// Tries the exact-length (substitution-only) window at every offset first,
/// same as originally; a real, single-base indel landing *before* the
/// anchor is already absorbed by sliding `from` forward. But a real indel
/// landing *inside* the 15bp comparison span itself misaligns every base
/// after it, so no substitution-only offset can find it regardless of where
/// the window starts -- real SBX errors are dominated by single-base
/// indels, not substitutions (see `probe_set.rs`'s
/// `single_deletion_variants`/`single_insertion_variants` doc comments and
/// the Flex design notes' real error-profile finding: 85.2% deletion, 9.2%
/// insertion of recoverable Hamming-miss shapes). Falls back to trying each
/// single-deletion/-insertion variant of the anchor prefix at every offset
/// too, mirroring that same fix. Real measured cost of the original,
/// substitution-only-only version on real CAT06 data: hard-gating on it
/// dropped `class_full` from 66.07% to 56.39% on the same raw shards
/// (23,185,519 fewer real, probe-matching reads) with *no* improvement in
/// hit-rate-among-survivors (96.4% -> 96.5%, unchanged) -- the signature of
/// a real, fixable recall gap in anchor detection, not a genuine precision
/// gain from being stricter.
///
/// Bounded to [`FLEX_R2_END_ANCHOR_MAX_SEARCH`] bp past `from` -- see its
/// doc comment; this both caps the wasted work on non-GEX reads that will
/// never match and avoids ever reporting a spurious, structurally-
/// impossible match far into the read.
fn find_flex_r2_end(seq: &[u8], from: usize) -> Option<usize> {
    let a = &FLEX_R2_END_ANCHOR[..FLEX_R2_END_ANCHOR_PREFIX_LEN];
    let m = a.len();
    let end = seq.len().min(from + FLEX_R2_END_ANCHOR_MAX_SEARCH);
    if end < from + m {
        return None;
    }
    if let Some(i) =
        (from..=end - m).find(|&i| hamming_leq(&seq[i..i + m], a, FLEX_R2_END_ANCHOR_MAX_MM))
    {
        return Some(i);
    }
    if end >= from + (m - 1) {
        for variant in single_deletion_variants(a) {
            let w = variant.len();
            if let Some(i) = (from..=end - w)
                .find(|&i| hamming_leq(&seq[i..i + w], &variant, FLEX_R2_END_ANCHOR_MAX_MM))
            {
                return Some(i);
            }
        }
    }
    if end >= from + (m + 1) {
        for variant in single_insertion_variants(a) {
            let w = variant.len();
            if let Some(i) = (from..=end - w)
                .find(|&i| hamming_leq(&seq[i..i + w], &variant, FLEX_R2_END_ANCHOR_MAX_MM))
            {
                return Some(i);
            }
        }
    }
    None
}

/// Bracket R1 (CB+UMI) and R2 (cDNA) on a canonical read.
///
/// `cb_umi` and `polyt_start` come from [`crate::orient::Oriented`].
///
/// `likely_non_gex`: routing hint for [`ExtractCfg::flex_r2_boundary`] --
/// true when the *caller* already has good reason (typically an early,
/// i5-only sample-sheet lookup -- see `lib.rs::process_one`) to believe this
/// read is NOT from a Flex GEX library, even though Flex mode is on for the
/// whole run (raw shards mix every library together). Ignored when
/// `cfg.flex_r2_boundary` is false. Always `false` is a safe default --
/// `--flex` assumes every read is Flex GEX unless told otherwise.
pub fn extract(
    seq: &[u8],
    cb_umi: Option<(usize, usize)>,
    polyt_start: Option<usize>,
    cfg: &ExtractCfg,
    likely_non_gex: bool,
) -> Extracted {
    let mut polyt_end = None;
    let mut polyt_run_len = 0u16;
    let mut r2 = None;
    let mut adapter_pos = None;
    let mut adapter_trimmed = false;
    let mut tso_pos = None;
    let mut tso_trimmed = false;

    if let Some(pts) = polyt_start {
        let cdna_start = dt_run_end(seq, pts);
        polyt_run_len = cdna_start.saturating_sub(pts).min(u16::MAX as usize) as u16;
        polyt_end = Some(cdna_start);

        if cfg.flex_r2_boundary && !likely_non_gex {
            // Verified boundary: the real probe-insert/constant boundary
            // (FLEX_R2_END_ANCHOR). `--flex` defaults to assuming every read
            // is Flex GEX; the caller is responsible for ruling that out via
            // `likely_non_gex` (an early, cheap i5-only sample-sheet lookup)
            // for libraries that don't carry this construct at all (CRISPR/
            // Feature Barcoding sharing the same invocation) -- see
            // `lib.rs::process_one`.
            if let Some(cdna_end) = find_flex_r2_end(seq, cdna_start) {
                adapter_pos = Some(cdna_end);
                adapter_trimmed = true;
                if cdna_end > cdna_start {
                    r2 = Some((cdna_start, cdna_end));
                }
            } else {
                // Anchor not confidently located, even with its own indel
                // tolerance -- real data shows this population's
                // hit-rate-among-survivors is statistically identical to the
                // verified population (96.4-96.5% either way), so most of it
                // isn't garbage; it just can't be verified via this one
                // landmark (didn't sequence far enough to reach it, or
                // carries more errors than even the indel-tolerant search's
                // budget). Fall back to a length floor sized to structurally
                // guarantee the real probe sequence fits, rather than
                // rejecting outright -- see `flex_min_cdna_len_no_anchor`'s
                // doc comment for the real derivation. `adapter_trimmed`
                // stays false here: this boundary is a length floor, not a
                // verified landmark, and downstream accounting (e.g.
                // `per_sample_early_boundary_found`) should reflect that
                // honestly.
                let cdna_end = seq.len().min(cdna_start + FLEX_R2_END_ANCHOR_MAX_SEARCH);
                if cdna_end >= cdna_start + cfg.flex_min_cdna_len_no_anchor {
                    r2 = Some((cdna_start, cdna_end));
                }
            }
        } else {
            // Standard path: non-Flex mode, or Flex mode but the caller
            // already ruled out this read carrying FLEX_R2_END_ANCHOR at
            // all -- skip straight to the generic-adapter/length-floor
            // logic rather than wasting a search that was never going to
            // match.
            adapter_pos = find_adapter(seq, cdna_start, cfg.adapter_max_mm);
            let mut cdna_end = adapter_pos.unwrap_or(seq.len());
            adapter_trimmed = adapter_pos.is_some();
            // The TSO sits between the cDNA and the adapter on any fragment that
            // reached the cDNA's 5' end. Trim before the length gate, so a fragment
            // that is *only* TSO is rejected rather than emitted as cDNA.
            if cfg.trim_tso && cdna_end > cdna_start {
                if let Some(pos) = find_tso(seq, cdna_start, cdna_end) {
                    tso_pos = Some(pos);
                    tso_trimmed = true;
                    cdna_end = pos;
                }
            }
            // Retain the WHOLE adapter-trimmed cDNA fragment here; any `--max-r2-len`
            // cap is applied at write time, after the fragment is reverse-complemented
            // into R2 orientation, so the cap keeps the bases the R2 primer would
            // actually read (see the emitter in lib.rs). See `ExtractCfg::max_r2_len`.
            if cdna_end > cdna_start && cdna_end - cdna_start >= cfg.min_cdna_len {
                r2 = Some((cdna_start, cdna_end));
            }
        }
    }

    Extracted {
        r1: cb_umi,
        r2,
        polyt_end,
        polyt_run_len,
        adapter_pos,
        adapter_trimmed,
        tso_pos,
        tso_trimmed,
    }
}
