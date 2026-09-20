//! R1/R2 bracket extraction and adapter/homopolymer trimming (plan §4 step 4, M2).
//!
//! Operates on the canonical (P5→P7) read produced by M1, reusing the coordinates
//! it already found:
//! - **R1** = the 28 bp CB+UMI block located by [`crate::orient`].
//! - **R2** = cDNA between the dT-primer run end and the 3' read-through adapter
//!   (`AGATCGGAAGAGC`), found by a cheap Hamming seed; if the read ends mid-cDNA
//!   (no adapter), R2 runs to the read end. **Standard 3' capture only** --
//!   see [`ExtractCfg::flex_r2_boundary`] for 10x Flex, which bounds R2 by a
//!   real structural landmark ([`FLEX_R2_END_ANCHOR`]) instead.
//!
//! All positions are returned as ranges into the canonical sequence — zero-copy;
//! the caller slices `seq`/`qual`.

use crate::anchors::{hamming_leq, ILLUMINA_3P_ADAPTER};

/// Reverse-complement of [`crate::anchors::TSO`], as the template-switch oligo
/// appears in the canonical frame: between the cDNA and the 3' read-through
/// adapter. A fragment carries it only when it reached the cDNA's 5' end, so it
/// is enriched in short fragments. Verified against `revcomp(TSO)` in tests.
const TSO_RC: &[u8] = b"CCCATGTACTCTGCGTTGATACCACTGCTT";

/// Shortest partial TSO trimmed at the 3' end of a read that never reached an
/// adapter. 12 bp with a 1-mismatch budget keeps false trims near zero.
const MIN_PARTIAL_TSO: usize = 12;

/// 10x Flex GEX/CRISPR construct: the constant sequence immediately following
/// `Ligated_Probe_insert` (the 50bp gene/guide probe), in canonical (pre-RC)
/// coordinates -- i.e. this marks the FAR end of the cDNA/R2 region for Flex
/// reads, the same role `ILLUMINA_3P_ADAPTER` plays for standard 3' capture.
/// Only a confident prefix is searched, not the full 33bp -- the boundary IS
/// the prefix's start position, so a short high-confidence match locates it
/// exactly as precisely as a full-length one while being far more tolerant
/// of per-base sequencing error accumulating over a long required span.
/// Empirically compared on real CAT06 Flex data (25-shard diagnostic,
/// 2026-09-17): a 10-15bp prefix at max_mm 1-2 is found in ~82% of reads,
/// vs 74.4% requiring the full 33bp at max_mm 2 -- and the resulting
/// immediate-downstream probe match rate is identical either way (~48.5%),
/// confirming the short prefix locates the same true boundary just as
/// reliably. 15bp/max_mm 2 measured marginally best of the lengths tried.
const FLEX_R2_END_ANCHOR: &[u8] = b"TGGAATTCTCGGGTGCCAAGGAACTCCAGTCAC";
const FLEX_R2_END_ANCHOR_PREFIX_LEN: usize = 15;
const FLEX_R2_END_ANCHOR_MAX_MM: u32 = 2;

/// Max bp past `cdna_start` searched for [`FLEX_R2_END_ANCHOR`]. The real
/// molecule structure fixes this distance: `BC10`(10) + a 29bp constant +
/// `Ligated_Probe_insert`(50) = 89bp past the probe barcode, which itself
/// starts ~10bp past `cdna_start` (`probe_barcode::DEFAULT_PROBE_OFFSET`,
/// with its own `probe_barcode_window` of jitter) -- ~109bp worst case, so
/// 150bp leaves real margin. This is every read's FIRST boundary check
/// (library identity isn't known until demux, which itself runs *after*
/// this), so an unbounded scan would waste a full-read pass on every
/// CRISPR/Feature-Barcoding read (SBX reads run to 1000+bp) before ever
/// falling back to the generic adapter search they actually need.
const FLEX_R2_END_ANCHOR_MAX_SEARCH: usize = 150;

/// Tuning for [`extract`].
#[derive(Debug, Clone, Copy)]
pub struct ExtractCfg {
    /// Minimum cDNA length to emit an R2.
    pub min_cdna_len: usize,
    /// Max mismatches for the 3' adapter trim seed.
    pub adapter_max_mm: u32,
    /// Trim the template-switch oligo when it sits between the cDNA and the 3'
    /// adapter. Off by default: enabling it changes emitted cDNA.
    pub trim_tso: bool,
    /// Cap emitted R2 length; 0 = uncapped. SBX cDNA can be >1 kb, but 10x reference
    /// hard-caps reads at 649 bp. NOT applied here — `extract` retains the whole
    /// fragment; the cap is applied by the emitter after reverse-complementing into
    /// R2 orientation, so it keeps the bases the R2 primer would read.
    pub max_r2_len: usize,
    /// 10x Flex mode: require [`FLEX_R2_END_ANCHOR`] to bound R2's far end,
    /// instead of the generic `ILLUMINA_3P_ADAPTER` search + `min_cdna_len`
    /// floor. Flex GEX's R2 always contains a fixed-length probe insert that
    /// a plain length floor can't guarantee is actually present -- finding
    /// the real boundary is a tighter cut (drops the uninformative tail past
    /// the probe insert that the generic adapter search usually never
    /// reaches -- real data: only 13.5% of reads hit it) that lands the
    /// probe insert at a fixed, predictable position in the emitted R2.
    ///
    /// Hard gate by default (`--flex` assumes every read is Flex GEX; no
    /// anchor found means no complete probe insert was captured, so R2 is
    /// correctly rejected). The one exception is `extract()`'s
    /// `likely_non_gex` parameter: a single invocation processes every
    /// library in a raw shard together (GEX, CRISPR/Feature Barcoding, ...),
    /// demuxed only downstream of this call, so this flag being on doesn't
    /// mean every read is Flex GEX -- callers rule that out per-read via a
    /// cheap early sample-sheet lookup (see `lib.rs::process_one`) rather
    /// than this flag falling back silently for everyone (confirmed on real
    /// data: a blanket fallback and a blanket hard gate are both wrong --
    /// the former lands the probe insert at a jittery offset for GEX reads,
    /// the latter drops virtually all non-GEX reads sharing the run). Leave
    /// `false` for non-Flex (standard 3' capture, real free cDNA, no such
    /// landmark to find).
    pub flex_r2_boundary: bool,
    /// When `flex_r2_boundary` is on and [`FLEX_R2_END_ANCHOR`] can't be
    /// confidently located (even with its own indel tolerance), fall back to
    /// accepting R2 up to [`FLEX_R2_END_ANCHOR_MAX_SEARCH`] anyway *if* that
    /// span is at least this long -- rather than rejecting the read
    /// outright. Real, measured rationale: on real CAT06 data, requiring the
    /// literal anchor discarded a large population of reads whose
    /// hit-rate-among-survivors was statistically identical to the reads
    /// that *did* verify (96.4-96.5% either way) -- i.e. most of what the
    /// hard gate excludes isn't garbage, it's real, otherwise-valid reads
    /// that simply can't be verified via this one specific landmark (didn't
    /// sequence far enough to reach it, or carry more errors than even the
    /// indel-tolerant search's budget allows). A length floor sized to
    /// *structurally guarantee the real probe sequence fits* gives these
    /// reads a fair shot at `probeset-align`/`probe_barcode` matching
    /// instead of being dropped before ever trying. Default derivation is
    /// real, not guessed: probe length is 50bp (confirmed directly from
    /// `Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv`)
    /// plus `ProbeMatcher`'s own real default search jitter (`--window 15`)
    /// = 65bp worst case.
    pub flex_min_cdna_len_no_anchor: usize,
}

/// Extraction result: ranges `[start, end)` into the canonical sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Extracted {
    /// CB+UMI range (R1), 28 bp, or `None` if not located.
    pub r1: Option<(usize, usize)>,
    /// cDNA range (R2) in the canonical (polyT) frame, or `None` if absent / below
    /// `min_cdna_len`. NOTE: this slice is on the same strand as R1; the emitter
    /// reverse-complements it (and reverses its quality) so R2 follows the Illumina
    /// paired-end convention and maps sense to the gene.
    pub r2: Option<(usize, usize)>,
    /// End of the dT-primer run = cDNA start.
    pub polyt_end: Option<usize>,
    /// Length of the dT-primer run actually observed (QC).
    pub polyt_run_len: u16,
    /// Position of the 3' read-through adapter (= start of `rc(R2 primer)`), when
    /// found. Marks the cDNA end and lets demux locate the downstream i7.
    pub adapter_pos: Option<usize>,
    /// The 3' read-through adapter was found and trimmed.
    pub adapter_trimmed: bool,
    /// Where the reverse-complemented TSO starts, when found and trimmed.
    pub tso_pos: Option<usize>,
    /// The TSO was found abutting the cDNA's 3' end and trimmed.
    pub tso_trimmed: bool,
}

/// Where the reverse-complemented TSO starts, if this read carries one at the
/// cDNA's 3' end.
///
/// **Anchored, never scanned.** The TSO can only sit immediately 5' of the
/// read-through adapter, or run off the 3' end of a read that never reached one.
/// Searching for it inside the cDNA would let a chance match truncate real
/// sequence — and at 30 bp a partial match is not rare. The mismatch budget
/// scales with the matched length (10%), mirroring [`find_adapter`].
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

/// Extend the dT-primer run rightward from `approx_start` (tolerating isolated
/// mismatches) and return the cDNA start = one past the run's last T. If no T is
/// near `approx_start`, returns `approx_start` (run length 0 — polyT absent).
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

/// Minimum 3' read-through adapter overlap to trim a partial adapter at the read
/// end. Short partials must match near-exactly (mismatch budget scales with
/// length, below), so 7 bp is safe — chance of a 7 bp exact adapter prefix at a
/// read's 3' end is ~6e-5.
const MIN_PARTIAL_ADAPTER: usize = 7;

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
        b"ACGT".iter().filter(move |&&b| b != existing).map(move |&b| {
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
/// FLEX_SUPPORT_PLAN.md's real error-profile finding: 85.2% deletion, 9.2%
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
    if let Some(i) = (from..=end - m).find(|&i| hamming_leq(&seq[i..i + m], a, FLEX_R2_END_ANCHOR_MAX_MM)) {
        return Some(i);
    }
    if end >= from + (m - 1) {
        for variant in single_deletion_variants(a) {
            let w = variant.len();
            if let Some(i) = (from..=end - w).find(|&i| hamming_leq(&seq[i..i + w], &variant, FLEX_R2_END_ANCHOR_MAX_MM)) {
                return Some(i);
            }
        }
    }
    if end >= from + (m + 1) {
        for variant in single_insertion_variants(a) {
            let w = variant.len();
            if let Some(i) = (from..=end - w).find(|&i| hamming_leq(&seq[i..i + w], &variant, FLEX_R2_END_ANCHOR_MAX_MM)) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchors::{revcomp, P5, R1_PRIMER};
    use crate::orient::{orient, FallbackMode, OrientCfg};
    use crate::pipeline::Orientation;

    #[test]
    fn adapter_exact_and_substitution() {
        let cdna = b"ACGTACGTACGTACGTACGT"; // 20bp, no adapter
        let mut s = cdna.to_vec();
        let astart = s.len();
        s.extend_from_slice(ILLUMINA_3P_ADAPTER); // exact adapter
        s.extend_from_slice(b"CACGTCTGAACTCC");
        assert_eq!(find_adapter(&s, 0, 2), Some(astart));
        // one substitution in the adapter still found
        let mut s2 = cdna.to_vec();
        s2.extend_from_slice(b"AGATCGGAAGTGC"); // A->T at pos 10
        assert_eq!(find_adapter(&s2, 0, 2), Some(cdna.len()));
    }

    #[test]
    fn adapter_homopolymer_insertion() {
        // The real compatibility-breaking case: AAGAGC -> AAAGAGC (extra A).
        let cdna = b"ACGTACGTACGTACGTACGT";
        let mut s = cdna.to_vec();
        let astart = s.len();
        s.extend_from_slice(b"AGATCGGAAAGAGCACACGTCTGAACTCC"); // indel adapter
                                                               // pure-Hamming (no indel handling) would miss this; find_adapter catches it
        assert_eq!(find_adapter(&s, 0, 2), Some(astart));
    }

    #[test]
    fn adapter_partial_at_end() {
        let cdna = b"ACGTACGTACGTACGTACGT";
        let mut s = cdna.to_vec();
        let astart = s.len();
        s.extend_from_slice(&ILLUMINA_3P_ADAPTER[..11]); // 11bp partial (>= MIN_PARTIAL)
        assert_eq!(find_adapter(&s, 0, 2), Some(astart));
        // a 7bp exact partial (the dominant real residual) IS called
        let mut s7 = cdna.to_vec();
        s7.extend_from_slice(&ILLUMINA_3P_ADAPTER[..7]); // "AGATCGG"
        assert_eq!(find_adapter(&s7, 0, 2), Some(cdna.len()));
        // a too-short tail (< MIN_PARTIAL) is not called
        let mut s2 = cdna.to_vec();
        s2.extend_from_slice(&ILLUMINA_3P_ADAPTER[..5]);
        assert_eq!(find_adapter(&s2, 0, 2), None);
    }

    #[test]
    fn adapter_absent_not_falsely_trimmed() {
        let s = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        assert_eq!(find_adapter(s, 0, 2), None);
    }

    /// Build a canonical read: CB+UMI | dT | rc(cDNA) | [rc(TSO)] | [adapter].
    fn canonical(cdna: &[u8], with_tso: bool, with_adapter: bool) -> Vec<u8> {
        // A cDNA ending in A makes revcomp start with T, which the dT-run walker
        // legitimately absorbs — a fixture trap, not a bug. Refuse to build one.
        assert_ne!(*cdna.last().unwrap(), b'A', "cDNA must not end in A");
        let mut v = vec![b'A'; 28];
        v.extend(std::iter::repeat_n(b'T', 30));
        v.extend_from_slice(&revcomp(cdna));
        if with_tso {
            v.extend_from_slice(TSO_RC);
        }
        if with_adapter {
            v.extend_from_slice(ILLUMINA_3P_ADAPTER);
        }
        v
    }

    fn tso_cfg() -> ExtractCfg {
        ExtractCfg {
            trim_tso: true,
            ..ecfg()
        }
    }

    #[test]
    fn tso_rc_matches_the_declared_oligo() {
        assert_eq!(TSO_RC, revcomp(crate::anchors::TSO).as_slice());
    }

    #[test]
    fn trims_a_full_tso_between_cdna_and_adapter() {
        let cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCG";
        let seq = canonical(cdna, true, true);
        let e = extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false);
        assert!(e.tso_trimmed);
        let (s, t) = e.r2.unwrap();
        assert_eq!(&seq[s..t], revcomp(cdna).as_slice());
        assert_eq!(e.tso_pos, Some(t));
    }

    #[test]
    fn leaves_the_tso_alone_by_default() {
        let cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCG";
        let seq = canonical(cdna, true, true);
        let e = extract(&seq, Some((0, 28)), Some(28), &ecfg(), false);
        assert!(!e.tso_trimmed);
        let (s, t) = e.r2.unwrap();
        // The TSO is still sitting on the end of the "cDNA".
        assert_eq!(t - s, cdna.len() + TSO_RC.len());
    }

    #[test]
    fn trims_a_partial_tso_running_off_the_read_end() {
        let cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCG";
        let mut seq = canonical(cdna, false, false);
        seq.extend_from_slice(&TSO_RC[..18]); // read stopped mid-TSO
        let e = extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false);
        assert!(e.tso_trimmed);
        let (s, t) = e.r2.unwrap();
        assert_eq!(&seq[s..t], revcomp(cdna).as_slice());
    }

    #[test]
    fn tolerates_mismatches_in_the_tso() {
        let cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCG";
        let mut seq = canonical(cdna, true, true);
        let tso_at = seq.len() - ILLUMINA_3P_ADAPTER.len() - TSO_RC.len();
        // 3 mismatches in 30 bp is inside the 10% budget.
        for offset in [4usize, 11, 19] {
            seq[tso_at + offset] = if seq[tso_at + offset] == b'A' {
                b'C'
            } else {
                b'A'
            };
        }
        assert!(extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false).tso_trimmed);

        // 6 is not.
        let mut seq = canonical(cdna, true, true);
        for offset in [2usize, 6, 10, 14, 18, 22] {
            seq[tso_at + offset] = if seq[tso_at + offset] == b'A' {
                b'C'
            } else {
                b'A'
            };
        }
        assert!(!extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false).tso_trimmed);
    }

    #[test]
    fn does_not_trim_a_tso_buried_inside_the_cdna() {
        // The oligo appears mid-cDNA, not abutting the adapter. Anchoring is what
        // stops a chance match from truncating real sequence.
        let mut cdna = b"GCTAGCTAGCATCGATCGTAGCTAGCATCG".to_vec();
        cdna.extend_from_slice(&revcomp(TSO_RC));
        cdna.extend_from_slice(b"ATCGATCGATGCTAGCTAGCATCGATCGTC");
        let seq = canonical(&cdna, false, true);
        let e = extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false);
        assert!(!e.tso_trimmed);
        let (s, t) = e.r2.unwrap();
        assert_eq!(&seq[s..t], revcomp(&cdna).as_slice());
    }

    #[test]
    fn a_fragment_that_is_only_tso_yields_no_cdna() {
        // Trimming happens before the length gate, so an adapter-dimer-like
        // fragment is rejected rather than emitted as 30 bp of "cDNA".
        let mut seq = vec![b'A'; 28];
        seq.extend(std::iter::repeat_n(b'T', 30));
        seq.extend_from_slice(TSO_RC);
        seq.extend_from_slice(ILLUMINA_3P_ADAPTER);
        let e = extract(&seq, Some((0, 28)), Some(28), &tso_cfg(), false);
        assert!(e.tso_trimmed);
        assert!(e.r2.is_none());
    }

    fn ecfg() -> ExtractCfg {
        ExtractCfg {
            min_cdna_len: 20,
            adapter_max_mm: 2,
            trim_tso: false,
            max_r2_len: 0,
            flex_r2_boundary: false,
            flex_min_cdna_len_no_anchor: 65,
        }
    }

    fn flex_ecfg() -> ExtractCfg {
        ExtractCfg {
            flex_r2_boundary: true,
            ..ecfg()
        }
    }

    fn ocfg() -> OrientCfg {
        OrientCfg {
            cb_len: 16,
            umi_len: 12,
            min_homopolymer: 12,
            r1_tail_max_mm: 2,
            anchor_max_edits: 6,
            fallback: FallbackMode::HammingTail,
            expected: None,
        }
    }

    /// Canonical 3' read with a real 3' read-through adapter after the cDNA.
    fn synthetic(cb: &[u8], umi: &[u8], cdna: &[u8], dt: usize) -> Vec<u8> {
        let mut v = Vec::new();
        v.extend_from_slice(b"AACAAT"); // leader
        v.extend_from_slice(P5);
        v.extend_from_slice(b"ACGTACGTAC"); // i5
        v.extend_from_slice(R1_PRIMER);
        v.extend_from_slice(cb);
        v.extend_from_slice(umi);
        v.extend(std::iter::repeat(b'T').take(dt)); // dT primer
        v.extend_from_slice(cdna);
        v.extend_from_slice(ILLUMINA_3P_ADAPTER);
        v.extend_from_slice(b"CACGTCTGAACTCCAGTCACTGCATGCATG"); // rest of rc(R2)+i7
        v
    }

    #[test]
    fn extracts_r1_and_trimmed_r2() {
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let cdna = b"GGCATTACTGGCATTACTGGCATTACTGGCATTACTG"; // 37 bp, no long T-run
        let seq = synthetic(cb, umi, cdna, 30);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        assert_eq!(o.orientation, Orientation::P5Led);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);

        // R1 == 28 bp and equals CB+UMI.
        let (r1s, r1e) = e.r1.unwrap();
        assert_eq!(r1e - r1s, 28);
        assert_eq!(&o.seq[r1s..r1s + 16], cb);
        assert_eq!(&o.seq[r1s + 16..r1e], umi);

        // R2 == cDNA, adapter trimmed.
        assert!(e.adapter_trimmed);
        let (r2s, r2e) = e.r2.unwrap();
        assert_eq!(&o.seq[r2s..r2e], cdna);
    }

    #[test]
    fn r2_runs_to_end_without_adapter() {
        let cb = b"AAACCCATCAAACCTG";
        let umi = b"TTTTGGGGCCCC";
        let cdna = b"GGCATTACTGGCATTACTGGCATTACTGGCATTACTG";
        // build a read that simply ends in cDNA (truncated, no adapter)
        let mut seq = Vec::new();
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(cdna);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);
        assert!(!e.adapter_trimmed);
        let (r2s, r2e) = e.r2.unwrap();
        assert_eq!(r2e, o.seq.len());
        assert_eq!(&o.seq[r2s..r2e], cdna);
    }

    /// Build a canonical Flex read: [leader][P5][i5][R1 primer][CB][UMI]
    /// [dT run][junk before the anchor, standing in for BC10+constant29]
    /// [FLEX_R2_END_ANCHOR][trailing junk, standing in for i7/P7]. Mirrors
    /// `synthetic()`, but with the Flex-specific tail instead of the
    /// standard-3'-capture `ILLUMINA_3P_ADAPTER` tail.
    fn synthetic_flex(cb: &[u8], umi: &[u8], pre_anchor: &[u8], dt: usize) -> Vec<u8> {
        let mut v = Vec::new();
        v.extend_from_slice(b"AACAAT");
        v.extend_from_slice(P5);
        v.extend_from_slice(b"ACGTACGTAC");
        v.extend_from_slice(R1_PRIMER);
        v.extend_from_slice(cb);
        v.extend_from_slice(umi);
        v.extend(std::iter::repeat(b'T').take(dt));
        v.extend_from_slice(pre_anchor);
        v.extend_from_slice(FLEX_R2_END_ANCHOR);
        v.extend_from_slice(b"ACGTACGTACACGTACGTAC"); // trailing i7/P7-ish junk
        v
    }

    #[test]
    fn flex_r2_boundary_stops_exactly_at_the_anchor() {
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        // Stand-in for BC10(10)+constant29(29): any 39bp, content doesn't matter.
        let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
        let seq = synthetic_flex(cb, umi, pre_anchor, 30);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);

        let (r2s, r2e) = e.r2.expect("flex R2 should be located");
        assert_eq!(&o.seq[r2s..r2e], pre_anchor.as_slice());
        assert!(e.adapter_trimmed);
        assert_eq!(e.adapter_pos, Some(r2e));
    }

    #[test]
    fn flex_r2_boundary_tolerates_mismatches_within_budget() {
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
        let mut anchor = FLEX_R2_END_ANCHOR.to_vec();
        // 2 mismatches inside the searched prefix (positions 3 and 9 of the
        // first FLEX_R2_END_ANCHOR_PREFIX_LEN bases) -- within max_mm budget.
        anchor[3] = if anchor[3] == b'A' { b'C' } else { b'A' };
        anchor[9] = if anchor[9] == b'A' { b'C' } else { b'A' };
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACAAT");
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(pre_anchor);
        seq.extend_from_slice(&anchor);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
        let (r2s, r2e) = e.r2.expect("anchor within mismatch budget should still be found");
        assert_eq!(&o.seq[r2s..r2e], pre_anchor.as_slice());
    }

    #[test]
    fn flex_r2_boundary_recovers_a_single_deletion_in_the_anchor() {
        // Real SBX error mode (85.2% of recoverable Hamming-miss shapes):
        // a genuine anchor with one base deleted from inside the searched
        // prefix. A pure substitution-only scan can never find this at any
        // offset, since every base after the deletion point is shifted by 1.
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
        let mut anchor = FLEX_R2_END_ANCHOR.to_vec();
        anchor.remove(7); // delete one base from inside the 15bp prefix
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACAAT");
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(pre_anchor);
        seq.extend_from_slice(&anchor);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
        let (r2s, _r2e) = e.r2.expect("single-deletion-shifted anchor should still be found");
        assert_eq!(&o.seq[r2s..r2s + pre_anchor.len()], pre_anchor.as_slice());
    }

    #[test]
    fn flex_r2_boundary_recovers_a_single_insertion_in_the_anchor() {
        // Rarer real shape (9.2% of recoverable Hamming-miss shapes) but
        // real: one extra base inserted inside the searched prefix.
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
        let mut anchor = FLEX_R2_END_ANCHOR.to_vec();
        anchor.insert(7, b'G'); // insert one extra base inside the 15bp prefix
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACAAT");
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(pre_anchor);
        seq.extend_from_slice(&anchor);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
        let (r2s, _r2e) = e.r2.expect("single-insertion-shifted anchor should still be found");
        assert_eq!(&o.seq[r2s..r2s + pre_anchor.len()], pre_anchor.as_slice());
    }

    #[test]
    fn flex_r2_boundary_rejects_when_anchor_absent_and_payload_too_short_for_a_probe() {
        // likely_non_gex=false (the default -- caller has no reason to think
        // otherwise): no FLEX_R2_END_ANCHOR anywhere, AND the payload (63bp)
        // is below flex_min_cdna_len_no_anchor's real 65bp floor (50bp probe
        // + 15bp real search jitter) -- too short to structurally contain a
        // real probe match either way, so correctly rejected even under the
        // length-floor fallback.
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let payload = b"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACAAT");
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(payload);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
        assert!(e.r2.is_none());
        assert!(!e.adapter_trimmed);
    }

    #[test]
    fn flex_r2_boundary_falls_back_to_a_length_floor_when_anchor_absent_but_payload_long_enough() {
        // Same shape as the too-short case, but 70bp of payload -- clears
        // flex_min_cdna_len_no_anchor's real 65bp floor even with no
        // FLEX_R2_END_ANCHOR anywhere. Real rationale: this read's
        // hit-rate-among-survivors on real CAT06 data is statistically the
        // same as the verified-anchor population (96.4-96.5% either way),
        // so it deserves a real shot at probeset-align/probe_barcode
        // matching instead of being dropped outright.
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let payload = b"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
        assert!(payload.len() >= 65);
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACAAT");
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(payload);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
        let (r2s, r2e) = e.r2.expect("long-enough payload should fall back to a length floor");
        assert_eq!(&o.seq[r2s..r2e], payload.as_slice());
        // Not a verified landmark -- honest accounting, unlike the real
        // anchor-found path.
        assert!(!e.adapter_trimmed);
    }

    #[test]
    fn flex_min_cdna_len_no_anchor_is_configurable() {
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let payload = b"GATTACAGATTACAGATTACA"; // 22bp -- below the real 65bp default
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACAAT");
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(payload);
        let qual = vec![b'I'; seq.len()];
        let o = orient(&seq, &qual, &ocfg());

        // Default (65) rejects a 22bp payload.
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
        assert!(e.r2.is_none());

        // Lowering the floor to <= payload.len() accepts it.
        let lenient = ExtractCfg { flex_min_cdna_len_no_anchor: payload.len(), ..flex_ecfg() };
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &lenient, false);
        assert!(e.r2.is_some());
    }

    #[test]
    fn flex_r2_boundary_skips_search_for_known_non_gex_reads() {
        // likely_non_gex=true: the caller's early i5-only sample-sheet
        // lookup (lib.rs::process_one) already identified this as e.g. a
        // CRISPR/Feature-Barcoding read sharing the same flex_r2_boundary=
        // true invocation as real Flex GEX reads. Must go straight to the
        // standard generic-adapter/read-end/min_cdna_len path -- same
        // outcome as non-Flex mode -- without ever attempting the (futile)
        // Flex anchor search.
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let payload = b"GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACAAT");
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(payload);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), true);
        let (r2s, r2e) = e.r2.expect("known-non-GEX reads must not be dropped");
        assert_eq!(r2e, o.seq.len());
        assert_eq!(&o.seq[r2s..r2e], payload.as_slice());
        assert!(!e.adapter_trimmed);
    }

    #[test]
    fn flex_r2_boundary_search_is_bounded_past_the_real_structural_range() {
        // The anchor sits real bp past cdna_start (FLEX_R2_END_ANCHOR_MAX_SEARCH,
        // ~150bp) -- a structurally-impossible position for the real
        // Ligated_Probe_insert/constant junction. Under the hard-gate default
        // (likely_non_gex=false) this must not be matched as the verified
        // landmark: confirms the search is actually bounded (an unbounded
        // scan would wrongly find it and mis-locate the boundary far too
        // late). The padding alone clears flex_min_cdna_len_no_anchor's real
        // 65bp floor though, so the length-floor fallback correctly still
        // emits an R2 -- bounded at the search limit, never at the real
        // (too-far-away) anchor position.
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let padding = vec![b'G'; FLEX_R2_END_ANCHOR_MAX_SEARCH + 20];
        let mut seq = Vec::new();
        seq.extend_from_slice(b"AACAAT");
        seq.extend_from_slice(P5);
        seq.extend_from_slice(b"ACGTACGTAC");
        seq.extend_from_slice(R1_PRIMER);
        seq.extend_from_slice(cb);
        seq.extend_from_slice(umi);
        seq.extend(std::iter::repeat(b'T').take(30));
        seq.extend_from_slice(&padding);
        seq.extend_from_slice(FLEX_R2_END_ANCHOR);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &flex_ecfg(), false);
        let (r2s, r2e) = e.r2.expect("length floor should still fire past the search bound");
        assert_eq!(r2e - r2s, FLEX_R2_END_ANCHOR_MAX_SEARCH, "bounded at the search limit, not the real anchor");
        assert!(!e.adapter_trimmed, "not a verified landmark match");
    }

    #[test]
    fn flex_r2_boundary_is_off_by_default_non_flex_path_unaffected() {
        // The same read that gets correctly bounded under flex_ecfg() should,
        // under the plain (non-Flex) ecfg(), fall through to the standard
        // generic-adapter/read-end behavior -- it has neither
        // ILLUMINA_3P_ADAPTER nor a length below min_cdna_len, so R2 just
        // runs to the read's own end, ignoring FLEX_R2_END_ANCHOR entirely.
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let pre_anchor = b"GATTACAGATTACAGATTACAGATTACAGATTACAGAT";
        let seq = synthetic_flex(cb, umi, pre_anchor, 30);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);
        let (r2s, r2e) = e.r2.expect("non-flex path should still emit R2 (long enough)");
        assert_eq!(r2e, o.seq.len());
        assert!(&o.seq[r2s..r2e] != pre_anchor.as_slice());
    }

    #[test]
    fn reversed_read_extracts_same() {
        let cb = b"AAACCCAAGAAACCTG";
        let umi = b"ACGTACGTACGT";
        let cdna = b"GGCATTACTGGCATTACTGGCATTACTGGCATTACTG";
        let canonical = synthetic(cb, umi, cdna, 30);
        let reversed = revcomp(&canonical);
        let qual = vec![b'I'; reversed.len()];

        let o = orient(&reversed, &qual, &ocfg());
        assert_eq!(o.orientation, Orientation::P7Led);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);
        let (r1s, r1e) = e.r1.unwrap();
        assert_eq!(r1e - r1s, 28);
        assert_eq!(&o.seq[r1s..r1s + 16], cb);
        let (r2s, r2e) = e.r2.unwrap();
        assert_eq!(&o.seq[r2s..r2e], cdna);
    }

    #[test]
    fn short_cdna_dropped() {
        let cb = b"AAACCCAAGAAACCCA";
        let umi = b"ACGTACGTACGT";
        let cdna = b"GGCATTAC"; // 8 bp < min_cdna_len
        let seq = synthetic(cb, umi, cdna, 30);
        let qual = vec![b'I'; seq.len()];

        let o = orient(&seq, &qual, &ocfg());
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg(), false);
        assert!(e.r2.is_none());
        assert!(e.r1.is_some());
    }
}
