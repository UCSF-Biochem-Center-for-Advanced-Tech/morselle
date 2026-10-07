use crate::anchors::{hamming_leq, ILLUMINA_3P_ADAPTER};
use super::types::{
    ExtractCfg, Extracted, FLEX_R2_END_ANCHOR, FLEX_R2_END_ANCHOR_MAX_MM,
    FLEX_R2_END_ANCHOR_MAX_SEARCH, FLEX_R2_END_ANCHOR_PREFIX_LEN, MIN_PARTIAL_ADAPTER,
    MIN_PARTIAL_TSO, TSO_RC,
};

/// Where the reverse-complemented TSO starts, if this read carries one at the cDNA's 3' end.
fn find_tso(seq: &[u8], cdna_start: usize, cdna_end: usize) -> Option<usize> {
    let m = TSO_RC.len();
    if cdna_end >= cdna_start + m {
        let start = cdna_end - m;
        if hamming_leq(&seq[start..cdna_end], TSO_RC, (m / 10) as u32) {
            return Some(start);
        }
    }
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
            break;
        }
        i += 1;
    }
    last + 1
}

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
            if seq[i] == a[0]
                && seq[i + 1] == a[1]
                && seq[i + 2] == a[2]
                && adapter_single_indel(&seq[i..], a, max_mm)
            {
                return Some(i);
            }
        } else if remain >= MIN_PARTIAL_ADAPTER {
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

fn adapter_single_indel(s: &[u8], a: &[u8], max_mm: u32) -> bool {
    let m = a.len();
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

fn single_deletion_variants(seq: &[u8]) -> impl Iterator<Item = Vec<u8>> + '_ {
    (0..seq.len()).map(move |i| {
        let mut v = Vec::with_capacity(seq.len() - 1);
        v.extend_from_slice(&seq[..i]);
        v.extend_from_slice(&seq[i + 1..]);
        v
    })
}

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
            if let Some(cdna_end) = find_flex_r2_end(seq, cdna_start) {
                adapter_pos = Some(cdna_end);
                adapter_trimmed = true;
                if cdna_end > cdna_start {
                    r2 = Some((cdna_start, cdna_end));
                }
            } else {
                let cdna_end = seq.len().min(cdna_start + FLEX_R2_END_ANCHOR_MAX_SEARCH);
                if cdna_end >= cdna_start + cfg.flex_min_cdna_len_no_anchor {
                    r2 = Some((cdna_start, cdna_end));
                }
            }
        } else {
            adapter_pos = find_adapter(seq, cdna_start, cfg.adapter_max_mm);
            let mut cdna_end = adapter_pos.unwrap_or(seq.len());
            adapter_trimmed = adapter_pos.is_some();
            if cfg.trim_tso && cdna_end > cdna_start {
                if let Some(pos) = find_tso(seq, cdna_start, cdna_end) {
                    tso_pos = Some(pos);
                    tso_trimmed = true;
                    cdna_end = pos;
                }
            }
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
