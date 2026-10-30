//! Where the UMI starts after a corrected cell barcode.
//!
//! The UMI is the `umi_len` bases right after the corrected CB *as it sits in
//! the read*. The CB end is where the whitelist CB aligns with the fewest
//! edits (route default ±1); when an end-of-CB indel and a substitution
//! explain the read equally well, the exact poly-T start breaks the tie.
//!
//! Why not anchor on the poly-T start (`Oriented::polyt_start`, as realign and
//! reframe did): that start tolerates isolated mismatches, so a UMI ending in
//! T (or T-rich, e.g. `..TGTA`) is partly absorbed into it and the UMI comes
//! out shifted. On CAT02 those UMIs never ended in T (0.0% / 2.6%, against
//! 26.6% for exact reads and 26.4% for Illumina UB); the positional rescues
//! ended in T too often (48-62%) when the CB end was taken from the route's
//! fixed frame.

/// Exact T run that marks the dT primer for the tie-break.
const MIN_T_RUN: usize = 6;

/// Start of the UMI in `seq` for a corrected `cb` that starts at `cb_start`;
/// `default_end` is the CB end the route's own reading implies.
pub(crate) fn umi_start(
    seq: &[u8],
    cb: &[u8],
    cb_start: usize,
    default_end: usize,
    umi_len: usize,
) -> usize {
    let lo = default_end.saturating_sub(1).max(cb_start + 1);
    let hi = (default_end + 1).min(seq.len());
    if hi < lo || cb_start >= seq.len() {
        return default_end;
    }
    let cost = end_costs(cb, &seq[cb_start..hi]);
    let q = t_run_start(seq, cb_start + cb.len() + umi_len / 2);
    (lo..=hi)
        .min_by_key(|&end| {
            let e = end + umi_len;
            let fit = match q {
                Some(q) if q <= e => e - q,
                Some(q) => 16 + q - e,
                None => 0,
            };
            (cost[end - cb_start], fit, end.abs_diff(default_end))
        })
        .unwrap_or(default_end)
}

/// Edit distance of all of `cb` against each prefix of `read`
/// (`out[j]` = distance to `read[..j]`).
fn end_costs(cb: &[u8], read: &[u8]) -> Vec<usize> {
    let mut row: Vec<usize> = (0..=read.len()).collect();
    for (i, &c) in cb.iter().enumerate() {
        let mut diag = row[0];
        row[0] = i + 1;
        for j in 1..=read.len() {
            let up = row[j];
            row[j] = (diag + usize::from(c != read[j - 1]))
                .min(up + 1)
                .min(row[j - 1] + 1);
            diag = up;
        }
    }
    row
}

/// First exact run of `MIN_T_RUN` T at or after `from`, within 24 bases.
fn t_run_start(seq: &[u8], from: usize) -> Option<usize> {
    let to = (from + 24).min(seq.len());
    (from..to.saturating_sub(MIN_T_RUN - 1))
        .find(|&i| seq[i..i + MIN_T_RUN].iter().all(|&b| b == b'T'))
}

#[cfg(test)]
#[path = "umi_anchor_tests.rs"]
mod tests;
