//! Per-file orientation pre-scan (`--orientation auto`), moved out of
//! `process.rs` to keep that file under the size ceiling.
use anyhow::Result;
use std::path::Path;

use crate::extract::{extract, ExtractCfg};
use crate::gzin::open_fastx_scan;
use crate::orient::{orient, OrientCfg};
use crate::pipeline::Orientation;

/// Detect a file's dominant orientation from a pre-scan of confident reads.
/// Returns `Some` only when ≥80% of confident reads agree (else `None` → per-read
/// auto). Each SBX run is orientation-pure, so this is normally unambiguous.
pub(crate) fn detect_orientation(
    input: &Path,
    base_cfg: &OrientCfg,
    extract_cfg: &ExtractCfg,
) -> Result<Option<Orientation>> {
    const SCAN: usize = 2000;
    let cfg = OrientCfg {
        expected: None,
        ..*base_cfg
    };
    let mut reader = open_fastx_scan(input)?;
    let (mut p5, mut p7) = (0u64, 0u64);
    let mut seen = 0;
    while let Some(rec) = reader.next() {
        if seen >= SCAN {
            break;
        }
        seen += 1;
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        let o = orient(&seq, qual, &cfg);
        // Only `polyt_run_len` below is read (set before the boundary branch)
        // -- this pre-scan never touches `e.r2`, so the routing hint is moot.
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, extract_cfg, false);
        if o.r1_tail_confirmed && e.polyt_run_len >= 10 {
            match o.orientation {
                Orientation::P5Led => p5 += 1,
                Orientation::P7Led => p7 += 1,
                Orientation::Unknown => {}
            }
        }
    }
    let total = p5 + p7;
    if total < 50 {
        return Ok(None); // too few confident reads to decide
    }
    let frac5 = p5 as f64 / total as f64;
    Ok(if frac5 >= 0.8 {
        Some(Orientation::P5Led)
    } else if frac5 <= 0.2 {
        Some(Orientation::P7Led)
    } else {
        None
    })
}
