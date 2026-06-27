//! Dev tool: distance-to-nearest-whitelist distribution + positional mismatch
//! profile for extracted cell barcodes, split by location confidence.
//! Usage: cargo run --release --example cb_distance -- <fastq> <whitelist> [max_reads]

use std::env;

use needletail::parse_fastx_file;
use rayon::prelude::*;
use sbx10x::barcode::{encode16, Whitelist};
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::Orientation;

/// Hamming distance (in bases) between two 2-bit-packed 16-mers.
#[inline]
fn dist16(a: u32, b: u32) -> u32 {
    let x = a ^ b;
    (((x | (x >> 1)) & 0x5555_5555).count_ones()) as u32
}

fn main() -> anyhow::Result<()> {
    let mut args = env::args().skip(1);
    let fq = args.next().expect("usage: cb_distance <fastq> <whitelist> [max]");
    let wlp = args.next().expect("need whitelist");
    let max: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(40000);

    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
    let codes = wl.codes();
    let ocfg = OrientCfg {
        cb_len: 16, umi_len: 12, min_homopolymer: 12, r1_tail_max_mm: 2,
        anchor_max_edits: 6, fallback: FallbackMode::HammingTail, expected: None,
    };
    let ecfg = ExtractCfg { min_cdna_len: 20, adapter_max_mm: 2, max_r2_len: 0 };

    // Collect (code, located, reframe_offset) for every read with a located CB.
    // reframe_offset = Some(k) if the offset-0 CB is off-whitelist but an EXACT
    // whitelist hit exists at a nearby frame shift k (|k| smallest first).
    let mut items: Vec<(u32, bool, Option<i64>)> = Vec::new();
    let mut reader = parse_fastx_file(&fq)?;
    let mut n = 0;
    while let Some(rec) = reader.next() {
        if n >= max { break; }
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        n += 1;
        let o = orient(&seq, qual, &ocfg);
        if o.orientation == Orientation::Unknown { continue; }
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);
        let Some((r1s, _)) = e.r1 else { continue; };
        let Some(code) = encode16(&o.seq[r1s..r1s + 16]) else { continue; };
        let located = o.r1_tail_confirmed && e.polyt_run_len >= 10;
        let mut reframe = None;
        if !wl.contains(&o.seq[r1s..r1s + 16]) {
            for k in [1i64, -1, 2, -2, 3, -3, 4, -4] {
                let s = r1s as i64 + k;
                if s < 0 { continue; }
                let s = s as usize;
                if s + 16 <= o.seq.len() && wl.contains(&o.seq[s..s + 16]) {
                    reframe = Some(k);
                    break;
                }
            }
        }
        items.push((code, located, reframe));
    }

    // Nearest whitelist distance + nearest code (brute force, parallel).
    let results: Vec<(u32, u32, bool)> = items
        .par_iter()
        .map(|&(code, located, _)| {
            let mut best = 99u32;
            let mut best_code = 0u32;
            for &w in codes {
                let d = dist16(code, w);
                if d < best { best = d; best_code = w; if d == 0 { break; } }
            }
            (best, code ^ best_code, located)
        })
        .collect();

    // Distance histogram split by confidence.
    let mut hist = [[0u64; 2]; 10]; // [dist capped 9][located]
    let mut pos_conf = [0u64; 16]; // positional mismatches, confident, d in 1..=3
    let mut pos_far = [0u64; 16];  // confident, d >= 4
    let (mut nconf, mut nunc) = (0u64, 0u64);
    for &(d, xor, located) in &results {
        let li = if located { 0 } else { 1 };
        hist[d.min(9) as usize][li] += 1;
        if located { nconf += 1; } else { nunc += 1; }
        if located && d >= 1 {
            for i in 0..16 {
                let shift = 2 * (15 - i);
                if (xor >> shift) & 0b11 != 0 {
                    if d <= 3 { pos_conf[i] += 1; } else { pos_far[i] += 1; }
                }
            }
        }
    }

    println!("reads scanned {n}, located CBs {} (confident {nconf}, uncertain {nunc})", items.len());
    println!("\ndistance-to-nearest-whitelist  (count and % within group)");
    println!("  d :  confident            uncertain");
    for d in 0..10 {
        let c = hist[d][0]; let u = hist[d][1];
        let label = if d == 9 { "9+".to_string() } else { d.to_string() };
        println!("  {label:<2}:  {c:>9} ({:>5.1}%)   {u:>9} ({:>5.1}%)",
            100.0 * c as f64 / nconf.max(1) as f64, 100.0 * u as f64 / nunc.max(1) as f64);
    }

    println!("\npositional mismatch profile, CONFIDENT reads (CB position 0=5' … 15=3' of canonical)");
    let pmax = *pos_conf.iter().max().unwrap_or(&1) as f64;
    for i in 0..16 {
        let bar = "#".repeat((40.0 * pos_conf[i] as f64 / pmax.max(1.0)) as usize);
        println!("  pos {i:>2} (d≤3): {:>8}  {bar}", pos_conf[i]);
    }
    println!("\n  (confident d≥4 positional totals: {pos_far:?})");

    // Reframe recovery: of off-whitelist CBs, how many get an EXACT hit at a
    // nearby frame shift? Split by location confidence + offset histogram.
    let mut off0 = [0u64; 2];      // off-whitelist at offset 0, [located]
    let mut reframed = [0u64; 2];  // recovered by a nonzero offset, [located]
    let mut offset_hist = std::collections::BTreeMap::<i64, u64>::new();
    for &(code, located, reframe) in &items {
        if wl.codes().binary_search(&code).is_ok() {
            continue; // already exact at offset 0
        }
        let li = if located { 0 } else { 1 };
        off0[li] += 1;
        if let Some(k) = reframe {
            reframed[li] += 1;
            *offset_hist.entry(k).or_default() += 1;
        }
    }
    println!("\nreframe recovery (exact whitelist hit at a nearby frame shift)");
    println!("  confident: {} of {} off-whitelist recovered ({:.1}%)",
        reframed[0], off0[0], 100.0 * reframed[0] as f64 / off0[0].max(1) as f64);
    println!("  uncertain: {} of {} off-whitelist recovered ({:.1}%)",
        reframed[1], off0[1], 100.0 * reframed[1] as f64 / off0[1].max(1) as f64);
    println!("  winning offset histogram: {offset_hist:?}");
    Ok(())
}
