//! Dev tool: is the CB whitelist still assignable after homopolymer collapse?
//! SBX indels concentrate in homopolymer runs. If we collapse runs (RLE) to absorb
//! them, we can only assign a CB when its collapsed form is UNIQUE in the whitelist.
//! Reports: homopolymer run-length stats, and the collision rate under two models:
//!   (1) full RLE collapse (most lossy upper bound)
//!   (2) length-bucketed RLE (run length capped at 1/2/3+ — tolerates +-1 in runs >=3)
//! Usage: cargo run --release --example hpc_whitelist -- <whitelist>

use std::collections::HashMap;
use std::env;

use sbx10x::barcode::{decode16, Whitelist};

// full RLE: "AAB CCC" -> bytes [A,B,C] (lengths dropped entirely)
fn rle_full(cb: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(cb.len());
    for &b in cb {
        if out.last() != Some(&b) {
            out.push(b);
        }
    }
    out
}

// bucketed RLE: keep base + capped run length (1,2, or 3 for >=3). Runs >=3 (where
// SBX +-1 indels live) all map to the same bucket, so a +-1 length error is absorbed
// only for long runs; short runs (1,2) stay distinct (less lossy than full collapse).
fn rle_bucketed(cb: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(cb.len() * 2);
    let mut i = 0;
    while i < cb.len() {
        let b = cb[i];
        let mut run = 0u8;
        while i < cb.len() && cb[i] == b {
            run += 1;
            i += 1;
        }
        out.push(b);
        out.push(run.min(3)); // 1,2,3+
    }
    out
}

fn collisions(forms: &[Vec<u8>]) -> (usize, usize) {
    let mut counts: HashMap<&Vec<u8>, u32> = HashMap::new();
    for f in forms {
        *counts.entry(f).or_default() += 1;
    }
    let distinct = counts.len();
    // barcodes whose collapsed form is shared by >=2 barcodes (not uniquely assignable)
    let colliding: usize = forms.iter().filter(|f| counts[*f] >= 2).count();
    (distinct, colliding)
}

fn main() -> anyhow::Result<()> {
    let wl = Whitelist::load(std::path::Path::new(&env::args().nth(1).unwrap()))?;
    let n = wl.codes().len();
    println!("whitelist barcodes: {n}");

    // homopolymer run-length distribution (max run per barcode)
    let mut max_run_hist = [0u64; 17];
    let mut any_run3 = 0u64;
    let mut any_run4 = 0u64;
    let mut full: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut bucket: Vec<Vec<u8>> = Vec::with_capacity(n);
    for &code in wl.codes() {
        let cb = decode16(code);
        let mut maxrun = 1usize;
        let mut i = 0;
        while i < 16 {
            let b = cb[i];
            let mut r = 0;
            while i < 16 && cb[i] == b {
                r += 1;
                i += 1;
            }
            maxrun = maxrun.max(r);
        }
        max_run_hist[maxrun.min(16)] += 1;
        if maxrun >= 3 { any_run3 += 1; }
        if maxrun >= 4 { any_run4 += 1; }
        full.push(rle_full(&cb));
        bucket.push(rle_bucketed(&cb));
    }

    println!("\nmax homopolymer run per barcode:");
    for (len, &c) in max_run_hist.iter().enumerate() {
        if c > 0 {
            println!("  run {len:2}: {c:8} ({:.1}%)", 100.0 * c as f64 / n as f64);
        }
    }
    println!("  barcodes with a run >=3: {any_run3} ({:.1}%)  ·  run >=4: {any_run4} ({:.2}%)",
        100.0 * any_run3 as f64 / n as f64, 100.0 * any_run4 as f64 / n as f64);

    let (df, cf) = collisions(&full);
    let (db, cb_) = collisions(&bucket);
    println!("\nCOLLAPSE FEASIBILITY (lower colliding = more assignable):");
    println!("  full RLE:      distinct forms {df:8} / {n}   colliding barcodes {cf} ({:.2}%)",
        100.0 * cf as f64 / n as f64);
    println!("  bucketed RLE:  distinct forms {db:8} / {n}   colliding barcodes {cb_} ({:.2}%)",
        100.0 * cb_ as f64 / n as f64);
    Ok(())
}
