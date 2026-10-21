//! Dev tool: are "no cDNA" reads just short molecules?
//! Among reads with R1 located, split by whether R2 (cDNA >= min_cdna_len) was
//! extracted, and compare raw read-length distributions.
//! Usage: cargo run --release --example nocdna_len -- <fastq> <p5-led|p7-led> [max]

use std::env;

use morselle::extract::{extract, ExtractCfg};
use morselle::orient::{orient, FallbackMode, OrientCfg};
use morselle::pipeline::Orientation;
use needletail::parse_fastx_file;

fn stats(name: &str, v: &mut [u64]) {
    if v.is_empty() {
        println!("  {name}: n=0");
        return;
    }
    v.sort_unstable();
    let n = v.len();
    let sum: u64 = v.iter().sum();
    let mean = sum as f64 / n as f64;
    let pct = |p: f64| v[((n as f64 * p) as usize).min(n - 1)];
    println!(
        "  {name:14} n={n:7}  mean={mean:6.0}  p10={:5}  p25={:5}  median={:5}  p75={:5}  p90={:5}  max={:5}",
        pct(0.10), pct(0.25), pct(0.50), pct(0.75), pct(0.90), v[n - 1]
    );
}

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a.next().unwrap();
    let orient_arg = a.next().unwrap();
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let expected = if orient_arg == "p7-led" {
        Orientation::P7Led
    } else {
        Orientation::P5Led
    };
    let ocfg = OrientCfg {
        cb_len: 16,
        umi_len: 12,
        min_homopolymer: 12,
        r1_tail_max_mm: 2,
        anchor_max_edits: 6,
        fallback: FallbackMode::HammingTail,
        expected: Some(expected),
    };
    let ecfg = ExtractCfg {
        min_cdna_len: 50,
        adapter_max_mm: 2,
        trim_tso: false,
        max_r2_len: 0,
        flex_r2_boundary: false,
        flex_min_cdna_len_no_anchor: 65,
    };

    let mut all: Vec<u64> = Vec::new();
    let mut r1ok_r2ok: Vec<u64> = Vec::new();
    let mut r1ok_nocdna: Vec<u64> = Vec::new();
    // for no-cDNA reads, how much sequence sits 3' of the polyT end (would-be cDNA window)?
    let mut tail_after_polyt: Vec<u64> = Vec::new();
    let mut adapter_present = 0u64;

    let mut reader = parse_fastx_file(&fq)?;
    let mut n = 0usize;
    while let Some(rec) = reader.next() {
        if n >= max {
            break;
        }
        n += 1;
        let rec = rec?;
        let qual = rec.qual().unwrap_or(&[]).to_vec();
        let raw = rec.seq();
        let len = raw.len() as u64;
        all.push(len);
        let o = orient(&raw, &qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg, false);
        if e.r1.is_none() {
            continue;
        }
        if e.r2.is_some() {
            r1ok_r2ok.push(len);
        } else {
            r1ok_nocdna.push(len);
            if let Some(pe) = e.polyt_end {
                tail_after_polyt.push((o.seq.len().saturating_sub(pe)) as u64);
            }
            if e.adapter_pos.is_some() {
                adapter_present += 1;
            }
        }
    }

    println!("reads {n}");
    stats("ALL", &mut all);
    stats("R1+R2 ok", &mut r1ok_r2ok);
    stats("R1, no cDNA", &mut r1ok_nocdna);
    println!(
        "\n  no-cDNA reads with a 3' adapter detected: {adapter_present} ({:.1}% of no-cDNA)",
        100.0 * adapter_present as f64 / r1ok_nocdna.len().max(1) as f64
    );
    stats("3'-of-polyT bp (no-cDNA)", &mut tail_after_polyt);
    Ok(())
}
