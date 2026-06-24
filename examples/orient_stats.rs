//! Dev tool: run M1 orientation over a FASTQ and print summary stats + timing.
//! Usage: cargo run --release --example orient_stats -- <fastq.gz> <off|tail|lev> [max_reads]

use std::env;
use std::time::Instant;

use needletail::parse_fastx_file;
use sbx10x::orient::{orient, FallbackMode, OrientCfg, OrientMethod};
use sbx10x::pipeline::Orientation;

fn main() -> anyhow::Result<()> {
    let mut args = env::args().skip(1);
    let path = args.next().expect("usage: orient_stats <fastq> <off|tail|lev> [max]");
    let mode = match args.next().as_deref() {
        Some("off") => FallbackMode::Off,
        Some("lev") => FallbackMode::Levenshtein,
        _ => FallbackMode::HammingTail,
    };
    let max: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let cfg = OrientCfg {
        cb_len: 16,
        umi_len: 12,
        min_homopolymer: 12,
        r1_tail_max_mm: 2,
        anchor_max_edits: 6,
        fallback: mode,
        expected: None,
    };

    let (mut total, mut p5, mut p7, mut unknown) = (0u64, 0u64, 0u64, 0u64);
    let (mut homop, mut tail, mut lev, mut confirmed) = (0u64, 0u64, 0u64, 0u64);
    let mut orient_nanos = 0u128;

    let mut reader = parse_fastx_file(path)?;
    while let Some(rec) = reader.next() {
        if total as usize >= max {
            break;
        }
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        let t0 = Instant::now();
        let o = orient(&seq, qual, &cfg);
        orient_nanos += t0.elapsed().as_nanos();
        total += 1;
        match o.orientation {
            Orientation::P5Led => p5 += 1,
            Orientation::P7Led => p7 += 1,
            Orientation::Unknown => unknown += 1,
        }
        match o.method {
            OrientMethod::Homopolymer => homop += 1,
            OrientMethod::TailFallback => tail += 1,
            OrientMethod::LevenshteinFallback => lev += 1,
            OrientMethod::None => {}
        }
        if o.r1_tail_confirmed {
            confirmed += 1;
        }
    }

    let pct = |x: u64| 100.0 * x as f64 / total as f64;
    println!("mode={mode:?}  reads={total}");
    println!("  P5-led {p5} ({:.1}%)  P7-led {p7} ({:.1}%)  Unknown {unknown} ({:.1}%)",
        pct(p5), pct(p7), pct(unknown));
    println!("  via homopolymer {:.1}%  via tail {:.1}%  via lev {:.1}%  confirmed {:.1}%",
        pct(homop), pct(tail), pct(lev), pct(confirmed));
    println!("  orient() time: {:.2} ms total, {:.0} ns/read, {:.2} M reads/s (single thread)",
        orient_nanos as f64 / 1e6,
        orient_nanos as f64 / total as f64,
        total as f64 / (orient_nanos as f64 / 1e9) / 1e6);
    Ok(())
}
