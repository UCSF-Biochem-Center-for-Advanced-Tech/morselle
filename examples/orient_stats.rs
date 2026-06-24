//! Dev tool: run M1 orientation over a FASTQ and print summary stats.
//! Usage: cargo run --example orient_stats -- <fastq.gz> [max_reads]

use std::env;

use needletail::parse_fastx_file;
use sbx10x::orient::{orient, OrientCfg, OrientMethod};
use sbx10x::pipeline::Orientation;

fn main() -> anyhow::Result<()> {
    let mut args = env::args().skip(1);
    let path = args.next().expect("usage: orient_stats <fastq> [max]");
    let max: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let cfg = OrientCfg {
        cb_len: 16,
        umi_len: 12,
        min_homopolymer: 12,
        r1_tail_max_mm: 2,
        anchor_max_edits: 6,
        expected: None,
    };

    let (mut total, mut p5, mut p7, mut unknown) = (0u64, 0u64, 0u64, 0u64);
    let (mut homop, mut fallback, mut confirmed, mut leader) = (0u64, 0u64, 0u64, 0u64);

    let mut reader = parse_fastx_file(path)?;
    while let Some(rec) = reader.next() {
        if total as usize >= max {
            break;
        }
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        let o = orient(&seq, qual, &cfg);
        total += 1;
        match o.orientation {
            Orientation::P5Led => p5 += 1,
            Orientation::P7Led => p7 += 1,
            Orientation::Unknown => unknown += 1,
        }
        match o.method {
            OrientMethod::Homopolymer => homop += 1,
            OrientMethod::R1Fallback => fallback += 1,
            OrientMethod::None => {}
        }
        if o.r1_tail_confirmed {
            confirmed += 1;
        }
        if o.leader_present {
            leader += 1;
        }
    }

    let pct = |x: u64| 100.0 * x as f64 / total as f64;
    println!("reads            {total}");
    println!("P5-led (no flip) {p5:>8}  {:.1}%", pct(p5));
    println!("P7-led (flipped) {p7:>8}  {:.1}%", pct(p7));
    println!("Unknown          {unknown:>8}  {:.1}%", pct(unknown));
    println!("via homopolymer  {homop:>8}  {:.1}%", pct(homop));
    println!("via R1 fallback  {fallback:>8}  {:.1}%", pct(fallback));
    println!("R1-tail confirmed{confirmed:>8}  {:.1}%", pct(confirmed));
    println!("leader present   {leader:>8}  {:.1}%", pct(leader));
    Ok(())
}
