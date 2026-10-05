//! Dev tool: is each run orientation-pure? Break down P5-led / P7-led by
//! location confidence (R1-tail confirmed + polyT flank). If confident reads are
//! ~one orientation per run, the wrong-orientation calls are false anchors.
//! Usage: cargo run --release --example orient_purity -- <fastq> [max]

use std::env;

use needletail::parse_fastx_file;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::Orientation;

fn main() -> anyhow::Result<()> {
    let mut args = env::args().skip(1);
    let fq = args.next().expect("usage: orient_purity <fastq> [max]");
    let max: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(40000);

    let ocfg = OrientCfg {
        cb_len: 16,
        umi_len: 12,
        min_homopolymer: 12,
        r1_tail_max_mm: 2,
        anchor_max_edits: 6,
        fallback: FallbackMode::HammingTail,
        expected: None,
    };
    let ecfg = ExtractCfg {
        min_cdna_len: 20,
        adapter_max_mm: 2,
        trim_tso: false,
        max_r2_len: 0,
        flex_r2_boundary: false,
        flex_min_cdna_len_no_anchor: 65,
    };

    // [confident?][orientation: 0=P5led,1=P7led,2=unknown]
    let mut tab = [[0u64; 3]; 2];
    let mut reader = parse_fastx_file(&fq)?;
    let mut n = 0;
    while let Some(rec) = reader.next() {
        if n >= max {
            break;
        }
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        n += 1;
        let o = orient(&seq, qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg, false);
        let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
        let oi = match o.orientation {
            Orientation::P5Led => 0,
            Orientation::P7Led => 1,
            Orientation::Unknown => 2,
        };
        tab[confident as usize][oi] += 1;
    }

    let row = |t: &[u64; 3], label: &str| {
        let s = (t[0] + t[1] + t[2]).max(1);
        println!("  {label:<22} P5-led {:>7} ({:>5.1}%)   P7-led {:>7} ({:>5.1}%)   Unknown {:>7} ({:>5.1}%)",
            t[0], 100.0 * t[0] as f64 / s as f64,
            t[1], 100.0 * t[1] as f64 / s as f64,
            t[2], 100.0 * t[2] as f64 / s as f64);
    };
    println!("reads {n}");
    row(&tab[1], "CONFIDENT location:");
    row(&tab[0], "unconfident:");
    Ok(())
}
