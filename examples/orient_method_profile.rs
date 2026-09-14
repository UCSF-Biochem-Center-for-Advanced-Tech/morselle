//! Dev tool: which orient() path is actually firing, and how long is the
//! detected dT/polyA run when it does? Prompted by a visual spot-check on
//! Flex reads not showing an obvious poly-T. Usage:
//!   cargo run --release --example orient_method_profile -- <fastq> [<fastq> ...]

use std::env;

use needletail::parse_fastx_file;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg, OrientMethod};

fn main() -> anyhow::Result<()> {
    let fastqs: Vec<String> = env::args().skip(1).collect();
    anyhow::ensure!(!fastqs.is_empty(), "usage: orient_method_profile <fastq> [...]");

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
    };

    let mut n = 0u64;
    let mut homopolymer = 0u64;
    let mut tail_fallback = 0u64;
    let mut levenshtein = 0u64;
    let mut none = 0u64;
    // polyt_run_len histogram, bucketed
    let mut run_len_0 = 0u64;
    let mut run_len_1_9 = 0u64;
    let mut run_len_10_19 = 0u64;
    let mut run_len_20_plus = 0u64;
    let mut run_len_sum = 0u64;
    let mut run_len_n = 0u64;

    for fq in &fastqs {
        let mut reader = parse_fastx_file(fq)?;
        while let Some(rec) = reader.next() {
            let rec = rec?;
            let seq = rec.seq();
            let qual = rec.qual().unwrap_or(&[]);
            n += 1;
            let o = orient(&seq, qual, &ocfg);
            match o.method {
                OrientMethod::Homopolymer => homopolymer += 1,
                OrientMethod::TailFallback => tail_fallback += 1,
                OrientMethod::LevenshteinFallback => levenshtein += 1,
                OrientMethod::None => {
                    none += 1;
                    continue;
                }
            }
            let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);
            let rl = e.polyt_run_len;
            run_len_sum += rl as u64;
            run_len_n += 1;
            match rl {
                0 => run_len_0 += 1,
                1..=9 => run_len_1_9 += 1,
                10..=19 => run_len_10_19 += 1,
                _ => run_len_20_plus += 1,
            }
        }
    }

    println!("reads: {n}");
    let p = |x: u64| 100.0 * x as f64 / n.max(1) as f64;
    println!("\n--- orient() method ---");
    println!("Homopolymer (real T/A run found): {homopolymer:>10} ({:.1}%)", p(homopolymer));
    println!("TailFallback (no run; R1-tail scan instead): {tail_fallback:>10} ({:.1}%)", p(tail_fallback));
    println!("LevenshteinFallback: {levenshtein:>10} ({:.1}%)", p(levenshtein));
    println!("None (unoriented): {none:>10} ({:.1}%)", p(none));

    println!("\n--- extract()'s measured polyt_run_len, among oriented reads ---");
    println!("run_len == 0:   {run_len_0:>10} ({:.1}%)", p(run_len_0));
    println!("run_len 1-9:    {run_len_1_9:>10} ({:.1}%)", p(run_len_1_9));
    println!("run_len 10-19:  {run_len_10_19:>10} ({:.1}%)", p(run_len_10_19));
    println!("run_len 20+:    {run_len_20_plus:>10} ({:.1}%)", p(run_len_20_plus));
    if run_len_n > 0 {
        println!("mean run_len (oriented reads): {:.2}", run_len_sum as f64 / run_len_n as f64);
    }

    Ok(())
}
