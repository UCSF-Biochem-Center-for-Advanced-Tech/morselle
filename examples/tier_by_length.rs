//! Dev tool: cross-tabulate per-read pipeline tier against RAW read length, and
//! report the true cDNA insert-length distribution for reads whose barcode was
//! located. Answers: is a short read a truncated long molecule (sequencing-side,
//! "Unparseable") or a complete short-insert molecule (library-side, "NoCdna")?
//! Usage: cargo run --release --example tier_by_length -- <fastq> <whitelist> <p5-led|p7-led> [max]

use std::env;

use needletail::parse_fastx_file;
use sbx10x::barcode::Whitelist;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};

#[derive(Clone, Copy, PartialEq)]
enum Tier {
    Unparseable = 0, // no CB bracket located (orient/anchor fail) -> truncated/junk
    NoCdna = 1,      // CB located but cDNA absent/too short -> short-insert molecule
    CbFail = 2,      // R1+R2 present, CB degraded/uncertain
    ValidPair = 3,   // R1+R2 present, CB ok
}

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a
        .next()
        .expect("usage: tier_by_length <fastq> <whitelist> <p5-led|p7-led> [max]");
    let wlp = a.next().expect("whitelist");
    let orient_arg = a.next().unwrap_or_else(|| "p5-led".into());
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
    let expected = match orient_arg.as_str() {
        "p7-led" => Orientation::P7Led,
        _ => Orientation::P5Led,
    };
    let (cb_len, umi_len) = (16usize, 12usize);
    let ocfg = OrientCfg {
        cb_len,
        umi_len,
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

    const NB: usize = 64; // 25bp bins, 0..1600
    const W: usize = 25;
    let mut grid = [[0u64; 4]; NB]; // [len_bin][tier]
    let mut tot = [0u64; 4];
    // true cDNA insert length histogram among reads with a located barcode
    let mut insert_hist = [0u64; NB];
    let mut insert_with_bc = 0u64;
    let mut total = 0u64;

    let mut reader = parse_fastx_file(&fq)?;
    while let Some(rec) = reader.next() {
        if total as usize >= max {
            break;
        }
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        let raw_len = seq.len();
        total += 1;

        let o = orient(&seq, qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg, false);

        let tier = classify(&o, &e, &wl, cb_len, umi_len);
        let b = (raw_len / W).min(NB - 1);
        grid[b][tier as usize] += 1;
        tot[tier as usize] += 1;

        // true insert length = cDNA span (polyt_end .. adapter_pos/seq_end), when barcode located
        if e.r1.is_some() {
            if let Some(pe) = e.polyt_end {
                let cend = e.adapter_pos.unwrap_or(o.seq.len());
                let insert = cend.saturating_sub(pe);
                insert_hist[(insert / W).min(NB - 1)] += 1;
                insert_with_bc += 1;
            }
        }
    }

    let pc = |x: u64, d: u64| {
        if d > 0 {
            100.0 * x as f64 / d as f64
        } else {
            0.0
        }
    };
    println!("reads scanned: {total}\n");
    println!("OVERALL TIERS");
    let tnames = [
        "Unparseable (truncated/no barcode)",
        "NoCdna (barcode ok, empty/short insert)",
        "CbFail (barcode degraded)",
        "ValidPair (recoverable)",
    ];
    for t in 0..4 {
        println!(
            "  {:<42} {:>12}  {:>5.1}%",
            tnames[t],
            tot[t],
            pc(tot[t], total)
        );
    }

    println!("\nTIER MIX BY RAW READ LENGTH (row % across tiers)");
    println!(
        "{:>9}  {:>9}  {:>7} {:>7} {:>7} {:>7}",
        "len_bin", "reads", "Unpars", "NoCdna", "CbFail", "Valid"
    );
    for b in 1..NB {
        let rb: u64 = (0..4).map(|t| grid[b][t]).sum();
        if rb < 200 {
            continue;
        }
        let lo = b * W;
        println!(
            "{:>4}-{:<4} {:>9}  {:>6.1}% {:>6.1}% {:>6.1}% {:>6.1}%",
            lo,
            lo + W - 1,
            rb,
            pc(grid[b][0], rb),
            pc(grid[b][1], rb),
            pc(grid[b][2], rb),
            pc(grid[b][3], rb)
        );
    }

    println!("\nTRUE cDNA INSERT LENGTH (reads with located barcode; n={insert_with_bc})");
    for b in 0..NB {
        if insert_hist[b] == 0 {
            continue;
        }
        let lo = b * W;
        if lo > 900 {
            break;
        }
        println!(
            "  {:>4}-{:<4} {:>10}  {:>5.1}%",
            lo,
            lo + W - 1,
            insert_hist[b],
            pc(insert_hist[b], insert_with_bc)
        );
    }
    Ok(())
}

fn classify(
    o: &sbx10x::orient::Oriented,
    e: &sbx10x::extract::Extracted,
    wl: &Whitelist,
    cb_len: usize,
    umi_len: usize,
) -> Tier {
    let Some((r1s, _)) = e.r1 else {
        return Tier::Unparseable;
    };
    if e.r2.is_none() {
        return Tier::NoCdna;
    }
    let cb = &o.seq[r1s..r1s + cb_len];
    let status = wl.classify(cb);
    let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
    let mut cb_ok = !matches!(status, CbStatus::Fail);
    if !cb_ok {
        if let Some(pts) = o.polyt_start {
            let delta = (pts as i64 - r1s as i64) - (cb_len + umi_len) as i64;
            if (delta == 1 || delta == -1)
                && pts >= umi_len
                && pts <= o.seq.len()
                && wl.realign_cb(&o.seq[r1s..], delta as i32).is_some()
            {
                cb_ok = true;
            }
        }
        if !cb_ok && confident && wl.reframe_cb(&o.seq, r1s, 2).is_some() {
            cb_ok = true;
        }
        if !cb_ok && confident && wl.correct_within2(cb).is_some() {
            cb_ok = true;
        }
    }
    if cb_ok {
        Tier::ValidPair
    } else {
        Tier::CbFail
    }
}
