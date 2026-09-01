//! Dev tool: profile the CB-fail bucket — are the failing 16 bp barcodes seen at
//! HIGH counts (incomplete whitelist / artifact) or flat (random degradation)?
//! And do they look like adapter / primer / homopolymer sequence?
//! Usage: cargo run --release --example cbfail_profile -- <fastq> <whitelist> <p5-led|p7-led> [max] [topN]

use std::collections::HashMap;
use std::env;

use needletail::parse_fastx_file;
use sbx10x::anchors::{
    revcomp, ILLUMINA_3P_ADAPTER, P5, P7, R1_PRIMER, R2_PRIMER, SBX_LEADER, TSO,
};
use sbx10x::barcode::Whitelist;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};

fn ham(a: &[u8], b: &[u8]) -> u32 {
    a.iter().zip(b).filter(|(x, y)| x != y).count() as u32
}

// min Hamming of `cb` (16) against any 16 bp window of `r` or its revcomp.
fn min_window(cb: &[u8], r: &[u8]) -> u32 {
    let rc = revcomp(r);
    let mut best = 99;
    for s in [r, &rc] {
        if s.len() >= 16 {
            for w in s.windows(16) {
                best = best.min(ham(cb, w));
            }
        }
    }
    best
}

fn max_run(cb: &[u8]) -> usize {
    let (mut best, mut cur) = (1, 1);
    for i in 1..cb.len() {
        if cb[i] == cb[i - 1] {
            cur += 1;
            best = best.max(cur);
        } else {
            cur = 1;
        }
    }
    best
}

fn adapter_tag(cb: &[u8]) -> String {
    let refs: [(&str, &[u8]); 6] = [
        ("R1primer", R1_PRIMER),
        ("R2primer", R2_PRIMER),
        ("P5", P5),
        ("P7", P7),
        ("3pAdapter", ILLUMINA_3P_ADAPTER),
        ("TSO", TSO),
    ];
    let mut hits = vec![];
    for (name, r) in refs {
        let d = min_window(cb, r);
        if d <= 2 {
            hits.push(format!("{name}~{d}"));
        }
    }
    if min_window(cb, SBX_LEADER) <= 2 {
        hits.push("leader".into());
    }
    let mr = max_run(cb);
    if mr >= 8 {
        hits.push(format!("homopolymer{mr}"));
    }
    if hits.is_empty() {
        "-".into()
    } else {
        hits.join(",")
    }
}

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a.next().unwrap();
    let wl = Whitelist::load(std::path::Path::new(&a.next().unwrap()))?;
    let orient_arg = a.next().unwrap();
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);
    let topn: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(25);

    let (cb_len, umi_len) = (16usize, 12usize);
    let expected = if orient_arg == "p7-led" {
        Orientation::P7Led
    } else {
        Orientation::P5Led
    };
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
    };

    let mut freq: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut total = 0u64;
    let mut n = 0usize;

    let mut reader = parse_fastx_file(&fq)?;
    while let Some(rec) = reader.next() {
        if n >= max {
            break;
        }
        n += 1;
        let rec = rec?;
        let q = rec.qual().unwrap_or(&[]).to_vec();
        let raw = rec.seq();
        let o = orient(&raw, &q, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);
        let Some((r1s, _)) = e.r1 else { continue };
        if e.r2.is_none() {
            continue;
        }
        let cb = &o.seq[r1s..r1s + cb_len];
        let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
        let mut cb_ok = !matches!(wl.classify(cb), CbStatus::Fail);
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
        if !cb_ok {
            *freq.entry(cb.to_vec()).or_default() += 1;
            total += 1;
        }
    }

    let distinct = freq.len();
    let singletons = freq.values().filter(|&&c| c == 1).count();
    let mut v: Vec<(Vec<u8>, u64)> = freq.into_iter().collect();
    v.sort_by(|a, b| b.1.cmp(&a.1));
    let cum = |k: usize| -> f64 {
        100.0 * v.iter().take(k).map(|x| x.1).sum::<u64>() as f64 / total.max(1) as f64
    };

    println!("CB-fail reads: {total}   distinct CBs: {distinct}   singletons: {singletons} ({:.1}% of distinct)",
        100.0 * singletons as f64 / distinct.max(1) as f64);
    println!(
        "concentration: top-1 {:.2}%  top-10 {:.2}%  top-100 {:.2}%  top-1000 {:.2}%",
        cum(1),
        cum(10),
        cum(100),
        cum(1000)
    );
    println!("  (flat/singleton-dominated => random degradation; top-heavy => recurring seq / artifact)\n");
    println!(
        "{:<18} {:>7} {:>6}   {:<8} {:>6}  {}",
        "rank/CB", "count", "%", "maxrun", "wl_d2", "adapter/motif"
    );
    for (i, (cb, c)) in v.iter().take(topn).enumerate() {
        let s = String::from_utf8_lossy(cb);
        let within2 = if wl.correct_within2(cb).is_some() {
            "yes"
        } else {
            "-"
        };
        println!(
            "{:>2} {:<15} {:>7} {:>5.2}%   run{:<5} {:>6}  {}",
            i + 1,
            s,
            c,
            100.0 * *c as f64 / total.max(1) as f64,
            max_run(cb),
            within2,
            adapter_tag(cb)
        );
    }
    Ok(())
}
