//! Dev tool: how many MORE reads could we recover at the CB gate?
//! For every CB-fail read that is otherwise emittable (R1 located + R2 present),
//! test the levers that mirror the index-reframe win:
//!   (A) CB reframe   — unique within-1 whitelist hit at a shifted frame (k != 0)
//!   (B) lenient w2    — unique within-2 hit at k=0 (current --lenient, confident only)
//!   (C) lenient w2 ALL — same but ignoring location confidence (precision risk)
//! Split by location confidence so genuine base errors are separated from framing.
//! Usage: cargo run --release --example cb_headroom -- <fastq> <wl> <p5-led|p7-led> [max]

use std::collections::BTreeMap;
use std::env;

use needletail::parse_fastx_file;
use sbx10x::barcode::Whitelist;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};

const FLANK: u16 = 10;

// unique within-`maxmm` whitelist hit for a 16-mer; returns true if exactly one.
fn unique_within(wl: &Whitelist, cb: &[u8], maxmm: u32) -> bool {
    if cb.len() != 16 {
        return false;
    }
    // exact?
    if wl.contains(cb) {
        return true;
    }
    // count whitelist members within maxmm by probing neighbors is expensive;
    // instead scan: classify gives exact/within1; for within2 use correct_within2.
    match maxmm {
        1 => matches!(wl.classify(cb), CbStatus::Within1),
        _ => wl.correct_within2(cb).is_some(),
    }
}

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a.next().unwrap();
    let wl = Whitelist::load(std::path::Path::new(&a.next().unwrap()))?;
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
        max_r2_len: 0,
    };

    let mut reads = 0u64;
    let mut emittable = 0u64; // R1 located + R2 present
    let mut cb_ok = 0u64; // already exact/within1 at k=0
    let mut fail = 0u64; // CB-fail among emittable
    let (mut fail_conf, mut fail_unconf) = (0u64, 0u64);

    // recovery counters (among CB-fail emittable)
    let mut reframe_w1 = 0u64; // unique within-1 at some k != 0
    let mut reframe_conf = 0u64;
    let mut lenient_w2_conf = 0u64; // unique within-2 at k=0, confident
    let mut lenient_w2_all = 0u64; // unique within-2 at k=0, any
    let mut off_hist: BTreeMap<i64, u64> = BTreeMap::new();

    let mut reader = parse_fastx_file(&fq)?;
    while let Some(rec) = reader.next() {
        if reads >= max as u64 {
            break;
        }
        let rec = rec?;
        let qual = rec.qual().unwrap_or(&[]).to_vec();
        let raw = rec.seq();
        reads += 1;
        let o = orient(&raw, &qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);
        let (Some((r1s, _)), true) = (e.r1, e.r2.is_some()) else {
            continue;
        };
        emittable += 1;
        let cb = &o.seq[r1s..r1s + 16];
        if !matches!(wl.classify(cb), CbStatus::Fail) {
            cb_ok += 1;
            continue;
        }
        fail += 1;
        let conf = o.r1_tail_confirmed && e.polyt_run_len >= FLANK;
        if conf {
            fail_conf += 1;
        } else {
            fail_unconf += 1;
        }

        // (A) CB reframe: unique within-1 at a shifted frame k != 0
        let mut framed = false;
        for k in [1i64, -1, 2, -2] {
            let s = r1s as i64 + k;
            if s < 0 {
                continue;
            }
            let s = s as usize;
            if s + 16 > o.seq.len() {
                continue;
            }
            if unique_within(&wl, &o.seq[s..s + 16], 1) {
                reframe_w1 += 1;
                if conf {
                    reframe_conf += 1;
                }
                *off_hist.entry(k).or_default() += 1;
                framed = true;
                break;
            }
        }
        let _ = framed;

        // (B,C) lenient within-2 at k=0
        if unique_within(&wl, cb, 2) {
            lenient_w2_all += 1;
            if conf {
                lenient_w2_conf += 1;
            }
        }
    }

    let pe = |x: u64| 100.0 * x as f64 / emittable.max(1) as f64;
    let pr = |x: u64| 100.0 * x as f64 / reads.max(1) as f64;
    println!("reads {reads}");
    println!(
        "emittable (R1 located + R2 present): {emittable} ({:.1}% of reads)",
        pr(emittable)
    );
    println!(
        "  CB ok (exact/within1):  {cb_ok} ({:.1}% of emittable)",
        pe(cb_ok)
    );
    println!(
        "  CB fail:                {fail} ({:.1}% of emittable)",
        pe(fail)
    );
    println!("    of fails — confident location: {fail_conf} ({:.1}%) · uncertain: {fail_unconf} ({:.1}%)",
        100.0 * fail_conf as f64 / fail.max(1) as f64, 100.0 * fail_unconf as f64 / fail.max(1) as f64);
    println!("\n  RECOVERY HEADROOM (unique whitelist hit; among CB-fail emittable):");
    println!("    (A) CB reframe (within-1 at shifted frame): {reframe_w1}  = {:.2}% of reads  [{} confident]",
        pr(reframe_w1), reframe_conf);
    println!("        winning offset histogram: {off_hist:?}");
    println!(
        "    (B) lenient within-2 @k0, confident only:   {lenient_w2_conf}  = {:.2}% of reads",
        pr(lenient_w2_conf)
    );
    println!(
        "    (C) lenient within-2 @k0, ALL (precision risk): {lenient_w2_all}  = {:.2}% of reads",
        pr(lenient_w2_all)
    );
    Ok(())
}
