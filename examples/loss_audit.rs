//! Dev tool: stage-by-stage loss funnel + CB-fail recoverability analysis.
//! Usage: cargo run --release --example loss_audit -- <fastq> <whitelist> [max_reads]

use std::env;

use needletail::parse_fastx_file;
use sbx10x::barcode::Whitelist;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};

fn main() -> anyhow::Result<()> {
    let mut args = env::args().skip(1);
    let fq = args.next().expect("usage: loss_audit <fastq> <whitelist> [max]");
    let wlp = args.next().expect("need whitelist");
    let max: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
    let ocfg = OrientCfg {
        cb_len: 16, umi_len: 12, min_homopolymer: 12, r1_tail_max_mm: 2,
        anchor_max_edits: 6, fallback: FallbackMode::HammingTail, expected: None,
    };
    let ecfg = ExtractCfg { min_cdna_len: 20, adapter_max_mm: 2, max_r2_len: 0 };

    let (mut total, mut unparseable, mut no_r2, mut emit_ok) = (0u64, 0u64, 0u64, 0u64);
    let (mut cb_fail_r2, mut tail_conf_among_fail) = (0u64, 0u64);
    // CB-fail breakdown
    let (mut d2_unique, mut d2_ambig, mut far) = (0u64, 0u64, 0u64);
    // location confidence among CB-fails
    let (mut fail_located, mut fail_located_d2unique) = (0u64, 0u64);
    let mut span_delta = std::collections::BTreeMap::<i64, u64>::new();

    let mut reader = parse_fastx_file(&fq)?;
    while let Some(rec) = reader.next() {
        if total as usize >= max { break; }
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        total += 1;

        let o = orient(&seq, qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);

        let Some((r1s, _)) = e.r1 else { unparseable += 1; continue; };
        if o.orientation == Orientation::Unknown { unparseable += 1; continue; }

        let cb = &o.seq[r1s..r1s + 16];
        let status = wl.classify(cb);
        let has_r2 = e.r2.is_some();

        match status {
            CbStatus::OnWhitelist | CbStatus::Within1 => {
                if has_r2 { emit_ok += 1; } else { no_r2 += 1; }
            }
            _ => {
                // CB fails the whitelist gate
                if !has_r2 { no_r2 += 1; continue; }
                cb_fail_r2 += 1;
                let located = o.r1_tail_confirmed && e.polyt_run_len >= 10;
                if o.r1_tail_confirmed { tail_conf_among_fail += 1; }
                if located { fail_located += 1; }
                if let Some(pts) = o.polyt_start {
                    let delta = pts as i64 - r1s as i64 - 28;
                    *span_delta.entry(delta.clamp(-5, 5)).or_default() += 1;
                }
                // Hamming-2 recoverability: unique whitelist hit within 2 subs?
                match unique_within2(&wl, cb) {
                    Rec::Unique => {
                        d2_unique += 1;
                        if located { fail_located_d2unique += 1; }
                    }
                    Rec::Ambiguous => d2_ambig += 1,
                    Rec::None => far += 1,
                }
            }
        }
    }

    let pc = |x: u64| 100.0 * x as f64 / total as f64;
    println!("=== funnel ({total} reads) ===");
    println!("  unparseable (no orient/CB)   {unparseable:>10}  {:.1}%", pc(unparseable));
    println!("  no R2 (cDNA missing/short)   {no_r2:>10}  {:.1}%", pc(no_r2));
    println!("  CB-fail (R1+R2, gated)       {cb_fail_r2:>10}  {:.1}%", pc(cb_fail_r2));
    println!("  EMITTED (full)               {emit_ok:>10}  {:.1}%", pc(emit_ok));
    println!("\n=== CB-fail breakdown ({cb_fail_r2}) ===");
    let f = |x: u64| 100.0 * x as f64 / cb_fail_r2.max(1) as f64;
    println!("  located-confident (both flanks) {fail_located:>10}  {:.1}%  <- genuine base errors", f(fail_located));
    println!("  R1-tail confirmed            {tail_conf_among_fail:>10}  {:.1}%", f(tail_conf_among_fail));
    println!("  unique within-2 (recoverable){d2_unique:>10}  {:.1}%  (of which located: {fail_located_d2unique})", f(d2_unique));
    println!("  ambiguous within-2           {d2_ambig:>10}  {:.1}%", f(d2_ambig));
    println!("  far (>2 edits / garbage)     {far:>10}  {:.1}%", f(far));
    println!("  span delta histogram (clamped +-5): {span_delta:?}");
    Ok(())
}

enum Rec { Unique, Ambiguous, None }

/// Is there a unique whitelist barcode within Hamming-2 of `cb`?
fn unique_within2(wl: &Whitelist, cb: &[u8]) -> Rec {
    let mut hit: Option<Vec<u8>> = None;
    let mut multiple = false;
    let bases = b"ACGT";
    let mut buf = cb.to_vec();
    // exact or distance-1 first (cheap), then distance-2
    for i in 0..16 {
        let oi = buf[i];
        for &bi in bases {
            if bi == oi { continue; }
            buf[i] = bi;
            check(wl, &buf, &mut hit, &mut multiple);
            for j in (i + 1)..16 {
                let oj = buf[j];
                for &bj in bases {
                    if bj == oj { continue; }
                    buf[j] = bj;
                    check(wl, &buf, &mut hit, &mut multiple);
                }
                buf[j] = oj;
            }
        }
        buf[i] = oi;
    }
    if multiple { Rec::Ambiguous } else if hit.is_some() { Rec::Unique } else { Rec::None }
}

fn check(wl: &Whitelist, cand: &[u8], hit: &mut Option<Vec<u8>>, multiple: &mut bool) {
    if *multiple || !wl.contains(cand) { return; }
    match hit {
        None => *hit = Some(cand.to_vec()),
        Some(h) if h != cand => *multiple = true,
        _ => {}
    }
}
