//! Dev tool: complete per-read disposition — every read lands in exactly one
//! terminal bucket (sums to 100%), attributed to the stage + reason it's lost.
//! Mirrors the full pipeline: orient(enforced) → extract → realign → lenient →
//! gate → demux.
//! Usage: cargo run --release --example disposition -- <fastq> <whitelist> <sheet> <p5-led|p7-led> [max]

use std::env;

use needletail::parse_fastx_file;
use sbx10x::barcode::Whitelist;
use sbx10x::demux::{index_anchor_starts, DemuxCfg, Demuxer};
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};
use sbx10x::sample_sheet::SampleSheet;

#[derive(Clone, Copy)]
enum Bucket {
    Unparseable,        // no CB bracket located (orient/anchor fail)
    NoCdna,             // CB located but no cDNA (R2 absent/short)
    CbFailConfident,    // R1+R2, CB confidently located but degraded (>=2 err)
    CbFailUncertain,    // R1+R2, CB not confidently located (no/false barcode)
    UndetNoIndex,       // valid pair but neither i5/i7 region present
    UndetNoMatch,       // valid pair, index present but unmatched/conflict
    Assigned,           // emitted to a sample
}

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a.next().expect("usage: disposition <fastq> <whitelist> <sheet> <p5-led|p7-led> [max]");
    let wlp = a.next().expect("whitelist");
    let ssp = a.next().expect("sheet");
    let orient_arg = a.next().expect("orientation p5-led|p7-led");
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
    let sheet = SampleSheet::load(std::path::Path::new(&ssp))?;
    let demux = Demuxer::new(&sheet, DemuxCfg { index_max_mm: 1, allow_single_index: true });
    let expected = match orient_arg.as_str() {
        "p7-led" => Orientation::P7Led,
        _ => Orientation::P5Led,
    };
    let (cb_len, umi_len) = (16usize, 12usize);
    let ocfg = OrientCfg {
        cb_len, umi_len, min_homopolymer: 12, r1_tail_max_mm: 2,
        anchor_max_edits: 6, fallback: FallbackMode::HammingTail, expected: Some(expected),
    };
    let ecfg = ExtractCfg { min_cdna_len: 50, adapter_max_mm: 2, max_r2_len: 0 };

    let mut counts = [0u64; 7];
    let mut total = 0u64;
    let mut reader = parse_fastx_file(&fq)?;
    while let Some(rec) = reader.next() {
        if total as usize >= max { break; }
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        total += 1;

        let o = orient(&seq, qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);

        let bucket = classify(&o, &e, &wl, &demux, cb_len, umi_len);
        counts[bucket as usize] += 1;
    }

    let names = [
        ("Unparseable — no CB bracket located", "no dT-primer + no R1-tail; read lacks barcode structure (truncated/junk)"),
        ("No cDNA — R2 absent/too short", "barcode region present but read ends before any transcript"),
        ("CB gated — confident, degraded", "CB correctly located but >=2 sequencing errors (genuinely unreadable)"),
        ("CB gated — uncertain location", "no confident anchor; likely no real barcode (truncated/false frame)"),
        ("Undetermined — no index", "valid pair but the i5/i7 region wasn't reached (3' truncation)"),
        ("Undetermined — index unmatched", "index present but no sample match / i5-i7 conflict"),
        ("ASSIGNED to a sample", "recovered"),
    ];
    println!("reads {total}\n");
    println!("{:<42} {:>10} {:>7}", "disposition", "reads", "%");
    let lost: u64 = (0..6).map(|i| counts[i]).sum();
    for i in [6usize, 0, 3, 2, 1, 4, 5] {
        let pct = 100.0 * counts[i] as f64 / total as f64;
        println!("{:<42} {:>10} {:>6.1}%", names[i].0, counts[i], pct);
        println!("    └ {}", names[i].1);
    }
    println!("\nTOTAL LOST (not in a sample): {} ({:.1}%)", lost, 100.0 * lost as f64 / total as f64);
    Ok(())
}

fn classify(
    o: &sbx10x::orient::Oriented,
    e: &sbx10x::extract::Extracted,
    wl: &Whitelist,
    demux: &Demuxer,
    cb_len: usize,
    umi_len: usize,
) -> Bucket {
    let Some((r1s, _)) = e.r1 else { return Bucket::Unparseable };
    if e.r2.is_none() {
        return Bucket::NoCdna;
    }
    // CB gate with realign + CB reframe + lenient rescue (matches the pipeline).
    let cb = &o.seq[r1s..r1s + cb_len];
    let status = wl.classify(cb);
    let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
    let mut cb_ok = !matches!(status, CbStatus::Fail);
    if !cb_ok {
        // realign (single indel)
        if let Some(pts) = o.polyt_start {
            let delta = (pts as i64 - r1s as i64) - (cb_len + umi_len) as i64;
            if (delta == 1 || delta == -1) && pts >= umi_len && pts <= o.seq.len()
                && wl.realign_cb(&o.seq[r1s..], delta as i32).is_some()
            {
                cb_ok = true;
            }
        }
        // CB reframe (homopolymer-indel frame shift), confident only
        if !cb_ok && confident && wl.reframe_cb(&o.seq, r1s, 2).is_some() {
            cb_ok = true;
        }
        // lenient within-2, confident only
        if !cb_ok && confident && wl.correct_within2(cb).is_some() {
            cb_ok = true;
        }
    }
    if !cb_ok {
        return if confident { Bucket::CbFailConfident } else { Bucket::CbFailUncertain };
    }
    // Demux (offset-tolerant, matches the pipeline's --index-window 3)
    let (i5s, i7s) = index_anchor_starts(r1s, e.adapter_pos, demux.index_len());
    let outcome = demux.assign_windowed(&o.seq, i5s, i7s, 3);
    if outcome.sample.is_some() {
        Bucket::Assigned
    } else if i5s.is_none() && i7s.is_none() {
        Bucket::UndetNoIndex
    } else {
        Bucket::UndetNoMatch
    }
}
