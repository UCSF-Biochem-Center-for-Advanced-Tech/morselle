//! Dev tool: the demux upper bound, decoupled from the CB *gate*.
//! For every structurally-located read we assign the sample from i5/i7
//! REGARDLESS of whether the cell barcode passes the whitelist, then cross-
//! tabulate with CB success. Answers: what fraction of reads are sample-
//! identifiable at all, and of those how often do we still lose the CB
//! (the actionable gap)? Plus: reads we can't structurally locate but whose
//! distal i7 alone still IDs the sample (extra headroom).
//! Usage: cargo run --release --example index_upper_bound -- <fastq> <wl> <sheet> <p5-led|p7-led> [max]

use std::env;

use morselle::anchors::{hamming, ILLUMINA_3P_ADAPTER};
use morselle::barcode::Whitelist;
use morselle::demux::{index_regions, DemuxCfg, Demuxer};
use morselle::extract::{extract, ExtractCfg};
use morselle::orient::{orient, FallbackMode, OrientCfg};
use morselle::pipeline::{CbStatus, Orientation};
use morselle::sample_sheet::SampleSheet;
use needletail::parse_fastx_file;

/// Independent i7: find the 3' adapter anywhere, i7 = +rc(R2 primer) (34 bp).
fn indep_i7<'a>(seq: &'a [u8], ilen: usize) -> Option<&'a [u8]> {
    let m = ILLUMINA_3P_ADAPTER.len();
    if seq.len() < m {
        return None;
    }
    let mut best: Option<(u32, usize)> = None;
    for i in 0..=seq.len() - m {
        let d = hamming(&seq[i..i + m], ILLUMINA_3P_ADAPTER);
        if d <= 2 && best.map_or(true, |(b, _)| d < b) {
            best = Some((d, i));
            if d == 0 {
                break;
            }
        }
    }
    best.and_then(|(_, i)| {
        let s = i + 34;
        (s + ilen <= seq.len()).then(|| &seq[s..s + ilen])
    })
}

fn cb_ok(o: &morselle::orient::Oriented, e: &morselle::extract::Extracted, wl: &Whitelist) -> bool {
    let Some((r1s, _)) = e.r1 else { return false };
    let cb = &o.seq[r1s..r1s + 16];
    if !matches!(wl.classify(cb), CbStatus::Fail) {
        return true;
    }
    if let Some(pts) = o.polyt_start {
        let delta = (pts as i64 - r1s as i64) - 28;
        if (delta == 1 || delta == -1)
            && pts >= 12
            && pts <= o.seq.len()
            && wl.realign_cb(&o.seq[r1s..], delta as i32).is_some()
        {
            return true;
        }
    }
    let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
    confident && wl.correct_within2(cb).is_some()
}

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a
        .next()
        .expect("usage: index_upper_bound <fastq> <wl> <sheet> <orient> [max]");
    let wl = Whitelist::load(std::path::Path::new(&a.next().expect("wl")))?;
    let sheet = SampleSheet::load(std::path::Path::new(&a.next().expect("sheet")))?;
    let orient_arg = a.next().expect("orient");
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let demux = Demuxer::new(
        &sheet,
        DemuxCfg {
            index_max_mm: 1,
            mode: None,
        },
    )?;
    let ilen = demux.index_len();
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

    let mut total = 0u64;
    let (mut located, mut has_i5, mut has_i7) = (0u64, 0u64, 0u64);
    let (mut assignable, mut a_cb_ok, mut a_cb_lost) = (0u64, 0u64, 0u64);
    let (mut unlocated, mut unloc_i7_id) = (0u64, 0u64);

    let mut reader = parse_fastx_file(&fq)?;
    while let Some(rec) = reader.next() {
        if total as usize >= max {
            break;
        }
        let rec = rec?;
        let qual = rec.qual().unwrap_or(&[]);
        total += 1;
        let raw = rec.seq();
        let o = orient(&raw, qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg, false);

        let Some((r1s, _)) = e.r1 else {
            unlocated += 1;
            if let Some(i7) = indep_i7(&o.seq, ilen) {
                if demux.assign(None, Some(i7)).sample.is_some() {
                    unloc_i7_id += 1;
                }
            }
            continue;
        };
        located += 1;
        let (i5, i7) = index_regions(&o.seq, r1s, e.adapter_pos, ilen);
        if i5.is_some() {
            has_i5 += 1;
        }
        if i7.is_some() {
            has_i7 += 1;
        }
        if demux.assign(i5, i7).sample.is_some() {
            assignable += 1;
            if cb_ok(&o, &e, &wl) {
                a_cb_ok += 1;
            } else {
                a_cb_lost += 1;
            }
        }
    }

    let pc = |x: u64| 100.0 * x as f64 / total as f64;
    println!("reads {total}");
    println!(
        "\nStructurally located (CB bracket found): {located} ({:.1}%)",
        pc(located)
    );
    println!(
        "  i5 present {has_i5} ({:.1}%) · i7 present {has_i7} ({:.1}%)",
        pc(has_i5),
        pc(has_i7)
    );
    println!(
        "\nSAMPLE-ASSIGNABLE regardless of CB gate: {assignable} ({:.1}%)  <- upper bound",
        pc(assignable)
    );
    let d = assignable.max(1);
    println!(
        "  CB also recovered                {a_cb_ok:>8}  {:.1}% of assignable",
        100.0 * a_cb_ok as f64 / d as f64
    );
    println!("  CB LOST despite known sample     {a_cb_lost:>8}  {:.1}% of assignable  (= {:.1}% of all reads)  <- ACTIONABLE GAP",
        100.0 * a_cb_lost as f64 / d as f64, pc(a_cb_lost));
    println!("\nUnlocated reads {unlocated} ({:.1}%); of these, distal i7 alone still IDs the sample: {unloc_i7_id} ({:.1}% of all reads extra headroom)",
        pc(unlocated), pc(unloc_i7_id));
    Ok(())
}
