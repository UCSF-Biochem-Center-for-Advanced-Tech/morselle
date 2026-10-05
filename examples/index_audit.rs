//! Dev tool: why do otherwise-valid reads (CB ok, R2 ok) fail sample assignment?
//! For undetermined reads, characterise the extracted i5/i7 — presence, distance
//! to the nearest sample index, and the most frequent unmatched sequences
//! (recurring = foreign/un-listed sample; diverse = degradation).
//! Usage: cargo run --release --example index_audit -- <fastq> <wl> <sheet> <p5-led|p7-led> [max]

use std::collections::HashMap;
use std::env;

use needletail::parse_fastx_file;
use sbx10x::anchors::revcomp;
use sbx10x::barcode::Whitelist;
use sbx10x::demux::{index_regions, DemuxCfg, Demuxer};
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};
use sbx10x::sample_sheet::SampleSheet;

fn ham(a: &[u8], b: &[u8]) -> u32 {
    a.iter().zip(b).filter(|(x, y)| x != y).count() as u32 + (a.len().abs_diff(b.len())) as u32
}

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a.next().unwrap();
    let wl = Whitelist::load(std::path::Path::new(&a.next().unwrap()))?;
    let sheet = SampleSheet::load(std::path::Path::new(&a.next().unwrap()))?;
    let orient_arg = a.next().unwrap();
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    // read-oriented sample indices: i5 = revcomp(sheet i5), i7 = sheet i7
    let i5_ref: Vec<Vec<u8>> = sheet.samples.iter().map(|s| revcomp(&s.i5)).collect();
    let i7_ref: Vec<Vec<u8>> = sheet.samples.iter().map(|s| s.i7.clone()).collect();
    let nearest = |r: &[u8], refs: &[Vec<u8>]| refs.iter().map(|x| ham(r, x)).min().unwrap_or(99);

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

    let mut undet = 0u64;
    let (mut both, mut i5only, mut i7only, mut neither) = (0u64, 0u64, 0u64, 0u64);
    let mut i5_dist = [0u64; 6]; // distance of present i5 to nearest sample (capped 5)
    let mut i5_freq: HashMap<Vec<u8>, u64> = HashMap::new();
    let mut reframe_ok = 0u64;
    let mut off_hist: std::collections::BTreeMap<i64, u64> = std::collections::BTreeMap::new();

    let mut reader = parse_fastx_file(&fq)?;
    let mut n = 0;
    while let Some(rec) = reader.next() {
        if n >= max {
            break;
        }
        let rec = rec?;
        let qual = rec.qual().unwrap_or(&[]);
        n += 1;
        let raw = rec.seq();
        let o = orient(&raw, qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg, false);
        let Some((r1s, _)) = e.r1 else { continue };
        if e.r2.is_none() {
            continue;
        }
        if matches!(wl.classify(&o.seq[r1s..r1s + 16]), CbStatus::Fail) {
            continue;
        } // CB must be ok
        let (i5, i7) = index_regions(&o.seq, r1s, e.adapter_pos, ilen);
        if demux.assign(i5, i7).sample.is_some() {
            continue;
        } // only undetermined
        undet += 1;
        match (i5.is_some(), i7.is_some()) {
            (true, true) => both += 1,
            (true, false) => i5only += 1,
            (false, true) => i7only += 1,
            (false, false) => neither += 1,
        }
        if let Some(r) = i5 {
            i5_dist[nearest(r, &i5_ref).min(5) as usize] += 1;
            *i5_freq.entry(r.to_vec()).or_default() += 1;
            // reframe: try i5 at cb_start-43+k; does a unique sample match appear?
            for k in [1i64, -1, 2, -2, 3, -3] {
                let s = r1s as i64 - 43 + k;
                if s < 0 {
                    continue;
                }
                let s = s as usize;
                if s + ilen > o.seq.len() {
                    continue;
                }
                let cand = &o.seq[s..s + ilen];
                let hits = i5_ref.iter().filter(|x| ham(cand, x) <= 1).count();
                if hits == 1 {
                    reframe_ok += 1;
                    *off_hist.entry(k).or_default() += 1;
                    break;
                }
            }
        }
        let _ = (i7, &i7_ref);
    }

    let pu = |x: u64| 100.0 * x as f64 / undet.max(1) as f64;
    println!("reads {n};  undetermined (CB ok, R2 ok, no sample): {undet}");
    println!("  index presence: both {} ({:.0}%) · i5-only {} ({:.0}%) · i7-only {} ({:.0}%) · neither {} ({:.0}%)",
        both, pu(both), i5only, pu(i5only), i7only, pu(i7only), neither, pu(neither));
    println!("  i5 distance to nearest sample index:");
    for d in 0..6 {
        let lbl = if d == 5 { "5+".into() } else { d.to_string() };
        println!(
            "    d{lbl}: {} ({:.1}%)",
            i5_dist[d],
            100.0 * i5_dist[d] as f64 / i5_dist.iter().sum::<u64>().max(1) as f64
        );
    }
    let mut top: Vec<_> = i5_freq.into_iter().collect();
    top.sort_by(|a, b| b.1.cmp(&a.1));
    println!("  top unmatched i5 sequences (read-oriented) — recurring = foreign sample:");
    for (seq, c) in top.iter().take(12) {
        let d = nearest(seq, &i5_ref);
        println!(
            "    {}  x{:<6}  (nearest sample i5: d{})",
            String::from_utf8_lossy(seq),
            c,
            d
        );
    }
    println!("\n  REFRAME: undetermined reads recovered by shifting the i5 window: {} ({:.1}% of undetermined)",
        reframe_ok, 100.0 * reframe_ok as f64 / undet.max(1) as f64);
    println!("  winning offset histogram: {off_hist:?}");
    Ok(())
}
