//! Dev tool: real-data precision proxy for the single-indel realign rescue.
//!
//! Roche gives only the SET of cell barcodes (Cell Ranger matrix), not per-read CB
//! calls, so a literal per-read "vs Roche" check isn't possible. Instead: build the
//! trusted CB abundance map from confident exact / within-1 reads (the real cells),
//! then ask where realign-recovered CBs land. Correct recoveries should pile onto
//! ABUNDANT trusted barcodes (real cells); false ones scatter to rare/absent
//! whitelist entries. Compare against a random-whitelist null to quantify the lift.
//!
//! Usage: cargo run --release --example realign_concordance -- <fastq> <whitelist> <p5-led|p7-led> [max]

use std::collections::HashMap;
use std::env;

use needletail::parse_fastx_file;
use sbx10x::barcode::{encode16, Whitelist};
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};

struct Rng(u64);
impl Rng { fn next(&mut self) -> u64 { self.0 = self.0.wrapping_add(0x9E3779B97F4A7C15); let mut z=self.0; z=(z^(z>>30)).wrapping_mul(0xBF58476D1CE4E5B9); z=(z^(z>>27)).wrapping_mul(0x94D049BB133111EB); z^(z>>31) } }

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a.next().expect("usage: realign_concordance <fastq> <whitelist> <p5-led|p7-led> [max]");
    let wlp = a.next().expect("whitelist");
    let orient_arg = a.next().expect("orientation");
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
    let expected = if orient_arg == "p7-led" { Orientation::P7Led } else { Orientation::P5Led };
    let ocfg = OrientCfg { cb_len: 16, umi_len: 12, min_homopolymer: 12, r1_tail_max_mm: 2, anchor_max_edits: 6, fallback: FallbackMode::HammingTail, expected: Some(expected) };
    let ecfg = ExtractCfg { min_cdna_len: 50, adapter_max_mm: 2, max_r2_len: 0 };

    let mut trusted: HashMap<u32, u32> = HashMap::new(); // exact/within-1 CB → count
    let mut recovered: Vec<u32> = Vec::new();            // realign-recovered CB codes

    let mut reader = parse_fastx_file(&fq)?;
    let mut seen = 0usize;
    while let Some(rec) = reader.next() {
        if seen >= max { break; }
        seen += 1;
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        let o = orient(&seq, qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);
        let Some((r1s, _)) = e.r1 else { continue };
        if r1s + 16 > o.seq.len() { continue; }
        let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
        let cb = &o.seq[r1s..r1s + 16];
        match wl.classify(cb) {
            CbStatus::OnWhitelist => { if let Some(c) = encode16(cb) { *trusted.entry(c).or_insert(0) += 1; } }
            CbStatus::Within1 => { if let Some(bc) = wl.correct_within1(cb) { if let Some(c) = encode16(&bc) { *trusted.entry(c).or_insert(0) += 1; } } }
            CbStatus::Fail => {
                if !confident { continue; }
                if let Some(pts) = o.polyt_start {
                    let delta = pts as i64 - r1s as i64 - 28;
                    if (delta == 1 || delta == -1) && pts >= 12 {
                        if let Some(bc) = wl.realign_cb(&o.seq[r1s..], delta as i32) {
                            if let Some(c) = encode16(&bc) { recovered.push(c); }
                        }
                    }
                }
            }
            CbStatus::Indel => {}
        }
    }

    let thresholds = [1u32, 5, 20];
    println!("trusted (exact/within-1) distinct CBs: {}  ({} reads total)", trusted.len(),
        trusted.values().map(|&v| v as u64).sum::<u64>());
    println!("realign-recovered reads: {}\n", recovered.len());
    println!("Where recovered CBs land vs a random-whitelist null (higher lift = more precise):");
    println!("{:>16} {:>14} {:>14} {:>8}", "trusted count >=", "recovered%", "random null%", "lift");

    let mut rng = Rng(0x1357_9BDF);
    let codes = wl.codes_sorted();
    let nrand = 200_000usize;
    for t in thresholds {
        let rec_hit = recovered.iter().filter(|c| trusted.get(c).copied().unwrap_or(0) >= t).count();
        let rec_pct = 100.0 * rec_hit as f64 / recovered.len().max(1) as f64;
        let mut null_hit = 0usize;
        for _ in 0..nrand {
            let c = codes[(rng.next() as usize) % codes.len()];
            if trusted.get(&c).copied().unwrap_or(0) >= t { null_hit += 1; }
        }
        let null_pct = 100.0 * null_hit as f64 / nrand as f64;
        println!("{:>16} {:>13.1}% {:>13.2}% {:>7.1}x", t, rec_pct, null_pct, rec_pct / null_pct.max(1e-6));
    }
    println!("\nIf recovered CBs land on abundant real-cell barcodes far above the random null,");
    println!("the realign rescues are landing on true cells (high precision), not scattering.");
    Ok(())
}
