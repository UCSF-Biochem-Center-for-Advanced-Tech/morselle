//! Dev tool: extract the ORIGINAL raw reads that fall into specific loss buckets,
//! for manual inspection. Mirrors the pipeline classifier (orient enforced →
//! extract → realign → CB reframe → lenient → gate). Writes the untouched input
//! record (id/seq/qual, original orientation) — NOT the oriented/trimmed version.
//!
//! Buckets written:
//!   unparseable  — no CB bracket located (e.r1 is None)
//!   cbfail       — R1+R2 present but CB fails the whitelist after all rescues
//!
//! Usage: cargo run --release --example extract_loss -- \
//!   <fastq> <whitelist> <sheet> <p5-led|p7-led> <out_unparseable.fastq.gz> <out_cbfail.fastq.gz> [max]

use std::env;
use std::fs::File;
use std::io::{BufWriter, Write};

use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use sbx10x::barcode::Whitelist;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};

type GzW = GzEncoder<BufWriter<File>>;

fn writer(path: &str) -> GzW {
    GzEncoder::new(BufWriter::new(File::create(path).unwrap()), Compression::default())
}

fn write_rec(w: &mut GzW, id: &[u8], seq: &[u8], qual: &[u8]) {
    w.write_all(b"@").unwrap();
    w.write_all(id).unwrap();
    w.write_all(b"\n").unwrap();
    w.write_all(seq).unwrap();
    w.write_all(b"\n+\n").unwrap();
    w.write_all(qual).unwrap();
    w.write_all(b"\n").unwrap();
}

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a.next().unwrap();
    let wl = Whitelist::load(std::path::Path::new(&a.next().unwrap()))?;
    let _sheet = a.next().unwrap(); // accepted for symmetry; CB buckets don't need it
    let orient_arg = a.next().unwrap();
    let out_unp = a.next().unwrap();
    let out_cbf = a.next().unwrap();
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let (cb_len, umi_len) = (16usize, 12usize);
    let expected = if orient_arg == "p7-led" { Orientation::P7Led } else { Orientation::P5Led };
    let ocfg = OrientCfg { cb_len, umi_len, min_homopolymer: 12, r1_tail_max_mm: 2,
        anchor_max_edits: 6, fallback: FallbackMode::HammingTail, expected: Some(expected) };
    let ecfg = ExtractCfg { min_cdna_len: 50, adapter_max_mm: 2 };

    let mut w_unp = writer(&out_unp);
    let mut w_cbf = writer(&out_cbf);
    let (mut n, mut n_unp, mut n_cbf) = (0usize, 0u64, 0u64);

    let mut reader = parse_fastx_file(&fq)?;
    while let Some(rec) = reader.next() {
        if n >= max { break; }
        n += 1;
        let rec = rec?;
        let id = rec.id().to_vec();
        let raw_seq = rec.seq().to_vec(); // ORIGINAL, untouched
        let raw_qual = rec.qual().unwrap_or(&[]).to_vec();

        let o = orient(&raw_seq, &raw_qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);

        let Some((r1s, _)) = e.r1 else {
            // Unparseable: no CB bracket located.
            write_rec(&mut w_unp, &id, &raw_seq, &raw_qual);
            n_unp += 1;
            continue;
        };
        if e.r2.is_none() {
            continue; // no-cDNA bucket — not requested
        }
        // CB gate with the full rescue stack (matches the pipeline).
        let cb = &o.seq[r1s..r1s + cb_len];
        let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
        let mut cb_ok = !matches!(wl.classify(cb), CbStatus::Fail);
        if !cb_ok {
            if let Some(pts) = o.polyt_start {
                let delta = (pts as i64 - r1s as i64) - (cb_len + umi_len) as i64;
                if (delta == 1 || delta == -1) && pts >= umi_len && pts <= o.seq.len()
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
            // CB fails the whitelist after realign/reframe/lenient.
            write_rec(&mut w_cbf, &id, &raw_seq, &raw_qual);
            n_cbf += 1;
        }
    }

    w_unp.finish()?.flush()?;
    w_cbf.finish()?.flush()?;
    let pct = |x: u64| 100.0 * x as f64 / n.max(1) as f64;
    println!("reads {n}");
    println!("  unparseable -> {out_unp}: {n_unp} ({:.1}%)", pct(n_unp));
    println!("  cbfail      -> {out_cbf}: {n_cbf} ({:.1}%)", pct(n_cbf));
    Ok(())
}
