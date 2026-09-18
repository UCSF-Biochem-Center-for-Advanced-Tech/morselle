//! Dev tool: characterize the CB single-substitution error by POSITION and TYPE.
//!
//! Among confidently-located CBs that are exactly Hamming-1 from a unique whitelist
//! barcode, the true base is known (the whitelist neighbor) and the observed base is
//! the read's — so each is a clean (position, truth→observed) substitution event.
//! Tabulates the position histogram and, for the hotspot, whether the errors are a
//! SPECIFIC substitution (systematic miscall → potentially correctable) or spread
//! across all types (uniform noise). Also prints whitelist base composition per
//! position (is the hotspot a low-complexity region of the barcode design?).
//!
//! Usage: cargo run --release --example cb_error_profile -- <fastq> <whitelist> <p5-led|p7-led> [max]

use std::env;

use needletail::parse_fastx_file;
use sbx10x::barcode::{decode16, Whitelist};
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::Orientation;

fn bi(b: u8) -> usize {
    match b {
        b'A' => 0,
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 4,
    }
}
const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn main() -> anyhow::Result<()> {
    let mut a = env::args().skip(1);
    let fq = a
        .next()
        .expect("usage: cb_error_profile <fastq> <whitelist> <p5-led|p7-led> [max]");
    let wlp = a.next().expect("whitelist");
    let orient_arg = a.next().expect("orientation");
    let max: usize = a.next().map(|s| s.parse().unwrap()).unwrap_or(usize::MAX);

    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
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
    };

    // whitelist base composition per position
    let mut wl_comp = [[0u64; 4]; 16];
    for code in wl.codes_sorted() {
        let bc = decode16(code);
        for i in 0..16 {
            wl_comp[i][bi(bc[i])] += 1;
        }
    }

    // subst[pos][truth][obs]
    let mut subst = [[[0u64; 4]; 4]; 16];
    let mut pos_total = [0u64; 16];
    let mut d1_reads = 0u64;
    // Per-position frame tests (hypothesis: 3' errors are UMI pulled in by an indel).
    let mut span_short = [0u64; 16]; // CB+UMI span < 28 (deletion)
    let mut span_ok = [0u64; 16]; // span == 28
    let mut span_long = [0u64; 16]; // span > 28 (insertion)
    let mut reframe_fix = [0u64; 16]; // a shifted frame yields a unique within-1 hit

    let mut reader = parse_fastx_file(&fq)?;
    let mut seen = 0usize;
    while let Some(rec) = reader.next() {
        if seen >= max {
            break;
        }
        seen += 1;
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        let o = orient(&seq, qual, &ocfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);
        let Some((r1s, _)) = e.r1 else { continue };
        let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
        if !confident || r1s + 16 > o.seq.len() {
            continue;
        }
        let cb = &o.seq[r1s..r1s + 16];
        // Unique within-1 neighbor = the presumed truth. Skip exact (d0).
        let Some(truth) = wl.correct_within1(cb) else {
            continue;
        };
        let mut diff = None;
        let mut ndiff = 0;
        for i in 0..16 {
            if cb[i] != truth[i] {
                ndiff += 1;
                diff = Some(i);
            }
        }
        if ndiff != 1 {
            continue;
        } // clean single-substitution only
        let p = diff.unwrap();
        let (t, ob) = (bi(truth[p]), bi(cb[p]));
        if t < 4 && ob < 4 {
            subst[p][t][ob] += 1;
            pos_total[p] += 1;
            d1_reads += 1;
            // Frame tests: measured CB+UMI span, and whether a shifted frame gives
            // a unique within-1 whitelist hit (⇒ the "d1" is really a misalignment).
            if let Some(pts) = o.polyt_start {
                match (pts as i64) - (r1s as i64) - 28 {
                    d if d < 0 => span_short[p] += 1,
                    0 => span_ok[p] += 1,
                    _ => span_long[p] += 1,
                }
            }
            if let Some((_, off)) = wl.reframe_cb(&o.seq, r1s, 2) {
                if off != 0 {
                    reframe_fix[p] += 1;
                }
            }
        }
    }

    let avg = d1_reads as f64 / 16.0;
    // transition = A<->G or C<->T
    let is_transition = |t: usize, o: usize| {
        (t == 0 && o == 2) || (t == 2 && o == 0) || (t == 1 && o == 3) || (t == 3 && o == 1)
    };
    println!(
        "clean single-substitution (d1, confident) events: {d1_reads}  (avg {:.0}/pos)",
        avg
    );
    println!("frame tests per position: span!=28 ⇒ CB±UMI indel; reframe-fixable ⇒ a shifted frame is a unique whitelist hit (⇒ misaligned, not a substitution)\n");
    println!(
        "{:>4} {:>7} {:>5}  {:>6}  {:>8}  {:>8}  {:>10}",
        "pos", "errors", "xavg", "trans%", "span<28%", "span>28%", "reframe-fix%"
    );
    println!("{}", "-".repeat(60));
    for p in 0..16 {
        let tot = pos_total[p].max(1);
        let mut trans = 0u64;
        for t in 0..4 {
            for o in 0..4 {
                if is_transition(t, o) {
                    trans += subst[p][t][o];
                }
            }
        }
        let mark = if pos_total[p] as f64 > 1.8 * avg {
            " <=="
        } else {
            ""
        };
        println!(
            "{:>4} {:>7} {:>4.1}x  {:>5.0}%  {:>7.1}%  {:>7.1}%  {:>9.1}%{}",
            p,
            pos_total[p],
            pos_total[p] as f64 / avg.max(1.0),
            100.0 * trans as f64 / tot as f64,
            100.0 * span_short[p] as f64 / tot as f64,
            100.0 * span_long[p] as f64 / tot as f64,
            100.0 * reframe_fix[p] as f64 / tot as f64,
            mark
        );
    }

    // Detail the hotspot positions: full substitution matrix + wl base composition.
    println!("\nHotspot detail (truth→obs counts) + whitelist base composition:");
    for p in [7usize, 8, 9] {
        println!("  pos {p}:");
        print!("    subst  ");
        for t in 0..4 {
            for ob in 0..4 {
                if t != ob && subst[p][t][ob] > 0 {
                    print!(
                        "{}→{}:{}  ",
                        BASES[t] as char, BASES[ob] as char, subst[p][t][ob]
                    );
                }
            }
        }
        let tot: u64 = wl_comp[p].iter().sum();
        print!("\n    wl comp  ");
        for b in 0..4 {
            print!(
                "{}:{:.0}%  ",
                BASES[b] as char,
                100.0 * wl_comp[p][b] as f64 / tot as f64
            );
        }
        println!();
    }
    // baseline position for comparison
    println!("  (baseline) pos 3 wl comp: ");
    let tot: u64 = wl_comp[3].iter().sum();
    print!("    ");
    for b in 0..4 {
        print!(
            "{}:{:.0}%  ",
            BASES[b] as char,
            100.0 * wl_comp[3][b] as f64 / tot as f64
        );
    }
    println!();
    Ok(())
}
