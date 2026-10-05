//! Dev tool: is the "uncertain" CB-fail group caused by false homopolymer
//! anchoring (an internal poly-A/T in cDNA) rather than absent barcodes?
//! For each read we independently scan the canonical sequence for the R1-primer
//! tail with a VALID downstream CB — i.e. the real barcode structure — and see
//! whether it exists at a position other than where we anchored.
//! Usage: cargo run --release --example false_anchor -- <fastq> <whitelist> [max]

use std::env;

use needletail::parse_fastx_file;
use sbx10x::anchors::{hamming, R1_PRIMER_TAIL};
use sbx10x::barcode::Whitelist;
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};
use sbx10x::pipeline::{CbStatus, Orientation};

/// Find an R1-primer tail (≤1 mm) anywhere in `seq` whose downstream 16 bp CB is
/// on/near the whitelist. Returns the CB start (tail end) of the first such hit.
fn real_structure(seq: &[u8], wl: &Whitelist) -> Option<usize> {
    let m = R1_PRIMER_TAIL.len();
    if seq.len() < m + 16 {
        return None;
    }
    for i in 0..=seq.len() - m - 16 {
        if hamming(&seq[i..i + m], R1_PRIMER_TAIL) <= 1 {
            let cb = &seq[i + m..i + m + 16];
            if !matches!(wl.classify(cb), CbStatus::Fail) {
                return Some(i + m);
            }
        }
    }
    None
}

fn main() -> anyhow::Result<()> {
    let mut args = env::args().skip(1);
    let fq = args
        .next()
        .expect("usage: false_anchor <fastq> <whitelist> [max]");
    let wlp = args.next().expect("need whitelist");
    let max: usize = args.next().map(|s| s.parse().unwrap()).unwrap_or(40000);

    let wl = Whitelist::load(std::path::Path::new(&wlp))?;
    let ocfg = OrientCfg {
        cb_len: 16,
        umi_len: 12,
        min_homopolymer: 12,
        r1_tail_max_mm: 2,
        anchor_max_edits: 6,
        fallback: FallbackMode::HammingTail,
        expected: None,
    };
    let ecfg = ExtractCfg {
        min_cdna_len: 20,
        adapter_max_mm: 2,
        trim_tso: false,
        max_r2_len: 0,
        flex_r2_boundary: false,
        flex_min_cdna_len_no_anchor: 65,
    };

    let mut unc_fail = 0u64; // uncertain + CB-fail reads
    let mut unc_has_real = 0u64; // ... with real structure elsewhere
    let mut unc_real_diff_pos = 0u64; // ... at a different position than we anchored
    let mut unc_len_sum = 0u64;
    let mut anchored_depth_sum = 0u64; // where we put cb_start
    let mut real_depth_sum = 0u64; // where the real CB actually is
                                   // control: confident + CB-fail
    let (mut conf_fail, mut conf_has_real) = (0u64, 0u64);

    let mut reader = parse_fastx_file(&fq)?;
    let mut n = 0;
    while let Some(rec) = reader.next() {
        if n >= max {
            break;
        }
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        n += 1;
        let o = orient(&seq, qual, &ocfg);
        if o.orientation == Orientation::Unknown {
            continue;
        }
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg, false);
        let Some((r1s, _)) = e.r1 else {
            continue;
        };
        if !matches!(wl.classify(&o.seq[r1s..r1s + 16]), CbStatus::Fail) {
            continue; // only CB-fails
        }
        let confident = o.r1_tail_confirmed && e.polyt_run_len >= 10;
        let real = real_structure(&o.seq, &wl);
        if confident {
            conf_fail += 1;
            if real.is_some() {
                conf_has_real += 1;
            }
        } else {
            unc_fail += 1;
            unc_len_sum += o.seq.len() as u64;
            anchored_depth_sum += r1s as u64;
            if let Some(rp) = real {
                unc_has_real += 1;
                real_depth_sum += rp as u64;
                if rp.abs_diff(r1s) > 2 {
                    unc_real_diff_pos += 1;
                }
            }
        }
    }

    let pu = |x: u64| 100.0 * x as f64 / unc_fail.max(1) as f64;
    println!("reads scanned {n}");
    println!("\nUNCERTAIN + CB-fail: {unc_fail}");
    println!(
        "  has real structure elsewhere (R1-tail + valid CB)  {unc_has_real}  ({:.1}%)",
        pu(unc_has_real)
    );
    println!(
        "    of those, at a DIFFERENT position (mis-anchored)  {unc_real_diff_pos}  ({:.1}%)",
        pu(unc_real_diff_pos)
    );
    println!(
        "  mean read length                {}",
        unc_len_sum / unc_fail.max(1)
    );
    println!(
        "  mean anchored cb_start (ours)   {}",
        anchored_depth_sum / unc_fail.max(1)
    );
    if unc_has_real > 0 {
        println!(
            "  mean TRUE cb_start (found)      {}",
            real_depth_sum / unc_has_real
        );
    }
    println!("\nCONFIDENT + CB-fail (control): {conf_fail}");
    println!(
        "  has real structure elsewhere    {conf_has_real}  ({:.1}%)",
        100.0 * conf_has_real as f64 / conf_fail.max(1) as f64
    );
    Ok(())
}
