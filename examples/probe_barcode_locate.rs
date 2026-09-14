//! Phase 0 of FLEX_SUPPORT_PLAN.md: confirm, with the real orient()/extract()
//! code path (not an ad hoc Python k-mer scan), where the 10x Flex probe
//! barcode and the CRISPR Guide Capture anchor+capture-barcode sit relative to
//! coordinates `extract()` already computes.
//!
//! Usage:
//!   cargo run --release --example probe_barcode_locate -- \
//!     --whitelist <flex-v2-384.txt> --features <feature_reference.csv> \
//!     [--max-per-file N] <fastq> [<fastq> ...]
//!
//! For each confidently-oriented read (same confidence gate as orient_purity):
//!   1. Probe-barcode scan: slide a 10bp window over a neighborhood around
//!      `Extracted::polyt_end` in *canonical* coordinates, both as-is and
//!      revcomp'd, and report the offset/orientation with the best hit rate
//!      against the 384-entry whitelist (<=2 mismatches).
//!   2. CRISPR anchor scan: reconstruct the emitted-R2 orientation (revcomp of
//!      the extract()-bracketed r2 span, uncapped -- identical bytes to what
//!      lib.rs's push_fastq_r2_rc would emit) and slide each feature
//!      reference's anchor over a small window near the start, reporting the
//!      offset with the best hit rate; on a hit, checks the following
//!      `capture_len` bases against each guide's sequence.
//!
//! This is a diagnostic dev tool, not part of the shipped pipeline -- nothing
//! here is hardcoded into `src/` until the offsets it reports are confirmed.

use std::collections::HashMap;
use std::env;
use std::fs;
use std::path::PathBuf;

use needletail::parse_fastx_file;
use sbx10x::anchors::{hamming, revcomp};
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::orient::{orient, FallbackMode, OrientCfg};

struct Whitelist {
    // (forward 10bp, id)
    entries: Vec<(Vec<u8>, String)>,
}

impl Whitelist {
    fn load(path: &str) -> anyhow::Result<Self> {
        let text = fs::read_to_string(path)?;
        let mut entries = Vec::new();
        for line in text.lines() {
            let mut f = line.split_whitespace();
            let seq = f.next().unwrap_or_default();
            let id = f.nth(1).unwrap_or_default(); // 3rd column
            if !seq.is_empty() {
                entries.push((seq.as_bytes().to_vec(), id.to_string()));
            }
        }
        anyhow::ensure!(!entries.is_empty(), "empty whitelist: {path}");
        Ok(Self { entries })
    }

    /// Best (mismatches, id) match of `query` against the whitelist.
    fn best_match(&self, query: &[u8]) -> Option<(u32, &str)> {
        self.entries
            .iter()
            .map(|(seq, id)| (hamming(query, seq), id.as_str()))
            .min_by_key(|(mm, _)| *mm)
    }

    fn min_pairwise_hamming(&self) -> u32 {
        let mut best = u32::MAX;
        for i in 0..self.entries.len() {
            for j in i + 1..self.entries.len() {
                let d = hamming(&self.entries[i].0, &self.entries[j].0);
                if d < best {
                    best = d;
                }
            }
        }
        best
    }
}

struct Feature {
    id: String,
    anchor: Vec<u8>,
    barcode: Vec<u8>,
}

fn load_features(path: &str) -> anyhow::Result<Vec<Feature>> {
    let text = fs::read_to_string(path)?;
    let mut lines = text.lines();
    let header = lines.next().unwrap_or_default();
    let cols: Vec<&str> = header.split(',').collect();
    let idx = |name: &str| cols.iter().position(|c| *c == name).unwrap();
    let (id_i, pattern_i, seq_i) = (idx("id"), idx("pattern"), idx("sequence"));
    let mut out = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split(',').collect();
        let pattern = f[pattern_i];
        anyhow::ensure!(
            pattern.ends_with("(BC)"),
            "unsupported feature pattern (expected <anchor>(BC)): {pattern}"
        );
        let anchor = pattern.trim_end_matches("(BC)").as_bytes().to_vec();
        out.push(Feature {
            id: f[id_i].to_string(),
            anchor,
            barcode: f[seq_i].as_bytes().to_vec(),
        });
    }
    Ok(out)
}

fn main() -> anyhow::Result<()> {
    let mut whitelist_path = None;
    let mut features_path = None;
    let mut max_per_file: usize = 200_000;
    let mut fastqs: Vec<PathBuf> = Vec::new();

    let mut args = env::args().skip(1).peekable();
    while let Some(a) = args.next() {
        match a.as_str() {
            "--whitelist" => whitelist_path = args.next(),
            "--features" => features_path = args.next(),
            "--max-per-file" => max_per_file = args.next().unwrap().parse().unwrap(),
            other => fastqs.push(PathBuf::from(other)),
        }
    }
    let whitelist_path = whitelist_path.expect("--whitelist <flex-v2-384.txt> required");
    let features_path = features_path.expect("--features <feature_reference.csv> required");
    anyhow::ensure!(!fastqs.is_empty(), "usage: probe_barcode_locate --whitelist <FILE> --features <CSV> [--max-per-file N] <fastq> [<fastq> ...]");

    let wl = Whitelist::load(&whitelist_path)?;
    println!(
        "loaded {} probe barcodes from {whitelist_path} (min pairwise Hamming = {})",
        wl.entries.len(),
        wl.min_pairwise_hamming()
    );
    let features = load_features(&features_path)?;
    println!(
        "loaded {} feature(s) from {features_path}: {}",
        features.len(),
        features
            .iter()
            .map(|f| format!("{}(anchor={}bp,bc={}bp)", f.id, f.anchor.len(), f.barcode.len()))
            .collect::<Vec<_>>()
            .join(", ")
    );

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
    };

    const PB_LEN: usize = 10;
    const PB_MAX_MM: u32 = 2;
    const PB_WINDOW: isize = 20; // offsets tried after polyt_end
    const ANCHOR_MAX_MM_FRAC: f64 = 0.15; // ~3/20 for a 20bp anchor
    const ANCHOR_WINDOW: usize = 60; // bp of emitted-R2 5' end searched

    let mut n_reads = 0u64;
    let mut n_confident = 0u64;
    let mut n_has_r2 = 0u64;

    // probe barcode: (offset, revcomp?) -> hit count (<= PB_MAX_MM)
    let mut pb_hits: HashMap<(isize, bool), u64> = HashMap::new();
    // best-offset mismatch histogram at the eventual winning (offset, revcomp)
    let mut pb_id_hits: HashMap<String, u64> = HashMap::new();

    // per-feature: offset -> hit count
    let mut anchor_hits: HashMap<(String, usize), u64> = HashMap::new();
    let mut guide_calls: HashMap<String, u64> = HashMap::new();
    let mut guide_capture_mm_hist: HashMap<u32, u64> = HashMap::new();

    for fq in &fastqs {
        println!("--- {} ---", fq.display());
        let mut reader = parse_fastx_file(fq)?;
        let mut n_file = 0usize;
        while let Some(rec) = reader.next() {
            if n_file >= max_per_file {
                break;
            }
            let rec = rec?;
            let seq = rec.seq();
            let qual = rec.qual().unwrap_or(&[]);
            n_file += 1;
            n_reads += 1;

            let o = orient(&seq, qual, &ocfg);
            let confident = o.r1_tail_confirmed;
            if !confident {
                continue;
            }
            n_confident += 1;
            let e = extract(&o.seq, o.cb_umi, o.polyt_start, &ecfg);
            let Some(polyt_end) = e.polyt_end else {
                continue;
            };
            let canon = o.seq.as_ref();

            // --- 1. probe barcode window scan ---
            for offset in 0..=PB_WINDOW {
                let start = polyt_end as isize + offset;
                if start < 0 {
                    continue;
                }
                let start = start as usize;
                let end = start + PB_LEN;
                if end > canon.len() {
                    break;
                }
                let window = &canon[start..end];
                if let Some((mm, id)) = wl.best_match(window) {
                    if mm <= PB_MAX_MM {
                        *pb_hits.entry((offset, false)).or_default() += 1;
                        if offset == 0 {
                            *pb_id_hits.entry(id.to_string()).or_default() += 1;
                        }
                    }
                }
                let window_rc = revcomp(window);
                if let Some((mm, _id)) = wl.best_match(&window_rc) {
                    if mm <= PB_MAX_MM {
                        *pb_hits.entry((offset, true)).or_default() += 1;
                    }
                }
            }

            // --- 2. CRISPR anchor scan, in emitted-R2 orientation ---
            if let Some((r2s, r2e)) = e.r2 {
                n_has_r2 += 1;
                let emitted_r2 = revcomp(&canon[r2s..r2e]);
                for feat in &features {
                    let alen = feat.anchor.len();
                    let max_mm = ((alen as f64) * ANCHOR_MAX_MM_FRAC).round() as u32;
                    let hay_end = emitted_r2.len().min(ANCHOR_WINDOW);
                    if hay_end < alen {
                        continue;
                    }
                    for offset in 0..=(hay_end - alen) {
                        let window = &emitted_r2[offset..offset + alen];
                        let mm = hamming(window, &feat.anchor);
                        if mm <= max_mm {
                            *anchor_hits.entry((feat.id.clone(), offset)).or_default() += 1;
                            let bc_start = offset + alen;
                            let bc_end = bc_start + feat.barcode.len();
                            if bc_end <= emitted_r2.len() {
                                let bc_window = &emitted_r2[bc_start..bc_end];
                                let bc_mm = hamming(bc_window, &feat.barcode);
                                *guide_capture_mm_hist.entry(bc_mm).or_default() += 1;
                                if bc_mm <= 2 {
                                    *guide_calls.entry(feat.id.clone()).or_default() += 1;
                                }
                            }
                            break; // first (leftmost) anchor hit for this read
                        }
                    }
                }
            }
        }
        println!("  read {n_file} records");
    }

    println!("\n=== summary ===");
    println!("reads scanned:      {n_reads}");
    println!("confidently oriented: {n_confident}");
    println!("with r2 (would-be cDNA/probe span): {n_has_r2}");

    println!("\n--- probe barcode: hit rate by (offset from polyt_end, revcomp) ---");
    let mut pb_rows: Vec<_> = pb_hits.into_iter().collect();
    pb_rows.sort_by_key(|(_, c)| std::cmp::Reverse(*c));
    for ((offset, rc), count) in pb_rows.iter().take(15) {
        println!(
            "  offset={offset:>3} revcomp={rc:<5} hits={count:>10} ({:.2}% of confident)",
            100.0 * *count as f64 / n_confident.max(1) as f64
        );
    }
    if !pb_id_hits.is_empty() {
        let mut ids: Vec<_> = pb_id_hits.into_iter().collect();
        ids.sort_by_key(|(_, c)| std::cmp::Reverse(*c));
        println!("\n--- probe barcode ids seen at offset=0, forward (top 15) ---");
        for (id, count) in ids.iter().take(15) {
            println!("  {id:<10} {count}");
        }
        println!("  distinct ids at offset=0: {}", ids.len());
    }

    println!("\n--- CRISPR anchor: hit rate by (feature id, offset in emitted R2) ---");
    let mut a_rows: Vec<_> = anchor_hits.into_iter().collect();
    a_rows.sort_by_key(|(_, c)| std::cmp::Reverse(*c));
    for ((id, offset), count) in a_rows.iter().take(15) {
        println!(
            "  {id:<10} offset={offset:>3} hits={count:>8} ({:.4}% of confident)",
            100.0 * *count as f64 / n_confident.max(1) as f64
        );
    }

    println!("\n--- guide capture-barcode mismatch histogram (at winning anchor offset) ---");
    let mut mm_rows: Vec<_> = guide_capture_mm_hist.into_iter().collect();
    mm_rows.sort_by_key(|(mm, _)| *mm);
    for (mm, count) in mm_rows {
        println!("  mm={mm:<3} count={count}");
    }

    println!("\n--- guide calls (capture-barcode mm<=2) ---");
    let mut g_rows: Vec<_> = guide_calls.into_iter().collect();
    g_rows.sort_by_key(|(_, c)| std::cmp::Reverse(*c));
    for (id, count) in g_rows {
        println!("  {id:<10} {count}");
    }

    Ok(())
}
