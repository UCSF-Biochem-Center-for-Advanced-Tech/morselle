//! Phase 6 prerequisite: where does the 10x Flex probe-set's 50bp
//! `probe_seq` actually appear in a real emitted R2 read for this SBX-adapted
//! chemistry? Not assumed to match either the probe-barcode's own canonical
//! offset or the CRISPR anchor's emitted-R2 offset -- confirmed empirically,
//! same as every other coordinate in this plan.
//!
//! Usage:
//!   cargo run --release --example probe_seq_locate -- \
//!     --probe-set <probe_set.csv> [--max-per-file N] <extracted_R2.fastq[.gz]> [...]
//!
//! Scans each read's full sequence for an *exact* 50bp substring match
//! against any real probe_seq (forward and revcomp), recording (offset,
//! orientation) to find the dominant position -- exact search is enough for
//! a first empirical pass given real per-base accuracy.

use std::collections::{HashMap, HashSet};
use std::env;
use std::fs;
use std::path::PathBuf;

use morselle::anchors::revcomp;
use needletail::parse_fastx_file;

fn main() -> anyhow::Result<()> {
    let mut probe_set_path = None;
    let mut max_per_file: usize = 100_000;
    let mut fastqs: Vec<PathBuf> = Vec::new();
    let mut args = env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--probe-set" => probe_set_path = args.next(),
            "--max-per-file" => max_per_file = args.next().unwrap().parse().unwrap(),
            other => fastqs.push(PathBuf::from(other)),
        }
    }
    let probe_set_path = probe_set_path.expect("--probe-set <CSV> required");
    anyhow::ensure!(
        !fastqs.is_empty(),
        "usage: probe_seq_locate --probe-set <CSV> [--max-per-file N] <fastq> [...]"
    );

    let text = fs::read_to_string(&probe_set_path)?;
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let header = lines.next().unwrap();
    let cols: Vec<&str> = header.split(',').collect();
    let idx = |name: &str| cols.iter().position(|c| *c == name).unwrap();
    let (seq_i, gene_i, included_i) = (idx("probe_seq"), idx("gene_name"), idx("included"));
    let mut probes: HashMap<Vec<u8>, String> = HashMap::new();
    let mut n_included = 0;
    for line in lines {
        let f: Vec<&str> = line.split(',').collect();
        if f[included_i] != "TRUE" {
            continue;
        }
        n_included += 1;
        probes.insert(f[seq_i].as_bytes().to_vec(), f[gene_i].to_string());
    }
    println!("loaded {n_included} included probes from {probe_set_path}");
    let probe_set: HashSet<&Vec<u8>> = probes.keys().collect();
    let probe_len = probes.keys().next().unwrap().len();

    let mut n_reads = 0u64;
    // (offset, revcomp?) -> hit count
    let mut hits: HashMap<(usize, bool), u64> = HashMap::new();
    let mut gene_hits: HashMap<String, u64> = HashMap::new();

    for fq in &fastqs {
        let mut reader = parse_fastx_file(fq)?;
        let mut n_file = 0usize;
        while let Some(rec) = reader.next() {
            if n_file >= max_per_file {
                break;
            }
            let rec = rec?;
            let seq = rec.seq();
            n_file += 1;
            n_reads += 1;
            if seq.len() < probe_len {
                continue;
            }
            for offset in 0..=(seq.len() - probe_len) {
                let window = &seq[offset..offset + probe_len];
                if probe_set.contains(&window.to_vec()) {
                    *hits.entry((offset, false)).or_default() += 1;
                    *gene_hits.entry(probes[window].clone()).or_default() += 1;
                    break; // leftmost hit for this read
                }
                let window_rc = revcomp(window);
                if probe_set.contains(&window_rc) {
                    *hits.entry((offset, true)).or_default() += 1;
                    *gene_hits.entry(probes[&window_rc].clone()).or_default() += 1;
                    break;
                }
            }
        }
        println!("  {}: {n_file} reads scanned", fq.display());
    }

    println!("\n=== summary ===\nreads: {n_reads}");
    let mut rows: Vec<_> = hits.into_iter().collect();
    rows.sort_by_key(|(_, c)| std::cmp::Reverse(*c));
    println!("\n--- exact probe_seq hits by (offset, revcomp) ---");
    for ((offset, rc), count) in rows.iter().take(15) {
        println!(
            "  offset={offset:>3} revcomp={rc:<5} hits={count:>8} ({:.3}% of reads)",
            100.0 * *count as f64 / n_reads.max(1) as f64
        );
    }
    println!("\n--- distinct genes hit: {} ---", gene_hits.len());
    let mut gene_rows: Vec<_> = gene_hits.into_iter().collect();
    gene_rows.sort_by_key(|(_, c)| std::cmp::Reverse(*c));
    for (gene, count) in gene_rows.iter().take(15) {
        println!("  {gene:<15} {count}");
    }
    Ok(())
}
