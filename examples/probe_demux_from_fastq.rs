//! Real-data exercise of the actual `probe_demux::assign`/`gem_overlap` logic
//! (the Flex design notes Phase 2), without needing a full alignment pipeline.
//!
//! `probe_demux::read_bam_tally` (already unit-tested against a synthetic BAM)
//! just turns a BAM into a `(CB, PB) -> molecule count` tally; this builds
//! that same tally directly from a `--minimap2`-preset R2 FASTQ's `CB`/`UB`/
//! `PB` comment tags, deduping by `(CB, UB, PB)` as a proxy for the real
//! `xf`-gated per-molecule count (no gene gate, since this predates
//! alignment/annotation — a slight overcount vs. real per-gene molecules, but
//! close enough to exercise the real assignment/overlap code against real
//! Flex data).
//!
//! Usage: cargo run --release --example probe_demux_from_fastq -- \
//!   <R2.fastq.gz> --probe-sample-sheet <config.csv> [--min-molecules N]

use std::collections::{BTreeMap, HashSet};
use std::env;
use std::path::PathBuf;

use morselle::probe_demux::{assign, gem_overlap, DemuxCfg};
use morselle::probe_sample_sheet::ProbeSampleSheet;
use needletail::parse_fastx_file;

fn main() -> anyhow::Result<()> {
    let mut fastq = None;
    let mut sheet_path = None;
    let mut min_molecules = 10u64;
    let mut args = env::args().skip(1);
    while let Some(a) = args.next() {
        match a.as_str() {
            "--probe-sample-sheet" => sheet_path = args.next().map(PathBuf::from),
            "--min-molecules" => min_molecules = args.next().unwrap().parse().unwrap(),
            other => fastq = Some(PathBuf::from(other)),
        }
    }
    let fastq =
        fastq.expect("usage: probe_demux_from_fastq <R2.fastq.gz> --probe-sample-sheet <CSV>");
    let sheet_path = sheet_path.expect("--probe-sample-sheet <CSV> required");

    let sheet = ProbeSampleSheet::load(&sheet_path)?;
    println!(
        "loaded {} samples from {}",
        sheet.samples.len(),
        sheet_path.display()
    );

    let mut seen: HashSet<(String, String, String)> = HashSet::new(); // (CB, UB, PB) dedup
    let mut tally: BTreeMap<(String, String), u64> = BTreeMap::new();
    let mut n = 0u64;
    let mut with_pb = 0u64;

    let mut reader = parse_fastx_file(&fastq)?;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        n += 1;
        let id = String::from_utf8_lossy(rec.id());
        let cb = extract_tag(&id, "CB:Z:");
        let ub = extract_tag(&id, "UB:Z:");
        let pb = extract_tag(&id, "PB:Z:");
        let (Some(cb), Some(ub), Some(pb)) = (cb, ub, pb) else {
            continue;
        };
        if pb.is_empty() {
            continue;
        }
        with_pb += 1;
        if seen.insert((cb.clone(), ub, pb.clone())) {
            *tally.entry((cb, pb)).or_insert(0) += 1;
        }
        if n % 20_000_000 == 0 {
            eprintln!("...{n} reads scanned");
        }
    }
    println!(
        "reads: {n}, with confident PB: {with_pb}, distinct (CB,UMI,PB) molecules: {}",
        seen.len()
    );

    let cfg = DemuxCfg { min_molecules };
    let assignments = assign(&sheet, &tally, &cfg);
    let declared: std::collections::HashSet<&str> = sheet
        .samples
        .iter()
        .flat_map(|s| s.probe_barcode_ids.iter().map(String::as_str))
        .collect();
    let overlap = gem_overlap(&tally, &declared);

    let mut assigned = 0u64;
    let mut unassigned = 0u64;
    let mut high_occ = 0u64;
    let mut per_sample: std::collections::HashMap<&str, u64> = std::collections::HashMap::new();
    for a in &assignments {
        use morselle::probe_demux::Status::*;
        match a.status {
            Assigned => {
                assigned += 1;
                *per_sample
                    .entry(a.sample_id.as_deref().unwrap())
                    .or_insert(0) += 1;
            }
            Unassigned => unassigned += 1,
            HighOccupancyGem => high_occ += 1,
        }
    }
    println!("\ncell barcodes total: {}", assignments.len());
    println!("Assigned: {assigned}  Unassigned: {unassigned}  HighOccupancyGem: {high_occ}");
    println!(
        "\nper-sample raw Assigned CB counts (pre-cell-calling, --min-molecules={min_molecules}):"
    );
    let mut rows: Vec<_> = sheet
        .samples
        .iter()
        .map(|s| {
            (
                s.id.as_str(),
                per_sample.get(s.id.as_str()).copied().unwrap_or(0),
            )
        })
        .collect();
    for (id, count) in &rows {
        println!("  {id:<10} {count}");
    }
    rows.sort_by_key(|&(_, c)| std::cmp::Reverse(c));

    println!("\ntop 10 gem_barcode_overlap.csv rows by overlap:");
    for r in overlap.iter().take(10) {
        println!(
            "  {} x {}  gems=({},{})  common={}  overlap={:.4}",
            r.barcode1, r.barcode2, r.barcode1_gems, r.barcode2_gems, r.common_gems, r.overlap
        );
    }

    Ok(())
}

fn extract_tag(header: &str, prefix: &str) -> Option<String> {
    let start = header.find(prefix)? + prefix.len();
    let rest = &header[start..];
    let end = rest.find('\t').unwrap_or(rest.len());
    Some(rest[..end].to_string())
}
