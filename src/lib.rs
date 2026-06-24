//! sbx10x — convert raw Roche SBX reads into per-sample 10x R1/R2 FASTQs.
//!
//! See `sbx10x_implementation_plan.md` for the milestone-phased design. M0 wires
//! up the crate, CLI, logging, and a FASTQ read->write passthrough; subsequent
//! milestones fill in orientation, demux, extraction, gating, and realignment.

pub mod accounting;
pub mod anchors;
pub mod barcode;
pub mod config;
pub mod demux;
pub mod extract;
pub mod homopolymer;
pub mod io;
pub mod orient;
pub mod pipeline;
pub mod sample_sheet;

use std::path::Path;

use anyhow::Result;
use needletail::parse_fastx_file;
use tracing::info;

use crate::config::Cli;
use crate::extract::extract;
use crate::io::{discover_fastqs, FastqGzWriter};
use crate::orient::orient;
use crate::sample_sheet::SampleSheet;

/// Running counts for the M2 pipeline summary.
struct Stats {
    reads: u64,
    oriented: u64,
    pairs_emitted: u64,
    r1_not_28: u64,
    adapter_trimmed: u64,
    r2_len_sum: u64,
    r2_len_min: u64,
    r2_len_max: u64,
}

impl Default for Stats {
    fn default() -> Self {
        Stats {
            reads: 0,
            oriented: 0,
            pairs_emitted: 0,
            r1_not_28: 0,
            adapter_trimmed: 0,
            r2_len_sum: 0,
            r2_len_min: u64::MAX,
            r2_len_max: 0,
        }
    }
}

/// Entry point invoked by the binary. M2: orient each read to canonical P5→P7,
/// extract R1 (CB+UMI) and R2 (cDNA), and emit combined `R1.fastq.gz` /
/// `R2.fastq.gz` (demux comes in M3). Reports validation stats.
pub fn run(cli: &Cli) -> Result<()> {
    if let Some(path) = &cli.sample_sheet {
        let sheet = SampleSheet::load(path)?;
        info!(
            samples = sheet.samples.len(),
            i5_revcomp = sheet.i5_reverse_complemented,
            chemistry = ?sheet.chemistry,
            "loaded sample sheet"
        );
        if sheet.source == crate::sample_sheet::SourceFormat::BclConvert {
            let flat_path = cli.out_dir.join("sample_sheet.csv");
            sheet.write_flat_csv(&flat_path)?;
            info!(out = %flat_path.display(), "wrote normalized flat sample sheet");
        }
    }

    let inputs = discover_fastqs(&cli.fastqs)?;
    info!(n_files = inputs.len(), "discovered FASTQ inputs");
    if inputs.is_empty() {
        anyhow::bail!("no FASTQ inputs found");
    }

    let orient_cfg = cli.params.orient_cfg();
    let extract_cfg = cli.params.extract_cfg();
    let mut r1w = FastqGzWriter::create(&cli.out_dir.join("R1.fastq.gz"))?;
    let mut r2w = FastqGzWriter::create(&cli.out_dir.join("R2.fastq.gz"))?;
    let mut stats = Stats::default();

    for input in &inputs {
        process_file(input, &orient_cfg, &extract_cfg, &mut r1w, &mut r2w, &mut stats)?;
    }
    r1w.finish()?;
    r2w.finish()?;

    let pct = |x: u64| 100.0 * x as f64 / stats.reads.max(1) as f64;
    info!(
        reads = stats.reads,
        oriented_pct = format!("{:.1}", pct(stats.oriented)),
        pairs = stats.pairs_emitted,
        pairs_pct = format!("{:.1}", pct(stats.pairs_emitted)),
        r1_not_28 = stats.r1_not_28,
        adapter_trimmed_pct = format!("{:.1}", pct(stats.adapter_trimmed)),
        r2_len_mean = stats.pairs_emitted.gt(&0).then(|| stats.r2_len_sum / stats.pairs_emitted),
        r2_len_min = (stats.r2_len_min != u64::MAX).then_some(stats.r2_len_min),
        r2_len_max = stats.r2_len_max,
        "M2 extract complete"
    );
    if stats.r1_not_28 > 0 {
        tracing::warn!(count = stats.r1_not_28, "R1 records not exactly 28 bp");
    }
    Ok(())
}

/// Orient + extract every read in one file, writing paired R1/R2 records.
fn process_file(
    input: &Path,
    orient_cfg: &crate::orient::OrientCfg,
    extract_cfg: &crate::extract::ExtractCfg,
    r1w: &mut FastqGzWriter,
    r2w: &mut FastqGzWriter,
    stats: &mut Stats,
) -> Result<()> {
    let mut reader = parse_fastx_file(input)?;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        let id = rec.id();
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        stats.reads += 1;

        let o = orient(&seq, qual, orient_cfg);
        if o.orientation != crate::pipeline::Orientation::Unknown {
            stats.oriented += 1;
        }
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, extract_cfg);

        // Emit only complete pairs (R1 and R2 both present) for M2.
        if let (Some((r1s, r1e)), Some((r2s, r2e))) = (e.r1, e.r2) {
            if r1e - r1s != 28 {
                stats.r1_not_28 += 1;
            }
            r1w.write_record(id, &o.seq[r1s..r1e], &o.qual[r1s..r1e])?;
            r2w.write_record(id, &o.seq[r2s..r2e], &o.qual[r2s..r2e])?;
            stats.pairs_emitted += 1;
            if e.adapter_trimmed {
                stats.adapter_trimmed += 1;
            }
            let len = (r2e - r2s) as u64;
            stats.r2_len_sum += len;
            stats.r2_len_min = stats.r2_len_min.min(len);
            stats.r2_len_max = stats.r2_len_max.max(len);
        }
    }
    info!(file = %input.display(), reads = stats.reads, "processed FASTQ");
    Ok(())
}
