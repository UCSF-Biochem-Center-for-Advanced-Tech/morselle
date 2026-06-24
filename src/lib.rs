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

use anyhow::Result;
use needletail::parse_fastx_file;
use tracing::info;

use crate::config::Cli;
use crate::demux::{index_regions, DemuxVia, Demuxer};
use crate::extract::extract;
use crate::io::{discover_fastqs, PairWriter};
use crate::orient::orient;
use crate::pipeline::Orientation;
use crate::sample_sheet::SampleSheet;

/// Running counts for the pipeline summary.
struct Stats {
    reads: u64,
    oriented: u64,
    pairs_emitted: u64,
    r1_not_28: u64,
    adapter_trimmed: u64,
    r2_len_sum: u64,
    r2_len_max: u64,
    // demux
    undetermined: u64,
    via_dual: u64,
    via_i5: u64,
    via_i7: u64,
    via_conflict: u64,
    per_sample: Vec<u64>,
}

impl Stats {
    fn new(n_samples: usize) -> Self {
        Stats {
            reads: 0,
            oriented: 0,
            pairs_emitted: 0,
            r1_not_28: 0,
            adapter_trimmed: 0,
            r2_len_sum: 0,
            r2_len_max: 0,
            undetermined: 0,
            via_dual: 0,
            via_i5: 0,
            via_i7: 0,
            via_conflict: 0,
            per_sample: vec![0; n_samples],
        }
    }
}

/// Entry point invoked by the binary. Orients each read to canonical P5→P7,
/// extracts R1 (CB+UMI) and R2 (cDNA), and — when a sample sheet is given —
/// demultiplexes by i5/i7 into per-sample FASTQs (plus `Undetermined`). Without a
/// sheet it writes a single combined `sbx10x_S1_*` pair.
pub fn run(cli: &Cli) -> Result<()> {
    let sheet = match &cli.sample_sheet {
        Some(path) => {
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
            Some(sheet)
        }
        None => None,
    };

    let inputs = discover_fastqs(&cli.fastqs)?;
    info!(n_files = inputs.len(), "discovered FASTQ inputs");
    if inputs.is_empty() {
        anyhow::bail!("no FASTQ inputs found");
    }

    let orient_cfg = cli.params.orient_cfg();
    let extract_cfg = cli.params.extract_cfg();

    // Set up writers: per-sample + Undetermined when demuxing, else one combined.
    let demuxer = sheet
        .as_ref()
        .map(|s| Demuxer::new(s, cli.params.demux_cfg()));
    let (mut sample_writers, mut undetermined) = match &demuxer {
        Some(d) => {
            let mut w = Vec::with_capacity(d.sample_ids().len());
            for (i, id) in d.sample_ids().iter().enumerate() {
                w.push(PairWriter::create(&cli.out_dir, id, i + 1)?);
            }
            (w, Some(PairWriter::create(&cli.out_dir, "Undetermined", 0)?))
        }
        None => (vec![PairWriter::create(&cli.out_dir, "sbx10x", 1)?], None),
    };
    let index_len = demuxer.as_ref().map_or(10, |d| d.index_len());
    let mut stats = Stats::new(sample_writers.len());

    for input in &inputs {
        let mut reader = parse_fastx_file(input)?;
        while let Some(rec) = reader.next() {
            let rec = rec?;
            let id = rec.id().to_vec();
            let seq = rec.seq();
            let qual = rec.qual().unwrap_or(&[]);
            stats.reads += 1;

            let o = orient(&seq, qual, &orient_cfg);
            if o.orientation != Orientation::Unknown {
                stats.oriented += 1;
            }
            let e = extract(&o.seq, o.cb_umi, o.polyt_start, &extract_cfg);

            let (Some((r1s, r1e)), Some((r2s, r2e))) = (e.r1, e.r2) else {
                continue; // need both R1 and R2 to emit a pair
            };
            if r1e - r1s != 28 {
                stats.r1_not_28 += 1;
            }

            // Route: demux to a sample (or Undetermined), or the combined writer.
            let target: &mut PairWriter = if let Some(d) = &demuxer {
                let (i5, i7) = index_regions(&o.seq, r1s, e.adapter_pos, index_len);
                let outcome = d.assign(i5, i7);
                match outcome.via {
                    DemuxVia::Dual => stats.via_dual += 1,
                    DemuxVia::I5Only => stats.via_i5 += 1,
                    DemuxVia::I7Only => stats.via_i7 += 1,
                    DemuxVia::Conflict => stats.via_conflict += 1,
                    DemuxVia::None => {}
                }
                match outcome.sample {
                    Some(i) => {
                        stats.per_sample[i] += 1;
                        &mut sample_writers[i]
                    }
                    None => {
                        stats.undetermined += 1;
                        undetermined.as_mut().expect("undetermined writer")
                    }
                }
            } else {
                &mut sample_writers[0]
            };

            target.write(
                &id,
                &o.seq[r1s..r1e],
                &o.qual[r1s..r1e],
                &o.seq[r2s..r2e],
                &o.qual[r2s..r2e],
            )?;
            stats.pairs_emitted += 1;
            if e.adapter_trimmed {
                stats.adapter_trimmed += 1;
            }
            let len = (r2e - r2s) as u64;
            stats.r2_len_sum += len;
            stats.r2_len_max = stats.r2_len_max.max(len);
        }
        info!(file = %input.display(), reads = stats.reads, "processed FASTQ");
    }

    for w in sample_writers {
        w.finish()?;
    }
    if let Some(w) = undetermined {
        w.finish()?;
    }

    report(&stats, demuxer.as_ref());
    Ok(())
}

/// Log the run summary.
fn report(stats: &Stats, demuxer: Option<&Demuxer>) {
    let pct = |x: u64| 100.0 * x as f64 / stats.reads.max(1) as f64;
    info!(
        reads = stats.reads,
        oriented_pct = format!("{:.1}", pct(stats.oriented)),
        pairs = stats.pairs_emitted,
        pairs_pct = format!("{:.1}", pct(stats.pairs_emitted)),
        r1_not_28 = stats.r1_not_28,
        adapter_trimmed_pct = format!("{:.1}", pct(stats.adapter_trimmed)),
        r2_len_mean = (stats.pairs_emitted > 0).then(|| stats.r2_len_sum / stats.pairs_emitted),
        r2_len_max = stats.r2_len_max,
        "pipeline complete"
    );
    if let Some(d) = demuxer {
        let assigned = stats.pairs_emitted - stats.undetermined;
        info!(
            assigned,
            assigned_pct = format!("{:.1}", pct(assigned)),
            undetermined = stats.undetermined,
            via_dual = stats.via_dual,
            via_i5_only = stats.via_i5,
            via_i7_only = stats.via_i7,
            via_conflict = stats.via_conflict,
            "demux summary"
        );
        for (i, id) in d.sample_ids().iter().enumerate() {
            info!(sample = %id, pairs = stats.per_sample[i], "  per-sample");
        }
    }
    if stats.r1_not_28 > 0 {
        tracing::warn!(count = stats.r1_not_28, "R1 records not exactly 28 bp");
    }
}
