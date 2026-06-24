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

use crate::accounting::{pct, DemuxMetrics, Metrics, SampleCount, WhitelistMetrics};
use crate::barcode::Whitelist;
use crate::config::Cli;
use crate::demux::{index_regions, DemuxVia, Demuxer};
use crate::extract::extract;
use crate::io::{discover_fastqs, PairWriter};
use crate::orient::orient;
use crate::pipeline::{CbStatus, Orientation, ReadClass};
use crate::sample_sheet::SampleSheet;

/// Running counts for the pipeline summary.
#[derive(Default)]
struct Stats {
    reads: u64,
    oriented: u64,
    pairs_emitted: u64,
    r1_eq_28: u64,
    r1_not_28: u64,
    adapter_trimmed: u64,
    r2_len_sum: u64,
    r2_len_max: u64,
    // read classes
    class_full: u64,
    class_partial: u64,
    class_unparseable: u64,
    // cb gate
    cb_on_whitelist: u64,
    cb_within1: u64,
    cb_fail: u64,
    gated_out: u64,
    realigned: u64,
    // demux
    undetermined: u64,
    via_dual: u64,
    via_i5: u64,
    via_i7: u64,
    via_conflict: u64,
    per_sample: Vec<u64>,
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
    let cb_len = cli.params.cb_len;
    let umi_len = cli.params.umi_len;

    // Load the cell-barcode whitelist (M4 gate), if provided.
    let whitelist = match &cli.whitelist {
        Some(path) => {
            let w = Whitelist::load(path)?;
            info!(barcodes = w.len(), path = %path.display(), "loaded CB whitelist");
            Some(w)
        }
        None => None,
    };

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
    let mut stats = Stats { per_sample: vec![0; sample_writers.len()], ..Default::default() };

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

            // Classify the cell barcode (CB = first cb_len bp of R1) if gating.
            let cb_status = match (&whitelist, e.r1) {
                (Some(w), Some((r1s, _))) => Some(w.classify(&o.seq[r1s..r1s + cb_len])),
                _ => None,
            };
            match cb_status {
                Some(CbStatus::OnWhitelist) => stats.cb_on_whitelist += 1,
                Some(CbStatus::Within1) => stats.cb_within1 += 1,
                Some(CbStatus::Fail) => stats.cb_fail += 1,
                _ => {}
            }

            // M5: when the CB fails and the measured CB+UMI span is one off 28,
            // try to rescue it via single-indel realignment (off unless --realign).
            // On success we synthesize R1 = corrected CB + UMI re-registered off the
            // dT anchor.
            let mut realigned_r1: Option<(Vec<u8>, Vec<u8>)> = None;
            if cli.params.realign && cb_status == Some(CbStatus::Fail) {
                if let (Some(w), Some((r1s, _)), Some(pts)) = (&whitelist, e.r1, o.polyt_start) {
                    let delta = (pts as i64 - r1s as i64) - (cb_len + umi_len) as i64;
                    if (delta == 1 || delta == -1) && pts >= umi_len && pts <= o.seq.len() {
                        if let Some(cb) = w.realign_cb(&o.seq[r1s..], delta as i32) {
                            let mut s = cb.to_vec();
                            s.extend_from_slice(&o.seq[pts - umi_len..pts]);
                            let mut q = o.qual[r1s..(r1s + cb_len).min(o.qual.len())].to_vec();
                            q.resize(cb_len, b'I');
                            q.extend_from_slice(&o.qual[pts - umi_len..pts]);
                            realigned_r1 = Some((s, q));
                            stats.realigned += 1;
                        }
                    }
                }
            }
            let cb_ok = realigned_r1.is_some()
                || !matches!(cb_status, Some(CbStatus::Fail));

            // Tier the read.
            let class = classify_read(o.orientation, e.r1.is_some(), e.r2.is_some(), cb_ok);
            match class {
                ReadClass::Full => stats.class_full += 1,
                ReadClass::Partial => stats.class_partial += 1,
                ReadClass::Unparseable => stats.class_unparseable += 1,
            }

            let (Some((r1s, r1e)), Some((r2s, r2e))) = (e.r1, e.r2) else {
                continue; // need both R1 and R2 to emit a pair
            };

            // Gate: drop pairs whose CB fails the whitelist (and wasn't rescued),
            // since Cell Ranger can't correct them either — unless gating is off.
            if cli.params.cb_gate && !cb_ok {
                stats.gated_out += 1;
                continue;
            }

            // R1 is either the synthesized realigned pair or the raw 28 bp slice.
            let (r1_seq, r1_qual): (&[u8], &[u8]) = match &realigned_r1 {
                Some((s, q)) => (s, q),
                None => (&o.seq[r1s..r1e], &o.qual[r1s..r1e]),
            };
            if r1_seq.len() == 28 {
                stats.r1_eq_28 += 1;
            } else {
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

            target.write(&id, r1_seq, r1_qual, &o.seq[r2s..r2e], &o.qual[r2s..r2e])?;
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
    let metrics = build_metrics(&stats, whitelist.as_ref(), demuxer.as_ref());
    let metrics_path = cli.out_dir.join("metrics_summary.json");
    metrics.write_json(&metrics_path)?;
    info!(out = %metrics_path.display(), "wrote metrics summary");
    Ok(())
}

/// Tier a read: `Unparseable` if it couldn't be oriented or no barcode located;
/// `Full` if R1+R2 present and the CB is usable (on-whitelist / within-1 /
/// realigned, or no whitelist given); `Partial` otherwise.
fn classify_read(orientation: Orientation, has_r1: bool, has_r2: bool, cb_ok: bool) -> ReadClass {
    if orientation == Orientation::Unknown || !has_r1 {
        ReadClass::Unparseable
    } else if has_r2 && cb_ok {
        ReadClass::Full
    } else {
        ReadClass::Partial
    }
}

/// Assemble the serializable metrics summary.
fn build_metrics(stats: &Stats, whitelist: Option<&Whitelist>, demuxer: Option<&Demuxer>) -> Metrics {
    let r2_len_mean = if stats.pairs_emitted > 0 {
        stats.r2_len_sum / stats.pairs_emitted
    } else {
        0
    };
    Metrics {
        total_reads: stats.reads,
        oriented: stats.oriented,
        oriented_pct: pct(stats.oriented, stats.reads),
        class_full: stats.class_full,
        class_partial: stats.class_partial,
        class_unparseable: stats.class_unparseable,
        pairs_emitted: stats.pairs_emitted,
        pairs_pct: pct(stats.pairs_emitted, stats.reads),
        r1_eq_28: stats.r1_eq_28,
        r1_not_28: stats.r1_not_28,
        adapter_trimmed: stats.adapter_trimmed,
        r2_len_mean,
        r2_len_max: stats.r2_len_max,
        whitelist: whitelist.map(|w| {
            let classified = stats.cb_on_whitelist + stats.cb_within1 + stats.cb_fail;
            WhitelistMetrics {
                size: w.len(),
                cb_on_whitelist: stats.cb_on_whitelist,
                cb_within1: stats.cb_within1,
                cb_fail: stats.cb_fail,
                cb_hit_rate: pct(stats.cb_on_whitelist + stats.cb_within1, classified),
                gated_out: stats.gated_out,
                realigned: stats.realigned,
            }
        }),
        demux: demuxer.map(|d| {
            let assigned = stats.pairs_emitted - stats.undetermined;
            DemuxMetrics {
                assigned,
                assigned_pct: pct(assigned, stats.pairs_emitted),
                undetermined: stats.undetermined,
                via_dual: stats.via_dual,
                via_i5_only: stats.via_i5,
                via_i7_only: stats.via_i7,
                via_conflict: stats.via_conflict,
                per_sample: d
                    .sample_ids()
                    .iter()
                    .enumerate()
                    .map(|(i, id)| SampleCount {
                        sample: id.clone(),
                        pairs: stats.per_sample[i],
                    })
                    .collect(),
            }
        }),
    }
}

/// Log the run summary.
fn report(stats: &Stats, demuxer: Option<&Demuxer>) {
    let p = |x: u64| 100.0 * x as f64 / stats.reads.max(1) as f64;
    info!(
        reads = stats.reads,
        oriented_pct = format!("{:.1}", p(stats.oriented)),
        full = stats.class_full,
        partial = stats.class_partial,
        unparseable = stats.class_unparseable,
        pairs = stats.pairs_emitted,
        pairs_pct = format!("{:.1}", p(stats.pairs_emitted)),
        r2_len_mean = (stats.pairs_emitted > 0).then(|| stats.r2_len_sum / stats.pairs_emitted),
        r2_len_max = stats.r2_len_max,
        "pipeline complete"
    );
    if stats.cb_on_whitelist + stats.cb_within1 + stats.cb_fail > 0 {
        let classified = stats.cb_on_whitelist + stats.cb_within1 + stats.cb_fail;
        info!(
            on_whitelist = stats.cb_on_whitelist,
            within1 = stats.cb_within1,
            fail = stats.cb_fail,
            hit_rate = format!(
                "{:.1}",
                100.0 * (stats.cb_on_whitelist + stats.cb_within1) as f64 / classified.max(1) as f64
            ),
            gated_out = stats.gated_out,
            "CB gate"
        );
    }
    if let Some(d) = demuxer {
        let assigned = stats.pairs_emitted - stats.undetermined;
        info!(
            assigned,
            assigned_pct = format!("{:.1}", p(assigned)),
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
