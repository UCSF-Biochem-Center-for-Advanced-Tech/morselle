//! sbx10x — convert raw Roche SBX reads into per-sample 10x R1/R2 FASTQs.
//!
//! See `sbx10x_implementation_plan.md` for the milestone-phased design. M0 wires
//! up the crate, CLI, logging, and a FASTQ read->write passthrough; subsequent
//! milestones fill in orientation, demux, extraction, gating, and realignment.

pub mod accounting;
pub mod anchors;
pub mod annotate_bam;
pub mod barcode;
pub mod cell_calling;
pub mod config;
pub mod cr_bam_to_mtx;
pub mod demux;
pub mod extract;
pub mod feature_capture;
pub mod gzin;
pub mod hamming_index;
pub mod homopolymer;
pub mod io;
pub mod orient;
pub mod overlap;
pub mod pipeline;
pub mod probe_barcode;
pub mod probe_demux;
pub mod probe_sample_sheet;
pub mod protospacer_calling;
pub mod sample_sheet;
pub mod yc;

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use anyhow::Result;
use crossbeam_channel::{bounded, Receiver, Sender};
use rayon::prelude::*;
use tracing::info;

use crate::accounting::{pct, DemuxMetrics, Metrics, ProbeBarcodeMetrics, SampleCount, WhitelistMetrics};
use crate::anchors::comp;
use crate::barcode::Whitelist;
use crate::config::Cli;
use crate::demux::{index_anchor_starts, DemuxVia, Demuxer};
use crate::extract::{extract, ExtractCfg};
use crate::gzin::{open_fastx, open_fastx_scan};
use crate::io::{discover_fastqs, PairShardWriter, ShardWriter};
use crate::orient::{orient, OrientCfg, OrientMode};
use crate::pipeline::{CbStatus, Orientation, ReadClass};
use crate::probe_barcode::ProbeWhitelist;
use crate::sample_sheet::SampleSheet;

/// Flush a target's batch buffer once its raw R1 bytes reach this size. Bounds
/// memory (≈ threads × targets × this) while keeping gzip members reasonably sized.
const FLUSH_BYTES: usize = 256 * 1024;

/// Minimum dT-primer run length to count as a confident 3' flank for the CB.
/// Together with a confirmed R1-tail (5' flank) this brackets the 16 bp CB exactly,
/// so a whitelist miss is a genuine base error rather than a mis-location.
const MIN_POLYT_FLANK: u16 = 10;

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
    reframed: u64,
    corrected2: u64,
    /// CBs whose location is dual-anchor confirmed (R1-tail flank + polyT flank).
    cb_located_confident: u64,
    /// CB-fails that are nonetheless confidently located (= genuine base errors).
    cb_fail_confident: u64,
    // demux
    undetermined: u64,
    /// Reads written to the `--failed-output` side channel.
    failed_emitted: u64,
    via_dual: u64,
    via_i5: u64,
    via_i7: u64,
    via_conflict: u64,
    per_sample: Vec<u64>,
    // probe barcode (FLEX_SUPPORT_PLAN.md Phase 1)
    probe_located: u64,
    probe_hit: u64,
}

impl Stats {
    /// A zeroed `Stats` sized for `n_samples` per-sample counters.
    fn new(n_samples: usize) -> Self {
        Stats {
            per_sample: vec![0; n_samples],
            ..Default::default()
        }
    }

    /// Accumulate a per-file `Stats` into this total.
    fn merge(&mut self, o: Stats) {
        self.reads += o.reads;
        self.oriented += o.oriented;
        self.pairs_emitted += o.pairs_emitted;
        self.r1_eq_28 += o.r1_eq_28;
        self.r1_not_28 += o.r1_not_28;
        self.adapter_trimmed += o.adapter_trimmed;
        self.r2_len_sum += o.r2_len_sum;
        self.r2_len_max = self.r2_len_max.max(o.r2_len_max);
        self.class_full += o.class_full;
        self.class_partial += o.class_partial;
        self.class_unparseable += o.class_unparseable;
        self.cb_on_whitelist += o.cb_on_whitelist;
        self.cb_within1 += o.cb_within1;
        self.cb_fail += o.cb_fail;
        self.gated_out += o.gated_out;
        self.realigned += o.realigned;
        self.reframed += o.reframed;
        self.corrected2 += o.corrected2;
        self.cb_located_confident += o.cb_located_confident;
        self.cb_fail_confident += o.cb_fail_confident;
        self.undetermined += o.undetermined;
        self.failed_emitted += o.failed_emitted;
        self.via_dual += o.via_dual;
        self.via_i5 += o.via_i5;
        self.via_i7 += o.via_i7;
        self.via_conflict += o.via_conflict;
        self.probe_located += o.probe_located;
        self.probe_hit += o.probe_hit;
        if self.per_sample.len() < o.per_sample.len() {
            self.per_sample.resize(o.per_sample.len(), 0);
        }
        for (a, b) in self.per_sample.iter_mut().zip(o.per_sample.iter()) {
            *a += b;
        }
    }
}

/// Wall-clock breakdown of a run, for the summary log and metrics JSON.
struct Timing {
    elapsed: f64,
    load_secs: f64,
    process_secs: f64,
    threads: usize,
}

/// Read-only context shared across the parallel per-file workers.
struct Shared<'a> {
    orient_cfg: OrientCfg,
    orient_mode: OrientMode,
    extract_cfg: ExtractCfg,
    demuxer: Option<&'a Demuxer>,
    whitelist: Option<&'a Whitelist>,
    /// 10x Flex probe-barcode whitelist (FLEX_SUPPORT_PLAN.md Phase 1). `None`
    /// unless `--probe-barcode-whitelist` was given — zero cost otherwise.
    probe_wl: Option<&'a ProbeWhitelist>,
    cb_len: usize,
    umi_len: usize,
    cb_gate: bool,
    realign: bool,
    cb_reframe: bool,
    reframe_window: usize,
    lenient: bool,
    index_len: usize,
    /// Max frame shift tried when matching i5/i7 (absorbs index-shifting indels).
    index_window: isize,
    /// Number of sample targets; the Undetermined target (when demuxing) is at
    /// index `n_samples`.
    n_samples: usize,
    /// Output chunks per logical target (`--split`). Physical writer for logical
    /// target `t`, chunk `c` is `targets[t * split + c]`.
    split: usize,
    /// Emit only R2 FASTQs with 10x-style barcode/UMI SAM tags in the
    /// FASTQ comment, for minimap2 `-y` pass-through.
    minimap2: bool,
    /// Add SBX processing provenance tags to R2 FASTQ comments.
    processing_tags: bool,
    targets: &'a [PairShardWriter],
    /// Side channel for reads that never produce a pair. See `--failed-output`.
    failed: Option<&'a ShardWriter>,
}

/// Entry point invoked by the binary. Orients each read to canonical P5→P7,
/// extracts R1 (CB+UMI) and R2 (cDNA), and — when a sample sheet is given —
/// demultiplexes by i5/i7 into per-sample FASTQs (plus `Undetermined`). Without a
/// sheet it writes a single combined `sbxpress_S1_*` pair.
pub fn run(cli: &Cli) -> Result<()> {
    if cli.params.cellranger && cli.params.minimap2 {
        anyhow::bail!("--cellranger and --minimap2 are mutually exclusive output presets");
    }
    let t_start = Instant::now();
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

    // Load the Flex probe-barcode whitelist (FLEX_SUPPORT_PLAN.md Phase 1), if
    // provided. `--minimap2` is required (clap `requires`) since PR/PB ride the
    // same FASTQ-comment -> SAM-tag mechanism as CR/CB.
    let probe_wl = match &cli.params.probe_barcode_whitelist {
        Some(path) => {
            let w = ProbeWhitelist::load(
                path,
                cli.params.probe_barcode_max_mm,
                cli.params.probe_barcode_window,
            )?;
            info!(barcodes = w.len(), path = %path.display(), "loaded probe-barcode whitelist");
            Some(w)
        }
        None => None,
    };

    // Output targets: per-sample then Undetermined (last) when demuxing, else one
    // combined target. Each file is a single thread-safe gzip member.
    let demuxer = sheet
        .as_ref()
        .map(|s| Demuxer::new(s, cli.params.demux_cfg()));
    let n_samples = demuxer.as_ref().map_or(1, |d| d.sample_ids().len());
    let glevel = cli.params.output_gzip_level();
    let split = cli.params.split.max(1);
    // Each logical target gets `split` chunk files (_001.._00N). Layout is
    // [target0 chunk0..N-1][target1 chunk0..N-1]...; physical = t*split + c.
    let mk = |prefix: &str, s_num: usize, out: &mut Vec<PairShardWriter>| -> Result<()> {
        for c in 0..split {
            out.push(PairShardWriter::create(
                &cli.out_dir,
                prefix,
                s_num,
                c + 1,
                glevel,
                !cli.params.minimap2,
            )?);
        }
        Ok(())
    };
    let targets: Vec<PairShardWriter> = match &demuxer {
        Some(d) => {
            let mut v = Vec::with_capacity((n_samples + 1) * split);
            for (i, id) in d.sample_ids().iter().enumerate() {
                mk(id, i + 1, &mut v)?;
            }
            mk("Undetermined", 0, &mut v)?;
            v
        }
        None => {
            let mut v = Vec::with_capacity(split);
            mk("sbxpress", 1, &mut v)?;
            v
        }
    };
    let index_len = demuxer.as_ref().map_or(10, |d| d.index_len());

    let failed_writer = cli
        .params
        .failed_output
        .as_ref()
        .map(|path| ShardWriter::create(path, glevel))
        .transpose()?;
    let shared = Shared {
        orient_cfg,
        orient_mode: cli.params.orientation,
        extract_cfg,
        demuxer: demuxer.as_ref(),
        whitelist: whitelist.as_ref(),
        probe_wl: probe_wl.as_ref(),
        cb_len,
        umi_len,
        cb_gate: cli.params.cb_gate,
        realign: cli.params.realign,
        cb_reframe: cli.params.cb_reframe,
        reframe_window: cli.params.reframe_window,
        lenient: cli.params.lenient,
        index_len,
        index_window: cli.params.index_window,
        n_samples,
        split,
        minimap2: cli.params.minimap2,
        processing_tags: cli.params.processing_tags,
        targets: &targets,
        failed: failed_writer.as_ref(),
    };

    // Process files in parallel — each plain-gzip file inflates serially on one
    // worker, so parallelism comes from running many SBX files at once. Sort
    // largest-first (LPT) for better load balance across the rayon pool.
    let mut inputs = inputs;
    inputs.sort_by_key(|p| std::cmp::Reverse(std::fs::metadata(p).map(|m| m.len()).unwrap_or(0)));

    // Everything up to here (sample sheet, file discovery, whitelist load, target
    // creation) is one-time setup; time it separately from per-read processing.
    let load_secs = t_start.elapsed().as_secs_f64();

    // Path selection. With many files (files ≥ workers) the one-file-per-thread
    // path is optimal: each file inflates+processes on its own core, giving up to
    // `workers` concurrent decodes. With few files (e.g. a single large one) that
    // path leaves cores idle — serial gzip inflate caps one file to one core — so
    // the pipelined path dedicates a reader per file and hands record batches to the
    // full worker pool. `--reader-threads N` forces the pipelined path with N readers.
    let n_workers = rayon::current_num_threads();
    let forced_readers = cli.params.reader_threads;
    let use_pipeline = forced_readers > 0 || inputs.len() < n_workers;
    info!(decoder = %crate::gzin::backend_description(), "gzip input decoder");

    let t_proc = Instant::now();
    let locals = if use_pipeline {
        let n_readers = if forced_readers > 0 {
            forced_readers.min(inputs.len())
        } else {
            inputs.len()
        }
        .max(1);
        info!(
            mode = "pipelined",
            readers = n_readers,
            workers = n_workers,
            n_files = inputs.len(),
            "starting parallel processing"
        );
        run_pipelined(&inputs, &shared, n_readers, n_workers)?
    } else {
        info!(
            mode = "per-file",
            rayon_threads = n_workers,
            n_files = inputs.len(),
            "starting parallel processing"
        );
        inputs
            .par_iter()
            .map(|input| process_file(input, &shared))
            .collect::<Result<Vec<_>>>()?
    };

    let mut stats = Stats::new(n_samples);
    for s in locals {
        stats.merge(s);
    }

    for t in &targets {
        t.finish()?;
    }
    if let Some(writer) = &failed_writer {
        writer.finish()?;
    }
    let process_secs = t_proc.elapsed().as_secs_f64();

    let elapsed = t_start.elapsed().as_secs_f64();
    let threads = rayon::current_num_threads();
    let timing = Timing {
        elapsed,
        load_secs,
        process_secs,
        threads,
    };
    report(&stats, demuxer.as_ref(), &timing);
    let metrics = build_metrics(
        &stats,
        whitelist.as_ref(),
        demuxer.as_ref(),
        probe_wl.as_ref(),
        &timing,
        cli.params.minimap2,
        cli.params.processing_tags,
    );
    let metrics_path = cli.out_dir.join("metrics_summary.json");
    metrics.write_json(&metrics_path)?;
    info!(out = %metrics_path.display(), "wrote metrics summary");
    Ok(())
}

/// Append one FASTQ record (id without '@') to a raw byte buffer.
fn push_fastq(buf: &mut Vec<u8>, id: &[u8], seq: &[u8], qual: &[u8]) {
    buf.push(b'@');
    buf.extend_from_slice(id);
    buf.push(b'\n');
    buf.extend_from_slice(seq);
    buf.extend_from_slice(b"\n+\n");
    buf.extend_from_slice(qual);
    buf.push(b'\n');
}

/// Append one FASTQ record whose SEQ is the reverse-complement of `seq[start..end]`
/// and whose QUAL is the reverse of `qual[start..end]`, optionally capped to
/// `max_len` bases (0 = uncapped). This is the R2 emit: the cap keeps the first
/// `max_len` bases of the reverse-complemented read (the R2 primer's first cycles),
/// so we only touch the last `max_len` bases of the cDNA — no `revcomp` allocation
/// and no work past the cap. Returns the emitted length. Identical bytes to
/// `revcomp(&seq[start..end])` reversed-qual, then truncated.
fn push_fastq_r2_rc(
    buf: &mut Vec<u8>,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    start: usize,
    end: usize,
    max_len: usize,
) -> usize {
    let full = end - start;
    let out = if max_len > 0 { full.min(max_len) } else { full };
    buf.reserve(id.len() + 2 * out + 6);
    buf.push(b'@');
    buf.extend_from_slice(id);
    buf.push(b'\n');
    for k in 0..out {
        buf.push(comp(seq[end - 1 - k]));
    }
    buf.extend_from_slice(b"\n+\n");
    for k in 0..out {
        buf.push(qual[end - 1 - k]);
    }
    buf.push(b'\n');
    out
}

/// Build a FASTQ comment that minimap2 `-y` copies verbatim to SAM/BAM as
/// auxiliary tags. UB is a pre-alignment (provisional) UMI; downstream UMI
/// correction after gene assignment may replace it.
fn minimap2_header(
    id: &[u8],
    cr: &[u8],
    cb: &[u8],
    ur: &[u8],
    ub: &[u8],
    cy: &[u8],
    uy: &[u8],
) -> Vec<u8> {
    let mut header = Vec::with_capacity(id.len() + 6 * 32 + 16);
    header.extend_from_slice(id);
    for (tag, value) in [
        (b"\tCR:Z:".as_slice(), cr),
        (b"\tCB:Z:".as_slice(), cb),
        (b"\tUR:Z:".as_slice(), ur),
        (b"\tUB:Z:".as_slice(), ub),
        (b"\tCY:Z:".as_slice(), cy),
        (b"\tUY:Z:".as_slice(), uy),
    ] {
        header.extend_from_slice(tag);
        header.extend_from_slice(value);
    }
    header
}

/// Compact, versioned SBX provenance carried from extraction into the alignment.
#[derive(Debug, Clone, Copy)]
struct ProcessingTags {
    orientation: u8,
    barcode_route: u8,
    index_route: u8,
    anchor_bits: u8,
    extracted_cdna_len: usize,
}

fn append_processing_tags(header: &mut Vec<u8>, tags: ProcessingTags) {
    header.extend_from_slice(b"\tXO:A:");
    header.push(tags.orientation);
    for (name, value) in [
        ("\tXB:i:", tags.barcode_route as usize),
        ("\tXI:i:", tags.index_route as usize),
        ("\tXA:i:", tags.anchor_bits as usize),
        ("\tXL:i:", tags.extracted_cdna_len),
    ] {
        header.extend_from_slice(name.as_bytes());
        header.extend_from_slice(value.to_string().as_bytes());
    }
}

/// Append the Flex probe-barcode tags (FLEX_SUPPORT_PLAN.md Phase 1): `PR` is
/// the raw observed 10bp probe barcode (whitelist orientation), `PB` the
/// corrected probe-barcode id (e.g. `D-A01`), empty when unmatched.
fn append_probe_tags(header: &mut Vec<u8>, raw: &[u8], id: Option<&str>) {
    header.extend_from_slice(b"\tPR:Z:");
    header.extend_from_slice(raw);
    header.extend_from_slice(b"\tPB:Z:");
    header.extend_from_slice(id.unwrap_or("").as_bytes());
}

/// Allocate the per-LOGICAL-target raw FASTQ buffers `(R1 bytes, R2 bytes, count)`.
/// Pre-sized to the flush threshold (plus slack for the record that crosses it) so
/// a hot target never reallocates mid-run.
fn make_target_bufs(n_targets: usize) -> Vec<(Vec<u8>, Vec<u8>, u64)> {
    (0..n_targets)
        .map(|_| {
            (
                Vec::with_capacity(FLUSH_BYTES + 4096),
                Vec::with_capacity(FLUSH_BYTES + 4096),
                0,
            )
        })
        .collect()
}

/// Flush this worker's buffered failed reads to the side channel.
fn flush_failed_buf(buf: &mut Vec<u8>, sh: &Shared) -> Result<()> {
    if !buf.is_empty() {
        if let Some(writer) = sh.failed {
            writer.append(buf)?;
        }
        buf.clear();
    }
    Ok(())
}

/// Flush any records still buffered per target to this worker's `chunk` shard.
fn flush_target_bufs(bufs: &[(Vec<u8>, Vec<u8>, u64)], sh: &Shared, chunk: usize) -> Result<()> {
    for (i, buf) in bufs.iter().enumerate() {
        if buf.2 > 0 {
            sh.targets[i * sh.split + chunk].append(&buf.0, &buf.1, buf.2)?;
        }
    }
    Ok(())
}

/// Process one read: orient → extract → CB classify/rescue → gate → demux → emit
/// into `bufs` (flushing the crossing batch to `chunk`), accumulating into `stats`.
/// `orient_cfg` carries the enforced per-file orientation. Shared by both the
/// one-file-per-thread path ([`process_file`]) and the pipelined path
/// ([`worker_loop`]).
#[allow(clippy::too_many_arguments)]
fn process_one(
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    orient_cfg: &OrientCfg,
    sh: &Shared,
    bufs: &mut [(Vec<u8>, Vec<u8>, u64)],
    failed_buf: &mut Vec<u8>,
    chunk: usize,
    stats: &mut Stats,
) -> Result<()> {
    stats.reads += 1;

    {
        let o = orient(seq, qual, orient_cfg);
        if o.orientation != Orientation::Unknown {
            stats.oriented += 1;
        }
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, &sh.extract_cfg);

        // Classify the cell barcode (CB = first cb_len bp of R1) if a whitelist.
        let cb_status = match (sh.whitelist, e.r1) {
            (Some(w), Some((r1s, _))) => Some(w.classify(&o.seq[r1s..r1s + sh.cb_len])),
            _ => None,
        };
        match cb_status {
            Some(CbStatus::OnWhitelist) => stats.cb_on_whitelist += 1,
            Some(CbStatus::Within1) => stats.cb_within1 += 1,
            Some(CbStatus::Fail) => stats.cb_fail += 1,
            _ => {}
        }

        // Location confidence: the 16 bp CB is definitively bracketed when both
        // flanks are present — the R1-primer tail (5') is lock-confirmed and the
        // dT-primer run (3') is long enough. Then a whitelist miss is a true base
        // error, not a mis-location.
        let located_confident =
            o.r1_tail_confirmed && e.polyt_run_len >= MIN_POLYT_FLANK && e.r1.is_some();
        if cb_status.is_some() && located_confident {
            stats.cb_located_confident += 1;
        }
        if cb_status == Some(CbStatus::Fail) && located_confident {
            stats.cb_fail_confident += 1;
        }

        // CB rescue (synthesizes a corrected R1). M5 indel realign (--realign);
        // then CB reframe (--cb-reframe); then lenient Hamming-2 (--lenient).
        // Reframe and lenient are gated on confident location.
        // A unique Hamming-1 match is a correction, not merely a classification:
        // emit the canonical whitelist sequence in R1 so downstream CB tags and
        // fixed barcode templates refer to the same cell identifier.
        let mut corrected_r1: Option<(Vec<u8>, Vec<u8>)> = match (sh.whitelist, e.r1, cb_status) {
            (Some(w), Some((r1s, r1e)), Some(CbStatus::Within1)) => {
                w.correct_within1(&o.seq[r1s..r1s + sh.cb_len]).map(|cb| {
                    let mut s = cb.to_vec();
                    s.extend_from_slice(&o.seq[r1s + sh.cb_len..r1e]);
                    (s, o.qual[r1s..r1e].to_vec())
                })
            }
            _ => None,
        };
        // XB: how the emitted barcode passed or was rescued. Values 0/1 are
        // whitelist exact/within-1; 2/3/4 are indel/reframe/Hamming-2 rescue;
        // 5 means no whitelist and 6 means an unrescued miss emitted with gating off.
        let mut barcode_route = match cb_status {
            Some(CbStatus::OnWhitelist) => 0,
            Some(CbStatus::Within1) => 1,
            Some(CbStatus::Fail) => 6,
            _ => 5,
        };
        if cb_status == Some(CbStatus::Fail) {
            if let (Some(w), Some((r1s, r1e))) = (sh.whitelist, e.r1) {
                if sh.realign {
                    if let Some(pts) = o.polyt_start {
                        let delta = (pts as i64 - r1s as i64) - (sh.cb_len + sh.umi_len) as i64;
                        if (delta == 1 || delta == -1) && pts >= sh.umi_len && pts <= o.seq.len() {
                            if let Some(cb) = w.realign_cb(&o.seq[r1s..], delta as i32) {
                                let mut s = cb.to_vec();
                                s.extend_from_slice(&o.seq[pts - sh.umi_len..pts]);
                                let mut q =
                                    o.qual[r1s..(r1s + sh.cb_len).min(o.qual.len())].to_vec();
                                q.resize(sh.cb_len, b'I');
                                q.extend_from_slice(&o.qual[pts - sh.umi_len..pts]);
                                corrected_r1 = Some((s, q));
                                barcode_route = 2;
                                stats.realigned += 1;
                            }
                        }
                    }
                }
                // CB reframe rescue: a homopolymer indel shifted the 16 bp CB grab;
                // try nearby frames (+1 first) for a unique within-1 whitelist hit,
                // re-registering the UMI off the dT anchor (so R1 stays 28 bp).
                // Gated on confident location to keep precision.
                if corrected_r1.is_none() && sh.cb_reframe && located_confident {
                    if let Some(pts) = o.polyt_start {
                        if pts >= sh.umi_len {
                            if let Some((cb, off)) = w.reframe_cb(&o.seq, r1s, sh.reframe_window) {
                                let cb_s = (r1s as i64 + off).max(0) as usize;
                                let mut s = cb.to_vec();
                                s.extend_from_slice(&o.seq[pts - sh.umi_len..pts]);
                                let mut q =
                                    o.qual[cb_s..(cb_s + sh.cb_len).min(o.qual.len())].to_vec();
                                q.resize(sh.cb_len, b'I');
                                q.extend_from_slice(&o.qual[pts - sh.umi_len..pts]);
                                corrected_r1 = Some((s, q));
                                barcode_route = 3;
                                stats.reframed += 1;
                            }
                        }
                    }
                }
                // Lenient within-2 substitution rescue: only when the CB is
                // confidently located (so the 2 mismatches are real base errors),
                // keeping the original (correctly-positioned) UMI.
                if corrected_r1.is_none() && sh.lenient && located_confident {
                    if let Some(cb2) = w.correct_within2(&o.seq[r1s..r1s + sh.cb_len]) {
                        let mut s = cb2.to_vec();
                        s.extend_from_slice(&o.seq[r1s + sh.cb_len..r1e]);
                        corrected_r1 = Some((s, o.qual[r1s..r1e].to_vec()));
                        barcode_route = 4;
                        stats.corrected2 += 1;
                    }
                }
            }
        }
        let cb_ok = corrected_r1.is_some() || !matches!(cb_status, Some(CbStatus::Fail));

        let class = classify_read(o.orientation, e.r1.is_some(), e.r2.is_some(), cb_ok);
        match class {
            ReadClass::Full => stats.class_full += 1,
            ReadClass::Partial => stats.class_partial += 1,
            ReadClass::Unparseable => stats.class_unparseable += 1,
        }

        // Emit the ORIGINAL read (pre-orientation, untrimmed) to the side
        // channel on every path that does not produce a pair, so no input read
        // is silently dropped.
        let mut push_failed = |reason: &str, stats: &mut Stats| -> Result<()> {
            if sh.failed.is_none() {
                return Ok(());
            }
            stats.failed_emitted += 1;
            failed_buf.push(b'@');
            failed_buf.extend_from_slice(id);
            failed_buf.extend_from_slice(format!("\tXF:Z:{reason}").as_bytes());
            failed_buf.push(b'\n');
            failed_buf.extend_from_slice(seq);
            failed_buf.extend_from_slice(b"\n+\n");
            failed_buf.extend_from_slice(qual);
            failed_buf.push(b'\n');
            if failed_buf.len() >= FLUSH_BYTES {
                flush_failed_buf(failed_buf, sh)?;
            }
            Ok(())
        };

        let (Some((r1s, r1e)), Some((r2s, r2e))) = (e.r1, e.r2) else {
            // need both R1 and R2 to emit a pair
            let reason = if o.orientation == Orientation::Unknown {
                "unoriented"
            } else if e.r1.is_none() {
                "no_barcode"
            } else {
                "no_cdna"
            };
            push_failed(reason, stats)?;
            return Ok(());
        };

        // Gate: drop CB-fail pairs that weren't rescued (unless gating is off).
        if sh.cb_gate && !cb_ok {
            stats.gated_out += 1;
            push_failed("cb_fail", stats)?;
            return Ok(());
        }

        let (r1_seq, r1_qual): (&[u8], &[u8]) = match &corrected_r1 {
            Some((s, q)) => (s, q),
            None => (&o.seq[r1s..r1e], &o.qual[r1s..r1e]),
        };
        if r1_seq.len() == 28 {
            stats.r1_eq_28 += 1;
        } else {
            stats.r1_not_28 += 1;
        }

        // Route to a target index: sample, Undetermined (= n_samples), or combined.
        let mut index_route = 5; // XI=5: demultiplexing was not requested.
        let target_idx = if let Some(d) = sh.demuxer {
            let (i5s, i7s) = index_anchor_starts(r1s, e.adapter_pos, sh.index_len);
            let outcome = d.assign_windowed(&o.seq, i5s, i7s, sh.index_window);
            index_route = match outcome.via {
                DemuxVia::Dual => 0,
                DemuxVia::I5Only => 1,
                DemuxVia::I7Only => 2,
                DemuxVia::Conflict => 3,
                DemuxVia::None => 4,
            };
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
                    i
                }
                None => {
                    stats.undetermined += 1;
                    sh.n_samples
                }
            }
        } else {
            0
        };

        // R1 (CB+UMI) and R2 (cDNA) sit on the SAME strand of the canonical read,
        // but Illumina paired-end mates come from OPPOSITE strands. The first-strand
        // cDNA following polyT is the reverse complement of the mRNA sense, so the
        // canonical R2 slice maps antisense. Emit revcomp(R2) with reversed quality
        // to match the Illumina convention (R1 as-is, R2 flipped) — otherwise Cell
        // Ranger reports ~50% "Reads Mapped Antisense to Gene". The `--max-r2-len`
        // cap is applied here (at write, on the reverse-complemented read): the R2
        // primer reads from the adapter-proximal end, so keeping the first N bases
        // mirrors a sequencer stopping after N cycles. `push_fastq_r2_rc` writes the
        // revcomp straight into the buffer, touching only the emitted (≤N) bases.
        let buf = &mut bufs[target_idx];
        if !sh.minimap2 {
            push_fastq(&mut buf.0, id, r1_seq, r1_qual);
        }
        let raw_cb = &o.seq[r1s..r1s + sh.cb_len];
        let raw_umi = &o.seq[r1s + sh.cb_len..r1e];
        let raw_cb_qual = &o.qual[r1s..r1s + sh.cb_len];
        let mut header;
        let processing = ProcessingTags {
            orientation: match o.orientation {
                Orientation::P5Led => b'5',
                Orientation::P7Led => b'7',
                Orientation::Unknown => b'?',
            },
            barcode_route,
            index_route,
            // XA bits: 1=R1-primer tail, 2=poly-T anchor, 4=cDNA adapter
            // trimmed, 8=file-level expected/detected orientation applied.
            anchor_bits: (o.r1_tail_confirmed as u8)
                | (((e.polyt_run_len >= MIN_POLYT_FLANK) as u8) << 1)
                | ((e.adapter_trimmed as u8) << 2)
                | ((orient_cfg.expected.is_some() as u8) << 3),
            extracted_cdna_len: r2e - r2s,
        };
        let r2_id = if sh.minimap2 {
            let mut corrected_cb = r1_seq[..sh.cb_len].to_vec();
            corrected_cb.extend_from_slice(b"-1");
            header = minimap2_header(
                id,
                raw_cb,
                &corrected_cb,
                raw_umi,
                &r1_seq[sh.cb_len..],
                raw_cb_qual,
                &r1_qual[sh.cb_len..],
            );
            if sh.processing_tags {
                append_processing_tags(&mut header, processing);
            }
            if let Some(wl) = sh.probe_wl {
                if let Some(hit) = e.polyt_end.and_then(|pts| wl.locate(&o.seq, pts)) {
                    stats.probe_located += 1;
                    if hit.id.is_some() {
                        stats.probe_hit += 1;
                    }
                    append_probe_tags(&mut header, &hit.raw, hit.id.as_deref());
                }
            }
            &header
        } else if sh.processing_tags {
            header = id.to_vec();
            append_processing_tags(&mut header, processing);
            &header
        } else {
            id
        };
        let emitted_r2 = push_fastq_r2_rc(
            &mut buf.1,
            r2_id,
            &o.seq,
            &o.qual,
            r2s,
            r2e,
            sh.extract_cfg.max_r2_len,
        );
        buf.2 += 1;
        if buf.0.len().max(buf.1.len()) >= FLUSH_BYTES {
            sh.targets[target_idx * sh.split + chunk].append(&buf.0, &buf.1, buf.2)?;
            buf.0.clear();
            buf.1.clear();
            buf.2 = 0;
        }

        stats.pairs_emitted += 1;
        if e.adapter_trimmed {
            stats.adapter_trimmed += 1;
        }
        let len = emitted_r2 as u64; // emitted (post-cap) R2 length
        stats.r2_len_sum += len;
        stats.r2_len_max = stats.r2_len_max.max(len);
    }
    Ok(())
}

/// Resolve the per-file orientation to enforce (detecting it for `Auto`).
fn file_orientation(input: &Path, sh: &Shared) -> Result<Option<Orientation>> {
    let expected = match sh.orient_mode {
        OrientMode::P5Led => Some(Orientation::P5Led),
        OrientMode::P7Led => Some(Orientation::P7Led),
        OrientMode::Auto => detect_orientation(input, &sh.orient_cfg, &sh.extract_cfg)?,
    };
    info!(file = %input.display(), orientation = ?expected, "orientation");
    Ok(expected)
}

/// One-file-per-thread path (used when files ≥ workers): orient → extract → demux →
/// gate → emit every read in one file on a single worker. Returns this file's
/// `Stats` (merged by the caller). Each file inflates serially on its own thread, so
/// parallelism comes from many files running at once.
fn process_file(input: &Path, sh: &Shared) -> Result<Stats> {
    let mut stats = Stats::new(sh.n_samples);
    let expected = file_orientation(input, sh)?;
    let orient_cfg = OrientCfg {
        expected,
        ..sh.orient_cfg
    };

    // This worker writes to one chunk per logical target (so each chunk file is
    // touched by few threads — parallel per-sample compression). The rayon worker
    // index is stable for the duration of this task.
    let chunk = rayon::current_thread_index().unwrap_or(0) % sh.split;
    let n_targets = sh.targets.len() / sh.split;
    let mut bufs = make_target_bufs(n_targets);
    let mut failed_buf = Vec::new();

    let mut reader = open_fastx(input)?;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        process_one(
            rec.id(),
            &rec.seq(),
            rec.qual().unwrap_or(&[]),
            &orient_cfg,
            sh,
            &mut bufs,
            &mut failed_buf,
            chunk,
            &mut stats,
        )?;
    }
    flush_target_bufs(&bufs, sh, chunk)?;
    flush_failed_buf(&mut failed_buf, sh)?;
    Ok(stats)
}

/// Records read from one file, flat-packed for cheap cross-thread transfer: all ids
/// / seqs / quals concatenated, with per-record lengths. One `Batch` carries the
/// enforced orientation of its source file. Avoids per-record heap allocations.
struct Batch {
    orientation: Option<Orientation>,
    ids: Vec<u8>,
    seqs: Vec<u8>,
    quals: Vec<u8>,
    /// `(id_len, seq_len)` per record; `qual_len == seq_len` for FASTQ.
    lens: Vec<(u32, u32)>,
}

/// Records per batch — big enough to amortize channel/sync overhead, small enough
/// to keep memory bounded and workers fed promptly.
const BATCH_READS: usize = 4096;

impl Batch {
    fn with_orientation(orientation: Option<Orientation>) -> Self {
        Batch {
            orientation,
            ids: Vec::new(),
            seqs: Vec::new(),
            quals: Vec::new(),
            lens: Vec::with_capacity(BATCH_READS),
        }
    }

    fn push(&mut self, id: &[u8], seq: &[u8], qual: &[u8]) {
        self.ids.extend_from_slice(id);
        self.seqs.extend_from_slice(seq);
        self.quals.extend_from_slice(qual);
        self.lens.push((id.len() as u32, seq.len() as u32));
    }

    fn len(&self) -> usize {
        self.lens.len()
    }
}

/// Reader thread body: pull files from the shared work queue, detect+enforce each
/// file's orientation, and stream record batches into `tx`. One reader per file
/// (serial gzip inflate can't be parallelized within a file), so many readers only
/// help across files; a single large file uses one reader feeding the whole pool.
fn reader_loop(
    files: &[PathBuf],
    next: &AtomicUsize,
    sh: &Shared,
    tx: Sender<Batch>,
) -> Result<()> {
    loop {
        let idx = next.fetch_add(1, Ordering::Relaxed);
        let Some(input) = files.get(idx) else { break };
        let expected = file_orientation(input, sh)?;
        let mut reader = open_fastx(input)?;
        let mut batch = Batch::with_orientation(expected);
        while let Some(rec) = reader.next() {
            let rec = rec?;
            batch.push(rec.id(), &rec.seq(), rec.qual().unwrap_or(&[]));
            if batch.len() >= BATCH_READS {
                let full = std::mem::replace(&mut batch, Batch::with_orientation(expected));
                if tx.send(full).is_err() {
                    return Ok(()); // all workers gone
                }
            }
        }
        if batch.len() > 0 && tx.send(batch).is_err() {
            return Ok(());
        }
    }
    Ok(())
}

/// Worker thread body: consume batches, run [`process_one`] per record, and flush
/// output to this worker's `chunk` shard. Returns the worker's accumulated `Stats`.
fn worker_loop(worker_id: usize, sh: &Shared, rx: Receiver<Batch>) -> Result<Stats> {
    let chunk = worker_id % sh.split;
    let n_targets = sh.targets.len() / sh.split;
    let mut bufs = make_target_bufs(n_targets);
    let mut failed_buf = Vec::new();
    let mut stats = Stats::new(sh.n_samples);

    while let Ok(batch) = rx.recv() {
        let orient_cfg = OrientCfg {
            expected: batch.orientation,
            ..sh.orient_cfg
        };
        let (mut io, mut so) = (0usize, 0usize);
        for &(id_len, seq_len) in &batch.lens {
            let (id_len, seq_len) = (id_len as usize, seq_len as usize);
            let id = &batch.ids[io..io + id_len];
            let seq = &batch.seqs[so..so + seq_len];
            let qual = &batch.quals[so..so + seq_len];
            process_one(
                id,
                seq,
                qual,
                &orient_cfg,
                sh,
                &mut bufs,
                &mut failed_buf,
                chunk,
                &mut stats,
            )?;
            io += id_len;
            so += seq_len;
        }
    }
    flush_target_bufs(&bufs, sh, chunk)?;
    flush_failed_buf(&mut failed_buf, sh)?;
    Ok(stats)
}

/// Pipelined path (used when files < workers): `n_readers` reader threads decouple
/// serial gzip decode from processing, streaming batches over a bounded channel to
/// `n_workers` workers — so a single large file still saturates the machine. Returns
/// one `Stats` per worker (merged by the caller).
fn run_pipelined(
    inputs: &[PathBuf],
    sh: &Shared,
    n_readers: usize,
    n_workers: usize,
) -> Result<Vec<Stats>> {
    // Bounded so readers apply backpressure instead of racing ahead into memory.
    let (tx, rx) = bounded::<Batch>(n_workers * 2);
    let next = AtomicUsize::new(0);

    std::thread::scope(|scope| -> Result<Vec<Stats>> {
        let reader_handles: Vec<_> = (0..n_readers)
            .map(|_| {
                let tx = tx.clone();
                scope.spawn(|| reader_loop(inputs, &next, sh, tx))
            })
            .collect();
        drop(tx); // workers see the channel close once every reader clone is dropped

        let worker_handles: Vec<_> = (0..n_workers)
            .map(|w| {
                let rx = rx.clone();
                scope.spawn(move || worker_loop(w, sh, rx))
            })
            .collect();
        drop(rx);

        // Surface any reader error, then collect the workers' stats.
        for h in reader_handles {
            h.join().expect("reader thread panicked")?;
        }
        let mut locals = Vec::with_capacity(n_workers);
        for h in worker_handles {
            locals.push(h.join().expect("worker thread panicked")?);
        }
        Ok(locals)
    })
}

/// Detect a file's dominant orientation from a pre-scan of confident reads.
/// Returns `Some` only when ≥80% of confident reads agree (else `None` → per-read
/// auto). Each SBX run is orientation-pure, so this is normally unambiguous.
fn detect_orientation(
    input: &Path,
    base_cfg: &OrientCfg,
    extract_cfg: &ExtractCfg,
) -> Result<Option<Orientation>> {
    const SCAN: usize = 2000;
    let cfg = OrientCfg {
        expected: None,
        ..*base_cfg
    };
    let mut reader = open_fastx_scan(input)?;
    let (mut p5, mut p7) = (0u64, 0u64);
    let mut seen = 0;
    while let Some(rec) = reader.next() {
        if seen >= SCAN {
            break;
        }
        seen += 1;
        let rec = rec?;
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        let o = orient(&seq, qual, &cfg);
        let e = extract(&o.seq, o.cb_umi, o.polyt_start, extract_cfg);
        if o.r1_tail_confirmed && e.polyt_run_len >= 10 {
            match o.orientation {
                Orientation::P5Led => p5 += 1,
                Orientation::P7Led => p7 += 1,
                Orientation::Unknown => {}
            }
        }
    }
    let total = p5 + p7;
    if total < 50 {
        return Ok(None); // too few confident reads to decide
    }
    let frac5 = p5 as f64 / total as f64;
    Ok(if frac5 >= 0.8 {
        Some(Orientation::P5Led)
    } else if frac5 <= 0.2 {
        Some(Orientation::P7Led)
    } else {
        None
    })
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
#[allow(clippy::too_many_arguments)]
fn build_metrics(
    stats: &Stats,
    whitelist: Option<&Whitelist>,
    demuxer: Option<&Demuxer>,
    probe_wl: Option<&ProbeWhitelist>,
    t: &Timing,
    minimap2: bool,
    processing_tags: bool,
) -> Metrics {
    let r2_len_mean = if stats.pairs_emitted > 0 {
        stats.r2_len_sum / stats.pairs_emitted
    } else {
        0
    };
    let reads_per_sec = if t.elapsed > 0.0 {
        (stats.reads as f64 / t.elapsed) as u64
    } else {
        0
    };
    Metrics {
        output_mode: if minimap2 {
            "minimap2_r2_only"
        } else {
            "paired_fastq"
        }
        .into(),
        processing_tags,
        total_reads: stats.reads,
        reads_failed: stats.failed_emitted,
        elapsed_secs: t.elapsed,
        load_secs: t.load_secs,
        process_secs: t.process_secs,
        reads_per_sec,
        threads: t.threads,
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
                cb_located_confident: stats.cb_located_confident,
                cb_located_confident_pct: pct(stats.cb_located_confident, classified),
                cb_fail_confident: stats.cb_fail_confident,
                cb_fail_confident_pct: pct(stats.cb_fail_confident, stats.cb_fail),
                gated_out: stats.gated_out,
                realigned: stats.realigned,
                reframed: stats.reframed,
                corrected_hamming2: stats.corrected2,
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
        probe_barcode: probe_wl.map(|w| ProbeBarcodeMetrics {
            size: w.len(),
            located: stats.probe_located,
            hit: stats.probe_hit,
            hit_rate: pct(stats.probe_hit, stats.probe_located),
        }),
    }
}

/// Log the run summary.
fn report(stats: &Stats, demuxer: Option<&Demuxer>, t: &Timing) {
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
                100.0 * (stats.cb_on_whitelist + stats.cb_within1) as f64
                    / classified.max(1) as f64
            ),
            located_confident_pct = format!(
                "{:.1}",
                100.0 * stats.cb_located_confident as f64 / classified.max(1) as f64
            ),
            fail_confident_pct = format!(
                "{:.1}",
                100.0 * stats.cb_fail_confident as f64 / stats.cb_fail.max(1) as f64
            ),
            realigned = stats.realigned,
            reframed = stats.reframed,
            corrected2 = stats.corrected2,
            gated_out = stats.gated_out,
            "CB gate"
        );
    }
    if stats.probe_located > 0 {
        info!(
            located = stats.probe_located,
            hit = stats.probe_hit,
            hit_rate = format!(
                "{:.1}",
                100.0 * stats.probe_hit as f64 / stats.probe_located.max(1) as f64
            ),
            "probe barcode"
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

    // Headline run summary: time, throughput, and where reads ended up.
    let has_demux = demuxer.is_some();
    let assigned = if has_demux {
        stats.pairs_emitted.saturating_sub(stats.undetermined)
    } else {
        stats.pairs_emitted
    };
    let final_ok = if has_demux {
        assigned
    } else {
        stats.pairs_emitted
    };
    let lost = stats.reads.saturating_sub(final_ok);
    let no_cdna = stats.class_partial.saturating_sub(stats.gated_out);
    let rps = if t.elapsed > 0.0 {
        (stats.reads as f64 / t.elapsed) as u64
    } else {
        0
    };
    info!(
        elapsed = format!("{:.1}s ({:.1} min)", t.elapsed, t.elapsed / 60.0),
        load_secs = format!("{:.1}", t.load_secs),
        process_secs = format!("{:.1}", t.process_secs),
        threads = t.threads,
        reads_per_sec = rps,
        processed = stats.reads,
        emitted = stats.pairs_emitted,
        assigned = if has_demux { Some(assigned) } else { None },
        lost,
        lost_pct = format!("{:.1}", p(lost)),
        "RUN SUMMARY"
    );
    info!(
        unparseable = stats.class_unparseable,
        no_cdna,
        cb_fail_gated = stats.gated_out,
        index_undetermined = has_demux.then_some(stats.undetermined),
        "  lost breakdown"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::anchors::revcomp;

    /// The optimized R2 emitter must produce byte-identical output to the original
    /// `revcomp(seq[start..end])` + reversed-qual + truncate approach it replaced.
    fn reference_r2(
        id: &[u8],
        seq: &[u8],
        qual: &[u8],
        start: usize,
        end: usize,
        cap: usize,
    ) -> Vec<u8> {
        let mut r2_seq = revcomp(&seq[start..end]);
        let mut r2_qual = qual[start..end].to_vec();
        r2_qual.reverse();
        if cap > 0 && r2_seq.len() > cap {
            r2_seq.truncate(cap);
            r2_qual.truncate(cap);
        }
        let mut buf = Vec::new();
        push_fastq(&mut buf, id, &r2_seq, &r2_qual);
        buf
    }

    #[test]
    fn r2_rc_matches_reference() {
        let seq = b"ACGTACGTNNGGCCTTAAACGTACGTTTGGCCAA";
        let qual: Vec<u8> = (0..seq.len() as u8).map(|i| b'!' + (i % 40)).collect();
        let id = b"read/1";
        // Sweep sub-ranges and caps (including 0 = uncapped and caps > length).
        for start in 0..6 {
            for end in (start + 1)..=seq.len() {
                for cap in [0usize, 1, 5, 10, end - start, end - start + 5] {
                    let mut got = Vec::new();
                    let out = push_fastq_r2_rc(&mut got, id, seq, &qual, start, end, cap);
                    let want = reference_r2(id, seq, &qual, start, end, cap);
                    assert_eq!(got, want, "start={start} end={end} cap={cap}");
                    let expected_len = if cap > 0 {
                        (end - start).min(cap)
                    } else {
                        end - start
                    };
                    assert_eq!(out, expected_len);
                }
            }
        }
    }

    #[test]
    fn minimap2_header_uses_10x_style_tags() {
        let got = minimap2_header(
            b"read_id_001",
            b"AACCGGTT",
            b"AACCGGTT-1",
            b"ACTG",
            b"ACTG",
            b"IIIIIIII",
            b"JJJJ",
        );
        assert_eq!(
            got,
            b"read_id_001\tCR:Z:AACCGGTT\tCB:Z:AACCGGTT-1\tUR:Z:ACTG\tUB:Z:ACTG\tCY:Z:IIIIIIII\tUY:Z:JJJJ"
        );
    }

    #[test]
    fn processing_tags_are_valid_sam_style_fastq_comments() {
        let mut got = b"read_id_001".to_vec();
        append_processing_tags(
            &mut got,
            ProcessingTags {
                orientation: b'5',
                barcode_route: 3,
                index_route: 1,
                anchor_bits: 15,
                extracted_cdna_len: 427,
            },
        );
        assert_eq!(
            got,
            b"read_id_001\tXO:A:5\tXB:i:3\tXI:i:1\tXA:i:15\tXL:i:427"
        );
    }
}
