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
use std::time::Instant;

use anyhow::Result;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use tracing::info;

use crate::accounting::{pct, DemuxMetrics, Metrics, SampleCount, WhitelistMetrics};
use crate::anchors::revcomp;
use crate::barcode::Whitelist;
use crate::config::Cli;
use crate::demux::{index_anchor_starts, DemuxVia, Demuxer};
use crate::extract::{extract, ExtractCfg};
use crate::io::{discover_fastqs, PairShardWriter};
use crate::orient::{orient, OrientCfg, OrientMode};
use crate::pipeline::{CbStatus, Orientation, ReadClass};
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
    via_dual: u64,
    via_i5: u64,
    via_i7: u64,
    via_conflict: u64,
    per_sample: Vec<u64>,
}

impl Stats {
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
        self.via_dual += o.via_dual;
        self.via_i5 += o.via_i5;
        self.via_i7 += o.via_i7;
        self.via_conflict += o.via_conflict;
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
    targets: &'a [PairShardWriter],
}

/// Entry point invoked by the binary. Orients each read to canonical P5→P7,
/// extracts R1 (CB+UMI) and R2 (cDNA), and — when a sample sheet is given —
/// demultiplexes by i5/i7 into per-sample FASTQs (plus `Undetermined`). Without a
/// sheet it writes a single combined `sbx10x_S1_*` pair.
pub fn run(cli: &Cli) -> Result<()> {
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

    // Output targets: per-sample then Undetermined (last) when demuxing, else one
    // combined target. Each file is a single thread-safe gzip member.
    let demuxer = sheet
        .as_ref()
        .map(|s| Demuxer::new(s, cli.params.demux_cfg()));
    let n_samples = demuxer.as_ref().map_or(1, |d| d.sample_ids().len());
    let glevel = cli.params.gzip_level;
    let split = cli.params.split.max(1);
    // Each logical target gets `split` chunk files (_001.._00N). Layout is
    // [target0 chunk0..N-1][target1 chunk0..N-1]...; physical = t*split + c.
    let mk = |prefix: &str, s_num: usize, out: &mut Vec<PairShardWriter>| -> Result<()> {
        for c in 0..split {
            out.push(PairShardWriter::create(&cli.out_dir, prefix, s_num, c + 1, glevel)?);
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
            mk("sbx10x", 1, &mut v)?;
            v
        }
    };
    let index_len = demuxer.as_ref().map_or(10, |d| d.index_len());

    let shared = Shared {
        orient_cfg,
        orient_mode: cli.params.orientation,
        extract_cfg,
        demuxer: demuxer.as_ref(),
        whitelist: whitelist.as_ref(),
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
        targets: &targets,
    };

    // Process files in parallel — each plain-gzip file inflates serially on one
    // worker, so parallelism comes from running many SBX files at once. Sort
    // largest-first (LPT) for better load balance across the rayon pool.
    let mut inputs = inputs;
    inputs.sort_by_key(|p| std::cmp::Reverse(std::fs::metadata(p).map(|m| m.len()).unwrap_or(0)));

    // Everything up to here (sample sheet, file discovery, whitelist load, target
    // creation) is one-time setup; time it separately from per-read processing.
    let load_secs = t_start.elapsed().as_secs_f64();

    info!(rayon_threads = rayon::current_num_threads(), n_files = inputs.len(), "starting parallel processing");
    let t_proc = Instant::now();
    let locals = inputs
        .par_iter()
        .map(|input| process_file(input, &shared))
        .collect::<Result<Vec<_>>>()?;

    let mut stats = Stats { per_sample: vec![0; n_samples], ..Default::default() };
    for s in locals {
        stats.merge(s);
    }

    for t in &targets {
        t.finish()?;
    }
    let process_secs = t_proc.elapsed().as_secs_f64();

    let elapsed = t_start.elapsed().as_secs_f64();
    let threads = rayon::current_num_threads();
    let timing = Timing { elapsed, load_secs, process_secs, threads };
    report(&stats, demuxer.as_ref(), &timing);
    let metrics = build_metrics(&stats, whitelist.as_ref(), demuxer.as_ref(), &timing);
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

/// Orient → extract → demux → gate → realign every read in one file, buffering
/// per-target output and flushing compressed batches to the shared writers.
/// Returns this file's `Stats` (merged by the caller).
fn process_file(input: &Path, sh: &Shared) -> Result<Stats> {
    let mut stats = Stats { per_sample: vec![0; sh.n_samples], ..Default::default() };

    // Determine the run orientation for this file and enforce it (each run is
    // orientation-pure). Auto = detect from a pre-scan of confident reads.
    let expected = match sh.orient_mode {
        OrientMode::P5Led => Some(Orientation::P5Led),
        OrientMode::P7Led => Some(Orientation::P7Led),
        OrientMode::Auto => detect_orientation(input, &sh.orient_cfg, &sh.extract_cfg)?,
    };
    info!(file = %input.display(), orientation = ?expected, "orientation");
    let orient_cfg = OrientCfg { expected, ..sh.orient_cfg };

    // This worker writes to one chunk per logical target (so each chunk file is
    // touched by few threads — parallel per-sample compression). The rayon worker
    // index is stable for the duration of this task.
    let chunk = rayon::current_thread_index().unwrap_or(0) % sh.split;

    // Per-LOGICAL-target raw FASTQ buffers: (R1 bytes, R2 bytes, record count).
    let n_targets = sh.targets.len() / sh.split;
    let mut bufs: Vec<(Vec<u8>, Vec<u8>, u64)> = (0..n_targets)
        .map(|_| (Vec::new(), Vec::new(), 0))
        .collect();

    let mut reader = parse_fastx_file(input)?;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        let id = rec.id();
        let seq = rec.seq();
        let qual = rec.qual().unwrap_or(&[]);
        stats.reads += 1;

        let o = orient(&seq, qual, &orient_cfg);
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
        let mut corrected_r1: Option<(Vec<u8>, Vec<u8>)> = None;
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
                            if let Some((cb, off)) =
                                w.reframe_cb(&o.seq, r1s, sh.reframe_window)
                            {
                                let cb_s = (r1s as i64 + off).max(0) as usize;
                                let mut s = cb.to_vec();
                                s.extend_from_slice(&o.seq[pts - sh.umi_len..pts]);
                                let mut q = o.qual
                                    [cb_s..(cb_s + sh.cb_len).min(o.qual.len())]
                                    .to_vec();
                                q.resize(sh.cb_len, b'I');
                                q.extend_from_slice(&o.qual[pts - sh.umi_len..pts]);
                                corrected_r1 = Some((s, q));
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

        let (Some((r1s, r1e)), Some((r2s, r2e))) = (e.r1, e.r2) else {
            continue; // need both R1 and R2 to emit a pair
        };

        // Gate: drop CB-fail pairs that weren't rescued (unless gating is off).
        if sh.cb_gate && !cb_ok {
            stats.gated_out += 1;
            continue;
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
        let target_idx = if let Some(d) = sh.demuxer {
            let (i5s, i7s) = index_anchor_starts(r1s, e.adapter_pos, sh.index_len);
            let outcome = d.assign_windowed(&o.seq, i5s, i7s, sh.index_window);
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
        // Ranger reports ~50% "Reads Mapped Antisense to Gene".
        let mut r2_seq = revcomp(&o.seq[r2s..r2e]);
        let mut r2_qual = o.qual[r2s..r2e].to_vec();
        r2_qual.reverse();
        // Apply the R2 length cap HERE, at write, on the reverse-complemented read:
        // the R2 primer reads from the (now leading) adapter-proximal end, so keeping
        // the first `max_r2_len` bases mirrors a real sequencer stopping after that
        // many cycles. The whole fragment was retained through extraction/trimming.
        if sh.extract_cfg.max_r2_len > 0 && r2_seq.len() > sh.extract_cfg.max_r2_len {
            r2_seq.truncate(sh.extract_cfg.max_r2_len);
            r2_qual.truncate(sh.extract_cfg.max_r2_len);
        }

        let buf = &mut bufs[target_idx];
        push_fastq(&mut buf.0, id, r1_seq, r1_qual);
        push_fastq(&mut buf.1, id, &r2_seq, &r2_qual);
        buf.2 += 1;
        if buf.0.len() >= FLUSH_BYTES {
            sh.targets[target_idx * sh.split + chunk].append(&buf.0, &buf.1, buf.2)?;
            buf.0.clear();
            buf.1.clear();
            buf.2 = 0;
        }

        stats.pairs_emitted += 1;
        if e.adapter_trimmed {
            stats.adapter_trimmed += 1;
        }
        let len = r2_seq.len() as u64; // emitted (post-cap) R2 length
        stats.r2_len_sum += len;
        stats.r2_len_max = stats.r2_len_max.max(len);
    }

    // Flush any remaining buffered records (to this worker's chunk).
    for (i, buf) in bufs.iter().enumerate() {
        if buf.2 > 0 {
            sh.targets[i * sh.split + chunk].append(&buf.0, &buf.1, buf.2)?;
        }
    }
    Ok(stats)
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
    let cfg = OrientCfg { expected: None, ..*base_cfg };
    let mut reader = parse_fastx_file(input)?;
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
fn build_metrics(
    stats: &Stats,
    whitelist: Option<&Whitelist>,
    demuxer: Option<&Demuxer>,
    t: &Timing,
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
        total_reads: stats.reads,
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
                100.0 * (stats.cb_on_whitelist + stats.cb_within1) as f64 / classified.max(1) as f64
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
    let final_ok = if has_demux { assigned } else { stats.pairs_emitted };
    let lost = stats.reads.saturating_sub(final_ok);
    let no_cdna = stats.class_partial.saturating_sub(stats.gated_out);
    let rps = if t.elapsed > 0.0 { (stats.reads as f64 / t.elapsed) as u64 } else { 0 };
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
