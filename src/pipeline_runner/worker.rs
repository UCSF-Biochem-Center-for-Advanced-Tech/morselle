use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use anyhow::Result;
use crossbeam_channel::{bounded, Receiver, Sender};
use tracing::info;

use crate::barcode::Whitelist;
use crate::demux::Demuxer;
use crate::extract::ExtractCfg;
use crate::gzin::open_fastx;
use crate::io::{PairShardWriter, ShardWriter};
use crate::orient::{OrientCfg, OrientMode};
use crate::pipeline::Orientation;
use crate::probe_barcode::ProbeWhitelist;

use super::formatting::{flush_failed_buf, flush_target_bufs, make_target_bufs};
use super::orientation_scan::detect_orientation;
use super::process::process_one;
use super::types::{Batch, Stats, BATCH_READS};

/// Read-only context shared across the parallel per-file workers.
pub(crate) struct Shared<'a> {
    pub orient_cfg: OrientCfg,
    pub orient_mode: OrientMode,
    pub extract_cfg: ExtractCfg,
    pub demuxer: Option<&'a Demuxer>,
    /// Per-sample (indexed like `demuxer`'s `sample_ids()`) "does this
    /// sample's declared name look like a Flex GEX library" -- the routing
    /// signal for `ExtractCfg::flex_r2_boundary`'s hard-gate default. Empty
    /// when `demuxer` is `None`.
    pub gex_by_sample: &'a [bool],
    pub whitelist: Option<&'a Whitelist>,
    /// 10x Flex probe-barcode whitelist (the Flex design notes Phase 1). `None`
    /// unless `--probe-barcode-whitelist` was given — zero cost otherwise.
    pub probe_wl: Option<&'a ProbeWhitelist>,
    pub cb_len: usize,
    pub umi_len: usize,
    pub cb_gate: bool,
    /// Per-barcode counters (`--barcode-metrics`); `None` skips them.
    pub barcodes: Option<&'a crate::barcode_metrics::SharedBarcodeTable<10>>,
    /// Low-complexity cDNA check (`--low-complexity-*`).
    pub complexity: crate::complexity::ComplexityCfg,
    pub realign: bool,
    pub cb_reframe: bool,
    pub reframe_window: usize,
    pub lenient: bool,
    pub index_len: usize,
    /// Max frame shift tried when matching i5/i7 (absorbs index-shifting indels).
    pub index_window: isize,
    /// Number of sample targets; the Undetermined target (when demuxing) is at
    /// index `n_samples`.
    pub n_samples: usize,
    /// Output chunks per logical target (`--split`). Physical writer for logical
    /// target `t`, chunk `c` is `targets[t * split + c]`.
    pub split: usize,
    /// Emit only R2 FASTQs with 10x-style barcode/UMI SAM tags in the
    /// FASTQ comment, for minimap2 `-y` pass-through.
    pub minimap2: bool,
    /// Add SBX processing provenance tags to R2 FASTQ comments.
    pub processing_tags: bool,
    pub targets: &'a [PairShardWriter],
    /// Side channel for reads that never produce a pair. See `--failed-output`.
    pub failed: Option<&'a ShardWriter>,
    /// Round-robin source for `process_file`'s per-file chunk assignment (the
    /// "many files" path). Shared and atomic, not per-thread: real per-file
    /// yield varies a lot (some raw shards pass far more reads than others),
    /// and rayon's work-stealing doesn't guarantee an even file count per
    /// thread either, so keying chunk-choice off which thread happened to
    /// claim a file let a handful of high-yield files landing on the same
    /// thread permanently skew that one chunk file to several times the
    /// median size. A shared counter, incremented once per file regardless
    /// of which thread claims it, decorrelates chunk assignment from any
    /// thread-level yield clustering.
    pub chunk_counter: &'a AtomicUsize,
}

/// Resolve the per-file orientation to enforce (detecting it for `Auto`).
pub(crate) fn file_orientation(input: &Path, sh: &Shared) -> Result<Option<Orientation>> {
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
pub(crate) fn process_file(input: &Path, sh: &Shared) -> Result<Stats> {
    let mut stats = Stats::new(sh.n_samples);
    let expected = file_orientation(input, sh)?;
    let orient_cfg = OrientCfg {
        expected,
        ..sh.orient_cfg
    };

    // One chunk per logical target for this whole file (parallel per-sample
    // compression). Chosen via a shared round-robin counter, not this
    // thread's identity: real per-file yield varies substantially, and
    // tying chunk choice to the processing thread let a handful of
    // high-yield files that happened to land on the same thread (via
    // rayon's work-stealing) permanently skew that thread's chunk to
    // several times the median size (real CAT06 full run: 2 of 128 chunks
    // came out ~4.4x the median, becoming the long-pole stragglers for
    // every downstream per-chunk consumer).
    let chunk = sh.chunk_counter.fetch_add(1, Ordering::Relaxed) % sh.split;
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
pub(crate) fn run_pipelined(
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
