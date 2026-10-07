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
use super::process::{detect_orientation, process_one};
use super::types::{Batch, Stats, BATCH_READS};

/// Read-only context shared across the parallel per-file workers.
pub(crate) struct Shared<'a> {
    pub orient_cfg: OrientCfg,
    pub orient_mode: OrientMode,
    pub extract_cfg: ExtractCfg,
    pub demuxer: Option<&'a Demuxer>,
    pub gex_by_sample: &'a [bool],
    pub whitelist: Option<&'a Whitelist>,
    pub probe_wl: Option<&'a ProbeWhitelist>,
    pub cb_len: usize,
    pub umi_len: usize,
    pub cb_gate: bool,
    pub realign: bool,
    pub cb_reframe: bool,
    pub reframe_window: usize,
    pub lenient: bool,
    pub index_len: usize,
    pub index_window: isize,
    pub n_samples: usize,
    pub split: usize,
    pub minimap2: bool,
    pub processing_tags: bool,
    pub targets: &'a [PairShardWriter],
    pub failed: Option<&'a ShardWriter>,
    pub chunk_counter: &'a AtomicUsize,
}

pub(crate) fn file_orientation(input: &Path, sh: &Shared) -> Result<Option<Orientation>> {
    let expected = match sh.orient_mode {
        OrientMode::P5Led => Some(Orientation::P5Led),
        OrientMode::P7Led => Some(Orientation::P7Led),
        OrientMode::Auto => detect_orientation(input, &sh.orient_cfg, &sh.extract_cfg)?,
    };
    info!(file = %input.display(), orientation = ?expected, "orientation");
    Ok(expected)
}

pub(crate) fn process_file(input: &Path, sh: &Shared) -> Result<Stats> {
    let mut stats = Stats::new(sh.n_samples);
    let expected = file_orientation(input, sh)?;
    let orient_cfg = OrientCfg {
        expected,
        ..sh.orient_cfg
    };

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
                    return Ok(());
                }
            }
        }
        if batch.len() > 0 && tx.send(batch).is_err() {
            return Ok(());
        }
    }
    Ok(())
}

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

pub(crate) fn run_pipelined(
    inputs: &[PathBuf],
    sh: &Shared,
    n_readers: usize,
    n_workers: usize,
) -> Result<Vec<Stats>> {
    let (tx, rx) = bounded::<Batch>(n_workers * 2);
    let next = AtomicUsize::new(0);

    std::thread::scope(|scope| -> Result<Vec<Stats>> {
        let reader_handles: Vec<_> = (0..n_readers)
            .map(|_| {
                let tx = tx.clone();
                scope.spawn(|| reader_loop(inputs, &next, sh, tx))
            })
            .collect();
        drop(tx);

        let worker_handles: Vec<_> = (0..n_workers)
            .map(|w| {
                let rx = rx.clone();
                scope.spawn(move || worker_loop(w, sh, rx))
            })
            .collect();
        drop(rx);

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
