//! BAM output for `sbx10x-molecules`: tag every eligible record with its
//! corrected UMI (`UB`), representative flag (`xf:25` counted / `xf:17`
//! not), collision degree (`XV`), outcome (`XW`), and, for fractional
//! molecules, the matrix weight (`ZW:f`).
//!
//! Restored from the pre-0.5.7 implementation (commit `4d42a18`). `ZG` is
//! never written here: it belongs to `annotate` (terminal gene call used to
//! resolve multi-gene `GX`) and is read back by `choose_gene`.
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use anyhow::{Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, Read};

use super::dedup::build_decisions;
use super::evidence::{collate, collect_evidence};
use super::helpers::{aux_string, choose_gene, record_identity};
use super::types::{merge_metrics, Decision, Dna, Metrics, RawKey, NO_SAMPLE};

/// Settings shared by the unsharded and sharded write paths.
#[derive(Clone, Copy, Debug)]
pub struct MoleculeOptions {
    pub threads: usize,
    pub evidence_workers: usize,
    pub flex: bool,
    pub tagseq: bool,
    pub fractional: bool,
}

/// Lookup tables produced by evidence collection, needed to rebuild a
/// record's `RawKey` at write time.
pub struct KeyTables<'a> {
    pub genes: &'a [String],
    pub gene_ids: &'a HashMap<String, u32>,
    pub sample_ids: &'a HashMap<String, u32>,
}

fn replace_string(record: &mut bam::Record, tag: &[u8; 2], value: &str) -> Result<()> {
    let _ = record.remove_aux(tag);
    record.push_aux(tag, Aux::String(value))?;
    Ok(())
}

fn record_key(record: &bam::Record, flex: bool, tables: &KeyTables) -> Option<RawKey> {
    if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
        return None;
    }
    let cell = aux_string(record, b"CB").and_then(|v| Dna::pack(&v))?;
    let umi = aux_string(record, b"UR")
        .or_else(|| aux_string(record, b"UB"))
        .and_then(|v| Dna::pack(&v))?;
    let gene = choose_gene(record).and_then(|(g, _)| tables.gene_ids.get(&g).copied())?;
    let sample = if flex {
        aux_string(record, b"PB")
            .filter(|v| !v.is_empty())
            .and_then(|pb| tables.sample_ids.get(&pb).copied())?
    } else {
        NO_SAMPLE
    };
    Some(RawKey { cell, sample, gene, umi })
}

/// Rewrite one record's molecule tags in place. `identity` is the value the
/// evidence pass stored as `best_ordinal` for this record (positional
/// ordinal in serial mode, `record_identity()` in indexed mode).
fn tag_record(
    record: &mut bam::Record,
    identity: impl Fn(&bam::Record) -> u64,
    flex: bool,
    decisions: &HashMap<RawKey, Decision>,
    tables: &KeyTables,
    metrics: &mut Metrics,
) -> Result<()> {
    let key = record_key(record, flex, tables);
    for tag in [b"xf", b"XV", b"XW", b"ZW"] {
        let _ = record.remove_aux(tag);
    }
    let Some((key, decision)) = key.and_then(|key| decisions.get(&key).map(|d| (key, d))) else {
        return Ok(());
    };
    replace_string(record, b"UB", &decision.corrected.unpack())?;
    // A unique ZG terminal call can resolve a multi-gene GX. Materialize
    // that decision so counted representatives satisfy the single-GX
    // invariant expected by cr-bam-to-mtx.
    if aux_string(record, b"GX").is_some_and(|gx| gx.contains(';')) {
        let selected = &tables.genes[key.gene as usize];
        let selected_name = choose_gene(record)
            .map(|(_, name)| name)
            .unwrap_or_else(|| selected.clone());
        replace_string(record, b"GX", selected)?;
        replace_string(record, b"GN", &selected_name)?;
    }
    let representative =
        decision.representative == identity(record) && decision.outcome.is_counted();
    record.push_aux(b"xf", Aux::I32(if representative { 25 } else { 17 }))?;
    record.push_aux(b"XV", Aux::U32(decision.collision_degree))?;
    record.push_aux(b"XW", Aux::Char(decision.outcome.xw_code()))?;
    if decision.outcome == super::types::Outcome::Fractional {
        record.push_aux(b"ZW", Aux::Float(decision.weight))?;
    }
    if representative {
        metrics.representatives_tagged += 1;
    }
    Ok(())
}

/// Unsharded write: re-read the single input BAM in order and tag it.
pub fn write_output(
    input: &Path,
    output: &Path,
    opts: &MoleculeOptions,
    decisions: &HashMap<RawKey, Decision>,
    tables: &KeyTables,
    metrics: &mut Metrics,
) -> Result<()> {
    eprintln!("[sbxpress molecules] phase=bam_write start");
    let mut reader = bam::Reader::from_path(input)?;
    reader.set_threads(opts.threads)?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;
    writer.set_threads(opts.threads)?;
    let indexed = opts.evidence_workers > 1;
    let mut written = 0_u64;
    for (ordinal, result) in reader.records().enumerate() {
        let mut record = result?;
        let identity = |r: &bam::Record| if indexed { record_identity(r) } else { ordinal as u64 };
        tag_record(&mut record, identity, opts.flex, decisions, tables, metrics)?;
        writer.write(&record)?;
        written += 1;
        if written % 25_000_000 == 0 {
            eprintln!(
                "[sbxpress molecules] phase=bam_write records={} representatives={}",
                written, metrics.representatives_tagged
            );
        }
    }
    eprintln!(
        "[sbxpress molecules] phase=bam_write complete records={} representatives={}",
        written, metrics.representatives_tagged
    );
    Ok(())
}

/// Sharded write: one shard's own temporary BAM, appended into the shared
/// output writer. Always ordinal-based: a shard is built via the serial
/// `collect_evidence()`, so ordinals are stable within this sequential read.
fn write_shard_output(
    shard_path: &Path,
    opts: &MoleculeOptions,
    decisions: &HashMap<RawKey, Decision>,
    tables: &KeyTables,
    writer: &mut bam::Writer,
    metrics: &mut Metrics,
) -> Result<()> {
    let mut reader = bam::Reader::from_path(shard_path)?;
    reader.set_threads(opts.threads)?;
    for (ordinal, result) in reader.records().enumerate() {
        let mut record = result?;
        tag_record(&mut record, |_| ordinal as u64, opts.flex, decisions, tables, metrics)?;
        writer.write(&record)?;
    }
    Ok(())
}

/// Full sharded (collate + quant) path: bounds peak memory to roughly
/// `total/shards` by processing one `(cell, sample)`-hashed shard fully
/// (evidence -> correct -> aggregate -> decide -> write) before discarding
/// it. Shard files live in `<output>.shards/` (same filesystem as the
/// output, not the system temp dir) and are removed afterwards.
pub fn run_sharded(
    inputs: &[PathBuf],
    output: &Path,
    shards: usize,
    shard_concurrency: usize,
    opts: &MoleculeOptions,
    metrics: &mut Metrics,
) -> Result<()> {
    let shard_dir = PathBuf::from(format!("{}.shards", output.display()));
    fs::create_dir_all(&shard_dir)
        .with_context(|| format!("creating shard directory {}", shard_dir.display()))?;
    let shard_paths: Vec<PathBuf> = (0..shards)
        .map(|i| shard_dir.join(format!("shard_{i:04}.bam")))
        .collect();
    eprintln!("[sbxpress molecules] phase=collate start shards={shards}");
    collate(inputs, opts.threads, &shard_paths)?;

    let header = {
        let reader = bam::Reader::from_path(&inputs[0])?;
        bam::Header::from_template(reader.header())
    };
    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;
    writer.set_threads(opts.threads)?;
    let writer = Mutex::new(writer);

    // Shards are fully independent until the final write (proof:
    // sharding_by_cell_sample_matches_global_build_decisions); only the
    // write holds the output lock.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(shard_concurrency.min(shard_paths.len()))
        .build()
        .context("creating shard-concurrency worker pool")?;
    let total_shards = shard_paths.len();
    let completed = Arc::new(AtomicUsize::new(0));
    let shard_metrics: Vec<Metrics> = pool.install(|| {
        shard_paths
            .par_iter()
            .map(|shard_path| -> Result<Metrics> {
                let mut local = Metrics::default();
                let (evidence, genes, gene_ids, _samples, sample_ids) =
                    collect_evidence(shard_path, opts.threads, opts.flex, &mut local)?;
                let (decisions, _) = build_decisions(
                    &evidence,
                    &mut local,
                    None,
                    opts.evidence_workers,
                    opts.tagseq,
                    opts.fractional,
                );
                drop(evidence);
                let tables = KeyTables { genes: &genes, gene_ids: &gene_ids, sample_ids: &sample_ids };
                {
                    let mut guard = writer.lock().expect("output writer mutex poisoned");
                    write_shard_output(shard_path, opts, &decisions, &tables, &mut guard, &mut local)?;
                }
                let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                eprintln!("[sbxpress molecules] phase=shard_quant complete shards={done}/{total_shards}");
                Ok(local)
            })
            .collect::<Result<Vec<_>>>()
    })?;
    for local in &shard_metrics {
        merge_metrics(metrics, local);
    }
    drop(writer);
    if let Err(err) = fs::remove_dir_all(&shard_dir) {
        eprintln!(
            "[sbxpress molecules] warning: failed to clean up shard directory {}: {err}",
            shard_dir.display()
        );
    }
    Ok(())
}
