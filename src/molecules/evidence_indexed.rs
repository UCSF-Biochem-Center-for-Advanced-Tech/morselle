//! Parallel evidence collection over indexed BAM regions
//! (`--evidence-workers` > 1). Split from `evidence.rs` for size.
// FxHash: these maps hold billions of fixed-size integer keys at real
// scale; SipHash's DoS resistance buys nothing here and cost ~6% of runtime.
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use rustc_hash::FxHashMap as HashMap;
use std::path::Path;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc,
};
use tracing::info;

use super::evidence::{
    collect_region_with, merge_evidence_maps, merge_item, RecordFilter, RegionEvidence,
};
use super::types::{Evidence, Metrics, RawKey, NO_SAMPLE};

pub fn indexed_regions(path: &Path, chunk_bases: u64) -> Result<Vec<(u32, u64, u64)>> {
    let reader = bam::IndexedReader::from_path(path)
        .with_context(|| format!("opening indexed {}", path.display()))?;
    let header = reader.header();
    let mut regions = Vec::new();
    for tid in 0..header.target_count() {
        let length = header.target_len(tid).unwrap_or(0);
        for start in (0..length).step_by(chunk_bases as usize) {
            regions.push((tid, start, (start + chunk_bases).min(length)));
        }
    }
    Ok(regions)
}

/// [`collect_evidence_indexed_with`] without the mapping-quality gate.
#[allow(clippy::type_complexity)]
pub fn collect_evidence_indexed(
    path: &Path,
    workers: usize,
    chunk_bases: u64,
    require_pb: bool,
    position_aware: bool,
    metrics: &mut Metrics,
) -> Result<(
    HashMap<RawKey, Evidence>,
    Vec<String>,
    HashMap<String, u32>,
    Vec<String>,
    HashMap<String, u32>,
)> {
    let filter = RecordFilter::unfiltered(require_pb, position_aware);
    collect_evidence_indexed_with(path, workers, chunk_bases, &filter, metrics)
}

#[allow(clippy::type_complexity)]
pub fn collect_evidence_indexed_with(
    path: &Path,
    workers: usize,
    chunk_bases: u64,
    filter: &RecordFilter,
    metrics: &mut Metrics,
) -> Result<(
    HashMap<RawKey, Evidence>,
    Vec<String>,
    HashMap<String, u32>,
    Vec<String>,
    HashMap<String, u32>,
)> {
    if chunk_bases == 0 {
        bail!("--evidence-chunk-bases must be greater than zero");
    }
    let regions = indexed_regions(path, chunk_bases)?;
    info!(
        workers,
        regions = regions.len(),
        chunk_bases,
        "dispatching indexed evidence regions"
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()?;
    eprintln!(
        "[sbxpress molecules] phase=evidence_collect mode=indexed start workers={} regions={} chunk_bases={}",
        workers,
        regions.len(),
        chunk_bases
    );
    let completed = Arc::new(AtomicUsize::new(0));
    let checkpoint = (regions.len() / 100).max(1);
    let total_regions = regions.len();
    let parts: Vec<RegionEvidence> = pool.install(|| {
        regions
            .par_iter()
            .map(|&(tid, start, end)| {
                let result = collect_region_with(path, tid, start, end, filter);
                let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                if done % checkpoint == 0 || done == total_regions {
                    eprintln!(
                        "[sbxpress molecules] phase=evidence_collect mode=indexed regions={done}/{total_regions}"
                    );
                }
                result
            })
            .collect::<Result<Vec<_>>>()
    })?;
    eprintln!("[sbxpress molecules] phase=evidence_collect mode=indexed complete");
    let mut genes = Vec::<String>::new();
    let mut gene_ids = HashMap::<String, u32>::default();
    let mut samples = Vec::<String>::new();
    let mut sample_ids = HashMap::<String, u32>::default();
    let mut remaps = Vec::with_capacity(parts.len());
    let mut sample_remaps = Vec::with_capacity(parts.len());
    for part in &parts {
        metrics.records += part.records;
        metrics.eligible_records += part.eligible_records;
        metrics.no_sample_records += part.no_sample_records;
        part.gate_counts.add_to(metrics);
        let remap: Vec<u32> = part
            .genes
            .iter()
            .map(|gene| {
                if let Some(&id) = gene_ids.get(gene) {
                    id
                } else {
                    let id = genes.len() as u32;
                    gene_ids.insert(gene.clone(), id);
                    genes.push(gene.clone());
                    id
                }
            })
            .collect();
        remaps.push(remap);
        let sample_remap: Vec<u32> = part
            .samples
            .iter()
            .map(|sample| {
                if let Some(&id) = sample_ids.get(sample) {
                    id
                } else {
                    let id = samples.len() as u32;
                    sample_ids.insert(sample.clone(), id);
                    samples.push(sample.clone());
                    id
                }
            })
            .collect();
        sample_remaps.push(sample_remap);
        super::barcode_metrics::merge_region(metrics, part, sample_remaps.last().unwrap());
    }
    eprintln!(
        "[sbxpress molecules] phase=evidence_remap start regions={} records={} eligible={}",
        parts.len(),
        metrics.records,
        metrics.eligible_records
    );
    let maps: Vec<HashMap<RawKey, Evidence>> = pool.install(|| {
        parts
            .into_par_iter()
            .zip(remaps.into_par_iter())
            .zip(sample_remaps.into_par_iter())
            .map(|((part, remap), sample_remap)| {
                let mut map =
                    HashMap::with_capacity_and_hasher(part.evidence.len(), Default::default());
                for (mut key, item) in part.evidence {
                    key.gene = remap[key.gene as usize];
                    if key.sample != NO_SAMPLE {
                        key.sample = sample_remap[key.sample as usize];
                    }
                    map.entry(key)
                        .and_modify(|current| merge_item(current, item))
                        .or_insert(item);
                }
                map
            })
            .collect()
    });
    eprintln!("[sbxpress molecules] phase=evidence_remap complete; phase=evidence_merge start");
    let evidence = pool.install(|| {
        maps.into_par_iter()
            .reduce_with(merge_evidence_maps)
            .unwrap_or_default()
    });
    eprintln!(
        "[sbxpress molecules] phase=evidence_merge complete raw_nodes={}",
        evidence.len()
    );
    Ok((evidence, genes, gene_ids, samples, sample_ids))
}
