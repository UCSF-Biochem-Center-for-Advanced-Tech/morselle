use std::collections::HashMap;
use std::path::Path;
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc,
};
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read};
use tracing::info;

use super::helpers::{
    aligned_bases, aux_string, better, choose_gene, intern_sample, read_hash,
    record_identity, shard_hash,
};
use super::types::{Evidence, Metrics, RawKey, Dna, NO_SAMPLE};

pub fn collate<P: AsRef<Path>>(inputs: &[P], threads: usize, shard_paths: &[P]) -> Result<()> {
    let first = inputs
        .first()
        .context("collate() requires at least one input BAM")?
        .as_ref();

    let header = {
        let reader = bam::Reader::from_path(first)
            .with_context(|| format!("opening {}", first.display()))?;
        bam::Header::from_template(reader.header())
    };
    // Uncompressed: shard files are pure ephemeral intermediates, deleted
    // right after quant reads each one exactly once -- same lesson as
    // extract's --no-gzip finding (PROFILING.md: 53.9% of extract's CPU was
    // zlib, --no-gzip gave an 8x win). Collate was previously the dominant
    // cost of the whole --shards run for exactly this reason: 8 writers all
    // doing default-level BGZF compression, serialized onto the single
    // collate loop thread.
    let mut writers: Vec<bam::Writer> = shard_paths
        .iter()
        .map(|path| {
            let path = path.as_ref();
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)
                .with_context(|| format!("creating shard {}", path.display()))?;
            writer.set_compression_level(bam::CompressionLevel::Uncompressed)?;
            Ok(writer)
        })
        .collect::<Result<Vec<_>>>()?;
    let shards = writers.len();
    let mut routed = 0_u64;
    for input in inputs {
        let input = input.as_ref();
        let mut reader = bam::Reader::from_path(input)
            .with_context(|| format!("opening {}", input.display()))?;

        reader.set_threads(threads)?;
        for result in reader.records() {
            let record = result
                .with_context(|| format!("reading BAM in collate pass ({})", input.display()))?;
            let cb = aux_string(&record, b"CB");
            let pb = aux_string(&record, b"PB");
            let idx = match &cb {
                Some(cb) => {
                    shard_hash(cb.as_bytes(), pb.as_deref().map(str::as_bytes)) as usize % shards
                }
                None => 0,
            };
            writers[idx].write(&record)?;
            routed += 1;
            if routed % 25_000_000 == 0 {
                eprintln!("[sbxpress molecules] phase=collate records={routed}");
            }
        }
    }
    eprintln!(
        "[sbxpress molecules] phase=collate complete records={routed} inputs={}",
        inputs.len()
    );
    Ok(())
}

#[allow(clippy::type_complexity)]
pub fn collect_evidence(
    path: &Path,
    threads: usize,
    require_pb: bool,
    metrics: &mut Metrics,
) -> Result<(
    HashMap<RawKey, Evidence>,
    Vec<String>,
    HashMap<String, u32>,
    Vec<String>,
    HashMap<String, u32>,
)> {
    let mut reader =
        bam::Reader::from_path(path).with_context(|| format!("opening {}", path.display()))?;
    reader.set_threads(threads)?;
    let mut evidence = HashMap::<RawKey, Evidence>::new();
    let mut genes = Vec::<String>::new();
    let mut gene_ids = HashMap::<String, u32>::new();
    let mut samples = Vec::<String>::new();
    let mut sample_ids = HashMap::<String, u32>::new();
    eprintln!("[sbxpress molecules] phase=evidence_collect mode=serial start");
    for (ordinal, result) in reader.records().enumerate() {
        let record = result.context("reading BAM in evidence pass")?;
        metrics.records += 1;
        if metrics.records % 25_000_000 == 0 {
            eprintln!(
                "[sbxpress molecules] phase=evidence_collect mode=serial records={} eligible={}",
                metrics.records, metrics.eligible_records
            );
        }
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }
        let Some(cell) = aux_string(&record, b"CB").and_then(|v| Dna::pack(&v)) else {
            continue;
        };
        let Some(umi) = aux_string(&record, b"UR")
            .or_else(|| aux_string(&record, b"UB"))
            .and_then(|v| Dna::pack(&v))
        else {
            continue;
        };
        let Some((gene, _name)) = choose_gene(&record) else {
            continue;
        };
        // PB is a Flex-only identity component. Tag-seq (and ordinary
        // non-Flex assays) deliberately use CB alone and must not be gated on
        // a PB tag that their chemistry never produces.
        let sample = if let Some(pb) = aux_string(&record, b"PB").filter(|v| !v.is_empty()) {
            intern_sample(&pb, &mut samples, &mut sample_ids)
        } else if !require_pb {
            NO_SAMPLE
        } else {
            metrics.no_sample_records += 1;
            continue;
        };
        let gene = if let Some(&id) = gene_ids.get(&gene) {
            id
        } else {
            let id = genes.len() as u32;
            genes.push(gene.clone());
            gene_ids.insert(gene, id);
            id
        };
        metrics.eligible_records += 1;
        let key = RawKey {
            cell,
            sample,
            gene,
            umi,
        };
        let ordinal = ordinal as u64;
        evidence
            .entry(key)
            .and_modify(|item| {
                item.reads += 1;
                item.min_read_hash = item.min_read_hash.min(read_hash(record.qname()));
                if better(&record, ordinal, item) {
                    item.best_ordinal = ordinal;
                    item.best_mapq = record.mapq();
                    item.best_aligned = aligned_bases(&record);
                }
            })
            .or_insert(Evidence {
                reads: 1,
                best_ordinal: ordinal,
                best_mapq: record.mapq(),
                best_aligned: aligned_bases(&record),
                min_read_hash: read_hash(record.qname()),
            });
    }
    eprintln!(
        "[sbxpress molecules] phase=evidence_collect mode=serial complete records={} eligible={} raw_nodes={}",
        metrics.records,
        metrics.eligible_records,
        evidence.len()
    );
    Ok((evidence, genes, gene_ids, samples, sample_ids))
}

pub struct RegionEvidence {
    pub evidence: HashMap<RawKey, Evidence>,
    pub genes: Vec<String>,
    pub samples: Vec<String>,
    pub records: u64,
    pub eligible_records: u64,
    pub no_sample_records: u64,
}

pub fn merge_item(target: &mut Evidence, source: Evidence) {
    target.reads += source.reads;
    target.min_read_hash = target.min_read_hash.min(source.min_read_hash);
    if source.best_mapq > target.best_mapq
        || (source.best_mapq == target.best_mapq
            && (source.best_aligned > target.best_aligned
                || (source.best_aligned == target.best_aligned
                    && source.best_ordinal < target.best_ordinal)))
    {
        target.best_ordinal = source.best_ordinal;
        target.best_mapq = source.best_mapq;
        target.best_aligned = source.best_aligned;
    }
}

pub fn merge_evidence_maps(
    mut left: HashMap<RawKey, Evidence>,
    mut right: HashMap<RawKey, Evidence>,
) -> HashMap<RawKey, Evidence> {
    // Keeping the larger allocation avoids repeated hash-table growth during
    // the parallel reduction tree.
    if left.len() < right.len() {
        std::mem::swap(&mut left, &mut right);
    }
    for (key, item) in right {
        left.entry(key)
            .and_modify(|current| merge_item(current, item))
            .or_insert(item);
    }
    left
}

pub fn collect_region(
    path: &Path,
    tid: u32,
    start: u64,
    end: u64,
    require_pb: bool,
) -> Result<RegionEvidence> {
    let mut reader = bam::IndexedReader::from_path(path)
        .with_context(|| format!("opening indexed {}", path.display()))?;
    reader
        .fetch((tid, start, end))
        .with_context(|| format!("fetching tid {tid}:{start}-{end}"))?;
    let mut evidence = HashMap::<RawKey, Evidence>::new();
    let mut genes = Vec::<String>::new();
    let mut gene_ids = HashMap::<String, u32>::new();
    let mut samples = Vec::<String>::new();
    let mut sample_ids = HashMap::<String, u32>::new();
    let mut records = 0_u64;
    let mut eligible_records = 0_u64;
    let mut no_sample_records = 0_u64;
    for result in reader.records() {
        let record = result.context("reading BAM region in evidence pass")?;
        // An indexed query can return an alignment overlapping its left edge.
        // Assign by alignment start so adjacent work items never double-count.
        if record.tid() != tid as i32 || record.pos() < start as i64 || record.pos() >= end as i64 {
            continue;
        }
        records += 1;
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }
        let Some(cell) = aux_string(&record, b"CB").and_then(|v| Dna::pack(&v)) else {
            continue;
        };
        let Some(umi) = aux_string(&record, b"UR")
            .or_else(|| aux_string(&record, b"UB"))
            .and_then(|v| Dna::pack(&v))
        else {
            continue;
        };
        let Some((gene, _name)) = choose_gene(&record) else {
            continue;
        };
        let sample = if let Some(pb) = aux_string(&record, b"PB").filter(|v| !v.is_empty()) {
            intern_sample(&pb, &mut samples, &mut sample_ids)
        } else if !require_pb {
            NO_SAMPLE
        } else {
            no_sample_records += 1;
            continue;
        };
        let gene = if let Some(&id) = gene_ids.get(&gene) {
            id
        } else {
            let id = genes.len() as u32;
            genes.push(gene.clone());
            gene_ids.insert(gene, id);
            id
        };
        eligible_records += 1;
        let key = RawKey {
            cell,
            sample,
            gene,
            umi,
        };
        let candidate = Evidence {
            reads: 1,
            best_ordinal: record_identity(&record),
            best_mapq: record.mapq(),
            best_aligned: aligned_bases(&record),
            min_read_hash: read_hash(record.qname()),
        };
        evidence
            .entry(key)
            .and_modify(|item| merge_item(item, candidate))
            .or_insert(candidate);
    }
    Ok(RegionEvidence {
        evidence,
        genes,
        samples,
        records,
        eligible_records,
        no_sample_records,
    })
}

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

#[allow(clippy::type_complexity)]
pub fn collect_evidence_indexed(
    path: &Path,
    workers: usize,
    chunk_bases: u64,
    require_pb: bool,
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
                let result = collect_region(path, tid, start, end, require_pb);
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
    let mut gene_ids = HashMap::<String, u32>::new();
    let mut samples = Vec::<String>::new();
    let mut sample_ids = HashMap::<String, u32>::new();
    let mut remaps = Vec::with_capacity(parts.len());
    let mut sample_remaps = Vec::with_capacity(parts.len());
    for part in &parts {
        metrics.records += part.records;
        metrics.eligible_records += part.eligible_records;
        metrics.no_sample_records += part.no_sample_records;
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
                let mut map = HashMap::with_capacity(part.evidence.len());
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
