use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};
use anyhow::Result;
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, Read};

use super::dedup::process_and_write_shard;
use super::evidence::collate;
use super::helpers::{aux_string, choose_gene, record_identity};
use super::types::{
    Decision, Dna, Metrics, Outcome, RawKey, NO_SAMPLE, merge_metrics,
};

pub fn run_sharded(
    inputs: &[PathBuf],
    output: &Path,
    threads: usize,
    shards: usize,
    shard_concurrency: usize,
    evidence_workers: usize,
    tagseq: bool,
    require_pb: bool,
    metrics: &mut Metrics,
) -> Result<()> {
    let shard_dir = tempfile::tempdir()?;
    let shard_in_paths: Vec<PathBuf> = (0..shards)
        .map(|i| shard_dir.path().join(format!("shard_in_{i:04}.bam")))
        .collect();
    let shard_out_paths: Vec<PathBuf> = (0..shards)
        .map(|i| shard_dir.path().join(format!("shard_out_{i:04}.bam")))
        .collect();

    let inputs_path: Vec<&Path> = inputs.iter().map(Path::new).collect();
    let shard_in_refs: Vec<&Path> = shard_in_paths.iter().map(|p| p.as_path()).collect();
    collate(&inputs_path, threads, &shard_in_refs)?;

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(shard_concurrency.max(1))
        .build()?;

    let shard_jobs: Vec<(&PathBuf, &PathBuf)> = shard_in_paths
        .iter()
        .zip(&shard_out_paths)
        .collect();

    let results = pool.install(|| {
        shard_jobs
            .par_iter()
            .map(|(in_path, out_path)| {
                process_and_write_shard(
                    in_path,
                    out_path,
                    1,
                    require_pb,
                    None,
                    evidence_workers,
                    tagseq,
                )
            })
            .collect::<Result<Vec<_>>>()
    })?;

    for (m, _sat) in results {
        merge_metrics(metrics, &m);
    }

    concatenate_bam_files(&shard_out_paths, output, threads)?;
    Ok(())
}

pub fn write_shard_output(
    input_path: &Path,
    output_path: &Path,
    threads: usize,
    require_pb: bool,
    decisions: &HashMap<RawKey, Decision>,
    gene_ids: &HashMap<String, u32>,
    sample_ids: &HashMap<String, u32>,
    metrics: &mut Metrics,
) -> Result<()> {
    if let Some(parent) = output_path.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut reader = bam::Reader::from_path(input_path)?;
    reader.set_threads(threads)?;

    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;
    writer.set_compression_level(bam::CompressionLevel::Uncompressed)?;
    writer.set_threads(threads)?;

    for rec in reader.records() {
        let mut record = rec?;
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            writer.write(&record)?;
            continue;
        }
        let Some(cell) = aux_string(&record, b"CB").and_then(|v| Dna::pack(&v)) else {
            writer.write(&record)?;
            continue;
        };
        let Some(umi) = aux_string(&record, b"UR")
            .or_else(|| aux_string(&record, b"UB"))
            .and_then(|v| Dna::pack(&v))
        else {
            writer.write(&record)?;
            continue;
        };
        let Some((gene_str, _)) = choose_gene(&record) else {
            writer.write(&record)?;
            continue;
        };
        let Some(&gene_id) = gene_ids.get(&gene_str) else {
            writer.write(&record)?;
            continue;
        };
        let sample_id = if let Some(pb) = aux_string(&record, b"PB").filter(|v| !v.is_empty()) {
            sample_ids.get(&pb).copied().unwrap_or(NO_SAMPLE)
        } else if !require_pb {
            NO_SAMPLE
        } else {
            writer.write(&record)?;
            continue;
        };

        let key = RawKey {
            cell,
            sample: sample_id,
            gene: gene_id,
            umi,
        };

        if let Some(dec) = decisions.get(&key) {
            let is_rep = record_identity(&record) == dec.representative;
            let tag_val = match (dec.outcome, is_rep) {
                (Outcome::Won | Outcome::Unique, true) => "won",
                (Outcome::Won | Outcome::Unique, false) => "lost",
                (Outcome::Lost, _) => "lost",
                (Outcome::Tie, _) => "tie",
                (Outcome::Shared, true) => "shared",
                (Outcome::Shared, false) => "lost",
            };
            record.push_aux(b"ZG", Aux::String(tag_val))?;
            let corr_str = dec.corrected.unpack();
            record.push_aux(b"UB", Aux::String(&corr_str))?;
            if is_rep && matches!(dec.outcome, Outcome::Won | Outcome::Unique | Outcome::Shared) {
                metrics.representatives_tagged += 1;
            }
        }
        writer.write(&record)?;
    }
    Ok(())
}

pub fn concatenate_bam_files(
    shard_paths: &[PathBuf],
    output: &Path,
    threads: usize,
) -> Result<()> {
    if shard_paths.is_empty() {
        anyhow::bail!("no shard files to concatenate");
    }
    let first = &shard_paths[0];
    let reader = bam::Reader::from_path(first)?;
    let header = bam::Header::from_template(reader.header());

    if let Some(parent) = output.parent() {
        fs::create_dir_all(parent)?;
    }
    let mut writer = bam::Writer::from_path(output, &header, bam::Format::Bam)?;
    writer.set_threads(threads)?;

    for path in shard_paths {
        let mut r = bam::Reader::from_path(path)?;
        r.set_threads(threads)?;
        for result in r.records() {
            let record = result?;
            writer.write(&record)?;
        }
    }
    Ok(())
}

pub fn write_output(
    inputs: &[PathBuf],
    output: &Path,
    threads: usize,
    _workers: usize,
    require_pb: bool,
    decisions: &HashMap<RawKey, Decision>,
    _genes: &[String],
    gene_ids: &HashMap<String, u32>,
    sample_ids: &HashMap<String, u32>,
    metrics: &mut Metrics,
) -> Result<()> {
    write_shard_output(
        &inputs[0],
        output,
        threads,
        require_pb,
        decisions,
        gene_ids,
        sample_ids,
        metrics,
    )
}
