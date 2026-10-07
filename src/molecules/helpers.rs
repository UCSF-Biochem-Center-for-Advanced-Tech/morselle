use std::collections::HashMap;
use rust_htslib::bam::{self, record::Aux};
use super::types::Evidence;

/// Stable FNV-1a followed by SplitMix64.
pub fn read_hash(name: &[u8]) -> u32 {
    let mut value = 0xcbf29ce484222325u64;
    for &byte in name {
        value ^= byte as u64;
        value = value.wrapping_mul(0x100000001b3);
    }
    let mut z = value.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    (z ^ (z >> 31)) as u32
}

pub fn aux_string(record: &bam::Record, tag: &[u8; 2]) -> Option<String> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Some(value.to_owned()),
        _ => None,
    }
}

pub fn choose_gene(record: &bam::Record) -> Option<(String, String)> {
    let gx = aux_string(record, b"GX")?;
    let gn = aux_string(record, b"GN").unwrap_or_else(|| gx.clone());
    let genes: Vec<&str> = gx.split(';').filter(|v| !v.is_empty()).collect();
    let names: Vec<&str> = gn.split(';').collect();
    if genes.len() == 1 {
        return Some((
            genes[0].to_owned(),
            names.first().copied().unwrap_or(genes[0]).to_owned(),
        ));
    }
    // For a 3' run, the annotator's unique terminal-supported gene resolves a
    // genomic overlap without discarding the alignment.
    let zg = aux_string(record, b"ZG")?;
    if zg.contains(';') || zg.is_empty() {
        return None;
    }
    let index = genes.iter().position(|gene| *gene == zg)?;
    Some((
        zg,
        names.get(index).copied().unwrap_or(genes[index]).to_owned(),
    ))
}

pub fn aligned_bases(record: &bam::Record) -> u32 {
    record
        .cigar()
        .iter()
        .map(|op| match op.char() {
            'M' | '=' | 'X' => op.len(),
            _ => 0,
        })
        .sum()
}

pub fn better(record: &bam::Record, ordinal: u64, evidence: &Evidence) -> bool {
    record.mapq() > evidence.best_mapq
        || (record.mapq() == evidence.best_mapq
            && (aligned_bases(record) > evidence.best_aligned
                || (aligned_bases(record) == evidence.best_aligned
                    && ordinal < evidence.best_ordinal)))
}

/// Intern a `PB` tag value into a compact id. Callers must have already
/// rejected records with a missing/empty `PB` -- see `RawKey::sample`'s doc
/// comment: a record with no confident sample-barcode call must not be
/// evidence-eligible at all, not merged into a shared placeholder bucket.
pub fn intern_sample(
    pb: &str,
    samples: &mut Vec<String>,
    sample_ids: &mut HashMap<String, u32>,
) -> u32 {
    if let Some(&id) = sample_ids.get(pb) {
        id
    } else {
        let id = samples.len() as u32;
        samples.push(pb.to_string());
        sample_ids.insert(pb.to_string(), id);
        id
    }
}

/// Stable hash of a record's raw (cell barcode, sample/probe barcode) pair,
/// used only to pick a shard during `collate()` -- deliberately independent
/// of any interned id, so collate never needs gene/sample interning tables
/// built up front. Every real grouping/collision key in this file is scoped
/// to `(cell, sample)`, so routing by this pair alone guarantees a shard
/// never needs to see another shard's data.
pub fn shard_hash(cb: &[u8], pb: Option<&[u8]>) -> u32 {
    let mut value = 0xcbf29ce484222325u64;
    for &byte in cb {
        value ^= byte as u64;
        value = value.wrapping_mul(0x100000001b3);
    }
    // Separator so e.g. cb="AB",pb="C" hashes differently from cb="A",pb="BC".
    value ^= 0xff;
    value = value.wrapping_mul(0x100000001b3);
    for &byte in pb.unwrap_or(&[]) {
        value ^= byte as u64;
        value = value.wrapping_mul(0x100000001b3);
    }
    let mut z = value.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    (z ^ (z >> 31)) as u32
}

pub fn record_identity(record: &bam::Record) -> u64 {
    // A deterministic representative tie-breaker for independently fetched
    // regions. QNAMEs are unique for this input; coordinates/flags make the
    // identity robust if a producer reuses names.
    let mut value = 0xcbf29ce484222325u64;
    for part in [
        record.qname(),
        &record.tid().to_le_bytes(),
        &record.pos().to_le_bytes(),
        &record.flags().to_le_bytes(),
    ] {
        for &byte in part {
            value ^= byte as u64;
            value = value.wrapping_mul(0x100000001b3);
        }
    }
    value
}
