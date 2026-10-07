use super::types::Evidence;
use crate::bam_helpers::for_each_aux;
use rust_htslib::bam::{self, record::Aux};
use std::collections::HashMap;

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
    let tags = MoleculeTags::read(record);
    tags.gene()
        .map(|(id, name)| (id.to_owned(), name.to_owned()))
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
pub fn intern_sample<S: std::hash::BuildHasher>(
    pb: &str,
    samples: &mut Vec<String>,
    sample_ids: &mut HashMap<String, u32, S>,
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

/// Borrowed views of every tag `molecules` reads, gathered in ONE pass over
/// the record's aux data. Real annotated records carry ~24 aux tags
/// (minimap2's, extract's comment tags, annotate's), and `record.aux(tag)`
/// is a linear scan from the start each call, so the 6-7 separate
/// `aux_string` lookups per record (each also allocating a `String`) were a
/// large share of the per-record cost in every BAM pass. Semantics match
/// `aux_string`: only `Z` (string) values count, and the first occurrence of
/// a tag wins.
#[derive(Default, Clone, Copy)]
pub struct MoleculeTags<'a> {
    pub cb: Option<&'a str>,
    pub ur: Option<&'a str>,
    pub ub: Option<&'a str>,
    pub gx: Option<&'a str>,
    pub gn: Option<&'a str>,
    pub zg: Option<&'a str>,
    pub pb: Option<&'a str>,
    /// Any tag `molecules` itself writes (`xf`/`XV`/`XW`/`ZW`) is already
    /// present, e.g. when re-running on a molecules BAM. Lets the writer skip
    /// four `remove_aux` calls -- each a full aux scan -- in the common case.
    pub has_output_tags: bool,
}

impl<'a> MoleculeTags<'a> {
    /// One bounds-checked walk over the raw aux bytes. Only the seven wanted
    /// tags are decoded (and UTF-8-checked); every other field is skipped by
    /// its type length, the way htslib's `bam_aux_get` walks. For each wanted
    /// tag the FIRST occurrence decides, matching `record.aux(tag)`: a first
    /// occurrence that is not a valid `Z`/`H` string yields `None` (later
    /// duplicates are not consulted). Malformed aux data ends the walk.
    pub fn read(record: &'a bam::Record) -> Self {
        let mut tags = Self::default();
        let mut seen = 0u8;
        for_each_aux(record, |tag, type_code, value| {
            let (index, slot) = match &tag {
                b"CB" => (0, &mut tags.cb),
                b"UR" => (1, &mut tags.ur),
                b"UB" => (2, &mut tags.ub),
                b"GX" => (3, &mut tags.gx),
                b"GN" => (4, &mut tags.gn),
                b"ZG" => (5, &mut tags.zg),
                b"PB" => (6, &mut tags.pb),
                b"xf" | b"XV" | b"XW" | b"ZW" => {
                    tags.has_output_tags = true;
                    return;
                }
                _ => return,
            };
            if seen & (1 << index) == 0 {
                seen |= 1 << index;
                // rust-htslib decodes both `Z` and `H` as `Aux::String`;
                // `value` includes the NUL terminator.
                if type_code == b'Z' || type_code == b'H' {
                    *slot = std::str::from_utf8(&value[..value.len() - 1]).ok();
                }
            }
        });
        tags
    }

    /// `UR`, falling back to `UB` (raw UMI preferred over corrected).
    pub fn umi(&self) -> Option<&'a str> {
        self.ur.or(self.ub)
    }

    /// Borrowing equivalent of [`choose_gene`].
    pub fn gene(&self) -> Option<(&'a str, &'a str)> {
        choose_gene_parts(self.gx?, self.gn, self.zg)
    }
}

/// Allocation-free core of [`choose_gene`]: `(gene id, gene name)` for a
/// single-gene `GX`, or the `ZG`-selected gene of a multi-gene `GX`.
pub fn choose_gene_parts<'a>(
    gx: &'a str,
    gn: Option<&'a str>,
    zg: Option<&'a str>,
) -> Option<(&'a str, &'a str)> {
    let gn = gn.unwrap_or(gx);
    let mut genes = gx.split(';').filter(|v| !v.is_empty());
    let first = genes.next();
    if let (Some(only), None) = (first, genes.clone().next()) {
        // `split` always yields at least one item, so the first name exists.
        return Some((only, gn.split(';').next().unwrap_or(only)));
    }
    let zg = zg?;
    if zg.contains(';') || zg.is_empty() {
        return None;
    }
    let index = gx
        .split(';')
        .filter(|v| !v.is_empty())
        .position(|gene| gene == zg)?;
    let name = gn.split(';').nth(index).unwrap_or(zg);
    Some((zg, name))
}
