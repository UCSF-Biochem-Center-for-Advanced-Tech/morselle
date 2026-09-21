use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, Mutex,
};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::{CommandFactory, FromArgMatches, Parser};
use flate2::{write::GzEncoder, Compression};
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, Read};
use serde::Serialize;
use tracing::info;

#[derive(Debug, Parser)]
#[command(
    name = "sbx10x-molecules",
    version,
    about = "Correct gene-aware UMIs, resolve gene collisions, and mark molecule representatives"
)]
struct Cli {
    /// Annotated BAM(s) containing CB, UB/UR, GX, and optionally ZG tags.
    /// Multiple files are only supported together with --shards > 1: they're
    /// read directly, in sequence, by collate() -- no pre-merge/sort/index
    /// step needed, since collate() only ever does a plain sequential scan
    /// (never indexed access) regardless of how many files it reads from.
    /// With --shards == 1 (the default), exactly one file is required.
    #[arg(value_name = "BAM", num_args = 1..)]
    bam: Vec<PathBuf>,
    /// Output BAM. Input order is preserved.
    #[arg(short, long)]
    output: PathBuf,
    /// JSON metrics output (default: <output>.metrics.json).
    #[arg(long)]
    metrics: Option<PathBuf>,
    /// BAM compression/decompression threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
    /// Number of indexed BAM-region workers for the evidence pass. One retains
    /// the single-reader implementation; values above one require a .bai/.csi.
    #[arg(long, default_value_t = 1)]
    evidence_workers: usize,
    /// Genomic width of a work item for parallel indexed evidence collection.
    #[arg(long, default_value_t = 10_000_000, hide = true)]
    evidence_chunk_bases: u64,
    /// Write an order-independent sequencing-saturation curve as a gzip TSV.
    /// Rows are emitted for the whole run and every observed cell barcode.
    #[arg(long)]
    saturation: Option<PathBuf>,
    /// Eligible-read interval between saturation observations.
    #[arg(long, default_value_t = 1_000_000)]
    saturation_interval: u64,
    /// Bound peak memory by partitioning evidence into N shards, hashed by
    /// (cell barcode, sample/probe barcode), and processing each shard fully
    /// independently -- collate (bucket records into N temporary BAMs),
    /// then quant (evidence -> correct -> aggregate -> write) one shard at a
    /// time, discarding memory before the next -- rather than materializing
    /// one global evidence/roots/decisions table. This is the same
    /// collate/quant split salmon alevin-fry uses to solve the identical
    /// memory-blowup problem. Real motivation: --evidence-workers 64 on a
    /// full-scale (~11B-read) library OOM-killed at ~1.49TB anon-rss on a
    /// 1.5TB box (see PIPELINE_NOTES.md). Safe by construction: every
    /// grouping/collision key sbx10x-molecules uses is scoped to
    /// (cell, sample), so no shard ever needs to see another shard's data.
    /// Default 1 = today's unsharded behavior, unchanged. Not yet supported
    /// together with --saturation (needs a global read-count known up
    /// front; deferred).
    #[arg(long, default_value_t = 1)]
    shards: usize,
    /// Number of shards to process concurrently (evidence -> correct ->
    /// aggregate -> decide) when --shards > 1. Each concurrent shard holds
    /// its own evidence/roots/molecules/decisions in memory simultaneously,
    /// so peak memory scales roughly linearly with this value -- real
    /// observed peak per shard on CAT06 L1_GEX (8 shards) was ~44-90GB, so
    /// e.g. concurrency=4 projects to ~180-360GB, well within a 1.5TB box's
    /// margin, but always check real per-shard peak memory for your own
    /// data before raising this. Default 1 preserves the original strictly
    /// sequential behavior (lowest memory, slowest). Must not exceed
    /// --shards. Only the final output write is serialized (one shard's
    /// write pass holds the output-file lock at a time); evidence
    /// collection and correction run fully concurrently.
    #[arg(long, default_value_t = 1)]
    shard_concurrency: usize,
    /// Pooled-sample mode: a barcode is a whole sample sequenced far deeper
    /// than a single cell, not a droplet. At that depth, cross-gene UMI
    /// collisions (the same (cell, UMI) independently drawn by two
    /// unrelated genes) become common rather than rare, so the default
    /// single-cell policy of keeping only the higher-support gene and
    /// discarding the rest silently throws away real molecules from
    /// genuinely distinct genes. `--tagseq` counts every gene's corrected
    /// UMI nodes independently instead, tagging shared ones as `Outcome::Shared`
    /// rather than adjudicating a winner. Leaves single-cell runs (the
    /// default) unchanged.
    #[arg(long)]
    tagseq: bool,
    /// 10x Flex mode. Require a non-empty probe/sample barcode (`PB`) on
    /// every counted molecule so identical GEM barcodes from distinct Flex
    /// samples cannot be merged.
    #[arg(long, conflicts_with = "tagseq")]
    flex: bool,
}

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, Ord, PartialOrd)]
struct Dna(u64);

impl Dna {
    fn pack(value: &str) -> Option<Self> {
        let value = value.split('-').next().unwrap_or(value);
        if value.len() > 28 || value.is_empty() {
            return None;
        }
        let mut bits = 0_u64;
        for base in value.bytes() {
            bits = (bits << 2)
                | match base {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => return None,
                };
        }
        Some(Self(bits | ((value.len() as u64) << 56)))
    }

    fn len(self) -> usize {
        (self.0 >> 56) as usize
    }

    fn hamming(self, other: Self) -> Option<u32> {
        if self.len() != other.len() {
            return None;
        }
        let mut xor = (self.0 ^ other.0) & ((1_u64 << (2 * self.len())) - 1);
        let mut distance = 0;
        while xor != 0 {
            distance += ((xor & 3) != 0) as u32;
            xor >>= 2;
        }
        Some(distance)
    }

    /// All sequences at Hamming distance exactly 1 (single-base substitution).
    fn one_mismatch_neighbors(self) -> impl Iterator<Item = Dna> {
        let len = self.len();
        let base = self.0;
        let len_tag = (len as u64) << 56;
        (0..len).flat_map(move |pos| {
            let shift = 2 * (len - 1 - pos);
            let orig = (base >> shift) & 3;
            (0..4u64).filter(move |&b| b != orig).map(move |b| {
                let cleared = base & !(3u64 << shift);
                Dna((cleared | (b << shift)) | len_tag)
            })
        })
    }

    fn unpack(self) -> String {
        let len = self.len();
        let mut out = vec![b'A'; len];
        let mut bits = self.0;
        for base in out.iter_mut().rev() {
            *base = b"ACGT"[(bits & 3) as usize];
            bits >>= 2;
        }
        String::from_utf8(out).unwrap()
    }
}

/// No `PB` tag on this record (non-Flex data, or a probe barcode that
/// wasn't confidently located) -- `sample` is genuinely absent, not merely
/// unknown-but-present, so this must never collide with a real interned id.
const NO_SAMPLE: u32 = u32::MAX;

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct RawKey {
    cell: Dna,
    /// Interned `PB` (probe/sample barcode) tag, or `NO_SAMPLE` if absent.
    /// Required for Flex data: the same 16bp `CB` whitelist is shared
    /// across every sample in a pooled Flex run (confirmed on real CAT06
    /// data -- 99.8% of observed cell barcodes are shared by more than one
    /// sample), so `cell` alone is not a unique cell identity there. Every
    /// place that groups or looks up by `(cell, ...)` for correction or
    /// collision-resolution purposes must include `sample` too, or reads
    /// from different samples sharing a `CB` get silently merged.
    sample: u32,
    gene: u32,
    umi: Dna,
}

#[derive(Clone, Copy, Debug)]
struct Evidence {
    reads: u32,
    best_ordinal: u64,
    best_mapq: u8,
    best_aligned: u32,
    // A 32-bit rank has sub-read resolution even for this 1.6B-read run, and
    // keeps Evidence at its pre-saturation 24-byte layout.
    min_read_hash: u32,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Outcome {
    Unique,
    Won,
    Lost,
    Tie,
    /// --tagseq only: this (cell, UMI) was independently claimed by more
    /// than one gene and every claimant was counted, rather than picking a
    /// single winner.
    Shared,
}

#[derive(Clone, Copy, Debug, PartialEq)]
struct Decision {
    corrected: Dna,
    representative: u64,
    collision_degree: u32,
    outcome: Outcome,
}

#[derive(Default, Serialize)]
struct Metrics {
    input_bam: String,
    output_bam: String,
    tagseq: bool,
    records: u64,
    eligible_records: u64,
    /// Otherwise-eligible (mapped, valid CB+UMI+gene) records dropped for
    /// missing/empty `PB` (no confident sample-barcode call). Never counted
    /// as a molecule: `probe-demux` requires the same non-empty `PB` and
    /// silently discards anything without it, so keeping these in `molecules`
    /// bought nothing downstream -- and worse, letting them share one
    /// `NO_SAMPLE` bucket risked exactly the cross-sample UMI collision
    /// `RawKey`'s `sample` field exists to prevent (see its doc comment: 99.8%
    /// of real CAT06 cell barcodes are shared by more than one sample).
    no_sample_records: u64,
    raw_cell_gene_umis: u64,
    corrected_umi_nodes: u64,
    molecules_before_gene_resolution: u64,
    molecules_counted: u64,
    umi_nodes_collapsed: u64,
    collision_groups: u64,
    collision_molecules_won: u64,
    collision_molecules_lost: u64,
    collision_ties: u64,
    /// --tagseq only: collision-group members counted independently instead
    /// of adjudicated. 0 under the default single-cell policy.
    collision_molecules_shared: u64,
    representatives_tagged: u64,
    saturation: Option<SaturationMetrics>,
    elapsed_secs: f64,
}

/// Adds one shard's local counters into the run-wide accumulator. Used by
/// the concurrent shard path, where each shard builds its own `Metrics`
/// (starting from `Metrics::default()`) rather than sharing one `&mut
/// Metrics` across threads. `input_bam`/`output_bam`/`tagseq`/`saturation`/
/// `elapsed_secs` are set once at the top level, not per-shard, so they're
/// deliberately left alone here.
fn merge_metrics(target: &mut Metrics, source: &Metrics) {
    target.records += source.records;
    target.eligible_records += source.eligible_records;
    target.no_sample_records += source.no_sample_records;
    target.raw_cell_gene_umis += source.raw_cell_gene_umis;
    target.corrected_umi_nodes += source.corrected_umi_nodes;
    target.molecules_before_gene_resolution += source.molecules_before_gene_resolution;
    target.molecules_counted += source.molecules_counted;
    target.umi_nodes_collapsed += source.umi_nodes_collapsed;
    target.collision_groups += source.collision_groups;
    target.collision_molecules_won += source.collision_molecules_won;
    target.collision_molecules_lost += source.collision_molecules_lost;
    target.collision_ties += source.collision_ties;
    target.collision_molecules_shared += source.collision_molecules_shared;
    target.representatives_tagged += source.representatives_tagged;
}

#[derive(Serialize)]
struct SaturationMetrics {
    output: String,
    interval_eligible_reads: u64,
    points: usize,
    observed_cells: usize,
}

/// Stable FNV-1a followed by SplitMix64.  `DefaultHasher` is intentionally not
/// used: saturation curves must be identical across machines and releases.
fn read_hash(name: &[u8]) -> u32 {
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

#[derive(Clone)]
struct CellSaturation {
    barcode: Dna,
    molecule_events: Vec<u64>,
    gene_events: Vec<u64>,
}

struct Saturation {
    total_reads: u64,
    interval: u64,
    cells: HashMap<Dna, usize>,
    per_cell: Vec<CellSaturation>,
    molecule_events: Vec<u64>,
    gene_min_hash: HashMap<u32, u32>,
}

impl Saturation {
    fn new(total_reads: u64, interval: u64) -> Self {
        let points = total_reads.div_ceil(interval) as usize;
        Self {
            total_reads,
            interval,
            cells: HashMap::new(),
            per_cell: Vec::new(),
            molecule_events: vec![0; points],
            gene_min_hash: HashMap::new(),
        }
    }

    fn bin(&self, hash: u32) -> usize {
        // `hash / 2^32 * total_reads` is the zero-based read depth at which a
        // uniformly selected read becomes available.  This avoids any reliance
        // on the coordinate order of the BAM.
        (((hash as u128 * self.total_reads as u128) >> 32) as u64 / self.interval) as usize
    }

    fn cell_index(&mut self, barcode: Dna) -> usize {
        if let Some(&index) = self.cells.get(&barcode) {
            return index;
        }
        let index = self.per_cell.len();
        let points = self.molecule_events.len();
        self.cells.insert(barcode, index);
        self.per_cell.push(CellSaturation {
            barcode,
            molecule_events: vec![0; points],
            gene_events: vec![0; points],
        });
        index
    }

    fn add_cell_gene(&mut self, cell: Dna, gene: u32, hash: u32) {
        let bin = self.bin(hash);
        let index = self.cell_index(cell);
        self.per_cell[index].gene_events[bin] += 1;
        self.gene_min_hash
            .entry(gene)
            .and_modify(|minimum| *minimum = (*minimum).min(hash))
            .or_insert(hash);
    }

    fn add_molecule(&mut self, cell: Dna, hash: u32) {
        let bin = self.bin(hash);
        let index = self.cell_index(cell);
        self.molecule_events[bin] += 1;
        self.per_cell[index].molecule_events[bin] += 1;
    }

    fn write(self, path: &Path) -> Result<SaturationMetrics> {
        let mut global_gene_events = vec![0_u64; self.molecule_events.len()];
        for hash in self.gene_min_hash.values() {
            global_gene_events[self.bin(*hash)] += 1;
        }
        let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
        let mut out = GzEncoder::new(BufWriter::new(file), Compression::new(4));
        writeln!(out, "scope\tcell_barcode\teligible_reads\tread_fraction\tmolecules\tgenes\tnew_molecules\tnew_genes")?;

        let mut global_molecules = 0_u64;
        let mut global_genes = 0_u64;
        for bin in 0..self.molecule_events.len() {
            global_molecules += self.molecule_events[bin];
            global_genes += global_gene_events[bin];
            let depth = ((bin as u64 + 1) * self.interval).min(self.total_reads);
            writeln!(
                out,
                "run\t.\t{depth}\t{}\t{global_molecules}\t{global_genes}\t{}\t{}",
                depth as f64 / self.total_reads as f64,
                self.molecule_events[bin],
                global_gene_events[bin]
            )?;
        }

        let mut cells = self.per_cell;
        cells.sort_by_key(|cell| cell.barcode);
        for cell in cells {
            let mut molecules = 0_u64;
            let mut genes = 0_u64;
            for bin in 0..cell.molecule_events.len() {
                molecules += cell.molecule_events[bin];
                genes += cell.gene_events[bin];
                let depth = ((bin as u64 + 1) * self.interval).min(self.total_reads);
                writeln!(
                    out,
                    "cell\t{}\t{depth}\t{}\t{molecules}\t{genes}\t{}\t{}",
                    cell.barcode.unpack(),
                    depth as f64 / self.total_reads as f64,
                    cell.molecule_events[bin],
                    cell.gene_events[bin]
                )?;
            }
        }
        out.finish()?;
        Ok(SaturationMetrics {
            output: path.display().to_string(),
            interval_eligible_reads: self.interval,
            points: self.molecule_events.len(),
            observed_cells: self.cells.len(),
        })
    }
}

fn aux_string(record: &bam::Record, tag: &[u8; 2]) -> Option<String> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Some(value.to_owned()),
        _ => None,
    }
}

fn choose_gene(record: &bam::Record) -> Option<(String, String)> {
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

fn aligned_bases(record: &bam::Record) -> u32 {
    record
        .cigar()
        .iter()
        .map(|op| match op.char() {
            'M' | '=' | 'X' => op.len(),
            _ => 0,
        })
        .sum()
}

fn better(record: &bam::Record, ordinal: u64, evidence: &Evidence) -> bool {
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
fn intern_sample(
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
fn shard_hash(cb: &[u8], pb: Option<&[u8]>) -> u32 {
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

/// Collate stage (salmon alevin-fry naming): a single sequential pass over
/// the input BAM(s), bucketing every record (eligible or not -- ineligible
/// records just pass through unmodified during quant) into one of
/// `shard_paths.len()` temporary BAM files by `shard_hash` of its `CB`+`PB`
/// tags. A record with no `CB` at all can never be correction-eligible
/// either way, so it's routed to shard 0 as a fixed fallback.
///
/// Reads `inputs` in sequence -- one, or many (e.g. the raw per-chunk
/// aligned BAMs straight from probeset-align/guides, unmerged, unsorted,
/// unindexed). This never needed indexed/sorted access in the first place
/// (only a plain sequential scan), so there was never a real reason to pay
/// for `samtools cat` + `sort` + `index` before this step -- that was
/// always just an artifact of the non-sharded --evidence-workers path's
/// real requirement (IndexedReader), not something collate() itself needs.
/// All input files must share the same header (verified once, from the
/// first file, and reused for every shard writer and every subsequent
/// input file).
fn collate(inputs: &[PathBuf], threads: usize, shard_paths: &[PathBuf]) -> Result<()> {
    let first = inputs
        .first()
        .context("collate() requires at least one input BAM")?;
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
            let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam)
                .with_context(|| format!("creating shard {}", path.display()))?;
            writer.set_compression_level(bam::CompressionLevel::Uncompressed)?;
            Ok(writer)
        })
        .collect::<Result<Vec<_>>>()?;
    let shards = writers.len();
    let mut routed = 0_u64;
    for input in inputs {
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
fn collect_evidence(
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

struct RegionEvidence {
    evidence: HashMap<RawKey, Evidence>,
    genes: Vec<String>,
    samples: Vec<String>,
    records: u64,
    eligible_records: u64,
    no_sample_records: u64,
}

fn merge_item(target: &mut Evidence, source: Evidence) {
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

fn merge_evidence_maps(
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

fn record_identity(record: &bam::Record) -> u64 {
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

fn collect_region(
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

fn indexed_regions(path: &Path, chunk_bases: u64) -> Result<Vec<(u32, u64, u64)>> {
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
fn collect_evidence_indexed(
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

fn correct_group(
    mut keys: Vec<RawKey>,
    evidence: &HashMap<RawKey, Evidence>,
) -> HashMap<RawKey, RawKey> {
    keys.sort_by_key(|key| (std::cmp::Reverse(evidence[key].reads), key.umi));
    // Every raw UMI in a (cell, gene) group is already distinct (RawKey
    // includes the UMI, and `evidence` is keyed on the full RawKey), so
    // indexing by UMI once and generating each child's 3*len single-mismatch
    // neighbors directly replaces an O(n^2) all-pairs Hamming scan with O(n)
    // hash lookups. This matters because raw UMI diversity per (cell, gene)
    // is heavily skewed toward highly expressed genes, where n^2 previously
    // dominated the whole molecules stage.
    let by_umi: HashMap<Dna, RawKey> = keys.iter().map(|&k| (k.umi, k)).collect();
    let mut roots = HashMap::with_capacity(keys.len());
    for &child in &keys {
        let mut parent = child;
        let threshold = evidence[&child].reads.saturating_mul(2).saturating_sub(1);
        for neighbor in child.umi.one_mismatch_neighbors() {
            let Some(&candidate) = by_umi.get(&neighbor) else {
                continue;
            };
            if evidence[&candidate].reads < threshold {
                continue;
            }
            if (evidence[&candidate].reads, std::cmp::Reverse(candidate.umi))
                > (evidence[&parent].reads, std::cmp::Reverse(parent.umi))
            {
                parent = candidate;
            }
        }
        let root = roots.get(&parent).copied().unwrap_or(parent);
        roots.insert(child, root);
    }
    roots
}

fn merge_roots(
    mut left: HashMap<RawKey, RawKey>,
    mut right: HashMap<RawKey, RawKey>,
) -> HashMap<RawKey, RawKey> {
    if left.len() < right.len() {
        std::mem::swap(&mut left, &mut right);
    }
    left.extend(right.drain());
    left
}

fn build_decisions(
    evidence: &HashMap<RawKey, Evidence>,
    metrics: &mut Metrics,
    saturation_interval: Option<u64>,
    workers: usize,
    tagseq: bool,
) -> (HashMap<RawKey, Decision>, Option<Saturation>) {
    // Grouped by (cell, sample, gene), not just (cell, gene): the same CB
    // whitelist is shared across every sample in a pooled Flex run, so
    // omitting `sample` here would let correct_group()'s internal by-UMI
    // map (keyed on UMI alone within a group) silently overwrite one
    // sample's evidence with another's whenever two samples' reads for the
    // same cell+gene happen to share a UMI value.
    let mut groups = HashMap::<(Dna, u32, u32), Vec<RawKey>>::new();
    for &key in evidence.keys() {
        groups
            .entry((key.cell, key.sample, key.gene))
            .or_default()
            .push(key);
    }
    eprintln!(
        "[sbxpress molecules] phase=umi_group_index complete groups={} raw_nodes={}",
        groups.len(),
        evidence.len()
    );
    let mut saturation =
        saturation_interval.map(|interval| Saturation::new(metrics.eligible_records, interval));
    if let Some(curve) = saturation.as_mut() {
        for (&(cell, _sample, gene), keys) in &groups {
            let min_hash = keys
                .iter()
                .map(|key| evidence[key].min_read_hash)
                .min()
                .expect("non-empty UMI group");
            curve.add_cell_gene(cell, gene, min_hash);
        }
    }
    let roots = if tagseq {
        // Pooled-sample mode: directional UMI correction assumes a
        // low-count 1-mismatch neighbor is a sequencing error of a nearby
        // high-count node. At pooled-sample UMI diversity that assumption
        // breaks down the same way cross-gene collision does -- two
        // genuinely distinct molecules can land 1 Hamming-distance apart
        // in a 12bp space purely by chance and get merged into one. Skip
        // correction entirely: every raw (already exact-duplicate-collapsed)
        // node is its own molecule. Trade-off: genuine 1bp sequencing
        // errors in the UMI itself are no longer corrected either, so
        // molecule counts trend slightly high rather than slightly low.
        eprintln!(
            "[sbxpress molecules] phase=umi_correct skipped (--tagseq) raw_nodes={}",
            evidence.len()
        );
        evidence.keys().map(|&key| (key, key)).collect()
    } else if workers > 1 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(workers)
            .build()
            .expect("creating UMI-correction worker pool");
        info!(
            workers,
            groups = groups.len(),
            "correcting UMI groups in parallel"
        );
        eprintln!(
            "[sbxpress molecules] phase=umi_correct start workers={} groups={}",
            workers,
            groups.len()
        );
        let total_groups = groups.len();
        let completed = Arc::new(AtomicUsize::new(0));
        let checkpoint = (total_groups / 100).max(1_000);
        pool.install(|| {
            groups
                .into_par_iter()
                .map(|(_, keys)| {
                    let roots = correct_group(keys, evidence);
                    let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                    if done % checkpoint == 0 || done == total_groups {
                        eprintln!(
                            "[sbxpress molecules] phase=umi_correct groups={done}/{total_groups}"
                        );
                    }
                    roots
                })
                .reduce(HashMap::new, merge_roots)
        })
    } else {
        eprintln!(
            "[sbxpress molecules] phase=umi_correct mode=serial start groups={}",
            groups.len()
        );
        groups
            .into_values()
            .map(|keys| correct_group(keys, evidence))
            .fold(HashMap::new(), merge_roots)
    };
    eprintln!(
        "[sbxpress molecules] phase=umi_correct complete roots={}",
        roots.len()
    );
    // `+=`, not `=`: safe both for the single global call (Metrics starts
    // zeroed, so the result is identical) and for the sharded path, which
    // calls build_decisions() once per shard against the same accumulating
    // Metrics.
    metrics.raw_cell_gene_umis += evidence.len() as u64;
    metrics.umi_nodes_collapsed += roots
        .iter()
        .filter(|(key, root)| key.umi != root.umi)
        .count() as u64;

    #[derive(Clone, Copy)]
    struct Molecule {
        support: u64,
        representative: u64,
        // MAPQ (high 8 bits) followed by aligned bases (low 24 bits).  This
        // preserves the representative ordering without growing this large,
        // temporary hash-map value.
        best_rank: u32,
        min_read_hash: u32,
    }
    eprintln!("[sbxpress molecules] phase=molecule_aggregate start");
    let mut molecules = HashMap::<RawKey, Molecule>::new();
    for (&key, item) in evidence {
        let root = roots[&key];
        molecules
            .entry(root)
            .and_modify(|m| {
                m.support += item.reads as u64;
                m.min_read_hash = m.min_read_hash.min(item.min_read_hash);
                let rank = (u32::from(item.best_mapq) << 24) | item.best_aligned.min(0x00ff_ffff);
                // On an exact rank tie (common for true PCR/optical
                // duplicates of the same molecule, which share mapq and
                // aligned_bases at the same locus), fall back to
                // `best_ordinal` as an explicit, deterministic tiebreak
                // instead of leaving the winner to HashMap iteration order.
                // `best_ordinal` is a positional ordinal in serial mode and
                // a content-based `record_identity()` hash in indexed mode
                // (see `collect_region`) -- either way it's a fixed
                // per-read value, so "smaller wins" gives the same answer
                // regardless of which order this loop visits entries in,
                // which HashMap iteration order does not otherwise
                // guarantee (confirmed by a real cross-run diff: identical
                // molecule/collision decisions, but an inconsistent choice
                // of which duplicate read got the representative flag).
                if rank > m.best_rank
                    || (rank == m.best_rank && item.best_ordinal < m.representative)
                {
                    m.representative = item.best_ordinal;
                    m.best_rank = rank;
                }
            })
            .or_insert(Molecule {
                support: item.reads as u64,
                representative: item.best_ordinal,
                best_rank: (u32::from(item.best_mapq) << 24) | item.best_aligned.min(0x00ff_ffff),
                min_read_hash: item.min_read_hash,
            });
    }
    metrics.corrected_umi_nodes += molecules.len() as u64;
    metrics.molecules_before_gene_resolution += molecules.len() as u64;
    eprintln!(
        "[sbxpress molecules] phase=molecule_aggregate complete corrected_nodes={}",
        molecules.len()
    );

    eprintln!("[sbxpress molecules] phase=gene_collision_index start");
    // (cell, sample, umi): same reasoning as the umi_group_index above --
    // without `sample`, two different samples' molecules sharing a cell+UMI
    // would look like one real molecule colliding across two genes.
    let mut collisions = HashMap::<(Dna, u32, Dna), Vec<RawKey>>::new();
    for &root in molecules.keys() {
        collisions
            .entry((root.cell, root.sample, root.umi))
            .or_default()
            .push(root);
    }
    eprintln!(
        "[sbxpress molecules] phase=gene_collision_index complete collision_keys={}",
        collisions.len()
    );
    eprintln!("[sbxpress molecules] phase=gene_collision_resolve start");
    let mut root_outcomes = HashMap::<RawKey, (Outcome, u32)>::new();
    for roots_for_umi in collisions.values() {
        let degree = roots_for_umi.len() as u32;
        if degree == 1 {
            root_outcomes.insert(roots_for_umi[0], (Outcome::Unique, 1));
            metrics.molecules_counted += 1;
            continue;
        }
        metrics.collision_groups += 1;
        if tagseq {
            // Pooled-sample mode: at this depth a shared UMI is usually two
            // genuinely distinct genes' molecules colliding by chance in a
            // 12bp UMI space, not a mapping artifact -- count every
            // claimant instead of discarding all but the highest-support
            // one.
            for &root in roots_for_umi {
                root_outcomes.insert(root, (Outcome::Shared, degree));
                metrics.collision_molecules_shared += 1;
                metrics.molecules_counted += 1;
            }
            continue;
        }
        let best = roots_for_umi
            .iter()
            .map(|root| molecules[root].support)
            .max()
            .unwrap();
        let winners: Vec<RawKey> = roots_for_umi
            .iter()
            .copied()
            .filter(|root| molecules[root].support == best)
            .collect();
        if winners.len() == 1 {
            for &root in roots_for_umi {
                if root == winners[0] {
                    root_outcomes.insert(root, (Outcome::Won, degree));
                    metrics.collision_molecules_won += 1;
                    metrics.molecules_counted += 1;
                } else {
                    root_outcomes.insert(root, (Outcome::Lost, degree));
                    metrics.collision_molecules_lost += 1;
                }
            }
        } else {
            metrics.collision_ties += 1;
            for &root in roots_for_umi {
                root_outcomes.insert(root, (Outcome::Tie, degree));
            }
        }
    }
    eprintln!(
        "[sbxpress molecules] phase=gene_collision_resolve complete molecules={} collision_groups={}",
        metrics.molecules_counted, metrics.collision_groups
    );
    if let Some(curve) = saturation.as_mut() {
        for (&root, molecule) in &molecules {
            let (outcome, _) = root_outcomes[&root];
            if matches!(outcome, Outcome::Unique | Outcome::Won | Outcome::Shared) {
                curve.add_molecule(root.cell, molecule.min_read_hash);
            }
        }
    }
    eprintln!("[sbxpress molecules] phase=decision_materialize start");
    let decisions: HashMap<RawKey, Decision> = evidence
        .keys()
        .map(|&key| {
            let root = roots[&key];
            let (outcome, collision_degree) = root_outcomes[&root];
            (
                key,
                Decision {
                    corrected: root.umi,
                    representative: molecules[&root].representative,
                    collision_degree,
                    outcome,
                },
            )
        })
        .collect();
    eprintln!(
        "[sbxpress molecules] phase=decision_materialize complete decisions={}",
        decisions.len()
    );
    (decisions, saturation)
}

fn replace_string(record: &mut bam::Record, tag: &[u8; 2], value: &str) -> Result<()> {
    let _ = record.remove_aux(tag);
    record.push_aux(tag, Aux::String(value))?;
    Ok(())
}

fn write_output(
    cli: &Cli,
    decisions: &HashMap<RawKey, Decision>,
    genes: &[String],
    gene_ids: &HashMap<String, u32>,
    sample_ids: &HashMap<String, u32>,
    metrics: &mut Metrics,
) -> Result<()> {
    eprintln!("[sbxpress molecules] phase=bam_write start");
    // Non-sharded path only: cli.shards <= 1 is validated (in main()) to
    // require exactly one input file.
    let mut reader = bam::Reader::from_path(&cli.bam[0])?;
    reader.set_threads(cli.threads)?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)?;
    writer.set_threads(cli.threads)?;
    let mut written_records = 0_u64;
    for (ordinal, result) in reader.records().enumerate() {
        let mut record = result?;
        let key = if !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary() {
            let cell = aux_string(&record, b"CB").and_then(|v| Dna::pack(&v));
            let umi = aux_string(&record, b"UR")
                .or_else(|| aux_string(&record, b"UB"))
                .and_then(|v| Dna::pack(&v));
            let gene = choose_gene(&record).and_then(|(g, _)| gene_ids.get(&g).copied());
            let sample = if cli.flex {
                aux_string(&record, b"PB")
                    .filter(|v| !v.is_empty())
                    .and_then(|pb| sample_ids.get(&pb).copied())
            } else {
                Some(NO_SAMPLE)
            };
            match (cell, sample, gene, umi) {
                (Some(cell), Some(sample), Some(gene), Some(umi)) => Some(RawKey {
                    cell,
                    sample,
                    gene,
                    umi,
                }),
                _ => None,
            }
        } else {
            None
        };
        let _ = record.remove_aux(b"xf");
        for tag in [b"XV", b"XW"] {
            let _ = record.remove_aux(tag);
        }
        if let Some((key, decision)) =
            key.and_then(|key| decisions.get(&key).map(|decision| (key, decision)))
        {
            replace_string(&mut record, b"UB", &decision.corrected.unpack())?;
            // A unique ZG terminal call can resolve a multi-gene GX. Materialize
            // that decision so counted representatives satisfy the single-GX
            // invariant expected by cr-bam-to-mtx.
            if aux_string(&record, b"GX").is_some_and(|gx| gx.contains(';')) {
                let selected = &genes[key.gene as usize];
                let selected_name = choose_gene(&record)
                    .map(|(_, name)| name)
                    .unwrap_or_else(|| selected.clone());
                replace_string(&mut record, b"GX", selected)?;
                replace_string(&mut record, b"GN", &selected_name)?;
            }
            let representative = (if cli.evidence_workers > 1 {
                decision.representative == record_identity(&record)
            } else {
                decision.representative == ordinal as u64
            }) && matches!(
                decision.outcome,
                Outcome::Unique | Outcome::Won | Outcome::Shared
            );
            record.push_aux(b"xf", Aux::I32(if representative { 25 } else { 17 }))?;
            record.push_aux(b"XV", Aux::U32(decision.collision_degree))?;
            record.push_aux(
                b"XW",
                Aux::Char(match decision.outcome {
                    Outcome::Unique => b'U',
                    Outcome::Won => b'W',
                    Outcome::Lost => b'L',
                    Outcome::Tie => b'T',
                    Outcome::Shared => b'S',
                }),
            )?;
            if representative {
                metrics.representatives_tagged += 1;
            }
        }
        writer.write(&record)?;
        written_records += 1;
        if written_records % 25_000_000 == 0 {
            eprintln!(
                "[sbxpress molecules] phase=bam_write records={} representatives={}",
                written_records, metrics.representatives_tagged
            );
        }
    }
    eprintln!(
        "[sbxpress molecules] phase=bam_write complete records={} representatives={}",
        written_records, metrics.representatives_tagged
    );
    Ok(())
}

/// Quant stage's write step for the sharded path: identical logic to
/// `write_output`, except it reads one shard's own temporary BAM (not
/// `cli.bam`) and appends into a shared `writer` opened once for the whole
/// run, rather than creating its own. Always uses ordinal-based
/// representative comparison (never `record_identity`): a shard is always
/// built via the plain serial `collect_evidence()`, so ordinals are stable
/// within this single sequential read of the shard file, matching the
/// ordinals `collect_evidence` assigned when building `decisions` for this
/// same shard.
#[allow(clippy::too_many_arguments)]
fn write_shard_output(
    cli: &Cli,
    shard_path: &Path,
    decisions: &HashMap<RawKey, Decision>,
    genes: &[String],
    gene_ids: &HashMap<String, u32>,
    sample_ids: &HashMap<String, u32>,
    writer: &mut bam::Writer,
    metrics: &mut Metrics,
) -> Result<()> {
    let mut reader = bam::Reader::from_path(shard_path)?;
    reader.set_threads(cli.threads)?;
    for (ordinal, result) in reader.records().enumerate() {
        let mut record = result?;
        let key = if !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary() {
            let cell = aux_string(&record, b"CB").and_then(|v| Dna::pack(&v));
            let umi = aux_string(&record, b"UR")
                .or_else(|| aux_string(&record, b"UB"))
                .and_then(|v| Dna::pack(&v));
            let gene = choose_gene(&record).and_then(|(g, _)| gene_ids.get(&g).copied());
            let sample = if cli.flex {
                aux_string(&record, b"PB")
                    .filter(|v| !v.is_empty())
                    .and_then(|pb| sample_ids.get(&pb).copied())
            } else {
                Some(NO_SAMPLE)
            };
            match (cell, sample, gene, umi) {
                (Some(cell), Some(sample), Some(gene), Some(umi)) => Some(RawKey {
                    cell,
                    sample,
                    gene,
                    umi,
                }),
                _ => None,
            }
        } else {
            None
        };
        let _ = record.remove_aux(b"xf");
        for tag in [b"XV", b"XW"] {
            let _ = record.remove_aux(tag);
        }
        if let Some((key, decision)) =
            key.and_then(|key| decisions.get(&key).map(|decision| (key, decision)))
        {
            replace_string(&mut record, b"UB", &decision.corrected.unpack())?;
            if aux_string(&record, b"GX").is_some_and(|gx| gx.contains(';')) {
                let selected = &genes[key.gene as usize];
                let selected_name = choose_gene(&record)
                    .map(|(_, name)| name)
                    .unwrap_or_else(|| selected.clone());
                replace_string(&mut record, b"GX", selected)?;
                replace_string(&mut record, b"GN", &selected_name)?;
            }
            let representative = decision.representative == ordinal as u64
                && matches!(
                    decision.outcome,
                    Outcome::Unique | Outcome::Won | Outcome::Shared
                );
            record.push_aux(b"xf", Aux::I32(if representative { 25 } else { 17 }))?;
            record.push_aux(b"XV", Aux::U32(decision.collision_degree))?;
            record.push_aux(
                b"XW",
                Aux::Char(match decision.outcome {
                    Outcome::Unique => b'U',
                    Outcome::Won => b'W',
                    Outcome::Lost => b'L',
                    Outcome::Tie => b'T',
                    Outcome::Shared => b'S',
                }),
            )?;
            if representative {
                metrics.representatives_tagged += 1;
            }
        }
        writer.write(&record)?;
    }
    Ok(())
}

/// Full sharded (collate + quant) path: bounds peak memory to roughly
/// `total/shards` by processing one `(cell, sample)`-hashed shard fully
/// (evidence -> correct -> aggregate -> decide -> write) before moving to
/// the next, discarding each shard's `evidence`/`roots`/`molecules`/
/// `collisions`/`decisions` before starting the next one. See the `--shards`
/// CLI doc comment and PIPELINE_NOTES.md for the real OOM this replaces.
fn run_sharded(cli: &Cli, metrics: &mut Metrics) -> Result<()> {
    let shard_dir = PathBuf::from(format!("{}.shards", cli.output.display()));
    fs::create_dir_all(&shard_dir)
        .with_context(|| format!("creating shard directory {}", shard_dir.display()))?;
    let shard_paths: Vec<PathBuf> = (0..cli.shards)
        .map(|i| shard_dir.join(format!("shard_{i:04}.bam")))
        .collect();

    eprintln!(
        "[sbxpress molecules] phase=collate start shards={}",
        cli.shards
    );
    collate(&cli.bam, cli.threads, &shard_paths)?;

    let header = {
        let reader = bam::Reader::from_path(&cli.bam[0])?;
        bam::Header::from_template(reader.header())
    };
    let mut writer = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)?;
    writer.set_threads(cli.threads)?;
    let writer = Mutex::new(writer);

    // Each shard runs evidence-collect -> correct -> aggregate -> decide
    // fully independently (real per-shard proof: sharding_by_cell_sample_
    // matches_global_build_decisions), so shard_concurrency shards can run
    // truly concurrently with no coordination needed until the final write.
    // Only that write is serialized (one shard holds the output-file lock
    // at a time) -- everything before it runs lock-free.
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(cli.shard_concurrency.min(shard_paths.len()))
        .build()
        .context("creating shard-concurrency worker pool")?;
    let total_shards = shard_paths.len();
    let completed = Arc::new(AtomicUsize::new(0));

    let shard_metrics: Vec<Metrics> = pool.install(|| {
        shard_paths
            .par_iter()
            .map(|shard_path| -> Result<Metrics> {
                let mut local_metrics = Metrics::default();
                let (evidence, genes, gene_ids, _samples, sample_ids) =
                    collect_evidence(shard_path, cli.threads, cli.flex, &mut local_metrics)?;
                let (decisions, _saturation) = build_decisions(
                    &evidence,
                    &mut local_metrics,
                    None,
                    cli.evidence_workers,
                    cli.tagseq,
                );
                {
                    let mut writer_guard = writer.lock().expect("output writer mutex poisoned");
                    write_shard_output(
                        cli,
                        shard_path,
                        &decisions,
                        &genes,
                        &gene_ids,
                        &sample_ids,
                        &mut writer_guard,
                        &mut local_metrics,
                    )?;
                }
                let done = completed.fetch_add(1, Ordering::Relaxed) + 1;
                eprintln!(
                    "[sbxpress molecules] phase=shard_quant complete shards={done}/{total_shards}"
                );
                Ok(local_metrics)
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

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();
    let cli = if std::env::var_os("SBXPRESS_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("sbxpress molecules").get_matches())?
    } else {
        Cli::parse()
    };
    if cli.bam.iter().any(|b| b == &cli.output) {
        bail!("input and output BAM must differ");
    }
    let started = Instant::now();
    let mut metrics = Metrics {
        input_bam: cli
            .bam
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join(","),
        output_bam: cli.output.display().to_string(),
        tagseq: cli.tagseq,
        ..Default::default()
    };
    if cli.saturation_interval == 0 {
        bail!("--saturation-interval must be greater than zero");
    }
    if cli.evidence_workers == 0 {
        bail!("--evidence-workers must be greater than zero");
    }
    if cli.shards == 0 {
        bail!("--shards must be greater than zero");
    }
    if cli.shard_concurrency == 0 {
        bail!("--shard-concurrency must be greater than zero");
    }
    if cli.shards <= 1 && cli.bam.len() != 1 {
        bail!(
            "multiple input BAMs are only supported with --shards > 1 (got {} inputs, --shards {}); \
             the non-sharded path needs a single indexed/re-scannable file",
            cli.bam.len(),
            cli.shards
        );
    }
    if cli.shards > 1 {
        if cli.saturation.is_some() {
            bail!("--saturation is not yet supported together with --shards > 1");
        }
        info!(shards = cli.shards, "running sharded (collate + quant)");
        run_sharded(&cli, &mut metrics)?;
        metrics.elapsed_secs = started.elapsed().as_secs_f64();
        let metrics_path = cli
            .metrics
            .clone()
            .unwrap_or_else(|| PathBuf::from(format!("{}.metrics.json", cli.output.display())));
        fs::write(&metrics_path, serde_json::to_string_pretty(&metrics)?)?;
        info!(
            records = metrics.records,
            molecules = metrics.molecules_counted,
            collisions = metrics.collision_groups,
            elapsed = format!("{:.1}s", metrics.elapsed_secs),
            "RUN SUMMARY"
        );
        return Ok(());
    }
    // Validated above: cli.shards <= 1 requires exactly one input file.
    let single_bam = &cli.bam[0];
    let (evidence, genes, gene_ids, _samples, sample_ids) = if cli.evidence_workers > 1 {
        info!(
            workers = cli.evidence_workers,
            "collecting evidence from indexed BAM regions"
        );
        collect_evidence_indexed(
            single_bam,
            cli.evidence_workers,
            cli.evidence_chunk_bases,
            cli.flex,
            &mut metrics,
        )?
    } else {
        collect_evidence(single_bam, cli.threads, cli.flex, &mut metrics)?
    };
    eprintln!(
        "[sbxpress molecules] phase=evidence_ready raw_nodes={} eligible={}",
        evidence.len(),
        metrics.eligible_records
    );
    let (decisions, saturation) = build_decisions(
        &evidence,
        &mut metrics,
        cli.saturation.as_ref().map(|_| cli.saturation_interval),
        cli.evidence_workers,
        cli.tagseq,
    );
    if let (Some(curve), Some(path)) = (saturation, cli.saturation.as_deref()) {
        eprintln!(
            "[sbxpress molecules] phase=saturation_write start output={}",
            path.display()
        );
        metrics.saturation = Some(curve.write(path)?);
        eprintln!("[sbxpress molecules] phase=saturation_write complete");
    }
    write_output(
        &cli,
        &decisions,
        &genes,
        &gene_ids,
        &sample_ids,
        &mut metrics,
    )?;
    metrics.elapsed_secs = started.elapsed().as_secs_f64();
    let metrics_path = cli
        .metrics
        .clone()
        .unwrap_or_else(|| PathBuf::from(format!("{}.metrics.json", cli.output.display())));
    fs::write(&metrics_path, serde_json::to_string_pretty(&metrics)?)?;
    info!(
        records = metrics.records,
        molecules = metrics.molecules_counted,
        collisions = metrics.collision_groups,
        elapsed = format!("{:.1}s", metrics.elapsed_secs),
        "RUN SUMMARY"
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_round_trip_and_hamming() {
        let a = Dna::pack("AACCGGTT").unwrap();
        let b = Dna::pack("AACCGGTA").unwrap();
        assert_eq!(a.unpack(), "AACCGGTT");
        assert_eq!(a.hamming(b), Some(1));
        assert_eq!(Dna::pack("AACCGGTT-1"), Some(a));
        assert_eq!(std::mem::size_of::<Evidence>(), 24);
    }

    #[test]
    fn one_mismatch_neighbors_match_brute_force_hamming() {
        use rand::{Rng, SeedableRng};
        let mut rng = rand_chacha::ChaCha8Rng::seed_from_u64(42);
        let bases = [b'A', b'C', b'G', b'T'];
        for _ in 0..200 {
            let seq: String = (0..12)
                .map(|_| bases[rng.gen_range(0..4)] as char)
                .collect();
            let target = Dna::pack(&seq).unwrap();

            let neighbors: std::collections::HashSet<Dna> =
                target.one_mismatch_neighbors().collect();
            assert_eq!(neighbors.len(), 36, "3*len distinct neighbors for {seq}");
            assert!(!neighbors.contains(&target));
            for &n in &neighbors {
                assert_eq!(target.hamming(n), Some(1));
            }

            // Cross-check against brute-force scan over random candidates:
            // membership in the neighbor set must agree with hamming == 1.
            for _ in 0..50 {
                let candidate_seq: String = (0..12)
                    .map(|_| bases[rng.gen_range(0..4)] as char)
                    .collect();
                let candidate = Dna::pack(&candidate_seq).unwrap();
                assert_eq!(
                    neighbors.contains(&candidate),
                    target.hamming(candidate) == Some(1)
                );
            }
        }
    }

    fn evidence(reads: u32, ordinal: u64) -> Evidence {
        Evidence {
            reads,
            best_ordinal: ordinal,
            best_mapq: 60,
            best_aligned: 100,
            min_read_hash: ordinal as u32,
        }
    }

    #[test]
    fn corrects_within_gene_then_resolves_cross_gene_collision() {
        let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
        let major = Dna::pack("AAAAAAAAAAAA").unwrap();
        let error = Dna::pack("AAAAAAAAAAAC").unwrap();
        let mut input = HashMap::new();
        let gene0_major = RawKey {
            cell,
            sample: NO_SAMPLE,
            gene: 0,
            umi: major,
        };
        let gene0_error = RawKey {
            cell,
            sample: NO_SAMPLE,
            gene: 0,
            umi: error,
        };
        let gene1_major = RawKey {
            cell,
            sample: NO_SAMPLE,
            gene: 1,
            umi: major,
        };
        input.insert(gene0_major, evidence(10, 0));
        input.insert(gene0_error, evidence(2, 1));
        input.insert(gene1_major, evidence(3, 2));

        let mut metrics = Metrics::default();
        let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, false);
        assert_eq!(decisions[&gene0_error].corrected, major);
        assert_eq!(decisions[&gene0_major].outcome, Outcome::Won);
        assert_eq!(decisions[&gene0_error].outcome, Outcome::Won);
        assert_eq!(decisions[&gene1_major].outcome, Outcome::Lost);
        assert_eq!(metrics.umi_nodes_collapsed, 1);
        assert_eq!(metrics.collision_groups, 1);
        assert_eq!(metrics.molecules_counted, 1);

        let mut parallel_metrics = Metrics::default();
        let (parallel, _) = build_decisions(&input, &mut parallel_metrics, None, 2, false);
        assert_eq!(parallel[&gene0_error].corrected, major);
        assert_eq!(parallel[&gene0_major].outcome, Outcome::Won);
        assert_eq!(parallel[&gene1_major].outcome, Outcome::Lost);
        assert_eq!(
            parallel_metrics.molecules_counted,
            metrics.molecules_counted
        );
    }

    /// Real-data-motivated regression: on real CAT06 Flex data, 99.8% of
    /// observed cell barcodes are shared by more than one sample (the same
    /// CB whitelist is reused across every pooled sample), so `cell` alone
    /// is not a unique identity there. Two different samples' reads sharing
    /// a cell+gene+UMI must never be merged into one evidence node, corrected
    /// against each other, or resolved as a cross-gene collision against
    /// each other -- each must come out `Outcome::Unique` on its own.
    #[test]
    fn different_samples_sharing_cell_gene_umi_are_never_merged() {
        let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
        let umi = Dna::pack("AAAAAAAAAAAA").unwrap();
        let mut input = HashMap::new();
        let sample_a = RawKey {
            cell,
            sample: 0,
            gene: 0,
            umi,
        };
        let sample_b = RawKey {
            cell,
            sample: 1,
            gene: 0,
            umi,
        };
        input.insert(sample_a, evidence(10, 0));
        input.insert(sample_b, evidence(7, 1));

        let mut metrics = Metrics::default();
        let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, false);
        // Neither read of `sample_a`'s evidence, nor a collision group of 2.
        assert_eq!(decisions[&sample_a].outcome, Outcome::Unique);
        assert_eq!(decisions[&sample_b].outcome, Outcome::Unique);
        assert_eq!(decisions[&sample_a].collision_degree, 1);
        assert_eq!(decisions[&sample_b].collision_degree, 1);
        assert_eq!(metrics.collision_groups, 0);
        assert_eq!(metrics.molecules_counted, 2);
        assert_eq!(metrics.umi_nodes_collapsed, 0);

        let mut parallel_metrics = Metrics::default();
        let (parallel, _) = build_decisions(&input, &mut parallel_metrics, None, 2, false);
        assert_eq!(parallel[&sample_a].outcome, Outcome::Unique);
        assert_eq!(parallel[&sample_b].outcome, Outcome::Unique);
        assert_eq!(parallel_metrics.molecules_counted, 2);
    }

    #[test]
    fn intern_sample_assigns_stable_ids_per_distinct_pb() {
        let mut samples = Vec::new();
        let mut sample_ids = HashMap::new();
        let a = intern_sample("D-A01", &mut samples, &mut sample_ids);
        let b = intern_sample("D-A02", &mut samples, &mut sample_ids);
        let a_again = intern_sample("D-A01", &mut samples, &mut sample_ids);
        assert_ne!(a, b);
        assert_eq!(a, a_again);
        assert_eq!(samples.len(), 2);
    }

    /// A record with valid CB+UMI+GX but a missing or empty PB must be
    /// dropped entirely, not folded into a shared placeholder bucket --
    /// see RawKey::sample's doc comment (99.8% of real CAT06 cell barcodes
    /// are shared by more than one sample, so a bucket that doesn't
    /// discriminate by sample risks merging two different samples' real
    /// molecules). probe-demux independently requires the same non-empty
    /// PB and discards anything without it, so such reads could never
    /// have become a real counted molecule downstream anyway.
    #[test]
    fn records_without_a_confident_sample_barcode_are_dropped_not_bucketed() {
        use rust_htslib::bam::{header::HeaderRecord, record::Cigar, record::CigarString};

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("input.bam");
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "ENSG1")
                .push_tag(b"LN", 1_000_000),
        );
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        let cigar = CigarString(vec![Cigar::Match(20)]);
        // (name, cb, umi, pb) -- pb: None = no tag at all, Some("") = tag
        // present but empty (a located-but-unresolved probe barcode).
        for (name, cb, umi, pb) in [
            ("no_pb_tag", "AAAACCCCGGGGTTTT", "AAAAAAAAAAAA", None),
            ("empty_pb", "AAAACCCCGGGGTTTT", "CCCCCCCCCCCC", Some("")),
            ("real_pb", "AAAACCCCGGGGTTTT", "GGGGGGGGGGGG", Some("D-A01")),
        ] {
            let mut record = bam::Record::new();
            record.set(
                name.as_bytes(),
                Some(&cigar),
                b"ACGTACGTACGTACGTACGT",
                &[30; 20],
            );
            record.set_tid(0);
            record.set_pos(1);
            record.set_flags(0);
            record.push_aux(b"CB", Aux::String(cb)).unwrap();
            record.push_aux(b"UR", Aux::String(umi)).unwrap();
            record.push_aux(b"GX", Aux::String("ENSG1")).unwrap();
            if let Some(pb) = pb {
                record.push_aux(b"PB", Aux::String(pb)).unwrap();
            }
            writer.write(&record).unwrap();
        }
        drop(writer);

        let mut metrics = Metrics::default();
        let (evidence, ..) = collect_evidence(&path, 1, true, &mut metrics).unwrap();
        assert_eq!(
            metrics.eligible_records, 1,
            "only the real-PB record is eligible"
        );
        assert_eq!(
            metrics.no_sample_records, 2,
            "both no-PB and empty-PB records are counted as dropped"
        );
        assert_eq!(evidence.len(), 1);

        let mut tagseq_metrics = Metrics::default();
        let (tagseq_evidence, ..) = collect_evidence(&path, 1, false, &mut tagseq_metrics).unwrap();
        assert_eq!(
            tagseq_metrics.eligible_records, 3,
            "tag-seq must not require PB"
        );
        assert_eq!(tagseq_metrics.no_sample_records, 0);
        assert_eq!(tagseq_evidence.len(), 3);
    }

    #[test]
    fn tagseq_counts_every_gene_in_a_collision_instead_of_picking_a_winner() {
        // Same raw evidence as corrects_within_gene_then_resolves_cross_gene_collision,
        // but under --tagseq: UMI correction is also skipped (pooled-sample UMI
        // diversity makes the 1-mismatch "sequencing error" assumption unsafe
        // the same way winner-take-all collision resolution is unsafe), so
        // gene0's "error" UMI no longer merges into "major" -- it stands as
        // its own uncollided molecule. gene0_major and gene1_major still
        // collide on the shared "major" UMI, and --tagseq keeps both instead
        // of picking a winner.
        let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
        let major = Dna::pack("AAAAAAAAAAAA").unwrap();
        let error = Dna::pack("AAAAAAAAAAAC").unwrap();
        let mut input = HashMap::new();
        let gene0_major = RawKey {
            cell,
            sample: NO_SAMPLE,
            gene: 0,
            umi: major,
        };
        let gene0_error = RawKey {
            cell,
            sample: NO_SAMPLE,
            gene: 0,
            umi: error,
        };
        let gene1_major = RawKey {
            cell,
            sample: NO_SAMPLE,
            gene: 1,
            umi: major,
        };
        input.insert(gene0_major, evidence(10, 0));
        input.insert(gene0_error, evidence(2, 1));
        input.insert(gene1_major, evidence(3, 2));

        let mut metrics = Metrics::default();
        let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, true);
        // No UMI correction under --tagseq: "error" stays uncorrected and
        // uncollided (no other gene/UMI shares it).
        assert_eq!(decisions[&gene0_error].corrected, error);
        assert_eq!(decisions[&gene0_error].outcome, Outcome::Unique);
        // "major" is shared by two genes -- both kept, neither discarded.
        assert_eq!(decisions[&gene0_major].outcome, Outcome::Shared);
        assert_eq!(decisions[&gene1_major].outcome, Outcome::Shared);
        assert_eq!(decisions[&gene0_major].collision_degree, 2);
        assert_eq!(decisions[&gene1_major].collision_degree, 2);
        assert_eq!(metrics.umi_nodes_collapsed, 0);
        // Diagnostics still reflect that a collision happened...
        assert_eq!(metrics.collision_groups, 1);
        // ...but nothing is discarded: all three raw nodes are counted.
        assert_eq!(metrics.collision_molecules_shared, 2);
        assert_eq!(metrics.collision_molecules_won, 0);
        assert_eq!(metrics.collision_molecules_lost, 0);
        assert_eq!(metrics.molecules_counted, 3);

        let mut parallel_metrics = Metrics::default();
        let (parallel, _) = build_decisions(&input, &mut parallel_metrics, None, 2, true);
        assert_eq!(parallel[&gene0_error].corrected, error);
        assert_eq!(parallel[&gene0_major].outcome, Outcome::Shared);
        assert_eq!(parallel[&gene1_major].outcome, Outcome::Shared);
        assert_eq!(
            parallel_metrics.molecules_counted,
            metrics.molecules_counted
        );
    }

    #[test]
    fn equal_cross_gene_support_is_not_counted() {
        let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
        let umi = Dna::pack("ACGTACGTACGT").unwrap();
        let a = RawKey {
            cell,
            sample: NO_SAMPLE,
            gene: 0,
            umi,
        };
        let b = RawKey {
            cell,
            sample: NO_SAMPLE,
            gene: 1,
            umi,
        };
        let input = HashMap::from([(a, evidence(2, 0)), (b, evidence(2, 1))]);
        let mut metrics = Metrics::default();
        let (decisions, _) = build_decisions(&input, &mut metrics, None, 1, false);
        assert_eq!(decisions[&a].outcome, Outcome::Tie);
        assert_eq!(decisions[&b].outcome, Outcome::Tie);
        assert_eq!(metrics.molecules_counted, 0);
        assert_eq!(metrics.collision_ties, 1);
    }

    #[test]
    fn saturation_is_nested_and_reports_run_and_cells() {
        let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
        let mut curve = Saturation::new(10, 5);
        // Hash 0 arrives in the first five-read sample; this molecule and its
        // gene must remain counted at all later depths.
        curve.add_cell_gene(cell, 7, 0);
        curve.add_molecule(cell, 0);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("curve.tsv.gz");
        let summary = curve.write(&path).unwrap();
        assert_eq!(summary.points, 2);
        assert_eq!(summary.observed_cells, 1);
        let contents = {
            use flate2::read::MultiGzDecoder;
            use std::io::Read;
            let mut decoder = MultiGzDecoder::new(File::open(path).unwrap());
            let mut text = String::new();
            decoder.read_to_string(&mut text).unwrap();
            text
        };
        assert!(contents.contains("run\t.\t5\t0.5\t1\t1\t1\t1"));
        assert!(contents.contains("run\t.\t10\t1\t1\t1\t0\t0"));
        assert_eq!(contents.lines().count(), 5); // header, two run, two cell rows
    }

    #[test]
    fn indexed_regions_collect_each_alignment_once() {
        use rust_htslib::bam::{header::HeaderRecord, record::Cigar, record::CigarString};

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("input.bam");
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 100),
        );
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        let cigar = CigarString(vec![Cigar::Match(20)]);
        for (name, pos, umi) in [
            ("a", 1, "AAAAAAAAAAAA"),
            ("b", 9, "AAAAAAAAAAAA"),
            ("c", 10, "CCCCCCCCCCCC"),
        ] {
            let mut record = bam::Record::new();
            record.set(
                name.as_bytes(),
                Some(&cigar),
                b"ACGTACGTACGTACGTACGT",
                &[30; 20],
            );
            record.set_tid(0);
            record.set_pos(pos);
            record.set_flags(0);
            record
                .push_aux(b"CB", Aux::String("AACCGGTTAACCGGTT"))
                .unwrap();
            record.push_aux(b"UR", Aux::String(umi)).unwrap();
            record.push_aux(b"GX", Aux::String("ENSG1")).unwrap();
            record.push_aux(b"PB", Aux::String("D-A01")).unwrap();
            writer.write(&record).unwrap();
        }
        drop(writer);
        bam::index::build(&path, None, bam::index::Type::Bai, 1).unwrap();
        assert_eq!(indexed_regions(&path, 10).unwrap().len(), 10);
        assert_eq!(
            collect_region(&path, 0, 0, 10, false)
                .unwrap()
                .eligible_records,
            2
        );

        let mut serial_metrics = Metrics::default();
        let (serial, _, _, _, _) = collect_evidence(&path, 1, false, &mut serial_metrics).unwrap();
        let mut parallel_metrics = Metrics::default();
        let (parallel, _, _, _, _) =
            collect_evidence_indexed(&path, 2, 10, false, &mut parallel_metrics).unwrap();
        assert_eq!(serial_metrics.eligible_records, 3);
        assert_eq!(parallel_metrics.eligible_records, 3);
        assert_eq!(serial.len(), 2);
        assert_eq!(parallel.len(), 2);
        assert_eq!(
            serial.values().map(|item| item.reads as u64).sum::<u64>(),
            3
        );
        assert_eq!(
            parallel.values().map(|item| item.reads as u64).sum::<u64>(),
            3
        );
        let mut serial_decision_metrics = Metrics::default();
        let mut parallel_decision_metrics = Metrics::default();
        let _ = build_decisions(&serial, &mut serial_decision_metrics, None, 2, false);
        let _ = build_decisions(&parallel, &mut parallel_decision_metrics, None, 2, false);
        assert_eq!(
            serial_decision_metrics.molecules_counted,
            parallel_decision_metrics.molecules_counted
        );
    }

    /// Core correctness guarantee behind `--shards`: every real grouping
    /// (`(cell, sample, gene)`) and collision-resolution (`(cell, sample,
    /// umi)`) key in `build_decisions` is scoped to `(cell, sample)`, so
    /// partitioning evidence by `(cell, sample)` into independent shards and
    /// merging each shard's own `build_decisions()` output must produce
    /// exactly the decisions (and exactly the accumulated metrics) that a
    /// single global `build_decisions()` call over the same evidence would
    /// -- this is what makes it safe to process shards one at a time,
    /// discarding memory in between, instead of materializing one global
    /// evidence/roots/decisions table.
    #[test]
    fn sharding_by_cell_sample_matches_global_build_decisions() {
        let cells: Vec<Dna> = [
            "AAAAAAAAAAAAAAAA",
            "CCCCCCCCCCCCCCCC",
            "GGGGGGGGGGGGGGGG",
            "ACGTACGTACGTACGT",
            "TGCATGCATGCATGCA",
            "AACCGGTTAACCGGTT",
        ]
        .iter()
        .map(|s| Dna::pack(s).unwrap())
        .collect();
        let umis: Vec<Dna> = [
            "AAAAAAAAAAAA",
            "AAAAAAAAAAAC",
            "GGGGGGGGGGGG",
            "TTTTTTTTTTTT",
        ]
        .iter()
        .map(|s| Dna::pack(s).unwrap())
        .collect();

        // A synthetic but real-shaped evidence set: several cells, some with
        // more than one sample (matching the pooled-Flex architecture this
        // whole fix targets), several genes, a within-group UMI-correction
        // pair (umis[0]/umis[1] differ by 1 base), and a cross-gene
        // collision (umis[2] shared by two genes in the same cell+sample).
        let mut evidence = HashMap::new();
        let mut next_ordinal = 0_u64;
        let mut push = |evidence: &mut HashMap<RawKey, Evidence>,
                        cell: Dna,
                        sample: u32,
                        gene: u32,
                        umi: Dna,
                        reads: u32| {
            let ordinal = next_ordinal;
            next_ordinal += 1;
            evidence.insert(
                RawKey {
                    cell,
                    sample,
                    gene,
                    umi,
                },
                Evidence {
                    reads,
                    best_ordinal: ordinal,
                    best_mapq: 60,
                    // Distinct per entry (real reads essentially never tie
                    // exactly): an exact rank tie makes `molecule_aggregate`'s
                    // representative choice depend on HashMap iteration
                    // order, which a shard's own map instance need not share
                    // with the global one -- a real but orthogonal,
                    // pre-existing property of tie-breaking, not something
                    // this test is meant to exercise.
                    best_aligned: 100 + ordinal as u32,
                    min_read_hash: ordinal as u32,
                },
            );
        };
        for (i, &cell) in cells.iter().enumerate() {
            let sample = (i % 3) as u32; // several cells reuse the same sample id, and vice versa
            push(&mut evidence, cell, sample, 0, umis[0], 10); // major
            push(&mut evidence, cell, sample, 0, umis[1], 2); // 1bp error of major -> corrects into it
            push(&mut evidence, cell, sample, 1, umis[2], 4); // gene 1, umi X
            push(&mut evidence, cell, sample, 2, umis[2], 3); // gene 2, same umi X, same (cell,sample) -> real cross-gene collision
                                                              // Same cell, a DIFFERENT sample, reusing the same UMI+a new gene:
                                                              // must NOT be pulled into the collision group above.
            push(&mut evidence, cell, (sample + 1) % 3, 3, umis[2], 5);
            push(&mut evidence, cell, sample, 4, umis[3], 1); // plain unique molecule
        }

        let mut global_metrics = Metrics::default();
        let (global_decisions, _) = build_decisions(&evidence, &mut global_metrics, None, 1, false);

        const SHARDS: u64 = 4;
        let mut shard_maps: Vec<HashMap<RawKey, Evidence>> =
            (0..SHARDS).map(|_| HashMap::new()).collect();
        for (&key, &item) in &evidence {
            // Mirrors shard_hash's real invariant (route by (cell, sample)
            // alone) without needing raw bytes -- any valid partition
            // function proves the same property.
            let idx = (key.cell.0 ^ (key.sample as u64).wrapping_mul(0x9E3779B97F4A7C15)) % SHARDS;
            shard_maps[idx as usize].insert(key, item);
        }

        let mut sharded_decisions = HashMap::new();
        let mut sharded_metrics = Metrics::default();
        for shard in &shard_maps {
            if shard.is_empty() {
                continue;
            }
            let (decisions, _) = build_decisions(shard, &mut sharded_metrics, None, 1, false);
            sharded_decisions.extend(decisions);
        }

        assert_eq!(sharded_decisions.len(), global_decisions.len());
        for (key, global_decision) in &global_decisions {
            assert_eq!(
                sharded_decisions.get(key),
                Some(global_decision),
                "decision for {key:?} differs between sharded and global runs"
            );
        }
        assert_eq!(
            sharded_metrics.molecules_counted,
            global_metrics.molecules_counted
        );
        assert_eq!(
            sharded_metrics.collision_groups,
            global_metrics.collision_groups
        );
        assert_eq!(
            sharded_metrics.umi_nodes_collapsed,
            global_metrics.umi_nodes_collapsed
        );
        assert_eq!(
            sharded_metrics.raw_cell_gene_umis,
            global_metrics.raw_cell_gene_umis
        );
        // Sanity: this synthetic input actually exercises correction and
        // collision resolution, not just a trivial all-unique case.
        assert!(global_metrics.umi_nodes_collapsed > 0);
        assert!(global_metrics.collision_groups > 0);
    }

    #[test]
    fn shard_hash_is_stable_and_distinguishes_cb_pb_boundary() {
        let a = shard_hash(b"AB", Some(b"C"));
        let b = shard_hash(b"A", Some(b"BC"));
        assert_ne!(
            a, b,
            "CB/PB boundary must affect the hash, not just concatenation"
        );
        assert_eq!(
            a,
            shard_hash(b"AB", Some(b"C")),
            "hash must be deterministic"
        );
        assert_eq!(
            shard_hash(b"AB", None),
            shard_hash(b"AB", Some(b"")),
            "no PB tag and an empty PB tag hash the same (both contribute zero bytes)"
        );
    }

    fn write_test_bam(path: &Path, records: &[(&str, i64, &str, &str, &str)]) {
        use rust_htslib::bam::{header::HeaderRecord, record::Cigar, record::CigarString};
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "ENSG1")
                .push_tag(b"LN", 1_000_000),
        );
        let mut writer = bam::Writer::from_path(path, &header, bam::Format::Bam).unwrap();
        let cigar = CigarString(vec![Cigar::Match(20)]);
        for (name, pos, cb, umi, gene) in records {
            let mut record = bam::Record::new();
            record.set(
                name.as_bytes(),
                Some(&cigar),
                b"ACGTACGTACGTACGTACGT",
                &[30; 20],
            );
            record.set_tid(0);
            record.set_pos(*pos);
            record.set_flags(0);
            record.push_aux(b"CB", Aux::String(cb)).unwrap();
            record.push_aux(b"UR", Aux::String(umi)).unwrap();
            record.push_aux(b"GX", Aux::String(gene)).unwrap();
            writer.write(&record).unwrap();
        }
    }

    fn read_tag_set(path: &Path) -> std::collections::BTreeSet<(String, String, String)> {
        let mut reader = bam::Reader::from_path(path).unwrap();
        reader
            .records()
            .map(|r| {
                let record = r.unwrap();
                (
                    String::from_utf8(record.qname().to_vec()).unwrap(),
                    aux_string(&record, b"CB").unwrap(),
                    aux_string(&record, b"UR").unwrap(),
                )
            })
            .collect()
    }

    /// Core correctness guarantee behind reading multiple un-merged input
    /// BAMs directly (the real motivation: skip samtools cat/sort/index
    /// entirely, since collate() only ever does a plain sequential scan --
    /// see the `bam` CLI field doc comment). Splitting the exact same
    /// records across two "chunk" files must route every record into the
    /// same shard, in the same multiset, as collating from one file holding
    /// all of them -- collate() doesn't care how many files or what order
    /// the records arrive in, only their CB+PB content.
    #[test]
    fn collate_from_multiple_inputs_matches_collate_from_one() {
        let dir = tempfile::tempdir().unwrap();
        let all_records = [
            ("r1", 1, "AAAACCCCGGGGTTTTAAAA", "AAAAAAAAAAAA", "ENSG1"),
            ("r2", 2, "AAAACCCCGGGGTTTTAAAA", "CCCCCCCCCCCC", "ENSG1"),
            ("r3", 3, "CCCCGGGGTTTTAAAACCCC", "GGGGGGGGGGGG", "ENSG1"),
            ("r4", 4, "GGGGTTTTAAAACCCCGGGG", "TTTTTTTTTTTT", "ENSG1"),
            ("r5", 5, "TTTTAAAACCCCGGGGTTTT", "AAAACCCCGGGG", "ENSG1"),
        ];

        // One pre-merged file, all records together.
        let merged_path = dir.path().join("merged.bam");
        write_test_bam(&merged_path, &all_records);
        let merged_shards: Vec<PathBuf> = (0..3)
            .map(|i| dir.path().join(format!("merged_shard_{i}.bam")))
            .collect();
        collate(&[merged_path], 1, &merged_shards).unwrap();

        // The same records, split across three un-merged "chunk" files
        // (simulating raw per-chunk probeset-align output, no cat/sort/index).
        let chunk_paths: Vec<PathBuf> = (0..3)
            .map(|i| dir.path().join(format!("chunk_{i}.bam")))
            .collect();
        write_test_bam(&chunk_paths[0], &all_records[0..2]);
        write_test_bam(&chunk_paths[1], &all_records[2..3]);
        write_test_bam(&chunk_paths[2], &all_records[3..5]);
        let split_shards: Vec<PathBuf> = (0..3)
            .map(|i| dir.path().join(format!("split_shard_{i}.bam")))
            .collect();
        collate(&chunk_paths, 1, &split_shards).unwrap();

        for i in 0..3 {
            assert_eq!(
                read_tag_set(&merged_shards[i]),
                read_tag_set(&split_shards[i]),
                "shard {i} differs between single-file and multi-file collate"
            );
        }
    }
}
