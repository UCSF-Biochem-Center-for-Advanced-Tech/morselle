use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::{CommandFactory, FromArgMatches, Parser};
use flate2::{write::GzEncoder, Compression};
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
    /// Annotated BAM containing CB, UB/UR, GX, and optionally ZG tags.
    #[arg(value_name = "BAM")]
    bam: PathBuf,
    /// Output BAM. Input order is preserved.
    #[arg(short, long)]
    output: PathBuf,
    /// JSON metrics output (default: <output>.metrics.json).
    #[arg(long)]
    metrics: Option<PathBuf>,
    /// BAM compression/decompression threads.
    #[arg(long, default_value_t = 4)]
    threads: usize,
    /// Write an order-independent sequencing-saturation curve as a gzip TSV.
    /// Rows are emitted for the whole run and every observed cell barcode.
    #[arg(long)]
    saturation: Option<PathBuf>,
    /// Eligible-read interval between saturation observations.
    #[arg(long, default_value_t = 1_000_000)]
    saturation_interval: u64,
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

#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
struct RawKey {
    cell: Dna,
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
}

#[derive(Clone, Copy, Debug)]
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
    records: u64,
    eligible_records: u64,
    raw_cell_gene_umis: u64,
    corrected_umi_nodes: u64,
    molecules_before_gene_resolution: u64,
    molecules_counted: u64,
    umi_nodes_collapsed: u64,
    collision_groups: u64,
    collision_molecules_won: u64,
    collision_molecules_lost: u64,
    collision_ties: u64,
    representatives_tagged: u64,
    saturation: Option<SaturationMetrics>,
    elapsed_secs: f64,
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

fn collect_evidence(
    path: &Path,
    threads: usize,
    metrics: &mut Metrics,
) -> Result<(HashMap<RawKey, Evidence>, Vec<String>, HashMap<String, u32>)> {
    let mut reader =
        bam::Reader::from_path(path).with_context(|| format!("opening {}", path.display()))?;
    reader.set_threads(threads)?;
    let mut evidence = HashMap::<RawKey, Evidence>::new();
    let mut genes = Vec::<String>::new();
    let mut gene_ids = HashMap::<String, u32>::new();
    for (ordinal, result) in reader.records().enumerate() {
        let record = result.context("reading BAM in evidence pass")?;
        metrics.records += 1;
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
        let gene = if let Some(&id) = gene_ids.get(&gene) {
            id
        } else {
            let id = genes.len() as u32;
            genes.push(gene.clone());
            gene_ids.insert(gene, id);
            id
        };
        metrics.eligible_records += 1;
        let key = RawKey { cell, gene, umi };
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
    Ok((evidence, genes, gene_ids))
}

fn build_decisions(
    evidence: &HashMap<RawKey, Evidence>,
    metrics: &mut Metrics,
    saturation_interval: Option<u64>,
) -> (HashMap<RawKey, Decision>, Option<Saturation>) {
    let mut groups = HashMap::<(Dna, u32), Vec<RawKey>>::new();
    for &key in evidence.keys() {
        groups.entry((key.cell, key.gene)).or_default().push(key);
    }
    let mut saturation =
        saturation_interval.map(|interval| Saturation::new(metrics.eligible_records, interval));
    if let Some(curve) = saturation.as_mut() {
        for (&(cell, gene), keys) in &groups {
            let min_hash = keys
                .iter()
                .map(|key| evidence[key].min_read_hash)
                .min()
                .expect("non-empty UMI group");
            curve.add_cell_gene(cell, gene, min_hash);
        }
    }
    let mut roots = HashMap::<RawKey, RawKey>::new();
    for keys in groups.values_mut() {
        keys.sort_by_key(|key| (std::cmp::Reverse(evidence[key].reads), key.umi));
        for &child in keys.iter() {
            let mut parent = child;
            for &candidate in keys.iter() {
                if candidate == child
                    || evidence[&candidate].reads
                        < evidence[&child].reads.saturating_mul(2).saturating_sub(1)
                {
                    continue;
                }
                if child.umi.hamming(candidate.umi) == Some(1)
                    && (evidence[&candidate].reads, std::cmp::Reverse(candidate.umi))
                        > (evidence[&parent].reads, std::cmp::Reverse(parent.umi))
                {
                    parent = candidate;
                }
            }
            let root = roots.get(&parent).copied().unwrap_or(parent);
            roots.insert(child, root);
        }
    }
    metrics.raw_cell_gene_umis = evidence.len() as u64;
    metrics.umi_nodes_collapsed = roots
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
    let mut molecules = HashMap::<RawKey, Molecule>::new();
    for (&key, item) in evidence {
        let root = roots[&key];
        molecules
            .entry(root)
            .and_modify(|m| {
                m.support += item.reads as u64;
                m.min_read_hash = m.min_read_hash.min(item.min_read_hash);
                let rank = (u32::from(item.best_mapq) << 24) | item.best_aligned.min(0x00ff_ffff);
                if rank > m.best_rank {
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
    metrics.corrected_umi_nodes = molecules.len() as u64;
    metrics.molecules_before_gene_resolution = molecules.len() as u64;

    let mut collisions = HashMap::<(Dna, Dna), Vec<RawKey>>::new();
    for &root in molecules.keys() {
        collisions
            .entry((root.cell, root.umi))
            .or_default()
            .push(root);
    }
    let mut root_outcomes = HashMap::<RawKey, (Outcome, u32)>::new();
    for roots_for_umi in collisions.values() {
        let degree = roots_for_umi.len() as u32;
        if degree == 1 {
            root_outcomes.insert(roots_for_umi[0], (Outcome::Unique, 1));
            metrics.molecules_counted += 1;
            continue;
        }
        metrics.collision_groups += 1;
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
    if let Some(curve) = saturation.as_mut() {
        for (&root, molecule) in &molecules {
            let (outcome, _) = root_outcomes[&root];
            if matches!(outcome, Outcome::Unique | Outcome::Won) {
                curve.add_molecule(root.cell, molecule.min_read_hash);
            }
        }
    }
    let decisions = evidence
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
    metrics: &mut Metrics,
) -> Result<()> {
    let mut reader = bam::Reader::from_path(&cli.bam)?;
    reader.set_threads(cli.threads)?;
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)?;
    writer.set_threads(cli.threads)?;
    for (ordinal, result) in reader.records().enumerate() {
        let mut record = result?;
        let key = if !record.is_unmapped() && !record.is_secondary() && !record.is_supplementary() {
            let cell = aux_string(&record, b"CB").and_then(|v| Dna::pack(&v));
            let umi = aux_string(&record, b"UR")
                .or_else(|| aux_string(&record, b"UB"))
                .and_then(|v| Dna::pack(&v));
            let gene = choose_gene(&record).and_then(|(g, _)| gene_ids.get(&g).copied());
            match (cell, gene, umi) {
                (Some(cell), Some(gene), Some(umi)) => Some(RawKey { cell, gene, umi }),
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
            let representative = decision.representative == ordinal as u64
                && matches!(decision.outcome, Outcome::Unique | Outcome::Won);
            record.push_aux(b"xf", Aux::I32(if representative { 25 } else { 17 }))?;
            record.push_aux(b"XV", Aux::U32(decision.collision_degree))?;
            record.push_aux(
                b"XW",
                Aux::Char(match decision.outcome {
                    Outcome::Unique => b'U',
                    Outcome::Won => b'W',
                    Outcome::Lost => b'L',
                    Outcome::Tie => b'T',
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

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(tracing_subscriber::EnvFilter::from_default_env())
        .init();
    let cli = if std::env::var_os("SBXPRESS_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("sbxpress molecules").get_matches())?
    } else {
        Cli::parse()
    };
    if cli.bam == cli.output {
        bail!("input and output BAM must differ");
    }
    let started = Instant::now();
    let mut metrics = Metrics {
        input_bam: cli.bam.display().to_string(),
        output_bam: cli.output.display().to_string(),
        ..Default::default()
    };
    let (evidence, genes, gene_ids) = collect_evidence(&cli.bam, cli.threads, &mut metrics)?;
    if cli.saturation_interval == 0 {
        bail!("--saturation-interval must be greater than zero");
    }
    let (decisions, saturation) = build_decisions(
        &evidence,
        &mut metrics,
        cli.saturation.as_ref().map(|_| cli.saturation_interval),
    );
    if let (Some(curve), Some(path)) = (saturation, cli.saturation.as_deref()) {
        metrics.saturation = Some(curve.write(path)?);
    }
    write_output(&cli, &decisions, &genes, &gene_ids, &mut metrics)?;
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
            gene: 0,
            umi: major,
        };
        let gene0_error = RawKey {
            cell,
            gene: 0,
            umi: error,
        };
        let gene1_major = RawKey {
            cell,
            gene: 1,
            umi: major,
        };
        input.insert(gene0_major, evidence(10, 0));
        input.insert(gene0_error, evidence(2, 1));
        input.insert(gene1_major, evidence(3, 2));

        let mut metrics = Metrics::default();
        let (decisions, _) = build_decisions(&input, &mut metrics, None);
        assert_eq!(decisions[&gene0_error].corrected, major);
        assert_eq!(decisions[&gene0_major].outcome, Outcome::Won);
        assert_eq!(decisions[&gene0_error].outcome, Outcome::Won);
        assert_eq!(decisions[&gene1_major].outcome, Outcome::Lost);
        assert_eq!(metrics.umi_nodes_collapsed, 1);
        assert_eq!(metrics.collision_groups, 1);
        assert_eq!(metrics.molecules_counted, 1);
    }

    #[test]
    fn equal_cross_gene_support_is_not_counted() {
        let cell = Dna::pack("AACCGGTTAACCGGTT").unwrap();
        let umi = Dna::pack("ACGTACGTACGT").unwrap();
        let a = RawKey { cell, gene: 0, umi };
        let b = RawKey { cell, gene: 1, umi };
        let input = HashMap::from([(a, evidence(2, 0)), (b, evidence(2, 1))]);
        let mut metrics = Metrics::default();
        let (decisions, _) = build_decisions(&input, &mut metrics, None);
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
}
