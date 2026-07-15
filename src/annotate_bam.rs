//! Add Cell Ranger-compatible GENCODE annotation tags to a genomic BAM.
//!
//! This stage performs feature assignment only. It deliberately does not correct
//! UMIs or set `xf`: those require evidence aggregated across all reads sharing a
//! cell/gene/UMI neighborhood and belong in the subsequent molecule-calling stage.

use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read as IoRead};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, Read};
use serde::Serialize;
use tracing::info;

const BIN_SIZE: i64 = 16_384;

#[derive(Debug, Parser)]
#[command(
    name = "sbx10x-annotate",
    version,
    about = "Add Cell Ranger-compatible GX/GN/TX/AN/RE tags to a genomic BAM"
)]
pub struct Cli {
    /// Genomic BAM/SAM produced by minimap2, or `-` for standard input.
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,

    /// GENCODE-compatible GTF, optionally gzip-compressed.
    #[arg(short, long, value_name = "GTF")]
    pub gtf: PathBuf,

    /// Tagged BAM output, or `-` for BAM on standard output.
    #[arg(short, long, value_name = "BAM")]
    pub output: PathBuf,

    /// Count same-strand gene-body overlaps as intronic assignments.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub include_introns: bool,

    /// BAM decompression/compression threads.
    #[arg(short, long, default_value_t = 4)]
    pub threads: usize,

    /// Annotate records in parallel and permit output reordering within batches.
    /// Intended for pipelines that immediately coordinate-sort the output.
    #[arg(long)]
    pub unordered: bool,

    /// Records per parallel annotation batch when --unordered is enabled.
    #[arg(long, default_value_t = 8192)]
    pub batch_size: usize,

    /// Optional JSON metrics path (default: OUTPUT.metrics.json).
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    /// Canonical gene feature manifest. Creates it from the filtered GTF when
    /// absent; when present, verifies byte-for-byte identity. Pass the same file
    /// to `sbxpress matrix --features` for stable rows across experiments.
    #[arg(long, value_name = "GENES.MTX")]
    pub genes_mtx: Option<PathBuf>,

    /// Apply the 10x GENCODE gene/transcript biotype, readthrough, and PAR rules.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub cellranger_filter: bool,

    /// Remove numeric version suffixes from Ensembl gene/transcript IDs.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub strip_id_versions: bool,

    /// Use terminal-exon support as a non-destructive preference signal for a 3' assay.
    #[arg(long)]
    pub three_prime: bool,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Strand {
    Plus,
    Minus,
}

impl Strand {
    fn parse(value: &str) -> Option<Self> {
        match value {
            "+" => Some(Self::Plus),
            "-" => Some(Self::Minus),
            _ => None,
        }
    }

    fn matches_record(self, reverse: bool) -> bool {
        matches!((self, reverse), (Self::Plus, false) | (Self::Minus, true))
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct Interval {
    start: i64,
    end: i64,
}

impl Interval {
    fn overlaps(self, other: Self) -> bool {
        self.start < other.end && other.start < self.end
    }

    fn contains(self, other: Self) -> bool {
        self.start <= other.start && self.end >= other.end
    }
}

#[derive(Debug)]
struct Gene {
    id: String,
    name: String,
    contig: usize,
    strand: Strand,
    span: Interval,
}

#[derive(Debug)]
struct Transcript {
    id: String,
    gene: usize,
    contig: usize,
    strand: Strand,
    exons: Vec<Interval>, // genomic ascending order
}

#[derive(Debug, Default)]
struct ContigIndex {
    exon_bins: Vec<Vec<usize>>, // transcript indexes
    gene_bins: Vec<Vec<usize>>,
}

#[derive(Debug, Default)]
struct Annotation {
    contig_names: Vec<String>,
    contig_ids: HashMap<String, usize>,
    genes: Vec<Gene>,
    transcripts: Vec<Transcript>,
    indexes: Vec<ContigIndex>,
    exons_seen: u64,
    exons_filtered_biotype: u64,
    exons_filtered_tag: u64,
}

#[derive(Debug)]
struct Call {
    region: u8,
    genes: Vec<usize>,
    tx: Vec<String>,
    antisense: Vec<String>,
    terminal_genes: Vec<usize>,
    terminal_tx: Vec<String>,
}

#[derive(Debug, Default, Serialize)]
struct Metrics {
    input_bam: String,
    output_bam: String,
    gtf: String,
    genes_mtx: Option<String>,
    include_introns: bool,
    cellranger_filter: bool,
    strip_id_versions: bool,
    three_prime: bool,
    unordered: bool,
    elapsed_secs: f64,
    records: u64,
    primary: u64,
    supplementary: u64,
    unmapped: u64,
    exonic: u64,
    intronic: u64,
    intergenic: u64,
    gene_assigned: u64,
    unique_gene: u64,
    multi_gene: u64,
    transcript_assigned: u64,
    antisense: u64,
    collision_records: u64,
    gene_degree_0: u64,
    gene_degree_1: u64,
    gene_degree_2: u64,
    gene_degree_3plus: u64,
    transcript_degree_0: u64,
    transcript_degree_1: u64,
    transcript_degree_2: u64,
    transcript_degree_3plus: u64,
    terminal_resolved_gene_collision: u64,
    terminal_unresolved_gene_collision: u64,
    terminal_unsupported_gene_collision: u64,
    genes_loaded: usize,
    transcripts_loaded: usize,
    gtf_exons_seen: u64,
    gtf_exons_filtered_biotype: u64,
    gtf_exons_filtered_tag: u64,
}

fn genes_mtx_bytes(annotation: &Annotation) -> Vec<u8> {
    let mut output = Vec::new();
    for gene in &annotation.genes {
        output.extend_from_slice(gene.id.as_bytes());
        output.push(b'\t');
        output.extend_from_slice(gene.name.as_bytes());
        output.extend_from_slice(b"\tGene Expression\n");
    }
    output
}

fn establish_genes_mtx(path: &Path, annotation: &Annotation) -> Result<()> {
    let expected = genes_mtx_bytes(annotation);
    if path.exists() {
        let observed = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        if observed != expected {
            bail!(
                "existing gene manifest {} does not match the filtered annotation; use a matching GTF/rule set or choose a new path",
                path.display()
            );
        }
        info!(path = %path.display(), genes = annotation.genes.len(), "verified canonical gene manifest");
        return Ok(());
    }
    if let Some(parent) = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
    {
        fs::create_dir_all(parent)?;
    }
    let temporary = PathBuf::from(format!("{}.tmp.{}", path.display(), std::process::id()));
    fs::write(&temporary, &expected).with_context(|| format!("writing {}", temporary.display()))?;
    fs::rename(&temporary, path).with_context(|| format!("creating {}", path.display()))?;
    info!(path = %path.display(), genes = annotation.genes.len(), "created canonical gene manifest");
    Ok(())
}

fn text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let file = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let input: Box<dyn IoRead> = if path.extension().is_some_and(|ext| ext == "gz") {
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    Ok(Box::new(BufReader::new(input)))
}

fn attribute<'a>(attrs: &'a str, key: &str) -> Option<&'a str> {
    attrs.split(';').find_map(|field| {
        let field = field.trim();
        let (k, value) = field.split_once(char::is_whitespace)?;
        (k == key).then(|| value.trim().trim_matches('"'))
    })
}

fn has_attribute_value(attrs: &str, key: &str, wanted: &str) -> bool {
    attrs.split(';').any(|field| {
        let field = field.trim();
        let Some((k, value)) = field.split_once(char::is_whitespace) else {
            return false;
        };
        k == key && value.trim().trim_matches('"') == wanted
    })
}

fn allowed_biotype(value: &str) -> bool {
    matches!(
        value,
        "protein_coding"
            | "protein_coding_LoF"
            | "lncRNA"
            | "antisense"
            | "IG_C_gene"
            | "IG_D_gene"
            | "IG_J_gene"
            | "IG_LV_gene"
            | "IG_V_gene"
            | "IG_V_pseudogene"
            | "IG_J_pseudogene"
            | "IG_C_pseudogene"
            | "TR_C_gene"
            | "TR_D_gene"
            | "TR_J_gene"
            | "TR_V_gene"
            | "TR_V_pseudogene"
            | "TR_J_pseudogene"
    )
}

fn normalize_ensembl_id(value: &str, strip_version: bool) -> String {
    if strip_version && (value.starts_with("ENSG") || value.starts_with("ENST")) {
        if let Some((base, version)) = value.rsplit_once('.') {
            if version.chars().all(|c| c.is_ascii_digit()) {
                return base.to_owned();
            }
        }
    }
    value.to_owned()
}

impl Annotation {
    fn load(path: &Path, cellranger_filter: bool, strip_id_versions: bool) -> Result<Self> {
        let mut a = Self::default();
        let mut genes_by_id: HashMap<String, usize> = HashMap::new();
        let mut tx_by_id: HashMap<String, usize> = HashMap::new();

        for (line_no, result) in text_reader(path)?.lines().enumerate() {
            let line = result.with_context(|| format!("reading GTF line {}", line_no + 1))?;
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != 9 || fields[2] != "exon" {
                continue;
            }
            a.exons_seen += 1;
            let Some(strand) = Strand::parse(fields[6]) else {
                continue;
            };
            let start = fields[3]
                .parse::<i64>()
                .with_context(|| format!("invalid GTF start at line {}", line_no + 1))?
                - 1;
            let end = fields[4]
                .parse::<i64>()
                .with_context(|| format!("invalid GTF end at line {}", line_no + 1))?;
            if start < 0 || end <= start {
                bail!("invalid GTF interval at line {}", line_no + 1);
            }
            let attrs = fields[8];
            if cellranger_filter {
                let gene_type =
                    attribute(attrs, "gene_type").or_else(|| attribute(attrs, "gene_biotype"));
                let transcript_type = attribute(attrs, "transcript_type")
                    .or_else(|| attribute(attrs, "transcript_biotype"));
                if !gene_type.is_some_and(allowed_biotype)
                    || !transcript_type.is_some_and(allowed_biotype)
                {
                    a.exons_filtered_biotype += 1;
                    continue;
                }
                if has_attribute_value(attrs, "tag", "readthrough_transcript")
                    || has_attribute_value(attrs, "tag", "PAR")
                {
                    a.exons_filtered_tag += 1;
                    continue;
                }
            }
            let gene_id = attribute(attrs, "gene_id")
                .with_context(|| format!("exon lacks gene_id at line {}", line_no + 1))?;
            let gene_name = attribute(attrs, "gene_name").unwrap_or(gene_id);
            let tx_id = attribute(attrs, "transcript_id")
                .with_context(|| format!("exon lacks transcript_id at line {}", line_no + 1))?;
            let gene_id = normalize_ensembl_id(gene_id, strip_id_versions);
            let tx_id = normalize_ensembl_id(tx_id, strip_id_versions);
            let contig = if let Some(&id) = a.contig_ids.get(fields[0]) {
                id
            } else {
                let id = a.contig_names.len();
                a.contig_names.push(fields[0].to_owned());
                a.contig_ids.insert(fields[0].to_owned(), id);
                a.indexes.push(ContigIndex::default());
                id
            };
            let exon = Interval { start, end };
            let gene = if let Some(&index) = genes_by_id.get(&gene_id) {
                let g = &mut a.genes[index];
                if g.contig != contig || g.strand != strand {
                    bail!("gene {gene_id} occurs on inconsistent contig/strand");
                }
                g.span.start = g.span.start.min(start);
                g.span.end = g.span.end.max(end);
                index
            } else {
                let index = a.genes.len();
                a.genes.push(Gene {
                    id: gene_id.clone(),
                    name: gene_name.to_owned(),
                    contig,
                    strand,
                    span: exon,
                });
                genes_by_id.insert(gene_id.clone(), index);
                index
            };
            let tx = if let Some(&index) = tx_by_id.get(&tx_id) {
                index
            } else {
                let index = a.transcripts.len();
                a.transcripts.push(Transcript {
                    id: tx_id.clone(),
                    gene,
                    contig,
                    strand,
                    exons: Vec::new(),
                });
                tx_by_id.insert(tx_id, index);
                index
            };
            a.transcripts[tx].exons.push(exon);
        }

        for tx in &mut a.transcripts {
            tx.exons.sort_unstable_by_key(|e| (e.start, e.end));
            tx.exons.dedup();
        }
        for (tx_index, tx) in a.transcripts.iter().enumerate() {
            for &exon in &tx.exons {
                add_to_bins(&mut a.indexes[tx.contig].exon_bins, exon, tx_index);
            }
        }
        for (gene_index, gene) in a.genes.iter().enumerate() {
            add_to_bins(&mut a.indexes[gene.contig].gene_bins, gene.span, gene_index);
        }
        Ok(a)
    }

    fn call(
        &self,
        contig_name: &str,
        blocks: &[Interval],
        reverse: bool,
        cigar: &str,
        include_introns: bool,
        three_prime: bool,
    ) -> Call {
        let Some(&contig) = self.contig_ids.get(contig_name) else {
            return Call {
                region: b'I',
                genes: vec![],
                tx: vec![],
                antisense: vec![],
                terminal_genes: vec![],
                terminal_tx: vec![],
            };
        };
        if blocks.is_empty() {
            return Call {
                region: b'I',
                genes: vec![],
                tx: vec![],
                antisense: vec![],
                terminal_genes: vec![],
                terminal_tx: vec![],
            };
        }
        let span = Interval {
            start: blocks[0].start,
            end: blocks.last().unwrap().end,
        };
        let mut tx_candidates = HashSet::new();
        collect_bins(&self.indexes[contig].exon_bins, span, &mut tx_candidates);
        let mut sense_exon_genes = HashSet::new();
        let mut any_exon_overlap = false;
        let mut sense_tx = Vec::new();
        let mut antisense_tx = Vec::new();
        let mut terminal_tx = Vec::new();
        let mut terminal_genes = HashSet::new();

        for tx_index in tx_candidates {
            let tx = &self.transcripts[tx_index];
            let sense = tx.strand.matches_record(reverse);
            if tx
                .exons
                .iter()
                .any(|&e| blocks.iter().any(|&b| e.overlaps(b)))
            {
                any_exon_overlap = true;
            }
            if let Some(pos) = compatible_position(tx, blocks) {
                let entry = format!(
                    "{},{}{},{}",
                    tx.id,
                    if sense { '+' } else { '-' },
                    pos,
                    transcript_cigar(cigar, tx.strand)
                );
                if sense {
                    sense_exon_genes.insert(tx.gene);
                    if three_prime && overlaps_terminal_exon(tx, blocks) {
                        terminal_genes.insert(tx.gene);
                        terminal_tx.push(entry.clone());
                    }
                    sense_tx.push(entry);
                } else {
                    antisense_tx.push(entry);
                }
            }
        }

        sense_tx.sort_unstable();
        sense_tx.dedup();
        antisense_tx.sort_unstable();
        antisense_tx.dedup();
        terminal_tx.sort_unstable();
        terminal_tx.dedup();

        let mut body_candidates = HashSet::new();
        collect_bins(&self.indexes[contig].gene_bins, span, &mut body_candidates);
        let mut sense_body = Vec::new();
        let mut antisense_body = Vec::new();
        for gene_index in body_candidates {
            let gene = &self.genes[gene_index];
            if !gene.span.overlaps(span) {
                continue;
            }
            if gene.strand.matches_record(reverse) {
                sense_body.push(gene_index);
            } else {
                antisense_body.push(gene_index);
            }
        }

        let exonic = any_exon_overlap;
        let intronic = !exonic && (!sense_body.is_empty() || !antisense_body.is_empty());
        let mut genes: Vec<usize> = if !sense_exon_genes.is_empty() {
            sense_exon_genes.into_iter().collect()
        } else if intronic && include_introns && !sense_body.is_empty() {
            sense_body.clone()
        } else {
            Vec::new()
        };
        genes.sort_unstable_by(|&x, &y| self.genes[x].id.cmp(&self.genes[y].id));
        genes.dedup();
        let mut terminal_genes: Vec<usize> = terminal_genes.into_iter().collect();
        terminal_genes.sort_unstable_by(|&x, &y| self.genes[x].id.cmp(&self.genes[y].id));
        terminal_genes.dedup();

        if sense_tx.is_empty() && intronic && include_introns {
            sense_tx = genes
                .iter()
                .map(|&g| format!("{},+", self.genes[g].id))
                .collect();
        }
        if antisense_tx.is_empty() && intronic && !antisense_body.is_empty() && include_introns {
            antisense_body.sort_unstable_by(|&x, &y| self.genes[x].id.cmp(&self.genes[y].id));
            antisense_body.dedup();
            antisense_tx = antisense_body
                .iter()
                .map(|&g| format!("{},-", self.genes[g].id))
                .collect();
        }
        Call {
            region: if exonic {
                b'E'
            } else if intronic {
                b'N'
            } else {
                b'I'
            },
            genes,
            tx: sense_tx,
            antisense: antisense_tx,
            terminal_genes,
            terminal_tx,
        }
    }
}

fn add_to_bins(bins: &mut Vec<Vec<usize>>, interval: Interval, value: usize) {
    let first = (interval.start / BIN_SIZE) as usize;
    let last = ((interval.end - 1) / BIN_SIZE) as usize;
    if bins.len() <= last {
        bins.resize_with(last + 1, Vec::new);
    }
    for bin in &mut bins[first..=last] {
        bin.push(value);
    }
}

fn collect_bins(bins: &[Vec<usize>], interval: Interval, output: &mut HashSet<usize>) {
    let first = (interval.start / BIN_SIZE).max(0) as usize;
    let last = ((interval.end - 1) / BIN_SIZE).max(0) as usize;
    for bin in first..=last.min(bins.len().saturating_sub(1)) {
        output.extend(bins[bin].iter().copied());
    }
}

fn compatible_position(tx: &Transcript, blocks: &[Interval]) -> Option<i64> {
    let mut exon_indexes = Vec::with_capacity(blocks.len());
    for &block in blocks {
        exon_indexes.push(tx.exons.iter().position(|&exon| exon.contains(block))?);
    }
    for (pair, exon_pair) in blocks.windows(2).zip(exon_indexes.windows(2)) {
        if exon_pair[1] != exon_pair[0] + 1
            || pair[0].end != tx.exons[exon_pair[0]].end
            || pair[1].start != tx.exons[exon_pair[1]].start
        {
            return None;
        }
    }
    let first = exon_indexes[0];
    let last = *exon_indexes.last()?;
    Some(match tx.strand {
        Strand::Plus => {
            tx.exons[..first]
                .iter()
                .map(|e| e.end - e.start)
                .sum::<i64>()
                + blocks[0].start
                - tx.exons[first].start
        }
        Strand::Minus => {
            tx.exons[last + 1..]
                .iter()
                .map(|e| e.end - e.start)
                .sum::<i64>()
                + tx.exons[last].end
                - blocks.last()?.end
        }
    })
}

fn overlaps_terminal_exon(tx: &Transcript, blocks: &[Interval]) -> bool {
    let terminal = match tx.strand {
        Strand::Plus => tx.exons.last(),
        Strand::Minus => tx.exons.first(),
    };
    terminal.is_some_and(|&exon| blocks.iter().any(|&block| exon.overlaps(block)))
}

fn transcript_cigar(cigar: &str, strand: Strand) -> String {
    let mut ops: Vec<(u64, char)> = Vec::new();
    let mut digits = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            digits.push(ch);
            continue;
        }
        if ch != 'N' {
            let length = digits.parse::<u64>().unwrap_or(0);
            if let Some((previous_length, previous_op)) = ops.last_mut() {
                if *previous_op == ch {
                    *previous_length += length;
                } else {
                    ops.push((length, ch));
                }
            } else {
                ops.push((length, ch));
            }
        }
        digits.clear();
    }
    if strand == Strand::Minus {
        ops.reverse();
    }
    ops.into_iter().map(|(n, op)| format!("{n}{op}")).collect()
}

fn alignment_blocks(record: &bam::Record) -> Vec<Interval> {
    use rust_htslib::bam::record::Cigar;
    let mut pos = record.pos();
    let mut start = pos;
    let mut blocks = Vec::new();
    for op in record.cigar().iter() {
        match *op {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) | Cigar::Del(n) => pos += n as i64,
            Cigar::RefSkip(n) => {
                if pos > start {
                    blocks.push(Interval { start, end: pos });
                }
                pos += n as i64;
                start = pos;
            }
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    if pos > start {
        blocks.push(Interval { start, end: pos });
    }
    blocks
}

fn replace_string_tag(record: &mut bam::Record, tag: &[u8; 2], value: Option<&str>) -> Result<()> {
    let _ = record.remove_aux(tag);
    if let Some(value) = value {
        record.push_aux(tag, Aux::String(value))?;
    }
    Ok(())
}

fn tag_record(
    record: &mut bam::Record,
    call: &Call,
    annotation: &Annotation,
    three_prime: bool,
) -> Result<()> {
    let _ = record.remove_aux(b"RE");
    record.push_aux(b"RE", Aux::Char(call.region))?;
    let gx = (!call.genes.is_empty()).then(|| {
        call.genes
            .iter()
            .map(|&g| annotation.genes[g].id.as_str())
            .collect::<Vec<_>>()
            .join(";")
    });
    let gn = (!call.genes.is_empty()).then(|| {
        call.genes
            .iter()
            .map(|&g| annotation.genes[g].name.as_str())
            .collect::<Vec<_>>()
            .join(";")
    });
    let tx = (!call.tx.is_empty()).then(|| call.tx.join(";"));
    let an = (!call.antisense.is_empty()).then(|| call.antisense.join(";"));
    let zg = (!call.terminal_genes.is_empty()).then(|| {
        call.terminal_genes
            .iter()
            .map(|&g| annotation.genes[g].id.as_str())
            .collect::<Vec<_>>()
            .join(";")
    });
    let zt = (!call.terminal_tx.is_empty()).then(|| call.terminal_tx.join(";"));
    replace_string_tag(record, b"GX", gx.as_deref())?;
    replace_string_tag(record, b"GN", gn.as_deref())?;
    replace_string_tag(record, b"TX", tx.as_deref())?;
    replace_string_tag(record, b"AN", an.as_deref())?;
    replace_string_tag(
        record,
        b"ZG",
        three_prime.then_some(zg.as_deref()).flatten(),
    )?;
    replace_string_tag(
        record,
        b"ZT",
        three_prime.then_some(zt.as_deref()).flatten(),
    )?;
    for tag in [b"XC", b"XT", b"XP"] {
        let _ = record.remove_aux(tag);
    }
    record.push_aux(b"XC", Aux::U32(call.genes.len() as u32))?;
    record.push_aux(b"XT", Aux::U32(call.tx.len() as u32))?;
    if three_prime {
        record.push_aux(b"XP", Aux::U32(call.terminal_genes.len() as u32))?;
    }
    Ok(())
}

fn increment_degree(
    value: usize,
    zero: &mut u64,
    one: &mut u64,
    two: &mut u64,
    three_plus: &mut u64,
) {
    match value {
        0 => *zero += 1,
        1 => *one += 1,
        2 => *two += 1,
        _ => *three_plus += 1,
    }
}

#[derive(Debug, Default)]
struct RecordMetrics {
    records: u64,
    primary: u64,
    supplementary: u64,
    unmapped: u64,
    exonic: u64,
    intronic: u64,
    intergenic: u64,
    gene_assigned: u64,
    unique_gene: u64,
    multi_gene: u64,
    transcript_assigned: u64,
    antisense: u64,
    collision_records: u64,
    gene_degree_0: u64,
    gene_degree_1: u64,
    gene_degree_2: u64,
    gene_degree_3plus: u64,
    transcript_degree_0: u64,
    transcript_degree_1: u64,
    transcript_degree_2: u64,
    transcript_degree_3plus: u64,
    terminal_resolved_gene_collision: u64,
    terminal_unresolved_gene_collision: u64,
    terminal_unsupported_gene_collision: u64,
}

impl Metrics {
    fn add_record_metrics(&mut self, value: RecordMetrics) {
        self.records += value.records;
        self.primary += value.primary;
        self.supplementary += value.supplementary;
        self.unmapped += value.unmapped;
        self.exonic += value.exonic;
        self.intronic += value.intronic;
        self.intergenic += value.intergenic;
        self.gene_assigned += value.gene_assigned;
        self.unique_gene += value.unique_gene;
        self.multi_gene += value.multi_gene;
        self.transcript_assigned += value.transcript_assigned;
        self.antisense += value.antisense;
        self.collision_records += value.collision_records;
        self.gene_degree_0 += value.gene_degree_0;
        self.gene_degree_1 += value.gene_degree_1;
        self.gene_degree_2 += value.gene_degree_2;
        self.gene_degree_3plus += value.gene_degree_3plus;
        self.transcript_degree_0 += value.transcript_degree_0;
        self.transcript_degree_1 += value.transcript_degree_1;
        self.transcript_degree_2 += value.transcript_degree_2;
        self.transcript_degree_3plus += value.transcript_degree_3plus;
        self.terminal_resolved_gene_collision += value.terminal_resolved_gene_collision;
        self.terminal_unresolved_gene_collision += value.terminal_unresolved_gene_collision;
        self.terminal_unsupported_gene_collision += value.terminal_unsupported_gene_collision;
    }
}

fn annotate_record(
    mut record: bam::Record,
    annotation: &Annotation,
    contig_names: &[Vec<u8>],
    include_introns: bool,
    three_prime: bool,
) -> Result<(bam::Record, RecordMetrics)> {
    let mut metrics = RecordMetrics {
        records: 1,
        ..Default::default()
    };
    if record.is_unmapped() {
        metrics.unmapped = 1;
        replace_string_tag(&mut record, b"GX", None)?;
        replace_string_tag(&mut record, b"GN", None)?;
        replace_string_tag(&mut record, b"TX", None)?;
        replace_string_tag(&mut record, b"AN", None)?;
        replace_string_tag(&mut record, b"ZG", None)?;
        replace_string_tag(&mut record, b"ZT", None)?;
        for tag in [b"XC", b"XT", b"XP"] {
            let _ = record.remove_aux(tag);
        }
        let _ = record.remove_aux(b"RE");
        return Ok((record, metrics));
    }

    if record.is_supplementary() {
        metrics.supplementary = 1;
    } else {
        metrics.primary = 1;
    }
    let tid = record.tid();
    let contig = contig_names
        .get(tid as usize)
        .with_context(|| format!("record has invalid target id {tid}"))?;
    let contig = std::str::from_utf8(contig)?;
    let blocks = alignment_blocks(&record);
    let call = annotation.call(
        contig,
        &blocks,
        record.is_reverse(),
        &record.cigar().to_string(),
        include_introns,
        three_prime,
    );
    match call.region {
        b'E' => metrics.exonic = 1,
        b'N' => metrics.intronic = 1,
        _ => metrics.intergenic = 1,
    }
    if !call.genes.is_empty() {
        metrics.gene_assigned = 1;
        if call.genes.len() == 1 {
            metrics.unique_gene = 1;
        } else {
            metrics.multi_gene = 1;
        }
    }
    metrics.transcript_assigned = (!call.tx.is_empty()) as u64;
    metrics.antisense = (!call.antisense.is_empty()) as u64;
    if !record.is_supplementary() {
        metrics.collision_records = 1;
        increment_degree(
            call.genes.len(),
            &mut metrics.gene_degree_0,
            &mut metrics.gene_degree_1,
            &mut metrics.gene_degree_2,
            &mut metrics.gene_degree_3plus,
        );
        increment_degree(
            call.tx.len(),
            &mut metrics.transcript_degree_0,
            &mut metrics.transcript_degree_1,
            &mut metrics.transcript_degree_2,
            &mut metrics.transcript_degree_3plus,
        );
        if three_prime && call.genes.len() > 1 {
            match call.terminal_genes.len() {
                0 => metrics.terminal_unsupported_gene_collision = 1,
                1 => metrics.terminal_resolved_gene_collision = 1,
                _ => metrics.terminal_unresolved_gene_collision = 1,
            }
        }
    }
    tag_record(&mut record, &call, annotation, three_prime)?;
    Ok((record, metrics))
}

fn write_unordered_batch(
    batch: Vec<bam::Record>,
    pool: &rayon::ThreadPool,
    annotation: &Annotation,
    contig_names: &[Vec<u8>],
    include_introns: bool,
    three_prime: bool,
    writer: &mut bam::Writer,
    metrics: &mut Metrics,
) -> Result<()> {
    let (sender, receiver) = crossbeam_channel::unbounded();
    pool.install(|| {
        batch
            .into_par_iter()
            .for_each_with(sender, |sender, record| {
                // The receiver is alive until every worker finishes this batch.
                let _ = sender.send(annotate_record(
                    record,
                    annotation,
                    contig_names,
                    include_introns,
                    three_prime,
                ));
            });
    });
    for result in receiver {
        let (record, record_metrics) = result?;
        writer.write(&record)?;
        metrics.add_record_metrics(record_metrics);
    }
    Ok(())
}

pub fn run(cli: &Cli) -> Result<()> {
    let stdin = cli.bam == Path::new("-");
    let stdout = cli.output == Path::new("-");
    if cli.bam == cli.output && !stdin {
        bail!("input and output BAM must differ");
    }
    if stdout && cli.metrics.is_none() {
        bail!("--metrics is required when --output - is used");
    }
    if cli.threads == 0 {
        bail!("--threads must be greater than zero");
    }
    if cli.batch_size == 0 {
        bail!("--batch-size must be greater than zero");
    }
    let started = Instant::now();
    info!(gtf = %cli.gtf.display(), "loading annotation");
    let annotation = Annotation::load(&cli.gtf, cli.cellranger_filter, cli.strip_id_versions)?;
    info!(
        genes = annotation.genes.len(),
        transcripts = annotation.transcripts.len(),
        "annotation loaded"
    );
    if let Some(path) = &cli.genes_mtx {
        establish_genes_mtx(path, &annotation)?;
    }

    let mut reader = if stdin {
        bam::Reader::from_stdin().context("opening BAM/SAM from standard input")?
    } else {
        bam::Reader::from_path(&cli.bam)
            .with_context(|| format!("opening {}", cli.bam.display()))?
    };
    reader.set_threads(cli.threads)?;
    let header_view = reader.header().clone();
    let contig_names = (0..header_view.target_count())
        .map(|tid| header_view.tid2name(tid).to_vec())
        .collect::<Vec<_>>();
    let header = bam::Header::from_template(&header_view);
    let mut writer = if stdout {
        bam::Writer::from_stdout(&header, bam::Format::Bam)
            .context("opening BAM on standard output")?
    } else {
        bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)
            .with_context(|| format!("creating {}", cli.output.display()))?
    };
    writer.set_threads(cli.threads)?;

    let mut metrics = Metrics {
        input_bam: cli.bam.display().to_string(),
        output_bam: cli.output.display().to_string(),
        gtf: cli.gtf.display().to_string(),
        genes_mtx: cli
            .genes_mtx
            .as_ref()
            .map(|path| path.display().to_string()),
        include_introns: cli.include_introns,
        cellranger_filter: cli.cellranger_filter,
        strip_id_versions: cli.strip_id_versions,
        three_prime: cli.three_prime,
        unordered: cli.unordered,
        genes_loaded: annotation.genes.len(),
        transcripts_loaded: annotation.transcripts.len(),
        gtf_exons_seen: annotation.exons_seen,
        gtf_exons_filtered_biotype: annotation.exons_filtered_biotype,
        gtf_exons_filtered_tag: annotation.exons_filtered_tag,
        ..Default::default()
    };
    if cli.unordered {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build()
            .context("creating annotation worker pool")?;
        let mut batch = Vec::with_capacity(cli.batch_size);
        for result in reader.records() {
            batch.push(result.context("reading BAM record")?);
            if batch.len() == cli.batch_size {
                write_unordered_batch(
                    std::mem::take(&mut batch),
                    &pool,
                    &annotation,
                    &contig_names,
                    cli.include_introns,
                    cli.three_prime,
                    &mut writer,
                    &mut metrics,
                )?;
                batch.reserve(cli.batch_size);
            }
        }
        if !batch.is_empty() {
            write_unordered_batch(
                batch,
                &pool,
                &annotation,
                &contig_names,
                cli.include_introns,
                cli.three_prime,
                &mut writer,
                &mut metrics,
            )?;
        }
    } else {
        for result in reader.records() {
            let record = result.context("reading BAM record")?;
            let (record, record_metrics) = annotate_record(
                record,
                &annotation,
                &contig_names,
                cli.include_introns,
                cli.three_prime,
            )?;
            writer.write(&record)?;
            metrics.add_record_metrics(record_metrics);
        }
    }
    drop(writer);
    metrics.elapsed_secs = started.elapsed().as_secs_f64();
    let metrics_path = cli
        .metrics
        .clone()
        .unwrap_or_else(|| PathBuf::from(format!("{}.metrics.json", cli.output.display())));
    fs::write(&metrics_path, serde_json::to_string_pretty(&metrics)?)?;
    info!(records = metrics.records, elapsed = format!("{:.1}s", metrics.elapsed_secs), output = %cli.output.display(), "RUN SUMMARY");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::header::HeaderRecord;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    fn spliced_transcript_compatibility_and_position() {
        let tx = Transcript {
            id: "tx".into(),
            gene: 0,
            contig: 0,
            strand: Strand::Plus,
            exons: vec![
                Interval {
                    start: 100,
                    end: 150,
                },
                Interval {
                    start: 200,
                    end: 260,
                },
            ],
        };
        assert_eq!(
            compatible_position(
                &tx,
                &[
                    Interval {
                        start: 120,
                        end: 150
                    },
                    Interval {
                        start: 200,
                        end: 230
                    }
                ]
            ),
            Some(20)
        );
        assert_eq!(
            compatible_position(
                &tx,
                &[
                    Interval {
                        start: 120,
                        end: 140
                    },
                    Interval {
                        start: 200,
                        end: 230
                    }
                ]
            ),
            None
        );
    }

    #[test]
    fn transcript_cigar_drops_introns_and_reverses_minus() {
        assert_eq!(
            transcript_cigar("5S20M100N30M2I4M", Strand::Plus),
            "5S50M2I4M"
        );
        assert_eq!(
            transcript_cigar("5S20M100N30M2I4M", Strand::Minus),
            "4M2I50M5S"
        );
    }

    #[test]
    fn minus_position_uses_high_genomic_block() {
        let tx = Transcript {
            id: "tx".into(),
            gene: 0,
            contig: 0,
            strand: Strand::Minus,
            exons: vec![
                Interval {
                    start: 100,
                    end: 150,
                },
                Interval {
                    start: 200,
                    end: 260,
                },
            ],
        };
        assert_eq!(
            compatible_position(
                &tx,
                &[
                    Interval {
                        start: 120,
                        end: 150
                    },
                    Interval {
                        start: 200,
                        end: 230
                    },
                ],
            ),
            Some(30),
        );
    }

    #[test]
    fn unordered_annotation_matches_ordered_tags_and_counts() {
        let temp = tempfile::tempdir().unwrap();
        let gtf = temp.path().join("genes.gtf");
        fs::write(
            &gtf,
            "chr1\ttest\texon\t101\t200\t.\t+\t.\tgene_id \"ENSG000001.1\"; transcript_id \"ENST000001.1\"; gene_name \"GENE1\"; gene_type \"protein_coding\"; transcript_type \"protein_coding\";\n",
        )
        .unwrap();

        let input = temp.path().join("input.bam");
        let mut header = bam::Header::new();
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 1000),
        );
        let mut writer = bam::Writer::from_path(&input, &header, bam::Format::Bam).unwrap();
        let cigar = CigarString(vec![Cigar::Match(20)]);
        for index in 0..128 {
            let mut record = bam::Record::new();
            let name = format!("read-{index:03}");
            record.set(
                name.as_bytes(),
                Some(&cigar),
                b"ACGTACGTACGTACGTACGT",
                &[30; 20],
            );
            record.set_tid(0);
            record.set_pos(if index % 2 == 0 { 110 } else { 300 });
            writer.write(&record).unwrap();
        }
        drop(writer);

        let run_mode = |unordered: bool, stem: &str| {
            let output = temp.path().join(format!("{stem}.bam"));
            let metrics = temp.path().join(format!("{stem}.json"));
            run(&Cli {
                bam: input.clone(),
                gtf: gtf.clone(),
                output: output.clone(),
                include_introns: true,
                threads: 4,
                unordered,
                batch_size: 17,
                metrics: Some(metrics.clone()),
                genes_mtx: None,
                cellranger_filter: true,
                strip_id_versions: true,
                three_prime: true,
            })
            .unwrap();
            let mut reader = bam::Reader::from_path(output).unwrap();
            let mut tags = reader
                .records()
                .map(|record| {
                    let record = record.unwrap();
                    let gx = match record.aux(b"GX").ok() {
                        Some(Aux::String(value)) => value.to_owned(),
                        _ => String::new(),
                    };
                    let re = match record.aux(b"RE").ok() {
                        Some(Aux::Char(value)) => value,
                        _ => b'.',
                    };
                    (record.qname().to_vec(), gx, re)
                })
                .collect::<Vec<_>>();
            tags.sort_unstable();
            let metrics: serde_json::Value =
                serde_json::from_slice(&fs::read(metrics).unwrap()).unwrap();
            (tags, metrics)
        };

        let (ordered_tags, ordered_metrics) = run_mode(false, "ordered");
        let (unordered_tags, unordered_metrics) = run_mode(true, "unordered");
        assert_eq!(ordered_tags, unordered_tags);
        for field in [
            "records",
            "primary",
            "exonic",
            "intergenic",
            "gene_assigned",
            "unique_gene",
            "collision_records",
        ] {
            assert_eq!(ordered_metrics[field], unordered_metrics[field], "{field}");
        }
        assert_eq!(unordered_metrics["unordered"], true);
    }
}
