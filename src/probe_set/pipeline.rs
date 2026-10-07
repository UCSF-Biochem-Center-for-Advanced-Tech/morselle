use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufWriter, Write as _};
use std::path::Path;
use std::time::Instant;

use anyhow::{ensure, Context, Result};
use needletail::parse_fastx_file;
use rayon::prelude::*;
use rust_htslib::bam::{
    self, header::HeaderRecord, record::Aux, record::Cigar, record::CigarString,
};
use tracing::info;

use super::aligner::ProbeMatcher;
use super::types::{Cli, Metrics, ProbeSetReference, SYNTHETIC_CONTIG_LEN};

pub fn parse_tag<'a>(comment: &'a str, tag: &str) -> Option<&'a str> {
    let start = comment.find(tag)? + tag.len();
    let rest = &comment[start..];
    let end = rest.find('\t').unwrap_or(rest.len());
    Some(&rest[..end])
}

fn genes_manifest_bytes(genes: &[(&str, &str)]) -> Vec<u8> {
    let mut out = Vec::new();
    for (id, name) in genes {
        out.extend_from_slice(id.as_bytes());
        out.push(b'\t');
        out.extend_from_slice(name.as_bytes());
        out.extend_from_slice(b"\tGene Expression\n");
    }
    out
}

pub fn establish_genes_manifest(path: &Path, genes: &[(&str, &str)]) -> Result<()> {
    let expected = genes_manifest_bytes(genes);
    if path.exists() {
        let observed = fs::read(path).with_context(|| format!("reading {}", path.display()))?;
        ensure!(
            observed == expected,
            "existing gene manifest {} does not match this probe set; use a matching CSV or choose a new path",
            path.display()
        );
        return Ok(());
    }
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, &expected).with_context(|| format!("writing {}", path.display()))
}

/// One FASTQ record -> its BAM record (mapped to the matched probe's gene
/// contig, or unmapped), plus the matched gene index for metrics.
fn build_record(
    header_line: &str,
    seq: &[u8],
    qual: &[u8],
    matcher: &ProbeMatcher,
    reference: &ProbeSetReference,
    gene_index: &HashMap<&str, usize>,
    gene_name: &HashMap<&str, &str>,
) -> Result<(bam::Record, Option<usize>)> {
    let (qname, comment) = header_line.split_once('\t').unwrap_or((header_line, ""));
    let mut record = bam::Record::new();
    let mut matched_gene_idx = None;
    if let Some(probe_idx) = matcher.locate(seq) {
        let gene_id = reference.probes[probe_idx].gene_id.as_str();
        let g_idx = gene_index[gene_id];
        matched_gene_idx = Some(g_idx);
        let cigar = CigarString(vec![Cigar::Match(seq.len() as u32)]);
        record.set(qname.as_bytes(), Some(&cigar), seq, qual);
        record.set_tid(g_idx as i32);
        record.set_pos(0);
        record.set_mapq(60);
        record.unset_unmapped();
        record.push_aux(b"GX", Aux::String(gene_id))?;
        record.push_aux(b"GN", Aux::String(gene_name[gene_id]))?;
    } else {
        record.set(qname.as_bytes(), None, seq, qual);
    }
    for tag in ["CR", "CB", "UR", "UB", "PR", "PB"] {
        if let Some(value) = parse_tag(comment, &format!("{tag}:Z:")) {
            if !value.is_empty() {
                record.push_aux(tag.as_bytes(), Aux::String(value))?;
            }
        }
    }
    Ok((record, matched_gene_idx))
}

pub fn run(cli: &Cli) -> Result<()> {
    let started = Instant::now();
    let reference = ProbeSetReference::load(&cli.probe_set)?;
    info!(probes = reference.probes.len(), path = %cli.probe_set.display(), "loaded probe set");
    let matcher = ProbeMatcher::build(&reference.probes, cli.max_mm, cli.window);

    // Distinct genes, first-seen order, for the synthetic contig header and
    // the (optional) canonical gene manifest.
    let mut gene_order: Vec<&str> = Vec::new();
    let mut gene_index: HashMap<&str, usize> = HashMap::new();
    let mut gene_name: HashMap<&str, &str> = HashMap::new();
    for p in &reference.probes {
        gene_index.entry(p.gene_id.as_str()).or_insert_with(|| {
            gene_order.push(p.gene_id.as_str());
            gene_order.len() - 1
        });
        gene_name.insert(p.gene_id.as_str(), p.gene_name.as_str());
    }

    if let Some(path) = &cli.genes_out {
        let genes: Vec<(&str, &str)> = gene_order.iter().map(|&id| (id, gene_name[id])).collect();
        establish_genes_manifest(path, &genes)?;
    }

    let mut header = bam::Header::new();
    for &gene_id in &gene_order {
        header.push_record(
            HeaderRecord::new(b"SQ")
                .push_tag(b"SN", gene_id)
                .push_tag(b"LN", SYNTHETIC_CONTIG_LEN),
        );
    }
    if let Some(parent) = cli.output.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    let mut writer = bam::Writer::from_path(&cli.output, &header, bam::Format::Bam)
        .with_context(|| format!("creating {}", cli.output.display()))?;

    let pool = if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build()
            .context("building thread pool")?
    } else {
        rayon::ThreadPoolBuilder::new()
            .build()
            .context("building thread pool")?
    };
    // BGZF compression was single-threaded and becomes the bottleneck once
    // matching is parallel; give htslib the same worker budget.
    writer.set_threads(pool.current_num_threads())?;

    let mut reads = 0u64;
    let mut matched = 0u64;
    let mut genes_hit: std::collections::HashSet<usize> = std::collections::HashSet::new();

    let mut fastq = parse_fastx_file(&cli.r2_fastq)
        .with_context(|| format!("opening {}", cli.r2_fastq.display()))?;

    const BATCH_SIZE: usize = 10_000;
    let mut batch: Vec<(String, Vec<u8>, Vec<u8>)> = Vec::with_capacity(BATCH_SIZE);
    let mut finished = false;

    while !finished {
        batch.clear();
        while batch.len() < BATCH_SIZE {
            if let Some(rec) = fastq.next() {
                let rec = rec.context("reading R2 FASTQ record")?;
                reads += 1;
                let header_line = String::from_utf8_lossy(rec.id()).to_string();
                let seq = rec.seq().to_vec();
                let qual = rec.qual().unwrap_or(&[]).to_vec();
                batch.push((header_line, seq, qual));
            } else {
                finished = true;
                break;
            }
        }

        if batch.is_empty() {
            break;
        }

        // `map` + `collect` (not a channel) keeps the output BAM in input
        // FASTQ order, so runs are reproducible and downstream ordinals
        // (molecules' representative tiebreak) don't depend on scheduling.
        let built: Vec<(bam::Record, Option<usize>)> = pool.install(|| {
            batch
                .par_iter()
                .map(|(header_line, seq, qual)| {
                    build_record(
                        header_line,
                        seq,
                        qual,
                        &matcher,
                        &reference,
                        &gene_index,
                        &gene_name,
                    )
                })
                .collect::<Result<Vec<_>>>()
        })?;
        for (record, matched_gene_idx) in &built {
            if let Some(g_idx) = matched_gene_idx {
                matched += 1;
                genes_hit.insert(*g_idx);
            }
            writer.write(record)?;
        }
    }
    drop(writer);

    let metrics = Metrics {
        reads,
        matched,
        hit_rate_pct: if reads > 0 {
            100.0 * matched as f64 / reads as f64
        } else {
            0.0
        },
        distinct_genes_hit: genes_hit.len(),
    };
    if let Some(path) = &cli.metrics {
        let mut w = BufWriter::new(File::create(path)?);
        w.write_all(serde_json::to_string_pretty(&metrics)?.as_bytes())?;
    }

    info!(
        reads = metrics.reads,
        matched = metrics.matched,
        hit_rate_pct = format!("{:.1}", metrics.hit_rate_pct),
        distinct_genes_hit = metrics.distinct_genes_hit,
        elapsed = format!("{:.1}s", started.elapsed().as_secs_f64()),
        out = %cli.output.display(),
        "RUN SUMMARY"
    );
    Ok(())
}
