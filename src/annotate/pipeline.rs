//! Add 10x-compatible GENCODE annotation tags to a genomic BAM.
//!
//! This stage performs feature assignment only. It deliberately does not correct
//! UMIs or set `xf`: those require evidence aggregated across all reads sharing a
//! cell/gene/UMI neighborhood and belong in the subsequent molecule-calling stage.
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use rust_htslib::bam::{self, record::Aux, Read};
use tracing::info;

use super::barcode_metrics::observe;
use super::call::alignment_blocks;
use super::index::establish_genes_mtx;
use super::multimap::write_multimap_group;
use super::types::{Annotation, Call, Cli, Metrics, RecordMetrics, RecordOptions};

fn replace_string_tag(record: &mut bam::Record, tag: &[u8; 2], value: Option<&str>) -> Result<()> {
    let _ = record.remove_aux(tag);
    if let Some(value) = value {
        record.push_aux(tag, Aux::String(value))?;
    }
    Ok(())
}

pub(super) fn replace_i32_tag(
    record: &mut bam::Record,
    tag: &[u8; 2],
    value: Option<i32>,
) -> Result<()> {
    let _ = record.remove_aux(tag);
    if let Some(value) = value {
        record.push_aux(tag, Aux::I32(value))?;
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

/// [`annotate_record_with`] without the `mm` rescue (pre-0.6.2 behaviour).
pub fn annotate_record(
    record: bam::Record,
    annotation: &Annotation,
    contig_names: &[Vec<u8>],
    include_introns: bool,
    three_prime: bool,
) -> Result<(bam::Record, RecordMetrics)> {
    let opts = RecordOptions {
        include_introns,
        three_prime,
        confident_mapq: 0,
    };
    annotate_record_with(record, annotation, contig_names, &opts)
}

/// Annotate one record. With `opts.confident_mapq > 0`, a primary read below
/// it whose genes are one exonic-compatible gene gets `mm:i:1` (Cell
/// Ranger's `rescue_alignments_se` for a read with a single reported hit).
pub fn annotate_record_with(
    mut record: bam::Record,
    annotation: &Annotation,
    contig_names: &[Vec<u8>],
    opts: &RecordOptions,
) -> Result<(bam::Record, RecordMetrics)> {
    let RecordOptions {
        include_introns,
        three_prime,
        confident_mapq,
    } = *opts;
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
    metrics.exonic_compatible = call.exonic_compatible;
    if confident_mapq > 0 {
        let rescued = !record.is_secondary()
            && !record.is_supplementary()
            && record.mapq() < confident_mapq
            && call.exonic_compatible
            && call.genes.len() == 1;
        replace_i32_tag(&mut record, b"mm", rescued.then_some(1))?;
        metrics.rescued = rescued as u64;
    }
    Ok((record, metrics))
}

fn write_unordered_batch(
    batch: Vec<bam::Record>,
    pool: &rayon::ThreadPool,
    annotation: &Annotation,
    contig_names: &[Vec<u8>],
    opts: &RecordOptions,
    writer: &mut bam::Writer,
    metrics: &mut Metrics,
) -> Result<()> {
    let (sender, receiver) = crossbeam_channel::unbounded();
    pool.install(|| {
        batch
            .into_par_iter()
            .for_each_with(sender, |sender, record| {
                // The receiver is alive until every worker finishes this batch.
                let _ = sender.send(annotate_record_with(record, annotation, contig_names, opts));
            });
    });
    for result in receiver {
        let (record, record_metrics) = result?;
        writer.write(&record)?;
        observe(metrics, &record, &record_metrics);
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
    if cli.unordered && cli.multimap_tags {
        bail!("--multimap-tags is incompatible with --unordered");
    }
    let started = Instant::now();
    info!(gtf = %cli.gtf.display(), "loading annotation");
    let annotation = Annotation::load(&cli.gtf, cli.reference_filter, cli.strip_id_versions)?;
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

    let opts = RecordOptions {
        include_introns: cli.include_introns,
        three_prime: cli.three_prime,
        confident_mapq: cli.confident_mapq,
    };
    let mut metrics = Metrics {
        input_bam: cli.bam.display().to_string(),
        output_bam: cli.output.display().to_string(),
        gtf: cli.gtf.display().to_string(),
        genes_mtx: cli
            .genes_mtx
            .as_ref()
            .map(|path| path.display().to_string()),
        include_introns: cli.include_introns,
        reference_filter: cli.reference_filter,
        strip_id_versions: cli.strip_id_versions,
        three_prime: cli.three_prime,
        unordered: cli.unordered,
        multimap_tags: cli.multimap_tags,
        confident_mapq: cli.confident_mapq,
        genes_loaded: annotation.genes.len(),
        transcripts_loaded: annotation.transcripts.len(),
        gtf_exons_seen: annotation.exons_seen,
        gtf_exons_filtered_biotype: annotation.exons_filtered_biotype,
        gtf_exons_filtered_tag: annotation.exons_filtered_tag,
        per_barcode: crate::barcode_metrics::BarcodeTable::new(cli.barcode_metrics),
        ..Default::default()
    };
    if cli.multimap_tags {
        let mut group = Vec::new();
        let mut qname = Vec::new();
        for result in reader.records() {
            let record = result.context("reading BAM record")?;
            if !group.is_empty() && record.qname() != qname.as_slice() {
                write_multimap_group(
                    std::mem::take(&mut group),
                    &annotation,
                    &contig_names,
                    &opts,
                    &mut writer,
                    &mut metrics,
                )?;
            }
            if group.is_empty() {
                qname = record.qname().to_vec();
            }
            group.push(record);
        }
        if !group.is_empty() {
            write_multimap_group(
                group,
                &annotation,
                &contig_names,
                &opts,
                &mut writer,
                &mut metrics,
            )?;
        }
    } else if cli.unordered {
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
                    &opts,
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
                &opts,
                &mut writer,
                &mut metrics,
            )?;
        }
    } else {
        for result in reader.records() {
            let record = result.context("reading BAM record")?;
            let (record, record_metrics) =
                annotate_record_with(record, &annotation, &contig_names, &opts)?;
            writer.write(&record)?;
            observe(&mut metrics, &record, &record_metrics);
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
    metrics.per_barcode.write(
        &metrics_path,
        &[
            "records",
            "primary",
            "unmapped",
            "exonic",
            "intronic",
            "intergenic",
            "antisense",
            "gene_assigned",
            "unique_gene",
            "multi_gene",
            "rescued",
        ],
    )?;
    info!(records = metrics.records, elapsed = format!("{:.1}s", metrics.elapsed_secs), output = %cli.output.display(), "RUN SUMMARY");
    Ok(())
}
