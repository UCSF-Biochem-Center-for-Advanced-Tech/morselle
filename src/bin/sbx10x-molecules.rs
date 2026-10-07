use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{bail, Result};
use clap::{CommandFactory, FromArgMatches, Parser};
use sbx10x::molecules::{
    build_decisions, collect_evidence, collect_evidence_indexed, run_sharded, write_output,
    Metrics,
};
use tracing::info;

#[derive(Debug, Parser)]
#[command(
    name = "sbx10x-molecules",
    version,
    about = "Correct gene-aware UMIs, resolve gene collisions, and mark molecule representatives"
)]
pub struct Cli {
    /// Annotated BAM(s) containing CB, UB/UR, GX, and optionally ZG tags.
    #[arg(value_name = "BAM", num_args = 1..)]
    pub bam: Vec<PathBuf>,
    /// Output BAM. Input order is preserved.
    #[arg(short, long)]
    pub output: PathBuf,
    /// JSON metrics output (default: <output>.metrics.json).
    #[arg(long)]
    pub metrics: Option<PathBuf>,
    /// BAM compression/decompression threads.
    #[arg(long, default_value_t = 4)]
    pub threads: usize,
    /// Number of indexed BAM-region workers for the evidence pass.
    #[arg(long, default_value_t = 1)]
    pub evidence_workers: usize,
    /// Genomic width of a work item for parallel indexed evidence collection.
    #[arg(long, default_value_t = 10_000_000, hide = true)]
    pub evidence_chunk_bases: u64,
    /// Write an order-independent sequencing-saturation curve as a gzip TSV.
    #[arg(long)]
    pub saturation: Option<PathBuf>,
    /// Eligible-read interval between saturation observations.
    #[arg(long, default_value_t = 1_000_000)]
    pub saturation_interval: u64,
    /// Bound peak memory by partitioning evidence into N shards.
    #[arg(long, default_value_t = 1)]
    pub shards: usize,
    /// Number of shards to process concurrently when --shards > 1.
    #[arg(long, default_value_t = 1)]
    pub shard_concurrency: usize,
    /// Pooled-sample mode: a barcode is a whole sample.
    #[arg(long)]
    pub tagseq: bool,
    /// 10x Flex mode. Require a non-empty probe/sample barcode (`PB`).
    #[arg(long, conflicts_with = "tagseq")]
    pub flex: bool,
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
        run_sharded(
            &cli.bam,
            &cli.output,
            cli.threads,
            cli.shards,
            cli.shard_concurrency,
            cli.evidence_workers,
            cli.tagseq,
            cli.flex,
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
        return Ok(());
    }
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
        &cli.bam,
        &cli.output,
        cli.threads,
        cli.evidence_workers,
        cli.flex,
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
