use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{bail, Result};
use clap::{CommandFactory, FromArgMatches, Parser};
use sbx10x::molecules::{
    build_decisions, collect_evidence, collect_evidence_indexed, run_sharded, write_output,
    KeyTables, Metrics, MoleculeOptions,
};
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
    /// UMI nodes instead of adjudicating a winner. It also skips 1-mismatch
    /// UMI correction (at pooled-sample UMI diversity, two distinct molecules
    /// 1bp apart are common). Collisions are split fractionally by default
    /// under --tagseq (see --fractional-collisions); pass
    /// --no-fractional-collisions to count every claimant in full
    /// (`XW:S`). Leaves single-cell runs (the default) unchanged.
    #[arg(long)]
    tagseq: bool,
    /// Split each cross-gene UMI collision across its genes in proportion to
    /// read support (weight = support / total support), instead of
    /// winner-take-all (single-cell default) or counting every claimant in
    /// full (--tagseq --no-fractional-collisions). Weights for one collision
    /// sum to 1.0 and are written as `ZW:f` with `XW:F`; cr-bam-to-mtx then
    /// writes a real-valued matrix. Default: on with --tagseq, off otherwise.
    #[arg(long)]
    fractional_collisions: bool,
    /// Disable fractional collision resolution (overrides the --tagseq default).
    #[arg(long, conflicts_with = "fractional_collisions")]
    no_fractional_collisions: bool,
    /// 10x Flex mode. Require a non-empty probe/sample barcode (`PB`) on
    /// every counted molecule so identical GEM barcodes from distinct Flex
    /// samples cannot be merged.
    #[arg(long, conflicts_with = "tagseq")]
    flex: bool,
}

fn write_metrics(cli: &Cli, metrics: &Metrics) -> Result<()> {
    let metrics_path = cli
        .metrics
        .clone()
        .unwrap_or_else(|| PathBuf::from(format!("{}.metrics.json", cli.output.display())));
    fs::write(&metrics_path, serde_json::to_string_pretty(metrics)?)?;
    info!(
        records = metrics.records,
        molecules = metrics.molecules_counted,
        collisions = metrics.collision_groups,
        elapsed = format!("{:.1}s", metrics.elapsed_secs),
        "RUN SUMMARY"
    );
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
    let opts = MoleculeOptions {
        threads: cli.threads,
        evidence_workers: cli.evidence_workers,
        flex: cli.flex,
        tagseq: cli.tagseq,
        fractional: (cli.tagseq || cli.fractional_collisions) && !cli.no_fractional_collisions,
    };
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
    if cli.shards > 1 {
        if cli.saturation.is_some() {
            bail!("--saturation is not yet supported together with --shards > 1");
        }
        info!(shards = cli.shards, fractional = opts.fractional, "running sharded (collate + quant)");
        run_sharded(&cli.bam, &cli.output, cli.shards, cli.shard_concurrency, &opts, &mut metrics)?;
        metrics.elapsed_secs = started.elapsed().as_secs_f64();
        return write_metrics(&cli, &metrics);
    }
    // Validated above: cli.shards <= 1 requires exactly one input file.
    let single_bam = &cli.bam[0];
    let (evidence, genes, gene_ids, _samples, sample_ids) = if cli.evidence_workers > 1 {
        info!(workers = cli.evidence_workers, "collecting evidence from indexed BAM regions");
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
        opts.tagseq,
        opts.fractional,
    );
    drop(evidence);
    if let (Some(curve), Some(path)) = (saturation, cli.saturation.as_deref()) {
        eprintln!("[sbxpress molecules] phase=saturation_write start output={}", path.display());
        metrics.saturation = Some(curve.write(path)?);
        eprintln!("[sbxpress molecules] phase=saturation_write complete");
    }
    let tables = KeyTables { genes: &genes, gene_ids: &gene_ids, sample_ids: &sample_ids };
    write_output(single_bam, &cli.output, &opts, &decisions, &tables, &mut metrics)?;
    metrics.elapsed_secs = started.elapsed().as_secs_f64();
    write_metrics(&cli, &metrics)
}
