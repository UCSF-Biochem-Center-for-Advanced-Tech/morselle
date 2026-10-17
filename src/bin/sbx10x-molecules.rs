use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use anyhow::{bail, Result};
use clap::{CommandFactory, FromArgMatches, Parser};
use sbx10x::molecules::{
    barcode_metrics, build_decisions_with, collect_evidence_indexed_with, collect_evidence_with,
    run_sharded, write_output, KeyTables, Metrics, MoleculeOptions, UmiRules,
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
    /// Output BAM. Input order is preserved unless --shards > 1, which
    /// emits unsorted hash buckets; sort that output before indexing.
    #[arg(short, long)]
    output: PathBuf,
    /// JSON metrics output (default: <output>.metrics.json).
    #[arg(long)]
    metrics: Option<PathBuf>,

    /// Also write per-barcode counters to `<metrics>.barcodes.tsv.gz` (MANUAL
    /// "Per-barcode read metrics"). Off by default; it adds run time and memory.
    #[arg(long)]
    barcode_metrics: bool,
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
    /// genuinely distinct genes. `--tagseq` counts every gene in a collision
    /// in full (`XW:S`) instead of adjudicating a winner, and skips 1-mismatch
    /// UMI correction (at pooled-sample UMI diversity, two distinct molecules
    /// 1bp apart are common) unless --position-aware-umis is also given.
    /// See also --fractional-collisions. Leaves single-cell runs unchanged.
    #[arg(long)]
    tagseq: bool,
    /// Split each cross-gene UMI collision across its genes in proportion to
    /// read support (weight = support / total support), instead of
    /// winner-take-all (single-cell default) or counting every claimant in
    /// full (--tagseq default). Weights for one collision sum to 1.0, so a
    /// collision counts as ONE molecule: right if collisions are mostly one
    /// molecule seen on two genes, an undercount if they are mostly distinct
    /// molecules sharing a UMI by chance (the usual case at tag-seq depth).
    /// Written as `ZW:f` with `XW:F`; matrix then writes a real-valued
    /// (non-integer) matrix. Off by default.
    #[arg(long)]
    fractional_collisions: bool,
    /// BGZF compression level (0-9) for the output BAM. Default: htslib's
    /// (6). Compression is a large share of a molecules run; on a 1M-read
    /// benchmark level 1 made the whole run ~24% faster for a ~12% larger
    /// file -- a good trade for an intermediate read once more by
    /// `matrix`/`audit`. Records are identical at every level.
    #[arg(long, value_parser = clap::value_parser!(u32).range(0..=9))]
    compression_level: Option<u32>,
    /// Accepted for compatibility with 0.5.6-0.5.9, where fractional was the
    /// --tagseq default. Now a no-op: fractional is opt-in.
    #[arg(long, conflicts_with = "fractional_collisions", hide = true)]
    no_fractional_collisions: bool,
    /// Treat reads as one molecule only if they also share their alignment
    /// start: (cell, gene, UMI, strand-aware unclipped 5' position) instead
    /// of (cell, gene, UMI). PCR duplicates share a start; distinct molecules
    /// that drew the same UMI almost never do, so chance UMI matches stay
    /// separate -- across genes and, unlike every other mode, within one
    /// gene. 1-mismatch UMI correction then runs among reads at the same
    /// position, even under --tagseq. The usual dedup key for bulk 3'
    /// tag-seq. Requires UMIs of at most 16 bases. Caveat: an indel or
    /// soft-clip difference at the read start splits duplicates (slight
    /// overcount); compare against an orthogonal count before relying on it.
    #[arg(long)]
    position_aware_umis: bool,
    /// 10x Flex mode. Require a non-empty probe/sample barcode (`PB`) on
    /// every counted molecule so identical GEM barcodes from distinct Flex
    /// samples cannot be merged.
    #[arg(long, conflicts_with = "tagseq")]
    flex: bool,
    /// Count only confidently mapped reads, as Cell Ranger does: MAPQ at
    /// least this (60 = minimap2's unique; STAR's 255 passes), or rescued by
    /// `sbxpress annotate` (`mm:i:1`: one compatible exonic gene). Dropped
    /// reads are reported as `low_mapq_records` / `multimapper_records`
    /// (`NH` > 1). 0 counts every MAPQ (the pre-0.6.2 behaviour).
    #[arg(long, default_value_t = 60)]
    min_mapq: u8,
    /// Also merge UMIs one insertion or deletion apart (same-length window,
    /// so the bases after the indel are shifted), under the same 2n-1
    /// read-support rule as 1-mismatch correction. SBX errors are mostly
    /// 1-base indels, which 1-mismatch correction cannot merge; on CAT06
    /// SBX CRISPR this cut molecules from 6.1x Cell Ranger's to 1.9x. Off by
    /// default; no effect where correction is skipped (--tagseq without
    /// --position-aware-umis). GEX impact not yet measured.
    #[arg(long)]
    umi_indel_correction: bool,
    /// After UMI correction, fold every 1-read molecule into its (cell,
    /// sample, gene) family's top molecule when that top has at least this
    /// many reads. For deeply sequenced feature-barcoding libraries, where a
    /// low-quality read with several UMI errors otherwise counts as a new
    /// molecule (CAT06 SBX CRISPR: with --umi-indel-correction and 10,
    /// molecules 0.99x Cell Ranger's). 0 (default) disables it. Not for GEX,
    /// where 1-read molecules are mostly real.
    #[arg(long, default_value_t = 0)]
    fold_singletons_min_top: u32,
}

fn write_metrics(cli: &Cli, metrics: &mut Metrics) -> Result<()> {
    metrics.artifacts = metrics.gene_tally.report();
    metrics.artifacts.warn();
    let metrics_path = cli
        .metrics
        .clone()
        .unwrap_or_else(|| PathBuf::from(format!("{}.metrics.json", cli.output.display())));
    fs::write(&metrics_path, serde_json::to_string_pretty(&*metrics)?)?;
    barcode_metrics::table(&metrics.barcode_rows, cli.barcode_metrics)
        .write(&metrics_path, &barcode_metrics::COLUMNS)?;
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
        fractional: cli.fractional_collisions,
        position_aware: cli.position_aware_umis,
        min_mapq: cli.min_mapq,
        compression_level: cli.compression_level,
        umi_rules: UmiRules {
            indels: cli.umi_indel_correction,
            fold_singletons_min_top: cli.fold_singletons_min_top,
        },
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
        fractional_collisions: cli.fractional_collisions,
        position_aware_umis: cli.position_aware_umis,
        umi_indel_correction: cli.umi_indel_correction,
        fold_singletons_min_top: cli.fold_singletons_min_top,
        min_mapq: cli.min_mapq,
        ..Default::default()
    };
    if cli.shards > 1 {
        if cli.saturation.is_some() {
            bail!("--saturation is not yet supported together with --shards > 1");
        }
        info!(
            shards = cli.shards,
            fractional = opts.fractional,
            "running sharded (collate + quant)"
        );
        run_sharded(
            &cli.bam,
            &cli.output,
            cli.shards,
            cli.shard_concurrency,
            &opts,
            &mut metrics,
        )?;
        metrics.elapsed_secs = started.elapsed().as_secs_f64();
        return write_metrics(&cli, &mut metrics);
    }
    // Validated above: cli.shards <= 1 requires exactly one input file.
    let single_bam = &cli.bam[0];
    let (evidence, genes, gene_ids, _samples, sample_ids) = if cli.evidence_workers > 1 {
        info!(
            workers = cli.evidence_workers,
            "collecting evidence from indexed BAM regions"
        );
        collect_evidence_indexed_with(
            single_bam,
            cli.evidence_workers,
            cli.evidence_chunk_bases,
            &opts.record_filter(),
            &mut metrics,
        )?
    } else {
        collect_evidence_with(single_bam, cli.threads, &opts.record_filter(), &mut metrics)?
    };
    eprintln!(
        "[sbxpress molecules] phase=evidence_ready raw_nodes={} eligible={}",
        evidence.len(),
        metrics.eligible_records
    );
    let (decisions, saturation) = build_decisions_with(
        &evidence,
        &mut metrics,
        cli.saturation.as_ref().map(|_| cli.saturation_interval),
        cli.evidence_workers,
        opts.dedup_policy(),
    );
    drop(evidence);
    if let (Some(curve), Some(path)) = (saturation, cli.saturation.as_deref()) {
        eprintln!(
            "[sbxpress molecules] phase=saturation_write start output={}",
            path.display()
        );
        metrics.saturation = Some(curve.write(path)?);
        eprintln!("[sbxpress molecules] phase=saturation_write complete");
    }
    let tables = KeyTables {
        genes: &genes,
        gene_ids: &gene_ids,
        sample_ids: &sample_ids,
    };
    write_output(
        single_bam,
        &cli.output,
        &opts,
        &decisions,
        &tables,
        &mut metrics,
    )?;
    metrics.elapsed_secs = started.elapsed().as_secs_f64();
    write_metrics(&cli, &mut metrics)
}
