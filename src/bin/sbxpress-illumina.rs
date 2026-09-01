//! Convert standard paired Illumina 10x GEX FASTQs into the tagged R2-only
//! FASTQ contract consumed by `minimap2 -y` in the sbxpress workflow.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::write::GzEncoder;
use flate2::Compression;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use sbx10x::anchors::revcomp;
use sbx10x::barcode::{encode16, Whitelist};
use sbx10x::extract::{extract, ExtractCfg};
use sbx10x::overlap::{overlap_in, InsertClass, Merged, OverlapCfg, Scratch};
use sbx10x::yc::YcTag;
use serde::Serialize;

const CB_LEN: usize = 16;
const UMI_LEN: usize = 12;

#[derive(Debug, Parser)]
#[command(about = "Convert paired Illumina 10x FASTQs to minimap2-tagged R2 FASTQ")]
struct Cli {
    /// Illumina Read 1 FASTQ. Repeat once per lane/shard, matching --r2 order.
    #[arg(long, required = true, action = clap::ArgAction::Append)]
    r1: Vec<PathBuf>,

    /// Illumina Read 2 FASTQ. Repeat once per lane/shard, matching --r1 order.
    #[arg(long, required = true, action = clap::ArgAction::Append)]
    r2: Vec<PathBuf>,

    /// 10x cell-barcode whitelist (plain text or .gz).
    #[arg(long)]
    whitelist: PathBuf,

    /// Output tagged R2 FASTQ (.gz recommended).
    #[arg(short, long)]
    output: PathBuf,

    /// JSON metrics path. Defaults to <output>.metrics.json.
    #[arg(long)]
    metrics: Option<PathBuf>,

    /// GEM-group suffix placed on corrected CB tags.
    #[arg(long, default_value_t = 1)]
    gem_group: u32,

    /// Maximum number of R1/R2 file pairs converted concurrently.
    #[arg(short, long, default_value_t = 0)]
    threads: usize,

    /// Merge overlapping mates into the full fragment before extraction. Pairs
    /// that merge are processed exactly like an SBX full-length fragment; pairs
    /// that do not fall back to the standard R1-barcode / R2-cDNA path.
    #[arg(long, default_value_t = false)]
    overlap: bool,

    /// Shortest mate overlap that may be merged.
    #[arg(long, default_value_t = 20)]
    overlap_min_overlap: usize,

    /// Seed length for the overlap lookup table.
    #[arg(long, default_value_t = 8)]
    overlap_kmer: usize,

    /// Stride between indexed overlap seeds.
    #[arg(long, default_value_t = 1)]
    overlap_kmer_offset: usize,

    /// Usable seeds probed inward from each end of the shorter mate.
    #[arg(long, default_value_t = 20)]
    overlap_check_lengths: usize,

    /// Absolute cap on mismatches within an overlap.
    #[arg(long, default_value_t = 100)]
    overlap_max_mismatch: usize,

    /// Cap on mismatches as a fraction of the overlap length.
    #[arg(long, default_value_t = 0.25)]
    overlap_max_mismatch_density: f64,

    /// Omit the YC tag from merged reads.
    #[arg(long, default_value_t = false)]
    no_emit_yc: bool,

    /// Minimum cDNA length to emit from a merged fragment.
    #[arg(long, default_value_t = 50)]
    min_cdna_len: usize,

    /// Max mismatches for the 3' adapter seed when extracting a merged fragment.
    #[arg(long, default_value_t = 2)]
    adapter_max_mm: u32,

    /// Emit every read that fails, as the original R1/R2 pair, to
    /// `<prefix>_R1.fastq.gz` and `<prefix>_R2.fastq.gz`. Each record keeps its
    /// original header and gains an `XF:Z:<reason>` tag. With this set,
    /// `records_read == records_emitted + records_failed` exactly.
    #[arg(long)]
    failed_output: Option<PathBuf>,

    /// Trim the template-switch oligo from merged fragments. The TSO survives
    /// only on fragments that reached the cDNA's 5' end — the same ones short
    /// enough to overlap — so merged reads are enriched for it. Off by default.
    #[arg(long, default_value_t = false)]
    trim_tso: bool,
}

#[derive(Debug, Default, Serialize)]
struct Metrics {
    r1_files: Vec<PathBuf>,
    r2_files: Vec<PathBuf>,
    whitelist: PathBuf,
    whitelist_size: usize,
    output: PathBuf,
    output_shards: Vec<PathBuf>,
    gem_group: u32,
    records_read: u64,
    records_emitted: u64,
    barcode_exact: u64,
    barcode_hamming1_corrected: u64,
    barcode_rejected: u64,
    short_r1: u64,
    overlap_enabled: bool,
    /// Merged, fragment shorter than a read: read-through was present on a mate.
    overlap_short_sins: u64,
    /// Merged, neither mate ran off the end of the molecule.
    overlap_medium_mins: u64,
    /// Not merged: the mates did not reach each other.
    overlap_long_lins: u64,
    /// Merged, but no cDNA survived extraction (below `--min-cdna-len`).
    merged_without_cdna: u64,
    merged_bases: u64,
    discordant_bases: u64,
    yc_n_columns_coerced: u64,
    yc_unencodable: u64,
    /// Merged fragments whose TSO was found and trimmed (needs `--trim-tso`).
    tso_trimmed: u64,
    /// Reads written to the failed-output side channel, by reason.
    records_failed: u64,
    failed_by_reason: BTreeMap<String, u64>,
    failed_outputs: Vec<PathBuf>,
    /// cDNA length of each merged fragment, excluding CB+UMI+dT.
    insert_length_histogram: BTreeMap<usize, u64>,
    elapsed_secs: f64,
}

#[derive(Debug, Default)]
struct PairMetrics {
    records_read: u64,
    records_emitted: u64,
    barcode_exact: u64,
    barcode_hamming1_corrected: u64,
    barcode_rejected: u64,
    short_r1: u64,
    overlap_short: u64,
    overlap_medium: u64,
    overlap_long: u64,
    merged_without_cdna: u64,
    merged_bases: u64,
    discordant_bases: u64,
    yc_n_columns_coerced: u64,
    yc_unencodable: u64,
    tso_trimmed: u64,
    records_failed: u64,
    failed_by_reason: BTreeMap<String, u64>,
    insert_length_histogram: BTreeMap<usize, u64>,
}

fn canonical_id(id: &[u8]) -> &[u8] {
    let id = id.split(|b| b.is_ascii_whitespace()).next().unwrap_or(id);
    if id.ends_with(b"/1") || id.ends_with(b"/2") {
        &id[..id.len() - 2]
    } else {
        id
    }
}

fn metrics_path(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_os_string();
    name.push(".metrics.json");
    PathBuf::from(name)
}

fn shard_path(output: &Path, index: usize, count: usize) -> PathBuf {
    if count == 1 {
        return output.to_path_buf();
    }
    let text = output.as_os_str().to_string_lossy();
    if let Some(prefix) = text.strip_suffix(".fastq.gz") {
        PathBuf::from(format!("{prefix}.shard-{:03}.fastq.gz", index + 1))
    } else if let Some(prefix) = text.strip_suffix(".fq.gz") {
        PathBuf::from(format!("{prefix}.shard-{:03}.fq.gz", index + 1))
    } else {
        PathBuf::from(format!("{text}.shard-{:03}", index + 1))
    }
}

/// Append one tagged FASTQ record. `seq`/`qual` are already in R2 (sense)
/// orientation; `yc` describes the merge in that same orientation.
#[allow(clippy::too_many_arguments)]
fn write_record(
    output: &mut dyn Write,
    id: &[u8],
    raw_cb: &[u8],
    raw_umi: &[u8],
    cb_qual: &[u8],
    umi_qual: &[u8],
    corrected: &[u8],
    gem_group: u32,
    seq: &[u8],
    qual: &[u8],
    yc: Option<&YcTag>,
) -> std::io::Result<()> {
    // Illumina FASTQ identifiers commonly include a whitespace-delimited
    // mate/control/index comment (for example `2:N:0:...`). Minimap2 -y treats
    // everything after the first whitespace as SAM-style auxiliary fields, so
    // emit only the canonical cluster name before our tags.
    output.write_all(b"@")?;
    output.write_all(canonical_id(id))?;
    output.write_all(b"\tCR:Z:")?;
    output.write_all(raw_cb)?;
    output.write_all(b"\tCB:Z:")?;
    output.write_all(corrected)?;
    write!(output, "-{gem_group}")?;
    output.write_all(b"\tUR:Z:")?;
    output.write_all(raw_umi)?;
    output.write_all(b"\tUB:Z:")?;
    output.write_all(raw_umi)?;
    output.write_all(b"\tCY:Z:")?;
    output.write_all(cb_qual)?;
    output.write_all(b"\tUY:Z:")?;
    output.write_all(umi_qual)?;
    if let Some(tag) = yc {
        write!(output, "\tYC:Z:{tag}")?;
    }
    output.write_all(b"\n")?;
    output.write_all(seq)?;
    output.write_all(b"\n+\n")?;
    output.write_all(qual)?;
    output.write_all(b"\n")
}

#[allow(clippy::too_many_arguments)]
/// Paths for a failed-pair shard: `<prefix>_R1.fastq.gz` / `<prefix>_R2.fastq.gz`,
/// carrying the same `.shard-NNN` infix as the main output when sharded.
fn failed_paths(prefix: &Path, index: usize, count: usize) -> (PathBuf, PathBuf) {
    let text = prefix.as_os_str().to_string_lossy().into_owned();
    let stem = if count == 1 {
        text
    } else {
        format!("{text}.shard-{:03}", index + 1)
    };
    (
        PathBuf::from(format!("{stem}_R1.fastq.gz")),
        PathBuf::from(format!("{stem}_R2.fastq.gz")),
    )
}

/// Both mates of a failed pair, written verbatim apart from an appended
/// `XF:Z:<reason>` tag. Emitting the untouched reads keeps the side channel
/// re-runnable through any other tool.
struct FailedWriter {
    r1: Box<dyn Write>,
    r2: Box<dyn Write>,
}

impl FailedWriter {
    fn create(r1_path: &Path, r2_path: &Path) -> Result<Self> {
        let open = |path: &Path| -> Result<Box<dyn Write>> {
            let file =
                File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
            Ok(Box::new(GzEncoder::new(
                BufWriter::new(file),
                Compression::default(),
            )))
        };
        Ok(FailedWriter {
            r1: open(r1_path)?,
            r2: open(r2_path)?,
        })
    }

    fn write_pair(
        &mut self,
        reason: &str,
        r1: (&[u8], &[u8], &[u8]),
        r2: (&[u8], &[u8], &[u8]),
    ) -> std::io::Result<()> {
        for (out, (id, seq, qual)) in [(&mut self.r1, r1), (&mut self.r2, r2)] {
            out.write_all(b"@")?;
            out.write_all(id)?;
            write!(out, "\tXF:Z:{reason}")?;
            out.write_all(b"\n")?;
            out.write_all(seq)?;
            out.write_all(b"\n+\n")?;
            out.write_all(qual)?;
            out.write_all(b"\n")?;
        }
        Ok(())
    }
}

#[allow(clippy::too_many_arguments)]
fn convert_pair(
    r1_path: &Path,
    r2_path: &Path,
    output_path: &Path,
    whitelist: &Whitelist,
    gem_group: u32,
    overlap_cfg: Option<&OverlapCfg>,
    extract_cfg: &ExtractCfg,
    emit_yc: bool,
    failed_paths: Option<(&Path, &Path)>,
) -> Result<PairMetrics> {
    let file = File::create(output_path)
        .with_context(|| format!("cannot create {}", output_path.display()))?;
    let gz = output_path.extension().is_some_and(|value| value == "gz");
    let mut output: Box<dyn Write> = if gz {
        Box::new(GzEncoder::new(BufWriter::new(file), Compression::default()))
    } else {
        Box::new(BufWriter::new(file))
    };
    let mut metrics = PairMetrics::default();
    let mut r1 = parse_fastx_file(r1_path)
        .with_context(|| format!("cannot open R1 {}", r1_path.display()))?;
    let mut r2 = parse_fastx_file(r2_path)
        .with_context(|| format!("cannot open R2 {}", r2_path.display()))?;
    let mut failed = failed_paths
        .map(|(r1, r2)| FailedWriter::create(r1, r2))
        .transpose()?;
    let mut pair_number = 0u64;
    let mut scratch = Scratch::default();
    loop {
        let a = r1.next();
        let b = r2.next();
        match (a, b) {
            (None, None) => break,
            (None, Some(_)) | (Some(_), None) => bail!(
                "R1/R2 record-count mismatch between {} and {} after {} pairs",
                r1_path.display(),
                r2_path.display(),
                pair_number
            ),
            (Some(a), Some(b)) => {
                let a = a.with_context(|| format!("invalid FASTQ in {}", r1_path.display()))?;
                let b = b.with_context(|| format!("invalid FASTQ in {}", r2_path.display()))?;
                pair_number += 1;
                metrics.records_read += 1;
                if canonical_id(a.id()) != canonical_id(b.id()) {
                    bail!(
                        "R1/R2 name mismatch at pair {}: '{}' versus '{}'",
                        pair_number,
                        String::from_utf8_lossy(a.id()),
                        String::from_utf8_lossy(b.id())
                    );
                }
                let r1_seq = a.seq();
                let r2_seq = b.seq();
                let Some(r1_qual) = a.qual() else {
                    bail!("R1 record {} has no qualities", pair_number)
                };
                let Some(r2_qual) = b.qual() else {
                    bail!("R2 record {} has no qualities", pair_number)
                };

                /// Record a failure, emit the untouched pair to the side channel,
                /// and move on. Every `continue` below goes through this, so
                /// `records_read == records_emitted + records_failed` holds.
                macro_rules! fail_pair {
                    ($reason:expr) => {{
                        let reason: &str = $reason;
                        metrics.records_failed += 1;
                        *metrics
                            .failed_by_reason
                            .entry(reason.to_string())
                            .or_insert(0) += 1;
                        if let Some(w) = failed.as_mut() {
                            w.write_pair(
                                reason,
                                (a.id(), &r1_seq[..], r1_qual),
                                (b.id(), &r2_seq[..], r2_qual),
                            )?;
                        }
                        continue;
                    }};
                }

                // Step 1: try to reconstruct the whole fragment from the RAW
                // mates. This must happen before anything is trimmed — the RT
                // primer's dT block is the polyA complement, so R1 and rc(R2)
                // describe the same span only while both are intact.
                let merged: Option<Merged> = overlap_cfg.and_then(|cfg| {
                    overlap_in(&mut scratch, &r1_seq, r1_qual, &r2_seq, r2_qual, cfg)
                });

                // Buffers for the merged path, kept alive across the branch.
                let mut merged_seq: Vec<u8> = Vec::new();
                let mut merged_qual: Vec<u8> = Vec::new();
                let mut merged_tag: Option<YcTag> = None;

                if let Some(m) = &merged {
                    // (a) sins / mins. The merged molecule is
                    // `CB | UMI | dT | rc(cDNA)` — the same canonical layout the
                    // SBX arm produces — so it goes through the *same* extractor
                    // rather than bespoke logic. Running identical code on both
                    // arms is what keeps the comparison honest.
                    match m.class {
                        InsertClass::Short => metrics.overlap_short += 1,
                        InsertClass::Medium => metrics.overlap_medium += 1,
                        InsertClass::Long => {}
                    }
                    metrics.merged_bases += m.seq.len() as u64;
                    metrics.discordant_bases += m.discordant as u64;
                    metrics.yc_n_columns_coerced += m.yc_stats.n_columns_coerced;
                    if emit_yc && m.yc.is_none() {
                        metrics.yc_unencodable += 1;
                    }
                    if m.seq.len() < CB_LEN + UMI_LEN {
                        metrics.short_r1 += 1;
                        metrics.barcode_rejected += 1;
                        fail_pair!("merged_too_short");
                    }

                    let e = extract(
                        &m.seq,
                        Some((0, CB_LEN + UMI_LEN)),
                        Some(CB_LEN + UMI_LEN),
                        extract_cfg,
                    );
                    let (cdna_start, cdna_end) = match e.r2 {
                        Some(range) => range,
                        None => {
                            metrics.merged_without_cdna += 1;
                            fail_pair!("no_cdna");
                        }
                    };
                    if e.tso_trimmed {
                        metrics.tso_trimmed += 1;
                    }
                    *metrics
                        .insert_length_histogram
                        .entry(cdna_end - cdna_start)
                        .or_insert(0) += 1;

                    // Emit the cDNA reverse-complemented into R2 orientation, so
                    // every downstream stage sees the convention it already
                    // expects. The tag is trimmed to the emitted window and
                    // flipped with it — reversing the string would be wrong.
                    merged_seq = revcomp(&m.seq[cdna_start..cdna_end]);
                    merged_qual = m.qual[cdna_start..cdna_end].iter().rev().copied().collect();
                    if emit_yc {
                        merged_tag = m.yc.clone().map(|mut tag| {
                            tag.trim(cdna_start as u64, (m.seq.len() - cdna_end) as u64);
                            tag.reverse_complement();
                            tag
                        });
                        debug_assert!(merged_tag
                            .as_ref()
                            .is_none_or(|t| t.sequence_length() as usize == merged_seq.len()));
                    }
                } else {
                    // (b) lins, or overlap disabled: the standard Illumina path.
                    if overlap_cfg.is_some() {
                        metrics.overlap_long += 1;
                    }
                    if r1_seq.len() < CB_LEN + UMI_LEN || r1_qual.len() < CB_LEN + UMI_LEN {
                        metrics.short_r1 += 1;
                        metrics.barcode_rejected += 1;
                        fail_pair!("short_r1");
                    }
                }

                // The barcode comes from the merged molecule when we have one —
                // on a short fragment R2 read back through it, so CB and UMI
                // carry summed duplex quality.
                let (raw_cb, raw_umi, cb_qual, umi_qual, out_seq, out_qual) = match &merged {
                    Some(m) => (
                        &m.seq[..CB_LEN],
                        &m.seq[CB_LEN..CB_LEN + UMI_LEN],
                        &m.qual[..CB_LEN],
                        &m.qual[CB_LEN..CB_LEN + UMI_LEN],
                        &merged_seq[..],
                        &merged_qual[..],
                    ),
                    None => (
                        &r1_seq[..CB_LEN],
                        &r1_seq[CB_LEN..CB_LEN + UMI_LEN],
                        &r1_qual[..CB_LEN],
                        &r1_qual[CB_LEN..CB_LEN + UMI_LEN],
                        &r2_seq[..],
                        r2_qual,
                    ),
                };

                let cb_invalid = encode16(raw_cb).is_none();
                let corrected = if cb_invalid {
                    None
                } else if whitelist.contains(raw_cb) {
                    metrics.barcode_exact += 1;
                    whitelist.correct_within1(raw_cb)
                } else {
                    let value = whitelist.correct_within1(raw_cb);
                    if value.is_some() {
                        metrics.barcode_hamming1_corrected += 1;
                    }
                    value
                };
                let corrected = match corrected {
                    Some(value) => value,
                    None => {
                        metrics.barcode_rejected += 1;
                        fail_pair!(if cb_invalid {
                            "cb_invalid"
                        } else {
                            "cb_unmatched"
                        });
                    }
                };

                write_record(
                    &mut output,
                    b.id(),
                    raw_cb,
                    raw_umi,
                    cb_qual,
                    umi_qual,
                    &corrected,
                    gem_group,
                    out_seq,
                    out_qual,
                    merged_tag.as_ref(),
                )?;
                metrics.records_emitted += 1;
            }
        }
    }
    output.flush()?;
    drop(output);
    Ok(metrics)
}

fn run(cli: Cli) -> Result<()> {
    if cli.r1.len() != cli.r2.len() {
        bail!("--r1 and --r2 must be supplied the same number of times");
    }
    if cli.gem_group == 0 {
        bail!("--gem-group must be greater than zero");
    }
    let started = Instant::now();
    let overlap_cfg = cli.overlap.then_some(OverlapCfg {
        kmer: cli.overlap_kmer,
        kmer_offset: cli.overlap_kmer_offset,
        check_lengths: cli.overlap_check_lengths,
        min_overlap: cli.overlap_min_overlap,
        max_mismatch: cli.overlap_max_mismatch,
        max_mismatch_density: cli.overlap_max_mismatch_density,
        qual_cap: 40,
    });
    let extract_cfg = ExtractCfg {
        min_cdna_len: cli.min_cdna_len,
        adapter_max_mm: cli.adapter_max_mm,
        trim_tso: cli.trim_tso,
        max_r2_len: 0,
    };
    let emit_yc = !cli.no_emit_yc;
    let cli_trim_tso = cli.trim_tso;
    let whitelist = Whitelist::load(&cli.whitelist)?;
    if whitelist.is_empty() {
        bail!("whitelist is empty: {}", cli.whitelist.display());
    }
    let output_shards = (0..cli.r1.len())
        .map(|index| shard_path(&cli.output, index, cli.r1.len()))
        .collect::<Vec<_>>();
    let failed_shards: Vec<(PathBuf, PathBuf)> = cli
        .failed_output
        .as_ref()
        .map(|prefix| {
            (0..cli.r1.len())
                .map(|index| failed_paths(prefix, index, cli.r1.len()))
                .collect()
        })
        .unwrap_or_default();
    let mut metrics = Metrics {
        r1_files: cli.r1.clone(),
        r2_files: cli.r2.clone(),
        whitelist: cli.whitelist.clone(),
        whitelist_size: whitelist.len(),
        output: cli.output.clone(),
        output_shards: output_shards.clone(),
        gem_group: cli.gem_group,
        overlap_enabled: cli.overlap,
        failed_outputs: failed_shards
            .iter()
            .flat_map(|(r1, r2)| [r1.clone(), r2.clone()])
            .collect(),
        ..Metrics::default()
    };
    let workers = if cli.threads == 0 {
        cli.r1.len()
    } else {
        cli.threads.min(cli.r1.len())
    }
    .max(1);
    eprintln!(
        "[sbxpress illumina] converting {} pair(s) with {} worker(s)",
        cli.r1.len(),
        workers
    );
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(workers)
        .build()
        .context("cannot create Illumina conversion thread pool")?;
    let pair_metrics = pool.install(|| {
        cli.r1
            .par_iter()
            .zip(cli.r2.par_iter())
            .zip(output_shards.par_iter())
            .enumerate()
            .map(|(index, ((r1, r2), output))| {
                let result = convert_pair(
                    r1,
                    r2,
                    output,
                    &whitelist,
                    cli.gem_group,
                    overlap_cfg.as_ref(),
                    &extract_cfg,
                    emit_yc,
                    failed_shards
                        .get(index)
                        .map(|(a, b)| (a.as_path(), b.as_path())),
                );
                if let Ok(value) = &result {
                    eprintln!(
                        "[sbxpress illumina] shard {}/{} complete: read={} emitted={} \
                         sins={} mins={} lins={} output={}",
                        index + 1,
                        cli.r1.len(),
                        value.records_read,
                        value.records_emitted,
                        value.overlap_short,
                        value.overlap_medium,
                        value.overlap_long,
                        output.display()
                    );
                }
                result
            })
            .collect::<Result<Vec<_>>>()
    })?;
    for value in pair_metrics {
        metrics.records_read += value.records_read;
        metrics.records_emitted += value.records_emitted;
        metrics.barcode_exact += value.barcode_exact;
        metrics.barcode_hamming1_corrected += value.barcode_hamming1_corrected;
        metrics.barcode_rejected += value.barcode_rejected;
        metrics.short_r1 += value.short_r1;
        metrics.overlap_short_sins += value.overlap_short;
        metrics.overlap_medium_mins += value.overlap_medium;
        metrics.overlap_long_lins += value.overlap_long;
        metrics.merged_without_cdna += value.merged_without_cdna;
        metrics.merged_bases += value.merged_bases;
        metrics.discordant_bases += value.discordant_bases;
        metrics.yc_n_columns_coerced += value.yc_n_columns_coerced;
        metrics.yc_unencodable += value.yc_unencodable;
        metrics.tso_trimmed += value.tso_trimmed;
        metrics.records_failed += value.records_failed;
        for (reason, count) in value.failed_by_reason {
            *metrics.failed_by_reason.entry(reason).or_insert(0) += count;
        }
        for (length, count) in value.insert_length_histogram {
            *metrics.insert_length_histogram.entry(length).or_insert(0) += count;
        }
    }
    metrics.elapsed_secs = started.elapsed().as_secs_f64();
    let metrics_path = cli.metrics.unwrap_or_else(|| metrics_path(&cli.output));
    let metrics_file = File::create(&metrics_path)
        .with_context(|| format!("cannot create {}", metrics_path.display()))?;
    serde_json::to_writer_pretty(BufWriter::new(metrics_file), &metrics)?;
    eprintln!(
        "[sbxpress illumina] read={} emitted={} rejected={} metrics={}",
        metrics.records_read,
        metrics.records_emitted,
        metrics.barcode_rejected,
        metrics_path.display()
    );
    if !metrics.failed_outputs.is_empty() {
        let reasons = metrics
            .failed_by_reason
            .iter()
            .map(|(reason, count)| format!("{reason}={count}"))
            .collect::<Vec<_>>()
            .join(" ");
        eprintln!(
            "[sbxpress illumina] failed: {} pair(s) [{}] -> {}",
            metrics.records_failed,
            reasons,
            metrics
                .failed_outputs
                .iter()
                .map(|p| p.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        // Every input record is now accounted for in exactly one place.
        debug_assert_eq!(
            metrics.records_read,
            metrics.records_emitted + metrics.records_failed
        );
    }
    if metrics.overlap_enabled {
        let merged = metrics.overlap_short_sins + metrics.overlap_medium_mins;
        let considered = merged + metrics.overlap_long_lins;
        let pct = |n: u64| {
            if considered > 0 {
                100.0 * n as f64 / considered as f64
            } else {
                0.0
            }
        };
        eprintln!(
            "[sbxpress illumina] overlap: sins(short)={} ({:.1}%) mins(medium)={} ({:.1}%) \
             lins(long)={} ({:.1}%)",
            metrics.overlap_short_sins,
            pct(metrics.overlap_short_sins),
            metrics.overlap_medium_mins,
            pct(metrics.overlap_medium_mins),
            metrics.overlap_long_lins,
            pct(metrics.overlap_long_lins),
        );
        let insert_mean = if merged > 0 {
            metrics
                .insert_length_histogram
                .iter()
                .map(|(len, n)| *len as u64 * n)
                .sum::<u64>() as f64
                / merged as f64
        } else {
            0.0
        };
        eprintln!(
            "[sbxpress illumina] merged: fragments={} mean_cdna={:.1}bp discordance={:.4}% \
             no_cdna={}",
            merged,
            insert_mean,
            if metrics.merged_bases > 0 {
                100.0 * metrics.discordant_bases as f64 / metrics.merged_bases as f64
            } else {
                0.0
            },
            metrics.merged_without_cdna,
        );
        if cli_trim_tso {
            eprintln!(
                "[sbxpress illumina] merged: tso_trimmed={} ({:.1}% of merged)",
                metrics.tso_trimmed,
                if merged > 0 {
                    100.0 * metrics.tso_trimmed as f64 / merged as f64
                } else {
                    0.0
                }
            );
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    run(Cli::parse())
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::read::MultiGzDecoder;
    use std::io::Read;

    #[test]
    fn normalizes_common_mate_suffixes() {
        assert_eq!(canonical_id(b"read/1"), b"read");
        assert_eq!(canonical_id(b"read/2 extra"), b"read");
        assert_eq!(canonical_id(b"read"), b"read");
    }

    /// Clap's `default_value_t` does not feed `Default`, so mirror the defaults
    /// here; tests flip only what they are exercising.
    fn base_cli(root: &Path, output: PathBuf) -> Cli {
        Cli {
            r1: vec![root.join("tests/fixtures/illumina_R1.fastq")],
            r2: vec![root.join("tests/fixtures/illumina_R2.fastq")],
            whitelist: root.join("tests/fixtures/illumina_whitelist.txt"),
            output,
            metrics: None,
            gem_group: 1,
            threads: 1,
            overlap: false,
            overlap_min_overlap: 20,
            overlap_kmer: 8,
            overlap_kmer_offset: 1,
            overlap_check_lengths: 20,
            overlap_max_mismatch: 100,
            overlap_max_mismatch_density: 0.25,
            no_emit_yc: false,
            min_cdna_len: 50,
            adapter_max_mm: 2,
            trim_tso: false,
            failed_output: None,
        }
    }

    fn read_gz(path: &Path) -> String {
        let mut text = String::new();
        MultiGzDecoder::new(File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        text
    }

    #[test]
    fn converts_exact_and_hamming1_and_gates_failure() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let output = temp.path().join("tagged.fastq.gz");
        run(base_cli(&root, output.clone())).unwrap();

        let text = read_gz(&output);
        assert!(text.starts_with("@read_exact\tCR:Z:"));
        assert!(!text.contains("read_exact/2"));
        assert_eq!(text.matches("\n+").count(), 2);
        assert!(text.contains("CR:Z:AAACCCAAGAAACACA"));
        assert!(text.contains("CB:Z:AAACCCAAGAAACACT-1"));
        assert!(!text.contains("read_rejected"));
        // Overlap is off by default, so nothing about the output may change.
        assert!(!text.contains("YC:Z:"));
    }

    #[test]
    fn converts_multiple_pairs_to_numbered_shards() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let output = temp.path().join("tagged.fastq.gz");
        let mut cli = base_cli(&root, output.clone());
        cli.r1 = vec![cli.r1[0].clone(), cli.r1[0].clone()];
        cli.r2 = vec![cli.r2[0].clone(), cli.r2[0].clone()];
        cli.threads = 2;
        run(cli).unwrap();

        assert!(!output.exists());
        assert!(temp.path().join("tagged.shard-001.fastq.gz").is_file());
        assert!(temp.path().join("tagged.shard-002.fastq.gz").is_file());
        let metrics = std::fs::read_to_string(metrics_path(&output)).unwrap();
        assert!(metrics.contains("tagged.shard-001.fastq.gz"));
        assert!(metrics.contains("\"records_read\": 6"));
    }

    /// A long-R1 pair that merges: R1 100 bp, R2 150 bp reading past the CB into
    /// the R1 primer, over an 82 bp cDNA.
    fn write_overlapping_fixture(dir: &Path) -> (PathBuf, PathBuf, Vec<u8>) {
        const CB: &[u8] = b"AAACCCAAGAAACACT"; // on the whitelist fixture
        const UMI: &[u8] = b"ACGTACGTACGT";
        let mut core: Vec<u8> =
            b"CTGATCGATCGTAGCTAGCATCGATCGTAGCTAGCATCGATCGATGCTAGCTAGCATCGATCGTAGCTAGCATCGATCGATCC"
                .to_vec();
        core.truncate(82);
        // End on C so rc(core) does not start with T and blur the dT boundary.
        *core.last_mut().unwrap() = b'C';
        assert_eq!(core.len(), 82);

        let mut molecule = Vec::new();
        molecule.extend_from_slice(CB);
        molecule.extend_from_slice(UMI);
        molecule.extend(std::iter::repeat_n(b'T', 30));
        molecule.extend_from_slice(&revcomp(&core));
        assert_eq!(molecule.len(), 140);

        let mut r1 = molecule.clone();
        r1.truncate(100);
        let mut r2 = revcomp(&molecule);
        r2.extend_from_slice(&revcomp(sbx10x::anchors::R1_PRIMER));
        r2.truncate(150);

        let write = |name: &str, seq: &[u8]| -> PathBuf {
            let path = dir.join(name);
            let text = format!(
                "@frag 1:N:0:1\n{}\n+\n{}\n",
                String::from_utf8_lossy(seq),
                "I".repeat(seq.len())
            );
            std::fs::write(&path, text).unwrap();
            path
        };
        (write("ov_R1.fastq", &r1), write("ov_R2.fastq", &r2), core)
    }

    #[test]
    fn overlap_merges_and_emits_the_full_cdna_with_a_yc_tag() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let (r1, r2, core) = write_overlapping_fixture(temp.path());
        let output = temp.path().join("tagged.fastq.gz");

        let mut cli = base_cli(&root, output.clone());
        cli.r1 = vec![r1];
        cli.r2 = vec![r2];
        cli.overlap = true;
        run(cli).unwrap();

        let text = read_gz(&output);
        let mut lines = text.lines();
        let header = lines.next().unwrap();
        let seq = lines.next().unwrap();

        // The barcode is still tagged exactly as on the unmerged path...
        assert!(header.contains("CR:Z:AAACCCAAGAAACACT"));
        assert!(header.contains("CB:Z:AAACCCAAGAAACACT-1"));
        // ...and the emitted read is now the whole 82 bp cDNA, sense strand,
        // rather than the 150 bp R2 with its barcode read-through.
        assert_eq!(seq.as_bytes(), &core[..]);
        // 5'-most 40 bp came from R2 alone; the remaining 42 bp are duplex.
        assert!(header.ends_with("\tYC:Z:40-42+"), "header was {header}");

        let metrics = std::fs::read_to_string(metrics_path(&output)).unwrap();
        assert!(metrics.contains("\"overlap_enabled\": true"));
        assert!(metrics.contains("\"overlap_short_sins\": 1"));
        assert!(metrics.contains("\"overlap_medium_mins\": 0"));
        assert!(metrics.contains("\"overlap_long_lins\": 0"));
        assert!(metrics.contains("\"82\": 1"), "insert histogram: {metrics}");
    }

    #[test]
    fn failed_output_accounts_for_every_input_read() {
        // The stock fixture has 3 pairs, one of which carries an unmatchable
        // barcode. With --failed-output nothing may be silently dropped.
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let output = temp.path().join("tagged.fastq.gz");
        let failed = temp.path().join("failed");
        let mut cli = base_cli(&root, output.clone());
        cli.failed_output = Some(failed.clone());
        run(cli).unwrap();

        let failed_r1 = temp.path().join("failed_R1.fastq.gz");
        let failed_r2 = temp.path().join("failed_R2.fastq.gz");
        assert!(failed_r1.is_file() && failed_r2.is_file());

        let r1_text = read_gz(&failed_r1);
        let r2_text = read_gz(&failed_r2);
        // Both mates of the failed pair, with the reason and the original header.
        assert!(r1_text.contains("@read_rejected"), "R1 was: {r1_text}");
        assert!(r2_text.contains("@read_rejected"), "R2 was: {r2_text}");
        assert!(r1_text.contains("XF:Z:cb_unmatched"));
        assert_eq!(r1_text.lines().count(), 4);
        assert_eq!(r2_text.lines().count(), 4);

        let metrics: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(metrics_path(&output)).unwrap()).unwrap();
        let read = metrics["records_read"].as_u64().unwrap();
        let emitted = metrics["records_emitted"].as_u64().unwrap();
        let failed_n = metrics["records_failed"].as_u64().unwrap();
        assert_eq!(read, emitted + failed_n, "every read accounted for");
        assert_eq!(failed_n, 1);
        assert_eq!(metrics["failed_by_reason"]["cb_unmatched"], 1);
    }

    #[test]
    fn failed_reads_are_emitted_verbatim() {
        // The side channel must be re-runnable, so the mates come back byte for
        // byte apart from the appended reason tag.
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let failed = temp.path().join("failed");
        let mut cli = base_cli(&root, temp.path().join("tagged.fastq.gz"));
        cli.failed_output = Some(failed);
        run(cli).unwrap();

        let source =
            std::fs::read_to_string(root.join("tests/fixtures/illumina_R1.fastq")).unwrap();
        let rejected: Vec<&str> = source
            .lines()
            .skip_while(|l| !l.starts_with("@read_rejected"))
            .take(4)
            .collect();
        let emitted: Vec<String> = read_gz(&temp.path().join("failed_R1.fastq.gz"))
            .lines()
            .map(str::to_string)
            .collect();
        assert_eq!(emitted[1], rejected[1], "sequence unchanged");
        assert_eq!(emitted[3], rejected[3], "qualities unchanged");
        assert!(
            emitted[0].starts_with(rejected[0]),
            "header preserved, tag appended: {} vs {}",
            emitted[0],
            rejected[0]
        );
    }

    #[test]
    fn both_paths_emit_the_same_orientation() {
        // The same fragment down both routes: once merged, once forced onto the
        // lins path by an unreachable --overlap-min-overlap. The merged read must
        // be the sense-strand cDNA, and R2 (emitted verbatim by the lins path)
        // begins with exactly those bases — so the two paths agree.
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let (r1, r2, core) = write_overlapping_fixture(temp.path());

        let merged_out = temp.path().join("merged.fastq.gz");
        let mut cli = base_cli(&root, merged_out.clone());
        cli.r1 = vec![r1.clone()];
        cli.r2 = vec![r2.clone()];
        cli.overlap = true;
        run(cli).unwrap();

        let lins_out = temp.path().join("lins.fastq.gz");
        let mut cli = base_cli(&root, lins_out.clone());
        cli.r1 = vec![r1];
        cli.r2 = vec![r2];
        cli.overlap = true;
        cli.overlap_min_overlap = 200; // longer than either read: cannot merge
        run(cli).unwrap();

        let merged_seq = read_gz(&merged_out).lines().nth(1).unwrap().to_string();
        let lins_seq = read_gz(&lins_out).lines().nth(1).unwrap().to_string();

        assert_eq!(
            merged_seq.as_bytes(),
            &core[..],
            "merged is sense-strand cDNA"
        );
        assert!(
            lins_seq.starts_with(&merged_seq),
            "the lins path emits R2 verbatim, which must begin with the same \
             bases in the same orientation\n  merged: {merged_seq}\n  lins:   {lins_seq}"
        );
        // The merge recovers the cDNA and drops the polyA/barcode read-through
        // that the lins path leaves on the end of R2.
        assert!(lins_seq.len() > merged_seq.len());
    }

    #[test]
    fn a_non_overlapping_pair_falls_back_to_the_standard_path() {
        // The stock fixtures are far too short to overlap, so with --overlap on
        // they must all be counted as lins and emitted exactly as before.
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let plain = temp.path().join("plain.fastq.gz");
        run(base_cli(&root, plain.clone())).unwrap();

        let overlapped = temp.path().join("overlapped.fastq.gz");
        let mut cli = base_cli(&root, overlapped.clone());
        cli.overlap = true;
        run(cli).unwrap();

        assert_eq!(read_gz(&plain), read_gz(&overlapped));
        let metrics = std::fs::read_to_string(metrics_path(&overlapped)).unwrap();
        assert!(metrics.contains("\"overlap_short_sins\": 0"));
        assert!(metrics.contains("\"overlap_medium_mins\": 0"));
        assert!(metrics.contains("\"overlap_long_lins\": 3"));
    }
}
