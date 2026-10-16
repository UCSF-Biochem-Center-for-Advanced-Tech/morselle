use crate::barcode_metrics::BarcodeTable;
use anyhow::{bail, Context, Result};
use clap::Parser;
use crossbeam_channel::{bounded, Sender};
use flate2::write::GzEncoder;
use flate2::Compression;
use serde::Serialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::thread::{self, JoinHandle};

pub const CB_LEN: usize = 16;
pub const UMI_LEN: usize = 12;
pub const WRITER_QUEUE_DEPTH: usize = 4;
pub const WRITER_BATCH_BYTES: usize = 256 * 1024;

#[derive(Debug, Parser)]
#[command(about = "Convert paired Illumina 10x FASTQs to minimap2-tagged R2 FASTQ")]
pub struct Cli {
    /// Illumina Read 1 FASTQ. Repeat once per lane/shard, matching --r2 order.
    #[arg(long, required = true, action = clap::ArgAction::Append)]
    pub r1: Vec<PathBuf>,

    /// Illumina Read 2 FASTQ. Repeat once per lane/shard, matching --r1 order.
    #[arg(long, required = true, action = clap::ArgAction::Append)]
    pub r2: Vec<PathBuf>,

    /// 10x cell-barcode whitelist (plain text or .gz).
    #[arg(long)]
    pub whitelist: PathBuf,

    /// Output tagged R2 FASTQ (.gz recommended).
    #[arg(short, long)]
    pub output: PathBuf,

    /// JSON metrics path. Defaults to <output>.metrics.json.
    #[arg(long)]
    pub metrics: Option<PathBuf>,

    /// Also write per-barcode counters to `<metrics>.barcodes.tsv.gz` (MANUAL
    /// "Per-barcode read metrics"). Off by default; it adds run time and memory.
    #[arg(long)]
    pub barcode_metrics: bool,

    /// GEM-group suffix placed on corrected CB tags.
    #[arg(long, default_value_t = 1)]
    pub gem_group: u32,

    /// Maximum number of R1/R2 file pairs converted concurrently.
    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,

    /// Total gzip shard budget, apportioned across all active output targets.
    #[arg(long, default_value_t = 1)]
    pub split: usize,

    /// Merge overlapping mates into the full fragment before extraction.
    #[arg(long, default_value_t = false)]
    pub overlap: bool,

    /// Shortest mate overlap that may be merged.
    #[arg(long, default_value_t = 20)]
    pub overlap_min_overlap: usize,

    /// Seed length for the overlap lookup table.
    #[arg(long, default_value_t = 8)]
    pub overlap_kmer: usize,

    /// Stride between indexed overlap seeds.
    #[arg(long, default_value_t = 1)]
    pub overlap_kmer_offset: usize,

    /// Usable seeds probed inward from each end of the shorter mate.
    #[arg(long, default_value_t = 20)]
    pub overlap_check_lengths: usize,

    /// Absolute cap on mismatches within an overlap.
    #[arg(long, default_value_t = 100)]
    pub overlap_max_mismatch: usize,

    /// Cap on mismatches as a fraction of the overlap length.
    #[arg(long, default_value_t = 0.25)]
    pub overlap_max_mismatch_density: f64,

    /// Omit the YC tag from merged reads.
    #[arg(long, default_value_t = false)]
    pub no_emit_yc: bool,

    /// Minimum cDNA length to emit.
    #[arg(long, default_value_t = 50)]
    pub min_cdna_len: usize,

    /// Max mismatches for the 3' adapter seed when extracting a merged fragment.
    #[arg(long, default_value_t = 2)]
    pub adapter_max_mm: u32,

    /// Emit every read that fails, as the original R1/R2 pair.
    #[arg(long)]
    pub failed_output: Option<PathBuf>,

    /// Trim the template-switch oligo: from merged fragments, and from the
    /// start of a non-merged R2 (3' GEX: R2 starts in the TSO on fragments
    /// that reach the transcript's 5' end). Off by default: minimap2 already
    /// soft-clips the TSO, and on CAT02 trimming it lowered MAPQ and cost
    /// 0.6% of countable reads (CHANGELOG 0.6.3).
    #[arg(long, default_value_t = false)]
    pub trim_tso: bool,

    /// Trim poly-A read-through in a non-merged R2. On by default (3' GEX
    /// and tag-seq both read into the tail); `--trim-r2-polya false` turns
    /// it off.
    #[arg(long, default_value_t = true, action = clap::ArgAction::Set,
          num_args = 0..=1, default_missing_value = "true")]
    pub trim_r2_polya: bool,

    /// Shortest poly-A run in R2 that counts as tail read-through (at least
    /// 90% A, isolated mismatches only).
    #[arg(long, default_value_t = 12)]
    pub r2_polya_min_run: usize,

    /// Low-complexity cDNA check (dropped by default).
    #[command(flatten)]
    pub complexity: crate::complexity::ComplexityCfg,
}

#[derive(Debug, Default, Serialize)]
pub struct Metrics {
    pub r1_files: Vec<PathBuf>,
    pub r2_files: Vec<PathBuf>,
    pub whitelist: PathBuf,
    pub whitelist_size: usize,
    pub output: PathBuf,
    pub output_shards: Vec<PathBuf>,
    pub split: usize,
    pub gem_group: u32,
    pub records_read: u64,
    pub records_emitted: u64,
    pub barcode_exact: u64,
    pub barcode_hamming1_corrected: u64,
    pub barcode_rejected: u64,
    pub short_r1: u64,
    pub overlap_enabled: bool,
    /// Merged, fragment shorter than a read: read-through was present on a mate.
    pub overlap_short_sins: u64,
    /// Merged, neither mate ran off the end of the molecule.
    pub overlap_medium_mins: u64,
    /// Not merged: the mates did not reach each other.
    pub overlap_long_lins: u64,
    /// Merged, but no cDNA survived extraction (below `--min-cdna-len`).
    pub merged_without_cdna: u64,
    pub merged_bases: u64,
    pub discordant_bases: u64,
    pub yc_n_columns_coerced: u64,
    pub yc_unencodable: u64,
    /// Merged fragments whose TSO was found and trimmed (needs `--trim-tso`).
    pub tso_trimmed: u64,
    /// Non-merged R2 reads trimmed at a detected poly-A read-through (needs
    /// `--trim-r2-polya`).
    pub r2_polya_trimmed: u64,
    /// Of those, the cDNA remaining before the poly-A run was below
    /// `--min-cdna-len`, so the pair was discarded instead of emitted (see
    /// `r2_polya_short` in `failed_by_reason`).
    pub r2_polya_discarded: u64,
    pub r2_polya_bases_trimmed: u64,
    /// Non-merged R2 reads that started with the TSO (whole or a suffix of at
    /// least 12 bases) and had it trimmed (needs `--trim-tso`).
    pub r2_tso_trimmed: u64,
    /// Of those, the cDNA left after the TSO was below `--min-cdna-len`, so
    /// the pair was discarded (`r2_tso_short` in `failed_by_reason`).
    pub r2_tso_discarded: u64,
    pub r2_tso_bases_trimmed: u64,
    /// Thresholds used for the low-complexity check.
    pub low_complexity_base_fraction: f64,
    pub low_complexity_run_fraction: f64,
    pub drop_low_complexity: bool,
    /// Barcode-passing pairs whose cDNA failed the low-complexity check
    /// (emitted tagged `YL:i:1`, or failed as `low_complexity` with
    /// `--drop-low-complexity`).
    pub low_complexity: u64,
    /// `low_complexity` by the read's most frequent base.
    pub low_complexity_by_base: crate::complexity::ByBase,
    /// Reads written to the failed-output side channel, by reason.
    pub records_failed: u64,
    pub failed_by_reason: BTreeMap<String, u64>,
    pub failed_outputs: Vec<PathBuf>,
    /// cDNA length of each merged fragment, excluding CB+UMI+dT.
    pub insert_length_histogram: BTreeMap<usize, u64>,
    pub elapsed_secs: f64,
}

#[derive(Debug, Default)]
pub struct PairMetrics {
    pub per_barcode: BarcodeTable<13>,
    pub records_read: u64,
    pub records_emitted: u64,
    pub barcode_exact: u64,
    pub barcode_hamming1_corrected: u64,
    pub barcode_rejected: u64,
    pub short_r1: u64,
    pub overlap_short: u64,
    pub overlap_medium: u64,
    pub overlap_long: u64,
    pub merged_without_cdna: u64,
    pub merged_bases: u64,
    pub discordant_bases: u64,
    pub yc_n_columns_coerced: u64,
    pub yc_unencodable: u64,
    pub tso_trimmed: u64,
    pub r2_polya_trimmed: u64,
    pub r2_polya_discarded: u64,
    pub r2_polya_bases_trimmed: u64,
    pub r2_tso_trimmed: u64,
    pub r2_tso_discarded: u64,
    pub r2_tso_bases_trimmed: u64,
    pub low_complexity: u64,
    pub low_complexity_by_base: crate::complexity::ByBase,
    pub records_failed: u64,
    pub failed_by_reason: BTreeMap<String, u64>,
    pub insert_length_histogram: BTreeMap<usize, u64>,
}

pub fn canonical_id(id: &[u8]) -> &[u8] {
    let id = id.split(|b| b.is_ascii_whitespace()).next().unwrap_or(id);
    if id.ends_with(b"/1") || id.ends_with(b"/2") {
        &id[..id.len() - 2]
    } else {
        id
    }
}

pub fn metrics_path(output: &Path) -> PathBuf {
    let mut name = output.as_os_str().to_os_string();
    name.push(".metrics.json");
    PathBuf::from(name)
}

pub fn shard_path(output: &Path, index: usize, count: usize) -> PathBuf {
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

/// Evenly apportion a total shard budget among logical output targets. Earlier
/// targets receive the one-slot remainder deterministically.
pub fn shard_allocation(total: usize, targets: usize) -> Result<Vec<usize>> {
    if total < targets {
        bail!("--split ({total}) must be at least the number of active output targets ({targets})");
    }
    Ok((0..targets)
        .map(|i| total / targets + usize::from(i < total % targets))
        .collect())
}

pub struct AsyncWriter {
    pub sender: Sender<Vec<u8>>,
    pub join: JoinHandle<Result<()>>,
}

impl AsyncWriter {
    pub fn create(path: &Path) -> Result<Self> {
        let file =
            File::create(path).with_context(|| format!("cannot create {}", path.display()))?;
        let (sender, receiver) = bounded::<Vec<u8>>(WRITER_QUEUE_DEPTH);
        let join = thread::spawn(move || -> Result<()> {
            let mut output = GzEncoder::new(BufWriter::new(file), Compression::new(4));
            for bytes in receiver {
                output.write_all(&bytes)?;
            }
            output.finish()?;
            Ok(())
        });
        Ok(Self { sender, join })
    }

    pub fn finish(self) -> Result<()> {
        drop(self.sender);
        self.join
            .join()
            .map_err(|_| anyhow::anyhow!("gzip writer panicked"))?
    }
}

pub struct AsyncFailedWriter {
    pub sender: Sender<(Vec<u8>, Vec<u8>)>,
    pub join: JoinHandle<Result<()>>,
}

impl AsyncFailedWriter {
    pub fn create(r1_path: &Path, r2_path: &Path) -> Result<Self> {
        let r1 = File::create(r1_path)
            .with_context(|| format!("cannot create {}", r1_path.display()))?;
        let r2 = File::create(r2_path)
            .with_context(|| format!("cannot create {}", r2_path.display()))?;
        let (sender, receiver) = bounded::<(Vec<u8>, Vec<u8>)>(WRITER_QUEUE_DEPTH);
        let join = thread::spawn(move || -> Result<()> {
            let mut r1 = GzEncoder::new(BufWriter::new(r1), Compression::new(4));
            let mut r2 = GzEncoder::new(BufWriter::new(r2), Compression::new(4));
            for (a, b) in receiver {
                r1.write_all(&a)?;
                r2.write_all(&b)?;
            }
            r1.finish()?;
            r2.finish()?;
            Ok(())
        });
        Ok(Self { sender, join })
    }

    pub fn finish(self) -> Result<()> {
        drop(self.sender);
        self.join
            .join()
            .map_err(|_| anyhow::anyhow!("failed writer panicked"))?
    }
}

/// Paths for a failed-pair shard: `<prefix>_R1.fastq.gz` / `<prefix>_R2.fastq.gz`,
/// carrying the same `.shard-NNN` infix as the main output when sharded.
pub fn failed_paths(prefix: &Path, index: usize, count: usize) -> (PathBuf, PathBuf) {
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
