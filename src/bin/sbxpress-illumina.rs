//! Convert standard paired Illumina 10x GEX FASTQs into the tagged R2-only
//! FASTQ contract consumed by `minimap2 -y` in the sbxpress workflow.

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
use sbx10x::barcode::{encode16, Whitelist};
use serde::Serialize;

const CB_LEN: usize = 16;
const UMI_LEN: usize = 12;

#[derive(Debug, Parser)]
#[command(about = "Convert paired Illumina 10x FASTQs to minimap2-tagged R2 FASTQ")]
struct Cli {
    /// Illumina Read 1 FASTQ. Repeat once per lane/shard, in R2 order.
    #[arg(long, required = true, action = clap::ArgAction::Append)]
    r1: Vec<PathBuf>,

    /// Illumina Read 2 FASTQ. Repeat once per lane/shard, in R1 order.
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

fn convert_pair(
    r1_path: &Path,
    r2_path: &Path,
    output_path: &Path,
    whitelist: &Whitelist,
    gem_group: u32,
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
    let mut pair_number = 0u64;
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
                let Some(r1_qual) = a.qual() else {
                    bail!("R1 record {} has no qualities", pair_number)
                };
                let Some(r2_qual) = b.qual() else {
                    bail!("R2 record {} has no qualities", pair_number)
                };
                if r1_seq.len() < CB_LEN + UMI_LEN || r1_qual.len() < CB_LEN + UMI_LEN {
                    metrics.short_r1 += 1;
                    metrics.barcode_rejected += 1;
                    continue;
                }
                let raw_cb = &r1_seq[..CB_LEN];
                let raw_umi = &r1_seq[CB_LEN..CB_LEN + UMI_LEN];
                let corrected = if encode16(raw_cb).is_none() {
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
                let Some(corrected) = corrected else {
                    metrics.barcode_rejected += 1;
                    continue;
                };
                // Illumina FASTQ identifiers commonly include a whitespace-delimited
                // mate/control/index comment (for example `2:N:0:...`). Minimap2 -y
                // treats everything after the first whitespace as SAM-style auxiliary
                // fields, so emit only the canonical cluster name before our tags.
                write!(output, "@{}", String::from_utf8_lossy(canonical_id(b.id())))?;
                write!(output, "\tCR:Z:{}", String::from_utf8_lossy(raw_cb))?;
                write!(
                    output,
                    "\tCB:Z:{}-{}",
                    String::from_utf8_lossy(&corrected),
                    gem_group
                )?;
                write!(
                    output,
                    "\tUR:Z:{}\tUB:Z:{}",
                    String::from_utf8_lossy(raw_umi),
                    String::from_utf8_lossy(raw_umi)
                )?;
                write!(
                    output,
                    "\tCY:Z:{}",
                    String::from_utf8_lossy(&r1_qual[..CB_LEN])
                )?;
                writeln!(
                    output,
                    "\tUY:Z:{}",
                    String::from_utf8_lossy(&r1_qual[CB_LEN..CB_LEN + UMI_LEN])
                )?;
                output.write_all(&b.seq())?;
                output.write_all(b"\n+\n")?;
                output.write_all(r2_qual)?;
                output.write_all(b"\n")?;
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
    let whitelist = Whitelist::load(&cli.whitelist)?;
    if whitelist.is_empty() {
        bail!("whitelist is empty: {}", cli.whitelist.display());
    }
    let output_shards = (0..cli.r1.len())
        .map(|index| shard_path(&cli.output, index, cli.r1.len()))
        .collect::<Vec<_>>();
    let mut metrics = Metrics {
        r1_files: cli.r1.clone(),
        r2_files: cli.r2.clone(),
        whitelist: cli.whitelist.clone(),
        whitelist_size: whitelist.len(),
        output: cli.output.clone(),
        output_shards: output_shards.clone(),
        gem_group: cli.gem_group,
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
                let result = convert_pair(r1, r2, output, &whitelist, cli.gem_group);
                if let Ok(value) = &result {
                    eprintln!(
                        "[sbxpress illumina] shard {}/{} complete: read={} emitted={} output={}",
                        index + 1,
                        cli.r1.len(),
                        value.records_read,
                        value.records_emitted,
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

    #[test]
    fn converts_exact_and_hamming1_and_gates_failure() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let output = temp.path().join("tagged.fastq.gz");
        run(Cli {
            r1: vec![root.join("tests/fixtures/illumina_R1.fastq")],
            r2: vec![root.join("tests/fixtures/illumina_R2.fastq")],
            whitelist: root.join("tests/fixtures/illumina_whitelist.txt"),
            output: output.clone(),
            metrics: None,
            gem_group: 1,
            threads: 1,
        })
        .unwrap();

        let mut text = String::new();
        MultiGzDecoder::new(File::open(output).unwrap())
            .read_to_string(&mut text)
            .unwrap();
        assert!(text.starts_with("@read_exact\tCR:Z:"));
        assert!(!text.contains("read_exact/2"));
        assert_eq!(text.matches("\n+").count(), 2);
        assert!(text.contains("CR:Z:AAACCCAAGAAACACA"));
        assert!(text.contains("CB:Z:AAACCCAAGAAACACT-1"));
        assert!(!text.contains("read_rejected"));
    }

    #[test]
    fn converts_multiple_pairs_to_numbered_shards() {
        let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let temp = tempfile::tempdir().unwrap();
        let output = temp.path().join("tagged.fastq.gz");
        let r1 = root.join("tests/fixtures/illumina_R1.fastq");
        let r2 = root.join("tests/fixtures/illumina_R2.fastq");
        run(Cli {
            r1: vec![r1.clone(), r1],
            r2: vec![r2.clone(), r2],
            whitelist: root.join("tests/fixtures/illumina_whitelist.txt"),
            output: output.clone(),
            metrics: None,
            gem_group: 1,
            threads: 2,
        })
        .unwrap();

        assert!(!output.exists());
        assert!(temp.path().join("tagged.shard-001.fastq.gz").is_file());
        assert!(temp.path().join("tagged.shard-002.fastq.gz").is_file());
        let metrics = std::fs::read_to_string(metrics_path(&output)).unwrap();
        assert!(metrics.contains("tagged.shard-001.fastq.gz"));
        assert!(metrics.contains("\"records_read\": 6"));
    }
}
