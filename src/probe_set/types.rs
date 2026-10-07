use std::fs;
use std::path::{Path, PathBuf};
use anyhow::{ensure, Context, Result};
use clap::Parser;
use serde::Serialize;

pub const SYNTHETIC_CONTIG_LEN: i64 = 1_000_000;
pub const DEFAULT_PROBE_SEQ_OFFSET: usize = 0;

/// A probe-set reference row.
#[derive(Debug, Clone)]
pub struct Probe {
    pub gene_id: String,
    pub gene_name: String,
    pub probe_id: String,
    pub region: String,
    pub seq: Vec<u8>,
}

#[derive(Debug)]
pub struct ProbeSetReference {
    pub probes: Vec<Probe>,
}

impl ProbeSetReference {
    pub fn load(path: &Path) -> Result<Self> {
        let text = fs::read_to_string(path)
            .with_context(|| format!("reading probe set {}", path.display()))?;
        let mut lines = text.lines().filter(|l| !l.starts_with('#'));
        let header = lines.next().context("empty probe set (no header)")?;
        let cols: Vec<&str> = header.split(',').collect();
        let col = |name: &str| {
            cols.iter()
                .position(|c| *c == name)
                .with_context(|| format!("probe set missing '{name}' column"))
        };
        let (gene_id_i, seq_i, probe_id_i, included_i, region_i, gene_name_i) = (
            col("gene_id")?,
            col("probe_seq")?,
            col("probe_id")?,
            col("included")?,
            col("region")?,
            col("gene_name")?,
        );

        let mut probes = Vec::new();
        let mut seq_len = None;
        for (lineno, line) in lines.enumerate() {
            if line.trim().is_empty() {
                continue;
            }
            let f: Vec<&str> = line.split(',').collect();
            let max = [
                gene_id_i,
                seq_i,
                probe_id_i,
                included_i,
                region_i,
                gene_name_i,
            ]
            .into_iter()
            .max()
            .unwrap();
            ensure!(
                f.len() > max,
                "{}:{}: row has too few columns: '{line}'",
                path.display(),
                lineno + 2
            );
            if f[included_i] != "TRUE" {
                continue;
            }
            let seq = f[seq_i].as_bytes().to_vec();
            match seq_len {
                None => seq_len = Some(seq.len()),
                Some(len) => ensure!(
                    seq.len() == len,
                    "{}:{}: probe_seq length {} != other probes' length {len} -- this parser assumes one fixed length",
                    path.display(),
                    lineno + 2,
                    seq.len()
                ),
            }
            probes.push(Probe {
                gene_id: f[gene_id_i].to_string(),
                gene_name: f[gene_name_i].to_string(),
                probe_id: f[probe_id_i].to_string(),
                region: f[region_i].to_string(),
                seq,
            });
        }
        ensure!(
            !probes.is_empty(),
            "probe set {} has no included=TRUE rows",
            path.display()
        );
        Ok(Self { probes })
    }
}

/// Assign Flex GEX reads to genes by probe-set lookup, not genome alignment.
#[derive(Debug, Parser)]
#[command(
    name = "sbxpress-probeset-align",
    version,
    about = "Assign 10x Flex GEX reads to genes via the real probe-set reference (not genome alignment)"
)]
pub struct Cli {
    #[arg(value_name = "R2_FASTQ")]
    pub r2_fastq: PathBuf,

    #[arg(long)]
    pub probe_set: PathBuf,

    #[arg(short, long, value_name = "BAM")]
    pub output: PathBuf,

    #[arg(long)]
    pub genes_out: Option<PathBuf>,

    #[arg(long)]
    pub metrics: Option<PathBuf>,

    #[arg(long, default_value_t = 1)]
    pub max_mm: u32,

    #[arg(long, default_value_t = 15)]
    pub window: isize,

    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,
}

#[derive(Debug, Serialize)]
pub struct Metrics {
    pub reads: u64,
    pub matched: u64,
    pub hit_rate_pct: f64,
    pub distinct_genes_hit: usize,
}
