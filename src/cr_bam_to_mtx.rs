//! Reconstruct a Cell Ranger raw gene-expression matrix from its barcoded BAM.

use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read as IoRead, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use rust_htslib::bam::{self, record::Aux, Read};
use serde::Serialize;
use tracing::info;

/// Reconstruct Cell Ranger's raw GEX matrix from a Cell Ranger BAM.
#[derive(Debug, Parser)]
#[command(
    name = "cr-bam-to-mtx",
    version,
    about = "Reconstruct a Cell Ranger raw gene-expression matrix from its barcoded BAM"
)]
pub struct Cli {
    /// Cell Ranger possorted_genome_bam.bam input.
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,

    /// Output directory for matrix files and metrics_summary.json.
    #[arg(short, long, default_value = "cr_bam_to_mtx_out")]
    pub output: PathBuf,

    /// Optional Cell Ranger features.tsv[.gz] template, including zero-count rows.
    #[arg(long)]
    pub features: Option<PathBuf>,

    /// Optional Cell Ranger barcodes.tsv[.gz] template, including absent columns.
    #[arg(long)]
    pub barcodes: Option<PathBuf>,

    /// Write plain .tsv/.mtx instead of gzip-compressed files.
    #[arg(long)]
    pub no_gzip: bool,
}

#[derive(Debug, Clone, Eq, PartialEq)]
struct Feature {
    id: String,
    name: String,
}

#[derive(Debug, Default)]
struct Matrix {
    features: Vec<Feature>,
    barcodes: Vec<String>,
    counts: BTreeMap<(usize, usize), u64>,
    records: u64,
    molecules: u64,
}

/// Machine-readable run accounting, following the main `sbx10x` app.
#[derive(Debug, Serialize)]
struct Metrics {
    input_bam: String,
    output_dir: String,
    elapsed_secs: f64,
    records_read: u64,
    records_per_sec: u64,
    molecules_counted: u64,
    features: usize,
    barcodes: usize,
    nonzero_entries: usize,
    used_feature_template: bool,
    used_barcode_template: bool,
}

/// Run the BAM-to-matrix application.
pub fn run(cli: &Cli) -> Result<()> {
    let started = Instant::now();
    let matrix = read_matrix(&cli.bam, cli.features.as_deref(), cli.barcodes.as_deref())?;
    write_matrix(&matrix, &cli.output, !cli.no_gzip)?;

    let elapsed = started.elapsed().as_secs_f64();
    let metrics = Metrics {
        input_bam: cli.bam.display().to_string(),
        output_dir: cli.output.display().to_string(),
        elapsed_secs: elapsed,
        records_read: matrix.records,
        records_per_sec: if elapsed > 0.0 {
            (matrix.records as f64 / elapsed) as u64
        } else {
            0
        },
        molecules_counted: matrix.molecules,
        features: matrix.features.len(),
        barcodes: matrix.barcodes.len(),
        nonzero_entries: matrix.counts.len(),
        used_feature_template: cli.features.is_some(),
        used_barcode_template: cli.barcodes.is_some(),
    };
    let metrics_path = cli.output.join("metrics_summary.json");
    fs::write(&metrics_path, serde_json::to_string_pretty(&metrics)?)
        .with_context(|| format!("writing {}", metrics_path.display()))?;

    info!(
        elapsed = format!("{elapsed:.1}s"),
        records = metrics.records_read,
        records_per_sec = metrics.records_per_sec,
        molecules = metrics.molecules_counted,
        features = metrics.features,
        barcodes = metrics.barcodes,
        nonzero_entries = metrics.nonzero_entries,
        output = %cli.output.display(),
        "RUN SUMMARY"
    );
    Ok(())
}

fn aux_string(record: &bam::Record, tag: &[u8; 2]) -> Result<Option<String>> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Ok(Some(value.to_owned())),
        Ok(_) => bail!("tag {}{} is not a string", tag[0] as char, tag[1] as char),
        Err(_) => Ok(None),
    }
}

fn aux_i64(record: &bam::Record, tag: &[u8; 2]) -> Result<Option<i64>> {
    let value = match record.aux(tag) {
        Ok(Aux::I8(v)) => v as i64,
        Ok(Aux::U8(v)) => v as i64,
        Ok(Aux::I16(v)) => v as i64,
        Ok(Aux::U16(v)) => v as i64,
        Ok(Aux::I32(v)) => v as i64,
        Ok(Aux::U32(v)) => v as i64,
        Ok(_) => bail!("tag {}{} is not an integer", tag[0] as char, tag[1] as char),
        Err(_) => return Ok(None),
    };
    Ok(Some(value))
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

fn read_features(path: &Path) -> Result<Vec<Feature>> {
    text_reader(path)?
        .lines()
        .map(|line| {
            let line = line?;
            let mut fields = line.split('\t');
            let id = fields.next().context("feature row has no ID")?.to_owned();
            let name = fields.next().unwrap_or(&id).to_owned();
            Ok(Feature { id, name })
        })
        .collect()
}

fn read_barcodes(path: &Path) -> Result<Vec<String>> {
    Ok(text_reader(path)?
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?)
}

fn read_matrix(
    path: &Path,
    feature_template: Option<&Path>,
    barcode_template: Option<&Path>,
) -> Result<Matrix> {
    let mut reader =
        bam::Reader::from_path(path).with_context(|| format!("opening BAM {}", path.display()))?;
    let mut feature_order = feature_template
        .map(read_features)
        .transpose()?
        .unwrap_or_default();
    let mut feature_indexes: HashMap<String, usize> = feature_order
        .iter()
        .enumerate()
        .map(|(index, feature)| (feature.id.clone(), index))
        .collect();
    let mut molecule_keys: Vec<(usize, String)> = Vec::new();
    let mut records = 0_u64;

    for result in reader.records() {
        let record = result.context("reading BAM record")?;
        records += 1;
        if aux_i64(&record, b"xf")?.unwrap_or(0) & 8 == 0 {
            continue;
        }
        let cb = aux_string(&record, b"CB")?
            .context("xf bit 8 record is missing corrected barcode tag CB")?;
        let _ub = aux_string(&record, b"UB")?
            .context("xf bit 8 record is missing corrected UMI tag UB")?;
        let gx =
            aux_string(&record, b"GX")?.context("xf bit 8 record is missing gene ID tag GX")?;
        let gn = aux_string(&record, b"GN")?.unwrap_or_else(|| gx.clone());
        if gx.contains(';') {
            bail!("xf bit 8 record has multiple genes in GX: {gx}");
        }
        let feature_index = if let Some(&index) = feature_indexes.get(&gx) {
            index
        } else {
            let index = feature_order.len();
            feature_indexes.insert(gx.clone(), index);
            feature_order.push(Feature { id: gx, name: gn });
            index
        };
        molecule_keys.push((feature_index, cb));
    }

    let barcodes = if let Some(path) = barcode_template {
        read_barcodes(path)?
    } else {
        let mut values: Vec<String> = molecule_keys.iter().map(|(_, cb)| cb.clone()).collect();
        values.sort_unstable();
        values.dedup();
        values
    };
    let barcode_indexes: HashMap<&str, usize> = barcodes
        .iter()
        .enumerate()
        .map(|(index, barcode)| (barcode.as_str(), index))
        .collect();
    let mut counts = BTreeMap::new();
    for (feature, barcode) in &molecule_keys {
        let Some(&barcode) = barcode_indexes.get(barcode.as_str()) else {
            bail!("counted barcode {barcode} is absent from barcode template");
        };
        *counts.entry((*feature, barcode)).or_insert(0) += 1;
    }

    Ok(Matrix {
        features: feature_order,
        barcodes,
        counts,
        records,
        molecules: molecule_keys.len() as u64,
    })
}

fn writer(path: &Path, gzip: bool) -> Result<Box<dyn Write>> {
    let file = File::create(path).with_context(|| format!("creating {}", path.display()))?;
    if gzip {
        Ok(Box::new(GzEncoder::new(file, Compression::default())))
    } else {
        Ok(Box::new(BufWriter::new(file)))
    }
}

fn write_matrix(matrix: &Matrix, output: &Path, gzip: bool) -> Result<()> {
    fs::create_dir_all(output).with_context(|| format!("creating {}", output.display()))?;
    let suffix = if gzip { ".gz" } else { "" };
    let mut features = writer(&output.join(format!("features.tsv{suffix}")), gzip)?;
    for feature in &matrix.features {
        writeln!(
            features,
            "{}\t{}\tGene Expression",
            feature.id, feature.name
        )?;
    }
    let mut barcodes = writer(&output.join(format!("barcodes.tsv{suffix}")), gzip)?;
    for barcode in &matrix.barcodes {
        writeln!(barcodes, "{barcode}")?;
    }
    let mut mtx = writer(&output.join(format!("matrix.mtx{suffix}")), gzip)?;
    writeln!(mtx, "%%MatrixMarket matrix coordinate integer general")?;
    writeln!(
        mtx,
        "%metadata_json: {{\"software_version\": \"cr-bam-to-mtx-{}\", \"format_version\": 2}}",
        env!("CARGO_PKG_VERSION")
    )?;
    writeln!(
        mtx,
        "{} {} {}",
        matrix.features.len(),
        matrix.barcodes.len(),
        matrix.counts.len()
    )?;
    for (&(feature, barcode), &count) in &matrix.counts {
        writeln!(mtx, "{} {} {}", feature + 1, barcode + 1, count)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_one_based_sparse_matrix() {
        let matrix = Matrix {
            features: vec![Feature {
                id: "g1".into(),
                name: "G1".into(),
            }],
            barcodes: vec!["bc-1".into()],
            counts: [((0, 0), 2)].into_iter().collect(),
            records: 3,
            molecules: 2,
        };
        let dir = tempfile::tempdir().unwrap();
        write_matrix(&matrix, dir.path(), false).unwrap();
        let text = fs::read_to_string(dir.path().join("matrix.mtx")).unwrap();
        assert!(text.ends_with("1 1 1\n1 1 2\n"));
    }
}
