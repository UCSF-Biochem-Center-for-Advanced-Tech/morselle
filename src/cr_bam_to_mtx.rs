//! Reconstruct a 10x reference raw gene-expression matrix from its barcoded BAM.

use std::collections::{BTreeMap, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Read as IoRead, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{bail, Context, Result};
use clap::Parser;
use flate2::{read::MultiGzDecoder, write::GzEncoder, Compression};
use hdf5::types::VarLenAscii;
use rust_htslib::bam::{self, record::Aux, Read};
use serde::Serialize;
use tracing::info;

/// Reconstruct 10x reference's raw GEX matrix from a 10x reference BAM.
#[derive(Debug, Parser)]
#[command(
    name = "cr-bam-to-mtx",
    version,
    about = "Reconstruct a 10x reference raw gene-expression matrix from its barcoded BAM"
)]
pub struct Cli {
    /// 10x reference possorted_genome_bam.bam input.
    #[arg(value_name = "BAM")]
    pub bam: PathBuf,

    /// 10x-style outs directory.
    #[arg(short, long, default_value = "outs")]
    pub output: PathBuf,

    /// Genome label stored for every feature in raw_feature_bc_matrix.h5.
    #[arg(long, default_value = "unknown")]
    pub genome: String,

    /// Optional 10x reference features.tsv[.gz] template, including zero-count rows.
    #[arg(long)]
    pub features: Option<PathBuf>,

    /// Optional 10x reference barcodes.tsv[.gz] template, including absent columns.
    #[arg(long, conflicts_with = "cell_barcodes")]
    pub barcodes: Option<PathBuf>,

    /// 10x reference barcode whitelist (plain or .gz). Preserves whitelist order,
    /// adds the gem-group suffix, and emits zero-count columns for absent cells.
    #[arg(long, value_name = "WHITELIST", conflicts_with = "barcodes")]
    pub cell_barcodes: Option<PathBuf>,

    /// Gem-group suffix added to entries read via --cell-barcodes.
    #[arg(long, default_value_t = 1)]
    pub gem_group: u32,

    /// With --barcodes/--cell-barcodes, silently skip a counted molecule
    /// whose barcode isn't in the template instead of erroring. Off by
    /// default (a barcode outside the template usually means a genuine
    /// whitelist/BAM mismatch worth failing loudly on) — turn on when the
    /// BAM is intentionally a pooled multi-sample input and the template is
    /// deliberately just one sample's slice of it (e.g. 10x Flex per-sample
    /// splitting via `probe-demux`'s `<sample>.barcodes.txt`).
    #[arg(long, default_value_t = false)]
    pub skip_unlisted_barcodes: bool,

    /// Write plain .tsv/.mtx instead of gzip-compressed files.
    #[arg(long)]
    pub no_gzip: bool,
}

#[derive(Debug, Clone, Eq, PartialEq)]
pub(crate) struct Feature {
    pub(crate) id: String,
    pub(crate) name: String,
}

#[derive(Debug, Default)]
pub(crate) struct Matrix {
    pub(crate) features: Vec<Feature>,
    pub(crate) barcodes: Vec<String>,
    /// (barcode index, feature index) -> molecule count. Keeping the barcode
    /// first makes Matrix Market output match 10x reference's column-major order.
    pub(crate) counts: BTreeMap<(usize, usize), u64>,
    pub(crate) records: u64,
    pub(crate) molecules: u64,
    /// Counted molecules skipped because their barcode wasn't in the
    /// template (only nonzero with `--skip-unlisted-barcodes`).
    pub(crate) barcodes_skipped: u64,
}

/// Machine-readable run accounting, following the main `sbx10x` app.
#[derive(Debug, Serialize)]
struct Metrics {
    input_bam: String,
    output_dir: String,
    mex_dir: String,
    h5_file: String,
    elapsed_secs: f64,
    records_read: u64,
    records_per_sec: u64,
    molecules_counted: u64,
    features: usize,
    barcodes: usize,
    nonzero_entries: usize,
    used_feature_template: bool,
    used_barcode_template: bool,
    used_cell_barcode_whitelist: bool,
    gem_group: Option<u32>,
    skip_unlisted_barcodes: bool,
    barcodes_skipped: u64,
}

/// Run the BAM-to-matrix application.
pub fn run(cli: &Cli) -> Result<()> {
    let started = Instant::now();
    if cli.gem_group == 0 {
        bail!("--gem-group must be greater than zero");
    }
    let matrix = read_matrix(
        &cli.bam,
        cli.features.as_deref(),
        cli.barcodes.as_deref(),
        cli.cell_barcodes.as_deref(),
        cli.gem_group,
        cli.skip_unlisted_barcodes,
    )?;
    fs::create_dir_all(&cli.output)
        .with_context(|| format!("creating {}", cli.output.display()))?;
    let mex_dir = cli.output.join("raw_feature_bc_matrix");
    let h5_file = cli.output.join("raw_feature_bc_matrix.h5");
    write_matrix(&matrix, &mex_dir, !cli.no_gzip)?;
    write_h5(&matrix, &h5_file, &cli.genome)?;

    let elapsed = started.elapsed().as_secs_f64();
    let metrics = Metrics {
        input_bam: cli.bam.display().to_string(),
        output_dir: cli.output.display().to_string(),
        mex_dir: mex_dir.display().to_string(),
        h5_file: h5_file.display().to_string(),
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
        used_cell_barcode_whitelist: cli.cell_barcodes.is_some(),
        gem_group: cli.cell_barcodes.as_ref().map(|_| cli.gem_group),
        skip_unlisted_barcodes: cli.skip_unlisted_barcodes,
        barcodes_skipped: matrix.barcodes_skipped,
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

fn read_cell_barcodes(path: &Path, gem_group: u32) -> Result<Vec<String>> {
    text_reader(path)?
        .lines()
        .enumerate()
        .map(|(index, line)| {
            let barcode = line?;
            let barcode = barcode.trim();
            if barcode.is_empty() {
                bail!("empty barcode at {}:{}", path.display(), index + 1);
            }
            if !barcode
                .bytes()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            {
                bail!(
                    "invalid 10x reference barcode at {}:{}: {barcode}",
                    path.display(),
                    index + 1
                );
            }
            Ok(format!("{barcode}-{gem_group}"))
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn read_matrix(
    path: &Path,
    feature_template: Option<&Path>,
    barcode_template: Option<&Path>,
    cell_barcode_whitelist: Option<&Path>,
    gem_group: u32,
    skip_unlisted_barcodes: bool,
) -> Result<Matrix> {
    let mut reader =
        bam::Reader::from_path(path).with_context(|| format!("opening BAM {}", path.display()))?;
    let strict_features = feature_template.is_some();
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
        } else if strict_features {
            bail!("counted gene {gx} is absent from the canonical feature template")
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
    } else if let Some(path) = cell_barcode_whitelist {
        read_cell_barcodes(path, gem_group)?
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
    let mut barcodes_skipped = 0u64;
    for (feature, barcode) in &molecule_keys {
        let Some(&barcode) = barcode_indexes.get(barcode.as_str()) else {
            if skip_unlisted_barcodes {
                barcodes_skipped += 1;
                continue;
            }
            bail!("counted barcode {barcode} is absent from barcode template");
        };
        *counts.entry((barcode, *feature)).or_insert(0) += 1;
    }

    Ok(Matrix {
        features: feature_order,
        barcodes,
        counts,
        records,
        molecules: molecule_keys.len() as u64 - barcodes_skipped,
        barcodes_skipped,
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

pub(crate) fn write_matrix(matrix: &Matrix, output: &Path, gzip: bool) -> Result<()> {
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
    for (&(barcode, feature), &count) in &matrix.counts {
        writeln!(mtx, "{} {} {}", feature + 1, barcode + 1, count)?;
    }
    Ok(())
}

fn ascii(value: &str, field: &str) -> Result<VarLenAscii> {
    VarLenAscii::from_ascii(value.as_bytes())
        .with_context(|| format!("{field} contains a non-ASCII or NUL character: {value:?}"))
}

/// Write the standard 10x CSC HDF5 feature-barcode matrix alongside the MEX
/// directory. The barcode-major BTreeMap is already in CSC traversal order.
pub(crate) fn write_h5(matrix: &Matrix, path: &Path, genome: &str) -> Result<()> {
    let file = hdf5::File::create(path).with_context(|| format!("creating {}", path.display()))?;
    let group = file.create_group("matrix")?;
    let features = group.create_group("features")?;

    let barcodes: Vec<_> = matrix
        .barcodes
        .iter()
        .map(|value| ascii(value, "barcode"))
        .collect::<Result<_>>()?;
    let ids: Vec<_> = matrix
        .features
        .iter()
        .map(|feature| ascii(&feature.id, "feature ID"))
        .collect::<Result<_>>()?;
    let names: Vec<_> = matrix
        .features
        .iter()
        .map(|feature| ascii(&feature.name, "feature name"))
        .collect::<Result<_>>()?;
    let gene_expression = ascii("Gene Expression", "feature type")?;
    let feature_types = vec![gene_expression; matrix.features.len()];
    let genome = ascii(genome, "genome")?;
    let genomes = vec![genome; matrix.features.len()];
    let all_tag_keys = vec![ascii("genome", "tag key")?];

    let mut data = Vec::<u32>::with_capacity(matrix.counts.len());
    let mut indices = Vec::<u32>::with_capacity(matrix.counts.len());
    let mut indptr = Vec::<u32>::with_capacity(matrix.barcodes.len() + 1);
    indptr.push(0);
    let mut entries = matrix.counts.iter().peekable();
    for barcode in 0..matrix.barcodes.len() {
        while let Some((&(entry_barcode, feature), &count)) = entries.peek() {
            if entry_barcode != barcode {
                break;
            }
            data.push(u32::try_from(count).context("UMI count exceeds HDF5 uint32 range")?);
            indices
                .push(u32::try_from(feature).context("feature index exceeds HDF5 uint32 range")?);
            entries.next();
        }
        indptr
            .push(u32::try_from(data.len()).context("matrix has more than 2^32 nonzero entries")?);
    }
    if entries.next().is_some() {
        bail!("matrix contains a barcode index outside the barcode table");
    }
    let shape = [matrix.features.len() as u64, matrix.barcodes.len() as u64];

    group
        .new_dataset_builder()
        .with_data(&barcodes)
        .create("barcodes")?;
    group
        .new_dataset_builder()
        .with_data(&data)
        .deflate(1)
        .create("data")?;
    group
        .new_dataset_builder()
        .with_data(&indices)
        .deflate(1)
        .create("indices")?;
    group
        .new_dataset_builder()
        .with_data(&indptr)
        .deflate(1)
        .create("indptr")?;
    group
        .new_dataset_builder()
        .with_data(&shape)
        .create("shape")?;
    features
        .new_dataset_builder()
        .with_data(&all_tag_keys)
        .create("_all_tag_keys")?;
    features
        .new_dataset_builder()
        .with_data(&feature_types)
        .create("feature_type")?;
    features
        .new_dataset_builder()
        .with_data(&genomes)
        .create("genome")?;
    features
        .new_dataset_builder()
        .with_data(&ids)
        .create("id")?;
    features
        .new_dataset_builder()
        .with_data(&names)
        .create("name")?;

    file.new_attr::<VarLenAscii>()
        .shape(())
        .create("filetype")?
        .write_scalar(&ascii("matrix", "filetype")?)?;
    file.new_attr::<u32>()
        .shape(())
        .create("version")?
        .write_scalar(&2)?;
    file.new_attr::<VarLenAscii>()
        .shape(())
        .create("software_version")?
        .write_scalar(&ascii(
            &format!("cr-bam-to-mtx-{}", env!("CARGO_PKG_VERSION")),
            "software version",
        )?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tenx_whitelist_gets_gem_group_suffix() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("barcodes.txt");
        fs::write(&path, "AACCGGTT\nTTGGCCAA\n").unwrap();
        assert_eq!(
            read_cell_barcodes(&path, 2).unwrap(),
            ["AACCGGTT-2", "TTGGCCAA-2"]
        );
    }

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
            barcodes_skipped: 0,
        };
        let dir = tempfile::tempdir().unwrap();
        write_matrix(&matrix, dir.path(), false).unwrap();
        let text = fs::read_to_string(dir.path().join("matrix.mtx")).unwrap();
        assert!(text.ends_with("1 1 1\n1 1 2\n"));
    }

    #[test]
    fn writes_compatibility_barcode_major_order() {
        let matrix = Matrix {
            features: vec![
                Feature {
                    id: "g1".into(),
                    name: "G1".into(),
                },
                Feature {
                    id: "g2".into(),
                    name: "G2".into(),
                },
            ],
            barcodes: vec!["bc1-1".into(), "bc2-1".into()],
            // Keys are (barcode, feature), deliberately inserted out of order.
            counts: [((1, 0), 3), ((0, 1), 2), ((0, 0), 1)]
                .into_iter()
                .collect(),
            records: 6,
            molecules: 6,
            barcodes_skipped: 0,
        };
        let dir = tempfile::tempdir().unwrap();
        write_matrix(&matrix, dir.path(), false).unwrap();
        let text = fs::read_to_string(dir.path().join("matrix.mtx")).unwrap();
        let data: Vec<_> = text.lines().filter(|line| !line.starts_with('%')).collect();
        assert_eq!(data, ["2 2 3", "1 1 1", "2 1 2", "1 2 3"]);
    }

    #[test]
    fn writes_compatibility_csc_h5() {
        let matrix = Matrix {
            features: vec![
                Feature {
                    id: "g1".into(),
                    name: "G1".into(),
                },
                Feature {
                    id: "g2".into(),
                    name: "G2".into(),
                },
            ],
            barcodes: vec!["bc1-1".into(), "bc2-1".into(), "bc3-1".into()],
            counts: [((0, 1), 2), ((0, 0), 1), ((2, 0), 3)]
                .into_iter()
                .collect(),
            records: 6,
            molecules: 6,
            barcodes_skipped: 0,
        };
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("raw_feature_bc_matrix.h5");
        write_h5(&matrix, &path, "GRCh38").unwrap();

        let file = hdf5::File::open(path).unwrap();
        let group = file.group("matrix").unwrap();
        assert_eq!(
            group.dataset("shape").unwrap().read_raw::<u64>().unwrap(),
            [2, 3]
        );
        assert_eq!(
            group.dataset("data").unwrap().read_raw::<u32>().unwrap(),
            [1, 2, 3]
        );
        assert_eq!(
            group.dataset("indices").unwrap().read_raw::<u32>().unwrap(),
            [0, 1, 0]
        );
        assert_eq!(
            group.dataset("indptr").unwrap().read_raw::<u32>().unwrap(),
            [0, 2, 2, 3]
        );
        let features = group.group("features").unwrap();
        let genomes = features
            .dataset("genome")
            .unwrap()
            .read_raw::<VarLenAscii>()
            .unwrap();
        assert_eq!(
            genomes.iter().map(|x| x.as_str()).collect::<Vec<_>>(),
            ["GRCh38", "GRCh38"]
        );
    }

    /// A tiny synthetic BAM with two distinct cell barcodes, each with one
    /// representative (xf bit 8 set) molecule on GENE1.
    fn two_barcode_bam(dir: &Path) -> PathBuf {
        let path = dir.join("in.bam");
        let mut header = bam::Header::new();
        header.push_record(
            bam::header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 1000),
        );
        let mut writer = bam::Writer::from_path(&path, &header, bam::Format::Bam).unwrap();
        let cigar = bam::record::CigarString(vec![bam::record::Cigar::Match(10)]);
        for (i, cb) in ["AAAAAAAAAAAAAAAA-1", "CCCCCCCCCCCCCCCC-1"].into_iter().enumerate() {
            let mut r = bam::Record::new();
            r.set(b"read", Some(&cigar), b"ACGTACGTAC", &[30; 10]);
            r.set_tid(0);
            r.set_pos(i as i64);
            r.push_aux(b"CB", Aux::String(cb)).unwrap();
            r.push_aux(b"UB", Aux::String("UMI0000AAAA")).unwrap();
            r.push_aux(b"GX", Aux::String("GENE1")).unwrap();
            r.push_aux(b"xf", Aux::I32(25)).unwrap();
            writer.write(&r).unwrap();
        }
        drop(writer);
        path
    }

    #[test]
    fn unlisted_barcode_errors_by_default() {
        let dir = tempfile::tempdir().unwrap();
        let bam = two_barcode_bam(dir.path());
        let barcodes_file = dir.path().join("only_one.txt");
        fs::write(&barcodes_file, "AAAAAAAAAAAAAAAA\n").unwrap();

        let err = read_matrix(&bam, None, None, Some(&barcodes_file), 1, false).unwrap_err();
        assert!(err.to_string().contains("absent from barcode template"));
    }

    #[test]
    fn skip_unlisted_barcodes_drops_them_instead_of_erroring() {
        let dir = tempfile::tempdir().unwrap();
        let bam = two_barcode_bam(dir.path());
        let barcodes_file = dir.path().join("only_one.txt");
        fs::write(&barcodes_file, "AAAAAAAAAAAAAAAA\n").unwrap();

        let matrix = read_matrix(&bam, None, None, Some(&barcodes_file), 1, true).unwrap();
        assert_eq!(matrix.barcodes, ["AAAAAAAAAAAAAAAA-1"]);
        assert_eq!(matrix.barcodes_skipped, 1);
        assert_eq!(matrix.molecules, 1);
        assert_eq!(matrix.counts.len(), 1);
    }
}
