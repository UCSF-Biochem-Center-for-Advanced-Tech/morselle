use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::Path;
use std::time::Instant;

use anyhow::{bail, Context, Result};
use flate2::{write::GzEncoder, Compression};
use hdf5::types::VarLenAscii;
use rust_htslib::bam::{self, Read};
use tracing::info;

use super::io_helpers::{read_barcodes, read_cell_barcodes, read_features, read_probe_barcodes};
use crate::bam_helpers::{aux_slot_f64, aux_slot_i64, aux_slot_string, aux_values};
use super::types::{Cli, Feature, Matrix, Metrics};

#[allow(clippy::too_many_arguments)]
pub(crate) fn read_matrix(
    path: &Path,
    feature_template: Option<&Path>,
    barcode_template: Option<&Path>,
    cell_barcode_whitelist: Option<&Path>,
    gem_group: u32,
    skip_unlisted_barcodes: bool,
    probe_barcode_allowlist: Option<&HashSet<String>>,
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
    // (feature index, interned CB id, weight). Interning keeps one String
    // per distinct barcode instead of one per counted molecule.
    let mut molecule_keys: Vec<(usize, u32, f64)> = Vec::new();
    let mut cb_names: Vec<String> = Vec::new();
    let mut cb_ids: HashMap<String, u32> = HashMap::new();
    let mut records = 0_u64;
    let mut probe_barcodes_skipped = 0_u64;

    for result in reader.records() {
        let record = result.context("reading BAM record")?;
        records += 1;
        // One aux scan per record; typing/errors identical to the old
        // per-tag aux_i64/aux_string/aux_f64 calls, and applied in the same
        // order (only xf is checked on non-representative records).
        let [xf, pb, cb, ub, gx, gn, zw] =
            aux_values(&record, [b"xf", b"PB", b"CB", b"UB", b"GX", b"GN", b"ZW"]);
        if aux_slot_i64(&xf, b"xf")?.unwrap_or(0) & 8 == 0 {
            continue;
        }
        if let Some(allowlist) = probe_barcode_allowlist {
            match aux_slot_string(&pb, b"PB")? {
                Some(pb) if allowlist.contains(pb) => {}
                _ => {
                    probe_barcodes_skipped += 1;
                    continue;
                }
            }
        }
        let cb = aux_slot_string(&cb, b"CB")?
            .context("xf bit 8 record is missing corrected barcode tag CB")?;
        let _ub = aux_slot_string(&ub, b"UB")?
            .context("xf bit 8 record is missing corrected UMI tag UB")?;
        let gx =
            aux_slot_string(&gx, b"GX")?.context("xf bit 8 record is missing gene ID tag GX")?;
        let gn = aux_slot_string(&gn, b"GN")?.unwrap_or(gx);
        if gx.contains(';') {
            bail!("xf bit 8 record has multiple genes in GX: {gx}");
        }
        let feature_index = if let Some(&index) = feature_indexes.get(gx) {
            index
        } else if strict_features {
            bail!("counted gene {gx} is absent from the canonical feature template")
        } else {
            let index = feature_order.len();
            feature_indexes.insert(gx.to_owned(), index);
            feature_order.push(Feature { id: gx.to_owned(), name: gn.to_owned() });
            index
        };
        let weight = aux_slot_f64(&zw, b"ZW")?.unwrap_or(1.0);
        let cb_id = match cb_ids.get(cb) {
            Some(&id) => id,
            None => {
                let id = cb_names.len() as u32;
                cb_names.push(cb.to_owned());
                cb_ids.insert(cb.to_owned(), id);
                id
            }
        };
        molecule_keys.push((feature_index, cb_id, weight));
    }
    drop(cb_ids);

    let barcodes = if let Some(path) = barcode_template {
        read_barcodes(path)?
    } else if let Some(path) = cell_barcode_whitelist {
        read_cell_barcodes(path, gem_group)?
    } else {
        let mut values = cb_names.clone();
        values.sort_unstable();
        values.dedup();
        values
    };
    let barcode_indexes: HashMap<&str, usize> = barcodes
        .iter()
        .enumerate()
        .map(|(index, barcode)| (barcode.as_str(), index))
        .collect();
    // Resolve each distinct CB once instead of once per molecule.
    let cb_to_column: Vec<Option<usize>> = cb_names
        .iter()
        .map(|cb| barcode_indexes.get(cb.as_str()).copied())
        .collect();
    let mut counts = BTreeMap::new();
    let mut barcodes_skipped = 0u64;
    for &(feature, cb_id, weight) in &molecule_keys {
        let Some(barcode) = cb_to_column[cb_id as usize] else {
            if skip_unlisted_barcodes {
                barcodes_skipped += 1;
                continue;
            }
            bail!(
                "counted barcode {} is absent from barcode template",
                cb_names[cb_id as usize]
            );
        };
        *counts.entry((barcode, feature)).or_insert(0.0) += weight;
    }

    Ok(Matrix {
        features: feature_order,
        barcodes,
        counts,
        records,
        molecules: molecule_keys.len() as u64 - barcodes_skipped,
        barcodes_skipped,
        probe_barcodes_skipped,
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
    let is_fractional = matrix.counts.values().any(|&c| (c - c.round()).abs() > 1e-5);
    let matrix_type = if is_fractional { "real" } else { "integer" };
    let mut mtx = writer(&output.join(format!("matrix.mtx{suffix}")), gzip)?;
    writeln!(mtx, "%%MatrixMarket matrix coordinate {matrix_type} general")?;
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
        if is_fractional {
            writeln!(mtx, "{} {} {:.4}", feature + 1, barcode + 1, count)?;
        } else {
            writeln!(mtx, "{} {} {}", feature + 1, barcode + 1, count.round() as u64)?;
        }
    }
    Ok(())
}

fn ascii(value: &str, field: &str) -> Result<VarLenAscii> {
    VarLenAscii::from_ascii(value.as_bytes())
        .with_context(|| format!("{field} contains a non-ASCII or NUL character: {value:?}"))
}

fn write_sparse_column<T: hdf5::H5Type>(group: &hdf5::Group, values: &[T], name: &str) -> Result<()> {
    if values.is_empty() {
        group.new_dataset_builder().with_data(values).create(name)?;
    } else {
        group
            .new_dataset_builder()
            .with_data(values)
            .chunk(values.len())
            .deflate(1)
            .create(name)?;
    }
    Ok(())
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

    let is_fractional = matrix.counts.values().any(|&c| (c - c.round()).abs() > 1e-5);
    let mut data_u32 = Vec::<u32>::new();
    let mut data_f32 = Vec::<f32>::new();
    let mut indices = Vec::<u32>::with_capacity(matrix.counts.len());
    let mut indptr = Vec::<u32>::with_capacity(matrix.barcodes.len() + 1);
    indptr.push(0);
    let mut entries = matrix.counts.iter().peekable();
    for barcode in 0..matrix.barcodes.len() {
        while let Some((&(entry_barcode, feature), &count)) = entries.peek() {
            if entry_barcode != barcode {
                break;
            }
            if is_fractional {
                data_f32.push(count as f32);
            } else {
                data_u32.push(u32::try_from(count.round() as u64).context("UMI count exceeds HDF5 uint32 range")?);
            }
            indices
                .push(u32::try_from(feature).context("feature index exceeds HDF5 uint32 range")?);
            entries.next();
        }
        let total_entries = if is_fractional { data_f32.len() } else { data_u32.len() };
        indptr
            .push(u32::try_from(total_entries).context("matrix has more than 2^32 nonzero entries")?);
    }
    if entries.next().is_some() {
        bail!("matrix contains a barcode index outside the barcode table");
    }
    let shape = [matrix.features.len() as u64, matrix.barcodes.len() as u64];

    group
        .new_dataset_builder()
        .with_data(&barcodes)
        .create("barcodes")?;
    if is_fractional {
        write_sparse_column(&group, &data_f32, "data")?;
    } else {
        write_sparse_column(&group, &data_u32, "data")?;
    }
    write_sparse_column(&group, &indices, "indices")?;
    write_sparse_column(&group, &indptr, "indptr")?;
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

/// Run the BAM-to-matrix application.
pub fn run(cli: &Cli) -> Result<()> {
    let started = Instant::now();
    if cli.gem_group == 0 {
        bail!("--gem-group must be greater than zero");
    }
    let probe_barcode_allowlist = cli
        .probe_barcodes
        .as_deref()
        .map(read_probe_barcodes)
        .transpose()?;
    let matrix = read_matrix(
        &cli.bam,
        cli.features.as_deref(),
        cli.barcodes.as_deref(),
        cli.cell_barcodes.as_deref(),
        cli.gem_group,
        cli.skip_unlisted_barcodes,
        probe_barcode_allowlist.as_ref(),
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
        used_probe_barcode_allowlist: cli.probe_barcodes.is_some(),
        probe_barcodes_skipped: matrix.probe_barcodes_skipped,
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
