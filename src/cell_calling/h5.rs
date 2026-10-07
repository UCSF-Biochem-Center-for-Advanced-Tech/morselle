use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};
use hdf5::types::VarLenAscii;

use crate::cr_bam_to_mtx::{Feature, Matrix};
use super::types::RawMatrix;

pub fn resolve_h5(path: &Path) -> PathBuf {
    if path.is_dir() {
        path.join("raw_feature_bc_matrix.h5")
    } else {
        path.to_owned()
    }
}

pub fn ascii_strings(dataset: hdf5::Dataset) -> Result<Vec<String>> {
    Ok(dataset
        .read_raw::<VarLenAscii>()?
        .into_iter()
        .map(|value| value.as_str().to_owned())
        .collect())
}

pub(crate) fn read_raw(path: &Path) -> Result<RawMatrix> {
    let file = hdf5::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    let matrix = file.group("matrix")?;
    let features_group = matrix.group("features")?;
    let ids = ascii_strings(features_group.dataset("id")?)?;
    let names = ascii_strings(features_group.dataset("name")?)?;
    let genomes = ascii_strings(features_group.dataset("genome")?)?;
    let barcodes = ascii_strings(matrix.dataset("barcodes")?)?;
    let data = matrix.dataset("data")?.read_raw::<u32>()?;
    let indices = matrix.dataset("indices")?.read_raw::<u32>()?;
    let indptr = matrix.dataset("indptr")?.read_raw::<u32>()?;
    if ids.len() != names.len() || data.len() != indices.len() || indptr.len() != barcodes.len() + 1
    {
        bail!("invalid 10x CSC matrix dimensions");
    }
    Ok(RawMatrix {
        features: ids
            .into_iter()
            .zip(names)
            .map(|(id, name)| Feature { id, name })
            .collect(),
        barcodes,
        genome: genomes.first().cloned().unwrap_or_else(|| "unknown".into()),
        data,
        indices,
        indptr,
    })
}

pub(crate) fn totals(raw: &RawMatrix) -> Vec<u64> {
    (0..raw.barcodes.len())
        .map(|bc| {
            let start = raw.indptr[bc] as usize;
            let end = raw.indptr[bc + 1] as usize;
            raw.data[start..end].iter().map(|&x| x as u64).sum()
        })
        .collect()
}

pub(crate) fn filtered_matrix(raw: &RawMatrix, selected: &[usize]) -> Matrix {
    let selected_set: BTreeSet<usize> = selected.iter().copied().collect();
    let ordered: Vec<usize> = (0..raw.barcodes.len())
        .filter(|i| selected_set.contains(i))
        .collect();
    let mut counts = BTreeMap::new();
    let mut molecules = 0u64;
    for (new_bc, &old_bc) in ordered.iter().enumerate() {
        for pos in raw.indptr[old_bc] as usize..raw.indptr[old_bc + 1] as usize {
            let count = raw.data[pos] as u64;
            counts.insert((new_bc, raw.indices[pos] as usize), count as f64);
            molecules += count;
        }
    }
    Matrix {
        features: raw.features.clone(),
        barcodes: ordered.iter().map(|&i| raw.barcodes[i].clone()).collect(),
        counts,
        records: 0,
        molecules,
        barcodes_skipped: 0,
        probe_barcodes_skipped: 0,
    }
}
