use std::collections::HashSet;
use std::io::BufRead;
use std::path::Path;

use anyhow::{bail, Context, Result};

use super::types::Feature;

pub use crate::bam_helpers::{aux_f64, aux_i64, aux_string};

pub fn text_reader(path: &Path) -> Result<Box<dyn BufRead>> {
    crate::io::open_text_reader(path)
}

pub(crate) fn read_features(path: &Path) -> Result<Vec<Feature>> {
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

pub fn read_barcodes(path: &Path) -> Result<Vec<String>> {
    Ok(text_reader(path)?
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?)
}

pub fn read_cell_barcodes(path: &Path, gem_group: u32) -> Result<Vec<String>> {
    text_reader(path)?
        .lines()
        .enumerate()
        .map(|(index, line)| {
            let barcode = line?;
            let barcode = barcode.trim();
            if barcode.is_empty() {
                bail!("empty barcode at {}:{}", path.display(), index + 1);
            }
            let bare = barcode.split('-').next().unwrap_or(barcode);
            if !bare
                .bytes()
                .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            {
                bail!(
                    "invalid 10x reference barcode at {}:{}: {barcode}",
                    path.display(),
                    index + 1
                );
            }
            Ok(format!("{bare}-{gem_group}"))
        })
        .collect()
}

pub fn read_probe_barcodes(path: &Path) -> Result<HashSet<String>> {
    let mut set = HashSet::new();
    for line in text_reader(path)?.lines() {
        let trimmed = line?.trim().to_owned();
        if !trimmed.is_empty() {
            set.insert(trimmed);
        }
    }
    Ok(set)
}
