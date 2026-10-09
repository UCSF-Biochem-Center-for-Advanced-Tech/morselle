//! Molecule counters use the evidence pass's CB/PB identity for later UMI groups.
//!
//! Rows are keyed by the packed `(cell, sample)` identity that molecule
//! grouping already uses, with the CB/PB strings stored once per row, so the
//! per-read cost is one lookup in a map the evidence pass builds anyway (a
//! string-keyed table added ~5% to `molecules` wall time on DZEC05).
use super::evidence::RegionEvidence;
use super::types::{Dna, Metrics, RawKey};
use crate::barcode_metrics::BarcodeTable;
use rustc_hash::FxHashMap;

pub const COLUMNS: [&str; 5] = [
    "eligible_reads",
    "raw_cell_gene_umis",
    "molecules_counted",
    "collision_groups",
    "representatives_tagged",
];

#[derive(Clone, Debug)]
pub struct BarcodeRow {
    pub barcode: String,
    pub sample: String,
    pub counts: [u64; 5],
}

pub type BarcodeRows = FxHashMap<(Dna, u32), BarcodeRow>;

pub(super) fn row<'a>(
    rows: &'a mut BarcodeRows,
    cell: Dna,
    sample: u32,
    cb: &str,
    pb: Option<&str>,
) -> &'a mut [u64; 5] {
    &mut rows
        .entry((cell, sample))
        .or_insert_with(|| BarcodeRow {
            barcode: cb.to_owned(),
            sample: pb.unwrap_or(".").to_owned(),
            counts: [0; 5],
        })
        .counts
}

pub(super) fn eligible(metrics: &mut Metrics, key: RawKey, cb: &str, pb: Option<&str>) {
    row(&mut metrics.barcode_rows, key.cell, key.sample, cb, pb)[0] += 1;
}

pub(super) fn bump(metrics: &mut Metrics, cell: Dna, sample: u32, column: usize) {
    if let Some(row) = metrics.barcode_rows.get_mut(&(cell, sample)) {
        row.counts[column] += 1;
    } else {
        // Direct decision tests construct evidence without an input BAM.
        let name = format!("{}-1", cell.unpack());
        row(&mut metrics.barcode_rows, cell, sample, &name, None)[column] += 1;
    }
}

/// Adds `source` rows into `target`; the first-seen CB/PB names win.
pub(super) fn merge_rows(
    target: &mut BarcodeRows,
    source: &BarcodeRows,
    remap: impl Fn(u32) -> u32,
) {
    for (&(cell, sample), src) in source {
        let dst = target
            .entry((cell, remap(sample)))
            .or_insert_with(|| BarcodeRow {
                counts: [0; 5],
                ..src.clone()
            });
        for (d, s) in dst.counts.iter_mut().zip(src.counts) {
            *d += s;
        }
    }
}

pub(super) fn merge_region(metrics: &mut Metrics, part: &RegionEvidence, sample_remap: &[u32]) {
    merge_rows(&mut metrics.barcode_rows, &part.barcode_rows, |sample| {
        if sample == super::types::NO_SAMPLE {
            sample
        } else {
            sample_remap[sample as usize]
        }
    });
}

/// The CB/PB-keyed table written with `--barcode-metrics`.
pub fn table(rows: &BarcodeRows, enabled: bool) -> BarcodeTable<5> {
    let mut table = BarcodeTable::new(enabled);
    if enabled {
        for row in rows.values() {
            let dst = table.row(&row.barcode, &row.sample);
            for (d, s) in dst.iter_mut().zip(row.counts) {
                *d += s;
            }
        }
    }
    table
}
