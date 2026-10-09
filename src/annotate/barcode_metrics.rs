//! Fold the same annotation counters into each input CB's row.
use super::types::{Metrics, RecordMetrics};
use crate::barcode_metrics::{NO_BARCODE, UNRESOLVED_BARCODE};
use rust_htslib::bam::{record::Aux, Record};

pub(super) fn observe(metrics: &mut Metrics, record: &Record, value: &RecordMetrics) {
    if !metrics.per_barcode.is_enabled() {
        return;
    }
    let barcode = match record.aux(b"CB") {
        Ok(Aux::String(cb)) if !cb.is_empty() => {
            if cb.split('-').next().is_some_and(|v| {
                v.len() == 16 && v.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
            }) {
                cb
            } else {
                UNRESOLVED_BARCODE
            }
        }
        _ => NO_BARCODE,
    };
    let row = metrics.per_barcode.row(barcode, ".");
    for (dst, src) in row.iter_mut().zip([
        value.records,
        value.primary,
        value.unmapped,
        value.exonic,
        value.intronic,
        value.intergenic,
        value.antisense,
        value.gene_assigned,
        value.unique_gene,
        value.multi_gene,
        value.rescued,
    ]) {
        *dst += src;
    }
}
