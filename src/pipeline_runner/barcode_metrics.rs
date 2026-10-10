//! Per-read accounting keyed by the emitted corrected CB and routing sample.
use crate::barcode_metrics::{SharedBarcodeTable, NO_BARCODE, UNRESOLVED_BARCODE};
use crate::pipeline::{CbStatus, ReadClass};

pub(super) enum BarcodeId {
    Special(&'static str),
    Stack { bytes: [u8; 64], len: usize },
    Long(String),
}

impl BarcodeId {
    fn from_bases(bases: &[u8]) -> Self {
        if bases.len() + 2 > 64 || std::str::from_utf8(bases).is_err() {
            return Self::Long(format!("{}-1", String::from_utf8_lossy(bases)));
        }
        let mut bytes = [0; 64];
        bytes[..bases.len()].copy_from_slice(bases);
        bytes[bases.len()..bases.len() + 2].copy_from_slice(b"-1");
        Self::Stack {
            bytes,
            len: bases.len() + 2,
        }
    }

    pub(super) fn as_str(&self) -> &str {
        match self {
            Self::Special(value) => value,
            Self::Stack { bytes, len } => std::str::from_utf8(&bytes[..*len]).unwrap(),
            Self::Long(value) => value,
        }
    }
}

pub(super) fn observe(
    table: Option<&SharedBarcodeTable<10>>,
    seq: &[u8],
    r1_start: Option<usize>,
    corrected_r1: Option<&[u8]>,
    cb_len: usize,
    cb_status: Option<CbStatus>,
    cb_gate: bool,
    class: ReadClass,
    sample: &str,
) -> BarcodeId {
    let barcode = if let Some(corrected) = corrected_r1 {
        BarcodeId::from_bases(&corrected[..cb_len])
    } else if let Some(start) = r1_start {
        if cb_gate && matches!(cb_status, Some(CbStatus::Fail)) {
            BarcodeId::Special(UNRESOLVED_BARCODE)
        } else {
            BarcodeId::from_bases(&seq[start..start + cb_len])
        }
    } else {
        BarcodeId::Special(NO_BARCODE)
    };
    let Some(table) = table else {
        return barcode;
    };
    table.update(barcode.as_str(), sample, |row| {
        row[0] += 1;
        row[1] += u64::from(matches!(cb_status, Some(CbStatus::OnWhitelist)));
        row[2] += u64::from(matches!(cb_status, Some(CbStatus::Within1)));
        match class {
            ReadClass::Full => row[4] += 1,
            ReadClass::Partial => row[5] += 1,
            ReadClass::Unparseable => row[6] += 1,
        }
    });
    barcode
}

pub(super) fn bump(
    table: Option<&SharedBarcodeTable<10>>,
    barcode: &BarcodeId,
    sample: &str,
    column: usize,
) {
    if let Some(table) = table {
        table.update(barcode.as_str(), sample, |row| row[column] += 1);
    }
}
