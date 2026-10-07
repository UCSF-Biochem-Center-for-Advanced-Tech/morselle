//! Shared helper functions for parsing BAM auxiliary tags.

use anyhow::{bail, Result};
use rust_htslib::bam::{self, record::Aux};

/// Extract an optional owned string aux tag from a BAM record.
pub fn aux_string(record: &bam::Record, tag: &[u8]) -> Result<Option<String>> {
    match record.aux(tag) {
        Ok(Aux::String(value)) => Ok(Some(value.to_owned())),
        Ok(_) => bail!("tag {} is not a string", String::from_utf8_lossy(tag)),
        Err(_) => Ok(None),
    }
}

/// Borrow an optional string slice aux tag from a BAM record.
pub fn aux_str<'a>(record: &'a bam::Record, tag: &[u8]) -> Option<&'a str> {
    match record.aux(tag) {
        Ok(Aux::String(s)) => Some(s),
        _ => None,
    }
}

/// Extract an optional integer aux tag from a BAM record.
pub fn aux_i64(record: &bam::Record, tag: &[u8]) -> Result<Option<i64>> {
    let value = match record.aux(tag) {
        Ok(Aux::I8(v)) => v as i64,
        Ok(Aux::U8(v)) => v as i64,
        Ok(Aux::I16(v)) => v as i64,
        Ok(Aux::U16(v)) => v as i64,
        Ok(Aux::I32(v)) => v as i64,
        Ok(Aux::U32(v)) => v as i64,
        Ok(_) => bail!("tag {} is not an integer", String::from_utf8_lossy(tag)),
        Err(_) => return Ok(None),
    };
    Ok(Some(value))
}

/// Borrow an optional integer aux tag value from a BAM record without throwing error.
pub fn aux_int(record: &bam::Record, tag: &[u8]) -> Option<i64> {
    aux_i64(record, tag).ok().flatten()
}

/// Extract an optional char aux tag from a BAM record.
pub fn aux_char(record: &bam::Record, tag: &[u8]) -> Option<u8> {
    match record.aux(tag) {
        Ok(Aux::Char(c)) => Some(c),
        _ => None,
    }
}

/// Extract an optional float aux tag from a BAM record.
pub fn aux_f64(record: &bam::Record, tag: &[u8]) -> Result<Option<f64>> {
    match record.aux(tag) {
        Ok(Aux::Float(v)) => Ok(Some(v as f64)),
        Ok(_) => bail!("tag {} is not a float", String::from_utf8_lossy(tag)),
        Err(_) => Ok(None),
    }
}

