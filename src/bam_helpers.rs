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

/// The record's variable-length data block (qname, cigar, seq, qual, aux).
/// rust-htslib keeps `Record::data` private; this is the same slice.
pub fn raw_data(record: &bam::Record) -> &[u8] {
    let inner = record.inner();
    if inner.data.is_null() || inner.l_data <= 0 {
        return &[];
    }
    // SAFETY: `data` points to `l_data` initialized bytes owned by `record`,
    // valid for the lifetime of the borrow.
    unsafe { std::slice::from_raw_parts(inner.data, inner.l_data as usize) }
}

/// Byte length of one aux value of `type_code` at the start of `body`
/// (SAM spec / htslib `skip_aux`), or `None` if malformed or truncated.
fn aux_value_len(type_code: u8, body: &[u8]) -> Option<usize> {
    let len = match type_code {
        b'A' | b'c' | b'C' => 1,
        b's' | b'S' => 2,
        b'i' | b'I' | b'f' => 4,
        b'd' => 8,
        b'Z' | b'H' => body.iter().position(|&b| b == 0)? + 1,
        b'B' => {
            let header = body.get(..5)?;
            let count = u32::from_le_bytes([header[1], header[2], header[3], header[4]]) as usize;
            let size = match header[0] {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return None,
            };
            5usize.checked_add(count.checked_mul(size)?)?
        }
        _ => return None,
    };
    (len <= body.len()).then_some(len)
}

/// Walk the raw aux fields once, calling `visit(tag, type_code, value)` for
/// each (`value` excludes the tag and type bytes; `Z`/`H` include the NUL).
/// Fields are skipped by their encoded length, like htslib's `bam_aux_get`,
/// so nothing is decoded or UTF-8-checked unless the caller does it. Stops
/// at malformed data. Much cheaper than `aux_iter()` or repeated `aux()`.
pub fn for_each_aux<'a>(record: &'a bam::Record, mut visit: impl FnMut([u8; 2], u8, &'a [u8])) {
    let core = &record.inner().core;
    let start = core.l_qname as usize
        + core.n_cigar as usize * 4
        + (core.l_qseq as usize).div_ceil(2)
        + core.l_qseq as usize;
    let mut aux = raw_data(record).get(start..).unwrap_or(&[]);
    while aux.len() >= 3 {
        let (tag, type_code, body) = ([aux[0], aux[1]], aux[2], &aux[3..]);
        let Some(len) = aux_value_len(type_code, body) else {
            break;
        };
        visit(tag, type_code, &body[..len]);
        aux = &body[len..];
    }
}

/// A decoded aux value, as far as the typed helpers below care: what
/// rust-htslib's `record.aux(tag)` would return, collapsed to the classes
/// `aux_string` / `aux_i64` / `aux_f64` distinguish.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AuxValue<'a> {
    /// `Z` or `H` (rust-htslib decodes both as `Aux::String`).
    Str(&'a str),
    /// `c C s S i I`.
    Int(i64),
    /// `f`.
    Float(f32),
    /// Any other type (`A`, `d`, `B` arrays): present but not one of the above.
    Other,
}

fn decode(type_code: u8, value: &[u8]) -> Option<AuxValue<'_>> {
    Some(match type_code {
        // `record.aux` fails on invalid UTF-8; the helpers treat that as absent.
        b'Z' | b'H' => AuxValue::Str(std::str::from_utf8(&value[..value.len() - 1]).ok()?),
        b'c' => AuxValue::Int(value[0] as i8 as i64),
        b'C' => AuxValue::Int(value[0] as i64),
        b's' => AuxValue::Int(i16::from_le_bytes([value[0], value[1]]) as i64),
        b'S' => AuxValue::Int(u16::from_le_bytes([value[0], value[1]]) as i64),
        b'i' => AuxValue::Int(i32::from_le_bytes([value[0], value[1], value[2], value[3]]) as i64),
        b'I' => AuxValue::Int(u32::from_le_bytes([value[0], value[1], value[2], value[3]]) as i64),
        b'f' => AuxValue::Float(f32::from_le_bytes([value[0], value[1], value[2], value[3]])),
        _ => AuxValue::Other,
    })
}

/// First occurrence of each of `tags`, gathered in ONE raw pass (see
/// [`for_each_aux`]). Feed the slots to [`aux_slot_string`] /
/// [`aux_slot_i64`] / [`aux_slot_f64`], which apply exactly the typing
/// rules and error messages of [`aux_string`] / [`aux_i64`] / [`aux_f64`].
pub fn aux_values<'a, const N: usize>(
    record: &'a bam::Record,
    tags: [&[u8; 2]; N],
) -> [Option<AuxValue<'a>>; N] {
    let mut out = [None; N];
    let mut seen = [false; N];
    for_each_aux(record, |tag, type_code, value| {
        if let Some(slot) = tags.iter().position(|t| **t == tag) {
            if !seen[slot] {
                seen[slot] = true;
                out[slot] = decode(type_code, value);
            }
        }
    });
    out
}

/// [`aux_string`] semantics for a slot from [`aux_values`], borrowed.
pub fn aux_slot_string<'a>(slot: &Option<AuxValue<'a>>, tag: &[u8]) -> Result<Option<&'a str>> {
    match slot {
        Some(AuxValue::Str(value)) => Ok(Some(value)),
        Some(_) => bail!("tag {} is not a string", String::from_utf8_lossy(tag)),
        None => Ok(None),
    }
}

/// [`aux_i64`] semantics for a slot from [`aux_values`].
pub fn aux_slot_i64(slot: &Option<AuxValue>, tag: &[u8]) -> Result<Option<i64>> {
    match slot {
        Some(AuxValue::Int(value)) => Ok(Some(*value)),
        Some(_) => bail!("tag {} is not an integer", String::from_utf8_lossy(tag)),
        None => Ok(None),
    }
}

/// [`aux_f64`] semantics for a slot from [`aux_values`].
pub fn aux_slot_f64(slot: &Option<AuxValue>, tag: &[u8]) -> Result<Option<f64>> {
    match slot {
        Some(AuxValue::Float(value)) => Ok(Some(*value as f64)),
        Some(_) => bail!("tag {} is not a float", String::from_utf8_lossy(tag)),
        None => Ok(None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn aux_values_matches_per_tag_helpers() {
        let mut record = bam::Record::new();
        record.set(b"q", None, b"ACGT", &[30; 4]);
        record.push_aux(b"CB", Aux::String("AAAC-1")).unwrap();
        record.push_aux(b"xf", Aux::U8(25)).unwrap();
        record.push_aux(b"GX", Aux::I32(7)).unwrap();
        record.push_aux(b"ZW", Aux::Float(0.25)).unwrap();
        let [cb, xf, gx, zw, missing] = aux_values(&record, [b"CB", b"xf", b"GX", b"ZW", b"PB"]);
        assert_eq!(aux_slot_string(&cb, b"CB").unwrap(), Some("AAAC-1"));
        assert_eq!(
            aux_slot_i64(&xf, b"xf").unwrap(),
            aux_i64(&record, b"xf").unwrap()
        );
        assert_eq!(
            aux_slot_f64(&zw, b"ZW").unwrap(),
            aux_f64(&record, b"ZW").unwrap()
        );
        assert_eq!(aux_slot_string(&missing, b"PB").unwrap(), None);
        assert_eq!(
            aux_slot_string(&gx, b"GX").unwrap_err().to_string(),
            aux_string(&record, b"GX").unwrap_err().to_string()
        );
        assert_eq!(
            aux_slot_i64(&cb, b"CB").unwrap_err().to_string(),
            aux_i64(&record, b"CB").unwrap_err().to_string()
        );
    }

    /// Randomized: every slot reader must return exactly what the per-tag
    /// helper returns (value, None, or the same error) for every aux type.
    #[test]
    fn aux_values_matches_per_tag_helpers_on_random_records() {
        use rust_htslib::bam::record::AuxArray;
        let pool: [&[u8; 2]; 5] = [b"xf", b"CB", b"ZW", b"GX", b"NM"];
        let mut state: u64 = 0x2545_F491_4F6C_DD1D;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };
        for case in 0..3_000 {
            let mut record = bam::Record::new();
            record.set(b"q", None, b"ACG", &[30; 3]);
            for _ in 0..(next() % 6) {
                let tag = pool[(next() % pool.len() as u64) as usize];
                if record.aux(tag).is_ok() {
                    continue;
                }
                let arr = [1u8, 2];
                let value = match next() % 11 {
                    0 => Aux::String("AAAC-1"),
                    1 => Aux::I8(-3),
                    2 => Aux::U8(25),
                    3 => Aux::I16(-300),
                    4 => Aux::U16(60000),
                    5 => Aux::I32(-7),
                    6 => Aux::U32(4_000_000_000),
                    7 => Aux::Float(0.25),
                    8 => Aux::Double(0.5),
                    9 => Aux::Char(b'F'),
                    _ => Aux::ArrayU8(AuxArray::from(&arr[..])),
                };
                record.push_aux(tag, value).unwrap();
            }
            let slots = aux_values(&record, pool);
            for (tag, slot) in pool.iter().zip(&slots) {
                let s = |r: Result<Option<String>>| r.map_err(|e| e.to_string());
                assert_eq!(
                    s(aux_slot_string(slot, *tag).map(|v| v.map(str::to_owned))),
                    s(aux_string(&record, *tag)),
                    "case {case}"
                );
                let e = |r: Result<Option<i64>>| r.map_err(|e| e.to_string());
                assert_eq!(
                    e(aux_slot_i64(slot, *tag)),
                    e(aux_i64(&record, *tag)),
                    "case {case}"
                );
                let f = |r: Result<Option<f64>>| r.map_err(|e| e.to_string());
                assert_eq!(
                    f(aux_slot_f64(slot, *tag)),
                    f(aux_f64(&record, *tag)),
                    "case {case}"
                );
            }
        }
    }
}
