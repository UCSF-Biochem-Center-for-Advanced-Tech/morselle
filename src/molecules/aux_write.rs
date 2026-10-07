//! Low-level aux-tag edits for the molecules write pass.
//!
//! `Record::push_aux` first calls `self.aux(tag)` -- a full linear scan of
//! the record's ~24 aux tags -- to refuse duplicates, and `remove_aux`
//! allocates a `CString` per call. On a 2M-read benchmark those scans were
//! ~17% of all instructions in `sbx10x-molecules`. The write pass already
//! knows (from one `MoleculeTags` pass, plus the removals it just did) that
//! the tags it appends are absent, so these helpers call htslib directly.
//! The resulting bytes are identical to `remove_aux` + `push_aux`: removal
//! uses `bam_aux_get`/`bam_aux_del`, appends use `bam_aux_append` with the
//! same type codes and lengths rust-htslib uses.
use std::os::raw::c_char;

use anyhow::{bail, Result};
use rust_htslib::{bam, htslib};

/// Remove `tag` if present; returns whether it was. One aux scan.
pub fn remove_tag(record: &mut bam::Record, tag: &[u8; 2]) -> bool {
    let b = record.inner_mut() as *mut htslib::bam1_t;
    // SAFETY: `b` is a valid bam1_t owned by `record`; htslib reads exactly
    // two bytes of `tag`; `s` points into this record's own aux data.
    unsafe {
        let s = htslib::bam_aux_get(b, tag.as_ptr() as *const c_char);
        if s.is_null() {
            return false;
        }
        htslib::bam_aux_del(b, s) == 0
    }
}

fn append(record: &mut bam::Record, tag: &[u8; 2], type_code: u8, bytes: &[u8]) -> Result<()> {
    debug_assert!(record.aux(tag).is_err(), "append_* requires an absent tag");
    let b = record.inner_mut() as *mut htslib::bam1_t;
    // SAFETY: valid bam1_t; htslib copies `bytes.len()` bytes from `bytes`.
    let ret = unsafe {
        htslib::bam_aux_append(
            b,
            tag.as_ptr() as *const c_char,
            type_code as c_char,
            bytes.len() as i32,
            bytes.as_ptr(),
        )
    };
    if ret < 0 {
        bail!("failed to append aux tag {}", String::from_utf8_lossy(tag));
    }
    Ok(())
}

/// Append a `Z` tag the caller knows is absent (NUL-terminated, as htslib).
pub fn append_str(record: &mut bam::Record, tag: &[u8; 2], value: &str) -> Result<()> {
    if value.as_bytes().contains(&0) {
        bail!(
            "aux string for {} contains NUL",
            String::from_utf8_lossy(tag)
        );
    }
    let mut bytes = Vec::with_capacity(value.len() + 1);
    bytes.extend_from_slice(value.as_bytes());
    bytes.push(0);
    append(record, tag, b'Z', &bytes)
}

pub fn append_i32(record: &mut bam::Record, tag: &[u8; 2], value: i32) -> Result<()> {
    append(record, tag, b'i', &value.to_le_bytes())
}

pub fn append_u32(record: &mut bam::Record, tag: &[u8; 2], value: u32) -> Result<()> {
    append(record, tag, b'I', &value.to_le_bytes())
}

pub fn append_char(record: &mut bam::Record, tag: &[u8; 2], value: u8) -> Result<()> {
    append(record, tag, b'A', &[value])
}

pub fn append_f32(record: &mut bam::Record, tag: &[u8; 2], value: f32) -> Result<()> {
    append(record, tag, b'f', &value.to_le_bytes())
}

/// `remove_aux` + `push_aux(String)`: the tag ends up last, exactly as before.
pub fn replace_str(record: &mut bam::Record, tag: &[u8; 2], value: &str) -> Result<()> {
    remove_tag(record, tag);
    append_str(record, tag, value)
}
