use std::io::Write;

use crate::yc::YcTag;

use super::types::canonical_id;

pub fn write_failed_record(
    out: &mut Vec<u8>,
    reason: &str,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
) -> std::io::Result<()> {
    out.write_all(b"@")?;
    out.write_all(id)?;
    write!(out, "\tXF:Z:{reason}\n")?;
    out.write_all(seq)?;
    out.write_all(b"\n+\n")?;
    out.write_all(qual)?;
    out.write_all(b"\n")
}

#[allow(clippy::too_many_arguments)]
pub fn write_record(
    output: &mut dyn Write,
    id: &[u8],
    raw_cb: &[u8],
    raw_umi: &[u8],
    cb_qual: &[u8],
    umi_qual: &[u8],
    corrected: &[u8],
    gem_group: u32,
    seq: &[u8],
    qual: &[u8],
    yc: Option<&YcTag>,
) -> std::io::Result<()> {
    output.write_all(b"@")?;
    output.write_all(canonical_id(id))?;
    output.write_all(b"\tCR:Z:")?;
    output.write_all(raw_cb)?;
    output.write_all(b"\tCB:Z:")?;
    output.write_all(corrected)?;
    write!(output, "-{gem_group}")?;
    output.write_all(b"\tUR:Z:")?;
    output.write_all(raw_umi)?;
    output.write_all(b"\tUB:Z:")?;
    output.write_all(raw_umi)?;
    output.write_all(b"\tCY:Z:")?;
    output.write_all(cb_qual)?;
    output.write_all(b"\tUY:Z:")?;
    output.write_all(umi_qual)?;
    if let Some(tag) = yc {
        write!(output, "\tYC:Z:{tag}")?;
    }
    output.write_all(b"\n")?;
    output.write_all(seq)?;
    output.write_all(b"\n+\n")?;
    output.write_all(qual)?;
    output.write_all(b"\n")
}
