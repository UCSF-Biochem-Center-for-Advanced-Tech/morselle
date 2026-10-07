use anyhow::Result;

use crate::anchors::comp;
use crate::pipeline_runner::types::ProcessingTags;
use crate::pipeline_runner::worker::Shared;

/// Flush a target's batch buffer once its raw R1 bytes reach this size.
pub(crate) const FLUSH_BYTES: usize = 256 * 1024;

/// Minimum dT-primer run length to count as a confident 3' flank for the CB.
pub(crate) const MIN_POLYT_FLANK: u16 = 10;

/// Append one FASTQ record (id without '@') to a raw byte buffer.
pub(crate) fn push_fastq(buf: &mut Vec<u8>, id: &[u8], seq: &[u8], qual: &[u8]) {
    buf.push(b'@');
    buf.extend_from_slice(id);
    buf.push(b'\n');
    buf.extend_from_slice(seq);
    buf.extend_from_slice(b"\n+\n");
    buf.extend_from_slice(qual);
    buf.push(b'\n');
}

/// Append one FASTQ record whose SEQ is the reverse-complement of `seq[start..end]`
/// and whose QUAL is the reverse of `qual[start..end]`, optionally capped to `max_len`.
pub(crate) fn push_fastq_r2_rc(
    buf: &mut Vec<u8>,
    id: &[u8],
    seq: &[u8],
    qual: &[u8],
    start: usize,
    end: usize,
    max_len: usize,
) -> usize {
    let full = end - start;
    let out = if max_len > 0 { full.min(max_len) } else { full };
    buf.reserve(id.len() + 2 * out + 6);
    buf.push(b'@');
    buf.extend_from_slice(id);
    buf.push(b'\n');
    for k in 0..out {
        buf.push(comp(seq[end - 1 - k]));
    }
    buf.extend_from_slice(b"\n+\n");
    for k in 0..out {
        buf.push(qual[end - 1 - k]);
    }
    buf.push(b'\n');
    out
}

/// Build a FASTQ comment that minimap2 `-y` copies verbatim to SAM/BAM as auxiliary tags.
pub(crate) fn minimap2_header(
    id: &[u8],
    cr: &[u8],
    cb: &[u8],
    ur: &[u8],
    ub: &[u8],
    cy: &[u8],
    uy: &[u8],
) -> Vec<u8> {
    let mut header = Vec::with_capacity(id.len() + 6 * 32 + 16);
    header.extend_from_slice(id);
    for (tag, value) in [
        (b"\tCR:Z:".as_slice(), cr),
        (b"\tCB:Z:".as_slice(), cb),
        (b"\tUR:Z:".as_slice(), ur),
        (b"\tUB:Z:".as_slice(), ub),
        (b"\tCY:Z:".as_slice(), cy),
        (b"\tUY:Z:".as_slice(), uy),
    ] {
        header.extend_from_slice(tag);
        header.extend_from_slice(value);
    }
    header
}

pub(crate) fn append_processing_tags(header: &mut Vec<u8>, tags: ProcessingTags) {
    header.extend_from_slice(b"\tXO:A:");
    header.push(tags.orientation);
    for (name, value) in [
        ("\tXB:i:", tags.barcode_route as usize),
        ("\tXI:i:", tags.index_route as usize),
        ("\tXA:i:", tags.anchor_bits as usize),
        ("\tXL:i:", tags.extracted_cdna_len),
    ] {
        header.extend_from_slice(name.as_bytes());
        header.extend_from_slice(value.to_string().as_bytes());
    }
}

/// Append the Flex probe-barcode tags (docs/archive/FLEX_SUPPORT_PLAN.md Phase 1).
pub(crate) fn append_probe_tags(header: &mut Vec<u8>, raw: &[u8], id: Option<&str>) {
    header.extend_from_slice(b"\tPR:Z:");
    header.extend_from_slice(raw);
    header.extend_from_slice(b"\tPB:Z:");
    header.extend_from_slice(id.unwrap_or("").as_bytes());
}

/// Allocate the per-LOGICAL-target raw FASTQ buffers `(R1 bytes, R2 bytes, count)`.
pub(crate) fn make_target_bufs(n_targets: usize) -> Vec<(Vec<u8>, Vec<u8>, u64)> {
    (0..n_targets)
        .map(|_| {
            (
                Vec::with_capacity(FLUSH_BYTES + 4096),
                Vec::with_capacity(FLUSH_BYTES + 4096),
                0,
            )
        })
        .collect()
}

/// Flush this worker's buffered failed reads to the side channel.
pub(crate) fn flush_failed_buf(buf: &mut Vec<u8>, sh: &Shared) -> Result<()> {
    if !buf.is_empty() {
        if let Some(writer) = sh.failed {
            writer.append(buf)?;
        }
        buf.clear();
    }
    Ok(())
}

/// Flush any records still buffered per target to this worker's `chunk` shard.
pub(crate) fn flush_target_bufs(
    bufs: &[(Vec<u8>, Vec<u8>, u64)],
    sh: &Shared,
    chunk: usize,
) -> Result<()> {
    for (i, buf) in bufs.iter().enumerate() {
        if buf.2 > 0 {
            sh.targets[i * sh.split + chunk].append(&buf.0, &buf.1, buf.2)?;
        }
    }
    Ok(())
}
