//! sbx10x — convert raw Roche SBX reads into per-sample 10x R1/R2 FASTQs.
//!
//! See `sbx10x_implementation_plan.md` for the milestone-phased design. M0 wires
//! up the crate, CLI, logging, and a FASTQ read->write passthrough; subsequent
//! milestones fill in orientation, demux, extraction, gating, and realignment.

pub mod accounting;
pub mod anchors;
pub mod barcode;
pub mod config;
pub mod demux;
pub mod extract;
pub mod homopolymer;
pub mod io;
pub mod orient;
pub mod pipeline;

use std::path::Path;

use anyhow::Result;
use needletail::parse_fastx_file;
use tracing::info;

use crate::config::Cli;
use crate::io::{discover_fastqs, FastqGzWriter};

/// Entry point invoked by the binary. M0: passthrough — read every input record
/// and re-emit it to a single FASTQ in the output dir, counting records.
pub fn run(cli: &Cli) -> Result<()> {
    let inputs = discover_fastqs(&cli.fastqs)?;
    info!(n_files = inputs.len(), "discovered FASTQ inputs");
    if inputs.is_empty() {
        anyhow::bail!("no FASTQ inputs found");
    }

    let out_path = cli.out_dir.join("passthrough.fastq.gz");
    let mut writer = FastqGzWriter::create(&out_path)?;

    let mut total: u64 = 0;
    for input in &inputs {
        total += passthrough_file(input, &mut writer)?;
    }
    writer.finish()?;

    info!(
        reads = total,
        out = %out_path.display(),
        "M0 passthrough complete"
    );
    Ok(())
}

/// Stream one FASTQ file's records to the writer, returning the count.
fn passthrough_file(input: &Path, writer: &mut FastqGzWriter) -> Result<u64> {
    let mut reader = parse_fastx_file(input)?;
    let mut n = 0u64;
    while let Some(rec) = reader.next() {
        let rec = rec?;
        let qual = rec.qual().unwrap_or(&[]);
        writer.write_record(rec.id(), &rec.seq(), qual)?;
        n += 1;
    }
    info!(file = %input.display(), reads = n, "read FASTQ");
    Ok(n)
}
