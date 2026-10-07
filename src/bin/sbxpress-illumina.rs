//! Convert standard paired Illumina 10x GEX FASTQs into the tagged R2-only
//! FASTQ contract consumed by `minimap2 -y` in the sbxpress workflow.

use anyhow::Result;
use clap::Parser;
use sbx10x::illumina::{run, Cli};

fn main() -> Result<()> {
    run(Cli::parse())
}
