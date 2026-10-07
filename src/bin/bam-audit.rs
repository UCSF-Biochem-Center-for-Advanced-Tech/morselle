//! `bam-audit` — single-pass, O(1)-memory audit of a 10x reference BAM.

use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use sbx10x::bam_audit::{run, Cli};

fn main() -> Result<()> {
    let cli = if std::env::var_os("SBXPRESS_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("sbxpress audit").get_matches())?
    } else {
        Cli::parse()
    };
    run(&cli)
}
