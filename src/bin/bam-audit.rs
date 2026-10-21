//! `bam-audit` — single-pass, O(1)-memory audit of a 10x reference BAM.

use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use morselle::bam_audit::{run, Cli};

fn main() -> Result<()> {
    let cli = if std::env::var_os("MORSELLE_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("morselle audit").get_matches())?
    } else {
        Cli::parse()
    };
    run(&cli)
}
