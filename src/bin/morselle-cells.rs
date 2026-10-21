use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use morselle::cell_calling::{self, Cli};

fn main() -> Result<()> {
    let cli = if std::env::var_os("MORSELLE_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("morselle cells").get_matches())?
    } else {
        Cli::parse()
    };
    cell_calling::run(&cli)
}
