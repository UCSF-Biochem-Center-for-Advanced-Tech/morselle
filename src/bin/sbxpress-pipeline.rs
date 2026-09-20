use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use sbx10x::orchestrator::{self, Cli};

fn main() -> Result<()> {
    let cli = if std::env::var_os("SBXPRESS_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("sbxpress pipeline").get_matches())?
    } else {
        Cli::parse()
    };
    orchestrator::run(&cli)
}
