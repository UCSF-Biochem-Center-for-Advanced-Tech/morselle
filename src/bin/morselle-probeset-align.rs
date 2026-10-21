use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use tracing_subscriber::EnvFilter;

use morselle::probe_set::{self, Cli};

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .with_target(false)
        .init();
    let cli = if std::env::var_os("MORSELLE_DISPATCH").is_some() {
        Cli::from_arg_matches(
            &Cli::command()
                .bin_name("morselle probeset-align")
                .get_matches(),
        )?
    } else {
        Cli::parse()
    };
    probe_set::run(&cli)
}
