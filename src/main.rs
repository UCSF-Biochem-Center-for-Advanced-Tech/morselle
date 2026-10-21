//! CLI wiring for the morselle-extract binary.

use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use tracing_subscriber::EnvFilter;

use morselle::config::Cli;

fn main() -> Result<()> {
    let mut cli = if std::env::var_os("MORSELLE_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("morselle extract").get_matches())?
    } else {
        Cli::parse()
    };
    cli.params.apply_presets();

    // Logging: respect RUST_LOG, default to info.
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .with_target(false)
        .init();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok();
    }

    morselle::run(&cli)
}
