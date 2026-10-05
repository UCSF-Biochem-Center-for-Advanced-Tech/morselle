use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use tracing_subscriber::EnvFilter;

use sbx10x::cr_bam_to_mtx::{self, Cli};

fn main() -> Result<()> {
    // Keep human-readable logging off stdout -- this binary's stdout is
    // reserved for its output contract; a log line mixed in silently
    // corrupts any consumer expecting clean output.
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .with_target(false)
        .init();
    let cli = if std::env::var_os("SBXPRESS_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("sbxpress matrix").get_matches())?
    } else {
        Cli::parse()
    };
    cr_bam_to_mtx::run(&cli)
}
