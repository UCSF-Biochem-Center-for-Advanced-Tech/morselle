use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use tracing_subscriber::EnvFilter;

use morselle::annotate::{self, Cli};

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        // Annotate may emit BAM on stdout. Keep all human-readable logging on
        // stderr so the BAM stream remains byte-for-byte valid for downstream
        // consumers such as `samtools sort`.
        .with_writer(std::io::stderr)
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .with_target(false)
        .init();
    let cli = if std::env::var_os("MORSELLE_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("morselle annotate").get_matches())?
    } else {
        Cli::parse()
    };
    annotate::run(&cli)
}
