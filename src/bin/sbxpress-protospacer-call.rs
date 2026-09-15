use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use tracing_subscriber::EnvFilter;

use sbx10x::protospacer_calling::{self, Cli};

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .with_target(false)
        .init();
    let cli = if std::env::var_os("SBXPRESS_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("sbxpress protospacer-call").get_matches())?
    } else {
        Cli::parse()
    };
    protospacer_calling::run(&cli)
}
