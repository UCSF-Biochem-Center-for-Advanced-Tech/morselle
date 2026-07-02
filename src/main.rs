//! CLI wiring for the sbx10x binary.

use anyhow::Result;
use clap::Parser;
use tracing_subscriber::EnvFilter;

use sbx10x::config::Cli;

fn main() -> Result<()> {
    let mut cli = Cli::parse();
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

    sbx10x::run(&cli)
}
