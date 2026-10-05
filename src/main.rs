//! CLI wiring for the sbx10x binary.

use anyhow::Result;
use clap::{CommandFactory, FromArgMatches, Parser};
use tracing_subscriber::EnvFilter;

use sbx10x::config::Cli;

fn main() -> Result<()> {
    let mut cli = if std::env::var_os("SBXPRESS_DISPATCH").is_some() {
        Cli::from_arg_matches(&Cli::command().bin_name("sbxpress extract").get_matches())?
    } else {
        Cli::parse()
    };
    cli.params.apply_presets();

    // Logging: respect RUST_LOG, default to info. Keep all human-readable
    // logging on stderr -- `--print-manifest` writes machine-readable JSON
    // to stdout, and any log line landing on the same stream silently
    // corrupts it for a consumer expecting clean JSON (real bug found this
    // session: `sbxpress pipeline`'s manifest fetch failed to parse because
    // a `demux mode resolved` INFO line was mixed into stdout ahead of the
    // JSON -- this file was the one binary in the tree missing this).
    tracing_subscriber::fmt()
        .with_writer(std::io::stderr)
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
