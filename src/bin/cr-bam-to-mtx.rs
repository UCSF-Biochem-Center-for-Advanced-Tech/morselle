use anyhow::Result;
use clap::Parser;
use tracing_subscriber::EnvFilter;

use sbx10x::cr_bam_to_mtx::{self, Cli};

fn main() -> Result<()> {
    tracing_subscriber::fmt()
        .with_env_filter(
            EnvFilter::try_from_default_env().unwrap_or_else(|_| EnvFilter::new("info")),
        )
        .with_target(false)
        .init();
    cr_bam_to_mtx::run(&Cli::parse())
}
