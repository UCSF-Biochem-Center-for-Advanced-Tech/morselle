//! Optional parallel gzip input decoding.
//!
//! Serial gzip inflate caps one input file to one core, which is the ceiling the
//! pipelined path ([`crate::run_pipelined`]) runs into on a single large FASTQ:
//! one reader thread feeds the whole worker pool, so above roughly two dozen
//! workers the reader — not the pool — sets the run's throughput.
//!
//! [`rapidgzip-core`] lifts that ceiling by decoding one gzip member across many
//! threads. It is decode-only and does not touch the output writers in
//! [`crate::io`]. Because the benefit depends entirely on how the input was
//! framed — BGZF and `pigz --independent` parallelize well, a plain `gzip` stream
//! moderately, and `pigz`'s default dependent-block framing not at all — the
//! backend is selected at *runtime* by environment variable so both can be
//! measured on the same binary and the same data:
//!
//! | Variable | Default | Meaning |
//! |---|---|---|
//! | `MORSELLE_RAPIDGZIP` | off | `1`/`true`/`on`/`yes` enables parallel decode |
//! | `MORSELLE_RAPIDGZIP_THREADS` | worker count | Per-file decoder worker budget |
//! | `MORSELLE_RAPIDGZIP_POOL` | worker count | Process-wide ceiling on decode slots |
//!
//! The pool matters on the one-file-per-thread path, where many files decode at
//! once: without an aggregate ceiling each file would claim its own worker budget
//! and oversubscribe the machine. Both paths share one pool, so total decode
//! concurrency stays bounded however the executor routes the run.
//!
//! Output framing is unaffected: each FASTQ this crate writes remains a single
//! gzip member, as 10x reference requires.

use std::path::Path;
use std::sync::OnceLock;

use anyhow::{Context, Result};
use needletail::{parse_fastx_file, parse_fastx_reader, FastxReader};
use rapidgzip_core::{Decoder, DecoderPool};

/// Resolved parallel-decode configuration, built once on first use.
struct Config {
    /// Pre-configured decoder; cloning is cheap and it is safe to `open` from
    /// many threads at once.
    decoder: Decoder,
    threads: usize,
    pool_workers: usize,
}

static CONFIG: OnceLock<Option<Config>> = OnceLock::new();

fn env_flag(key: &str) -> bool {
    match std::env::var(key) {
        Ok(v) => matches!(
            v.trim().to_ascii_lowercase().as_str(),
            "1" | "true" | "on" | "yes"
        ),
        Err(_) => false,
    }
}

fn env_usize(key: &str) -> Option<usize> {
    std::env::var(key)
        .ok()
        .and_then(|v| v.trim().parse::<usize>().ok())
        .filter(|n| *n > 0)
}

/// Build the decoder configuration, or `None` when parallel decode is disabled.
///
/// A decoder that fails to build (an invalid worker budget) disables the feature
/// rather than failing the run: the flate2 path stays correct, so a bad tuning
/// value should not abort a long job.
fn config() -> Option<&'static Config> {
    CONFIG
        .get_or_init(|| {
            if !env_flag("MORSELLE_RAPIDGZIP") {
                return None;
            }
            let workers = rayon::current_num_threads().max(1);
            let threads = env_usize("MORSELLE_RAPIDGZIP_THREADS").unwrap_or(workers);
            let pool_workers = env_usize("MORSELLE_RAPIDGZIP_POOL").unwrap_or(workers);

            let pool = DecoderPool::builder().workers(pool_workers).build().ok()?;
            let decoder = Decoder::builder()
                .decoder_threads(threads)
                .decoder_pool(pool)
                .build()
                .ok()?;
            Some(Config {
                decoder,
                threads,
                pool_workers,
            })
        })
        .as_ref()
}

/// One-line description of the active input decoder, for the run log.
pub fn backend_description() -> String {
    match config() {
        Some(c) => format!(
            "rapidgzip (parallel; {} decoder threads/file, {} pool slots)",
            c.threads, c.pool_workers
        ),
        None => "flate2/zlib (serial per file)".to_string(),
    }
}

/// Whether parallel decode is active. Callers use this only for reporting.
pub fn is_parallel() -> bool {
    config().is_some()
}

fn is_gzip(path: &Path) -> bool {
    path.file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.ends_with(".gz"))
}

/// Open a FASTQ for a full pass over every record.
///
/// Uses the parallel decoder when enabled and the input is gzip; plain FASTQ and
/// the disabled case both fall through to needletail's own file handling, which
/// keeps the flate2 behaviour byte-for-byte.
pub fn open_fastx(path: &Path) -> Result<Box<dyn FastxReader>> {
    let Some(cfg) = config().filter(|_| is_gzip(path)) else {
        return Ok(parse_fastx_file(path)?);
    };
    let reader = cfg
        .decoder
        .open(path)
        .with_context(|| format!("cannot open {} for parallel gzip decode", path.display()))?;
    Ok(parse_fastx_reader(reader)?)
}

/// Open a FASTQ for a short head scan that stops well before EOF.
///
/// Always the serial decoder. Orientation detection reads a couple of thousand
/// records and drops the reader, so spinning up a parallel decode — and the
/// speculative work it does to get ahead — would be spent and then discarded.
pub fn open_fastx_scan(path: &Path) -> Result<Box<dyn FastxReader>> {
    Ok(parse_fastx_file(path)?)
}
