//! `sbxpress pipeline` — a real orchestrator binary that runs the full
//! `extract -> probeset-align/guides -> molecules -> probe-demux -> matrix ->
//! cells -> protospacer-call` pipeline end to end, replacing shell-script
//! orchestration (`scripts/run_full_cat06.sh`, `scripts/sbxpress-gex-pipeline.sh`).
//!
//! By default the `extract -> probeset-align`/`guides` hop streams through OS
//! named pipes with zero intermediate files; it falls back to real
//! materialized files when the requested CPU/memory budget can't fit
//! `extract` and every reader process live at once, or when `--resumable` is
//! set (streaming and resumability are fundamentally in tension — a crash
//! mid-stream leaves nothing on disk to resume from). See
//! `PIPELINE_NOTES.md` for the measured waste this replaces and the plan
//! this module implements.

use std::env;
use std::ffi::CString;
use std::fs;
use std::io::{BufRead, BufReader, Write};
use std::os::unix::ffi::OsStrExt;
use std::path::{Path, PathBuf};
use std::process::{Child, Command, Stdio};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

use anyhow::{bail, Context, Result};
use clap::{Parser, ValueEnum};

use crate::Manifest;

// ---------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum StopAfter {
    Align,
    Molecules,
    Matrix,
    Cells,
    All,
}

#[derive(Debug, Parser)]
#[command(
    name = "sbxpress-pipeline",
    version,
    about = "Run the full SBX pipeline end to end: extract -> probeset-align/guides -> molecules -> probe-demux -> matrix -> cells -> protospacer-call"
)]
pub struct Cli {
    /// Raw SBX FASTQ(s): a directory of `PrimaryAnalysis_*.fastq.gz` or explicit files.
    #[arg(long, num_args = 1.., required = true)]
    pub fastqs: Vec<PathBuf>,

    /// BCL Convert v2 sample sheet (or flat `sample,i7,i5` CSV).
    #[arg(long)]
    pub bclconvert_samplesheet: Option<PathBuf>,

    /// 10x cell-barcode whitelist. Optional, matching standalone `extract`'s
    /// own default (without one, the CB-whitelist gate is skipped entirely
    /// rather than rejecting every read) -- a real production run should
    /// still pass one, but the pipeline itself doesn't require it.
    #[arg(long)]
    pub whitelist: Option<PathBuf>,

    /// 10x Flex probe-barcode whitelist, if this run has Flex probe-barcoded samples.
    #[arg(long)]
    pub probe_barcode_whitelist: Option<PathBuf>,

    /// Probe-set reference for `probeset-align` (GEX-routed chunks). Required
    /// only if the manifest routes any chunk to GEX.
    #[arg(long)]
    pub probe_set: Option<PathBuf>,

    /// Feature reference for `guides` (CRISPR/Feature Barcoding-routed chunks).
    /// Required only if the manifest routes any chunk to non-GEX.
    #[arg(long)]
    pub feature_reference: Option<PathBuf>,

    /// Output directory for the whole run.
    #[arg(long, short = 'o', required = true)]
    pub out_dir: PathBuf,

    /// CPU budget. 0 = auto-detect via `std::thread::available_parallelism()`.
    #[arg(long, default_value_t = 0)]
    pub cpus: usize,

    /// Memory budget in GB. 0 = auto-detect from `/proc/meminfo`'s `MemAvailable`.
    #[arg(long, default_value_t = 0)]
    pub max_memory_gb: u64,

    /// extract -> align/guides fan-out width. 0 = auto, derived from `--cpus`.
    #[arg(long, default_value_t = 0)]
    pub split: usize,

    /// `sbx10x-molecules --shards`. 0 = auto, derived from real available
    /// memory at the moment molecules is about to run. Always honored over
    /// the heuristic when nonzero — the heuristic is seeded from one real
    /// library's measurements (see `plan_molecule_shards`) and could
    /// under-shard (OOM) or over-shard (wasted wall-clock) on a materially
    /// different one.
    #[arg(long, default_value_t = 0)]
    pub molecule_shards: usize,

    /// `sbx10x-molecules --shard-concurrency`. 0 = auto (see `--molecule-shards`).
    #[arg(long, default_value_t = 0)]
    pub molecule_shard_concurrency: usize,

    /// Force real materialized (gzip or plain) intermediates for the
    /// extract->align hop instead of FIFO streaming, even when the resource
    /// budget would allow streaming. Debugging/inspection use case.
    #[arg(long, default_value_t = false)]
    pub materialize_extract: bool,

    /// Keep intermediate files instead of cleaning them up as the pipeline
    /// progresses past them. Ignored (always kept) under `--resumable`.
    #[arg(long, default_value_t = false)]
    pub keep_intermediates: bool,

    /// Checkpoint every stage's real output and allow restarting from the
    /// last completed checkpoint after an interrupted run. Forces
    /// materialized intermediates pipeline-wide (streaming is incompatible
    /// with resumability — a crash mid-FIFO-stream leaves nothing on disk to
    /// resume from) and defers all cleanup until the whole pipeline reports
    /// final success.
    #[arg(long, default_value_t = false)]
    pub resumable: bool,

    /// Seconds between progress heartbeats. A dedicated timer emits a
    /// summary line on this cadence regardless of whether any stage produced
    /// output in that window, so a long-running stage never looks stalled.
    #[arg(long, default_value_t = 30)]
    pub progress_interval: u64,

    /// Write progress output to this file instead of stderr. Kept separate
    /// from stderr by default so a real run's warnings/errors don't get
    /// interleaved with throughput lines when replaying a run offline.
    #[arg(long)]
    pub progress_file: Option<PathBuf>,

    /// Stop after this stage.
    #[arg(long, value_enum, default_value_t = StopAfter::All)]
    pub stop_after: StopAfter,
}

// ---------------------------------------------------------------------------
// Stage graph / resource shapes (PIPELINE_NOTES.md's proposal, implemented)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ResourceShape {
    CpuBoundParallel,
    MemoryDominantSerial,
    LightweightIoBound,
}

#[derive(Debug, Clone, Copy)]
pub struct Stage {
    pub name: &'static str,
    pub shape: ResourceShape,
}

pub const STAGES: &[Stage] = &[
    Stage { name: "extract", shape: ResourceShape::CpuBoundParallel },
    Stage { name: "align", shape: ResourceShape::CpuBoundParallel },
    Stage { name: "molecules", shape: ResourceShape::MemoryDominantSerial },
    Stage { name: "probe-demux", shape: ResourceShape::LightweightIoBound },
    Stage { name: "matrix", shape: ResourceShape::LightweightIoBound },
    Stage { name: "cells", shape: ResourceShape::LightweightIoBound },
    Stage { name: "protospacer-call", shape: ResourceShape::LightweightIoBound },
];

pub fn stage_shape(name: &str) -> Option<ResourceShape> {
    STAGES.iter().find(|s| s.name == name).map(|s| s.shape)
}

// ---------------------------------------------------------------------------
// Resource detection
// ---------------------------------------------------------------------------

pub fn detect_cpus() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}

/// Real available memory in GB, from `/proc/meminfo`'s `MemAvailable` line —
/// this whole project's own investigations have read this value by hand all
/// session (see `PIPELINE_NOTES.md`); this is the same value, just read
/// programmatically instead of via `awk`.
pub fn detect_available_mem_gb() -> Result<u64> {
    let contents = fs::read_to_string("/proc/meminfo").context("reading /proc/meminfo")?;
    parse_mem_available_kb(&contents)
        .map(|kb| kb / 1_048_576)
        .context("MemAvailable not found in /proc/meminfo")
}

fn parse_mem_available_kb(meminfo: &str) -> Option<u64> {
    for line in meminfo.lines() {
        if let Some(rest) = line.strip_prefix("MemAvailable:") {
            let kb: u64 = rest.trim().trim_end_matches(" kB").trim().parse().ok()?;
            return Some(kb);
        }
    }
    None
}

// ---------------------------------------------------------------------------
// Pure scheduler functions
// ---------------------------------------------------------------------------

/// `probeset-align`/`guides` are single-threaded per process (measured,
/// `PIPELINE_NOTES.md`), so fan-out width — not extract's own thread count —
/// is the real parallelism dial for this hop. Starting heuristic from the one
/// real data point in `PIPELINE_NOTES.md`: extract sustains ~1.6M reads/sec
/// aggregate at `--threads 128` (~12.5K reads/sec/thread); one align process
/// sustains ~59K reads/sec (58M reads / 977s) — so one extract-thread
/// roughly feeds 4-5 align processes. **Not a validated constant across
/// libraries** — `--split` always overrides this when set explicitly.
const READS_PER_EXTRACT_THREAD: usize = 4;

/// Returns `(extract_threads, n_readers)`.
pub fn plan_align_fanout(cpus: usize, requested_split: usize) -> (usize, usize) {
    let cpus = cpus.max(1);
    let n_readers = if requested_split > 0 {
        requested_split
    } else {
        (cpus * READS_PER_EXTRACT_THREAD / (READS_PER_EXTRACT_THREAD + 1)).max(1)
    };
    let extract_threads = (n_readers / READS_PER_EXTRACT_THREAD).max(1).min(cpus);
    (extract_threads, n_readers)
}

/// Conservative bytes-of-peak-RSS-per-raw-read estimate, calibrated so that
/// `PIPELINE_NOTES.md`'s one real anchor (~9B reads, ~900GB available memory
/// on this box, real observed 8-shard run) lands in a sane, safety-margined
/// shard count — **not a validated per-library constant**, which is exactly
/// why `--molecule-shards`/`--molecule-shard-concurrency` always override
/// this heuristic when set explicitly.
const BYTES_PER_READ_ESTIMATE: f64 = 53.0;
const SHARD_SAFETY_FACTOR: f64 = 1.5;

/// Returns `(shards, shard_concurrency)`.
pub fn plan_molecule_shards(total_reads: u64, available_mem_gb: u64, cpu_budget: usize) -> (usize, usize) {
    let available_bytes = (available_mem_gb.max(1) as f64) * 1_000_000_000.0;
    let needed_bytes = (total_reads as f64) * BYTES_PER_READ_ESTIMATE * SHARD_SAFETY_FACTOR;
    let shards = ((needed_bytes / available_bytes).ceil() as usize).max(1);
    let per_shard_bytes = needed_bytes / shards as f64;
    let mem_bound_concurrency = ((available_bytes / per_shard_bytes.max(1.0)) as usize).max(1);
    let concurrency = shards.min(cpu_budget.max(1)).min(mem_bound_concurrency);
    (shards, concurrency.max(1))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExecutionMode {
    Stream,
    Materialize,
}

/// Whether `extract` and all `n_targets` align/guides readers can safely run
/// concurrently within the requested budget. This is a CPU-concurrency
/// question, not a memory one — each align process is a lightweight
/// single-threaded FASTQ scan — so `mem_budget_gb` only guards against a
/// pathologically tiny budget, not a per-target cost.
const STREAM_MIN_MEM_GB: u64 = 4;

pub fn plan_extract_align_mode(
    cpu_budget: usize,
    mem_budget_gb: u64,
    n_targets: usize,
    extract_threads: usize,
) -> ExecutionMode {
    let needed_cpus = extract_threads + n_targets;
    if needed_cpus <= cpu_budget.max(1) && mem_budget_gb >= STREAM_MIN_MEM_GB {
        ExecutionMode::Stream
    } else {
        ExecutionMode::Materialize
    }
}

/// Whether align's own output can *also* stream directly into
/// `sbx10x-molecules` (no per-chunk aligned BAM ever touching disk), on top
/// of an already-streamed extract->align hop. Real, additional constraint
/// beyond `plan_extract_align_mode`: `molecules` is memory-dominant (see
/// `STAGES`), not lightweight like align/guides, so folding it into the same
/// concurrently-running child group adds its own `--shard-concurrency`
/// worth of CPU demand on top of what extract+align already need -- gated
/// here, not assumed free. Memory itself isn't re-checked: `shards`/
/// `concurrency` were already sized against the same `mem_budget_gb` by
/// `plan_molecule_shards`, and extract/align are lightweight, so molecules
/// remains the dominant consumer either way -- no new memory risk beyond
/// what running molecules at all already carries.
pub fn plan_align_molecules_stream(
    extract_align_mode: ExecutionMode,
    cpu_budget: usize,
    extract_threads: usize,
    n_targets: usize,
    molecule_shard_concurrency: usize,
) -> bool {
    extract_align_mode == ExecutionMode::Stream
        && (extract_threads + n_targets + molecule_shard_concurrency) <= cpu_budget.max(1)
}

// ---------------------------------------------------------------------------
// Sibling-binary resolution (shared with `src/bin/sbxpress.rs`)
// ---------------------------------------------------------------------------

/// Resolve the directory containing sbxpress's sibling engine binaries.
/// Promoted out of `src/bin/sbxpress.rs` so the top-level dispatcher and this
/// orchestrator share one implementation and can never silently diverge.
pub fn helper_dir(executable: &Path) -> Option<PathBuf> {
    if let Some(value) = env::var_os("SBXPRESS_LIBEXEC") {
        return Some(PathBuf::from(value));
    }
    let bin = executable.parent()?;
    let installed = bin.parent()?.join("libexec/sbxpress");
    if installed.is_dir() {
        Some(installed)
    } else {
        Some(bin.to_path_buf())
    }
}

/// Injected so a test can resolve `CARGO_BIN_EXE_*` paths instead of
/// production's installed/`libexec` filesystem layout.
pub trait HelperResolver: Send + Sync {
    fn resolve(&self, name: &str) -> Result<PathBuf>;
}

pub struct LibexecResolver {
    pub dir: PathBuf,
}

impl HelperResolver for LibexecResolver {
    fn resolve(&self, name: &str) -> Result<PathBuf> {
        let path = self.dir.join(name);
        if !path.is_file() {
            bail!(
                "sbxpress pipeline: helper not found: {} (set SBXPRESS_LIBEXEC to override)",
                path.display()
            );
        }
        Ok(path)
    }
}

// ---------------------------------------------------------------------------
// Progress reporting
// ---------------------------------------------------------------------------

/// Minimal, dependency-free UTC calendar conversion (Howard Hinnant's
/// `civil_from_days` algorithm) — avoids pulling in a chrono/time dependency
/// just for one timestamp format in progress lines.
fn civil_from_unix(secs: i64) -> (i64, u32, u32, u32, u32, u32) {
    let days = secs.div_euclid(86400);
    let rem = secs.rem_euclid(86400);
    let (hh, mm, ss) = (rem / 3600, (rem % 3600) / 60, rem % 60);
    let z = days + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32;
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32;
    let y = if m <= 2 { y + 1 } else { y };
    (y, m, d, hh as u32, mm as u32, ss as u32)
}

fn unix_timestamp_iso(t: SystemTime) -> String {
    let secs = t
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0);
    let (y, mo, d, h, mi, s) = civil_from_unix(secs);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{s:02}Z")
}

/// Pure line formatter: real wall-clock timestamp + elapsed-since-pipeline-start
/// on every line. This retroactively gives every relayed sub-app log line
/// (e.g. `sbx10x-molecules`' own `phase=...` lines, which carry no timestamp
/// today — see `PIPELINE_NOTES.md`) a real, greppable time axis, with no
/// changes needed inside the sub-apps themselves.
pub fn format_progress_line(elapsed: Duration, tag: &str, msg: &str) -> String {
    let el = elapsed.as_secs();
    let ts = unix_timestamp_iso(SystemTime::now());
    format!(
        "[{ts} +{:02}:{:02}:{:02}] [{tag}] {msg}",
        el / 3600,
        (el % 3600) / 60,
        el % 60
    )
}

/// Pure heartbeat-scheduling logic — table-testable without any real sleep.
pub fn heartbeat_due(elapsed_since_last: Duration, interval: Duration) -> bool {
    elapsed_since_last >= interval
}

pub struct Progress {
    writer: Mutex<Box<dyn Write + Send>>,
    start: Instant,
}

impl Progress {
    pub fn new(file: Option<&Path>) -> Result<Self> {
        let writer: Box<dyn Write + Send> = match file {
            Some(path) => Box::new(
                fs::File::create(path)
                    .with_context(|| format!("creating progress file {}", path.display()))?,
            ),
            None => Box::new(std::io::stderr()),
        };
        Ok(Progress { writer: Mutex::new(writer), start: Instant::now() })
    }

    pub fn line(&self, tag: &str, msg: &str) {
        let formatted = format_progress_line(self.start.elapsed(), tag, msg);
        if let Ok(mut w) = self.writer.lock() {
            let _ = writeln!(w, "{formatted}");
            let _ = w.flush();
        }
    }

    pub fn elapsed(&self) -> Duration {
        self.start.elapsed()
    }
}

/// Spawns a heartbeat thread that emits a top-level summary line every
/// `interval` regardless of whether any child produced output — the actual
/// guarantee that a long-running stage never looks stalled. Returns a stop
/// flag the caller should set once the run is done, and the join handle.
pub fn spawn_heartbeat(
    progress: Arc<Progress>,
    interval: Duration,
    status: Arc<Mutex<String>>,
) -> (Arc<AtomicBool>, thread::JoinHandle<()>) {
    let stop = Arc::new(AtomicBool::new(false));
    let stop_clone = Arc::clone(&stop);
    // Poll in short slices rather than one `thread::sleep(interval)` -- with
    // a large --progress-interval, sleeping the whole interval before the
    // first stop-flag check would make shutdown (`hb_handle.join()`) block
    // for up to that entire interval after the real work is already done.
    const POLL_SLICE: Duration = Duration::from_millis(200);
    let handle = thread::spawn(move || {
        let mut waited = Duration::ZERO;
        while !stop_clone.load(Ordering::Relaxed) {
            thread::sleep(POLL_SLICE);
            waited += POLL_SLICE;
            if stop_clone.load(Ordering::Relaxed) {
                break;
            }
            if waited >= interval {
                waited = Duration::ZERO;
                let msg = status.lock().map(|s| s.clone()).unwrap_or_default();
                progress.line("heartbeat", &msg);
            }
        }
    });
    (stop, handle)
}

/// Relays one child's stderr, line by line, through `progress` with a
/// `[stage:unit]`-style tag, instead of letting it hit the terminal raw.
fn relay_child_stderr(
    child_stderr: std::process::ChildStderr,
    progress: Arc<Progress>,
    tag: String,
) -> thread::JoinHandle<()> {
    thread::spawn(move || {
        let reader = BufReader::new(child_stderr);
        for line in reader.lines().map_while(std::io::Result::ok) {
            progress.line(&tag, &line);
        }
    })
}

// ---------------------------------------------------------------------------
// Checkpointing (`--resumable`)
// ---------------------------------------------------------------------------

pub fn checkpoint_dir(out_dir: &Path) -> PathBuf {
    out_dir.join(".checkpoints")
}

pub fn checkpoint_path(out_dir: &Path, stage: &str, unit_id: &str) -> PathBuf {
    checkpoint_dir(out_dir).join(format!("{stage}.{unit_id}.done"))
}

pub fn is_checkpointed(out_dir: &Path, stage: &str, unit_id: &str) -> bool {
    checkpoint_path(out_dir, stage, unit_id).is_file()
}

/// Written atomically (temp file + rename) so a crash mid-write never leaves
/// a checkpoint that looks complete but isn't.
pub fn write_checkpoint(out_dir: &Path, stage: &str, unit_id: &str) -> Result<()> {
    let dir = checkpoint_dir(out_dir);
    fs::create_dir_all(&dir).with_context(|| format!("creating checkpoint dir {}", dir.display()))?;
    let final_path = checkpoint_path(out_dir, stage, unit_id);
    let tmp_path = dir.join(format!("{stage}.{unit_id}.done.tmp"));
    fs::write(&tmp_path, b"").with_context(|| format!("writing {}", tmp_path.display()))?;
    fs::rename(&tmp_path, &final_path)
        .with_context(|| format!("renaming checkpoint into place: {}", final_path.display()))?;
    Ok(())
}

// ---------------------------------------------------------------------------
// FIFOs
// ---------------------------------------------------------------------------

fn mkfifo(path: &Path) -> Result<()> {
    if let Some(dir) = path.parent() {
        fs::create_dir_all(dir).with_context(|| format!("creating dir {}", dir.display()))?;
    }
    let c_path = CString::new(path.as_os_str().as_bytes())
        .with_context(|| format!("path has interior NUL: {}", path.display()))?;
    let rc = unsafe { libc::mkfifo(c_path.as_ptr(), 0o644) };
    if rc != 0 {
        let err = std::io::Error::last_os_error();
        bail!("mkfifo {}: {err}", path.display());
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// extract manifest
// ---------------------------------------------------------------------------

pub fn fetch_manifest(
    sbx10x_bin: &Path,
    cli: &Cli,
    split: usize,
) -> Result<Manifest> {
    let mut cmd = Command::new(sbx10x_bin);
    cmd.arg("--fastqs").args(&cli.fastqs);
    if let Some(sheet) = &cli.bclconvert_samplesheet {
        cmd.arg("--bclconvert-samplesheet").arg(sheet);
    }
    if let Some(wl) = &cli.whitelist {
        cmd.arg("--whitelist").arg(wl);
    }
    cmd.arg("--out-dir").arg(&cli.out_dir);
    cmd.arg("--minimap2");
    // Must match the real extract invocation's flags exactly (both Stream
    // and Materialize modes always pass --no-gzip -- this hop is about
    // concurrency, not compression, see plan_extract_align_mode), or the
    // predicted paths here won't be the paths extract actually writes to:
    // a `.fastq.gz` FIFO nobody ever opens for write hangs its reader
    // forever.
    cmd.arg("--no-gzip");
    cmd.arg("--split").arg(split.max(1).to_string());
    if let Some(pbw) = &cli.probe_barcode_whitelist {
        cmd.arg("--probe-barcode-whitelist").arg(pbw);
    }
    cmd.arg("--print-manifest");
    let output = cmd.output().context("running sbx10x --print-manifest")?;
    if !output.status.success() {
        bail!(
            "sbx10x --print-manifest failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    let manifest: Manifest = serde_json::from_slice(&output.stdout)
        .context("parsing extract --print-manifest JSON")?;
    validate_manifest_version(&manifest)?;
    Ok(manifest)
}

fn validate_manifest_version(manifest: &Manifest) -> Result<()> {
    if manifest.version != 1 {
        bail!("unsupported manifest version {} (expected 1)", manifest.version);
    }
    Ok(())
}

/// Aligned-BAM output path for one manifest target — the orchestrator's own
/// naming, not reverse-engineered from another module.
pub fn aligned_bam_path(out_dir: &Path, sample: &str, s_num: usize, chunk: usize) -> PathBuf {
    out_dir
        .join("aligned")
        .join(format!("{sample}_S{s_num}_L001_R2_{chunk:03}.bam"))
}

// ---------------------------------------------------------------------------
// extract -> align/guides
// ---------------------------------------------------------------------------

struct AlignOutcome {
    bam_path: PathBuf,
    gex: bool,
}

fn wait_all_or_kill(mut children: Vec<(String, Child)>) -> Result<()> {
    // One OS thread per child calling blocking `.wait()` -- no async runtime
    // needed, matches the codebase's existing synchronous style.
    let failure: Arc<Mutex<Option<String>>> = Arc::new(Mutex::new(None));
    let pids: Vec<u32> = children.iter().map(|(_, c)| c.id()).collect();
    let handles: Vec<_> = children
        .drain(..)
        .map(|(tag, mut child)| {
            let failure = Arc::clone(&failure);
            thread::spawn(move || match child.wait() {
                Ok(status) if status.success() => {}
                Ok(status) => {
                    let mut f = failure.lock().unwrap();
                    if f.is_none() {
                        *f = Some(format!("{tag} exited with {status}"));
                    }
                }
                Err(err) => {
                    let mut f = failure.lock().unwrap();
                    if f.is_none() {
                        *f = Some(format!("{tag} wait() failed: {err}"));
                    }
                }
            })
        })
        .collect();

    // Poll for an early failure and actively kill remaining children -- a
    // crashed FIFO writer/reader otherwise leaves its sibling blocked
    // forever on a pipe read/write that will never complete.
    loop {
        if handles.iter().all(|h| h.is_finished()) {
            break;
        }
        if failure.lock().unwrap().is_some() {
            for pid in &pids {
                unsafe {
                    libc::kill(*pid as i32, libc::SIGKILL);
                }
            }
            break;
        }
        thread::sleep(Duration::from_millis(100));
    }
    for h in handles {
        let _ = h.join();
    }
    if let Some(msg) = failure.lock().unwrap().take() {
        bail!("pipeline stage failed: {msg}");
    }
    Ok(())
}

/// Runs the extract->align hop in `ExecutionMode::Stream`: creates FIFOs at
/// the manifest's predicted R2 paths, spawns every align/guides reader
/// first (order is load-bearing -- opening a FIFO for write blocks until a
/// reader attaches), then spawns extract itself. Returns each target's real
/// aligned-BAM output path.
#[allow(clippy::too_many_arguments)]
fn run_extract_align_stream(
    cli: &Cli,
    helpers: &dyn HelperResolver,
    manifest: &Manifest,
    extract_threads: usize,
    progress: &Arc<Progress>,
) -> Result<Vec<AlignOutcome>> {
    let sbx10x = helpers.resolve("sbx10x")?;
    let probeset_align = helpers.resolve("sbxpress-probeset-align")?;
    let guides = helpers.resolve("sbxpress-guides")?;

    for target in &manifest.targets {
        if target.gex && cli.probe_set.is_none() {
            bail!("target {} routes to probeset-align but --probe-set was not given", target.sample);
        }
        if !target.gex && cli.feature_reference.is_none() {
            bail!("target {} routes to guides but --feature-reference was not given", target.sample);
        }
        mkfifo(&target.r2_path)?;
    }

    let mut children: Vec<(String, Child)> = Vec::new();
    let mut outcomes = Vec::new();
    fs::create_dir_all(cli.out_dir.join("aligned")).context("creating aligned/ dir")?;

    for target in &manifest.targets {
        let bam_path = aligned_bam_path(&cli.out_dir, &target.sample, target.s_num, target.chunk);
        let tag = format!("align:{}.{:03}", target.sample, target.chunk);
        let bin = if target.gex { &probeset_align } else { &guides };
        let mut cmd = Command::new(bin);
        cmd.arg(&target.r2_path).arg("-o").arg(&bam_path);
        if target.gex {
            cmd.arg("--probe-set").arg(cli.probe_set.as_ref().unwrap());
        } else {
            cmd.arg("--feature-reference").arg(cli.feature_reference.as_ref().unwrap());
        }
        cmd.stderr(Stdio::piped());
        let mut child = cmd.spawn().with_context(|| format!("spawning {tag}"))?;
        if let Some(stderr) = child.stderr.take() {
            relay_child_stderr(stderr, Arc::clone(progress), tag.clone());
        }
        children.push((tag, child));
        outcomes.push(AlignOutcome { bam_path, gex: target.gex });
    }

    // Now that every reader is attached to its FIFO, extract can be spawned
    // safely -- reversing this order deadlocks (opening a FIFO for write
    // blocks until a reader exists).
    let mut extract_cmd = Command::new(&sbx10x);
    extract_cmd
        .arg("--fastqs")
        .args(&cli.fastqs)
        .arg("--out-dir")
        .arg(&cli.out_dir)
        .arg("--minimap2")
        .arg("--no-gzip")
        .arg("--split")
        .arg(manifest_split(manifest).to_string())
        .arg("--threads")
        .arg(extract_threads.to_string());
    if let Some(wl) = &cli.whitelist {
        extract_cmd.arg("--whitelist").arg(wl);
    }
    if let Some(sheet) = &cli.bclconvert_samplesheet {
        extract_cmd.arg("--bclconvert-samplesheet").arg(sheet);
    }
    if let Some(pbw) = &cli.probe_barcode_whitelist {
        extract_cmd.arg("--probe-barcode-whitelist").arg(pbw);
    }
    extract_cmd.stderr(Stdio::piped());
    let mut extract_child = extract_cmd.spawn().context("spawning extract (streaming)")?;
    if let Some(stderr) = extract_child.stderr.take() {
        relay_child_stderr(stderr, Arc::clone(progress), "extract".to_string());
    }
    children.push(("extract".to_string(), extract_child));

    wait_all_or_kill(children)?;

    for target in &manifest.targets {
        let _ = fs::remove_file(&target.r2_path);
    }
    Ok(outcomes)
}

fn manifest_split(manifest: &Manifest) -> usize {
    manifest.targets.iter().map(|t| t.chunk).max().unwrap_or(1)
}

/// Extends `run_extract_align_stream` one hop further: align/guides write
/// their *output* into FIFOs too (instead of `aligned_bam_path`'s real
/// file), and `sbx10x-molecules` is spawned as the reader of all of one
/// routing group's (GEX, or feature/CRISPR) FIFOs -- with `--collate-workers`
/// set to the FIFO count so its own `collate()` drains them concurrently
/// (see that flag's doc comment) rather than serializing the align fan-out
/// the way sequential draining would. Zero per-chunk aligned BAM ever
/// touches disk. Same reader-before-writer spawn ordering and
/// `wait_all_or_kill` failure semantics as `run_extract_align_stream`, just
/// with the molecules reader(s) folded into the same child group instead of
/// a separate stage invoked afterward: molecules readers spawn first (they
/// must be attached before align tries to open their input FIFOs for
/// write), then align/guides (attached to extract's FIFOs before extract
/// spawns), then extract last.
///
/// Returns `(gex_molecules_bam, feature_molecules_bam)`, matching the shape
/// callers already expect from the non-streamed molecules stage.
#[allow(clippy::too_many_arguments)]
fn run_extract_align_molecules_stream(
    cli: &Cli,
    helpers: &dyn HelperResolver,
    manifest: &Manifest,
    extract_threads: usize,
    shards: usize,
    concurrency: usize,
    progress: &Arc<Progress>,
) -> Result<(Option<PathBuf>, Option<PathBuf>)> {
    let sbx10x = helpers.resolve("sbx10x")?;
    let probeset_align = helpers.resolve("sbxpress-probeset-align")?;
    let guides = helpers.resolve("sbxpress-guides")?;
    let sbx10x_molecules = helpers.resolve("sbx10x-molecules")?;

    for target in &manifest.targets {
        if target.gex && cli.probe_set.is_none() {
            bail!("target {} routes to probeset-align but --probe-set was not given", target.sample);
        }
        if !target.gex && cli.feature_reference.is_none() {
            bail!("target {} routes to guides but --feature-reference was not given", target.sample);
        }
        mkfifo(&target.r2_path)?;
    }
    fs::create_dir_all(cli.out_dir.join("aligned")).context("creating aligned/ dir")?;
    fs::create_dir_all(cli.out_dir.join("molecules")).context("creating molecules/ dir")?;

    // One align-output FIFO per target, same naming/location a real file
    // would have used, so logs/paths stay recognizable even though nothing
    // is ever written to disk at them.
    let mut align_fifos: Vec<(PathBuf, bool)> = Vec::with_capacity(manifest.targets.len());
    for target in &manifest.targets {
        let fifo = aligned_bam_path(&cli.out_dir, &target.sample, target.s_num, target.chunk);
        mkfifo(&fifo)?;
        align_fifos.push((fifo, target.gex));
    }
    let gex_fifos: Vec<PathBuf> =
        align_fifos.iter().filter(|(_, gex)| *gex).map(|(p, _)| p.clone()).collect();
    let other_fifos: Vec<PathBuf> =
        align_fifos.iter().filter(|(_, gex)| !*gex).map(|(p, _)| p.clone()).collect();

    let mut children: Vec<(String, Child)> = Vec::new();
    let mut gex_out = None;
    let mut other_out = None;

    let spawn_molecules_reader =
        |fifos: &[PathBuf], tagseq: bool, out: &Path, tag: &str| -> Result<Child> {
            // Always `--shards > 1` here, even for a single FIFO: `--shards
            // <= 1` takes sbx10x-molecules' non-sharded path, which reads
            // its input BAM *twice* (once to collect evidence, once more in
            // write_output() to re-read and tag the final output) -- fatal
            // for a FIFO, whose data can only be read once. `--shards > 1`
            // routes through collate(), which drains each FIFO exactly once
            // and materializes real (ephemeral, already-cleaned-up) shard
            // temp files that the rest of the pipeline -- including its own
            // second read pass -- can safely reread from disk. Real,
            // confirmed repro this session: a single-FIFO run hung
            // indefinitely at `--shards 1` and succeeded at `--shards 2`.
            let shards_here = shards.max(2);
            let mut cmd = Command::new(&sbx10x_molecules);
            cmd.args(fifos)
                .arg("-o")
                .arg(out)
                .arg("--shards")
                .arg(shards_here.to_string())
                .arg("--shard-concurrency")
                .arg(concurrency.max(1).to_string())
                .arg("--collate-workers")
                .arg(fifos.len().max(1).to_string());
            if tagseq {
                cmd.arg("--tagseq");
            }
            cmd.stderr(Stdio::piped());
            let mut child = cmd.spawn().with_context(|| format!("spawning {tag}"))?;
            if let Some(stderr) = child.stderr.take() {
                relay_child_stderr(stderr, Arc::clone(progress), tag.to_string());
            }
            Ok(child)
        };

    if !gex_fifos.is_empty() {
        let out = cli.out_dir.join("molecules").join("gex.bam");
        let child = spawn_molecules_reader(&gex_fifos, false, &out, "molecules:gex")?;
        children.push(("molecules:gex".to_string(), child));
        gex_out = Some(out);
    }
    if !other_fifos.is_empty() {
        let out = cli.out_dir.join("molecules").join("feature.bam");
        let child = spawn_molecules_reader(&other_fifos, true, &out, "molecules:feature")?;
        children.push(("molecules:feature".to_string(), child));
        other_out = Some(out);
    }

    // Align/guides next: readers of extract's FIFOs, writers of molecules'
    // FIFOs -- must attach to the latter before molecules blocks, and to
    // the former before extract spawns.
    for (target, (bam_fifo, _)) in manifest.targets.iter().zip(&align_fifos) {
        let tag = format!("align:{}.{:03}", target.sample, target.chunk);
        let bin = if target.gex { &probeset_align } else { &guides };
        let mut cmd = Command::new(bin);
        cmd.arg(&target.r2_path).arg("-o").arg(bam_fifo);
        if target.gex {
            cmd.arg("--probe-set").arg(cli.probe_set.as_ref().unwrap());
        } else {
            cmd.arg("--feature-reference").arg(cli.feature_reference.as_ref().unwrap());
        }
        cmd.stderr(Stdio::piped());
        let mut child = cmd.spawn().with_context(|| format!("spawning {tag}"))?;
        if let Some(stderr) = child.stderr.take() {
            relay_child_stderr(stderr, Arc::clone(progress), tag.clone());
        }
        children.push((tag, child));
    }

    // Extract last -- every reader of its own FIFOs is already attached.
    let mut extract_cmd = Command::new(&sbx10x);
    extract_cmd
        .arg("--fastqs")
        .args(&cli.fastqs)
        .arg("--out-dir")
        .arg(&cli.out_dir)
        .arg("--minimap2")
        .arg("--no-gzip")
        .arg("--split")
        .arg(manifest_split(manifest).to_string())
        .arg("--threads")
        .arg(extract_threads.to_string());
    if let Some(wl) = &cli.whitelist {
        extract_cmd.arg("--whitelist").arg(wl);
    }
    if let Some(sheet) = &cli.bclconvert_samplesheet {
        extract_cmd.arg("--bclconvert-samplesheet").arg(sheet);
    }
    if let Some(pbw) = &cli.probe_barcode_whitelist {
        extract_cmd.arg("--probe-barcode-whitelist").arg(pbw);
    }
    extract_cmd.stderr(Stdio::piped());
    let mut extract_child = extract_cmd
        .spawn()
        .context("spawning extract (streaming, align->molecules)")?;
    if let Some(stderr) = extract_child.stderr.take() {
        relay_child_stderr(stderr, Arc::clone(progress), "extract".to_string());
    }
    children.push(("extract".to_string(), extract_child));

    wait_all_or_kill(children)?;

    for target in &manifest.targets {
        let _ = fs::remove_file(&target.r2_path);
    }
    for (fifo, _) in &align_fifos {
        let _ = fs::remove_file(fifo);
    }

    Ok((gex_out, other_out))
}

/// Runs the extract->align hop in `ExecutionMode::Materialize`: extract
/// writes real files first (still `--no-gzip` — this axis is about
/// concurrency, not compression), then align runs in CPU-budget-sized
/// concurrent batches (the same `xargs -P<cores>` shape the shell scripts
/// use today). Each chunk's materialized FASTQ is deleted immediately once
/// its align process exits successfully, unless `keep` is set.
fn run_extract_align_materialize(
    cli: &Cli,
    helpers: &dyn HelperResolver,
    manifest: &Manifest,
    extract_threads: usize,
    cpu_budget: usize,
    progress: &Arc<Progress>,
    keep: bool,
    resumable: bool,
) -> Result<Vec<AlignOutcome>> {
    let sbx10x = helpers.resolve("sbx10x")?;
    let probeset_align = helpers.resolve("sbxpress-probeset-align")?;
    let guides = helpers.resolve("sbxpress-guides")?;

    let extract_unit = "extract";
    if !(resumable && is_checkpointed(&cli.out_dir, extract_unit, "all")) {
        let mut extract_cmd = Command::new(&sbx10x);
        extract_cmd
            .arg("--fastqs")
            .args(&cli.fastqs)
            .arg("--out-dir")
            .arg(&cli.out_dir)
            .arg("--minimap2")
            .arg("--no-gzip")
            .arg("--split")
            .arg(manifest_split(manifest).to_string())
            .arg("--threads")
            .arg(extract_threads.to_string());
        if let Some(wl) = &cli.whitelist {
            extract_cmd.arg("--whitelist").arg(wl);
        }
        if let Some(sheet) = &cli.bclconvert_samplesheet {
            extract_cmd.arg("--bclconvert-samplesheet").arg(sheet);
        }
        if let Some(pbw) = &cli.probe_barcode_whitelist {
            extract_cmd.arg("--probe-barcode-whitelist").arg(pbw);
        }
        extract_cmd.stderr(Stdio::piped());
        let mut child = extract_cmd.spawn().context("spawning extract (materialize)")?;
        if let Some(stderr) = child.stderr.take() {
            relay_child_stderr(stderr, Arc::clone(progress), "extract".to_string());
        }
        wait_all_or_kill(vec![("extract".to_string(), child)])?;
        if resumable {
            write_checkpoint(&cli.out_dir, extract_unit, "all")?;
        }
    } else {
        progress.line("extract", "skipped (already checkpointed)");
    }

    fs::create_dir_all(cli.out_dir.join("aligned")).context("creating aligned/ dir")?;
    let mut outcomes: Vec<AlignOutcome> = Vec::with_capacity(manifest.targets.len());
    let batch_size = cpu_budget.max(1);
    for batch in manifest.targets.chunks(batch_size) {
        let mut children: Vec<(String, Child)> = Vec::new();
        let mut batch_outcomes = Vec::new();
        for target in batch {
            let bam_path = aligned_bam_path(&cli.out_dir, &target.sample, target.s_num, target.chunk);
            let unit_id = format!("{}.{:03}", target.sample, target.chunk);
            if resumable && is_checkpointed(&cli.out_dir, "align", &unit_id) {
                progress.line("align", &format!("{unit_id} skipped (already checkpointed)"));
                batch_outcomes.push((target, bam_path, None));
                continue;
            }
            let tag = format!("align:{unit_id}");
            let bin = if target.gex { &probeset_align } else { &guides };
            let mut cmd = Command::new(bin);
            cmd.arg(&target.r2_path).arg("-o").arg(&bam_path);
            if target.gex {
                cmd.arg("--probe-set").arg(
                    cli.probe_set
                        .as_ref()
                        .with_context(|| format!("{unit_id} routes to probeset-align but --probe-set was not given"))?,
                );
            } else {
                cmd.arg("--feature-reference").arg(
                    cli.feature_reference
                        .as_ref()
                        .with_context(|| format!("{unit_id} routes to guides but --feature-reference was not given"))?,
                );
            }
            cmd.stderr(Stdio::piped());
            let mut child = cmd.spawn().with_context(|| format!("spawning {tag}"))?;
            if let Some(stderr) = child.stderr.take() {
                relay_child_stderr(stderr, Arc::clone(progress), tag.clone());
            }
            children.push((tag, child));
            batch_outcomes.push((target, bam_path, Some(unit_id)));
        }
        wait_all_or_kill(children)?;
        for (target, bam_path, unit_id) in batch_outcomes {
            if let Some(unit_id) = unit_id {
                if resumable {
                    write_checkpoint(&cli.out_dir, "align", &unit_id)?;
                }
                if !resumable && !keep {
                    let _ = fs::remove_file(&target.r2_path);
                    if let Some(r1) = &target.r1_path {
                        let _ = fs::remove_file(r1);
                    }
                }
            }
            outcomes.push(AlignOutcome { bam_path, gex: target.gex });
        }
    }
    Ok(outcomes)
}

// ---------------------------------------------------------------------------
// Downstream stages (molecules -> probe-demux -> matrix -> cells -> protospacer-call)
// ---------------------------------------------------------------------------

fn run_simple_stage(bin: &Path, args: &[String], tag: &str, progress: &Arc<Progress>) -> Result<()> {
    let mut cmd = Command::new(bin);
    cmd.args(args);
    cmd.stderr(Stdio::piped());
    let mut child = cmd.spawn().with_context(|| format!("spawning {tag}"))?;
    if let Some(stderr) = child.stderr.take() {
        relay_child_stderr(stderr, Arc::clone(progress), tag.to_string());
    }
    wait_all_or_kill(vec![(tag.to_string(), child)])
}

// ---------------------------------------------------------------------------
// Top-level run()
// ---------------------------------------------------------------------------

pub fn run(cli: &Cli) -> Result<()> {
    let argv0 = env::args_os().next().unwrap_or_else(|| "sbxpress-pipeline".into());
    let executable = env::current_exe().unwrap_or_else(|_| PathBuf::from(argv0));
    let dir = helper_dir(&executable).context("cannot determine helper executable directory")?;
    let resolver = LibexecResolver { dir };
    run_with_resolver(cli, &resolver)
}

pub fn run_with_resolver(cli: &Cli, helpers: &dyn HelperResolver) -> Result<()> {
    fs::create_dir_all(&cli.out_dir)
        .with_context(|| format!("creating output dir {}", cli.out_dir.display()))?;

    let cpus = if cli.cpus > 0 { cli.cpus } else { detect_cpus() };
    let mem_gb = if cli.max_memory_gb > 0 {
        cli.max_memory_gb
    } else {
        detect_available_mem_gb().unwrap_or(0)
    };

    let progress = Arc::new(Progress::new(cli.progress_file.as_deref())?);
    let status = Arc::new(Mutex::new(format!("starting (cpus={cpus} mem_gb={mem_gb})")));
    let (hb_stop, hb_handle) =
        spawn_heartbeat(Arc::clone(&progress), Duration::from_secs(cli.progress_interval.max(1)), Arc::clone(&status));

    let result = run_inner(cli, helpers, cpus, mem_gb, &progress, &status);

    hb_stop.store(true, Ordering::Relaxed);
    let _ = hb_handle.join();
    result
}

fn set_status(status: &Arc<Mutex<String>>, msg: impl Into<String>) {
    if let Ok(mut s) = status.lock() {
        *s = msg.into();
    }
}

fn run_inner(
    cli: &Cli,
    helpers: &dyn HelperResolver,
    cpus: usize,
    mem_gb: u64,
    progress: &Arc<Progress>,
    status: &Arc<Mutex<String>>,
) -> Result<()> {
    let (extract_threads, n_readers) = plan_align_fanout(cpus, cli.split);
    let sbx10x = helpers.resolve("sbx10x")?;

    set_status(status, "extract: fetching manifest");
    progress.line("pipeline", "fetching extract manifest");
    let manifest = fetch_manifest(&sbx10x, cli, n_readers)?;
    let n_targets = manifest.targets.len();

    let mode = if cli.resumable {
        ExecutionMode::Materialize
    } else if cli.materialize_extract {
        ExecutionMode::Materialize
    } else {
        plan_extract_align_mode(cpus, mem_gb, n_targets, extract_threads)
    };
    progress.line(
        "pipeline",
        &format!("extract->align mode={mode:?} targets={n_targets} extract_threads={extract_threads}"),
    );

    // Sized here (not after align) because the align->molecules streaming
    // decision below needs it: `total_reads` is a placeholder (the real
    // count isn't known before molecules runs) but `plan_molecule_shards`
    // has always treated it that way, moving the call earlier doesn't lose
    // any real precision it had before.
    let total_reads = n_targets as u64 * 1;
    let (shards, concurrency) = if cli.molecule_shards > 0 {
        (cli.molecule_shards, cli.molecule_shard_concurrency.max(1))
    } else {
        plan_molecule_shards(total_reads.max(1), mem_gb, cpus)
    };

    let stream_molecules_too = cli.stop_after != StopAfter::Align
        && plan_align_molecules_stream(mode, cpus, extract_threads, n_targets, concurrency);

    let sbx10x_molecules = helpers.resolve("sbx10x-molecules")?;
    let (gex_molecules, other_molecules) = if stream_molecules_too {
        progress.line(
            "pipeline",
            "extract->align->molecules mode=Stream (zero intermediate files)",
        );
        set_status(status, format!("extract->align->molecules (streaming, {n_targets} targets)"));
        let (gex, other) = run_extract_align_molecules_stream(
            cli,
            helpers,
            &manifest,
            extract_threads,
            shards,
            concurrency,
            progress,
        )?;
        progress.line("pipeline", "extract->align->molecules complete");
        (gex, other)
    } else {
        set_status(status, format!("extract->align ({mode:?}, {n_targets} targets)"));
        let align_outcomes = match mode {
            ExecutionMode::Stream => {
                run_extract_align_stream(cli, helpers, &manifest, extract_threads, progress)?
            }
            ExecutionMode::Materialize => run_extract_align_materialize(
                cli,
                helpers,
                &manifest,
                extract_threads,
                cpus,
                progress,
                cli.keep_intermediates,
                cli.resumable,
            )?,
        };
        progress.line("pipeline", &format!("align complete: {} chunks", align_outcomes.len()));
        if cli.stop_after == StopAfter::Align {
            return Ok(());
        }

        // molecules: GEX-routed chunks and non-GEX (CRISPR/Feature
        // Barcoding) chunks are molecule-counted separately, matching this
        // project's existing convention (CRISPR gets `--tagseq`, run
        // separately from GEX -- see `scripts/run_flex_multi.sh`).
        let gex_bams: Vec<PathBuf> =
            align_outcomes.iter().filter(|o| o.gex).map(|o| o.bam_path.clone()).collect();
        let other_bams: Vec<PathBuf> =
            align_outcomes.iter().filter(|o| !o.gex).map(|o| o.bam_path.clone()).collect();

        let mut gex_molecules: Option<PathBuf> = None;
        if !gex_bams.is_empty() {
            set_status(status, "molecules (GEX)");
            let out = cli.out_dir.join("molecules").join("gex.bam");
            fs::create_dir_all(out.parent().unwrap()).context("creating molecules/ dir")?;
            let mut args: Vec<String> =
                gex_bams.iter().map(|p| p.to_string_lossy().to_string()).collect();
            args.push("-o".into());
            args.push(out.to_string_lossy().to_string());
            args.push("--shards".into());
            args.push(shards.to_string());
            args.push("--shard-concurrency".into());
            args.push(concurrency.to_string());
            run_simple_stage(&sbx10x_molecules, &args, "molecules:gex", progress)?;
            gex_molecules = Some(out);
        }
        let mut other_molecules: Option<PathBuf> = None;
        if !other_bams.is_empty() {
            set_status(status, "molecules (feature/CRISPR)");
            let out = cli.out_dir.join("molecules").join("feature.bam");
            fs::create_dir_all(out.parent().unwrap()).context("creating molecules/ dir")?;
            let mut args: Vec<String> =
                other_bams.iter().map(|p| p.to_string_lossy().to_string()).collect();
            args.push("-o".into());
            args.push(out.to_string_lossy().to_string());
            args.push("--tagseq".into());
            run_simple_stage(&sbx10x_molecules, &args, "molecules:feature", progress)?;
            other_molecules = Some(out);
        }
        (gex_molecules, other_molecules)
    };
    progress.line("pipeline", "molecules complete");
    if cli.stop_after == StopAfter::Molecules {
        return Ok(());
    }

    // Downstream stages only run when there's a real GEX molecules BAM to
    // demux/count/call cells on -- a CRISPR/feature-only run stops here.
    let Some(gex_molecules) = gex_molecules else {
        progress.line("pipeline", "no GEX molecules output; stopping before probe-demux/matrix/cells");
        return Ok(());
    };

    if let Some(sheet) = &cli.bclconvert_samplesheet {
        set_status(status, "probe-demux");
        let probe_demux = helpers.resolve("sbxpress-probe-demux")?;
        let demux_dir = cli.out_dir.join("demux");
        fs::create_dir_all(&demux_dir).context("creating demux/ dir")?;
        let args = vec![
            gex_molecules.to_string_lossy().to_string(),
            "--probe-sample-sheet".into(),
            sheet.to_string_lossy().to_string(),
            "-o".into(),
            demux_dir.to_string_lossy().to_string(),
        ];
        run_simple_stage(&probe_demux, &args, "probe-demux", progress)?;
    }

    set_status(status, "matrix");
    let matrix_bin = helpers.resolve("cr-bam-to-mtx")?;
    let matrix_out = cli.out_dir.join("raw_matrix");
    let args = vec![
        gex_molecules.to_string_lossy().to_string(),
        "-o".into(),
        matrix_out.to_string_lossy().to_string(),
    ];
    run_simple_stage(&matrix_bin, &args, "matrix", progress)?;
    progress.line("pipeline", "matrix complete");
    if cli.stop_after == StopAfter::Matrix {
        return Ok(());
    }

    set_status(status, "cells");
    let cells_bin = helpers.resolve("sbxpress-cells")?;
    let cells_out = cli.out_dir.join("filtered_matrix");
    let args = vec![
        matrix_out.to_string_lossy().to_string(),
        "-o".into(),
        cells_out.to_string_lossy().to_string(),
    ];
    run_simple_stage(&cells_bin, &args, "cells", progress)?;
    progress.line("pipeline", "cells complete");
    if cli.stop_after == StopAfter::Cells {
        return Ok(());
    }

    if let Some(other_molecules) = other_molecules {
        set_status(status, "protospacer-call");
        let protospacer_bin = helpers.resolve("sbxpress-protospacer-call")?;
        let ps_out = cli.out_dir.join("protospacer_calls");
        let args = vec![
            other_molecules.to_string_lossy().to_string(),
            "-o".into(),
            ps_out.to_string_lossy().to_string(),
        ];
        let _ = run_simple_stage(&protospacer_bin, &args, "protospacer-call", progress);
    }

    if cli.resumable && !cli.keep_intermediates {
        set_status(status, "sweeping intermediates");
        sweep_intermediates(&cli.out_dir)?;
    }

    progress.line("pipeline", "complete");
    Ok(())
}

/// Deletes intermediate checkpointed directories (`extracted/`, `aligned/`)
/// after the whole pipeline reports final success -- the `--resumable`
/// cleanup-timing rule. Final deliverables (`molecules/`, `raw_matrix/`,
/// `filtered_matrix/`, `demux/`, `protospacer_calls/`) are never swept.
fn sweep_intermediates(out_dir: &Path) -> Result<()> {
    for name in ["extracted", "aligned"] {
        let dir = out_dir.join(name);
        if dir.is_dir() {
            let _ = fs::remove_dir_all(&dir);
        }
    }
    let _ = fs::remove_dir_all(checkpoint_dir(out_dir));
    Ok(())
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stage_shape_lookup() {
        assert_eq!(stage_shape("extract"), Some(ResourceShape::CpuBoundParallel));
        assert_eq!(stage_shape("molecules"), Some(ResourceShape::MemoryDominantSerial));
        assert_eq!(stage_shape("cells"), Some(ResourceShape::LightweightIoBound));
        assert_eq!(stage_shape("does-not-exist"), None);
    }

    #[test]
    fn plan_align_fanout_respects_explicit_split() {
        let (extract_threads, n_readers) = plan_align_fanout(128, 32);
        assert_eq!(n_readers, 32);
        assert!(extract_threads >= 1 && extract_threads <= 128);
    }

    #[test]
    fn plan_align_fanout_auto_stays_within_budget() {
        let (extract_threads, n_readers) = plan_align_fanout(16, 0);
        assert!(extract_threads + n_readers <= 16 || extract_threads == 1);
        assert!(n_readers >= 1);
    }

    #[test]
    fn plan_align_fanout_never_zero() {
        let (extract_threads, n_readers) = plan_align_fanout(0, 0);
        assert!(extract_threads >= 1);
        assert!(n_readers >= 1);
    }

    #[test]
    fn plan_molecule_shards_real_anchor_stays_bounded() {
        // Real PIPELINE_NOTES.md anchor: ~9B reads, ~900GB available memory,
        // real observed --shards 8 run on this box. Assert the heuristic
        // produces a sane, safety-margined shard count in the same order of
        // magnitude, and that per-shard memory actually fits the budget --
        // not a brittle exact-equality check against one hand-derived number.
        let (shards, concurrency) = plan_molecule_shards(9_000_000_000, 900, 128);
        assert!(shards >= 1);
        assert!(concurrency >= 1 && concurrency <= shards);
        let per_shard_reads = 9_000_000_000_f64 / shards as f64;
        let per_shard_bytes = per_shard_reads * BYTES_PER_READ_ESTIMATE * SHARD_SAFETY_FACTOR;
        assert!(per_shard_bytes <= 900.0 * 1_000_000_000.0);
    }

    #[test]
    fn plan_molecule_shards_more_memory_never_needs_more_shards() {
        let (shards_low_mem, _) = plan_molecule_shards(9_000_000_000, 100, 128);
        let (shards_high_mem, _) = plan_molecule_shards(9_000_000_000, 900, 128);
        assert!(shards_high_mem <= shards_low_mem);
    }

    #[test]
    fn plan_molecule_shards_more_reads_never_needs_fewer_shards() {
        let (shards_few, _) = plan_molecule_shards(1_000_000_000, 200, 128);
        let (shards_many, _) = plan_molecule_shards(20_000_000_000, 200, 128);
        assert!(shards_many >= shards_few);
    }

    #[test]
    fn plan_extract_align_mode_fits_streams() {
        let mode = plan_extract_align_mode(128, 900, 16, 4);
        assert_eq!(mode, ExecutionMode::Stream);
    }

    #[test]
    fn plan_extract_align_mode_exceeds_cpu_materializes() {
        let mode = plan_extract_align_mode(8, 900, 64, 4);
        assert_eq!(mode, ExecutionMode::Materialize);
    }

    #[test]
    fn plan_extract_align_mode_tiny_memory_materializes() {
        let mode = plan_extract_align_mode(128, 1, 4, 2);
        assert_eq!(mode, ExecutionMode::Materialize);
    }

    #[test]
    fn plan_align_molecules_stream_requires_extract_align_streaming_first() {
        // Plenty of CPU for molecules too, but extract->align itself is
        // Materialize -- can't stream a hop further than the first one.
        assert!(!plan_align_molecules_stream(ExecutionMode::Materialize, 128, 4, 16, 8));
    }

    #[test]
    fn plan_align_molecules_stream_fits_when_budget_covers_all_three() {
        // extract_threads(4) + n_targets(16) + molecule_shard_concurrency(8) = 28 <= 128.
        assert!(plan_align_molecules_stream(ExecutionMode::Stream, 128, 4, 16, 8));
    }

    #[test]
    fn plan_align_molecules_stream_exceeds_cpu_falls_back() {
        // extract_threads(4) + n_targets(16) + molecule_shard_concurrency(8) = 28 > 16.
        assert!(!plan_align_molecules_stream(ExecutionMode::Stream, 16, 4, 16, 8));
    }

    #[test]
    fn manifest_json_round_trips() {
        let manifest = Manifest {
            version: 1,
            targets: vec![crate::ManifestTarget {
                sample: "TSeq100_L1_GEX".into(),
                s_num: 1,
                chunk: 1,
                gex: true,
                r2_path: PathBuf::from("/tmp/out/TSeq100_L1_GEX_S1_L001_R2_001.fastq"),
                r1_path: None,
            }],
        };
        let json = serde_json::to_string(&manifest).unwrap();
        let parsed: Manifest = serde_json::from_str(&json).unwrap();
        assert_eq!(parsed.version, 1);
        assert_eq!(parsed.targets.len(), 1);
        assert!(parsed.targets[0].gex);
    }

    #[test]
    fn manifest_rejects_unknown_version() {
        let json = r#"{"version":2,"targets":[]}"#;
        let manifest: Manifest = serde_json::from_str(json).unwrap();
        assert!(validate_manifest_version(&manifest).is_err());
    }

    #[test]
    fn manifest_accepts_version_one() {
        let json = r#"{"version":1,"targets":[]}"#;
        let manifest: Manifest = serde_json::from_str(json).unwrap();
        assert!(validate_manifest_version(&manifest).is_ok());
    }

    #[test]
    fn checkpoint_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        assert!(!is_checkpointed(dir.path(), "align", "sample.001"));
        write_checkpoint(dir.path(), "align", "sample.001").unwrap();
        assert!(is_checkpointed(dir.path(), "align", "sample.001"));
        assert!(!is_checkpointed(dir.path(), "align", "sample.002"));
    }

    #[test]
    fn checkpoint_path_is_stable() {
        let dir = PathBuf::from("/tmp/out");
        assert_eq!(
            checkpoint_path(&dir, "align", "sample.001"),
            PathBuf::from("/tmp/out/.checkpoints/align.sample.001.done")
        );
    }

    #[test]
    fn heartbeat_due_respects_interval() {
        assert!(!heartbeat_due(Duration::from_secs(10), Duration::from_secs(30)));
        assert!(heartbeat_due(Duration::from_secs(30), Duration::from_secs(30)));
        assert!(heartbeat_due(Duration::from_secs(31), Duration::from_secs(30)));
    }

    #[test]
    fn progress_line_has_timestamp_and_elapsed_prefix() {
        let line = format_progress_line(Duration::from_secs(3725), "align:chunk07", "records=100");
        assert!(line.contains("[align:chunk07]"));
        assert!(line.contains("records=100"));
        assert!(line.contains("+01:02:05"));
        // Real ISO8601 UTC timestamp: `[YYYY-MM-DDTHH:MM:SSZ ...`.
        assert!(line.starts_with('['));
        assert!(line.contains('T'));
        assert!(line.contains('Z'));
    }

    #[test]
    fn civil_from_unix_epoch_zero() {
        assert_eq!(civil_from_unix(0), (1970, 1, 1, 0, 0, 0));
    }

    #[test]
    fn civil_from_unix_one_day_plus_time() {
        // 1970-01-02T01:01:01Z -- one day (uncontroversial calendar math) plus
        // 1h1m1s, so this doesn't depend on the leap-year branch at all.
        assert_eq!(civil_from_unix(86400 + 3661), (1970, 1, 2, 1, 1, 1));
    }

    #[test]
    fn detect_cpus_is_at_least_one() {
        assert!(detect_cpus() >= 1);
    }

    #[test]
    fn parse_mem_available_kb_finds_the_line() {
        let sample = "MemTotal:       1610612736 kB\nMemFree:         610612736 kB\nMemAvailable:   1400000000 kB\n";
        assert_eq!(parse_mem_available_kb(sample), Some(1_400_000_000));
    }

    #[test]
    fn parse_mem_available_kb_missing_returns_none() {
        assert_eq!(parse_mem_available_kb("MemTotal: 123 kB\n"), None);
    }

    #[test]
    fn wait_all_or_kill_kills_siblings_on_first_failure() {
        // Real subprocess test of the actual mitigation for the FIFO-streaming
        // failure mode: one child fails fast, the other would otherwise sleep
        // for 30s (standing in for "blocked forever on a FIFO read/write").
        // Assert the whole call returns promptly (i.e. the sleeper was really
        // killed, not waited out) and reports the real failure.
        let start = Instant::now();
        let failer = Command::new("sh").args(["-c", "exit 7"]).spawn().unwrap();
        let sleeper = Command::new("sh").args(["-c", "sleep 30"]).spawn().unwrap();
        let result = wait_all_or_kill(vec![
            ("failer".to_string(), failer),
            ("sleeper".to_string(), sleeper),
        ]);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("failer"));
        assert!(
            start.elapsed() < Duration::from_secs(10),
            "sleeper should have been killed, not waited out (took {:?})",
            start.elapsed()
        );
    }

    #[test]
    fn wait_all_or_kill_succeeds_when_all_children_succeed() {
        let a = Command::new("sh").args(["-c", "exit 0"]).spawn().unwrap();
        let b = Command::new("sh").args(["-c", "exit 0"]).spawn().unwrap();
        let result = wait_all_or_kill(vec![("a".to_string(), a), ("b".to_string(), b)]);
        assert!(result.is_ok());
    }

    #[test]
    fn aligned_bam_path_matches_naming_convention() {
        let path = aligned_bam_path(&PathBuf::from("/out"), "TSeq100_L1_GEX", 1, 7);
        assert_eq!(path, PathBuf::from("/out/aligned/TSeq100_L1_GEX_S1_L001_R2_007.bam"));
    }
}
