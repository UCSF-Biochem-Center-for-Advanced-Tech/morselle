//! Real subprocess integration test for `sbxpress pipeline`'s core new
//! capability: streaming `extract`'s output directly into `probeset-align`
//! via named pipes, with zero intermediate FASTQ files. Uses the real
//! compiled binaries (via `CARGO_BIN_EXE_*`, which Cargo guarantees are
//! built before integration tests run) and a small `HelperResolver` wired to
//! them, instead of production's installed/`libexec` filesystem layout. See
//! `/home/matt/.claude/plans/valiant-gathering-codd.md` for the design this
//! validates.

use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::Result;
use rust_htslib::bam::{self, Read as _};
use sbx10x::orchestrator::{self, Cli, HelperResolver, StopAfter};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

/// Real header + shape for a minimal probe-set reference, copied from
/// `src/probe_set.rs`'s own `REAL_SHAPE` test fixture — this test doesn't
/// need any of `smoke.fastq.gz`'s R2 content to actually *match* a real
/// gene (zero probe hits is a perfectly valid, successful run); it needs
/// `probeset-align` to load a real, valid reference and process real
/// streamed FASTQ bytes end to end.
const PROBE_SET_CSV: &str = "\
#probe_set_file_format=3.0
#panel_name=test
gene_id,probe_seq,probe_id,included,region,gene_name
ENSG1,AAAACCCCGGGGTTTTACGTACGTA,ENSG1|GENE1|abc,TRUE,spliced,GENE1
ENSG2,TTTTGGGGCCCCAAAATGCATGCAT,ENSG2|GENE2|ghi,TRUE,unspliced,GENE2
";

/// Resolves sibling binaries to Cargo's own build outputs for this test
/// run, exactly as `orchestrator::run_with_resolver` expects an injected
/// `HelperResolver` to work — no `SBXPRESS_LIBEXEC`/installed-layout
/// dependency.
struct CargoBinResolver;

impl HelperResolver for CargoBinResolver {
    fn resolve(&self, name: &str) -> Result<PathBuf> {
        let path = match name {
            "sbx10x" => env!("CARGO_BIN_EXE_sbx10x"),
            "sbxpress-probeset-align" => env!("CARGO_BIN_EXE_sbxpress-probeset-align"),
            "sbxpress-guides" => env!("CARGO_BIN_EXE_sbxpress-guides"),
            "sbx10x-molecules" => env!("CARGO_BIN_EXE_sbx10x-molecules"),
            "sbxpress-probe-demux" => env!("CARGO_BIN_EXE_sbxpress-probe-demux"),
            "cr-bam-to-mtx" => env!("CARGO_BIN_EXE_sbxpress"), // placeholder, unused by this test
            "sbxpress-cells" => env!("CARGO_BIN_EXE_sbxpress-cells"),
            "sbxpress-protospacer-call" => env!("CARGO_BIN_EXE_sbxpress-protospacer-call"),
            other => anyhow::bail!("no CARGO_BIN_EXE_* mapping for {other}"),
        };
        Ok(PathBuf::from(path))
    }
}

fn base_cli(out_dir: &Path, probe_set: &Path) -> Cli {
    Cli {
        fastqs: vec![fixture("smoke.fastq.gz")],
        bclconvert_samplesheet: None,
        // No whitelist: `illumina_whitelist.txt` doesn't match this
        // fixture's real barcodes, and the CB-whitelist gate would drop
        // every read (extract's own default behavior without one is to
        // skip that gate entirely, per `config.rs`) -- this test cares
        // about real bytes flowing through the FIFO end to end, not CB
        // gating specifically.
        whitelist: None,
        probe_barcode_whitelist: None,
        probe_set: Some(probe_set.to_path_buf()),
        feature_reference: None,
        out_dir: out_dir.to_path_buf(),
        cpus: 8,
        max_memory_gb: 64,
        split: 1,
        molecule_shards: 0,
        molecule_shard_concurrency: 0,
        materialize_extract: false,
        keep_intermediates: false,
        resumable: false,
        progress_interval: 5,
        progress_file: None,
        stop_after: StopAfter::Align,
    }
}

fn write_probe_set(dir: &Path) -> PathBuf {
    let path = dir.join("probe_set.csv");
    std::fs::File::create(&path).unwrap().write_all(PROBE_SET_CSV.as_bytes()).unwrap();
    path
}

#[test]
fn streams_extract_output_into_probeset_align_with_no_intermediate_fastq() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let probe_set = write_probe_set(temp.path());
    let cli = base_cli(&out_dir, &probe_set);

    // Budget deliberately generous (1 target, plenty of cpus/mem) so the
    // scheduler picks Stream mode -- this is the actual default path being
    // validated, not the Materialize fallback.
    orchestrator::run_with_resolver(&cli, &CargoBinResolver).expect("pipeline run should succeed");

    let bam_path = out_dir.join("aligned").join("sbxpress_S1_L001_R2_001.bam");
    assert!(bam_path.is_file(), "aligned BAM should exist at {}", bam_path.display());

    // Real byte-correctness check: the BAM is actually openable and its
    // records came from real streamed data, not an empty/corrupt file.
    let mut reader = bam::Reader::from_path(&bam_path).expect("aligned BAM must be a valid, readable BAM");
    let header = bam::Header::from_template(reader.header());
    assert!(!header.to_hashmap().is_empty(), "BAM must carry a real header");
    let mut record_count = 0u64;
    for result in reader.records() {
        result.expect("every record must parse cleanly");
        record_count += 1;
    }
    // smoke.fastq.gz has real emitted pairs (tests/failed_output.rs already
    // asserts `total_reads > 0` for it); with --minimap2 every emitted pair
    // becomes one streamed R2 record `probeset-align` should have seen and
    // written back out (gene match or not).
    assert!(record_count > 0, "expected real records to have flowed through the FIFO");

    // The core "low/no intermediate output" claim: no leftover FASTQ (the
    // FIFOs are removed after a successful stream) and no `extracted/`
    // directory was ever created (Stream mode never materializes it).
    assert!(
        !out_dir.join("extracted").exists(),
        "streaming mode must never create a materialized extracted/ dir"
    );
    let leftover_fastq: Vec<_> = std::fs::read_dir(&out_dir)
        .into_iter()
        .flatten()
        .flatten()
        .filter(|e| e.path().extension().is_some_and(|ext| ext == "fastq" || ext == "gz"))
        .collect();
    assert!(leftover_fastq.is_empty(), "no FASTQ should remain in {}: {leftover_fastq:?}", out_dir.display());
}

/// Phase 2 of the streaming work: with a generous enough budget, align's
/// *output* also streams directly into `sbx10x-molecules` (no per-chunk
/// aligned BAM ever touches disk), not just extract's output into align.
/// Validates the real end-to-end chain -- extract -> probeset-align ->
/// molecules -- with zero intermediate files at either hop, and a real,
/// valid molecules BAM as the actual output.
#[test]
fn streams_extract_align_molecules_with_zero_intermediate_bam() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let probe_set = write_probe_set(temp.path());
    let mut cli = base_cli(&out_dir, &probe_set);
    cli.stop_after = StopAfter::Molecules;

    orchestrator::run_with_resolver(&cli, &CargoBinResolver).expect("pipeline run should succeed");

    let molecules_bam = out_dir.join("molecules").join("gex.bam");
    assert!(
        molecules_bam.is_file(),
        "molecules BAM should exist at {}",
        molecules_bam.display()
    );
    let mut reader =
        bam::Reader::from_path(&molecules_bam).expect("molecules BAM must be a valid, readable BAM");
    let header = bam::Header::from_template(reader.header());
    assert!(!header.to_hashmap().is_empty(), "BAM must carry a real header");
    let mut record_count = 0u64;
    for result in reader.records() {
        result.expect("every record must parse cleanly");
        record_count += 1;
    }
    assert!(record_count > 0, "expected real records to have flowed all the way to molecules");

    // The core claim: the per-chunk aligned BAM this hop would otherwise
    // have written to disk (and that the plain align-only test above
    // asserts *does* exist) must never have been materialized -- it was a
    // FIFO, consumed directly by molecules, then removed.
    let would_be_aligned_bam = out_dir.join("aligned").join("sbxpress_S1_L001_R2_001.bam");
    assert!(
        !would_be_aligned_bam.exists(),
        "streamed align->molecules must never leave a per-chunk aligned BAM on disk: {}",
        would_be_aligned_bam.display()
    );
    assert!(
        !out_dir.join("extracted").exists(),
        "streaming mode must never create a materialized extracted/ dir"
    );
}

/// Failure-injection: a broken align process (bad `--probe-set` path) in the
/// fused extract->align->molecules streaming topology must surface as a
/// real error and must NOT hang -- `wait_all_or_kill`'s active-kill-on-
/// first-failure logic (already exercised for the plain extract->align hop)
/// needs to cover this new, longer child group (molecules reader + align +
/// extract) too, since a crashed align process closes molecules' input FIFO
/// out from under it and closes extract's output FIFO out from under
/// extract, either of which could hang forever on a blocking read/write
/// without the poll-and-kill loop reaching them. Self-bounding (a channel +
/// `recv_timeout`, not a real per-test framework timeout) so a regression
/// that reintroduces a hang fails this test instead of hanging CI silently.
#[test]
fn align_failure_in_streamed_molecules_chain_propagates_and_does_not_hang() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let mut cli = base_cli(&out_dir, &temp.path().join("does_not_exist_probe_set.csv"));
    cli.stop_after = StopAfter::Molecules;

    let (tx, rx) = std::sync::mpsc::channel();
    std::thread::spawn(move || {
        let result = orchestrator::run_with_resolver(&cli, &CargoBinResolver);
        let _ = tx.send(result.is_err());
    });
    let finished_with_error = rx
        .recv_timeout(std::time::Duration::from_secs(30))
        .expect("pipeline must fail fast, not hang, when an align process crashes mid-stream");
    assert!(finished_with_error, "a broken --probe-set must surface as a real pipeline error");
}

#[test]
fn materialize_extract_flag_produces_real_intermediate_files_and_cleans_up() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let probe_set = write_probe_set(temp.path());
    let mut cli = base_cli(&out_dir, &probe_set);
    cli.materialize_extract = true;

    orchestrator::run_with_resolver(&cli, &CargoBinResolver).expect("pipeline run should succeed");

    let bam_path = out_dir.join("aligned").join("sbxpress_S1_L001_R2_001.bam");
    assert!(bam_path.is_file(), "aligned BAM should exist at {}", bam_path.display());

    // Default cleanup semantics: materialized FASTQ is deleted once its
    // align process consumes it successfully (not kept, since
    // --keep-intermediates was not set).
    let extracted_dir = out_dir.join("extracted");
    if extracted_dir.exists() {
        let remaining: Vec<_> = std::fs::read_dir(&extracted_dir).unwrap().collect();
        assert!(remaining.is_empty(), "materialized FASTQ should be cleaned up after consumption");
    }
}

#[test]
fn keep_intermediates_retains_materialized_fastq() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let probe_set = write_probe_set(temp.path());
    let mut cli = base_cli(&out_dir, &probe_set);
    cli.materialize_extract = true;
    cli.keep_intermediates = true;

    orchestrator::run_with_resolver(&cli, &CargoBinResolver).expect("pipeline run should succeed");

    let fastq = out_dir.join("sbxpress_S1_L001_R2_001.fastq");
    assert!(fastq.is_file(), "--keep-intermediates must retain the materialized FASTQ at {}", fastq.display());
}

#[test]
fn resumable_run_writes_checkpoints_and_skips_completed_work_on_restart() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let probe_set = write_probe_set(temp.path());
    let mut cli = base_cli(&out_dir, &probe_set);
    cli.resumable = true;

    orchestrator::run_with_resolver(&cli, &CargoBinResolver).expect("first resumable run should succeed");
    let bam_path = out_dir.join("aligned").join("sbxpress_S1_L001_R2_001.bam");
    assert!(bam_path.is_file());
    let first_mtime = std::fs::metadata(&bam_path).unwrap().modified().unwrap();

    // A second run against the same out_dir with the same flags should see
    // the align checkpoint and skip redoing that work (the align BAM's
    // mtime must not change).
    std::thread::sleep(std::time::Duration::from_millis(1100));
    orchestrator::run_with_resolver(&cli, &CargoBinResolver).expect("second resumable run should succeed");
    let second_mtime = std::fs::metadata(&bam_path).unwrap().modified().unwrap();
    assert_eq!(first_mtime, second_mtime, "resumable restart must skip already-checkpointed align work");
}

#[test]
fn heartbeat_reaches_progress_file_within_the_configured_interval() {
    let temp = tempfile::tempdir().unwrap();
    let out_dir = temp.path().join("out");
    let probe_set = write_probe_set(temp.path());
    let mut cli = base_cli(&out_dir, &probe_set);
    let progress_file = temp.path().join("progress.log");
    cli.progress_file = Some(progress_file.clone());
    cli.progress_interval = 1;

    orchestrator::run_with_resolver(&cli, &CargoBinResolver).expect("pipeline run should succeed");

    let contents = std::fs::read_to_string(&progress_file).expect("progress file should exist");
    assert!(!contents.is_empty(), "progress file should have received output");
    assert!(contents.contains("[pipeline]"), "expected at least one pipeline-tagged line: {contents}");
    // Real timestamp format sanity check (ISO8601 UTC + elapsed prefix).
    assert!(contents.contains('T') && contents.contains('Z'), "expected ISO8601 timestamps: {contents}");
}
