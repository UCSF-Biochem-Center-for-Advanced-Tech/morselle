# Changelog

All notable changes to `sbxpress` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.1] — 2026-07-13

### Fixed
- Unique Hamming-1 cell barcodes are corrected to the canonical whitelist
  sequence before minimap2/Cell Ranger `CB` emission; ambiguous Hamming-1
  matches are rejected rather than assigned arbitrarily.

### Added
- `sbxpress matrix --cell-barcodes` reads a plain/gzipped Cell Ranger whitelist,
  applies a configurable gem-group suffix, and retains every zero-count cell column.
- `sbxpress annotate --genes-mtx` creates or verifies a canonical filtered gene
  manifest; `sbxpress matrix --features` now enforces it strictly for stable rows.
- Formal `resources/cellranger_barcodes/` staging and installation of the full
  available 10x barcode-resource collection under `share/sbxpress/barcodes/`.
- Idempotent `uninstall.sh` with custom-prefix, dry-run, data-preservation, and
  optional minimap2-companion removal controls.
- Unified `sbxpress` CLI with `extract`, `minimap`, `annotate`, `molecules`,
  `audit`, and `matrix` subcommands. `minimap` transparently wraps minimap2;
  processing engines are installed privately under `libexec/sbxpress`.
- `sbx10x-molecules`: two-pass, coordinate-order-preserving gene-aware UMI
  correction, read-support collision resolution, corrected `UB`, Cell Ranger
  `xf=17/25` eligibility/representative tagging, and `XV`/`XW` collision audit tags.
- `--processing-tags` emits `XO`, `XB`, `XI`, `XA`, and `XL` per-read SBX
  provenance in R2 FASTQ comments for minimap2 `-y` propagation. The
  `--minimap2` preset enables it automatically.
- `sbx10x-annotate`: strand-aware GTF feature assignment that preserves existing
  barcode/UMI tags and adds Cell Ranger-compatible `GX`, `GN`, `TX`, `AN`, and
  `RE` tags to coordinate-sorted minimap2 BAMs, with default 10x GENCODE filtering,
  optional intron inclusion, collision-degree tags (`XC`/`XT`), non-destructive
  terminal-exon preference tags (`ZG`/`ZT`/`XP`) behind an explicit
  `--three-prime` assay flag, BAM/SAM stdin and BAM stdout streaming, and JSON
  run metrics.
- `--minimap2` R2-only output preset. It puts raw/corrected Cell Barcode and UMI
  tags (`CR`, `CB`, `UR`, provisional `UB`, `CY`, `UY`) in FASTQ comments for
  minimap2 `-y` propagation to SAM/BAM, and records `output_mode` in run metrics.
- The installer now builds and installs the `k8` runtime and a version-matched
  `paftools.js` alongside minimap2, enabling annotated-junction BED preparation.

## [0.3.0] — 2026-07-12

### Changed
- **`--realign` (single-indel CB rescue) is now ON by default** (was off). Analysis of
  run_5p showed the residual 3′-end CB "errors" are largely misaligned barcodes — a
  homopolymer indel frame-shifts the 16 bp grab so a random base lands in the last CB
  slot (transition bias collapses 81%→43% toward pos 15; insertion rate rises 5×).
  These are recovered by single-indel realignment, gated on the structural CB+UMI span
  signal (±1). Measured ~+1% of reads recovered at high precision: genuine single-indels
  are correct-or-abstain in simulation (0% wrong), and recoveries land on real-cell
  barcodes ~4–10× above chance. The flag is now settable (`--realign true|false`).

### Added
- `cr-bam-to-mtx`: reconstructs Cell Ranger raw GEX matrices from corrected BAM
  molecule tags, including Cell Ranger-compatible MEX coordinate ordering, the
  `outs/raw_feature_bc_matrix/` layout, and its CSC HDF5 neighbor
  `outs/raw_feature_bc_matrix.h5`.
- Structured `metrics_summary.json`, a `RUN SUMMARY`, and fixture validation for
  `cr-bam-to-mtx`.
- The installer now builds, installs, and verifies both command-line apps.
- The installer installs and verifies `minimap2` and `samtools` through
  `apt-get` on Ubuntu/Debian and Homebrew on macOS.
- Diagnostic examples: `lenient_precision` and `realign_precision` (recall/precision of
  the within-2 and single-indel CB rescues via error injection), `cb_error_profile`
  (CB substitution position/type + frame tests), and `realign_concordance` (real-data
  precision proxy — do realign recoveries land on abundant real-cell barcodes?).
- Finding (no code change): lenient within-2 rescue stays OFF — ~51% precision on the
  dense 736K ARC whitelist (within-2 neighborhoods are usually ambiguous).

## [0.2.0] — 2026-07-07

### Added
- `bench_whitelist` example: microbenchmark comparing the whitelist membership
  backends (sorted `Vec` vs `FxHashSet`) across sizes 96 → ~2M on hit/fail streams.
- `--cellranger` compatibility preset that applies the settings Cell Ranger
  requires (currently `--max-r2-len 649`); explicit flags override it.
- `--reader-threads` (default `0` = auto) selecting the I/O strategy.
- Pipelined I/O path for runs with few input files: one reader thread per file
  decouples serial gzip decode from processing and streams record batches over a
  bounded channel to the full worker pool, so a single large file now uses every
  core. Runs with many files keep the one-file-per-thread path. Auto-selected by
  file count; mate-order (Nth R1 == Nth R2) is preserved in both paths.
- README with overview, pipeline description, install, usage, options, and
  architecture.

### Changed
- **R2 is now emitted reverse-complemented** (with reversed quality) to match the
  Illumina paired-end convention, fixing the ~46% "Reads Mapped Antisense to Gene"
  reported by Cell Ranger.
- `--max-r2-len` is applied at write time (after reverse-complementing) instead of
  during extraction, so the whole cDNA fragment is retained through adapter trimming
  and the cap keeps the bases the R2 primer would actually read.
- Default `--gzip-level` lowered from 6 to 4 (output compression is the wall-clock
  floor at scale; level 4 favors throughput for a small size increase).
- Hot-path speedups (byte-identical output): early-exit Hamming for the 3' adapter
  scan, and a zero-allocation R2 emitter that reverse-complements directly into the
  output buffer and only touches the emitted (≤ cap) bases.
- Whitelist membership now uses `FxHashSet<u32>` (O(1)) instead of a sorted `Vec`
  binary search. Benchmarked (`bench_whitelist`) at ~10× faster CB hits and ~2–5×
  faster CB-fails (the 48-neighbor path) on 736k–2M whitelists; ~2.5× the barcode
  memory (<20 MB at 2M). `Whitelist::codes()` → `codes_sorted()` (diagnostics only).

### Fixed
- Cell Ranger `quality string length is not equal to sequence length` on long reads:
  Cell Ranger hard-caps reads at 649 bp (fails at 650). Use `--cellranger` (or
  `--max-r2-len 649`) to cap R2 within the limit.

## [0.1.0] — milestone build (M0–M7)

Initial end-to-end pipeline converting raw Roche SBX long reads into per-sample 10x
R1/R2 FASTQs for Cell Ranger ARC.

### Added
- **M0** — crate scaffold: CLI, logging, module layout, FASTQ read/write passthrough.
- **M1** — per-read orientation to a canonical P5→P7 frame: homopolymer-first
  (poly-dT) anchoring with a windowed Hamming R1-primer-tail confirm, and per-run
  orientation detection/enforcement that re-anchors past false cDNA homopolymers.
- **M2** — R1 (CB+UMI) and R2 (cDNA) bracket extraction with 3' read-through adapter
  trimming, including indel-tolerant matching for SBX homopolymer variants
  (`AAGAGC`→`AAAGAGC`) and 3'-end partial adapters.
- **M3** — demultiplex by i5/i7 against the sample sheet into per-sample FASTQs, with
  offset-tolerant (windowed) index matching. BCL Convert v2 sample sheet parser
  (flat `sample,i7,i5` CSV fallback); normalized flat `sample_sheet.csv` emitted.
- **M4** — cell-barcode whitelist gate (`737K-arc-v1`, 2-bit packed), read-class
  tiers, and `metrics_summary.json` run accounting.
- **M5** — single-indel CB realignment behind `--realign`.
- **M6** — parallel processing across files (rayon); loss analysis with location
  confidence and `--lenient` rescue.
- CB reframe rescue for homopolymer-indel frame shifts (`--cb-reframe`, on by default).
- Single-member gzip output (Cell Ranger mis-frames records at multi-member gzip
  boundaries) written during processing.
- `--split` to shard each sample across N gzip files for parallel compression.
- End-of-run RUN SUMMARY report: timing, throughput, and lost breakdown.
- Ubuntu install script (`install.sh`) and installation manual (`INSTALL.md`).

[Unreleased]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbx10x/compare/v0.3.1...HEAD
[0.3.1]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbx10x/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbx10x/releases/tag/v0.3.0
[0.2.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbx10x/releases/tag/v0.2.0
