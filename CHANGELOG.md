# Changelog

All notable changes to `sbxpress` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.7] — 2026-09-03

### Added
- `sbx10x::yc` encodes and decodes the Roche SBX `YC` duplex-consensus tag,
  written from the published XOOS specification (not derived from Roche's
  reference implementation, whose licence restricts derivatives to Roche
  instrument data). Covers parse/serialise, reverse-complement, trim, and
  per-position provenance, with the four published worked examples as fixtures.
- `sbx10x::overlap` merges overlapping Illumina mates into the full fragment — a
  port of HTStream's `hts_Overlapper` adapted to 10x geometry. Merges the *raw*
  mates before any trimming, which is what makes the alignment gapless, and
  removes adapter read-through by geometry rather than sequence matching.
- `sbxpress illumina --overlap` (default off) reconstructs the fragment before
  extraction. Merged pairs (`sins`/`mins`) carry the same canonical layout as an
  SBX read and go through the same `extract`, gaining a `YC:Z:` tag; unmerged
  pairs (`lins`) fall back to the existing R1-barcode / R2-cDNA path unchanged.
  Insert-class counts and a cDNA-length histogram are logged and recorded in the
  metrics JSON.
- `--failed-output` (both arms) emits every read that fails to produce a pair,
  as the original untouched record with an `XF:Z:<reason>` tag appended, so no
  input read is silently dropped. The SBX arm writes one gzipped FASTQ; the
  Illumina arm writes `<prefix>_R1.fastq.gz` / `<prefix>_R2.fastq.gz` so the
  failed pair stays re-runnable. With it set, `total_reads` equals emitted plus
  failed exactly, and the per-reason breakdown is recorded in the metrics JSON.
- `--trim-tso` (default off, both arms) trims the template-switch oligo where it
  sits between the cDNA and the 3' read-through adapter. Anchored rather than
  scanned, so a chance match inside the cDNA cannot truncate real sequence.
- `SBXPRESS_RAPIDGZIP=1` enables parallel `rapidgzip` input decoding for
  `sbxpress illumina`; `SBXPRESS_RAPIDGZIP_THREADS` and
  `SBXPRESS_RAPIDGZIP_POOL` tune decoder workers and pooled input slots.
- `sbxpress illumina --split N` uses a total logical gzip-writer budget across
  all main and failed-output classes, producing disjoint numbered shards for
  parallel compression.
- `sbxpress molecules --saturation <path>` emits an order-independent,
  one-evidence-pass gzip TSV saturation curve at configurable eligible-read
  intervals. It includes run-wide and per-cell detected genes and resolved,
  gene-aware molecules.
- `sbxpress molecules --evidence-workers N` uses a BAM index to collect
  disjoint genomic regions concurrently, then merges evidence before global UMI
  correction and gene-collision resolution. Independent `(CB, gene)` UMI groups
  and regional-map reduction are parallelized as well.
- Molecule-stage stderr checkpoints now report evidence, merge, UMI correction,
  collision resolution, saturation, and BAM-writing progress without requiring
  `RUST_LOG` configuration.

### Changed
- Illumina main and failed FASTQ gzip writers use compression level 4 for higher
  throughput.

## [0.3.5] — 2026-07-16

### Added
- `sbxpress annotate --multimap-tags` processes query-grouped minimap2 primary
  and secondary alignments in one stream, adding reported-hit `NH`, hit-index
  `HI`, and an `mm:i:1` single-compatible-gene multimap-rescue signal while
  retaining minimap2's native MAPQ.
- `scripts/sbxpress-gex-pipeline.sh --multimap-tags` retains up to ten minimap2
  secondary hits (`-N 10`), yielding top-coded `NH=11` values that match the
  audit's `NH=11+` ambiguity bucket without an intermediate name-sorted BAM.

## [0.3.4] — 2026-07-15

### Added
- `sbxpress annotate --unordered` performs feature assignment across parallel
  record workers for pipelines that immediately coordinate-sort the BAM.
- `scripts/sbxpress-gex-pipeline.sh` provides a restart-safe raw-SBX-to-counts
  workflow with independent mapping, annotation, and sorting thread controls.
- `sbxpress audit --cells` accounts for representative molecules retained in a
  final called-cell barcode set.
- Audit JSON now includes `NH`/MAPQ/soft-clip/score-gap ambiguity diagnostics,
  aligned-length bins, top ambiguous genes, a raw-UMI-to-representative molecule
  funnel, per-molecule read-support distribution, and a read-depth saturation
  curve. The audit wrapper and local report accept and display these fields.


## [0.3.3] — 2026-07-14

### Added
- `sbxpress cells` implements 10x-style GEX cell calling with bootstrapped
  OrdMag, Simple Good-Turing ambient profiling, multinomial EmptyDrops Monte Carlo
  testing, BH FDR correction, and filtered MEX/HDF5 outputs.
- `sbxpress illumina` converts multiple R1/R2 lanes concurrently into numbered
  minimap2 FASTQ shards instead of serializing every lane through one
  gzip stream.

### Fixed
- `sbxpress annotate --output -` sends diagnostics to stderr so logging cannot
  corrupt the streamed BAM consumed by `samtools`.
- `sbxpress illumina` removes Illumina mate/control/index comments and `/1` or
  `/2` suffixes before appending minimap2-compatible SAM tags to read names.

## [0.3.2] — 2026-07-13

### Added
- `sbxpress illumina` converts paired 10x Illumina R1/R2 FASTQs into the same
  whitelist-corrected, minimap2-tagged R2-only handoff produced by SBX extraction.
- `sbxpress audit` now reports explicit mapped/unique-gene counts, 10x reference
  `xf` values and bits, UMI correction impact, SBX collision outcome records,
  detected genes/cells, and molecules-per-cell distributions. The comparison
  report adds platform-only molecules and per-cell/per-gene concordance.

## [0.3.1] — 2026-07-13

### Fixed
- Unique Hamming-1 cell barcodes are corrected to the canonical whitelist
  sequence before minimap2/10x reference `CB` emission; ambiguous Hamming-1
  matches are rejected rather than assigned arbitrarily.

### Added
- `sbxpress matrix --cell-barcodes` reads a plain/gzipped 10x reference whitelist,
  applies a configurable gem-group suffix, and retains every zero-count cell column.
- `sbxpress annotate --genes-mtx` creates or verifies a canonical filtered gene
  manifest; `sbxpress matrix --features` now enforces it strictly for stable rows.
- Formal `resources/10x_barcodes/` staging and installation of the full
  available 10x barcode-resource collection under `share/sbxpress/barcodes/`.
- Idempotent `uninstall.sh` with custom-prefix, dry-run, data-preservation, and
  optional minimap2-companion removal controls.
- Unified `sbxpress` CLI with `extract`, `minimap`, `annotate`, `molecules`,
  `audit`, and `matrix` subcommands. `minimap` transparently wraps minimap2;
  processing engines are installed privately under `libexec/sbxpress`.
- `sbx10x-molecules`: two-pass, coordinate-order-preserving gene-aware UMI
  correction, read-support collision resolution, corrected `UB`, 10x reference
  `xf=17/25` eligibility/representative tagging, and `XV`/`XW` collision audit tags.
- `--processing-tags` emits `XO`, `XB`, `XI`, `XA`, and `XL` per-read SBX
  provenance in R2 FASTQ comments for minimap2 `-y` propagation. The
  `--minimap2` preset enables it automatically.
- `sbx10x-annotate`: strand-aware GTF feature assignment that preserves existing
  barcode/UMI tags and adds 10x-compatible `GX`, `GN`, `TX`, `AN`, and
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
- `cr-bam-to-mtx`: reconstructs 10x reference raw GEX matrices from corrected BAM
  molecule tags, including 10x-compatible MEX coordinate ordering, the
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
- `--cellranger` compatibility preset that applies the settings 10x reference
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
  reported by 10x reference.
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
- 10x reference `quality string length is not equal to sequence length` on long reads:
  10x reference hard-caps reads at 649 bp (fails at 650). Use `--cellranger` (or
  `--max-r2-len 649`) to cap R2 within the limit.

## [0.1.0] — milestone build (M0–M7)

Initial end-to-end pipeline converting raw Roche SBX long reads into per-sample 10x
R1/R2 FASTQs for 10x reference ARC.

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
- Single-member gzip output (10x reference mis-frames records at multi-member gzip
  boundaries) written during processing.
- `--split` to shard each sample across N gzip files for parallel compression.
- End-of-run RUN SUMMARY report: timing, throughput, and lost breakdown.
- Ubuntu install script (`install.sh`) and installation manual (`INSTALL.md`).

[Unreleased]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbxpress/compare/v0.3.3...HEAD
[0.3.3]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbxpress/compare/v0.3.2...v0.3.3
[0.3.2]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbxpress/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbxpress/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbxpress/releases/tag/v0.3.0
[0.2.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbxpress/releases/tag/v0.2.0
