# Changelog

All notable changes to `sbxpress` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.2] — 2026-10-06

### Added
- Comprehensive integration test suite for Illumina + 10x Flex and Illumina + CAT Tagseq modes (`tests/illumina_modes_integration.rs`).
- Detailed documentation for 10x Flex (`sbxpress probe-demux`, `sbxpress probeset-align`) and CAT Tagseq (`sbxpress molecules --tagseq`, `sbxpress cells --tagseq-barcodes`) in `MANUAL.md` and `README.md`.

### Changed
- Exported `resolve_tagseq_barcodes` as `pub fn` in `src/cell_calling.rs` for integration testing.

## [0.5.1] — 2026-10-05

### Fixed
- `cargo bench` (`benches/pipeline_hot_paths.rs`) failed to compile:
  `ExtractCfg::flex_min_cdna_len_no_anchor` was added in `e36ee39`; 0.5.0's
  fix covered all 16 affected `examples/*.rs` but missed this bench file,
  since `cargo test` doesn't build benches by default and the gap went
  unnoticed. Now set to `65`, matching every other fixed call site.

## [0.5.0] — 2026-10-05

### Added
- Indel-tolerant anchor matching for CRISPR/Feature Barcoding guide capture
  (`feature_capture.rs`), mirroring the single-deletion/single-insertion
  recovery already used by GEX probe matching (`probe_set.rs`): at each read
  offset, try the exact-length anchor window, then a (len-1)-shifted window
  against a zero-tolerance deletion-variant index, then a (len+1)-shifted
  window against a zero-tolerance insertion-variant index. Real validation
  on a CAT06 L1_CRI fastq chunk: 73.4% → 81.3% hit rate.
- Windowed, indel-tolerant i5/i7 index matching in `demux.rs` (`find_windowed`),
  replacing fixed-position-only index lookup.
- GMM/EM-based 2-component protospacer UMI-count threshold fitting
  (`protospacer_calling.rs`): fits a 2-component 1D Gaussian mixture (several
  k-means/percentile initializations, keeping the best by log-likelihood)
  over log1p(UMI count) and picks the log-likelihood-ratio decision boundary
  between components, replacing a simpler fixed/heuristic threshold.

### Changed
- `probe-demux` composite-barcode handling refined (`probe_demux.rs`):
  real high-occupancy-GEM filtering and tightened `(GEM, ProbeBarcode)`
  assignment logic.

### Fixed
- 16 example binaries (`examples/*.rs`) failed to compile:
  `ExtractCfg::flex_min_cdna_len_no_anchor` was added in `e36ee39` but only
  `sbxpress-illumina.rs` was updated for it at the time (`4a7d760`). All
  affected examples now set it to `65` (matching Flex's own default;
  `flex_r2_boundary` is `false` in every one of these examples, so the value
  is never consulted).

## [0.4.1] — 2026-09-21

### Added
- `sbxpress molecules --flex` explicitly enables Flex pooled-sample handling
  and requires a non-empty `PB` probe/sample barcode on counted records.

### Changed
- `sbxpress molecules --tagseq` no longer requires `PB`. Pooled tag-seq
  chemistry does not produce Flex probe-barcode tags, so tag-seq molecule
  counting now uses the cell barcode without imposing a Flex-only gate.
- `--flex` and `--tagseq` are mutually exclusive to make the molecule
  identity policy explicit.

## [0.4.0] — 2026-09-18

### Added
- `probe-demux` writes `<sample>.probe_barcodes.txt` (that sample's declared
  probe-barcode ids, straight from the sample sheet) alongside
  `<sample>.barcodes.txt`.
- `matrix --probe-barcodes <file>`: an optional per-sample probe-barcode
  allow-list. Needed by the composite-cell-barcode change below — once a GEM
  can legitimately appear in more than one sample's `--cell-barcodes`, each
  sample's `matrix` run must also gate on `PB` matching its own declared
  set, or a shared GEM's other-sample molecules leak into both matrices.
  New `probe_barcodes_skipped` metric records how many counted molecules
  the allow-list excluded.
- `sbxpress extract --no-gzip` writes plain, uncompressed FASTQ instead of
  gzip. Profiling `extract` (see `PROFILING.md`) found gzip output
  compression consuming 54% of its CPU time — worse, since every worker
  thread shares one mutex-protected encoder per output target, it's a
  serialization point, not just added work. Measured on a 2M-read run:
  19.1s → 2.4s (~8x), far more than the 54% CPU share alone would predict,
  confirming the lock-contention effect. Meant for the common case where
  `extract`'s output is about to be read once by a co-located process
  (typically a mapper) rather than stored — compression is pure overhead
  there. Verified end-to-end: `extract --no-gzip` writing into a named pipe,
  read directly by `minimap2` with no intermediate file on either end, all
  1,649,018 records and their `CR`/`CB`/`UR`/`UB` tags surviving intact.
- `matrix --cell-barcodes`/`--barcodes` accepts an unlisted counted barcode
  without erroring via `--skip-unlisted-barcodes` (default off, existing
  strict behavior unchanged) — for splitting one pooled multi-sample BAM
  per sample, where every per-sample `matrix` call legitimately sees other
  samples' barcodes too.
- `matrix --cell-barcodes` now accepts a barcode list whose entries already
  carry a gem-group suffix (e.g. copied from a real `barcodes.tsv.gz` or
  `filtered_barcodes.csv`, not just a bare 10x whitelist), matching
  `cells --tagseq-barcodes`'s existing tolerance — `--gem-group` is always
  what's actually applied, so an existing suffix is optional, never
  required, and never silently wins over the flag.

### Changed
- **Breaking**: `extract --sample-sheet` renamed to `extract
  --bclconvert-samplesheet`. The old name didn't make clear that this flag
  is scoped specifically to reads that passed through an Illumina BCL
  Convert library-conversion step and therefore carry i5/i7 indices at all
  — a future direct-Roche-library input (no Illumina conversion step) has
  no i5/i7 present, needs its own demux mechanism entirely, and must never
  be pointed at this flag. Demux stays fully opt-in either way: unset, no
  index-based splitting is attempted regardless of what's in the raw reads.
- **Breaking**: `probe-demux`'s cell barcode is now the composite (GEM
  barcode, Probe Barcode) pair, matching 10x's own "Cell Ranger's Flex
  Algorithm" support doc, rather than treating the GEM barcode alone as the
  cell. A GEM that independently clears `--min-molecules` for more than one
  sample's declared probe-barcode set now produces one `Assigned` row per
  clearing sample (each carrying that sample's own molecule count and
  purity), instead of the previous all-or-nothing reject. Real CAT06 L1
  data showed that reject firing on 60.5% of barcodes with enough signal to
  assign — far too high a rate to be real cross-sample doublet biology —
  discarding the majority of otherwise-real, correctly multiplexed cells.
  `Status::HighOccupancy` is renamed `Status::HighOccupancyGem` and now
  means something narrower and rarer: a GEM's occupancy (number of
  independently-clearing samples) exceeds a negative-binomial-fit (`r` from
  the observed mean, `p=0.5`) 99.9th-percentile cutoff on the run's whole
  occupancy distribution — mirroring 10x's own "High occupancy GEM
  filtering" step, which targets chip-loading artifacts (partial clogs,
  wetting failures), not ordinary multi-sample pooling. Measured end-to-end
  impact on a real 25-shard L1 GEX run: demux assigned-barcode count
  26,667 → 424,711 (~16x), final called cells 18,825 → 332,192 (~17.6x),
  closing the gap to Illumina's real L1 total (427,331) from 22.7x under to
  1.29x under.
- `sample_assignment.csv` gains an `occupancy` column, and `molecules`/
  `purity` are now per-assigned-sample (that sample's own sum/fraction of
  the GEM total) rather than the whole GEM's total — see above.
- `probe_demux_metrics.json`'s `high_occupancy` field is renamed
  `high_occupancy_gems` and now counts only the rare NB-tail-filtered GEMs;
  new `occupancy_r`/`occupancy_cutoff` fields record the fitted parameters.

## [0.3.8] — 2026-09-14

### Added
- `sbx10x-molecules --tagseq` opts a run out of two single-cell-shaped
  assumptions that stop holding at pooled-sample depth (tens of millions of
  reads per barcode rather than a single cell's tens of thousands):
  - **Cross-gene UMI collision adjudication is skipped.** At this depth the
    same `(cell, UMI)` value is routinely drawn independently by two or more
    genuinely distinct, unrelated genes — a 12bp-UMI-space birthday-paradox
    effect, not a mapping artifact — and the default policy of keeping only
    the higher-support gene silently discards real molecules from the rest.
    Every gene sharing a collided UMI is counted independently and tagged
    `Outcome::Shared` (`XW:A:S`) instead of `Won`/`Lost`/`Tie`;
    `collision_groups`/`collision_degree` diagnostics are still recorded so
    collision pressure remains visible.
  - **Directional UMI error-correction is skipped.** The same UMI-space
    exhaustion that causes cross-gene collisions also makes the
    1-mismatch/2x-read-count "sequencing error" heuristic unsafe for
    high-diversity (cell, gene) groups: two genuinely distinct real
    molecules can land a single Hamming distance apart by chance and get
    incorrectly merged. Every raw, exact-duplicate-collapsed node is kept as
    its own molecule. Trade-off: genuine single-base UMI sequencing errors
    are no longer corrected either, so `--tagseq` molecule counts trend
    slightly high rather than the default policy's slightly low.
  - `Metrics.tagseq` records which policy a run used. Default (single-cell)
    runs are unaffected by either change — `Outcome::Shared` is never
    produced and UMI correction still runs unless the flag is passed.
- `sbxpress cells --tagseq-barcodes <FILE>` calls exactly the listed
  barcodes as cells/samples, skipping OrdMag and EmptyDrops entirely.
  Both algorithms are shaped around droplet single-cell data (a few real
  cells against a large empty-droplet background with a shared ambient-RNA
  noise profile); a pooled tag-seq run has a different noise mechanism
  (index hopping into a small closed barcode whitelist) and the real sample
  barcodes are already known from the experiment design, so there is
  nothing to infer statistically. One barcode per line, matching the
  extraction whitelist format (a gem-group `-1` suffix, if present, is
  ignored); an unrecognized barcode is a hard error rather than a silent
  skip.
- `sbxpress illumina --trim-r2-polya` (default off) catches poly-A
  read-through in R2 on the non-merged path — either `--overlap` is off
  entirely, or this particular pair failed to merge. On a fragment shorter
  than R2's read length, R2 reads past the cDNA insert into the poly-A tail
  and, past that, into the reverse-complemented UMI/CB/primer; `--overlap`
  already removes this by reconstructing the fragment from both mates, but a
  short-R1 run (28 bp CB+UMI, no cDNA read-through possible on R1 at all)
  never has anything to merge against, so the same read-through went
  undetected. Finds the leftmost poly-A run of at least
  `--r2-polya-min-run` bp (default 12, same isolated-mismatch tolerance as
  the existing dT-anchor scan) via `homopolymer::first_run_at_least`, and
  trims R2 back to just before it; the pair is discarded (`r2_polya_short`)
  instead of emitted if what remains is below `--min-cdna-len`.
  `Metrics.r2_polya_trimmed` / `r2_polya_discarded` / `r2_polya_bases_trimmed`
  record the effect. On a 200k-pair DMTS60 sample (R1 28 bp / R2 90 bp, no
  `--overlap`), 12.5% of barcode-valid reads carried a poly-A read-through
  that `--overlap` alone could not have caught.

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
