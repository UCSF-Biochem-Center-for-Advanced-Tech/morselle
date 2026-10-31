# Changelog

All notable changes to `morselle` are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- `extract --failed-output` no longer serialises the run. Every worker
  gzipped its failed reads into the one shared stream while holding its
  lock, so when failures were common the whole extract ran on about one
  core: on the full CAT02 3p run (53% failed), 4 of 64 threads were busy,
  at ~65k reads/s against 1.3M reads/s without the option. Each worker now
  compresses its batch into its own gzip member before taking the lock, so
  only the file write is serialised. The file is now multi-member gzip, which
  gzip/zcat, Python `gzip` and other standard readers decode as one stream.
  The failed reads are unchanged; only their order varies, as before. On one
  CAT06 shard with 14 threads: 38.9 s before, 29.1 s after, 27.0 s without
  `--failed-output`.

## [0.8.1] — 2026-10-30

### Changed
- `scripts/morselle-pipeline.sh` and `scripts/morselle-gex-pipeline.sh` no
  longer pass `extract --lenient` (Hamming-2 barcode rescue), and the MANUAL
  examples drop it. Scored against Illumina Cell Ranger molecules from the
  same library (CAT02), only about half of Hamming-2 barcodes were right
  (Hamming-1 78-83%); USAGE already said to leave it off. In 0.8.0 the
  re-lock and one-extra-base rescues run first and recover many of the same
  reads at Hamming-1 precision. Pass `--lenient` to `extract` directly to
  get the old behaviour.

## [0.8.0] — 2026-10-30

Headline, default settings against 0.7.0 on the same reads, barcodes
scored against Illumina Cell Ranger molecules from the same libraries:

| Data | Assigned reads | Reads on a right barcode | Reads on a wrong barcode |
|---|---|---|---|
| CAT02 Multiome, barcode-first 4M | 56.14% → 57.98% | +4.7% (M2_GEX) | −26% |
| CAT02 Multiome, barcode-last 2M | 45.35% → 45.69% | +2.1% (M2_GEX) | −30% |
| CAT06 Flex v2, one 43.9M shard | 68.46% → 70.84% | re-lock 91%, insertion 93% right (Hamming-1 91%) | |

Throughput is unchanged.

### Added
- `extract --flex`: Flex v2 reads have no dT run after the UMI, so almost
  no barcode was "confidently located" (both flanks found) and the rescues
  that require it (`--cb-insertion`, `--cb-reframe`, `--lenient`) never ran
  (CAT06: 0.1% of located barcodes). With `--flex`, the Flex v2 fixed
  sequence `TTGCTAGGACCG`, found with at most one edit 0-3 bp (plus 2 bp
  slack) after the UMI, now serves as the 3' flank: 93.8% confidently
  located. On a CAT06 shard with `--cb-insertion`, 818,380 reads (1.9%) are
  rescued, assigned reads go from 69.15% to 70.84%, and 93.3% of rescued
  barcodes are right against Illumina molecules (Hamming-1 91.6%). Runs
  without `--flex` are unchanged (CAT02 output byte-identical).
- `extract --index-indels`: match an i5/i7 that carries one insertion or
  deletion inside the index (precomputed one-indel lookup tables, floated
  over the same `--index-window` frames). It runs only after the
  mismatch-only search finds nothing, so reads that match today are
  unchanged. Each index is enabled only when the sample sheet is safe for it:
  minimum pairwise edit distance ≥ 3, and no one-edit variant of any index,
  placed in its real adapter flanks, is sent to another sample by the indel
  pass. On by default since 0.8.0 (see Changed). On the CAT02 Multiome slices: assigned reads +0.74%
  (barcode-first) / +0.20% (barcode-last) of all reads, dual-index
  assignments +2.2 / +1.9 points, no read moved between samples (reads whose
  two indexes now name different samples become conflicts, 0.04% / 0.03%),
  throughput unchanged.
- `extract --r1-tail-relock` (`XB:i:7`): an insertion or deletion in the
  R1-primer tail can leave the mismatch-only lock a base or two off, so the
  16 bp barcode grab fails the whitelist. For such barcode failures (after
  realign and reframe), re-find the tail's last 15 bp with up to two edits
  (bit-parallel Myers search over the first 120 bp) and read the barcode and
  UMI right after it; accepted only when exactly one tied candidate start
  gives a whitelist or unique within-1 barcode. On by default since 0.8.0 (see Changed). CAT02
  Multiome slices: 14,412 / 2,473 reads rescued, assigned reads +0.54% /
  +0.23% of all reads, rescued barcodes on cell-like barcodes at 58.2% /
  29.6% (exact matches 60.3% / 30.1%, reframe 38.4% / 19.6%), no
  previously emitted read changed, throughput unchanged.
  `whitelist.relocked` counts them. The index demux ran on the original
  lock, so its i5 anchor was off by the same base or two: 15.8% of re-locked
  reads were undetermined on the barcode-first slice (exact: 3.4%). A
  re-locked read the early demux left unassigned is now demultiplexed again
  from the re-located barcode start (`whitelist.relock_reanchored`):
  undetermined re-locked reads 15.8% → 5.3% / 4.1% → 2.6%, +1,346 / +35
  reads assigned, no previously assigned read changed sample.
- `extract --cb-insertion` (`XB:i:8`): for a confidently located barcode
  failure the earlier rescues leave, accept the unique whitelist barcode
  that removes one extra base inside the 17 bp at the barcode lock; no
  one-deletion reading may name another barcode, and the UMI is the 12 bp
  after it. Insertions only: re-measured after `--r1-tail-relock`, unique
  one-edit hits among barcode failures are 1.6% / 0.95% of reads, but
  deletion readings land on cell-like barcodes at 25.8% / 7.0% (near
  chance) against 50.0% / 20.8% for insertions. Gated on confident location
  like reframe (ungated: +1.22% / +0.80% assigned but 45.8% / 17.6%
  cell-like). On by default since 0.8.0 (see Changed). CAT02 slices: assigned +0.89% / +0.53% of all
  reads, cell-like 48.6% / 19.0% (Hamming-1 47.2% / 24.2%, reframe 37.8% /
  19.5%), no previously emitted read changed. `whitelist.cb_insertion`
  counts them.
- `tests/extract_relock_e2e.rs`: `extract --whitelist` end to end on
  generated reads with a generated whitelist (10x whitelists are not
  committed, so the extract golden runs without one). Every one-base
  insertion and deletion in the R1 primer's last 15 bp is planted; every
  emitted barcode and UMI must equal the planted one, the exact, realign,
  reframe and re-lock routes must all occur, and turning `--r1-tail-relock`
  on must not change any read emitted without it.
- `metrics_summary.json` `demux.index_evidence` and two new log lines: per
  R1-located read, `full` (both indexes seen and identified to one sample)
  versus `partial`, which index regions were seen, which were identified,
  and how many only through `--index-indels`; `demux.index_indels_i5/_i7`
  record whether it was active.
- `tests/extract_golden.rs` pins `extract`'s output (FASTQ digests and core
  counters) on the smoke fixture, as `golden_outputs.rs` does for
  `molecules`/`matrix`.

### Changed
- `extract --cb-insertion` is now on by default (`--cb-insertion false`
  turns it off; the bare flag still works). Against Illumina Cell Ranger
  molecules from the same libraries, its barcodes are right about as often
  as Hamming-1's: CAT02 Multiome 89% / 78% (Hamming-1 78% / 83%), CAT06 Flex
  93% (Hamming-1 92%, with the Flex fixed-sequence flank). It adds 0.9% /
  0.5% of reads on the CAT02 slices and 1.7 points on a CAT06 shard.
- `extract --r1-tail-relock` and `--index-indels` are now on by default
  (`--r1-tail-relock false` / `--index-indels false` turn them off; the bare
  flags still work). Against Illumina Cell Ranger molecules from the same
  libraries, re-locked barcodes are right as often as Hamming-1's (CAT02
  Multiome 91% / 79% vs 78% / 83%; CAT06 Flex 91% vs 91%). With both on and
  reframe off, assigned reads go from 56.14% / 45.35% (0.7.0) to 57.98% /
  45.69% on CAT02 slices and from 68.46% to 69.15% on a CAT06 shard, with
  about 30% fewer wrong-barcode reads on CAT02, at unchanged throughput.
  The extract golden is regenerated: the default cases now use one-indel
  index matching (smoke fixture: 2 more dual-index calls, 2 more reads
  assigned with `--demux-index both`), and a new `--index-indels false`
  case reproduces the previous default output exactly.
- `extract --cb-reframe` is now off by default. Scored against Cell Ranger
  molecules from the Illumina run of the same library (CAT02 M2_GEX: a read
  counts when an Illumina molecule on its barcode has a UMI within one edit,
  less a random-UMI control, relative to exact-match reads), only 51% / 59%
  of reframed barcodes are right (barcode-first / barcode-last), against
  Hamming-1 78% / 83%, re-lock 91% / 79%, insertion 89% / 78% and realign
  97% / 99%. Reframe had rescued about 1.6% / 1.3% of reads. Pass
  `--cb-reframe true` for the previous behaviour.
- Demultiplexing decides per read which indexes it needs: both i5 and i7
  are looked up, one alone assigns when it names exactly one sample, and
  otherwise the read goes to the single sample both indexes agree on (no
  common sample = conflict, several = undetermined). An index key several
  samples share -- a combinatorial sheet reusing an i7, or two indexes whose
  mismatch neighbourhoods overlap -- used to match nothing; it now keeps its
  candidate samples for the other index to choose among. Auto
  `--demux-index` therefore always resolves to `dual`; before, it fell back
  to a single index or `both` when an index was not standalone-safe across
  the whole sheet (a shared-i7 sheet could never assign those samples under
  `both`). Explicit `both`/`i5`/`i7` keep their meaning. Unchanged on
  CE15204 at the default `--index-max-mm 1` (no shared keys; both CAT02
  slices identical). At `--index-max-mm 2`, where auto used to pick `both`:
  assigned reads +24.5% (barcode-first) / +7.7% (barcode-last) of all
  reads, no previously assigned read changed sample, and against the
  default `--index-max-mm 1` run 23 / 20 reads (≤ 0.004%) differ.
  `demux.index_evidence.identified_i5_shared` / `_i7_shared` count shared
  matches.
- Sample-index matching uses one table per index, built from the sample
  sheet at startup, holding every allowed variant (≤`--index-max-mm`
  mismatches, plus one-indel variants with `--index-indels`). Keys are
  2-bit packed and looked up in a directly addressed array, so each frame is
  one array read instead of a byte-slice hash. Per read: 36 → 21 ns when
  both indexes match exactly, 125–395 → 39–86 ns when nothing matches; 384
  samples cost the same as 14. Output is identical (extract golden and both
  CAT02 slices, with and without `--index-indels`); end-to-end throughput is
  unchanged because demux was already a small share of `extract`.
- `benches/pipeline_hot_paths.rs`: `demux/*` benchmarks (CE15204 and a
  generated 384-sample sheet; exact, indel, and no-match reads).

### Fixed
- Rescued barcodes now take their UMI from right after the corrected
  barcode as it sits in the read, for every rescue route (Hamming-1,
  realign, reframe, re-lock, insertion, Hamming-2). The barcode end is
  where the whitelist barcode aligns with the fewest edits (route default
  ±1); when an indel at the barcode's last bases and a substitution fit
  equally, the exact poly-T start decides. Realign and reframe used to
  anchor the UMI on the poly-T start, which tolerates isolated mismatches,
  so a UMI ending in T or T-rich (`..TTATTG`) was partly absorbed and came
  out shifted, sometimes by several bases. Compared with Illumina UB on the
  same library (UMIs ending in T: 26.4%), on CAT02 slices UMIs from the
  rescue routes ended in T 17.8% / 14.9% (barcode-first / barcode-last)
  before and 28.5% / 30.3% after (exact reads 26.6% / 28.1%); their largest
  per-position deviation from exact reads fell from 9.6 / 13.7 to 1.9 / 2.2
  points. Changed UMIs: reframe 78% (old ones were mostly shifted),
  Hamming-1 6.7% / 4.6% (an indel at the barcode end read as a mismatch),
  realign, re-lock and insertion 0.5-7%. Barcodes, samples, routes and
  exact-match reads are unchanged. An insertion before the barcode's last
  base followed by a UMI ending in T still reads as a substitution (the
  read cannot tell them apart).
- The `--index-window` float no longer takes a one-mismatch hit in a
  shifted frame when another sample competes -- another frame naming a
  different sample, or the read being one insertion or deletion from
  another sample's index. Such a hit is two edits from the nominal index,
  so the index now abstains instead (`demux.index_evidence.i5_float_tie` /
  `_i7_float_tie`). In CE15204, V4_GEX and M2_GEX i5s are one-base shifts of
  each other, and a V4 i5 with one insertion used to land on M2: of the
  1,176 one-edit i5 variants, 0 now reach the wrong sample (was 1), and 2
  abstain. The one-indel tables are now always built for this check
  (`--index-indels` still decides whether they assign reads). CAT02 slices:
  66 / 6 index abstentions per 2M / 1M reads, no read changed sample, 26 /
  0 reads became undetermined, and 3 / 1 reads whose indexes used to
  conflict are now assigned on the other index. Demux bench unchanged.

## [0.7.0] — 2026-10-21

First public (alpha) release, under the name **MORSELLE** (MOlecular Read
Structure, Extraction, Labeling & Enumeration). Counting behaviour is
unchanged from 0.6.4: the golden outputs of `molecules` and `matrix` are
identical in every mode.

### Changed
- Renamed the project, crate, library, and every command from its internal
  pre-release name to `morselle`. The public command is `morselle`; its
  private engines are `morselle-extract`, `morselle-annotate`,
  `morselle-molecules`, `morselle-illumina`, `morselle-cells`,
  `morselle-probe-demux`, `morselle-guides`, `morselle-protospacer-call`,
  and `morselle-probeset-align` (plus `bam-audit` and `cr-bam-to-mtx`).
- Pipeline drivers are `scripts/morselle-pipeline.sh` and
  `scripts/morselle-gex-pipeline.sh`.
- Environment variables use the `MORSELLE_` prefix (`MORSELLE_LIBEXEC`,
  `MORSELLE_RAPIDGZIP`, ...). `install.sh` installs to
  `<prefix>/bin/morselle`, `<prefix>/libexec/morselle/` and
  `<prefix>/share/morselle/`. `extract`'s default output directory is
  `morselle_out`.
- Repository moved to
  <https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle>.

### Removed
- Presentation decks, presentation-building scripts, internal design and
  planning notes, and dataset-specific analysis scripts are no longer part
  of the repository.
- `uninstall.sh --remove-legacy` (it only removed pre-release binaries).

## [0.6.4] — 2026-10-16

### Fixed
- SBX CRISPR guide capture no longer overcounts molecules about 6.5x.
  On CAT06 L1, compared with Cell Ranger on full-depth Illumina in the cells
  both called, morselle had 6–7x more UMIs wherever a guide was really
  present. That pushed cells with 1–2 true guide UMIs over the 3-UMI
  calling floor: RZ211 CD45 81% of cells called vs 58%, Neg 4.1% vs 0.9%.
  On 3,000 RZ211 cells with depth-matched reads (10.1 M each), the cause
  was UMI errors in a saturated library (~450 reads per molecule).
  1-mismatch correction cannot merge SBX's 1-base UMI indels (6.1x Cell
  Ranger's molecules), and rare low-quality reads with several UMI errors
  each added one more. It is not an extract problem: these reads have normal
  barcode routes, anchors and UMI position, but 50% have two or more
  sub-Q20 UMI bases, vs 2.3% of real-molecule reads.
  `scripts/run_flex_multi.sh` now runs CRISPR `molecules` with the two new
  rules below. Recounting CP1 L1 CRISPR from its guide BAM with them:
  - Molecules: 6.11 M → 2.08 M; 814 k singletons folded.
  - Guide UMIs on shared cells: 1.01–1.12x Cell Ranger in every sample
    (was 1.45–6.92x); per-cell median ratio 1.00.
  - Cells called CD45 (Cell Ranger / before → after): Neg 0.5 / 4.1 → 0.8%,
    Pos 1.6 / 11.1 → 2.2%, RZ213RD 25.7 / 43.5 → 26.9%.
  - Cells called inactive: RZ212 3.3 / 17.4 → 4.8%.

  RZ211 and RZ214 still call more cells (52.8% and 53.7% vs 30.8% and
  30.1%), but not because of counts. On shared cells at Cell Ranger's own
  threshold, morselle gives 38.9% vs 36.7% and 42.0% vs 37.9%. The rest
  comes from two things, tracked in FUTURE_ROADMAP Priority 7:
  - `protospacer-call` fits a lower threshold than Cell Ranger: 4 vs 6
    and 3 vs 4.
  - Cell Ranger calls 5–8 k more cells, few of them guide-positive.

  Evidence: `CAT06/read_depth_comparison_out/diff_investigation/FINDINGS.md`.

### Added
- `molecules --umi-indel-correction`: also merge UMIs one insertion or
  deletion apart (same-length window), under the same 2n−1 read-support
  rule as 1-mismatch correction.
- `molecules --fold-singletons-min-top K`: after correction, fold each
  1-read molecule into its (cell, sample, gene) family's top molecule when
  that top has at least `K` reads. New metric `umi_singletons_folded`; the
  metrics JSON also records both settings.
- Both are off by default, so GEX output and the golden files are unchanged.
  They are meant for deeply sequenced feature-barcoding libraries; GEX
  impact has not been measured. New `tests/crispr_umi_rules_e2e.rs`
  (guide BAM fixture shaped like real `guides` output) checks that each
  rule removes only its own error class and that the matrix agrees with
  `molecules`.

## [0.6.3] — 2026-10-16

### Added
- `morselle illumina --trim-tso` now also trims the template-switch oligo
  from the start of a non-merged R2 (before, it acted only on `--overlap`
  merged fragments, so unmerged R2 was never trimmed). In 3' GEX a fragment
  that reached the transcript's 5' end starts R2 in the TSO: 8.5% of 1 M raw
  CAT02 M2_GEX R2 reads, 75% of those with all 30 bases and the rest with a
  suffix. `illumina::tso::r2_tso_prefix` finds the whole TSO or a suffix of
  at least 12 bases (one mismatch per 10 bases). Still off by default,
  and measured not worth turning on: on 2 M CAT02 M2_GEX pairs mapped with
  the 0.6.2 preset it trimmed 10.7% of emitted reads and **lowered**
  countable reads by 0.59% (739,622 -> 735,277). minimap2 already
  soft-clips the TSO (median aligned length 60 bases with or without it),
  and without it the same alignments get lower MAPQ (trimmed reads at
  MAPQ 60: 52,392 -> 38,985; 4,582 of the 5,775 lost reads stay at the same
  locus). SBX `extract --trim-tso` trimmed 4.5% (3') and 3.2% (5') of
  reads with no change in countable reads (+/- 0.01%). The `r2_tso_*`
  counters remain useful as a TSO-rate QC metric.
- `illumina` metrics gain `r2_tso_trimmed`, `r2_tso_discarded` and
  `r2_tso_bases_trimmed`; a pair whose cDNA after the TSO is below
  `--min-cdna-len` fails as `r2_tso_short` (new per-barcode column
  `failed_r2_tso_short`). `illumina::Metrics` and `illumina::PairMetrics`
  gain the fields; struct literals must set them.

## [0.6.2] — 2026-10-16

### Changed
- **Counts change:** `morselle-pipeline.sh` and `morselle-gex-pipeline.sh`
  map with `minimap2 -x splice -m 20 -s 30` instead of `-x splice:sr`.
  minimap2 2.28 accepts `-x splice:sr` but runs the long-read `splice`
  preset (2.29 added a real short-read preset; 2.30's `-x splice` gives the
  same output as 2.28's `-x splice:sr`), so every earlier run and
  validation used `splice` scoring. Naming it keeps the mapping the same on
  any minimap2 version. Under it, 90-base Illumina R2 often fell below the
  minimum chain and alignment scores (19.5% unmapped); `-m 20 -s 30` lowers
  them. CAT02 M2_GEX at matched depth (Cell Ranger Illumina 190.9 M UMIs in
  cells): Illumina 122.8 M -> 151.4 M UMIs in cells (64% -> 79% of Cell
  Ranger; per-cell r 0.997), SBX 3' 187.3 M -> 190.8 M, SBX 5' 185.9 M ->
  193.2 M (Cell Ranger on the same SBX 5' reads: 192.3 M). Reads already
  countable keep their gene, and the repeat-artifact genes do not rise.
  `-f 0.00002` (keep more repeat seeds) added little for 4.3x the CPU and
  is not used.
- **Counts change:** `morselle illumina` now trims poly-A read-through in
  unmerged R2 by default (`--trim-r2-polya`; 3' GEX and tag-seq both read
  into the tail). `--trim-r2-polya false` restores the old behaviour. The
  bare flag still works.
- **Counts change (fix):** the `--trim-r2-polya` tail must now be at least
  `--r2-polya-min-run` bases and at least 90% A
  (`homopolymer::first_tail_run`). The old rule tolerated a mismatch at
  every other base and cut at A-rich transcript sequence: on DZEC04 (5 M
  pairs, no `--overlap`) it discarded 123,568 counted reads (5.8%), 68% of
  them aligned over all 150 bases at MAPQ 60, losing 37% of ALB's reads and
  13-35% of several MT and ribosomal genes. With the fix, 8,831 counted
  reads (0.4%) are discarded, and molecules counted go from 1,989,952
  (trim off) to 2,004,230 (trim on, previously 1,809,922).
- **Counts change:** low-complexity cDNA (below) is dropped by default at
  conversion (reason `low_complexity`): a dead read 2 is a library fault,
  not biology. `--drop-low-complexity false` keeps such reads, tagged
  `YL:i:1`. `complexity::ComplexityCfg::default()` now drops, so
  `illumina::pair::convert_pair` / `convert_pair_with` do too.
- **Counts change:** `molecules` now counts only confidently mapped reads,
  as Cell Ranger does: MAPQ >= `--min-mapq` (default 60) or rescued by
  `annotate` (`mm:i:1`). On DZEC05 this removes the poly-G LINC00486
  artifact (22% of reads -> 0.01%); on a normal library (DZEC04 subset) it
  drops ~3% of single-gene reads, all low-MAPQ intronic placements.
  `--min-mapq 0` restores the old counts. See MANUAL "Confident mapping and
  multimapper rescue".
- `annotate` writes `mm:i:1` on a primary read below `--confident-mapq`
  (default 60) whose genes are one splice-compatible sense exonic gene
  (Cell Ranger's `rescue_alignments_se`), keeping minimap2's MAPQ. Intronic
  hits are never rescued. With `--multimap-tags` the rule now spans the
  query group like Cell Ranger's (no rescue if any hit is confident; only
  exonic-compatible hits count; `mm` on the primary only). The old looser
  group rule (every hit tagged when all hits' genes were one gene, intronic
  included) remains under `--confident-mapq 0`.
- **Counts change:** `annotate --max-end-overhang` (default 3) lets a read's
  first and last aligned block extend up to 3 bases past their exon and
  still be transcript-compatible (soft-clipped; `TX` position clamped;
  internal junctions stay exact). minimap2 extends rather than splices a
  1-3 base end anchor, which left such reads exonic with no gene: about 4%
  of DZEC04 records (55% of the 7% exonic-without-gene reads), 85% of them
  MAPQ 60. New `end_overhang` counter and per-barcode column. `0` restores
  exact ends. `annotate::Call` gains `end_overhang`, `RecordOptions` gains
  `max_end_overhang`; `Annotation::call_with` and
  `call::compatible_position_with` are added, and `call` /
  `compatible_position` keep exact ends.
- Before this, `molecules` never read `NH` or `mm`, so `--multimap-tags`
  had no effect on counts. The earlier CAT02 match to Cell Ranger came from
  `cr-bam-to-mtx` on Cell Ranger's own BAM, not from `molecules`.
- `annotate::Cli` gains `confident_mapq`, `annotate::Call` gains
  `exonic_compatible`, and `molecules::MoleculeOptions` gains `min_mapq`;
  struct literals must set them.

### Fixed
- `morselle-pipeline.sh --input-mode illumina` with more than one `--r1`/
  `--r2` pair (e.g. one per lane) failed after conversion: `illumina` writes
  one `.shard-NNN` output per pair, and the script looked for a single
  file. It now converts into a scratch directory, moves every shard out,
  marks the stage done, and maps all shards. Found on CAT02 M2_GEX (9 pairs).

### Added
- Low-complexity cDNA check at conversion (`illumina`, `extract`): a read
  where one base is >= 80% of the cDNA or the longest single-base run
  covers >= 50% (`--low-complexity-base-fraction`,
  `--low-complexity-run-fraction`) is dropped, or with
  `--drop-low-complexity false` emitted tagged `YL:i:1`, and counted
  (`low_complexity` and `low_complexity_by_base` in the metrics JSON,
  `low_complexity` in the per-barcode tables), with a warning above 5% of
  reads that are poly-G/N (all flagged reads for `illumina` without
  `--trim-r2-polya false`, where poly-A read-through is expected). Normal
  reads are written exactly as before. `morselle-probeset-align` copies `YL` into its BAM.
- `molecules` metrics gain an `artifacts` report: gene-assigned reads, how
  many are low-complexity (`YL:i:1`) or MAPQ < 60 without `mm:i:1`, and the
  genes where at least half of >= 100 reads are flagged; a listed gene with
  >= 1% of reads prints a warning. Nothing is dropped by it.
- `illumina::pair::convert_pair_full` (takes the low-complexity settings);
  `convert_pair` / `convert_pair_with` use the defaults. `illumina::Cli`
  and `config::Params` gain `complexity` (`complexity::ComplexityCfg`),
  `accounting::Metrics` and `molecules::Metrics` gain the new fields, and
  `MoleculeTags` gains `yl`; struct literals must set them.
- Counters: `annotate` `rescued`; `molecules` `rescued_records`,
  `low_mapq_records`, `multimapper_records` and `min_mapq`, with matching
  per-barcode columns (`rescued`; `rescued_reads`, `low_mapq_reads`,
  `multimapper_reads`).
- `annotate_record_with` / `RecordOptions`, and `collect_evidence_with`,
  `collect_region_with`, `collect_evidence_indexed_with` / `RecordFilter`.
  The existing functions keep their signatures and apply no MAPQ rule.
- Tests: SAM -> annotate -> molecules -> matrix chain tests for both routes
  (`tests/regressions/mapping_quality.rs`) and unit tests for the gate and
  `NH`/`mm` parsing.

### Notes
- Golden outputs regenerated: the fixtures give 1 in 5 reads MAPQ 3, and
  those are no longer eligible (`eligible_records` 11,491 -> 9,182 in the
  Flex mode, for example). The tag-seq ground-truth tests now count only
  MAPQ >= 60 reads as true molecules, and the indexed-representative
  regression runs with `--min-mapq 0` because it tests representative
  choice across MAPQs.

## [0.6.1] — 2026-10-08

### Added
- `--barcode-metrics` on `illumina`, `extract`, `annotate` and `molecules`
  (and on `morselle-pipeline.sh` / `morselle-gex-pipeline.sh`, which pass it
  to every stage) writes a per-barcode gzip TSV next to the stage's metrics
  JSON (`<metrics>.barcodes.tsv.gz`). Each count column sums exactly to the
  matching run-wide counter; `__no_barcode__` and `__unresolved_barcode__`
  collect reads without a usable barcode. See MANUAL "Per-barcode read
  metrics" for the columns and how to join the tables into a per-cell funnel.
- Off by default. Measured on real-data subsets, a run without the flag is
  unchanged in time and memory versus 0.6.0. With it: `annotate` +7% and SBX
  `extract` +12% wall time; `extract` also holds one row per observed
  barcode/sample pair (+375 MB on one CAT06 Flex chunk).
- `illumina::pair::convert_pair_with` (takes the flag); `convert_pair` keeps
  its signature and writes no table.
- Tests: an Illumina -> annotate -> molecules chain test reconciling every
  stage's table with its pool metrics and with the next stage
  (`tests/annotate_chain_e2e.rs`), a Flex probe-route chain test
  (`tests/per_barcode_chain_e2e.rs`), and a check that no table is written
  without the flag.

### Changed
- `illumina::Cli` and `annotate::Cli` gain a public `barcode_metrics` field;
  code constructing them with struct literals must set it.

### Notes
- Count matrices, BAMs and FASTQs are unchanged (golden outputs untouched).
- A first implementation kept one table per worker thread; on 5 CAT06 chunks
  with 64 threads `extract` peaked at 11.8 GB (0.83 GB without). `extract`
  now shares one sharded table across workers, and `molecules` keys its
  counters by the packed (cell, sample) identity it already groups on
  (a string-keyed table had cost it ~5% wall time).

## [0.6.0] — 2026-10-07

Correctness review of the 0.5.6–0.5.9 refactoring and repairs, with
regression coverage and a reproducible server validation procedure.
Server real-data acceptance (docs/SERVER_TESTING.md) against 0.5.6 passed
on 2026-10-08: CAT02 GEX (SBX, Illumina), CAT06 Flex multi GEX + CRISPR
(SBX, Illumina) and DZEC04 tag-seq gave identical count matrices, called
cells and protospacer calls in every run; molecule BAMs differ only in the
`xf` representative read. `molecules` was up to 3.1x faster with up to 40% less
peak memory (DZEC04 full depth: 27:20 -> 15:54, 45.9 -> 27.8 GB).

### Fixed
- `install.sh` installs all 11 private engines. Since 0.5.2 it installed
  only 7, omitting `morselle-probeset-align`, `morselle-probe-demux`,
  `morselle-guides` and `morselle-protospacer-call`, so an installed
  `morselle` failed every Flex and CRISPR command with "helper not found"
  (found by the server acceptance run on CAT06). The helper list is
  now one `HELPERS` array, and `tests/pipeline_scripts.rs` checks it
  against the dispatch table in `src/bin/morselle.rs`.
- Molecule output now uses the same present `PB` identity as evidence collection
  in every assay; `--flex` still requires it. Indexed representatives use unique
  physical BAM offsets, preventing duplicate representatives and selecting the
  first physical record on equal quality. Multiple input BAMs with incompatible
  reference dictionaries are rejected; sharded output declares unsorted order.
- Fractional MEX output retains full precision, including weights below 1e-5;
  non-finite/negative weights and malformed CSC matrices are rejected. Fractional
  saturation counts each collision group once. Protospacer calling preserves
  fractional values through fitting, threshold comparison, and summaries.
- Automatic EmptyDrops rejects fractional inputs because its integer multinomial
  null is not calibrated for them; explicit barcode/forced selection is supported.
- Pipeline drivers work with macOS Bash empty arrays, bind restarts to input and
  analysis configuration, and publish matrix/cell outputs only after successful
  completion. Existing output directories without a restart manifest require a
  fresh output directory. Flex worker arguments preserve whitespace.
- Restored the public `annotate_bam::{Cli, run}` compatibility path. Reporting
  consumes molecule/cell metrics and real pipeline filenames, labels historical
  slides, and is split into modules below the file ceiling. UMI benchmarks now
  use valid, deterministic identity groups.
- Added 16 regression tests without removing existing tests, plus Python restart
  and reporting cases; CI now checks Linux and macOS. Golden updates are limited
  to the intentional indexed tie-break and fractional MEX precision changes.

### Changed
- **`morselle molecules --tagseq` counts every gene in a UMI collision in
  full again (`XW:S`, the 0.5.6 behaviour); fractional splitting is opt-in
  (`--fractional-collisions`).** At tag-seq depth collisions are mostly
  distinct molecules that drew the same UMI (DZEC03: 81% of `(cell, UMI)`
  pairs, genomically unrelated genes), and fractional splitting counts each
  collision as one molecule — an undercount concentrated in the most
  expressed genes, plus non-integer matrices. `--no-fractional-collisions`
  is still accepted (no-op).

### Added
- `docs/SERVER_TESTING.md`: release installation, fresh-output real-data runs,
  indexed/sharded parity checks, tag-seq mode comparison, and acceptance criteria.
- `morselle molecules --position-aware-umis`: molecules keyed on
  `(cell, gene, UMI, strand-aware unclipped 5′ start)`, so chance UMI
  matches stay separate across and within genes; UMI correction then runs
  per start, even with `--tagseq`. No memory cost (the start is a 24-bit
  tag in spare bits of the packed UMI). On an oversubscribed simulation
  with known truth: 0.02% error vs 10% (shared) and 74% (fractional).
- `scripts/compare_tagseq_modes.sh` / `.py`: run a BAM through each mode
  and compare with a reference matrix (e.g. Illumina) on normalised
  correlation, per-gene log2 ratio, and expression-dependent bias.
- `scripts/morselle-pipeline.sh --fractional-collisions`,
  `--position-aware-umis`.
- molecules metrics JSON records `fractional_collisions` and
  `position_aware_umis`.

## [0.5.9] — 2026-10-07

Review and repair of the 0.5.2–0.5.8 changes made on 2026-10-06. Replaces
0.5.7 and 0.5.8, which are **BROKEN** (see their entries). Output of
`molecules` + `matrix` is identical to 0.5.6 in every mode 0.5.6 had
(verified by `tests/golden_outputs.rs` against the 0.5.6 binaries).

### Fixed
- `morselle molecules`: restored the 0.5.6 UMI-correction, collision and BAM
  tagging logic that the 0.5.7 modularization rewrote (crash on existing `UB`;
  missing `xf`/`XV`/`XW`; serial-mode representatives never matched; `ZG`
  overwritten with won/lost; `--tagseq` ignored; O(n²) decision step;
  uncompressed output; shard scratch in the system temp dir). Output is
  record-for-record identical to 0.5.6 in every non-fractional mode.
- Fractional collisions re-implemented on top: `XW:A:F` plus `ZW:f` weight,
  `xf:25` on the representative, `ZG` untouched. Under `--tagseq` UMI
  correction is skipped again (as in 0.5.6). A fractional collision group adds
  1 to `molecules_counted`.
- `morselle cells` (and `protospacer-call`) read `f32` fractional matrices
  without truncating (2.5 was read as 2, 0.25 as 0); filtered matrices keep
  fractional values.
- `scripts/morselle-pipeline.sh` never parsed (syntax error since 0.5.4) and
  passed flags no subcommand accepts; rewritten so the GEX path matches the
  production GEX script and Flex follows `run_full_cat06.sh`.
  `scripts/morselle-gex-pipeline.sh` is restored to the standalone 0.5.6 script.
- `probeset-align --threads`: output BAM is in input order again (was
  scheduling order); record-construction errors are no longer ignored; BAM
  compression uses the worker pool.
- Restored ~1,650 comment lines (including `--help` texts) and two test
  assertions dropped by the 0.5.7 modularization; code unchanged.
- MANUAL Flex section used non-existent `probe-demux` flags; fractional
  section described the old `ZG:Z:fractional` tagging.
- References to plan documents moved to `docs/archive/`; binary count in
  TECHNICAL_DESCRIPTION.md (11, not 12). The `codebase_gap_analysis.md`
  mentioned under 0.5.7 was never committed.

### Performance
Measured on a synthetic 2M-read annotated BAM with real tag layout (24 aux
tags per record); outputs byte-identical to before (SAM text incl. tag
order; MEX/H5 equal):
- `morselle molecules`: 16.2 s → 7.6 s wall (2.1×), 18.0 G → 7.6 G
  instructions; peak RSS −24% (−29% with `--tagseq`). Single raw aux-tag
  pass per record instead of ~6 `aux()` scans + String allocations;
  FxHash tables; tag appends without `push_aux`'s duplicate pre-scan;
  sorted-run grouping instead of `HashMap<_, Vec<_>>` indexes; collision
  outcomes stored inline.
- `morselle matrix`: 4.2 s → 2.7 s (1.6×); one aux pass per record and
  interned barcodes (~16 B instead of ~70 B per counted molecule).
- All BAM I/O: rust-htslib `libdeflate` feature (BGZF via libdeflate).
- `scripts/morselle-pipeline.sh --extract-no-gzip` (opt-in; PROFILING.md:
  7.3× extract, 2.2× extract+map, ~4× intermediate disk).

### Added
- `morselle molecules --compression-level 0-9` (default unchanged) and
  `scripts/morselle-pipeline.sh --molecule-compression-level`; level 1 was
  ~24% faster overall for a ~12% larger BAM on a 1M-read benchmark.
- Test gaps that let 0.5.7/0.5.8 pass every test while broken are closed:
  `tests/pipeline_chain_e2e.rs` (each stage on the previous stage's real
  output, incl. Flex probe-demux), `tests/golden_outputs.rs` (exact
  output digests for 9 modes), `tests/docs_commands.rs` (documented flags
  exist), `scripts/check_test_inventory.sh` (no silent test deletion),
  `scripts/check.sh` (single gate), GitHub Actions CI, a pre-push hook
  (`git config core.hooksPath scripts/hooks`), and testing rules in
  AGENTS.md/CLAUDE.md. All nine new chain/golden cases fail on 0.5.8.
- `tests/molecules_matrix_e2e.rs`: runs the real `morselle-molecules` and
  `cr-bam-to-mtx` binaries on a fixture with `UB`/`UR`/`ZG` and multi-gene
  `GX`, in serial, indexed and sharded modes and every collision policy.
- `tests/pipeline_scripts.rs`: `bash -n` on every script, shellcheck when
  installed, and both pipeline scripts run end to end in all input × analysis
  modes against stubs that reject any flag the real binary does not accept.

## [0.5.8] — 2026-10-06

> **BROKEN — do not use; use 0.5.9.** `morselle molecules` aborts on any
> BAM carrying `UB` (every extracted library) and otherwise writes no `xf:25`
> representatives, so GEX/tagseq matrices are empty; `morselle cells`
> truncates fractional counts; `scripts/morselle-pipeline.sh` (and the
> `morselle-gex-pipeline.sh` forwarder in 0.5.7+) does not parse.

### Documentation
- Updated `MANUAL.md` Section 7 with detailed technical documentation for fractional UMI collision resolution (`--fractional-collisions` in `morselle molecules --tagseq`) and real Matrix Market MTX / `f32` HDF5 count matrix materialization in `cr-bam-to-mtx`.
- Updated `TECHNICAL_DESCRIPTION.md` IP metadata specification to v0.5.8 (~19,400 lines of Rust across 25 library modules and 12 binary engines).
- Archived completed historic plan markdown files (`FLEX_SUPPORT_PLAN.md`, `GENE_COLLISION_RESOLUTION_PLAN.md`, `OVERLAP_PLAN.md`, `PRESENTATION_PLAN.md`, `mapping.md`) into `docs/archive/`.

## [0.5.7] — 2026-10-06

> **BROKEN — do not use; use 0.5.9.** `morselle molecules` aborts on any
> BAM carrying `UB` (every extracted library) and otherwise writes no `xf:25`
> representatives, so GEX/tagseq matrices are empty; `morselle cells`
> truncates fractional counts; `scripts/morselle-pipeline.sh` (and the
> `morselle-gex-pipeline.sh` forwarder in 0.5.7+) does not parse.

### Added
- Expanded hot-path microbenchmarks in `benches/pipeline_hot_paths.rs` (`bench_molecules_dedup`) benchmarking UMI directional graph collapse (~23 µs per 50 UMIs) and fractional vs winner-take-all decision matrix building (~128 µs per 1,000 UMIs).
- Comprehensive codebase gap analysis report (`codebase_gap_analysis.md`) evaluating input format and output table equivalence across SBX and Illumina input modes and GEX, TagSeq, and Flex analysis modes.

### Refactored
- Consolidated stub module `src/annotate_bam.rs` into `src/annotate/mod.rs` and updated `src/bin/morselle-annotate.rs`.
- Converted `scripts/morselle-gex-pipeline.sh` into a clean forwarder script invoking the unified multi-mode driver script `scripts/morselle-pipeline.sh`.

## [0.5.6] — 2026-10-06

> 0.5.4–0.5.6: the new `scripts/morselle-pipeline.sh` does not parse; the
> Rust tools are fine. Fixed in 0.5.9.

### Added
- Automated presentation & insert-size distribution reporting (`scripts/build_presentation.py`). Ingests pipeline metrics JSON files (`illumina_metrics.json`, `molecules_metrics.json`, `cell_calling_metrics.json`), exports cDNA insert-size distribution TSVs via `--export-insert-sizes`, and builds 10-slide PowerPoint presentations (`.pptx`).
- Fractional / proportional UMI collision resolution (`--fractional-collisions`) in `morselle-molecules` (`molecules` stage). Designed for high-depth bulk/pooled samples (e.g. DZEC03, 2–25M reads/sample) where 12bp UMI space is oversubscribed (~2.1x, ~81% collision rate). Instead of deterministic winner-take-all discards, fractional resolution distributes proportional molecule weights $w_i = r_i / \sum r$ across competing genes based on read support, tagging records with `ZG:Z:fractional` and emitting `ZW:f:<weight>` tags. `cr-bam-to-mtx` parses `ZW:f` weight tags and materializes `real` Matrix Market MTX files (`%%MatrixMarket matrix coordinate real general`) and `f32` CSC HDF5 matrices (`raw_feature_bc_matrix.h5`).
- Default fractional collision resolution in `--tagseq` mode (`morselle-molecules --tagseq`), with explicit control via `--fractional-collisions` and `--no-fractional-collisions`.
- Consolidated master tracking roadmap `FUTURE_ROADMAP.md` organizing completed v0.3.3–v0.5.5 milestones and open future priorities (sharded RAM bounding for `molecules`, fractional UMI assignment, automated slide deck reporting).
- Updated status headers across `FLEX_SUPPORT_PLAN.md`, `GENE_COLLISION_RESOLUTION_PLAN.md`, `OVERLAP_PLAN.md`, and `PIPELINE_NOTES.md` linking to `FUTURE_ROADMAP.md`.

## [0.5.5] — 2026-10-06

### Added
- Expanded microbenchmark suite (`benches/pipeline_hot_paths.rs`) to cover `ProbeMatcher::locate` probe sequence alignment performance (740 ns/read throughput).

## [0.5.4] — 2026-10-06

### Added
- Unified multi-mode pipeline script `scripts/morselle-pipeline.sh` supporting all combinations of `--input-mode {sbx|illumina}` and `--analysis-mode {gex|flex|tagseq}` with stage restartability and thread management.

## [0.5.3] — 2026-10-06

### Added
- Native multi-threading (`--threads N`) and chunk/batch processing via Rayon and crossbeam channels in `morselle probeset-align` (`src/probe_set.rs`), enabling multi-core scaling without requiring shell FASTQ splitting.

## [0.5.2] — 2026-10-06

### Added
- Comprehensive integration test suite for Illumina + 10x Flex and Illumina + CAT Tagseq modes (`tests/illumina_modes_integration.rs`).
- Detailed documentation for 10x Flex (`morselle probe-demux`, `morselle probeset-align`) and CAT Tagseq (`morselle molecules --tagseq`, `morselle cells --tagseq-barcodes`) in `MANUAL.md` and `README.md`.

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
  `morselle-illumina.rs` was updated for it at the time (`4a7d760`). All
  affected examples now set it to `65` (matching Flex's own default;
  `flex_r2_boundary` is `false` in every one of these examples, so the value
  is never consulted).

## [0.4.1] — 2026-09-21

### Added
- `morselle molecules --flex` explicitly enables Flex pooled-sample handling
  and requires a non-empty `PB` probe/sample barcode on counted records.

### Changed
- `morselle molecules --tagseq` no longer requires `PB`. Pooled tag-seq
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
- `morselle extract --no-gzip` writes plain, uncompressed FASTQ instead of
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
- `morselle-molecules --tagseq` opts a run out of two single-cell-shaped
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
- `morselle cells --tagseq-barcodes <FILE>` calls exactly the listed
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
- `morselle illumina --trim-r2-polya` (default off) catches poly-A
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
- `morselle::yc` encodes and decodes the Roche SBX `YC` duplex-consensus tag,
  written from the published XOOS specification (not derived from Roche's
  reference implementation, whose licence restricts derivatives to Roche
  instrument data). Covers parse/serialise, reverse-complement, trim, and
  per-position provenance, with the four published worked examples as fixtures.
- `morselle::overlap` merges overlapping Illumina mates into the full fragment — a
  port of HTStream's `hts_Overlapper` adapted to 10x geometry. Merges the *raw*
  mates before any trimming, which is what makes the alignment gapless, and
  removes adapter read-through by geometry rather than sequence matching.
- `morselle illumina --overlap` (default off) reconstructs the fragment before
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
- `MORSELLE_RAPIDGZIP=1` enables parallel `rapidgzip` input decoding for
  `morselle illumina`; `MORSELLE_RAPIDGZIP_THREADS` and
  `MORSELLE_RAPIDGZIP_POOL` tune decoder workers and pooled input slots.
- `morselle illumina --split N` uses a total logical gzip-writer budget across
  all main and failed-output classes, producing disjoint numbered shards for
  parallel compression.
- `morselle molecules --saturation <path>` emits an order-independent,
  one-evidence-pass gzip TSV saturation curve at configurable eligible-read
  intervals. It includes run-wide and per-cell detected genes and resolved,
  gene-aware molecules.
- `morselle molecules --evidence-workers N` uses a BAM index to collect
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
- `morselle annotate --multimap-tags` processes query-grouped minimap2 primary
  and secondary alignments in one stream, adding reported-hit `NH`, hit-index
  `HI`, and an `mm:i:1` single-compatible-gene multimap-rescue signal while
  retaining minimap2's native MAPQ.
- `scripts/morselle-gex-pipeline.sh --multimap-tags` retains up to ten minimap2
  secondary hits (`-N 10`), yielding top-coded `NH=11` values that match the
  audit's `NH=11+` ambiguity bucket without an intermediate name-sorted BAM.

## [0.3.4] — 2026-07-15

### Added
- `morselle annotate --unordered` performs feature assignment across parallel
  record workers for pipelines that immediately coordinate-sort the BAM.
- `scripts/morselle-gex-pipeline.sh` provides a restart-safe raw-SBX-to-counts
  workflow with independent mapping, annotation, and sorting thread controls.
- `morselle audit --cells` accounts for representative molecules retained in a
  final called-cell barcode set.
- Audit JSON now includes `NH`/MAPQ/soft-clip/score-gap ambiguity diagnostics,
  aligned-length bins, top ambiguous genes, a raw-UMI-to-representative molecule
  funnel, per-molecule read-support distribution, and a read-depth saturation
  curve. The audit wrapper and local report accept and display these fields.


## [0.3.3] — 2026-07-14

### Added
- `morselle cells` implements 10x-style GEX cell calling with bootstrapped
  OrdMag, Simple Good-Turing ambient profiling, multinomial EmptyDrops Monte Carlo
  testing, BH FDR correction, and filtered MEX/HDF5 outputs.
- `morselle illumina` converts multiple R1/R2 lanes concurrently into numbered
  minimap2 FASTQ shards instead of serializing every lane through one
  gzip stream.

### Fixed
- `morselle annotate --output -` sends diagnostics to stderr so logging cannot
  corrupt the streamed BAM consumed by `samtools`.
- `morselle illumina` removes Illumina mate/control/index comments and `/1` or
  `/2` suffixes before appending minimap2-compatible SAM tags to read names.

## [0.3.2] — 2026-07-13

### Added
- `morselle illumina` converts paired 10x Illumina R1/R2 FASTQs into the same
  whitelist-corrected, minimap2-tagged R2-only handoff produced by SBX extraction.
- `morselle audit` now reports explicit mapped/unique-gene counts, 10x reference
  `xf` values and bits, UMI correction impact, SBX collision outcome records,
  detected genes/cells, and molecules-per-cell distributions. The comparison
  report adds platform-only molecules and per-cell/per-gene concordance.

## [0.3.1] — 2026-07-13

### Fixed
- Unique Hamming-1 cell barcodes are corrected to the canonical whitelist
  sequence before minimap2/10x reference `CB` emission; ambiguous Hamming-1
  matches are rejected rather than assigned arbitrarily.

### Added
- `morselle matrix --cell-barcodes` reads a plain/gzipped 10x reference whitelist,
  applies a configurable gem-group suffix, and retains every zero-count cell column.
- `morselle annotate --genes-mtx` creates or verifies a canonical filtered gene
  manifest; `morselle matrix --features` now enforces it strictly for stable rows.
- Formal `resources/10x_barcodes/` staging and installation of the full
  available 10x barcode-resource collection under `share/morselle/barcodes/`.
- Idempotent `uninstall.sh` with custom-prefix, dry-run, data-preservation, and
  optional minimap2-companion removal controls.
- Unified `morselle` CLI with `extract`, `minimap`, `annotate`, `molecules`,
  `audit`, and `matrix` subcommands. `minimap` transparently wraps minimap2;
  processing engines are installed privately under `libexec/morselle`.
- `morselle-molecules`: two-pass, coordinate-order-preserving gene-aware UMI
  correction, read-support collision resolution, corrected `UB`, 10x reference
  `xf=17/25` eligibility/representative tagging, and `XV`/`XW` collision audit tags.
- `--processing-tags` emits `XO`, `XB`, `XI`, `XA`, and `XL` per-read SBX
  provenance in R2 FASTQ comments for minimap2 `-y` propagation. The
  `--minimap2` preset enables it automatically.
- `morselle-annotate`: strand-aware GTF feature assignment that preserves existing
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

[Unreleased]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/compare/v0.8.1...HEAD
[0.8.1]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/compare/v0.8.0...v0.8.1
[0.8.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/compare/v0.6.4...v0.7.0
[0.6.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/compare/v0.5.9...v0.6.0
[0.3.3]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/compare/v0.3.2...v0.3.3
[0.3.2]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/compare/v0.3.1...v0.3.2
[0.3.1]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/releases/tag/v0.3.0
[0.2.0]: https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle/releases/tag/v0.2.0
