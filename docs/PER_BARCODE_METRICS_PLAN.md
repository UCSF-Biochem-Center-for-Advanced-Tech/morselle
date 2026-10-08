# sbxpress 0.6.1 — per-barcode metrics and read-quality filtering (implementation brief)

Two tasks. Do **Part 1 first**: Part 2 reports its filter counts in Part 1's
per-barcode tables, and Part 1's reconciliation tests guard Part 2's changes.

# Part 1: per-barcode metrics in every stage

Hand this file to a Claude Code session started in
`/mnt/analysis/software_shared/sbxpress` on branch `fix/gemini-review`
(at or after `c52d940`, which is tagged `v0.6.0`). Do the work on a new branch
`feat/per-barcode-metrics`. Read `CLAUDE.md` first: its testing and change
rules apply to everything below.

## Why

Every stage reports read accounting for the whole run only. Per barcode,
sbxpress writes molecules and genes (the matrix) and `cells/cell_calls.csv`,
and nothing about reads. Answering "how many reads did barcode X get, and
where were they lost" means re-streaming every BAM record afterwards. On
DZEC05 (4 libraries, ~620M records) that took 30+ minutes of single-threaded
awk. Each stage already holds every read's barcode, so it can count per
barcode as it goes, at negligible cost.

## Goal

When a stage writes its metrics JSON, it also writes a per-barcode table next
to it. Joining the tables gives each barcode's full funnel:
read pairs → barcode valid → cDNA emitted → mapped → gene-assigned →
eligible → molecules.

Out of scope: poly-G detection (a one-off for DZEC05), new filtering, any
change to BAM, FASTQ or matrix output. Molecule counts, matrices and golden
outputs must not change.

## Output format (all stages)

- Path: `<metrics path without .json>.barcodes.tsv.gz`, written whenever the
  stage writes metrics (`--metrics`, or the stage's default metrics path).
  No new required flag. If you add an opt-out (e.g. `--no-barcode-metrics`),
  document it in MANUAL.md; `tests/docs_commands.rs` checks documented flags.
- Gzip TSV with a header. One row per observed barcode, sorted by barcode.
  First columns: `barcode` (as written to `CB`, including `-1` where the
  stage adds it), then `sample` (sample-sheet sample for SBX `extract`, `PB`
  for Flex in `molecules`, `.` otherwise), then the stage's counters.
- Reads that never got a barcode go in one aggregate row with
  `barcode = __no_barcode__`, so every column sums to the pool-level metric.
- Integers only; no derived ratios (reports compute those).

### Columns per stage

| Stage (binary) | Columns | Must sum to (pool metric) |
|---|---|---|
| `illumina` (`sbxpress-illumina`) | `reads`, `barcode_exact`, `barcode_corrected`, `emitted`, one `failed_<reason>` column per `failed_by_reason` key | `records_read`, `barcode_exact`, `barcode_hamming1_corrected`, `records_emitted`, `failed_by_reason` |
| `extract` (`sbx10x`, SBX) | `reads`, `cb_on_whitelist`, `cb_within1`, `emitted`, `class_full`, `class_partial`, `class_unparseable`; Flex adds `pb_located`, `pb_hit` | the matching `DemuxMetrics` / `WhitelistMetrics` / `ProbeBarcodeMetrics` fields in `metrics_summary.json` |
| `annotate` (`sbx10x-annotate`) | `records`, `primary`, `unmapped`, `exonic`, `intronic`, `intergenic`, `antisense`, `gene_assigned`, `unique_gene`, `multi_gene` | the same fields of annotate `Metrics` |
| `molecules` (`sbx10x-molecules`) | `eligible_reads`, `raw_cell_gene_umis`, `molecules_counted`, `collision_groups`, `representatives_tagged` | the same fields of molecules `Metrics` |

`matrix` and `cells` already produce per-barcode data; leave them alone.

## Where to hook in

Each stage already accumulates counters on worker-local structs that are
merged at the end. Add a barcode-keyed map to that same struct and merge it
in the same place. Don't add locks, channels or a second pass.

- **Key.** Pack the corrected 16-bp barcode to `u32` (2 bits/base; reuse
  `Dna::pack` or the barcode module's packing if present). Use
  `FxHashMap<u32, Counters>`; the crate already uses FxHash. `Counters` is a
  small `#[derive(Default, Clone, Copy)]` struct of `u64`s. SBX needs
  `(sample_index, cb)`; Flex in molecules needs `(cell, sample)`, the same
  scoping as `RawKey`.
- **`illumina`:** `src/illumina/pair.rs::convert_pair` increments
  `PairMetrics` per pair (`records_read`, `records_failed`,
  `failed_by_reason`). Add `per_barcode` to `PairMetrics`
  (`src/illumina/types.rs`) and merge it in `src/illumina/converter.rs` (~line
  189), where the other `PairMetrics` fields are summed. Write the TSV beside
  the metrics JSON.
- **`extract`:** `src/pipeline_runner/process.rs::process_one` updates
  `stats` per read (`cb_on_whitelist`, `class_*`, demux). Add the per-barcode
  map to that stats type, merge with the other worker stats
  (`src/pipeline_runner/engine.rs` / `worker.rs`), and write it next to
  `metrics_summary.json` in `src/pipeline_runner/report.rs`.
- **`annotate`:** `src/annotate/pipeline.rs::annotate_record` returns a
  per-record `RecordMetrics` (0/1 fields, `src/annotate/types.rs:229`) that
  is folded into `Metrics`. Read `CB` once per record (unordered and streaming
  paths both) and fold the same `RecordMetrics` into the barcode's counters.
  Records without `CB` go to `__no_barcode__`.
- **`molecules`:** count `eligible_reads` where `metrics.eligible_records += 1`
  in `src/molecules/evidence.rs` (serial ~line 168, indexed `collect_region`,
  and the merge at ~419). Count `molecules_counted` and `collision_groups` in
  `src/molecules/collisions.rs`, where the pool counters increment.
  `representatives_tagged` goes in `src/molecules/pipeline.rs` (~138). Sharded
  runs (`run_sharded`) merge each shard's map like the existing
  `Metrics::merge` (`src/molecules/types.rs` ~217). Write from `write_metrics`
  in `src/bin/sbx10x-molecules.rs`.

Keep files under the 500-line ceiling. If a module would cross it, put the
counter type and its TSV writer in a small new module (e.g.
`src/barcode_metrics.rs`) shared by all four stages.

## Pipeline scripts

`scripts/sbxpress-pipeline.sh` and `scripts/sbxpress-gex-pipeline.sh` need no
new flags, since the tables follow `--metrics`. Make sure each stage's
metrics path is set so the tables land in the run directory, and add the
`.barcodes.tsv.gz` files to the smoke test's expected outputs
(`tests/scripts/pipeline_smoke.sh`).

## Tests (required; see CLAUDE.md)

1. **Reconciliation, per stage.** On the existing fixtures, every
   `.barcodes.tsv.gz` column summed over rows equals the matching pool metric.
   This is the main correctness guard.
2. **Exact per-barcode values** on a tiny hand-checkable fixture for each
   stage (two or three barcodes, known reads per barcode, one barcode-less
   read).
3. **Parallel = serial.** `illumina --threads 1` vs `--threads 8`, annotate
   unordered vs streaming, molecules `--evidence-workers 1` vs `4` vs
   `--shards 4`: identical tables.
4. **Chain test.** In `tests/pipeline_chain_e2e.rs`, joining the four tables
   on the barcode gives a monotone non-increasing funnel per barcode
   (`reads ≥ emitted ≥ … ≥ molecules_counted`), and `molecules_counted` per
   barcode equals that barcode's column sum in the raw matrix.
5. **Golden outputs unchanged.** `tests/golden_outputs.rs` must pass
   untouched; this change adds files and changes no records.
6. Add new tests to `tests/test_inventory.txt` by inserting lines only.
   `scripts/check_test_inventory.sh --update` reorders the file under a
   different locale, so don't use it for this.

Run `scripts/check.sh` before every commit.

## Overhead check (acceptance)

Measure on the server with release builds, before and after, same inputs, in
a fresh output directory each time:

- `illumina` and `annotate`: DZEC05 30 min library inputs,
  `/mnt/data-ready/illumina/20261005_LH00826_0193_A2522CWLT4/DZEC05/DZEC05_SSS_30min_S3_L005_R{1,2}_001.fastq.gz`
  (111M pairs). Use the commands in
  `/mnt/analysis/active_projects/sbxpress_validation/dzec05_dedup/run_dzec05.sh`
  (`convert` and `annotate` functions).
- `molecules`: `/mnt/analysis/active_projects/sbxpress_validation/dzec05_dedup/30min/annotated.bam`
  with `--tagseq --evidence-workers 16`.
- `extract`: the CAT06 5-chunk SBX subset in
  `/mnt/analysis/active_projects/sbxpress_validation/v0.5.6_vs_v0.6.0/inputs/cat06_sbx`
  (Flex; see `docs/SERVER_TESTING.md` for the full command).

Run each twice and record `/usr/bin/time -v` wall time and peak RSS.
Acceptance: wall-time increase ≤ 2% and peak RSS increase ≤ 5% per stage
(well under these is expected). Cross-check one barcode's numbers against
the DZEC05 awk results in
`/mnt/analysis/active_projects/cat_tagseq/DZEC05/analysis/work/per_barcode_30min.tsv`.
Its columns are barcode, tagged reads, poly-G, mapped, gene-assigned,
LINC00486; tagged reads and mapped should match `annotate`'s `records`
(primary + unmapped) and mapped counts.

## Docs and release

- MANUAL.md: one short section listing each table, its columns and the
  sum-to-pool guarantee, plus one example join.
- CHANGELOG: `## [0.6.1]` → Added.
- `docs/SERVER_TESTING.md`: in the stage comparison, add the new tables to
  what is compared between builds.
- One version bump to 0.6.1 after `check.sh` passes, with tag `v0.6.1`.
  Don't move `v0.6.0` again.

## Part 1 done means

- Four stages write `.barcodes.tsv.gz` tables that sum exactly to their pool
  metrics, identical across thread and shard settings.
- Golden outputs unchanged; `check.sh` passes; overhead within the limits
  above, with the measurements recorded in the PR description.
- A one-paragraph note in the PR on how to join the tables into a
  per-barcode funnel.

# Part 2: read-quality filtering (bad-read handling)

## Why

sbxpress has no read-quality or mapping-uniqueness filter on the counting
path. MAPQ is used only to pick a molecule's representative read
(`src/molecules/helpers.rs:44`), and conversion checks only cDNA length and
uses base quality only when merging read pairs. DZEC05 showed what gets
through: 28–56% of read pairs had a dead read 2 (poly-G from cycle 1). Most
did not map, but 6.7M aligned to G runs of up to 81 bp inside LINC00486,
giving that gene 18–21% of all molecules (normally 0.03–0.14%). Those reads
are 99% MAPQ < 30, 53% carry spurious splice gaps, and all are intronic
(`RE:A:N`).

The fix should catch this and the wider class of bad reads (low-complexity
inserts, adapter or primer dimers, chimeras, multi-mapping artifacts), not
only poly-G.

## Measurements that constrain the design (2026-10-08)

- **Base quality does not identify dead reads.** NovaSeq X calls a no-signal
  G at fairly high quality: poly-G R2 reads average Q25, and only 20% have a
  mean below Q20 (6% of normal R2 reads do). A mean-quality filter removes a
  fifth of the problem.
- **A plain MAPQ threshold is too blunt.** A MAPQ ≥ 30 cut removes
  98.8–99.6% of LINC00486 reads, but also 18% of all other gene-assigned
  reads in DZEC04 (a normal library) and 28% in DZEC05. The losses are
  concentrated in highly expressed genes with processed pseudogenes or
  paralogs: EEF1A1 100%, RPL41 99%, STEAP1B 99%, RPLP0 98%, RPSA2 94%,
  PTMA 83%, FTH1 67%, MT-ATP6 38–49% (NUMTs), HBG2 38%, HBG1 28%. minimap2
  runs with `--secondary=no`, so sbxpress cannot tell "two copies of one
  gene's sequence" from "a G run that fits anywhere".
- **Cell Ranger's approach.** It counts only confidently mapped reads, with a
  MAPQ rescue: a multi-mapped read whose alignments include exactly one hit
  to an annotated exon is counted for that gene. Processed-pseudogene copies
  are unannotated or filtered, so EEF1A1 and RPL41 survive. LINC00486's reads
  are intronic multi-mappers, so they would not be rescued.

## Recommended design (layered; confirm decision points before coding)

**Layer A: sequence filters at conversion** (`illumina`, `extract`). Cheap,
works before mapping, catches the cause, and the counts are visible in
metrics.
1. Trim poly-G tails from the 3' end of cDNA reads (two-colour chemistry),
   then re-apply `--min-cdna-len`.
2. Low-complexity cDNA filter: drop a read if, after trimming, one base makes
   up ≥ 80% of it, or its longest homopolymer covers ≥ 50% of it (thresholds
   are flags). This catches poly-G, leftover poly-A/T and many primer-dimer
   or TSO concatemer products. Report it as `failed_by_reason.low_complexity`.
3. UMI sanity, as Cell Ranger does: drop reads whose UMI contains N or is a
   homopolymer; report `failed_by_reason.umi_invalid`. Leave a UMI
   base-quality threshold as an option, off by default.

**Layer B: gene-aware mapping uniqueness at counting** (`annotate` → `molecules`).
Replaces "MAPQ ≥ X" with the Cell Ranger rule.
1. Keep secondary alignments (minimap2 `-N`; `annotate --multimap-tags`
   already handles query-grouped multi-hits and NH/HI tags).
2. A read is countable if its best-scoring alignments resolve to exactly one
   gene, giving exonic hits priority over intronic and intronic over
   intergenic. If they tie across genes, or the only hits are intronic on
   more than one locus, mark it ambiguous and don't count it.
3. Keep `--min-mapq N` as an explicit option (default off) for users who
   want the simple rule.
4. Counters: `ambiguous_multigene`, `rescued_multimapper`, `low_mapq` (when
   the option is on).

**Layer C: alignment sanity, cheap, per read** (`annotate`).
1. Minimum aligned fraction of the read (e.g. ≥ 50%; soft-clipped
   chimeras and adapter read-through fail).
2. Splice support: an `N` gap must match an annotated junction or the
   `--junctions` BED, or the read is flagged. 53% of the LINC00486 reads had
   unsupported gaps.
3. Optional mismatch-density cap (NM / aligned bases).

**Layer D: monitoring** (no reads dropped).
1. Every filter's count appears in the stage metrics and in Part 1's
   per-barcode tables.
2. A per-gene artifact list in molecules metrics: genes where ≥ 50% of reads
   were low-complexity, ambiguous or unsupported-splice, with counts. On
   DZEC05 this would have named LINC00486 on day one.

## Decision points to confirm with the user

- Defaults on or off? Recommendation: Layer A and Layer C.1 on by default;
  Layer B on by default for tag-seq and GEX once validated (it changes counts
  for paralogous genes such as HBG1/HBG2, as Cell Ranger does); C.2 and C.3
  flag-only at first.
- Thresholds for low complexity and aligned fraction.
- Layer B needs secondaries from minimap2, which makes annotation slower and
  the BAM larger. Measure both before making it the default.

## Tests

- Unit tests per filter with hand-built reads: poly-G, poly-A tail, primer
  dimer, N-UMI, homopolymer UMI, soft-clipped chimera, unsupported splice,
  pseudogene multimapper (rescued), paralog multimapper (ambiguous).
- Golden outputs **will change** for any default-on filter. Regenerate in a
  separate commit that names the filter and the reason (CLAUDE.md rule).
- Chain test: per-barcode tables reconcile, with the new failure columns, to
  the stage metrics.

## Validation on real data (acceptance)

Run 0.6.0 and the new build on DZEC05 30 min, DZEC04 (subset and full depth),
DZEC03, CAT02 GEX and CAT06 Flex (procedures in `docs/SERVER_TESTING.md`)
and report:
- **DZEC05:** LINC00486 back to background (< 0.2% of molecules); dead R2
  pairs counted under `low_complexity`.
- **Normal libraries** (DZEC03, DZEC04, CAT02, CAT06): Layer A removes
  < 1% of reads; total molecules change < 2%; per-gene log2 ratio new/old
  within ±0.1 for ≥ 99% of genes with ≥ 100 molecules, with every larger
  change listed and explained (expected: paralogs under Layer B).
- **CAT02 against Cell Ranger:** concordance unchanged or better
  (`compare_tagseq_modes.py`).
- Runtime and peak memory per stage, before and after.

## Part 2 done means

Bad-read classes are removed or flagged with counts in every run's metrics
and per-barcode tables, LINC00486-type artifacts can't dominate a library,
normal libraries change only where the rules intend, and the defaults chosen
at the decision points are documented in MANUAL.md and CHANGELOG 0.6.1.
