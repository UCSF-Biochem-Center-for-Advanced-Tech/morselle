# sbxpress 0.6.1 per-barcode metrics; 0.6.2 read-quality filtering

## Implementation checkpoint (2026-10-08)

Part 1 is committed on `feat/per-barcode-metrics` as 0.6.1 (not tagged or
pushed). Part 2 has not been started. Its Layer B design was revised after
reading Cell Ranger's rescue code and measuring the rule on DZEC04/DZEC05,
and a read-end overhang tolerance was added (Layer E); see Part 2 below.

Changed from the brief below: the tables are **opt-in** (`--barcode-metrics`
on every stage and both pipeline scripts), because the first implementation
cost `extract` 11.8 GB peak RSS on Flex (one table per worker thread) and
`molecules` ~5% wall time. `extract` now shares one sharded table and
`molecules` keys counters by packed (cell, sample). Measured on real-data
subsets (3 alternating runs): flag off, no measurable change from 0.6.0;
flag on, `annotate` +7%, `extract` +12% and +375 MB on one CAT06 chunk.

Not yet validated on real data (run before tagging): table sums against pool
metrics for every stage; identical tables across repeat runs and across
thread/shard settings; flag-off outputs byte-identical to 0.6.0; the DZEC05
one-barcode cross-check; full 5-chunk Flex `extract` memory with the flag;
one real pipeline-script run with `--barcode-metrics`.

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
- Reads whose barcode could not be located go in `__no_barcode__`. Reads
  whose raw barcode was located but cannot be uniquely corrected go in
  `__unresolved_barcode__`. Both are aggregate rows; neither denotes a cell.
  A stage may use `__no_barcode__` for records with no usable `CB` tag.
  Every input read contributes to exactly one row's input counter.
- Join identified rows on `(barcode, sample)` when the sample identity is
  consistent across stages. For SBX/Flex, sample-sheet assignment and `PB`
  identify different things: retain both as distinct fields where needed,
  and join within a single routed sample/library using `CB` plus `PB` for
  Flex molecules. Never merge Flex `PB` groups sharing the same `CB`.
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

- **Key.** Preserve the emitted `CB` spelling, including the GEM suffix, and
  scope SBX by the routed sample and Flex molecules by `PB`, as `RawKey` does.
  A packed 16-bp barcode can be used internally where it preserves those
  identities. Avoid per-read string allocations in the hot path.
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
4. **Chain tests.** Test the Illumina → annotate → molecules and SBX extract
   → annotate → molecules routes separately. Flex uses probeset-align instead
   of annotate, so test extract → molecules there. Join using the identity
   rules above. On the default `--secondary=no` route, compare primary read
   counts: `reads ≥ emitted ≥ primary ≥ eligible_reads ≥ molecules_counted`.
   In multimap mode, `annotate` currently includes secondary hits in its
   `primary` metric, so test its table against its pool metrics without that
   cross-stage inequality. Supplementary alignments are additional BAM
   records, not additional input reads. The per-barcode molecule count
   must agree with the raw matrix column sum within the matrix's existing
   floating-point tolerance when fractional collisions are enabled.
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
- `extract`: the CAT06 5-chunk SBX subset, raw chunks
  `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT06/raw_fastqs/20260901_UCSF_Cat06_843_PT370_C3-000300.fastq.gz`
  through `-000304.fastq.gz` (symlink them into one directory for `--fastqs`).
  It's Flex; the full command is in `docs/SERVER_TESTING.md`. The 0.5.6/0.6.0
  validation tree was deleted after acceptance; its evidence is in
  `/mnt/analysis/active_projects/sbxpress_validation/v0.5.6_vs_v0.6.0_evidence.tar.gz`.
- Install the before/after builds into their own prefixes (`install.sh --prefix
  … --skip-packages --skip-rust`). `dzec05_dedup/run_dzec05.sh` points at a
  0.6.0 install that no longer exists, so take only its commands from it.

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

**Layer B: Cell Ranger MAPQ rule with exonic rescue** (`annotate` → `molecules`).
Agreed design (2026-10-09), replacing "MAPQ ≥ X" and the loose 0.3.5
`--multimap-tags` rule.

What Cell Ranger does (10.0.0, `tx_annotation/src/read.rs`): a read is
counted if it is not secondary, has CB and UB, is confidently mapped
(MAPQ == 255 for STAR, 60 for minimap2) and GX names exactly one gene.
`rescue_alignments_se` only runs when no alignment is confidently mapped:
of the alignments that are splice-compatible, sense and exonic (`tx_align`),
if they name exactly one gene the first is made primary with MAPQ = high
confidence and `mm:1`; the rest become secondary with MAPQ 0. Intronic
hits are never rescued, even with `include-introns`. NH is not read; HI is
a STAR pass-through.

1. **Default route first** (`--secondary=no`, one record per read; no
   sorting or grouping change). `annotate` sets `mm:i:1` on a primary read
   with MAPQ below the cutoff whose GX names one gene from a compatible
   sense exonic transcript (TX of the form `ENST…,+pos,cigar`). Annotate's
   GX is already splice-compatible (`compatible_position`), equivalent to
   Cell Ranger's `tx_align`, so no new check is needed. Intronic GX is never
   rescued, regardless of `--include-introns`.
2. **Keep minimap2's MAPQ** and flags unchanged; `mm:1` marks the rescue
   (Cell Ranger rewrites MAPQ to 255/60; we keep the original visible).
3. **Grouped route later** (`-N` with `--multimap-tags`, query-grouped):
   one record per read is written. Secondaries are dropped (or written to a
   `--secondary-output` file); the rescued hit is promoted to primary with
   `mm:1`; `NH` records the number of hits as a diagnostic. HI is not
   written (not meaningful outside STAR).
4. **`molecules` uses the tags** (next edit, see below).
5. Counters: `rescued_multimapper` (mm:1), `low_mapq_intronic` /
   `low_mapq_unrescued` (dropped), `multimapper_ambiguous` (NH > 1 without
   mm, grouped route), in stage metrics and Part 1's per-barcode tables.

**`molecules` uses NH and mm** (next edit). Today `molecules` and the matrix
never read NH or mm, so the 0.3.5 tags have no effect on counts (the CAT02
"exact copy of Cell Ranger" check ran `cr-bam-to-mtx` on Cell Ranger's own
BAM, using its `xf` decisions, not our molecules path).
- A record is countable if it is primary, has CB and UB, GX names one gene,
  and **MAPQ ≥ `--min-mapq` (default 60) or `mm:i:1`**.
- `NH` names the reason a low-MAPQ read is dropped (`NH` > 1:
  multimapper) but never vetoes a confidently mapped read: Cell Ranger
  counts MAPQ 60 whatever its other hits (decided 2026-10-09).
- `--min-mapq 0` restores today's behaviour (flag for comparison runs).
- New counters in molecules metrics and per-barcode table:
  `rescued_records`, `low_mapq_records`, `multimapper_records`.
- Golden outputs change; regenerate in a separate commit. Chain test in
  `tests/pipeline_chain_e2e.rs`: annotate a minimap2-style BAM with a
  MAPQ-60 read, a low-MAPQ compatible exonic read (rescued), a low-MAPQ
  intronic read (dropped) and an NH:2 read without mm (dropped), through
  molecules and matrix. **Implemented** (annotate flag `--confident-mapq`,
  molecules `--min-mapq`; `tests/regressions/mapping_quality.rs`).

Measured with the default-route rule (DZEC04 = normal library, DZEC05 30 min
= poly-G library; chr1/2/11/17, 25% sample, primary reads with CB/UB):

| | DZEC04 | DZEC05 |
|---|---|---|
| single-gene reads | 45.9 M | 1.18 M |
| MAPQ 60, counted | 78.3% | 26.5% |
| MAPQ < 60, compatible exonic, rescued | 19.0% | 22.5% |
| MAPQ < 60, intronic, dropped | 2.8% | 51.0% |
| LINC00486 share before → after | 0.245% → 0.000% | 48.5% → 0.05% |

The rescue is what makes a 60 cutoff usable (without it ~22% of DZEC04
would be lost; cf. MAPQ ≥ 30 above). Reads dropped in DZEC04 are intronic
and mostly not low-complexity (0–7% fail the Layer A test) but carry longer
homopolymers (longest run 10–15% of the read vs 3–5% kept): intronic
repeat / internal-priming reads, also dropped by Cell Ranger. Largest real
genes affected: AMBRA1 34%, NBEAL1 88%, PLEKHA7 64%, ZMYM4 49%; highly
expressed exonic genes lose almost nothing (SFTPB 2 reads). The 2.8% read
loss is near the < 2% molecule target; measure molecules, which Layer E
more than offsets.

**Layer E: read-end exon overhang tolerance** (`annotate`, separate commit).
7.0% of DZEC04 records (3.7 M in the subset) are `RE:E` with no gene: they
overlap an exon but fail `compatible_position`. Split by cause:

| cause | share | MAPQ 60 |
|---|---|---|
| straddles exon/intron, unspliced | 42.5% | 85% |
| straddles exon/intron, spliced | 20.4% | 80% |
| runs past gene 3' end | 11.4% | 69% |
| runs past gene 5' end | 0.6% | 52% |
| inside exons, junction mismatch | 6.1% | 84% |
| exon of an antisense gene only | 16.5% | 73% |
| antisense exon inside a sense gene body | 2.4% | 67% |

55% overhang an exon by only 1–3 bp. In a 3% sample (157 K reads with
≤ 10 bp overhang) 98.8% of overhangs are at a read end, not an internal
junction; at ≤ 3 bp the overhanging bases match the intron 100% and also
match the adjacent exon in 98.1% (minimap2 extends rather than splice a
1–3 bp anchor). These are mature-mRNA reads. At 4–10 bp only 68% match the
adjacent exon and 54% have no adjacent exon (transcript-end extensions).

1. In `compatible_position`, let the **first and last block** overhang their
   exon by ≤ `--max-end-overhang` bases (default 3; 0 = today's behaviour),
   treating the overhang as soft-clipped. Clamp the transcript position to
   the exon start so TX positions never go negative. No junction tolerance
   (junction wobble is 1% of cases).
2. Expected gain: ~2.0 M reads in the DZEC04 subset (~4% of records, 85%
   MAPQ 60), more than offsetting Layer B's loss.
3. Counter `end_overhang_rescued`. Golden outputs change (separate commit);
   unit tests in `src/annotate/tests.rs` (1, 3 and 4 bp overhangs, plus and
   minus strand, position clamp) and a chain case.
4. Later, optional: Cell Ranger-style intronic fallback for sense reads that
   overlap an exon but stay incompatible (> 3 bp straddles, 3' extensions,
   junction mismatches, antisense exon inside a sense gene body) when
   `--include-introns` is on. Antisense-exon-only reads stay uncounted, as in
   Cell Ranger.

**Build order:** Layer B on the default route (annotate `mm` + molecules
NH/mm rule) → Layer E → Layer A/C/D → grouped route.

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
- Layer B (default route) needs no secondaries; decided: keep MAPQ, add
  `mm`, cutoff 60. The grouped route (`-N`) needs secondaries, which makes
  annotation slower and the BAM larger; measure before adding it.
- Layer E default overhang (3 bp proposed) and whether to add the intronic
  fallback.

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

Part 2 is a separate 0.6.2 release. Its behavior-changing defaults and
golden-output update must not be included in the Part 1 metrics release.
Multi-mapper rescue must be validated specifically on intronic artifacts;
the per-gene artifact list can only attribute reads after they reach gene
annotation, unless a separate diagnostic mapping pass is specified.

## Part 2 done means

Bad-read classes are removed or flagged with counts in every run's metrics
and per-barcode tables, LINC00486-type artifacts can't dominate a library,
normal libraries change only where the rules intend, and the defaults chosen
at the decision points are documented in MANUAL.md and CHANGELOG 0.6.2.
