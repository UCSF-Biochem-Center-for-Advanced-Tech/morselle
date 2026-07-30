# Mapping and annotation analysis: long SBX cDNA reads

## Question

In the M2 comparison, why does strict genomic uniqueness decrease as SBX cDNA
reads become longer? At first glance this is counterintuitive: more sequence
should generally make an alignment more specific.

The important distinction is between **mapping somewhere**, **one reported
genomic locus**, and **one compatible gene**. The observed effect is primarily a
decrease in `NH=1` (one reported alignment), not evidence that longer reads fail
to map or broadly lose gene assignment.

## What the current audit shows

The audit stratifies mapped primary BAM records by raw read length and separately
by aligned length. Across the SBX P5 and P7 M2 runs:

- The fraction with `NH=1` falls with increasing read length.
- Gene assignment (`GX`) is substantially more stable than `NH=1`.
- The same uniqueness trend remains in aligned-length bins. It is therefore not
  explained solely by longer unaligned sequence.
- SBX has markedly fewer heavily soft-clipped reads than Illumina in this
  comparison. Residual read-through adapter or poly(A/T) sequence is unlikely to
  be the principal explanation.
- Multi-gene assignments remain uncommon. In particular, many reads with two
  reported genomic hits still have exactly one compatible gene.

This combination supports the interpretation that the added sequence exposes an
additional genomic placement more often than it creates an additional plausible
gene assignment.

## Interpretation of relevant tags

| Quantity | Meaning in this analysis | What it does *not* mean |
|---|---|---|
| mapped | minimap2 produced a genomic alignment | the alignment is unique or countable |
| `NH=1` | minimap2 reported one alignment under the selected secondary-hit policy | one biologically unique gene |
| `GX` | annotation found one or more compatible genes | one genomic locus |
| `GX` degree = 1 | all compatible evidence supports one gene | minimap2 found only one locus |
| `mm=1` | multiple reported hits but one compatible gene in sbxpress multimap mode | a MAPQ promotion or a molecule decision |

The multimap run uses minimap2 `-N 10`: one primary plus up to ten secondary
alignments are retained. `sbxpress annotate --multimap-tags` processes each
query-group together and emits `NH`, `HI`, and `mm`. The audit top-codes the last
ambiguity bucket as `NH=11+`.

Thus `NH=1` is a property of the mapper's reported-hit policy as well as the
sequence and reference. It should not be interpreted as a direct measure of
gene-level identifiability.

## Plausible mechanisms

1. **Additional homologous genomic evidence.** Longer cDNA can extend into
   shared UTR, intronic, paralogous, pseudogene, or repeat-associated sequence.
   A second sufficiently good chain can then be reported even when the annotated
   gene is unchanged.
2. **Transcript/gene structure.** A longer read may span sequence that is less
   locus-specific than its initial exonic portion. Conversely, a splice junction
   may improve gene assignment without eliminating another genomic placement.
3. **Mapper reporting policy.** minimap2 emits secondary chains that pass its
   score and chaining thresholds. The frequency of `NH>1` depends on the index,
   splice preset, `-N`, and secondary-hit threshold; it is not interchangeable
   with STAR's confidence conventions.
4. **Composition confounding.** Long-read bins are not random extensions of
   short-read bins. They can be enriched for different transcript classes,
   intronic molecules, gene families, or molecule-quality subpopulations.

The current results do **not** support the simpler claim that longer reads map
worse. They support: longer reads retain or improve gene-identifying evidence,
while more often revealing additional reported genomic homology.

## Read-level churn: baseline vs. multimap, matched by read name (compact fixture)

The aggregate molecule-count delta between the baseline and multimap arms
(+153 of 5,798,640, 0.0026%; see `mapping_evaluation/`) could in principle mask
a much larger pair of opposing effects — reads newly excluded because a
secondary alignment introduced ambiguity, offset by reads newly rescued via
`mm`. This was checked directly: primary/unmapped records from both arms
(9,419,913 in each, confirming primary-alignment selection itself is unaffected
by `--secondary=no` vs `-N 10`) were joined by read name and their `GX` degree
compared.

| baseline → multimap | count | % of baseline pool |
|---|---|---|
| single-gene (`GX` degree 1) → lost (0 or 2+) | 501 | 0.009% of 5,892,858 |
| non-single (0 or 2+) → gained single-gene | 652 | 0.018% of 3,527,055 |
| net | +151 | matches the aggregate molecule delta |

Both pools are genuinely small — not two large effects cancelling. Two
sub-findings complicate a clean "multimap rescue" story, though:

- Of the 501 lost reads, only 26 became genuinely multi-gene (`GX` degree ≥2).
  The other 475 dropped to zero compatible genes. That is a different failure
  mode than "excluded for being multimapped" and should be checked in the
  annotator's `--multimap-tags` code path rather than attributed to secondary
  alignment ambiguity by default.
- Of the 652 gained reads, only 115 (17.6%) carry the `mm` tag (the actual
  multimap-rescue signature: multiple reported loci, one compatible gene). The
  remaining 537 gained single-gene status with `NH=1` in both arms — the
  change came from the `--multimap-tags` vs. `--unordered` annotation code
  path itself, not from seeing any secondary alignment.

**Limitation (partially addressed 2026-07-30):** this fixture's own
diagnostic-arm `NH=1` fraction is ~91-92% and roughly flat across length bins
— it barely exercises the length-driven ambiguity this document is about, so
a churn rate measured only here would not bound the churn rate on real SBX
reads. The full-scale P5/P7 diagnostic run (`scripts/run_full_evaluation.sh`,
launched 2026-07-30 06:38) is still in progress and will give the
production-scale answer. In the meantime this churn comparison was repeated on
a matched-scale (~9.5M read) subsample of real raw P5 and P7 data (first
10M read pairs of shard 001 from `p5_out_7-7`/`p7_out_7-7`, output in
`mapping_evaluation/p5_prelim/` and `p7_prelim/`) — not a full-scale or
randomly-drawn sample, but confirmed representative in read length (mean 233bp
P5, 254bp P7 vs. the full-scale medians of 226bp/244bp established elsewhere).

| dataset | reads | lost (single→non-single) | lost, →0 genes / →2+ genes | gained (non-single→single) | gained, `mm`-tagged | net Δ molecules |
|---|---|---|---|---|---|---|
| Illumina (fixture) | 9,419,913 | 501 (0.009%) | 475 / 26 | 652 (0.018%) | 115/652 (17.6%) | +153 (0.0026%) |
| SBX P5 (prelim) | 9,497,915 | 379 (0.007%) | 351 / 28 | 624 (0.014%) | 89/624 (14.3%) | +232 (0.0046%) |
| SBX P7 (prelim) | 9,513,521 | 284 (0.005%) | 263 / 21 | 591 (0.015%) | 87/591 (14.7%) | +303 (0.0055%) |

All three datasets agree on two points that go beyond "multimap changes almost
nothing":

1. **The net effect stays small and positive at this scale for all three.**
   The magnitude is small enough (hundreds of reads out of ~9.5M) that each
   individual estimate carries substantial noise. The net Δ ordering
   (Illumina fixture +153 < P5 +232 < P7 +303) does *not* track mean read
   length — the Illumina fixture has the longest tagged-R2 reads in this
   comparison (mean 316 bp vs. 233 bp P5 and 254 bp P7), yet produces the
   smallest net gain. The ordering may reflect orientation-driven differences
   in cDNA quality (P7 > P5 is consistent with the CB/UMI position effect
   documented elsewhere), or may simply be noise at these event counts. A
   proper significance test and the full-scale P5/P7 run are needed before
   treating any trend here as real.
2. **Most of the churn on both sides is not actually about secondary
   alignments.** In every dataset, ~93-95% of "losses" drop to *zero*
   compatible genes rather than becoming genuinely multi-gene, and only
   ~14-18% of "gains" carry the `mm` tag (the actual multimap-rescue
   signature: multiple reported loci, one compatible gene). The remainder on
   both sides changes with `NH=1` unchanged in both arms — i.e., it comes from
   minimap2 choosing a slightly different primary alignment between
   `--secondary=no` and `-N 10` runs (documented mapper behaviour when the
   candidate chain set changes), not from a divergence in the annotation logic.
   Both code paths call the same `annotate_record` function for individual
   records. Read another way: on any of these three datasets, the *primary
   alignment selection variance* between the two minimap2 runs is larger than
   the *genuine multimap-rescue signal* (21-28 true ambiguity losses and 87-115
   true rescues per ~9.5M reads).

Still to do: re-run this same comparison on the full-scale P5/P7 diagnostic
output once `run_full_evaluation.sh` completes, both to get a proper random
(not shard-prefix) sample and to see whether the length-vs-net-effect trend
holds at full statistical power.

## Current controls and limitations

The aligned-length analysis and soft-clip diagnostics reduce the likelihood that
adapter/poly(A/T) carryover explains the pattern. They do not identify the source
of each secondary alignment.

The currently available BAMs do not consistently provide a primary-secondary
alignment-score comparison suitable for the report. Therefore the present audit
cannot yet distinguish a nearly tied secondary hit from a weak alternative chain.
Also, aggregate length bins cannot separate a true length effect from changes in
the genes represented at different lengths.

## Evaluation data and reference set

Use CAT02 M2 GEX as the first controlled cohort. The compact test fixtures are
useful for matrix regression but cannot answer mapper-ambiguity questions: their
BAMs contain no secondary alignments.

| Purpose | Location |
|---|---|
| Controlled M2 input | `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/example_fastqs/example_CAT02_fq/M2_GEX/` |
| Sample sheet | `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/p5_out_7-7/sample_sheet.csv` |
| ARC whitelist | `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/737K-arc-v1.txt.gz` |
| Splice index | `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi` |
| GTF | `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/gencode.v44.primary_assembly.basic.annotation.gtf.gz` |
| Junction BED | `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/gencode.v44.primary_assembly.junctions.bed` |
| Canonical features | `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/genes.mtx` |
| Existing comparators | `CAT02/illumina_data/cellranger_outputs/M2_GEX/outs/`, `CAT02/p5_out_7-7/cellranger_outputs/M2_GEX/outs/`, and `CAT02/p7_out_7-7/cellranger_outputs/M2_GEX/outs/` |
| Matrix regression fixtures | `/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/test_fixtures/` |

The M2 example directory has paired, extracted R1/R2 data across three lanes
(about 2.7 million pairs total). It is appropriate for the shared tagged-cDNA
mapping path and controlled minimap2 comparison. It is not raw SBX input. Use
`raw_fastqs_run_5p/` or `raw_fastqs_run1/` only after the small evaluation has
passed; their individual input files are several gigabytes.

## Controlled evaluation plan

These analyses should precede any SBX-specific mapping policy change.

### 1. Freeze the cohort and baseline

- Start with M2 GEX and retain the same read pairs in every arm. If downsampling
  is needed, use a deterministic hash or recorded random seed and preserve lane
  balance.
- Record exact input paths, read-pair count, sbxpress/minimap2 versions, index
  checksum, GTF checksum, junction BED checksum, and all commands in the result
  directory.
- Define length bins before analysis, and make the primary mapped alignment the
  read-level observation. Use deduplicated `(CB, UMI, gene)` records for
  molecule-level analyses so repeated reads do not inflate inference.
- First run `tests/validate_cr_bam_to_mtx.sh` against the CAT02 fixture directory
  to confirm the current binary still reconstructs the expected matrices.

### 2. Produce two minimap2 arms

Both arms must use the reference set above, `-x splice:sr`, `--junc-bed`, the
same tagged cDNA input, and the same annotation and molecule configuration.

| Arm | Mapper policy | Annotation policy | Valid endpoint |
|---|---|---|---|
| Production baseline | `--secondary=no` | `--unordered` | Throughput and downstream gene, molecule, matrix, and cell outputs. Do not interpret `NH` here. |
| Diagnostic multimap | `-N 10` | `--multimap-tags` | Reported-hit accounting with `NH`, `HI`, and `mm`, plus downstream impact. |

In the diagnostic arm, `NH` is the number of reported mapped primary plus
secondary records. It is capped by the primary plus ten secondaries; `NH=11`
therefore means 11 or more reportable hits under this policy. It is not a count
of all possible genomic placements. The fast arm suppresses secondary records,
so its apparent `NH=1` rate is not comparable to the diagnostic arm.

Keep the annotated multimap BAM as the evidence artifact. The molecule stage
intentionally ignores secondary and supplementary records, so primary-secondary
relationships cannot be reconstructed from the molecule BAM.

### 3. Analyze mapping ambiguity

For primary records in the diagnostic BAM:

1. **Within-gene stratification**
   - Compare `NH=1` by length within `GX`, transcript biotype, exon/intron class,
     and gene-length strata.
   - Report both unadjusted bin summaries and a molecule-aware adjusted analysis.
     This separates read length from transcript composition and read depth.
2. **Classify secondary hits**
   - For every `NH>1` query, classify each secondary record as compatible with the
     primary gene, another gene/paralog or pseudogene, repeat/decoy sequence, or
     unannotated/intergenic sequence.
   - Retain query name, primary and secondary contig/interval/strand, `GX`, `NH`,
     `HI`, and `mm`; summarize the classes by primary read-length bin.
3. **Compare alignment evidence**
   - Record available primary/secondary alignment scores, query-aligned bases,
     query coverage, identity, clipping, MAPQ, and genomic overlap.
   - Plot score-gap and coverage-gap distributions by length and secondary-hit
     class. Explicitly label unavailable mapper tags rather than treating them as
     zero.

### 4. Measure downstream impact

Run annotation, molecules, matrix construction, and cell calling for both arms.
Compare:

- mapped and annotated primary reads, `GX` degree, and `mm` frequency;
- retained and counted molecules, including changes for the same `(CB, UMI)`;
- raw and filtered cell counts, recovered molecules, and called-cell overlap;
- raw and filtered matrix totals, per-cell and per-gene correlations, and genes
  whose counts materially change; and
- wall time, peak memory, BAM size, and storage required by each arm.

Use the existing Illumina, P5, and P7 M2 outputs as external concordance
comparators for genes, molecules, and matrices. Do not compare their mapper
specific `NH` or MAPQ conventions directly with minimap2.

### 5. Decision criteria and scale-up

Keep the production fast path unless the diagnostic arm produces a defined,
reproducible improvement in gene- or molecule-level recovery. A higher `NH=1`
fraction alone is not a success criterion. Before changing policy, document:

- the fraction and identity of secondary-hit classes;
- whether length remains associated with ambiguity after composition adjustment;
- the magnitude and direction of molecule, cell, and matrix differences; and
- the runtime and storage cost of retaining secondary alignments.

Only then repeat the selected protocol on a larger deterministic CAT02 subset or
the raw P5/P7 SBX inputs and confirm that the conclusion holds across both
orientations.

## Implication for sbxpress

The default production path may continue to suppress secondary alignments and
use unordered annotation for throughput. The multimap path is a diagnostic and
validation mode: it preserves query-level hit context before molecule calling.

Any SBX-specific change—such as gene-aware multimap rescue, a revised mapping
filter, or a length-aware policy—should be evaluated only after the baseline
pipeline is validated for annotation tags, molecule representatives, matrix
construction, and cell calling. The relevant success criterion is improved
gene- and molecule-level recovery with a measured, explainable effect on
ambiguity, not an increase in `NH=1` alone.
