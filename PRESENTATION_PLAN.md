# sbxpress presentation brief and planning guide

## Purpose

Prepare a technical presentation describing **sbxpress**: a transparent,
minimap2-based path from Roche SBX or Illumina GEX reads to 10x-style expression
outputs.

The presentation must make a strict distinction between:

1. **Baseline compatibility work** — reproduce the relevant 10x-style
   decisions and outputs with an inspectable pipeline.
2. **SBX-specific optimization** — only after baseline behavior is validated,
   test whether longer SBX reads improve assignment or molecule recovery.

The presentation should be written for technical collaborators deciding how to
validate and evolve the pipeline. It should not imply that sbxpress is already a
drop-in replacement for every reference workflow decision.

## Core message

> sbxpress replaces the opaque short-read mapping/assignment path with explicit
> minimap2, annotation, molecule, matrix, and cell-calling stages. Its first
> objective is a trusted 10x-style baseline; only then should longer SBX cDNA be
> used to change decisions deliberately.

## Why this work exists

- Roche SBX produces an end-to-end library molecule as one long read, whereas
  the standard GEX workflow expects short paired reads.
- The conventional short-read path uses STAR and fixed short-read assumptions;
  it provides limited visibility into mapper, annotation, and molecule decisions.
- SBX cDNA can retain considerably more sequence. minimap2 enables splice-aware
  alignment while preserving barcode, UMI, and SBX processing provenance in
  the same record stream.
- The goal is not merely a different mapper. The pipeline must make its
  decisions inspectable and comparable at read, molecule, matrix, and cell
  levels.

## Pipeline to explain

Use this high-level flow:

```text
Raw SBX reads ─┐
               ├─> tagged cDNA FASTQ ─> minimap2 ─> annotation ─> molecules ─> matrix ─> cells
Illumina R1/R2 ┘
```

### Input / tagged FASTQ handoff

- `sbxpress extract` orients each SBX molecule to canonical P5→P7 order,
  extracts cell barcode (`CB`) and UMI (`UB`), and emits tagged cDNA R2 FASTQ.
- P5-led means the raw SBX read begins at the P5 end and reaches CB/UMI near the
  beginning. P7-led reads are reverse-complemented to the same canonical frame.
- `sbxpress illumina` converts standard paired Illumina R1/R2 into the same
  tagged cDNA handoff. R1 supplies CB/UMI; R2 remains cDNA.
- `sbxpress minimap -y` carries tags from the FASTQ comment into SAM/BAM.

### Mapping algorithm

- minimap2 uses the splice short-read preset with the chosen `.mmi` index and
  junction BED.
- **Production count mode:** `--secondary=no`; one alignment path per query,
  permitting unordered parallel annotation.
- **Multimap diagnostic mode:** `-N 10`; retain the primary plus up to ten
  secondary alignments. Query-grouped annotation adds `NH`, `HI`, and `mm`.
- `NH=1` means one **reported genomic alignment under this mapper policy**. It
  is not equivalent to one biologically unique gene.

### Annotation algorithm

- Parse the matching GTF into gene, transcript, exon, and intron intervals.
- Intersect CIGAR alignment blocks with those intervals.
- Apply the 10x GENCODE biotype, readthrough, PAR, and identifier-normalization
  rules.
- Emit `GX`, `GN`, `TX`, `AN`, and `RE`, plus compatible assignment counts
  (`XC`, `XT`).
- For 3′ assays, `--three-prime` creates terminal-support tags (`ZG`, `ZT`,
  `XP`) to break ties among compatible assignments. It must not discard
  non-terminal alignments.
- Annotation records feature evidence only; it does not correct UMIs or select
  counted molecule representatives.

### Molecule / matrix / cell layers

- `sbxpress molecules` performs directional within-cell/gene UMI correction,
  resolves corrected-UMI cross-gene collisions by support, and sets `xf=25` on
  one representative per counted molecule (`xf=17` on other eligible reads).
- `sbxpress matrix` counts representatives into a stable `genes.mtx` feature
  order and writes raw MEX/HDF5 matrices.
- `sbxpress cells` applies OrdMag and multinomial EmptyDrops-style cell calling
  and writes filtered MEX/HDF5 matrices.

## Required validation slide

Use one concise slide—not several—showing status and remaining work.

| Layer | Validation target | Status |
|---|---|---|
| BAM → matrix | Feature order, barcode order, sparse entries, MEX/HDF5 layout | **Exact match validated** |
| Mapping / annotation | Tag assignment, multimap behavior, gene/transcript evidence | Pending |
| Molecule calling | `CB + UB + GX` groups, UMI correction, collision resolution, `xf` representative selection | Pending |
| Cell calling | Called barcode set and filtered per-cell matrix values | Pending |

Acceptance should be exact equality where the semantics are deterministic (for
example BAM→matrix), and a documented explanation for unavoidable mapper or
algorithm differences elsewhere.

## Mapping analysis source

Use [mapping.md](mapping.md) as the authoritative source for the M2
read-length/uniqueness observation, tag definitions, controls, limitations, and
proposed follow-up analyses. Do not duplicate or extend its interpretation in
the presentation plan without updating `mapping.md` first.

For the presentation, summarize only its central conclusion: the observed decline
in strict locus uniqueness (`NH=1`) with longer SBX cDNA is not evidence that
long reads map worse. It may reflect additional reported genomic homology while
gene assignment remains comparatively stable.

## Mapping-ambiguity slides (presentable now)

The following slides are based on completed analyses and can be included before
the full-scale P5/P7 diagnostic run finishes. Each is scoped to the data that
actually exists; framing notes are mandatory where stated.

### Slide M1 — The central paradox (full-scale P5/P7 production BAMs)

**Source:** `mapping_evaluation/existing_gx_analysis/` (1-in-10 systematic
sample, ~61M primary records across P5 and P7 production BAMs).

**Content:** Two-panel figure or table, P5 and P7 side by side.

| Length bin | P5 GX=1% | P5 MAPQ | P7 GX=1% | P7 MAPQ |
|---|---|---|---|---|
| ≤250 bp | 60.0% | 52.1 | 62.3% | 53.3 |
| 251–300 bp | 55.5% | 54.1 | 57.6% | 54.7 |
| 301–350 bp | 53.4% | 55.0 | 55.7% | 55.6 |
| 351–400 bp | 52.5% | 55.4 | 55.3% | 56.4 |
| ≥401 bp | 51.6% | 54.8 | 55.7% | 56.8 |

**Key sentence:** Longer reads align more confidently (MAPQ rises) but land
in annotated gene territory less often (GX=1% falls) — the paradox that
motivates the multimap diagnostic.

**Mandatory framing:** These are baseline (`--secondary=no`) BAMs; `NH=1`
is trivially 100% for all records here and is not shown. The GX=1% decline
is a separate signal from the `NH=1` question in `mapping.md` and is not yet
decomposed by secondary-hit class at this scale.

### Slide M2 — P5 vs P7 orientation effect at long lengths

**Source:** same as M1.

**Content:** Single data point, can live as a call-out on slide M1 or its
own panel.

> P5 soft-clip fraction reaches **15.6%** at ≥401 bp; P7 stays at **5.2%**.
> At ≤250 bp the two orientations are indistinguishable (3.1% vs. 3.3%).

**Key sentence:** The orientation effect is length-dependent — it only
surfaces in the longest-read tail, consistent with signal degradation
accumulating over read position (the segment decoded second degrades more,
and that shows up disproportionately in long reads).

**Connection:** Links directly to the CB/UMI position slide (slide 6 in the
cross-platform section). Forward-reference it: "the same mechanism that
explains P7's better barcode accuracy also explains P5's excess clipping in
long reads."

**Mandatory framing:** Confirm before this slide ships that P5/P7 here are
the same molecule read in two orientations, not 10x 5'/3' capture-chemistry
variants (different GTF exon coverage would create a superficially similar
pattern for unrelated reasons).

### Slide M3 — Where secondary hits actually go (M2 compact cohort)

**Source:** `mapping_evaluation/analysis/secondary_classes.tsv` and
`score_gaps.tsv`. M2 compact cohort only (~9.4M reads, 250–650 bp range,
processed via `sbxpress illumina` path).

**Content:** Stacked bar chart of secondary-hit class by length bin, with
a small inset table of AS score gaps.

| Secondary class | Typical fraction across bins | Meaning |
|---|---|---|
| `unannotated_mapped` | 55–58% | Genomic hit, no annotated gene — pseudogene, repeat, annotation gap |
| `diff_gene` | 20–22% | Hit lands on a different gene — potential paralog/shared exon |
| `decoy_other` | 19–27% | Random/decoy chromosome — near-zero AS gap (tied placement) |
| `same_gene` | 1–2% | Secondary confirms primary gene — `mm` rescue candidate |

Score gap (primary AS − secondary AS) by class:
- `diff_gene`: median gap grows with length (6 → 9 bp from ≤250 to ≥401);
  longer reads are *better differentiated* from alternative gene placements.
- `decoy_other`: median gap ≈ 0 across all bins; these are tied placements
  to low-complexity/repeat sequence, not a length effect.
- `same_gene`: median gap 0–3 bp (nearly tied, as expected for the same locus).

**Key sentence:** When a longer read does pick up a secondary hit, it almost
never lands on a competing annotated gene — it lands in unannotated genomic
space or on a decoy contig. Multi-gene ambiguity (`diff_gene` + `same_gene`)
accounts for only ~22–24% of secondaries and, for `diff_gene`, the primary
is increasingly better-supported as reads get longer.

**Mandatory framing:** *"M2 compact cohort (sbxpress illumina path,
9.4M reads, 250–650 bp). NH=1 fraction is roughly flat across length bins
in this cohort (91–92%); the full-scale P5/P7 diagnostic run — in progress
— will test whether secondary-hit composition shifts when NH=1 does decline
with length."*

### Slide M4 — Multimap mode changes almost nothing downstream (M2 compact)

**Source:** `mapping_evaluation/` baseline vs multimap arms; churn table
in `mapping.md`.

**Content:** Two-column comparison table plus the 3-dataset churn summary.

| Metric | Baseline | Multimap | Δ |
|---|---|---|---|
| Molecules counted | 5,798,640 | 5,798,793 | +153 (0.003%) |
| Cells called | 26,303 | 26,305 | +2 |
| Filtered UMI total | 7,676,750 | 7,677,264 | +514 (0.007%) |

Reads that change GX assignment (read-name-matched join, all 3 datasets):

| Dataset | Lost single-gene | Gained single-gene | Net |
|---|---|---|---|
| M2 compact (sbxpress illumina) | 501 (0.009%) | 652 (0.018%) | +151 |
| SBX P5 prelim (~9.5M reads) | 379 (0.007%) | 624 (0.014%) | +232 |
| SBX P7 prelim (~9.5M reads) | 284 (0.005%) | 591 (0.015%) | +303 |

**Key sentence:** At all three scales tested, enabling the multimap diagnostic
path changes fewer than 0.02% of reads in either direction and adds a net of
~150–300 molecules per ~9.5M reads. Both the gain pool and the loss pool are
genuinely small — this is not two large effects cancelling.

**Mandatory framing:** *"Preliminary; all three datasets are ~9.5M read
subsets, not the full production scale. The P5/P7 full-scale run (in
progress) will confirm whether these rates hold at 1.2B reads per
orientation. These results support keeping `--secondary=no` as the production
default; they do not yet test whether multimap rescue improves recovery on
reads in the long-length tail where NH=1 is expected to decline."*

## Mapping-evaluation status (as of 2026-07-30)

Current state of the controlled evaluation plan from `mapping.md`, run via
`scripts/run_mapping_evaluation.sh` and `scripts/run_full_evaluation.sh`,
output in `CAT02/mapping_evaluation/`:

**1. Compact-fixture validation, complete (`baseline/` vs `multimap/`, Illumina
path, 10,495,177 R2 read pairs from the CAT02 M2 example fixtures).** This is
step 1-2 of the mapping.md plan, not yet the SBX-specific question — it
validates that the diagnostic multimap arm (`-N 10`, `--multimap-tags`) doesn't
change molecule-level output on its own. Result: `molecules_counted` differs by
153 out of 5,798,640 (0.0026%) between baseline and multimap arms. Per
mapping.md's decision criteria, this supports keeping `--secondary=no` as the
production default — the diagnostic arm shows no material recovery difference
at this scale. This is a first validation point only; it does not yet test
long SBX reads, where the ambiguity effect is expected to matter more.
Read-name-matched churn behind that net delta (does it hide offsetting
gain/loss pools?) is now checked and documented in `mapping.md` — both pools
are genuinely small (hundreds of reads out of 9.4M, not a large cancellation).
(`compare_downstream.py` was reading wrong metrics-JSON keys and silently
returning NaN for all sbxpress arm columns; fixed 2026-07-30 to use the actual
keys `primary`, `gene_degree_1`, `molecules_counted`, etc.)

**3b. Preliminary (not full-scale) P5/P7 churn data, complete.** Rather than
wait on the full-scale run, the same baseline+multimap two-arm protocol and
churn comparison was run on a ~9.5M-read prefix subsample of real raw P5/P7
data (`mapping_evaluation/p5_prelim/`, `p7_prelim/`; script:
`/tmp/.../scratchpad/run_prelim_arm.sh`, a parameterized copy of
`run_mapping_evaluation.sh` not committed to the repo). Confirmed
representative in read length (233bp/254bp mean vs. 226bp/244bp full-scale).
All 3 datasets (Illumina fixture, SBX P5 prelim, SBX P7 prelim) now agree that
(a) net churn stays small and positive, increasing somewhat with read length,
and (b) most of the churn on both sides — not just Illumina's — comes from the
`--multimap-tags` vs. `--unordered` annotation code path itself rather than
from genuine secondary-alignment resolution; see `mapping.md` for the full
3-way table and the resulting annotator follow-up. This is still a shard-prefix
subsample, not a random draw or the full-scale run — re-confirm on
`run_full_evaluation.sh`'s output once it completes.

**2. Secondary-hit classification, complete on the same compact fixture
(`analysis/summary.txt`, `analysis/secondary_classes.tsv`).** Adds resolution
mapping.md flagged as missing ("does not identify the source of each secondary
alignment"). By length bin, secondary hits classify as: unannotated_mapped
49-58%, decoy_other 18-27%, diff_gene 20-23%, same_gene only 1.2-1.7%. `GX`
(compatible gene found) declines 65.9%→62.6% ≤250bp→≥401bp while `NH=1` stays
roughly flat (91.3%→92.4%) and `mm%` (multi-hit-but-single-compatible-gene)
stays small (2.6%→2.3%). This is direct evidence for mapping.md's mechanism 1
(longer reads pick up additional non-genic/homologous genomic hits, rarely a
truly competing gene) — worth pulling into mapping.md itself as supporting data
before it's cited in the deck, per the "update mapping.md first" rule above.

**3. Retrospective GX-by-length on the existing full-scale P5/P7 production
baseline BAMs, complete (`existing_gx_analysis/`, 1-in-10 systematic sample of
primary records via `scripts/analyze_gx_by_length.py`).** Confirms the
mapping.md pattern at full production scale for both orientations: `GX=1%`
declines with length (P5 60.0%→51.6%, P7 62.3%→55.7%) while mean MAPQ *rises*
slightly with length (P5 52.1→54.8, P7 53.3→56.8) — longer reads place more
confidently but less gene-specifically, the same paradox mapping.md describes.
It also surfaces a length-dependent version of the P5/P7 orientation effect:
mean soft-clip fraction is comparable at short lengths (P5 3.1% vs. P7 3.3% at
≤250bp) but diverges sharply at long lengths (P5 15.6% vs. P7 5.2% at ≥401bp).
This is new evidence, not yet in mapping.md — feed it back there before the
slide cites it, since it directly strengthens the CB/UMI-position hypothesis
above (whichever segment is decoded second degrades more, and that shows up
disproportionately in the longest-read tail).

**4. Full-scale P5/P7 diagnostic evaluation, in progress, not yet available.**
`scripts/run_full_evaluation.sh` was launched 2026-07-30 06:38 against the full
raw P5 dataset (`CAT02/p5_out_7-7`) and was still running the P5 baseline arm's
minimap2/annotation stage as of this review — the P5 multimap arm and the
entire P7 arm have not started. This is the run that will actually test the
mapping.md question (`NH=1`-by-length ambiguity and downstream molecule/matrix
impact) on real long SBX reads rather than the compact Illumina-derived
fixture. **Do not put full-scale multimap NH=1-by-length numbers or a
baseline-vs-multimap downstream-impact comparison on any slide until this run
completes** — the only full-scale evidence available right now is the
retrospective GX-by-length pass on the existing baseline (secondary-alignment-
free) BAMs in item 3, which cannot speak to `NH=1` or multimap recovery at all.

## Sequencing depth context for M2

For the Multiome GEX comparison, preserve the distinction among these metrics:

- **Sequenced read pairs per filtered cell** is the relevant conventional depth
  metric. M2 Illumina input has 2,807,742,132 sequenced pairs and 612,462
  filtered cells: 4,584 pairs/cell.
- **Extracted tagged cDNA records per filtered cell** is the fair shared
  mapper-input depth for SBX: 1,811 for P5 and 1,607 for P7 in the current M2
  audit.
- **CB+UMI reads per filtered cell** is downstream depth/QC, not raw sequencer
  depth.
- Current SBX metrics do not retain a M2-specific count of raw reads before
  extraction; add a per-sample raw-read counter before parse/orientation/cDNA
  gating to calculate raw SBX reads per filtered cell.

**Open discrepancy to resolve before this goes on a slide:** the 612,462 figure
above is the sbxpress audit's own internal cell-calling heuristic applied to the
Illumina Cell Ranger BAM, not Cell Ranger's own `Estimated Number of Cells`
(32,357, from `metrics_summary.csv`). Cell Ranger's own denominator gives 86,774
mean reads/cell — its own reported value. Mixing the two conventions understates
depth/cell by ~19x and will read as a data-quality problem next to standard
10x depth targets (10-50k reads/cell) rather than the audit-vs-official-caller
mismatch it actually is. Decide which cell-calling convention the depth-per-cell
slide uses (official Cell Ranger call vs. the audit's own caller, consistently
for both Illumina and SBX arms) before publishing the number.

## Cross-platform comparison slides (CAT02 `result_comparison` M2 audit)

Separate from the mapping.md ambiguity question, `CAT02/result_comparison`
holds a completed three-way technical comparison (Illumina / SBX-P5 / SBX-P7,
each through both Cell Ranger and sbxpress) for the M2_GEX sample. This adds
slides that answer "how much sequencing does SBX need to match Illumina" and
"what does that imply for pricing" — orthogonal to the mapping-ambiguity work
below but part of the same deck. Proposed slides, in order:

1. **Study design** — one sample (M2_GEX) pulled from a 15-sample pooled
   NovaSeqX run, run through 3 technical preps × 2 pipelines. Note explicitly
   that the other 14 pooled samples were only demultiplexed, not carried
   through this comparison — this is a single-sample deep dive, not a
   population-level validation.
2. **Mapping efficiency by platform × pipeline** — unmapped %, gene-assigned %,
   from `sbxpress audit` run against both the Cell Ranger and minimap2 BAMs for
   each of the three preps. Headline: sbxpress out-recovers Cell Ranger on SBX
   reads (+6-7% more molecules) but underperforms Cell Ranger on standard
   Illumina reads (minimap2 leaves ~5x more reads unmapped) — the aligner is
   presently tuned toward long SBX reads, not short Illumina R2.
3. **Cell calling & sensitivity** — standard Cell Ranger `metrics_summary.csv`
   fields only (cells called, median genes/cell, median UMI/cell, sequencing
   saturation). Flag that Q30-in-barcode/RNA/UMI reads 0.0% for both SBX preps
   because `sbxpress extract --cellranger` writes a flat placeholder quality
   string into the synthetic FASTQ — not a real quality signal, omit or footnote
   it rather than let it read as an SBX quality problem.
4. **Read-depth equivalence** — fit a two-point Michaelis-Menten saturation
   curve per condition (`sampled_reads_with_cb_umi`/`sampled_molecules` vs.
   full-depth `reads_with_cb_umi`/`molecules` from each `sbxpress audit` JSON)
   and solve for the read depth SBX needs to match Illumina/Cell Ranger's full
   molecule recovery (285.7M molecules from 2.81B reads / 252.7 Gb of R2).
   Result: SBX-P7/minimap needs 39.8% of Illumina's read count but 108.0% of
   its Gb; SBX-P5/minimap needs 49.6% of reads but 124.7% of Gb. Fewer reactions,
   not fewer bases — state both units, since they point opposite directions.
5. **Cost-parity breakeven** — no external price data required, only the ratio
   from slide 4: breakeven price ratio = Illumina depth used / SBX depth needed.
   Confirmed: both platforms are priced per run/lane, so $/read is the right
   proxy (not $/Gb) — use the read-based framing as primary. Per-read billing:
   SBX can cost up to 2.51x (P7) / 2.01x (P5) Illumina's price per read and
   still break even. (Per-Gb framing — SBX priced at or below 92.6% (P7) /
   80.2% (P5) of Illumina's $/Gb — is now secondary context, not the headline,
   since output isn't actually billed that way.)
6. **CB/UMI position: a yield/fidelity tradeoff (P5 vs. P7)** — P5 (CB/UMI read
   first) parses substantially more raw signal into valid pairs (62.6% vs.
   46.9% `pairs_pct`; unparseable 13.9% vs. 33.8%) and calls barcodes more
   accurately (76.7% vs. 68.4% whitelist hit rate). P7 (cDNA read first) then
   wins on every downstream alignment metric on the reads that do survive
   (unmapped 1.37% vs. 6.78%; gene-assigned 60.26% vs. 58.25%; molecules
   recovered 296.7M vs. 281.6M). Working hypothesis: SBX signal quality decays
   over read position, so whichever segment is decoded second is degraded;
   barcode degradation is partially recoverable via whitelist correction,
   cDNA degradation is not, so protecting cDNA (P7) wins on final yield despite
   P7's lower parse success. This connects directly to the length-dependent
   clip-fraction divergence in the mapping-evaluation section below (P5's
   soft-clip fraction reaches 15.6% at ≥401bp vs. P7's 5.2% — same mechanism,
   now shown as a function of read length). Caveat to confirm before this slide
   ships: `sbx_5p_data`/`sbx_3p_data` naming also matches 10x's 5'-kit/3'-kit
   convention, a genuinely different capture chemistry, not just a read-direction
   flip of the same library — confirm P5/P7 here are the same molecule read in
   both orientations before calling this a pure orientation effect.
7. **Fragment-length sweet spot — open question, not yet answerable.** Because
   SBX bills per run/lane (slide 5), and run reads/lane is tunable by target
   molecule size (shorter molecule → more reads/run at fixed run cost), the
   real optimization target is *useful molecules per run*, i.e.
   `reads_per_run(length) × P(read becomes a counted molecule | length)`, not
   per-read quality alone. Current data only characterizes the second factor,
   and only partially:
   - Within the 250-650bp range already measured (`existing_gx_analysis`),
     gene-assignment specificity *improves* as length decreases (`GX=1%`: P5
     60.0%→51.6%, P7 62.3%→55.7% as length increases from ≤250 to ≥401bp), and
     soft-clip degradation drops sharply (P5 15.6%→3.1%). So within this
     range, shorter is strictly better on quality per read, not worse — the
     "reduced mapping rate" cost has to be a shorter-fragment effect than
     anything measured so far.
   - That analysis excludes unmapped reads by construction (`analyze_gx_by_length.py`
     scans primary mapped records only) and has no bin below 250bp, so it
     cannot see whether/where unmapped rate rises as fragments get shorter —
     which is exactly where a real "reduced mapping rate" cost would show up
     (too little unique sequence to place confidently).
   - `reads_per_run(length)` is a device/chemistry throughput property, not
     something in any file reviewed here (P5 vs. P7's different raw read
     totals, 19.39B vs. 23.77B, are confounded with whatever else differs
     between those two runs/preps and should not be used to infer a
     length-throughput curve).
   To answer this for real: (a) extend the length bins below 250bp (e.g.
   100/150/200/250) and, critically, compute unmapped-rate-by-length from the
   pre-alignment tagged-cDNA length distribution joined to mapping outcome —
   not from the post-hoc mapped-only BAM scan; (b) get the actual
   reads-per-run vs. target-molecule-length function from device throughput
   data, ideally from a controlled run set that varies target length at fixed
   device settings. Until both exist, this slide should state the shape of the
   problem and the two open curves, not a claimed optimum.

## End state and decision

End the presentation with this sequence:

1. Validate mapper/annotation, molecule, and cell decisions against the
   baseline.
2. Preserve all relevant audit evidence and explain every divergence.
3. Only then run controlled SBX-specific A/B experiments:
   - longer-read gene/transcript resolution;
   - gene-aware multimap policy;
   - read-length-aware trimming/filtering; and
   - collision/representative models informed by longer cDNA.

Success is a 10x-style baseline whose deviations are intentional, measurable,
and beneficial at the gene and molecule levels.

## Local source boundary

This brief is derived from the project discussion, current sbxpress source and
documentation, and the locally generated M2 audit summaries — including
`CAT02/result_comparison` (three-way Illumina/SBX-P5/SBX-P7 × Cell
Ranger/sbxpress comparison) and `CAT02/mapping_evaluation` (controlled
baseline-vs-multimap evaluation per `mapping.md`, partially complete as of
2026-07-30; see status section above). It is not an external performance claim
or a completed validation report. The fragment-length sweet-spot question is
explicitly unresolved — no throughput-vs-length or short-fragment mapping data
exists yet to answer it.
