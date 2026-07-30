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

## Mapping observation to include

The M2 audit observed lower strict locus uniqueness (`NH=1`) at longer SBX cDNA
lengths, while gene assignment is more stable. This should **not** be described
as “long reads map worse.”

Current interpretation:

- Longer sequence can reveal an additional homologous genomic placement without
  creating another compatible gene.
- `NH` depends on minimap2 secondary-hit reporting (`-N`, index, preset, score
  and chaining thresholds), not only on sequence information.
- The effect persists in aligned-length bins and SBX has less heavy soft clipping
  than Illumina, reducing the likelihood that adapter or poly(A/T) carryover is
  the main explanation.
- Long and short read bins may differ in transcript class, intronic content,
  gene family composition, or molecule quality; this is not yet causal evidence
  for a read-length effect.

Refer readers to [mapping.md](mapping.md) for full evidence, tag definitions,
limitations, and discriminating analyses.

## Follow-up mapping analyses

Before proposing an SBX-specific mapping policy, plan:

1. Compare `NH=1` by length **within the same gene**, biotype, exon/intron
   class, and gene-length strata.
2. Classify secondary hits as same-gene, paralog/pseudogene, repeat/decoy, or
   different-gene alternatives.
3. Record primary/secondary alignment score gap, coverage, identity, clipping,
   and overlap where available.
4. Re-map a stratified subset under `--secondary=no` and `-N 10` to separate
   reporting policy from intrinsic sequence ambiguity.
5. Measure molecule and matrix consequences; `NH>1` alone is not a failure if
   the gene-aware molecule decision remains stable.

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
documentation, and the locally generated M2 audit summaries. It is not an
external performance claim or a completed validation report.
