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

## Current controls and limitations

The aligned-length analysis and soft-clip diagnostics reduce the likelihood that
adapter/poly(A/T) carryover explains the pattern. They do not identify the source
of each secondary alignment.

The currently available BAMs do not consistently provide a primary-secondary
alignment-score comparison suitable for the report. Therefore the present audit
cannot yet distinguish a nearly tied secondary hit from a weak alternative chain.
Also, aggregate length bins cannot separate a true length effect from changes in
the genes represented at different lengths.

## Discriminating follow-up analyses

These analyses should precede any SBX-specific mapping policy change.

1. **Within-gene stratification**
   - Compare `NH=1` by length within the same `GX`, transcript biotype,
     exon/intron class, and gene-length strata.
   - This tests whether length itself is associated with ambiguity after removing
     gene composition.
2. **Classify secondary hits**
   - For each `NH>1` query, record whether secondary hits are compatible with the
     primary gene, a paralog/pseudogene, a repeat/decoy region, or a different
     gene.
   - Summarize the distribution by read-length bin.
3. **Compare alignment evidence**
   - Record primary and secondary alignment scores when present, aligned query
     coverage, identity, clipping, and genomic overlap.
   - Plot score-gap and coverage-gap distributions by length and ambiguity class.
4. **Controlled remapping**
   - Re-map a stratified subset under the same index/preset with controlled
     secondary-hit settings, including the production fast mode
     (`--secondary=no`) and multimap mode (`-N 10`).
   - This separates minimap2 reporting behavior from sequence-intrinsic
     ambiguity.
5. **Gene-level impact**
   - Compare molecule and matrix counts under the alternatives, not only read
     `NH`. An increase in `NH>1` is acceptable if it does not change the
     gene-aware molecule decision or if `mm` makes that decision explicit.

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
