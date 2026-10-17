# sbxpress: raw reads to 10x-style expression outputs

This manual covers the two supported GEX input paths:

- **Raw Roche SBX reads**: demultiplex, extract tagged cDNA, align, annotate,
  call molecules, make raw counts, then call cells.
- **Paired Illumina R1/R2 reads**: retain the existing barcode/UMI read, convert
  the cDNA read to the same tagged handoff, then run the shared downstream path.

The final layout is intentionally familiar:

```text
outs/
  raw_feature_bc_matrix/                 # MEX: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
  raw_feature_bc_matrix.h5
  filtered_feature_bc_matrix/            # called-cell MEX
  filtered_feature_bc_matrix.h5
  metrics_summary.json
```

`sbxpress` is responsible for these outputs; the input paths differ only before
alignment.

## Per-barcode read metrics

With `--barcode-metrics`, each read-processing stage (`illumina`, `extract`,
`annotate`, `molecules`) also writes a gzip TSV beside its metrics JSON; the
pipeline scripts accept the same flag and pass it to every stage. It is off by
default because it costs run time and memory (see below). Replace
`.json` in the metrics path with `.barcodes.tsv.gz` to find it. Rows are sorted
by `barcode`, then `sample`. The first columns are `barcode` and `sample`;
all other columns are integer counts. Summing each count column over all rows
reproduces the corresponding run-wide JSON counter.

| Stage | Count columns |
|---|---|
| `illumina` | `reads`, `barcode_exact`, `barcode_corrected`, `emitted`, `failed_merged_too_short`, `failed_no_cdna`, `failed_short_r1`, `failed_r2_polya_short`, `failed_cb_invalid`, `failed_cb_unmatched`, `low_complexity`, `failed_low_complexity`, `failed_r2_tso_short` |
| `extract` | `reads`, `cb_on_whitelist`, `cb_within1`, `emitted`, `class_full`, `class_partial`, `class_unparseable`, `pb_located`, `pb_hit`, `low_complexity` |
| `annotate` | `records`, `primary`, `unmapped`, `exonic`, `intronic`, `intergenic`, `antisense`, `gene_assigned`, `unique_gene`, `multi_gene`, `rescued`, `end_overhang` |
| `molecules` | `eligible_reads`, `raw_cell_gene_umis`, `molecules_counted`, `collision_groups`, `representatives_tagged`, `rescued_reads`, `low_mapq_reads`, `multimapper_reads` |

`__no_barcode__` groups reads without a located or usable barcode;
`__unresolved_barcode__` groups reads with a raw barcode that did not resolve
to a cell (including reads rejected before correction). Neither row is a
cell. For an Illumina GEX run, join the `illumina`, `annotate`, and `molecules`
tables by `barcode` to inspect one cell's read funnel. An SBX GEX run uses
`extract` in place of `illumina`. Flex molecules use `PB` in `sample` to keep
shared cell barcodes from different probe samples separate; `extract` uses
the outer sample-sheet assignment, so these two sample fields have different
meanings. Annotation counts BAM records, including extra alignments when
secondary hits are retained; its historical `primary` counter includes mapped
secondary records. Molecule eligibility counts primary reads only.

Cost of `--barcode-metrics`, measured against the same build without it on
subsets of real data (3 runs each): `illumina` and `molecules` within 1%,
`annotate` +7% wall time, SBX `extract` +12% wall time. Memory is unchanged
except in `extract`, where the table holds one row per observed
barcode/sample pair: +375 MB on one CAT06 Flex chunk (~43 M reads), growing
with the number of distinct barcodes in the run. Without the flag no stage
does any per-barcode work.

## Confident mapping and multimapper rescue

`sbxpress molecules` counts a read only if it is confidently mapped, the rule
Cell Ranger uses: MAPQ of at least `--min-mapq` (default 60, minimap2's unique
placement; STAR's 255 also passes), or rescued by `sbxpress annotate`. Without
it, reads that fit many places in the genome are counted wherever minimap2
happened to put them; on DZEC05 poly-G reads placed in LINC00486 gave that gene
22% of all reads, and the rule takes it to 0.01%.

`sbxpress annotate` marks a rescue with `mm:i:1` and keeps minimap2's MAPQ. A
primary read below `--confident-mapq` (default 60) is rescued when its genes
come from splice-compatible, sense, exonic transcripts of exactly one gene.
Intronic placements are never rescued, even with `--include-introns`, because
they carry no transcript evidence. With `--multimap-tags`, the rule spans the
read's reported hits as in Cell Ranger: no rescue if any hit is confidently
mapped; otherwise the exonic-compatible hits must name one gene and the
primary must be one of them. `NH` (reported hits) is written only with
`--multimap-tags`.

The molecules metrics count each outcome: `rescued_records` (counted through
`mm:i:1`), and, for dropped reads, `low_mapq_records` or `multimapper_records`
(`NH` > 1). A confidently mapped read is counted even when `NH` > 1, as in
Cell Ranger. Dropped reads are not in `eligible_records` and get no molecule
tags. `annotate` reports `rescued`.

On a normal library (DZEC04, four chromosomes sampled) 78% of single-gene
reads have MAPQ 60, 19% are rescued and 3% are dropped. The dropped reads are
intronic and mostly inside repeats. `molecules --min-mapq 0` and
`annotate --confident-mapq 0` restore the pre-0.6.2 counts.

### Read-end exon overhang

`sbxpress annotate` gives a read a gene only if its aligned blocks fit a
transcript's exons with exact splice junctions. minimap2 often extends a
read's alignment 1-3 bases past an exon end instead of splicing such a short
anchor, and the read is then exonic (`RE:A:E`) with no gene. On DZEC04 this
was 7% of records; 55% of them overhung by 1-3 bases, almost all at a read
end, and those bases also matched the next exon in 98% of cases.
`--max-end-overhang` (default 3) lets the first and last aligned block extend
that many bases past their exon. The overhang is treated as soft-clipped, and
the `TX` position is clamped to the exon. Internal junctions must still match
exactly. The `end_overhang` counter reports records that got genes this way;
`--max-end-overhang 0` restores exact ends.

## Low-complexity reads and the artifact report

A library or sequencing fault can fill reads with one base. On DZEC05,
28-56% of read 2s were dead (a two-colour sequencer calls "no signal" as G,
at about Q25, so quality filters miss it), and the poly-G reads that mapped
gave LINC00486 a fifth of all molecules. Such reads are a library fault,
not biology, so conversion drops them by default, and two reports make the
fault visible:

- **At conversion** (`sbxpress illumina`, `sbxpress extract`), each
  barcode-passing read's cDNA is checked: it is low-complexity if one base
  makes up at least `--low-complexity-base-fraction` of it (default 0.8) or
  its longest single-base run covers at least `--low-complexity-run-fraction`
  (default 0.5). Either test is off at 0. Flagged reads fail with reason
  `low_complexity` (written to `--failed-output` when set).
  `--drop-low-complexity false` keeps them instead, with `YL:i:1` in the
  FASTQ comment (minimap2 `-y` and `sbxpress-probeset-align` carry it into
  the BAM), so `molecules` can name the genes they land on. Either way the
  metrics JSON reports `low_complexity` and `low_complexity_by_base` (by the
  read's most frequent base), and the per-barcode table has `low_complexity`
  (and, for `illumina`, `failed_low_complexity`). The stage warns when more
  than 5% of reads are poly-G or poly-N (the dead-read signal).
  `illumina` trims poly-A read-through from unmerged R2 by default
  (`--trim-r2-polya`; 3' GEX and tag-seq both read into the tail); with
  `--trim-r2-polya false` that read-through is expected to look
  low-complexity, so the warning then counts all flagged reads (DZEC04
  flagged 4.6%, 90% poly-A, untrimmed; 0.48%, 85% poly-G, trimmed). The
  check is one pass over bases already in memory.
- **In `sbxpress molecules`**, the metrics JSON has an `artifacts` section:
  gene-assigned primary reads, how many carry `YL:i:1`, how many are below
  MAPQ 60 without `mm:i:1` (counted whatever `--min-mapq` is), and up to 25
  genes with at least 100 reads of which at least half are flagged. A listed
  gene holding at least 1% of gene-assigned reads is printed as a warning.
  In a normal library a few intronic-repeat genes are listed, each a small
  share of reads; DZEC05's LINC00486 would be listed with a fifth of them.
  With the default drop, `YL` reads never reach `molecules`, so the report
  then shows only the low-MAPQ side.

## 1. Inputs shared by both paths

Prepare one consistent reference set:

- a chemistry-matched 10x barcode whitelist (`.txt` or `.txt.gz`);
- a minimap2 splice index (`.mmi`);
- the matching annotation GTF (`.gtf` or `.gtf.gz`); and
- a junction BED used with the splice index.

The first annotation run writes a canonical `genes.mtx` feature manifest. Reuse
that exact file for every sample using the same reference; it fixes feature order
and preserves comparable raw and filtered matrices across an experiment.

For a 3′ GEX library, use `--three-prime` during annotation. This affects the
tie-breaking and terminal-support tags (`ZG`, `ZT`, `XP`); it does not discard
non-terminal alignments. The SBX P5/P7 orientation changes read extraction, not
the biological 3′ assay preference.

## 2. Raw SBX reads: complete, restartable pipeline

Both pipeline drivers record input paths, sizes, modification times, and analysis
settings for safe restarts. Changed inputs/settings or pre-existing outputs with
no restart manifest require a new output directory. Matrix and cell-calling
stages use temporary directories and completion markers, so interrupted writes
are rebuilt. Keep each driver's shared helper files beside the driver.


Use the repository workflow for production SBX processing. You can use the unified pipeline runner `scripts/sbxpress-pipeline.sh` which supports both `--input-mode {sbx|illumina}` and `--analysis-mode {gex|flex|tagseq}`:

```bash
REPO=/path/to/sbxpress
OUT=/analysis/run_sbx
REF=/analysis/resources

"$REPO/scripts/sbxpress-pipeline.sh" \
  --input-mode sbx \
  --analysis-mode gex \
  --fastqs /data/device_0/out /data/device_1/out /data/device_2/out /data/device_3/out \
  --bclconvert-samplesheet /data/samples.csv \
  --sample M2_GEX \
  --whitelist "$REF/737K-arc-v1.txt.gz" \
  --mmi "$REF/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi" \
  --gtf "$REF/gencode.v44.primary_assembly.basic.annotation.gtf.gz" \
  --junctions "$REF/gencode.v44.primary_assembly.junctions.bed" \
  --out-dir "$OUT" \
  --extract-threads 128 \
  --map-threads 32 \
  --annotate-threads 8 \
  --sort-threads 8 \
  --molecule-threads 8 \
  --split 16 \
  --audit
```

The script uses `sbxpress extract --minimap2`, which writes one tagged cDNA R2
FASTQ per shard under `OUT/extract/`. The read comment carries `CR`, `CB`, `UR`,
`UB`, `CY`, `UY`, and SBX provenance tags; `sbxpress minimap -y` transfers them
to the alignment records.

Use `--multimap-tags` when the downstream comparison needs reported-hit tags:

```bash
  --multimap-tags
```

This retains the primary plus up to ten secondary alignments (`NH` therefore
ranges through the top-coded `11+` audit bucket) and emits `NH`, `HI`, and the
single-gene multimap rescue tag `mm` (see "Confident mapping and multimapper
rescue"). It uses query-grouped streaming annotation,
so it trades the default unordered annotation parallelism for complete per-read
hit accounting. Omit it for the fastest production count path.

Key outputs for `M2_GEX` are:

```text
OUT/extract/metrics_summary.json
OUT/reference/genes.mtx
OUT/M2_GEX/annotation/M2_GEX.annotated.coordinate.bam
OUT/M2_GEX/molecules/M2_GEX.molecules.coordinate.bam
OUT/M2_GEX/outs/raw_feature_bc_matrix.h5
OUT/M2_GEX/outs/filtered_feature_bc_matrix.h5
OUT/M2_GEX/audit/M2_GEX.audit.json                 # when --audit is used
```

## 3. Raw SBX reads: paired R1/R2 compatibility output only

When only paired short-read-style FASTQs are required, use extraction directly
and omit `--minimap2`:

```bash
sbxpress extract \
  --fastqs /data/device_0/out /data/device_1/out /data/device_2/out /data/device_3/out \
  --bclconvert-samplesheet /data/samples.csv \
  --whitelist /references/737K-arc-v1.txt.gz \
  --out-dir /analysis/sbx_fastqs \
  --threads 128 --split 16 --lenient \
  --cellranger
```

The `--cellranger` flag is the paired-FASTQ compatibility preset. It writes
per-sample R1 (cell barcode + UMI) and R2 (cDNA) files, and caps R2 to 649 bases.
It is intentionally mutually exclusive with `--minimap2`.

## 4. Illumina R1/R2 reads: shared minimap2 path

Illumina input is already paired R1/R2. `sbxpress illumina` reads the barcode and
UMI from R1, checks/corrects the barcode against the whitelist, and emits a tagged
cDNA-only FASTQ from R2. Repeat `--r1` and `--r2` once per lane in the same order;
the lanes are converted concurrently with `--threads`.

```bash
OUT=/analysis/illumina_minimap/M2_GEX
REF=/analysis/resources
mkdir -p "$OUT" "$OUT/annotation" "$OUT/molecules" "$OUT/outs"

sbxpress illumina \
  --r1 /data/M2_GEX_S8_L003_R1_001.fastq.gz \
  --r2 /data/M2_GEX_S8_L003_R2_001.fastq.gz \
  --r1 /data/M2_GEX_S8_L004_R1_001.fastq.gz \
  --r2 /data/M2_GEX_S8_L004_R2_001.fastq.gz \
  --whitelist "$REF/737K-arc-v1.txt.gz" \
  --output "$OUT/M2_GEX.tagged.R2.fastq.gz" \
  --threads 6
```

Run the same mapping-to-cells path. The standard fast mode suppresses secondary
alignments and allows unordered annotation. `-x splice -m 20 -s 30` is the
pipeline's mapping preset: minimap2 2.28 runs `-x splice:sr` as `-x splice`
anyway, and the lower chain and alignment score minimums keep short reads that
would otherwise go unmapped (CHANGELOG 0.6.2):

```bash
sbxpress minimap -a -x splice -m 20 -s 30 \
  --junc-bed "$REF/gencode.v44.primary_assembly.junctions.bed" \
  --secondary=no -y -t 32 \
  "$REF/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi" \
  "$OUT/M2_GEX.tagged.R2.fastq.gz" \
| sbxpress annotate - \
    --gtf "$REF/gencode.v44.primary_assembly.basic.annotation.gtf.gz" \
    --output - --metrics "$OUT/annotation/M2_GEX.annotation.metrics.json" \
    --genes-mtx "$REF/genes.mtx" --three-prime --unordered --threads 8 \
| samtools sort -@ 8 -m 1G -T "$OUT/M2_GEX.annotated" \
    -o "$OUT/annotation/M2_GEX.annotated.coordinate.bam" -

samtools index -@ 8 "$OUT/annotation/M2_GEX.annotated.coordinate.bam"

sbxpress molecules "$OUT/annotation/M2_GEX.annotated.coordinate.bam" \
  --output "$OUT/molecules/M2_GEX.molecules.coordinate.bam" --threads 8
samtools index -@ 8 "$OUT/molecules/M2_GEX.molecules.coordinate.bam"

sbxpress matrix "$OUT/molecules/M2_GEX.molecules.coordinate.bam" \
  --output "$OUT/outs" --genome GRCh38 \
  --features "$REF/genes.mtx" --cell-barcodes "$REF/737K-arc-v1.txt.gz"

sbxpress cells "$OUT/outs/raw_feature_bc_matrix.h5" --output "$OUT/outs"
```

For the multimap-tag mode, replace `--secondary=no` with `-N 10` and replace
`--unordered` with `--multimap-tags`. Do not use both annotation modes together.

## 5. Audit and comparison

Audit operates on the molecule BAM, after UMI correction, collision resolution,
and `xf` representative tagging. Supply filtered barcodes to measure the share
that reaches the filtered expression matrix:

```bash
sbxpress audit "$OUT/molecules/M2_GEX.molecules.coordinate.bam" \
  --label ILLUMINA_MINIMAP \
  --cells "$OUT/outs/filtered_feature_bc_matrix/barcodes.tsv.gz" \
  --output "$OUT/audit/M2_GEX.audit.json" \
  --molecules "$OUT/audit/M2_GEX.molecules.tsv.gz" \
  --molecule-sample 64 --threads 8
```

Use the same `--molecule-sample` value for every compared BAM. The JSON is an
all-record summary; the TSV is a deterministic hash sample used for joinable
cross-platform molecule comparisons.

## 6. 10x Flex (Fixed RNA Profiling) pipeline

10x Flex assays use fixed synthetic probe pairs rather than polyA-primed cDNA. Reads are assigned directly to probe sequences in the probe-set reference (`.csv`) rather than aligned to the genome.

The simplest route is the unified runner (single sample, or add
`--probe-sample-sheet config.csv` for multiplexed samples):

```bash
"$REPO/scripts/sbxpress-pipeline.sh" --analysis-mode flex \
  --fastqs /data/device_*/out --bclconvert-samplesheet /data/samples.csv \
  --sample TSeq100_L1_GEX \
  --whitelist "$REF/737K-flex-v2.txt.gz" \
  --probe-barcode-whitelist "$REF/flex-v2-384.txt" \
  --probe-set "$REF/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv" \
  --probe-sample-sheet "$REF/config.csv" \
  --out-dir "$OUT"
```

### Workflow overview

1. **Extraction**: `sbxpress extract --flex --probe-barcode-whitelist` locates the probe barcode and writes it as `PR`/`PB` alongside `CB`/`UB`.
2. **Probe-set alignment**: `sbxpress probeset-align` assigns each extracted R2 chunk to the probe set, emitting a BAM on synthetic per-gene contigs with `GX`/`GN` tags (`--threads` parallelises matching; output order follows the FASTQ).
3. **Molecules**: `sbxpress molecules --flex` requires a non-empty `PB` so identical GEM barcodes from different samples are never merged.
4. **Counting**: a single sample goes straight to `matrix` + `cells`; multiplexed samples first run `sbxpress probe-demux` on the molecule BAM, then `scripts/run_flex_multi.sh` counts each sample.

### Example command sequence

```bash
# 1. Extract (probe-barcode aware)
sbxpress extract --fastqs /data/device_*/out \
  --bclconvert-samplesheet /data/samples.csv \
  --whitelist "$REF/737K-flex-v2.txt.gz" \
  --probe-barcode-whitelist "$REF/flex-v2-384.txt" --flex \
  --minimap2 --lenient --out-dir "$OUT/extract" --threads 128 --split 16

# 2. Probe-set alignment per extracted chunk, then one sorted BAM
mkdir -p "$OUT/aligned"
for fq in "$OUT"/extract/TSeq100_L1_GEX_S*_R2_*.fastq.gz; do
  sbxpress probeset-align "$fq" \
    --probe-set "$REF/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv" \
    --genes-out "$OUT/reference/probe_set_genes.tsv" --threads 16 \
    --output "$OUT/aligned/$(basename "$fq" .fastq.gz).bam"
done
samtools cat -o - "$OUT"/aligned/*.bam | samtools sort -@ 16 -o "$OUT/L1_GEX.probeset.bam" -
samtools index "$OUT/L1_GEX.probeset.bam"

# 3. Molecules (PB required)
sbxpress molecules "$OUT/L1_GEX.probeset.bam" --flex --output "$OUT/L1_GEX.molecules.bam"
samtools index "$OUT/L1_GEX.molecules.bam"

# 4a. Single sample
sbxpress matrix "$OUT/L1_GEX.molecules.bam" --output "$OUT/outs" \
  --features "$OUT/reference/probe_set_genes.tsv" --cell-barcodes "$REF/737K-flex-v2.txt.gz"
sbxpress cells "$OUT/outs/raw_feature_bc_matrix.h5" --output "$OUT/outs"

# 4b. Multiplexed samples (config.csv [samples] section)
sbxpress probe-demux "$OUT/L1_GEX.molecules.bam" \
  --probe-sample-sheet "$REF/config.csv" --output "$OUT/demux"
scripts/run_flex_multi.sh "$OUT/L1_GEX.molecules.bam" "$OUT/demux" "$OUT/results"
```

Illumina Flex input (`sbxpress illumina`) carries no probe barcode, so it supports
single-sample counting only. For the full multi-library CAT06 run (GEX + CRISPR),
see `scripts/run_full_cat06.sh`.

## 7. CAT Tagseq & Targeted Expression Mode

Targeted assays, high-depth bulk/pooled samples, and CAT Tagseq datasets require specialized molecule deduplication and cell-calling rules due to oversubscribed UMI space and pre-determined cell lists. (CRISPR guide capture at single-cell depth does *not* use `--tagseq`; see `scripts/run_flex_multi.sh`.)

### UMI collision policies (`sbxpress molecules`)

In high-depth bulk or targeted TagSeq samples (e.g., 2–25M reads/sample into a 12bp UMI space), distinct molecules often share a UMI by chance — DZEC03 measured ~2.1x oversubscription, with 81% of `(cell, UMI)` pairs claimed by more than one, genomically unrelated, gene. `--tagseq` therefore (a) skips 1-mismatch UMI correction, since two real UMIs 1bp apart are common at this depth, and (b) does not discard the minor gene in a collision:

| Flags | Collision outcome (`XW`) | Matrix contribution |
|---|---|---|
| default (single-cell) | `W` winner / `L` loser, `T` on exact tie | 1 for the winner, 0 otherwise |
| `--tagseq` | `S` for every gene | 1 for every gene (the collision is treated as distinct molecules) |
| `--fractional-collisions` (with or without `--tagseq`) | `F` for every gene | `w_i = r_i / Σr` (read support); weights of one collision sum to 1 |

Fractional splitting counts a collision as **one** molecule. That is right if collisions are mostly one molecule seen on two genes, and an undercount — largest for the most expressed genes — if they are mostly distinct molecules that drew the same UMI, which is what DZEC03 shows. It is therefore opt-in (it was the `--tagseq` default in 0.5.6–0.5.9). Fractional molecules carry the weight as `ZW:f`; the representative read is flagged `xf:25` like any counted molecule; `ZG` is never touched. In the metrics JSON a fractional collision adds 1 to `molecules_counted`.

```bash
sbxpress molecules "$OUT/annotation/sample.annotated.bam" \
  --tagseq \
  --output "$OUT/molecules/sample.molecules.bam"
```

Fractional matrices retain small positive weights and MEX precision; HDF5 stores
fractional values as `f32`. Automatic EmptyDrops requires integer counts and
rejects fractional input. For fractional matrices, select known samples with
`--tagseq-barcodes` or use `--force-cells` for an explicit ranked selection.
Fractional saturation counts a collision group once, at its first supporting
read; it uses full-depth molecule decisions rather than refitting at each depth.

### Position-aware UMIs (`--position-aware-umis`)

Every policy above still merges two distinct molecules of the **same** gene that drew the same UMI — invisible, and frequent for highly expressed genes at tag-seq depth. `--position-aware-umis` keys molecules on `(cell, gene, UMI, read 5′ start)` (strand-aware, unclipped), the usual deduplication key for bulk 3′ tag-seq: PCR duplicates share a start, chance UMI matches almost never do. Because 1-mismatch neighbours are then only compared at the same start, where they are almost always sequencing errors, UMI correction runs in this mode even with `--tagseq`. It works with any collision policy and with `--shards`/`--evidence-workers`; UMIs must be ≤16 bases.

```bash
sbxpress molecules "$OUT/annotation/sample.annotated.bam" \
  --tagseq --position-aware-umis \
  --output "$OUT/molecules/sample.molecules.bam"
```

Caveat: duplicates whose start differs (an indel or soft-clip difference at the read's 5′ end) are counted separately, so this mode can overcount slightly. Validate on your data before switching.

### SBX UMI rules for CRISPR guide capture (`--umi-indel-correction`, `--fold-singletons-min-top`)

Default UMI correction merges a UMI into a 1-mismatch neighbour with at least 2n−1 reads. SBX errors are mostly 1-base indels; inside the 12-bp UMI an indel shifts every later base, so its copy is several mismatches from the parent and never merges. In a deeply sequenced library this adds up. In CAT06 SBX CRISPR (~450 reads per guide molecule), sbxpress counted 6.1× Cell Ranger's molecules on depth-matched Illumina, which pushed cells with 1–2 true guide UMIs over the 3-UMI call floor (RZ211 CD45: 80% of cells called vs Cell Ranger's 60%; Neg: 4.5% vs 0.9%).

- `--umi-indel-correction` also merges one-indel neighbours under the same 2n−1 rule (CAT06: 6.1× → 1.9×).
- `--fold-singletons-min-top K` then folds every 1-read molecule into its (cell, sample, gene) family's top molecule when that top has at least `K` reads. The remaining excess was low-quality reads with several UMI errors (50% had two or more sub-Q20 UMI bases, vs 2.3% of real-molecule reads). With both rules and `K = 10`, CAT06 matched Cell Ranger: 0.99× the molecules, RZ211 59.4% called, Neg 0.80%. Any `K` from 5 to 50 gave results within about a point. Folded molecules are reported as `umi_singletons_folded` in the metrics JSON.

`scripts/run_flex_multi.sh` turns both on for the CRISPR stream only. Both are off by default and are not meant for GEX: most GEX molecules have few reads, so a 1-read UMI there is usually a real molecule. Their effect on GEX has not been measured.

```bash
sbxpress molecules "$OUT/crispr/guides.bam" --flex \
  --umi-indel-correction --fold-singletons-min-top 10 \
  --output "$OUT/crispr/crispr_molecules.bam"
```

### Choosing a mode: compare against a reference

`scripts/compare_tagseq_modes.sh` runs one annotated BAM through each mode and compares every matrix with an independent reference — typically the same library sequenced on Illumina and counted separately:

```bash
scripts/compare_tagseq_modes.sh \
  --bam "$OUT/DZEC03/annotation/DZEC03.annotated.coordinate.bam" \
  --reference /path/to/illumina/outs/raw_feature_bc_matrix.h5 \
  --features "$OUT/reference/genes.mtx" --out-dir "$OUT/mode_comparison" \
  --threads 16 --molecule-args "--shards 8"
```

It writes `summary.tsv` (per mode: per-barcode correlation of library-normalised counts, median per-gene |log2 ratio|, and the slope of log2 ratio against expression) and `per_gene.tsv`. Depth differs between platforms, so judge modes on the normalised metrics: collision undercounting shows up as a **negative slope** (the most expressed genes lose the most); a mode without collision bias is flat. Needs `numpy` (`h5py` for an `.h5` reference).

### Real & Float32 Expression Matrices (`sbxpress matrix`)

When `sbxpress matrix` (`cr-bam-to-mtx`) encounters `ZW:f` floating-point weight tags in the molecules BAM, it automatically materializes fractional count matrices:
- **Matrix Market MEX**: Emits real-valued coordinate headers (`%%MatrixMarket matrix coordinate real general`).
- **HDF5 Matrix**: Emits 32-bit floating-point datasets (`f32`) for `data` in `raw_feature_bc_matrix.h5`.

Standard integer runs without fractional weights continue emitting exact integer MTX and 32-bit unsigned integer (`u32`) HDF5 datasets.

`sbxpress cells` reads either type. Fractional values are kept in the filtered matrix; OrdMag and EmptyDrops, which need whole-UMI totals, use per-barcode and per-feature sums rounded to the nearest integer.

### Whitelist Cell Calling (`sbxpress cells --tagseq-barcodes`)

When exact cell barcode identities are known beforehand (e.g. called by a coupled 3′ GEX run or specified via an experimental whitelist), statistical knee-point cell calling (OrdMag / EmptyDrops) can be bypassed.

Pass `--tagseq-barcodes` to `sbxpress cells` to restrict matrix filtering strictly to the provided barcode file:

```bash
sbxpress cells "$OUT/outs/raw_feature_bc_matrix.h5" \
  --tagseq-barcodes "$REF/known_gex_cells.txt" \
  --output "$OUT/outs"
```
