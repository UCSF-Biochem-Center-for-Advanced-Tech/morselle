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

Use the repository workflow for production SBX processing. It runs extraction,
minimap2, annotation, one coordinate sort, molecules, raw matrix creation, and
cell calling. It reuses a completed stage on rerun; delete that stage's output if
you intentionally need to regenerate it.

```bash
REPO=/path/to/sbxpress
OUT=/analysis/run_sbx
REF=/analysis/resources

"$REPO/scripts/sbxpress-gex-pipeline.sh" \
  --fastqs /data/device_0/out /data/device_1/out /data/device_2/out /data/device_3/out \
  --sample-sheet /data/samples.csv \
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
single-gene multimap rescue tag `mm`. It uses query-grouped streaming annotation,
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
  --sample-sheet /data/samples.csv \
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
alignments and allows unordered annotation:

```bash
sbxpress minimap -a -x splice:sr \
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
