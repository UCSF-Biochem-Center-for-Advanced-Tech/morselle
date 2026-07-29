# Preparing a sbxpress minimap reference for SBX 3′ GEX counting

This manual prepares a human reference for `sbxpress extract --minimap2` output: accurate,
roughly 50–649 bp cDNA reads derived from a 3′ poly(dT)-primed GEX assay. It is
intended for gene-level molecule counting, not isoform discovery.

The recommended strategy is **spliced alignment to the genome**, using the GENCODE
annotation to prefer known splice junctions. Do not align only to a transcriptome:
genomic splice-aware alignments give the coordinate, intron/exon context, and
multimapping information needed for robust feature assignment and 3′-end QC.

## Inputs

The local reference directory used in these examples is:

```text
/Users/biomsettles/projects/src/ucsf_projects/roche_abrf/working_with_raw_data/reference
```

It contains the matched GENCODE/GRCh38 inputs:

```text
GRCh38.primary_assembly.genome.fa.gz
gencode.v50.primary_assembly.annotation.gtf.gz
```

Keep the FASTA and GTF from the same release/build. Their contig names must also
match (for example, both use `chr1`, not one `chr1` and the other `1`).

## 1. Expand the inputs

```bash
ref=/Users/biomsettles/projects/src/ucsf_projects/roche_abrf/working_with_raw_data/reference
cd "$ref"

gunzip -k GRCh38.primary_assembly.genome.fa.gz
gunzip -k gencode.v50.primary_assembly.annotation.gtf.gz
```

The `-k` keeps the downloaded compressed files. The uncompressed FASTA is needed
while building the index; sbxpress minimap must retain target sequence in the index because
we will request SAM alignments with CIGAR strings.

## 2. Create the sbxpress minimap index

Build an index with the seed parameters used by the `splice:sr` preset:

```bash
sbxpress minimap -k15 -w5 \
  -d GRCh38.primary_assembly.gencode-v50.splice-sr.mmi \
  GRCh38.primary_assembly.genome.fa
```

Index seed settings are stored in `.mmi`; mapping-time `-k`/`-w` values cannot
override them. Therefore, build this separate index rather than reusing an index
made with the default window size.

Do **not** use `--idx-no-seq`: that makes a smaller index but prevents base-level
SAM/CIGAR alignment (`-a`), which the counting workflow requires.

## 3. Make an annotated-junction BED file

Minimap2 accepts GENCODE splice junctions as BED12 via `--junc-bed`. This is
recommended: it favors known junctions but does not prohibit genuine novel ones.

The sbxpress minimap distribution supplies `paftools.js`; use it if available:

```bash
paftools.js gff2bed gencode.v50.primary_assembly.annotation.gtf \
  > gencode.v50.primary_assembly.junctions.bed
```

If `paftools.js` is not on `PATH`, obtain it from the same sbxpress minimap release as the
installed executable, or use an equivalent GTF/GFF-to-BED12 converter. Confirm the
result has 12 tab-delimited columns and contig names matching the FASTA:

```bash
awk -F '\t' 'NF != 12 { print "non-BED12 line:", NR; exit 1 }' \
  gencode.v50.primary_assembly.junctions.bed
```

## 4. Generate minimap2-ready SBX FASTQ

Use the R2-only mode. It emits reverse-complemented cDNA in standard R2 files and
puts barcode/UMI fields in the FASTQ comment:

```bash
sbxpress extract \
  --fastqs /data/sbx/run \
  --sample-sheet /data/sample_sheet.csv \
  --whitelist /path/to/737K-arc-v1.txt.gz \
  --out-dir /data/sbxpress_minimap2 \
  --minimap2
```

Each R2 header carries:

```text
CR:Z:<raw-cell-barcode>
CB:Z:<corrected-cell-barcode>-1
UR:Z:<raw-UMI>
UB:Z:<provisional-UMI>
CY:Z:<raw-cell-barcode-qualities>
UY:Z:<UMI-qualities>
```

`UB` is not a final molecule UMI: it is the extracted, pre-gene-assignment value.
It may change during gene-aware UMI correction.

## 5. Align reads

Use the accurate short-RNA preset, not ONT/PacBio presets. Although SBX cDNA reads
can be longer than conventional Illumina reads, they are still short, high-accuracy
fragments—not noisy multi-kilobase long reads.

```bash
ref=/Users/biomsettles/projects/src/ucsf_projects/roche_abrf/working_with_raw_data/reference
fastq=/data/sbxpress_minimap2/sample_S1_L001_R2_001.fastq.gz

sbxpress minimap -a -x splice:sr \
  --junc-bed "$ref/gencode.v50.primary_assembly.junctions.bed" \
  --secondary=no \
  -y \
  -t 16 \
  "$ref/GRCh38.primary_assembly.gencode-v50.splice-sr.mmi" \
  "$fastq" \
| samtools sort -@ 8 -o sample.coordinate.bam

samtools index sample.coordinate.bam
```

Important options:

- `-x splice:sr` selects spliced alignment and settings for accurate short RNA
  reads. It remains suitable for several-hundred-base SBX cDNA fragments.
- `--junc-bed` favors GENCODE splice junctions.
- `--secondary=no` prevents secondary alignments from being emitted as additional
  countable records. Retain MAPQ and apply an explicit multimapping policy later.
- `-y` is required: it copies FASTQ comments to SAM/BAM, preserving the CB/UMI tags.
- `samtools sort` without `-t` makes the coordinate-sorted BAM used for inspection,
  indexing, and QC.

Use each sample's R2 shards as inputs. For split output, pass all matching R2 files
to minimap2, or concatenate/decompress them as a stream; do not include R1 because
`--minimap2` intentionally does not emit it.

## 6. 3′-tag handling and feature assignment

Do not force the mapper to terminate at an annotated transcript end. The poly(dT)
primer identifies a 3′-captured molecule, but real data include alternate
polyadenylation, incomplete cDNA, internal priming, and annotation differences.

Instead, after alignment:

1. Assign each primary alignment to an exon/gene using the same GENCODE GTF.
2. Apply mapping-quality, alignment-length, and gene-overlap rules.
3. Measure the read end adjacent to the poly(A) capture site against annotated
   transcript 3′ ends as a QC/optional filter, with strand awareness.
4. Retain the gene assignment needed to constrain UMI correction.

Create or verify the canonical feature manifest while annotating. It is derived
from the same filtered GTF gene set used for `GX`, so reference construction and
feature assignment cannot drift apart:

```bash
sbxpress minimap ... \
| sbxpress annotate - \
    --gtf "$ref/gencode.v50.primary_assembly.annotation.gtf.gz" \
    --genes-mtx "$ref/genes.mtx" \
    --output - --metrics sample.annotation.metrics.json --three-prime \
| samtools sort -@ 8 -o sample.annotated.coordinate.bam -
```

On the first sample, `--genes-mtx` creates a three-column manifest containing
gene ID, gene name, and `Gene Expression`. Later samples verify that the existing
file is byte-identical to their filtered annotation and fail on a mismatch.
Always pass it to matrix generation:

```bash
sbxpress matrix sample.molecules.coordinate.bam \
  --features "$ref/genes.mtx" \
  --cell-barcodes /path/to/737K-arc-v1.txt.gz \
  --output sample.outs --genome GRCh38
```

With a feature manifest, an unlisted counted `GX` is an error. This prevents a
sample-specific gene from being silently appended and changing matrix row order.
`--cell-barcodes` similarly expands the 10x reference whitelist into canonical
`-1` suffixed matrix columns, retaining cells with zero molecules. Use the same
gene manifest and barcode whitelist for every comparable experiment.

The correct genomic end to compare is strand-dependent. Derive it from the mapped
alignment and the known library orientation; validate the convention on a small
set of reads before making a hard 3′-distance filter.

## 7. Molecule grouping and counting

The molecule count is one per corrected UMI per cell per gene:

```text
alignment → gene assignment → group (CB, UMI, gene) → UMI correction → count
```

Keep the standard tags (`CB`, `CR`, `UR`, `UB`) even if a composite sort key is
introduced. They preserve provenance and make the correction rules auditable.

`samtools sort -t TAG` sorts on exactly one auxiliary tag. For an efficient
pre-correction external sort, create a custom composite tag such as:

```text
ZU:Z:<CB>|<UR>
```

then run:

```bash
samtools sort -t ZU -@ 8 -o sample.cb_ur.bam sample.coordinate.bam
```

`ZU` is a transport/sort key, not a final molecule identity. After gene assignment,
the relevant grouping is `(CB, UB, gene)` (or the corrected descendant of `UR`),
not barcode+UMI alone. Avoid using `MI` for the provisional key because `MI` is
normally reserved for a final molecule identifier.

## Reference deliverables

Archive and record versions/checksums for these files with the analysis:

```text
GRCh38.primary_assembly.genome.fa
GRCh38.primary_assembly.gencode-v50.splice-sr.mmi
gencode.v50.primary_assembly.annotation.gtf
gencode.v50.primary_assembly.junctions.bed
genes.mtx
```

The index, annotation, sbxpress minimap version, exact command line, and feature-assignment
rules together define the counting reference. Changing any of them can change the
gene-expression matrix.
