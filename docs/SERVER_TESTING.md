# v0.6.0 real-data acceptance testing

Local regression checks are complete; server real-data acceptance is pending.
Use a fresh output directory for every version, mode, or changed configuration.
Do not reuse 0.5.x intermediates or the legacy hardcoded evaluation scripts for
release acceptance: their existence-only restart checks can reuse stale results.

## Install and verify on the Linux server

Install the prerequisites in INSTALL.md first, including Rust, CMake, clang,
samtools, minimap2, k8 and paftools.js. The pipeline also requires Python 3.
Keep the checkout and installation separate from the current production build.

```bash
git clone --branch v0.6.0 https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/sbxpress.git sbxpress-v0.6.0
cd sbxpress-v0.6.0
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install numpy python-pptx h5py
scripts/check.sh > server-check.log 2>&1
./install.sh --prefix "$HOME/.local/sbxpress-v0.6.0" --skip-packages --skip-rust
export PATH="$HOME/.local/sbxpress-v0.6.0/bin:$PATH"
sbxpress --version
sbxpress molecules --version
sbxpress matrix --version
git rev-parse HEAD > release-commit.txt
```

Confirm all three version commands report 0.6.0. Archive server-check.log and
release-commit.txt with the results. Record CPU/RAM, thread counts, reference
checksums, input file identities and versions of samtools/minimap2. Do not copy
the workstation's macOS binaries to Linux.

## Start with a representative sample

Choose a complete small sample or a representative raw FASTQ subset. Retain
duplicate reads, multiple genes in GX, existing UB/ZG and PB tags where present;
avoid selecting only idealized records. Start with moderate worker counts and
monitor peak memory before scaling to a full library. Sort memory is per worker.

Replace the example paths and sample name below with server inputs. The whitelist
must match the assay; the ARC whitelist is appropriate only for ARC-derived GEX.
The splice index, junction BED and annotation must use the same genome build.

```bash
SAMPLE=M2_GEX
FASTQS=/data/run/device_0/out
SHEET=/data/run/bclconvert_input.csv
WHITELIST=/data/reference/737K-arc-v1.txt.gz
MMI=/data/reference/genome.splice-sr.mmi
GTF=/data/reference/annotation.gtf.gz
JUNCTIONS=/data/reference/junctions.bed
RUN=/data/validation/sbxpress-v0.6.0-gex
mkdir -p "$RUN"
/usr/bin/time -v -o "$RUN/resources.txt" bash scripts/sbxpress-pipeline.sh \
  --input-mode sbx --analysis-mode gex --sample "$SAMPLE" \
  --fastqs "$FASTQS" --bclconvert-samplesheet "$SHEET" \
  --whitelist "$WHITELIST" --mmi "$MMI" --gtf "$GTF" \
  --junctions "$JUNCTIONS" --out-dir "$RUN" --threads 8 \
  --sort-mem 1G --audit > "$RUN/console.log" 2>&1
```

Inspect extraction, annotation, molecule and cell metrics; raw/filtered H5 and
MEX outputs should exist under RUN/SAMPLE/outs. Check logs for errors and BAMs
with samtools quickcheck. Rerun the identical command: completed stages should
be reused. A changed input or analysis configuration must require a new output
directory. Keep the pipeline-state manifests and completion markers for audit.

For Illumina, use the same driver with --input-mode illumina and paired --r1/--r2
instead of SBX input options. For Flex, use --analysis-mode flex and the actual
probe reference/whitelists and sample sheet documented in MANUAL.md. Include a
multi-sample PB case and verify identical CB/UMI combinations remain independent
across samples. Each assay needs a separate output directory.

## Compare molecule execution strategies

Use the same coordinate-sorted, indexed annotated BAM and feature manifest for
both runs. Select a BAM small enough for the unsharded evidence pass. Set BAM,
FEATURES and PARITY below; append identical assay flags to both molecule commands
when testing tag-seq or Flex. The sharded BAM is unsorted: sort before indexing.

```bash
BAM=/data/validation/sample.annotated.coordinate.bam
FEATURES=/data/validation/genes.mtx
PARITY=/data/validation/sbxpress-v0.6.0-parity
mkdir -p "$PARITY"
sbxpress molecules "$BAM" --evidence-workers 4 --threads 8 \
  --output "$PARITY/indexed.bam" --metrics "$PARITY/indexed.metrics.json"
sbxpress molecules "$BAM" --shards 8 --shard-concurrency 1 --threads 8 \
  --output "$PARITY/sharded.bam" --metrics "$PARITY/sharded.metrics.json"
samtools sort -@ 8 -o "$PARITY/sharded.sorted.bam" "$PARITY/sharded.bam"
samtools index "$PARITY/sharded.sorted.bam"
samtools quickcheck -v "$PARITY/indexed.bam" "$PARITY/sharded.sorted.bam"
sbxpress matrix "$PARITY/indexed.bam" --features "$FEATURES" \
  --output "$PARITY/indexed" --no-gzip
sbxpress matrix "$PARITY/sharded.sorted.bam" --features "$FEATURES" \
  --output "$PARITY/sharded" --no-gzip
for file in matrix.mtx features.tsv barcodes.tsv; do
  cmp "$PARITY/indexed/raw_feature_bc_matrix/$file" \
      "$PARITY/sharded/raw_feature_bc_matrix/$file" || exit 1
done
```

Compare biological totals in the metrics, excluding execution configuration and
timing. BAM byte equality is not expected across strategies because order and
headers differ. Investigate matrix differences before accepting the release.

## Evaluate tag-seq counting policy

Use an independently counted Illumina matrix for the same library. Choose shared
(the default), fractional and position-aware modes using count concordance and
expression-dependent bias, rather than assuming one mode is universally correct.
Use a new directory: the comparison helper reuses existing matrices by existence.

```bash
bash scripts/compare_tagseq_modes.sh \
  --bam /data/tagseq/sample.annotated.coordinate.bam \
  --reference /data/illumina/outs/raw_feature_bc_matrix.h5 \
  --features /data/reference/genes.mtx \
  --out-dir /data/validation/sbxpress-v0.6.0-tagseq --threads 8 \
  --modes "shared fractional position position_frac"
```

Review summary.tsv, per_gene.tsv and molecule metrics for every mode. Fractional
counts require explicit barcode or forced cell selection; automatic EmptyDrops
rejects them. MEX retains full precision; H5 stores float32. Fractional saturation
is a full-depth approximation, and saturation is not available with sharding.

## Acceptance and evidence to return

- Required server checks pass and installed binaries report 0.6.0.
- Every tested assay finishes, produces readable BAM/H5/MEX outputs, and retains
  correct sample identities and reference names.
- Indexed/sharded matrix comparisons pass; duplicates yield one representative
  per selected molecule and representative quality choices remain correct.
- Restart behavior passes; an interrupted run does not publish partial matrix
  or cell outputs as completed stages.
- Compare the same inputs/references against v0.5.9 in a separate installation
  and output directory. Explain differences using the v0.6.0 changelog: tag-seq
  shared defaults, optional position-aware counting, corrected PB identities,
  duplicate representatives and fractional precision can change results.
- Archive commands, commit/version, logs, manifests, metrics, comparisons and
  elapsed time/peak RSS. Unexpected count changes, crashes, sample mixing or
  unbounded memory growth block real-data acceptance.

After the representative runs pass, scale to one full library per assay and
repeat the count and resource checks before replacing the production build.
