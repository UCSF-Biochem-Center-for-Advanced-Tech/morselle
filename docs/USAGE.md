# MORSELLE command reference

Detailed reference for each `morselle` stage: what `extract` does to an SBX
read, every stage command, `extract` options and SAM tags, outputs, and the
`extract` architecture. For complete workflows see [MANUAL.md](MANUAL.md); for
installation see [INSTALL.md](INSTALL.md).

## What it does

Each SBX read holds a complete library molecule:

```
leader – P5 – i5 – R1primer – CB(16) – UMI(12) – polyT – cDNA – adapter – i7 – rc(P7)
```

For every read, `morselle`:

1. **Orients** it to a canonical P5→P7 frame (auto-detecting per file; each SBX run
   is orientation-pure but the two runs of a library are opposite strands).
2. **Extracts** R1 (the 16 bp cell barcode + 12 bp UMI) and R2 (the cDNA between the
   poly(dT) primer and the 3′ read-through adapter).
3. **Trims** the 3′ Illumina read-through adapter from the cDNA — including the
   homopolymer-indel variants SBX produces (`AAGAGC`→`AAAGAGC`) and 3′-end partials.
4. **Gates / rescues** the cell barcode against the 10x whitelist (exact, Hamming-1,
   and optional homopolymer-reframe / indel-realign rescues).
5. **Demultiplexes** by i5/i7 against the sample sheet (offset-tolerant matching).
6. **Emits** per-sample R1/R2 pairs in 10x reference naming, plus a full metrics report.

R2 is written **reverse-complemented** (with reversed quality) so it follows the
Illumina paired-end convention and maps sense to the gene.

---

## Stage commands

`morselle matrix` reconstructs a 10x reference
raw gene-expression matrix from the corrected molecule annotations in a complete
10x reference BAM:

```bash
morselle matrix possorted_genome_bam.bam \
  --features raw_feature_bc_matrix/features.tsv.gz \
  --cell-barcodes resources/10x_barcodes/737K-arc-v1.txt.gz \
  --genome GRCh38 \
  --output outs
```

It counts the `xf` bit-8 molecule representatives using the corrected `CB`, `UB`,
`GX`, and `GN` tags. The feature/barcode templates preserve zero-count rows and
columns and 10x reference ordering; without them, only observed features and
barcodes can be recovered. The 10x-style output layout is:

```text
outs/
├── raw_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── raw_feature_bc_matrix.h5
└── metrics_summary.json
```

The H5 file stores the same matrix in 10x CSC format. `--genome` supplies its
feature genome label because that value cannot be recovered reliably from a BAM.
Every run also prints an stderr `RUN SUMMARY`.

Call cells from the raw HDF5 matrix and create 10x-style filtered
matrices with the two-stage 10x GEX procedure (OrdMag followed by multinomial
EmptyDrops):

```bash
morselle cells outs/raw_feature_bc_matrix.h5 --output outs
```

For NextGEM/5′ chemistry, defaults follow 10x reference: ambient barcode ranks
45,000–90,000, at least 500 UMIs (raised above the ambient maximum), 100,000
Monte Carlo simulations, and 0.001 Benjamini-Hochberg FDR. Use
`--expect-cells N` when recovery is known, or `--force-cells N` to select exactly
the top N nonzero barcodes without EmptyDrops. Outputs include
`filtered_feature_bc_matrix/`, `filtered_feature_bc_matrix.h5`,
`filtered_barcodes.csv`, `cell_calls.csv`, and `cell_calling_metrics.json`.

### Raw SBX to counts pipeline

`scripts/morselle-gex-pipeline.sh` is the restart-safe production wrapper for
one demultiplexed sample. It performs extraction once for all supplied device
outputs, then maps, annotates in unordered parallel mode, coordinate-sorts,
resolves molecules, writes raw counts, and calls cells. The comparison audit is
optional so normal count generation does not make an additional full-BAM pass.

```bash
scripts/morselle-gex-pipeline.sh \
  --fastqs raw/device_0/out raw/device_1/out raw/device_2/out raw/device_3/out \
  --bclconvert-samplesheet samples.csv --sample M2_GEX \
  --whitelist 737K-arc-v1.txt.gz \
  --mmi GRCh38.splice-sr.mmi \
  --gtf gencode.v44.primary_assembly.basic.annotation.gtf.gz \
  --junctions gencode.v44.primary_assembly.junctions.bed \
  --out-dir results/sbx_p5 \
  --map-threads 32 --annotate-threads 8 --sort-threads 8
```

Use separate thread values for concurrent runs; the defaults are intentionally
conservative. The script writes its canonical `genes.mtx` to
`OUT/reference/genes.mtx` unless `--genes-mtx` is supplied.

The installed `morselle annotate` stage assigns alignments against a GENCODE-style
GTF and adds 10x-compatible `GX`, `GN`, `TX`, `AN`, and `RE` tags:

```bash
morselle annotate sample.coordinate.bam \
  --gtf gencode.v50.primary_assembly.annotation.gtf.gz \
  --output sample.annotated.bam \
  --three-prime \
  --threads 8
samtools index -@ 8 sample.annotated.bam
```

For a single-sort streaming workflow, use `-` for annotator stdin/stdout and
provide an explicit metrics path:

```bash
morselle minimap ... \
| morselle annotate - --gtf genes.gtf.gz --output - \
    --metrics sample.annotation.metrics.json --three-prime --unordered \
| samtools sort -@ 8 -o sample.annotated.coordinate.bam -
```

By default annotation preserves input order. `--unordered` annotates record
batches across `--threads` workers and emits records as workers finish; use it
when a downstream `samtools sort` will establish coordinate order. Existing
barcode/UMI tags are preserved in either mode. UMI correction,
molecule selection, and `xf` tagging remain a separate downstream stage.
By default it recreates the 10x GENCODE gene/transcript-biotype,
readthrough-transcript, PAR, and Ensembl ID-version rules while loading the GTF.

For mapper-ambiguity diagnostics, use `--multimap-tags` with minimap2 secondary
alignments retained (the production wrapper enables `-N 10`, plus the primary,
to match the audit's top-coded `NH=11+` bucket). It
streams query-grouped SAM directly through annotation and writes `NH` as the
number of **reported** primary/secondary hits, `HI` as the hit index, and `mm:i:1`
when multiple reported hits resolve to exactly one compatible gene. It preserves
minimap2's native MAPQ rather than applying 10x reference's MAPQ=255 promotion, so
the `mm` tag records the equivalent gene-rescue condition rather than an identical
STAR/10x reference mapping action. This mode is intentionally sequential and cannot
use `--unordered`; it avoids a name-sorted intermediate BAM.
The standard tags retain all compatible calls. Local-use `XC`/`XT` tags always
record gene/transcript collision degree. With the explicit `--three-prime` assay
flag, `ZG`/`ZT` and `XP` additionally record the subset supported by a
strand-aware terminal exon. These preference tags never remove an alignment or
rewrite `GX`/`TX`, and are absent unless the flag is given.

Run the molecule stage after annotation. It corrects UMIs within each `(CB, GX)`
group, resolves the same corrected `(CB, UB)` landing on multiple genes by total
read support, and marks exactly one counted representative per resolved molecule:

```bash
morselle molecules sample.annotated.coordinate.bam \
  --output sample.molecules.coordinate.bam --threads 8 \
  --evidence-workers 16 \
  --saturation sample.saturation.tsv.gz --saturation-interval 1000000
samtools index -@ 8 sample.molecules.coordinate.bam
```

Eligible reads receive `xf=17`; the selected representative receives `xf=25`
(10x reference bit 8). `XV` records the molecule-level gene collision degree and
`XW` records `U`nique, `W`inner, `L`oser, or unresolved `T`ie. Ties are not counted.

`--saturation` is collected during the molecule evidence pass (it does not add a
BAM pass). Its gzip TSV has one row per requested eligible-read depth for the
whole run (`scope=run`) and each observed cell barcode (`scope=cell`). The
`molecules` values are gene-aware, UMI-corrected, collision-resolved molecules;
`genes` is the number of detected genes. Read names are stably hashed before
downsampling, so coordinate-sorted input cannot bias the nested curve.

For an indexed coordinate BAM, `--evidence-workers 16` partitions the first
evidence pass into disjoint genomic ranges. The raw UMI evidence is merged before
correction, so molecule counts remain global and exact. Start with 16 workers;
raising this to 32 requires substantially more transient RAM during the merge.

Audit a molecule BAM with comparison-ready mapping, annotation, collision, UMI,
and sensitivity metrics:

```bash
morselle audit sample.molecules.coordinate.bam \
  --label SBX --output sample.audit.json \
  --molecules sample.molecules.tsv.gz --molecule-sample 64 \
  --cells filtered_feature_bc_matrix/barcodes.tsv.gz
```

The JSON includes mapping and region counts, `NH`/MAPQ/soft-clip ambiguity
diagnostics, gene/transcript compatibility degrees, exact `xf` bit/value counts,
UMI-corrected records, SBX `XW` outcome record counts, and a sampled molecule
funnel (raw and corrected UMI nodes, gene candidates, eligible nodes, and
representatives). It also reports the read-support histogram and a read-depth
saturation curve. With `--cells`, representatives are additionally split into
all barcodes versus final called cells. Consistent molecule subsamples support
platform-only molecule estimates, gene agreement, and per-cell/per-gene
concordance in `scripts/bam_audit_report.py`.

The three coordinate-sliced 10x reference fixtures can be checked with:

```bash
./tests/validate_cr_bam_to_mtx.sh /path/to/test_fixtures
```

See **[INSTALL.md](INSTALL.md)** for system requirements, the barcode whitelist,
throughput tuning, and troubleshooting.

### Barcode whitelist (required)

This is **Single Cell Multiome (ATAC + GEX)** data, GEX side, so the cell barcodes come
from the **ARC** gel-bead list **`737K-arc-v1.txt.gz`** (736,320 barcodes). It is 10x
reference data staged under the gitignored `resources/10x_barcodes/`
directory. Obtain it from a 10x reference ARC distribution:

```
$TENX_INSTALL/lib/python/barcodes/737K-arc-v1.txt.gz
```

and pass it at runtime with `--whitelist`.

For stable feature rows across samples, build one canonical `genes.mtx` from the
reference GTF with `morselle annotate --genes-mtx`, then pass that same file to
every `morselle matrix --features` invocation. Existing manifests are verified
against the filtered annotation rather than overwritten.

---

## Usage

For the end-to-end genomic minimap2 reference and counting workflow, see
[Preparing a minimap2 reference for SBX 3′ GEX counting](REFERENCE_PREPARATION.md).

Single directory of run files, all cores, ready for 10x reference:

```bash
morselle extract \
  --fastqs   /data/CAT02/raw_fastqs_run1/device_0/out \
  --bclconvert-samplesheet /data/CE15204_bclconvert_input.csv \
  --whitelist    /usr/local/share/morselle/737K-arc-v1.txt.gz \
  --out-dir  /data/out/run1 \
  --threads  0 \
  --cellranger
```

Explicit files across multiple devices:

```bash
morselle extract --fastqs /data/.../device_{0,1,2,3}/out/PrimaryAnalysis_*.fastq.gz \
  --bclconvert-samplesheet sheet.csv --whitelist 737K-arc-v1.txt.gz -o /data/out/run1 --cellranger
```

Process each SBX run separately (the two runs are the same library in opposite
orientations; orientation is auto-detected and enforced per file).

### Options

| Flag | Default | Purpose |
|---|---|---|
| `--fastqs` | (required) | Input SBX FASTQs: a directory of `PrimaryAnalysis_*.fastq.gz`, or explicit files. |
| `--bclconvert-samplesheet` | — | BCL Convert v2 CSV (primary) or a flat `sample,i7,i5` CSV. |
| `--whitelist` | — | 10x cell-barcode whitelist (`737K-arc-v1.txt.gz`). |
| `--out-dir` / `-o` | `morselle_out` | Output directory. |
| `--threads` | 0 (all cores) | Worker pool size. |
| `--cellranger` | off | **10x reference compatibility preset** — applies its required settings (currently `--max-r2-len 649`). Use this for 10x reference input. |
| `--minimap2` | off | Emit only R2 FASTQs. Adds 10x reference barcode/UMI SAM tags to each FASTQ comment; align with `morselle minimap -y`. Cannot be combined with `--cellranger`. |
| `--processing-tags` | off | Add SBX processing-provenance tags (`XO`, `XB`, `XI`, `XA`, `XL`) to R2 comments. Implied by `--minimap2`. |
| `--orientation` | auto | Force `p5-led` / `p7-led` if auto-detect is wrong. |
| `--cb-gate` | on | Drop pairs whose CB is neither on-whitelist nor within Hamming-1. |
| `--cb-reframe` | off | Rescue homopolymer-indel frame-shifted barcodes (`XB:i:3`). Off by default: on CAT02 Multiome, scored against Illumina molecules from the same library, only about half of reframed barcodes are right (Hamming-1 78–83%, realign 97–99%). |
| `--r1-tail-relock` | off | For barcode failures, re-find the R1-primer tail allowing up to two insertions/deletions/mismatches and read the barcode + UMI after it; accepted only on a whitelist (or unique within-1) barcode (`XB:i:7`). On CAT02 Multiome: +0.54% / +0.23% of reads assigned; rescued barcodes land on cell-like barcodes at the exact-match rate. A re-locked read the index demux left unassigned is demultiplexed again from the re-located barcode (`whitelist.relock_reanchored`). |
| `--index-window` | 3 | Offset-tolerant i5/i7 matching (absorbs indels shifting the index). |
| `--index-indels` | off | Also match an i5/i7 with one insertion or deletion inside the index, over the same `--index-window` frames, only when mismatch-only matching finds nothing. Enabled per index only if the sample sheet is safe for it (logged). On CAT02 Multiome: +0.2–0.7% of reads assigned and +2 points of dual-index assignment, no read moved between samples. |
| `--index-max-mm` | 1 | Max mismatches when matching sample-sheet indices. |
| `--demux-index` | auto | `dual`, `both`, `i5` or `i7`. Auto = `dual` with the per-read rule: an index alone assigns only when it names exactly one sample; an index shared by several samples (combinatorial sheets, or overlapping mismatch neighbourhoods) is narrowed to one by the other index, or the read stays undetermined. |
| `--realign` | on | Single-indel CB rescue — recovers ~+1% of reads at high precision (`--realign false` to disable). |
| `--cb-insertion` | off | For confidently located barcode failures, accept the unique whitelist barcode that removes one extra base inside the 17 bp at the barcode (`XB:i:8`); the UMI is the 12 bp after it. Insertions only (deletion readings land near chance). On CAT02 Multiome: +0.89% / +0.53% of reads assigned; rescued barcodes land on cell-like barcodes at 48.6% / 19.0% (Hamming-1 47.2% / 24.2%, reframe 37.8% / 19.5%). |
| `--lenient` | off | Within-2 CB rescue — leave off (~51% precise on the dense ARC whitelist). |
| `--min-cdna-len` | 50 | Shortest cDNA emitted as R2. |
| `--max-r2-len` | 0 (uncapped) | Cap emitted R2 length (seq+qual together). 10x reference hard-caps at 649 bp; `--cellranger` sets this. |
| `--gzip-level` | 4 | Output gzip compression level (lower = faster/larger). |
| `--split` | 1 | Split each sample across N gzip files (`_001`..`_00N`) so N threads compress one sample in parallel. |
| `--reader-threads` | 0 (auto) | I/O parallelism. Auto uses every core even for a single input file (see below); `N>0` forces the pipelined path with N readers. |

**Parallel gzip input decode (experimental, opt-in).** `MORSELLE_RAPIDGZIP=1`
decodes each gzip input across many threads via
[rapidgzip-core](https://crates.io/crates/rapidgzip-core) instead of the serial
flate2 reader, lifting the single-file ceiling described below. See
[Input parallelism](#architecture) for when it helps and what it costs.

| Variable | Default | Meaning |
|---|---|---|
| `MORSELLE_RAPIDGZIP` | off | `1`/`true`/`on`/`yes` enables parallel input decode |
| `MORSELLE_RAPIDGZIP_THREADS` | worker count | Decoder threads per input file |
| `MORSELLE_RAPIDGZIP_POOL` | worker count | Process-wide ceiling on decode slots |

Run `morselle extract --help` for extraction options. Set `RUST_LOG=debug` for verbose logging.

### 10x reference compatibility

Two quirks of 10x reference are handled for you:

- **Hard 649 bp read cap.** A read ≥ 650 bp fails with STAR `quality string length
  is not equal to sequence length`. `--cellranger` caps R2 to 649 bp — applied at
  write time, after the full-length fragment has been adapter-trimmed and
  reverse-complemented, so the retained bases are the ones the R2 primer would read.
- **Single-member gzip.** 10x reference mis-frames records at multi-member gzip
  boundaries, so each FASTQ is written as a single gzip member (compressed once,
  during processing).

### Minimap2-ready output

Use `--minimap2` when R2 will be aligned directly with minimap2. It emits only
`*_R2_*.fastq.gz`; there is deliberately no R1 file. Each R2 FASTQ comment contains
the following SAM-style tags, which `morselle minimap` preserves when invoked with `-y`:

- `CR` / `UR`: raw cell barcode / raw UMI.
- `CB`: whitelist-corrected cell barcode, with the 10x reference `-1` gem-group suffix.
- `UB`: provisional corrected UMI (the extracted UMI before gene-aware UMI correction).
- `CY` / `UY`: raw barcode and UMI Phred qualities.
- `XO`: original SBX read orientation (`5` = P5-led, `7` = P7-led).
- `XB`: barcode route: 0 exact, 1 Hamming-1, 2 indel realignment, 3 reframe,
  4 Hamming-2, 5 no whitelist, 6 unrescued miss emitted with gating disabled,
  7 R1-tail re-lock (`--r1-tail-relock`), 8 CB insertion (`--cb-insertion`).
- `XI`: index route: 0 dual, 1 i5-only, 2 i7-only, 3 conflict, 4 no usable
  index, 5 demultiplexing not requested.
- `XA`: anchor bitmask: 1 R1-primer tail, 2 poly-T, 4 cDNA adapter trimmed,
  8 file-level expected/detected orientation applied.
- `XL`: extracted cDNA length before `--max-r2-len` capping.

The five SBX tags are also available in paired-FASTQ mode with
`--processing-tags`; they are placed on R2, the record that is aligned.

For example:

```bash
morselle extract --fastqs /data/sbx --whitelist 737K-arc-v1.txt.gz -o /data/out --minimap2
morselle minimap -a -x splice -uf -k14 --secondary=no -y reference.mmi /data/out/*_R2_*.fastq.gz \
  | samtools sort -o aligned.bam
```

Do not create a combined barcode+UMI tag at this stage. Group records by `(CB, UR)`
for UMI correction, then by `(CB, UB, gene)` after gene assignment; an `MI` tag, if
needed, should be generated only at that downstream molecule-calling step.

### Illumina input for matched-platform comparisons

`morselle illumina` converts standard paired 10x Illumina FASTQs into the same
R2-only minimap2 handoff used by `morselle extract --minimap2`. R1 supplies the
16 bp cell barcode and 12 bp UMI; R2 is emitted unchanged because Illumina R2 is
already in mapping orientation. Exact and uniquely resolvable Hamming-1 cell
barcodes are retained, while invalid, absent, and ambiguous barcodes are gated.

Repeat `--r1` and `--r2` in matching order for multiple lanes or shards. File
pairs are converted concurrently (`--threads` limits the lane workers) and a
separate numbered output shard is created for each pair:

```bash
morselle illumina \
  --r1 sample_L001_R1_001.fastq.gz --r2 sample_L001_R2_001.fastq.gz \
  --r1 sample_L002_R1_001.fastq.gz --r2 sample_L002_R2_001.fastq.gz \
  --whitelist resources/10x_barcodes/3M-february-2018.txt.gz \
  --output sample.illumina.minimap.fastq.gz \
  --threads 2 --split 4

morselle minimap -a -x splice -m 20 -s 30 --secondary=no -y reference.mmi \
  sample.illumina.minimap.shard-*.fastq.gz \
  | morselle annotate - --gtf genes.gtf --output -
```

For one input pair with the default `--split 1`, the requested output filename
is retained. Otherwise, outputs are named `*.shard-001.fastq.gz`,
`*.shard-002.fastq.gz`, and so on. `--split N` routes records to N independent
level-4 gzip writer threads per input pair, removing the single-encoder
bottleneck for a large single lane. The files can all be passed directly to
minimap2, as above. For maximum single-pair throughput, combine this with
parallel gzip decode:

```bash
MORSELLE_RAPIDGZIP=1 MORSELLE_RAPIDGZIP_THREADS=32 MORSELLE_RAPIDGZIP_POOL=32 \
morselle illumina ... --threads 1 --split 8
```
Their comments carry `CR`, `CB`, `UR`, `UB`, `CY`, and `UY`. A combined JSON
report lists every shard and is written to `<output>.metrics.json` unless
`--metrics` specifies another path.

---

## Outputs

In `--out-dir`:

- `<sample>_S<n>_L001_R1_001.fastq.gz` / `_R2_001.fastq.gz` — per-sample pairs
  (R1 = 28 bp CB+UMI, R2 = cDNA).
- With `--minimap2`, only `<sample>_S<n>_L001_R2_001.fastq.gz` is produced; its
  FASTQ comments carry the barcode/UMI SAM tags described above.
- `Undetermined_S0_L001_R{1,2}_001.fastq.gz` — pairs with no matchable index.
- `metrics_summary.json` — full run accounting: totals, timing, read classes,
  CB-gate stats, demux breakdown, per-sample counts. `demux.index_evidence`
  counts, for every read whose R1 primer was found: `full` (i5 and i7 both
  seen and identified to the same sample) versus `partial`; which index
  regions are in the read (`seen_both` / `seen_i5_only` / `seen_i7_only` /
  `seen_none`); which matched a sample (`identified_both` / `_i5_only` /
  `_i7_only` / `_none`; `both` includes conflicts; a match to an index shared
  by several samples counts as identified); how many matched only through
  `--index-indels`; how many matched a shared index
  (`identified_i5_shared` / `_i7_shared`); and how many abstained because a
  shifted-frame mismatch hit had a competing sample (`i5_float_tie` /
  `i7_float_tie`).
- `sample_sheet.csv` — the normalized flat sample sheet.

Every run also prints a **RUN SUMMARY** to stderr (elapsed, threads, reads/sec,
`processed / emitted / assigned / lost`, and a lost breakdown). Capture it with
`2> run1.log`.

---

## Architecture

`morselle extract` modules in the library crate (`src/`); the other stages live in
`annotate/`, `molecules/`, `cr_bam_to_mtx/` (matrix), `cell_calling/`, `illumina/`,
`probe_demux/`, `probe_set/`, `feature_capture/`, and `protospacer_calling/`:

| Module | Role |
|---|---|
| `orient/` | Per-read orientation → canonical P5→P7 frame. |
| `homopolymer.rs` | poly(dT)/poly(A) run detection (the barcode/cDNA boundary). |
| `extract/` | R1/R2 bracketing + 3′ adapter trimming (indel-tolerant). |
| `barcode.rs` | Whitelist load + gate/rescue of the cell barcode. |
| `demux/` | i5/i7 extraction + offset-tolerant sample-sheet matching. |
| `sample_sheet.rs` | Sample-sheet parsing (BCL Convert v2 / flat CSV). |
| `anchors.rs` | Adapter / primer / leader constants — single source of truth. |
| `io.rs` | FASTQ input discovery and single-member gzip output writers. |
| `accounting.rs` | Metrics aggregation → `metrics_summary.json`. |
| `pipeline.rs` | Core per-read types and state. |
| `config.rs` | CLI definition and runtime parameters. |
| `pipeline_runner/` | `extract` orchestration: parallel per-file processing and emission. |

FASTQ parsing uses [needletail](https://crates.io/crates/needletail) and gzip uses
[flate2](https://crates.io/crates/flate2) with the zlib backend.

**Input parallelism.** Serial gzip inflate caps one file to one core, so the executor
adapts to the file count (override with `--reader-threads`):

- **Many files (files ≥ workers):** one file per worker thread via
  [rayon](https://crates.io/crates/rayon) — up to `workers` concurrent decodes, each
  thread decoding and processing its own file. Optimal when files are plentiful.
- **Few files (files < workers, e.g. one large file):** decode is decoupled from
  processing — one reader thread per file streams record batches over a bounded
  channel to the full worker pool, so even a single-file run saturates every core.
  A single huge file is then bounded only by serial gzip inflate throughput.

That last ceiling is real: one flate2 reader sustains roughly 0.9M reads/s, and
per-read processing runs about 39k reads/s per core, so a single-file run stops
scaling near two dozen workers no matter how many cores the node has. Setting
`MORSELLE_RAPIDGZIP=1` decodes one gzip member across many threads and moves that
ceiling out, but the gain depends entirely on how the input was framed — BGZF and
`pigz --independent` parallelize very well, a plain `gzip` stream moderately, and
`pigz`'s default dependent-block framing not at all. Measure your own inputs
before relying on it; on a small node, or when files ≥ workers (where the per-file
path already decodes many files at once), it changes little.

The feature is decode-only — output writers are untouched, and every emitted FASTQ
remains a single gzip member. Enabling it does not change which records are
emitted, and on the per-file path output stays byte-identical; on the pipelined
path record order varies between runs with either decoder (see below).

Output stays mate-ordered (the Nth R1 record is always the mate of the Nth R2) in
both modes; the pipelined path does not preserve input record order across the
file, and because record order there depends on worker scheduling, its output is
not byte-reproducible between runs. The per-file path is byte-reproducible.
