# sbxpress

**Convert raw Roche SBX long reads into per-sample 10x R1/R2 FASTQs for Cell Ranger.**

Roche SBX sequences each library molecule end-to-end as one long read (P5→P7),
so a 10x Single Cell Multiome **GEX** library comes off the instrument as full
molecules rather than the paired short reads Cell Ranger expects. `sbxpress` parses
each molecule, recovers the cell barcode / UMI / cDNA, demultiplexes by sample
index, and writes standard 10x-style **R1** (barcode+UMI) and **R2** (cDNA) FASTQ
pairs that feed directly into `cellranger-arc count`.

It is built for scale: a single multi-core node processes tens of billions of reads
in well under two hours, and output is deterministic (identical inputs → byte-identical
FASTQs).

---

## What it does

Each SBX read holds a complete library molecule:

```
leader – P5 – i5 – R1primer – CB(16) – UMI(12) – polyT – cDNA – adapter – i7 – rc(P7)
```

For every read, `sbxpress`:

1. **Orients** it to a canonical P5→P7 frame (auto-detecting per file; each SBX run
   is orientation-pure but the two runs of a library are opposite strands).
2. **Extracts** R1 (the 16 bp cell barcode + 12 bp UMI) and R2 (the cDNA between the
   poly(dT) primer and the 3′ read-through adapter).
3. **Trims** the 3′ Illumina read-through adapter from the cDNA — including the
   homopolymer-indel variants SBX produces (`AAGAGC`→`AAAGAGC`) and 3′-end partials.
4. **Gates / rescues** the cell barcode against the 10x whitelist (exact, Hamming-1,
   and optional homopolymer-reframe / indel-realign rescues).
5. **Demultiplexes** by i5/i7 against the sample sheet (offset-tolerant matching).
6. **Emits** per-sample R1/R2 pairs in Cell Ranger naming, plus a full metrics report.

R2 is written **reverse-complemented** (with reversed quality) so it follows the
Illumina paired-end convention and maps sense to the gene.

---

## Install

From the source directory on Ubuntu/Debian or macOS:

```bash
./install.sh
```

This installs build prerequisites, `minimap2`, `samtools`, and the `k8`/
`paftools.js` sbxpress minimap companions, plus the Rust toolchain (via rustup if absent).
It downloads the official `k8` runtime, builds the optimized binaries, and installs
them under `/usr/local/bin`. Linux uses `apt-get`; macOS uses Homebrew. User-local
install without sudo:

```bash
./install.sh --prefix "$HOME/.local"
```

Verify:

```bash
sbxpress --version
sbxpress --help
sbxpress extract --help
sbxpress minimap --version
sbxpress annotate --version
sbxpress annotate --help
sbxpress molecules --help
sbxpress audit --help
sbxpress matrix --version
sbxpress matrix --help
```

Uninstall using the same prefix supplied at installation:

```bash
./uninstall.sh
./uninstall.sh --prefix "$HOME/.local"
```

The uninstaller preserves the potentially shared `k8` and `paftools.js` tools by
default. Pass `--remove-companions` to remove them, `--keep-data` to retain the
installed whitelist, or `--dry-run` to preview all actions.

The repository also installs `sbxpress matrix`, which reconstructs a Cell Ranger
raw gene-expression matrix from the corrected molecule annotations in a complete
Cell Ranger BAM:

```bash
sbxpress matrix possorted_genome_bam.bam \
  --features raw_feature_bc_matrix/features.tsv.gz \
  --cell-barcodes resources/cellranger_barcodes/737K-arc-v1.txt.gz \
  --genome GRCh38 \
  --output outs
```

It counts the `xf` bit-8 molecule representatives using the corrected `CB`, `UB`,
`GX`, and `GN` tags. The feature/barcode templates preserve zero-count rows and
columns and Cell Ranger ordering; without them, only observed features and
barcodes can be recovered. The Cell Ranger-style output layout is:

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

Call cells from the raw HDF5 matrix and create Cell Ranger-style filtered
matrices with the two-stage 10x GEX procedure (OrdMag followed by multinomial
EmptyDrops):

```bash
sbxpress cells outs/raw_feature_bc_matrix.h5 --output outs
```

For NextGEM/5′ chemistry, defaults follow Cell Ranger: ambient barcode ranks
45,000–90,000, at least 500 UMIs (raised above the ambient maximum), 100,000
Monte Carlo simulations, and 0.001 Benjamini-Hochberg FDR. Use
`--expect-cells N` when recovery is known, or `--force-cells N` to select exactly
the top N nonzero barcodes without EmptyDrops. Outputs include
`filtered_feature_bc_matrix/`, `filtered_feature_bc_matrix.h5`,
`filtered_barcodes.csv`, `cell_calls.csv`, and `cell_calling_metrics.json`.

### Raw SBX to counts pipeline

`scripts/sbxpress-gex-pipeline.sh` is the restart-safe production wrapper for
one demultiplexed sample. It performs extraction once for all supplied device
outputs, then maps, annotates in unordered parallel mode, coordinate-sorts,
resolves molecules, writes raw counts, and calls cells. The comparison audit is
optional so normal count generation does not make an additional full-BAM pass.

```bash
scripts/sbxpress-gex-pipeline.sh \
  --fastqs raw/device_0/out raw/device_1/out raw/device_2/out raw/device_3/out \
  --sample-sheet samples.csv --sample M2_GEX \
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

The installed `sbxpress annotate` stage assigns alignments against a GENCODE-style
GTF and adds Cell Ranger-compatible `GX`, `GN`, `TX`, `AN`, and `RE` tags:

```bash
sbxpress annotate sample.coordinate.bam \
  --gtf gencode.v50.primary_assembly.annotation.gtf.gz \
  --output sample.annotated.bam \
  --three-prime \
  --threads 8
samtools index -@ 8 sample.annotated.bam
```

For a single-sort streaming workflow, use `-` for annotator stdin/stdout and
provide an explicit metrics path:

```bash
sbxpress minimap ... \
| sbxpress annotate - --gtf genes.gtf.gz --output - \
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
minimap2's native MAPQ rather than applying Cell Ranger's MAPQ=255 promotion, so
the `mm` tag records the equivalent gene-rescue condition rather than an identical
STAR/Cell Ranger mapping action. This mode is intentionally sequential and cannot
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
sbxpress molecules sample.annotated.coordinate.bam \
  --output sample.molecules.coordinate.bam --threads 8
samtools index -@ 8 sample.molecules.coordinate.bam
```

Eligible reads receive `xf=17`; the selected representative receives `xf=25`
(Cell Ranger bit 8). `XV` records the molecule-level gene collision degree and
`XW` records `U`nique, `W`inner, `L`oser, or unresolved `T`ie. Ties are not counted.

Audit a molecule BAM with comparison-ready mapping, annotation, collision, UMI,
and sensitivity metrics:

```bash
sbxpress audit sample.molecules.coordinate.bam \
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

The three coordinate-sliced Cell Ranger fixtures can be checked with:

```bash
./tests/validate_cr_bam_to_mtx.sh /path/to/test_fixtures
```

See **[INSTALL.md](INSTALL.md)** for the full operations manual: system requirements,
getting the source onto a server, the barcode whitelist, throughput tuning, and
troubleshooting.

### Barcode whitelist (required)

This is **Single Cell Multiome (ATAC + GEX)** data, GEX side, so the cell barcodes come
from the **ARC** gel-bead list **`737K-arc-v1.txt.gz`** (736,320 barcodes). It is 10x
reference data staged under the gitignored `resources/cellranger_barcodes/`
directory. Obtain it from a Cell Ranger ARC distribution:

```
cellranger-arc-<version>/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
```

and pass it with `--whitelist` (or place it in `resources/cellranger_barcodes/` before
running `install.sh`).

For stable feature rows across samples, build one canonical `genes.mtx` from the
reference GTF with `sbxpress annotate --genes-mtx`, then pass that same file to
every `sbxpress matrix --features` invocation. Existing manifests are verified
against the filtered annotation rather than overwritten.

---

## Usage

For the end-to-end genomic sbxpress minimap reference and counting workflow, see
[Preparing a sbxpress minimap reference for SBX 3′ GEX counting](REFERENCE_PREPARATION.md).

Single directory of run files, all cores, ready for Cell Ranger:

```bash
sbxpress extract \
  --fastqs   /data/CAT02/raw_fastqs_run1/device_0/out \
  --sample-sheet /data/CE15204_bclconvert_input.csv \
  --whitelist    /usr/local/share/sbxpress/737K-arc-v1.txt.gz \
  --out-dir  /data/out/run1 \
  --threads  0 \
  --cellranger
```

Explicit files across multiple devices:

```bash
sbxpress extract --fastqs /data/.../device_{0,1,2,3}/out/PrimaryAnalysis_*.fastq.gz \
  --sample-sheet sheet.csv --whitelist 737K-arc-v1.txt.gz -o /data/out/run1 --cellranger
```

Process each SBX run separately (the two runs are the same library in opposite
orientations; orientation is auto-detected and enforced per file).

### Options

| Flag | Default | Purpose |
|---|---|---|
| `--fastqs` | (required) | Input SBX FASTQs: a directory of `PrimaryAnalysis_*.fastq.gz`, or explicit files. |
| `--sample-sheet` | — | BCL Convert v2 CSV (primary) or a flat `sample,i7,i5` CSV. |
| `--whitelist` | — | 10x cell-barcode whitelist (`737K-arc-v1.txt.gz`). |
| `--out-dir` / `-o` | `sbxpress_out` | Output directory. |
| `--threads` | 0 (all cores) | Worker pool size. |
| `--cellranger` | off | **Cell Ranger compatibility preset** — applies its required settings (currently `--max-r2-len 649`). Use this for Cell Ranger input. |
| `--minimap2` | off | Emit only R2 FASTQs. Adds Cell Ranger barcode/UMI SAM tags to each FASTQ comment; align with sbxpress minimap `-y`. Cannot be combined with `--cellranger`. |
| `--processing-tags` | off | Add SBX processing-provenance tags (`XO`, `XB`, `XI`, `XA`, `XL`) to R2 comments. Implied by `--minimap2`. |
| `--orientation` | auto | Force `p5-led` / `p7-led` if auto-detect is wrong. |
| `--cb-gate` | on | Drop pairs whose CB is neither on-whitelist nor within Hamming-1. |
| `--cb-reframe` | on | Rescue homopolymer-indel frame-shifted barcodes. |
| `--index-window` | 3 | Offset-tolerant i5/i7 matching (absorbs indels shifting the index). |
| `--index-max-mm` | 1 | Max mismatches when matching sample-sheet indices. |
| `--realign` | on | Single-indel CB rescue — recovers ~+1% of reads at high precision (`--realign false` to disable). |
| `--lenient` | off | Within-2 CB rescue — leave off (~51% precise on the dense ARC whitelist). |
| `--min-cdna-len` | 50 | Shortest cDNA emitted as R2. |
| `--max-r2-len` | 0 (uncapped) | Cap emitted R2 length (seq+qual together). Cell Ranger hard-caps at 649 bp; `--cellranger` sets this. |
| `--gzip-level` | 4 | Output gzip compression level (lower = faster/larger). |
| `--split` | 1 | Split each sample across N gzip files (`_001`..`_00N`) so N threads compress one sample in parallel. |
| `--reader-threads` | 0 (auto) | I/O parallelism. Auto uses every core even for a single input file (see below); `N>0` forces the pipelined path with N readers. |

Run `sbxpress extract --help` for extraction options. Set `RUST_LOG=debug` for verbose logging.

### Cell Ranger compatibility

Two quirks of Cell Ranger are handled for you:

- **Hard 649 bp read cap.** A read ≥ 650 bp fails with STAR `quality string length
  is not equal to sequence length`. `--cellranger` caps R2 to 649 bp — applied at
  write time, after the full-length fragment has been adapter-trimmed and
  reverse-complemented, so the retained bases are the ones the R2 primer would read.
- **Single-member gzip.** Cell Ranger mis-frames records at multi-member gzip
  boundaries, so each FASTQ is written as a single gzip member (compressed once,
  during processing).

### Minimap2-ready output

Use `--minimap2` when R2 will be aligned directly with minimap2. It emits only
`*_R2_*.fastq.gz`; there is deliberately no R1 file. Each R2 FASTQ comment contains
the following SAM-style tags, which sbxpress minimap preserves when invoked with `-y`:

- `CR` / `UR`: raw cell barcode / raw UMI.
- `CB`: whitelist-corrected cell barcode, with the Cell Ranger `-1` gem-group suffix.
- `UB`: provisional corrected UMI (the extracted UMI before gene-aware UMI correction).
- `CY` / `UY`: raw barcode and UMI Phred qualities.
- `XO`: original SBX read orientation (`5` = P5-led, `7` = P7-led).
- `XB`: barcode route: 0 exact, 1 Hamming-1, 2 indel realignment, 3 reframe,
  4 Hamming-2, 5 no whitelist, 6 unrescued miss emitted with gating disabled.
- `XI`: index route: 0 dual, 1 i5-only, 2 i7-only, 3 conflict, 4 no usable
  index, 5 demultiplexing not requested.
- `XA`: anchor bitmask: 1 R1-primer tail, 2 poly-T, 4 cDNA adapter trimmed,
  8 file-level expected/detected orientation applied.
- `XL`: extracted cDNA length before `--max-r2-len` capping.

The five SBX tags are also available in paired-FASTQ mode with
`--processing-tags`; they are placed on R2, the record that is aligned.

For example:

```bash
sbxpress extract --fastqs /data/sbx --whitelist 737K-arc-v1.txt.gz -o /data/out --minimap2
sbxpress minimap -a -x splice -uf -k14 --secondary=no -y reference.mmi /data/out/*_R2_*.fastq.gz \
  | samtools sort -o aligned.bam
```

Do not create a combined barcode+UMI tag at this stage. Group records by `(CB, UR)`
for UMI correction, then by `(CB, UB, gene)` after gene assignment; an `MI` tag, if
needed, should be generated only at that downstream molecule-calling step.

### Illumina input for matched-platform comparisons

`sbxpress illumina` converts standard paired 10x Illumina FASTQs into the same
R2-only minimap2 handoff used by `sbxpress extract --minimap2`. R1 supplies the
16 bp cell barcode and 12 bp UMI; R2 is emitted unchanged because Illumina R2 is
already in mapping orientation. Exact and uniquely resolvable Hamming-1 cell
barcodes are retained, while invalid, absent, and ambiguous barcodes are gated.

Repeat `--r1` and `--r2` in matching order for multiple lanes or shards. File
pairs are converted concurrently (`--threads` limits the lane workers) and a
separate numbered output shard is created for each pair:

```bash
sbxpress illumina \
  --r1 sample_L001_R1_001.fastq.gz --r2 sample_L001_R2_001.fastq.gz \
  --r1 sample_L002_R1_001.fastq.gz --r2 sample_L002_R2_001.fastq.gz \
  --whitelist resources/cellranger_barcodes/3M-february-2018.txt.gz \
  --output sample.illumina.minimap.fastq.gz \
  --threads 2

sbxpress minimap -a -x splice:sr --secondary=no -y reference.mmi \
  sample.illumina.minimap.shard-*.fastq.gz \
  | sbxpress annotate - --gtf genes.gtf --output -
```

For one input pair, the requested output filename is retained. For multiple
pairs, outputs are named `*.shard-001.fastq.gz`, `*.shard-002.fastq.gz`, and so
on. Their comments carry `CR`, `CB`, `UR`, `UB`, `CY`, and `UY`. A combined JSON
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
  CB-gate stats, demux breakdown, per-sample counts.
- `sample_sheet.csv` — the normalized flat sample sheet.

Every run also prints a **RUN SUMMARY** to stderr (elapsed, threads, reads/sec,
`processed / emitted / assigned / lost`, and a lost breakdown). Capture it with
`2> run1.log`.

---

## Architecture

Library crate (`src/`), each module a pipeline stage:

| Module | Role |
|---|---|
| `orient.rs` | Per-read orientation → canonical P5→P7 frame. |
| `homopolymer.rs` | poly(dT)/poly(A) run detection (the barcode/cDNA boundary). |
| `extract.rs` | R1/R2 bracketing + 3′ adapter trimming (indel-tolerant). |
| `barcode.rs` | Whitelist load + gate/rescue of the cell barcode. |
| `demux.rs` | i5/i7 extraction + offset-tolerant sample-sheet matching. |
| `sample_sheet.rs` | Sample-sheet parsing (BCL Convert v2 / flat CSV). |
| `anchors.rs` | Adapter / primer / leader constants — single source of truth. |
| `io.rs` | FASTQ input discovery and single-member gzip output writers. |
| `accounting.rs` | Metrics aggregation → `metrics_summary.json`. |
| `pipeline.rs` | Core per-read types and state. |
| `config.rs` | CLI definition and runtime parameters. |
| `lib.rs` | Orchestration: parallel per-file processing and emission. |

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

Output stays mate-ordered (the Nth R1 record is always the mate of the Nth R2) in
both modes; the pipelined path does not preserve input record order across the file.

---

## Development

```bash
cargo build --release     # optimized binary
cargo test  --release     # unit tests
```

The `examples/` directory holds diagnostic/analysis tools (per-read disposition,
loss taxonomy, CB-distance and index audits, etc.), run with:

```bash
cargo run --release --example disposition -- <fastq> <whitelist> <sheet> <p5-led|p7-led> [max]
```

---

## License

See [LICENSE](LICENSE).
