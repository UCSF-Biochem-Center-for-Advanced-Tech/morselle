# MORSELLE — Installation & Operations Manual (Ubuntu and macOS)

MORSELLE (MOlecular Read Structure, Extraction, Labeling & Enumeration)
converts raw Roche SBX long reads (one full library molecule per read, P5→P7)
and paired Illumina 10x reads into per-sample 10x-style FASTQs, annotated BAMs,
and feature-barcode matrices. This guide covers installing and running it on
Ubuntu/Debian and macOS.

MORSELLE is alpha software; contact Matt Settles
([matt.settles@ucsf.edu](mailto:matt.settles@ucsf.edu)) before using it for
your own data.

---

## 1. System requirements

| | Recommended |
|---|---|
| OS | Ubuntu 20.04+ / Debian-family Linux, or macOS with Homebrew |
| CPU | 8+ cores; scales across files (a 64–128-core node clears a full run in minutes) |
| RAM | 4 GB minimum; ~0.5 GB peak per ~14 worker threads. The barcode whitelist adds ~50 MB |
| Disk | Free space ≥ ~1× the size of your input FASTQs for the output |
| Toolchain | Rust 1.74+ (the installer fetches the latest stable via rustup if absent) |
| Build deps | Linux: `build-essential`, `pkg-config`, `cmake`, `clang`, `libclang-dev`; macOS: Xcode Command Line Tools plus Homebrew `cmake` and `pkg-config` |
| Alignment tools | `minimap2`, `samtools`, and `paftools.js` (with its `k8` runtime; installed automatically) |

Throughput reference (Apple M4 Pro, 14 cores): ~142K reads/s single-thread,
~1.3M reads/s at 14 threads → **~13 min per billion reads**. Server CPUs with
more physical cores scale further.

---

## 2. Quick install

From the `morselle` source directory:

```bash
./install.sh
```

This will:

1. Install build prerequisites, `minimap2`, and `samtools` using `apt-get` on
   Ubuntu/Debian or Homebrew on macOS, then install the official `k8` runtime and matching
   `paftools.js` minimap2 companion script.
2. Install the Rust toolchain via rustup if `cargo` is not already on `PATH`.
3. Build the optimized release binaries (`cargo build --release --locked`).
4. Install the unified `morselle` command to `/usr/local/bin/`; its processing
   engines are installed privately under `/usr/local/libexec/morselle/`.
5. Run a built-in smoke test.

Common variants:

```bash
./install.sh --prefix "$HOME/.local"          # user-local install (no sudo)
./install.sh --skip-packages                    # packages already present / no sudo
./install.sh --skip-apt                         # legacy alias for --skip-packages
./install.sh --skip-rust                        # use the cargo already installed
```

If you installed to `$HOME/.local`, ensure it is on your `PATH`:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
```

Verify:

```bash
morselle --version
morselle --help
morselle annotate --version
morselle annotate --help
morselle molecules --version
morselle molecules --help
morselle audit --version
morselle audit --help
morselle matrix --version
morselle matrix --help
morselle minimap --version
samtools --version
paftools.js
```

### Uninstall

Use the same prefix used during installation:

```bash
./uninstall.sh --dry-run
./uninstall.sh
./uninstall.sh --prefix "$HOME/.local"
```

This removes `morselle`, its private engines, and `share/morselle`. Because `k8`
and `paftools.js` may be shared with other minimap2 workflows, they are preserved
unless `--remove-companions` is supplied. Use `--keep-data` to retain the
installed whitelist.

### BAM gene/transcript annotation

`morselle annotate` adds 10x-compatible `GX`, `GN`, `TX`, `AN`, and `RE`
tags to a coordinate-sorted genomic BAM using the same GTF used to build the
morselle minimap reference. Coordinate order is preserved; index the output afterward.

```bash
morselle annotate sample.coordinate.bam \
  --gtf gencode.v50.primary_assembly.annotation.gtf.gz \
  --output sample.annotated.bam \
  --three-prime \
  --threads 8
samtools index -@ 8 sample.annotated.bam
```

The annotator also supports streaming BAM/SAM input and BAM output with `-`. This
allows annotation before the only coordinate sort:

```bash
morselle minimap -a -x splice -m 20 -s 30 --secondary=no -y reference.mmi reads.fastq.gz \
| morselle annotate - --gtf genes.gtf.gz --output - \
    --metrics sample.annotation.metrics.json --three-prime \
| samtools sort -@ 8 -o sample.annotated.coordinate.bam -
samtools index -@ 8 sample.annotated.coordinate.bam
```

Intronic gene-body assignment is enabled by default, matching current 10x reference
GEX behavior; use `--include-introns false` for exon-only gene assignment. The
10x GENCODE biotype/readthrough/PAR rules and Ensembl ID-version stripping are
also enabled by default; each can be disabled explicitly for custom references.
This stage does not correct UMIs or set `xf` molecule flags.

Follow annotation with `morselle molecules`. It performs directional Hamming-1
UMI correction within each cell/gene, resolves cross-gene corrected-UMI
collisions by read support, writes corrected `UB`, and sets `xf=25` on exactly
one representative of each counted molecule (`xf=17` on other eligible reads).

```bash
morselle molecules sample.annotated.bam --output sample.molecules.bam --threads 8
samtools index -@ 8 sample.molecules.bam
```

The standard tags preserve all compatible assignments. Additional local-use tags
`XC` and `XT` always record compatible-gene and compatible-transcript counts. For
3' assays, pass `--three-prime` to additionally emit `ZG`/`ZT` terminal-supported
gene/transcript subsets and `XP`, the terminal-supported gene count. Those tags
and their JSON resolution metrics are not emitted for other assay types.

### Convert matched Illumina FASTQs

For a platform-controlled comparison, convert standard 10x Illumina R1/R2 files
to the same tagged R2-only handoff produced by `morselle extract --minimap2`:

```bash
morselle illumina \
  --r1 sample_L001_R1_001.fastq.gz \
  --r2 sample_L001_R2_001.fastq.gz \
  --whitelist /usr/local/share/morselle/barcodes/3M-february-2018.txt.gz \
  --output sample.illumina.minimap.fastq.gz
```

Repeat `--r1` and `--r2` in matching order for additional lanes. The converter
checks mate names and record counts, reads the 16 bp barcode plus 12 bp UMI from
R1, preserves R2 orientation, and emits `CR`, `CB`, `UR`, `UB`, `CY`, and `UY`.
Its JSON metrics default to `<output>.metrics.json`.

### BAM to raw gene-expression matrix

The second installed app reconstructs 10x reference molecule counts from a complete
10x reference BAM. For an exactly shaped raw matrix, provide the matching template
files so zero-count features/barcodes and original ordering are retained:

```bash
morselle matrix /data/outs/possorted_genome_bam.bam \
  --features /data/outs/raw_feature_bc_matrix/features.tsv.gz \
  --cell-barcodes /usr/local/share/morselle/barcodes/737K-arc-v1.txt.gz \
  --genome GRCh38 \
  --output /data/reconstructed/outs
```

The output follows 10x reference's layout: the three compressed Matrix Market files
are under `outs/raw_feature_bc_matrix/`, with `outs/raw_feature_bc_matrix.h5` as
their HDF5 neighbor and `outs/metrics_summary.json` for run accounting. The app
prints a `RUN SUMMARY` with elapsed time, throughput, records, molecules, features,
barcodes, and nonzero entries. Complete, matching BAM/matrix pairs should agree
exactly. Coordinate-sliced BAMs can under-count boundary genes when their
molecule-representative alignment lies outside the retained interval.

---

## 3. Getting the source onto the server

Clone a tagged release from GitHub, or transfer a checkout from your workstation:

```bash
# A) clone the v0.8.2 release (recommended)
git clone --branch v0.8.2 https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle.git
# NOTE: clone does NOT include the whitelist

# B) rsync a working tree (includes a locally staged whitelist)
rsync -av --exclude target/ /path/to/morselle/ user@server:~/morselle/

# C) tarball
tar czf morselle.tgz --exclude=target morselle/ && scp morselle.tgz user@server:~/
```

`rsync`/`tar` bring a locally staged whitelist along; `git clone` does **not** (§4).

---

## 4. The barcode whitelist (required for the cell-barcode gate)

This data is **10x Single Cell Multiome (ATAC + GEX), GEX side**, so the cell
barcodes come from the **ARC** gel-bead list **`737K-arc-v1.txt.gz`** (736,320
barcodes), *not* any 3′ GEX list.

These files are staged under `resources/10x_barcodes/` but remain ignored
by Git because they are external 10x reference reference data. Obtain
it from a 10x reference ARC distribution:

```
$TENX_INSTALL/lib/python/barcodes/737K-arc-v1.txt.gz
```

Pass it directly at runtime with `--whitelist /path/to/737K-arc-v1.txt.gz`.

Without a whitelist the pipeline still runs but cannot gate/rescue cell barcodes
(quality drops sharply) — always supply it for production.

---

## 5. Inputs

1. **Raw SBX FASTQs** — `--fastqs`. Either a directory of
   `PrimaryAnalysis_*.fastq.gz` or an explicit list of files. Plain gzip is
   expected (each file inflates serially; parallelism is across files).
2. **Sample sheet** — `--bclconvert-samplesheet`. The BCL Convert v2 CSV
   (e.g. `CE15204_bclconvert_input.csv`) is the primary format; a flat
   `sample,i7,i5` CSV is the fallback. A flat CSV is also written to the output.
3. **Whitelist** — `--whitelist 737K-arc-v1.txt.gz` (see §4).

---

## 6. Running

Single directory of run files, all cores:

```bash
morselle extract \
  --fastqs /data/CAT02/raw_fastqs_run1/device_0/out \
  --bclconvert-samplesheet /data/CE15204_bclconvert_input.csv \
  --whitelist /usr/local/share/morselle/737K-arc-v1.txt.gz \
  --out-dir /data/out/run1 \
  --threads 0          # 0 = use all cores
```

Explicit files / multiple devices:

```bash
morselle extract --fastqs /data/.../device_{0,1,2,3}/out/PrimaryAnalysis_*.fastq.gz \
  --bclconvert-samplesheet sheet.csv --whitelist 737K-arc-v1.txt.gz -o /data/out/run1
```

Process each SBX run separately (the two runs are the same library in opposite
orientations; orientation is auto-detected and enforced per file). Useful flags:

| Flag | Default | Purpose |
|---|---|---|
| `--threads` | 0 (all cores) | worker pool size |
| `--orientation` | auto | force `p5-led` / `p7-led` if auto-detect is wrong |
| `--cb-reframe` | off | rescue homopolymer-indel frame-shifted barcodes (about half right on CAT02; off by default) |
| `--r1-tail-relock` | on | rescue barcodes mislocated by an indel in the R1-primer tail |
| `--index-window` | 3 | offset-tolerant i5/i7 matching |
| `--index-indels` | on | also match an i5/i7 with one insertion or deletion (sheet-checked for safety) |
| `--demux-index` | auto | auto = per read: one index when it names exactly one sample, otherwise both must agree |
| `--realign` | on | single-indel CB rescue (recovers ~+1% of reads at high precision); disable with `--realign false` |
| `--cb-insertion` | on | rescue barcodes carrying one extra base (insertion) |
| `--lenient` | off | within-2 CB rescue — **leave off** (only ~51% precise on the dense ARC whitelist) |
| `--cb-gate` | on | drop pairs whose CB isn't on/near the whitelist |
| `--cellranger` | off | 10x reference compatibility preset — applies the settings it requires (currently `--max-r2-len 649`; explicit flags override). **Use this for 10x reference input.** |
| `--minimap2` | off | R2-only minimap2 preset. Puts `CR`, `CB`, `UR`, `UB`, `CY`, and `UY` SAM-style tags in the R2 FASTQ comment; use `morselle minimap` `-y`. Mutually exclusive with `--cellranger`. |
| `--processing-tags` | off | Add SBX provenance tags (`XO`, `XB`, `XI`, `XA`, `XL`) to R2 comments; automatically enabled by `--minimap2`. |
| `--max-r2-len` | 0 (uncapped) | cap R2 (cDNA) length, seq+qual together; 0 = full-length. 10x reference hard-caps at 649 bp (`--cellranger` sets this) |
| `--gzip-level` | 4 | output gzip compression level (lower = faster/larger); output compression is the wall-clock floor at scale |
| `--split` | 1 | split each sample across N gzip files (`_001`..`_00N`); N>1 parallelizes per-sample compression |
| `--reader-threads` | 0 (auto) | I/O strategy. Auto: with few input files (e.g. a single large one) decode is decoupled from processing so every core is used; with many files the classic one-file-per-thread path is used. `N>0` forces the pipelined path with N reader threads |

Tune logging with `RUST_LOG` (e.g. `RUST_LOG=debug morselle ...`).

---

## 7. Outputs

In `--out-dir`:

- `<sample>_S<n>_L001_R1_001.fastq.gz` / `_R2_001.fastq.gz` — per-sample pairs in
  10x reference naming (R1 = 28 bp CB+UMI, R2 = cDNA).
- With `--minimap2`, only the R2 files are written. Their comments carry raw and
  corrected barcode/UMI tags (`CR`/`CB`, `UR`/provisional `UB`) and qualities
  (`CY`, `UY`). Align with `morselle minimap -a -y`; `-y` copies these comments to SAM/BAM.
  `--minimap2` also adds SBX provenance tags for orientation (`XO`), barcode
  processing route (`XB`), index route (`XI`), anchor evidence (`XA` bitmask),
  and extracted pre-cap cDNA length (`XL`).
- `Undetermined_S0_L001_R{1,2}_001.fastq.gz` — pairs with no matchable index.
- `metrics_summary.json` — full run accounting: `total_reads`, `elapsed_secs`,
  `reads_per_sec`, `threads`, read classes, CB-gate stats, demux breakdown, and
  per-sample counts.
- `sample_sheet.csv` — the flat sample sheet derived from the input.

At the end of every run the console (stderr) prints a **RUN SUMMARY** —
elapsed time, threads, reads/sec, and `processed / emitted / assigned / lost`,
followed by a `lost breakdown` (unparseable · no_cdna · cb_fail_gated ·
index_undetermined). Capture it with `2> run1.log` (or `|& tee run1.log`).

These feed directly into the 10x reference count command (use the ARC reference). Each FASTQ
is written as a **single gzip member** (compressed once, during processing, via the
zlib backend). This matters: 10x reference mis-frames records at multi-member gzip
boundaries, so single-member output is required — it is always produced. (The
`--failed-output` file is the exception: it is multi-member gzip, which `zcat` and
other standard readers decode as one stream, and is not a 10x input.)

**Minimap2 workflow:** `--minimap2` is a distinct output mode, so it cannot be used
with `--cellranger`. It preserves the corrected cell barcode as `CB:Z:<barcode>-1`,
the raw barcode/UMI as `CR`/`UR`, and a provisional `UB` for later gene-aware UMI
correction. Keep barcode and UMI separate: group `(CB, UR)` before correction and
`(CB, UB, gene)` after assignment. A combined molecule identifier belongs downstream.

**Long R2 + 10x reference:** 10x reference has a **hard 649 bp read-length cap** — a read of
650 bp or longer fails with STAR `quality string length is not equal to sequence length`.
SBX cDNA routinely exceeds this. Pass **`--cellranger`** to apply the 10x reference
compatibility preset (currently `--max-r2-len 649`), which caps R2 to 649 bp (sequence +
quality together; read counts unchanged). morselle also trims the 3′ read-through adapter
itself — including SBX homopolymer-indel variants (`AAGAGC`→`AAAGAGC`) and 3′-end
partials — so the emitted cDNA is adapter-free; the cap simply keeps every read inside
10x reference's length limit. STAR uses ~90 bp for gene mapping regardless of R2 length.

**Throughput tip (`--split`):** compression for one sample goes through one writer,
so a hot sample can bottleneck on a single core. `--split N` writes that sample as N
files (`…_R1_001.fastq.gz` … `_00N`), which 10x reference reads together, letting N
threads compress it in parallel. Set N up to `--threads` and keep
`N × (samples+1) × 2` under your open-file limit (`ulimit -n`); e.g. `--split 8`.
Within each file R1/R2 stay mate-ordered.

---

## 8. Troubleshooting

| Symptom | Fix |
|---|---|
| `cargo: command not found` after install | `source "$HOME/.cargo/env"` or open a new shell |
| Permission denied installing to `/usr/local` | re-run with `sudo`, or `--prefix "$HOME/.local"` |
| `morselle: command not found` | add the prefix `bin/` to `PATH` (§2) |
| Very low yield / everything Undetermined | missing/wrong `--whitelist`, or wrong `--orientation` |
| Build fails on old Ubuntu | update Rust: `rustup update stable` |
| `is cmake not installed?` / `libz-sys` build error | you're on an old checkout — update the source (the `gzp`/zlib-ng dependency was removed; gzip is now pure-Rust and needs no cmake). As a one-off workaround: `sudo apt-get install -y cmake` |
| OOM with many threads | lower `--threads` (peak RAM ≈ threads × targets × 256 KB + whitelist) |

Validate a build at any time:

```bash
cargo test --release      # unit tests
./tests/validate_cr_bam_to_mtx.sh /path/to/test_fixtures
./install.sh              # re-run; includes the smoke test
```
