# sbx10x

**Convert raw Roche SBX long reads into per-sample 10x R1/R2 FASTQs for Cell Ranger.**

Roche SBX sequences each library molecule end-to-end as one long read (P5→P7),
so a 10x Single Cell Multiome **GEX** library comes off the instrument as full
molecules rather than the paired short reads Cell Ranger expects. `sbx10x` parses
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

For every read, `sbx10x`:

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

From the source directory on an Ubuntu/Debian host:

```bash
./install.sh
```

This installs build prerequisites, the Rust toolchain (via rustup if absent), builds
the optimized binary, and installs it to `/usr/local/bin/sbx10x`. User-local install
without sudo:

```bash
./install.sh --prefix "$HOME/.local"
```

Verify:

```bash
sbx10x --version
sbx10x --help
```

See **[INSTALL.md](INSTALL.md)** for the full operations manual: system requirements,
getting the source onto a server, the barcode whitelist, throughput tuning, and
troubleshooting.

### Barcode whitelist (required)

This is **Single Cell Multiome (ATAC + GEX)** data, GEX side, so the cell barcodes come
from the **ARC** gel-bead list **`737K-arc-v1.txt.gz`** (736,320 barcodes). It is 10x
reference data and is **not committed** to this repo. Obtain it from a Cell Ranger ARC
distribution:

```
cellranger-arc-<version>/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
```

and pass it with `--whitelist` (or place it in `resources/cellranger_barcodes/` before
running `install.sh`).

---

## Usage

Single directory of run files, all cores, ready for Cell Ranger:

```bash
sbx10x \
  --fastqs   /data/CAT02/raw_fastqs_run1/device_0/out \
  --sample-sheet /data/CE15204_bclconvert_input.csv \
  --whitelist    /usr/local/share/sbx10x/737K-arc-v1.txt.gz \
  --out-dir  /data/out/run1 \
  --threads  0 \
  --cellranger
```

Explicit files across multiple devices:

```bash
sbx10x --fastqs /data/.../device_{0,1,2,3}/out/PrimaryAnalysis_*.fastq.gz \
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
| `--out-dir` / `-o` | `sbx10x_out` | Output directory. |
| `--threads` | 0 (all cores) | Worker pool size. |
| `--cellranger` | off | **Cell Ranger compatibility preset** — applies its required settings (currently `--max-r2-len 649`). Use this for Cell Ranger input. |
| `--orientation` | auto | Force `p5-led` / `p7-led` if auto-detect is wrong. |
| `--cb-gate` | on | Drop pairs whose CB is neither on-whitelist nor within Hamming-1. |
| `--cb-reframe` | on | Rescue homopolymer-indel frame-shifted barcodes. |
| `--index-window` | 3 | Offset-tolerant i5/i7 matching (absorbs indels shifting the index). |
| `--index-max-mm` | 1 | Max mismatches when matching sample-sheet indices. |
| `--realign` / `--lenient` | off | Extra opt-in CB rescues (small marginal gain). |
| `--min-cdna-len` | 50 | Shortest cDNA emitted as R2. |
| `--max-r2-len` | 0 (uncapped) | Cap emitted R2 length (seq+qual together). Cell Ranger hard-caps at 649 bp; `--cellranger` sets this. |
| `--gzip-level` | 4 | Output gzip compression level (lower = faster/larger). |
| `--split` | 1 | Split each sample across N gzip files (`_001`..`_00N`) so N threads compress one sample in parallel. |
| `--reader-threads` | 0 (auto) | I/O parallelism. Auto uses every core even for a single input file (see below); `N>0` forces the pipelined path with N readers. |

Run `sbx10x --help` for the complete list. Set `RUST_LOG=debug` for verbose logging.

### Cell Ranger compatibility

Two quirks of Cell Ranger are handled for you:

- **Hard 649 bp read cap.** A read ≥ 650 bp fails with STAR `quality string length
  is not equal to sequence length`. `--cellranger` caps R2 to 649 bp — applied at
  write time, after the full-length fragment has been adapter-trimmed and
  reverse-complemented, so the retained bases are the ones the R2 primer would read.
- **Single-member gzip.** Cell Ranger mis-frames records at multi-member gzip
  boundaries, so each FASTQ is written as a single gzip member (compressed once,
  during processing).

---

## Outputs

In `--out-dir`:

- `<sample>_S<n>_L001_R1_001.fastq.gz` / `_R2_001.fastq.gz` — per-sample pairs
  (R1 = 28 bp CB+UMI, R2 = cDNA).
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
