# sbx10x — Installation & Operations Manual (Ubuntu)

`sbx10x` converts raw Roche SBX long reads (one full library molecule per read,
P5→P7) into per-sample 10x R1/R2 FASTQs ready for Cell Ranger ARC. This guide
covers installing and running it on an Ubuntu server.

---

## 1. System requirements

| | Recommended |
|---|---|
| OS | Ubuntu 20.04 LTS or newer (any Debian-family distro works) |
| CPU | 8+ cores; scales across files (a 64–128-core node clears a full run in minutes) |
| RAM | 4 GB minimum; ~0.5 GB peak per ~14 worker threads. The barcode whitelist adds ~50 MB |
| Disk | Free space ≥ ~1× the size of your input FASTQs for the output |
| Toolchain | Rust 1.74+ (the installer fetches the latest stable via rustup if absent) |
| Build deps | A C compiler (`build-essential`) — used to compile the bundled zlib for output gzip. **No cmake needed** (zlib1g-dev optional; the build vendors zlib if absent) |

Throughput reference (Apple M4 Pro, 14 cores): ~142K reads/s single-thread,
~1.3M reads/s at 14 threads → **~13 min per billion reads**. Server CPUs with
more physical cores scale further.

---

## 2. Quick install

From the `sbx10x` source directory:

```bash
./install.sh
```

This will:

1. `apt-get install` build prerequisites (`build-essential curl pkg-config git`).
2. Install the Rust toolchain via rustup if `cargo` is not already on `PATH`.
3. Build the optimized release binary (`cargo build --release --locked`).
4. Install the binary to `/usr/local/bin/sbx10x`.
5. Install the ARC whitelist to `/usr/local/share/sbx10x/` **if it is present**
   in `resources/cellranger_barcodes/` (see §4 — it is not committed to git).
6. Run a built-in smoke test.

Common variants:

```bash
./install.sh --prefix "$HOME/.local"          # user-local install (no sudo)
./install.sh --whitelist /data/737K-arc-v1.txt.gz
./install.sh --skip-apt                         # packages already present / no sudo
./install.sh --skip-rust                        # use the cargo already installed
```

If you installed to `$HOME/.local`, ensure it is on your `PATH`:

```bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
```

Verify:

```bash
sbx10x --version
sbx10x --help
```

---

## 3. Getting the source onto the server

There is **no public git remote** for this project, so use one of:

```bash
# A) rsync the working tree (recommended — includes the whitelist if you have it)
rsync -av --exclude target/ /path/to/sbx10x/ user@server:~/sbx10x/

# B) tarball
tar czf sbx10x.tgz --exclude=target sbx10x/ && scp sbx10x.tgz user@server:~/

# C) if you push it to a git host first
git clone <your-remote-url> sbx10x      # NOTE: clone does NOT include the whitelist
```

`rsync`/`tar` bring the bundled whitelist along; `git clone` does **not** (§4).

---

## 4. The barcode whitelist (required for the cell-barcode gate)

This data is **10x Single Cell Multiome (ATAC + GEX), GEX side**, so the cell
barcodes come from the **ARC** gel-bead list **`737K-arc-v1.txt.gz`** (736,320
barcodes), *not* any 3′ GEX list.

This file is **not committed to the repository** (it is external Cell Ranger
reference data, and `resources/cellranger_barcodes/*.gz` is gitignored). Obtain
it from a Cell Ranger ARC distribution:

```
cellranger-arc-<version>/lib/python/cellranger/barcodes/737K-arc-v1.txt.gz
```

Then either:

- place it at `resources/cellranger_barcodes/737K-arc-v1.txt.gz` **before**
  running `install.sh` (it will be installed to `share/sbx10x/`), or
- pass it directly at runtime with `--whitelist /path/to/737K-arc-v1.txt.gz`.

Without a whitelist the pipeline still runs but cannot gate/rescue cell barcodes
(quality drops sharply) — always supply it for production.

---

## 5. Inputs

1. **Raw SBX FASTQs** — `--fastqs`. Either a directory of
   `PrimaryAnalysis_*.fastq.gz` or an explicit list of files. Plain gzip is
   expected (each file inflates serially; parallelism is across files).
2. **Sample sheet** — `--sample-sheet`. The BCL Convert v2 CSV
   (e.g. `CE15204_bclconvert_input.csv`) is the primary format; a flat
   `sample,i7,i5` CSV is the fallback. A flat CSV is also written to the output.
3. **Whitelist** — `--whitelist 737K-arc-v1.txt.gz` (see §4).

---

## 6. Running

Single directory of run files, all cores:

```bash
sbx10x \
  --fastqs /data/CAT02/raw_fastqs_run1/device_0/out \
  --sample-sheet /data/CE15204_bclconvert_input.csv \
  --whitelist /usr/local/share/sbx10x/737K-arc-v1.txt.gz \
  --out-dir /data/out/run1 \
  --threads 0          # 0 = use all cores
```

Explicit files / multiple devices:

```bash
sbx10x --fastqs /data/.../device_{0,1,2,3}/out/PrimaryAnalysis_*.fastq.gz \
  --sample-sheet sheet.csv --whitelist 737K-arc-v1.txt.gz -o /data/out/run1
```

Process each SBX run separately (the two runs are the same library in opposite
orientations; orientation is auto-detected and enforced per file). Useful flags:

| Flag | Default | Purpose |
|---|---|---|
| `--threads` | 0 (all cores) | worker pool size |
| `--orientation` | auto | force `p5-led` / `p7-led` if auto-detect is wrong |
| `--cb-reframe` | on | rescue homopolymer-indel frame-shifted barcodes |
| `--index-window` | 3 | offset-tolerant i5/i7 matching |
| `--realign` / `--lenient` | off | extra opt-in CB rescues (small marginal gain) |
| `--cb-gate` | on | drop pairs whose CB isn't on/near the whitelist |
| `--cellranger` | off | Cell Ranger compatibility preset — applies the settings it requires (currently `--max-r2-len 649`; explicit flags override). **Use this for Cell Ranger input.** |
| `--max-r2-len` | 0 (uncapped) | cap R2 (cDNA) length, seq+qual together; 0 = full-length. Cell Ranger hard-caps at 649 bp (`--cellranger` sets this) |
| `--gzip-level` | 4 | output gzip compression level (lower = faster/larger); output compression is the wall-clock floor at scale |
| `--split` | 1 | split each sample across N gzip files (`_001`..`_00N`); N>1 parallelizes per-sample compression |

Tune logging with `RUST_LOG` (e.g. `RUST_LOG=debug sbx10x ...`).

---

## 7. Outputs

In `--out-dir`:

- `<sample>_S<n>_L001_R1_001.fastq.gz` / `_R2_001.fastq.gz` — per-sample pairs in
  Cell Ranger naming (R1 = 28 bp CB+UMI, R2 = cDNA).
- `Undetermined_S0_L001_R{1,2}_001.fastq.gz` — pairs with no matchable index.
- `metrics_summary.json` — full run accounting: `total_reads`, `elapsed_secs`,
  `reads_per_sec`, `threads`, read classes, CB-gate stats, demux breakdown, and
  per-sample counts.
- `sample_sheet.csv` — the flat sample sheet derived from the input.

At the end of every run the console (stderr) prints a **RUN SUMMARY** —
elapsed time, threads, reads/sec, and `processed / emitted / assigned / lost`,
followed by a `lost breakdown` (unparseable · no_cdna · cb_fail_gated ·
index_undetermined). Capture it with `2> run1.log` (or `|& tee run1.log`).

These feed directly into `cellranger-arc count` (use the ARC reference). Each FASTQ
is written as a **single gzip member** (compressed once, during processing, via the
zlib backend). This matters: Cell Ranger mis-frames records at multi-member gzip
boundaries, so single-member output is required — it is always produced.

**Long R2 + Cell Ranger:** Cell Ranger has a **hard 649 bp read-length cap** — a read of
650 bp or longer fails with STAR `quality string length is not equal to sequence length`.
SBX cDNA routinely exceeds this. Pass **`--cellranger`** to apply the Cell Ranger
compatibility preset (currently `--max-r2-len 649`), which caps R2 to 649 bp (sequence +
quality together; read counts unchanged). sbx10x also trims the 3′ read-through adapter
itself — including SBX homopolymer-indel variants (`AAGAGC`→`AAAGAGC`) and 3′-end
partials — so the emitted cDNA is adapter-free; the cap simply keeps every read inside
Cell Ranger's length limit. STAR uses ~90 bp for gene mapping regardless of R2 length.

**Throughput tip (`--split`):** compression for one sample goes through one writer,
so a hot sample can bottleneck on a single core. `--split N` writes that sample as N
files (`…_R1_001.fastq.gz` … `_00N`), which Cell Ranger reads together, letting N
threads compress it in parallel. Set N up to `--threads` and keep
`N × (samples+1) × 2` under your open-file limit (`ulimit -n`); e.g. `--split 8`.
Within each file R1/R2 stay mate-ordered.

---

## 8. Troubleshooting

| Symptom | Fix |
|---|---|
| `cargo: command not found` after install | `source "$HOME/.cargo/env"` or open a new shell |
| Permission denied installing to `/usr/local` | re-run with `sudo`, or `--prefix "$HOME/.local"` |
| `sbx10x: command not found` | add the prefix `bin/` to `PATH` (§2) |
| Very low yield / everything Undetermined | missing/wrong `--whitelist`, or wrong `--orientation` |
| Build fails on old Ubuntu | update Rust: `rustup update stable` |
| `is cmake not installed?` / `libz-sys` build error | you're on an old checkout — update the source (the `gzp`/zlib-ng dependency was removed; gzip is now pure-Rust and needs no cmake). As a one-off workaround: `sudo apt-get install -y cmake` |
| OOM with many threads | lower `--threads` (peak RAM ≈ threads × targets × 256 KB + whitelist) |

Validate a build at any time:

```bash
cargo test --release      # unit tests
./install.sh              # re-run; includes the smoke test
```
