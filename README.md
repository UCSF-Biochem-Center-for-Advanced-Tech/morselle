# MORSELLE

**MOlecular Read Structure, Extraction, Labeling & Enumeration**

MORSELLE turns raw Roche SBX long reads, or paired Illumina 10x reads, into
per-sample 10x-style FASTQs, annotated BAMs, corrected molecule counts, and
raw and filtered feature-barcode matrices. It supports 10x 3′ GEX, 10x Flex
(fixed RNA profiling, including multiplexed probe barcodes and CRISPR guide
capture), and pooled CAT Tag-seq.

> [!WARNING]
> **MORSELLE is alpha software.** Its commands, options, output formats, and
> default settings may change between releases without notice. Results have
> been checked against Cell Ranger on a limited set of in-house SBX and
> Illumina libraries only.
>
> **If you plan to use MORSELLE for its current capabilities, please contact
> Matt Settles at [matt.settles@ucsf.edu](mailto:matt.settles@ucsf.edu) before
> you start** for help choosing settings and interpreting results, and to hear
> about known limitations for your assay.

---

## How it works

Roche SBX sequences each library molecule end-to-end as one long read (P5→P7):

```text
leader – P5 – i5 – R1primer – CB(16) – UMI(12) – polyT – cDNA – adapter – i7 – rc(P7)
```

MORSELLE reads that structure from every read and carries the barcodes through
alignment to counts:

1. **Extract** (`morselle extract`): orient each SBX read, demultiplex by i5/i7,
   gate and rescue the cell barcode against the 10x whitelist, recover the UMI,
   and trim the cDNA. Illumina reads take the same path through
   `morselle illumina`.
2. **Align**: minimap2 (`morselle minimap`) splice-aligns the tagged cDNA to
   the genome (or `morselle probeset-align` for Flex probe sets).
3. **Label** (`morselle annotate`): assign genes, transcripts, and regions from
   a GENCODE GTF and add 10x-compatible `GX`/`GN`/`TX`/`RE` tags.
4. **Enumerate** (`morselle molecules`, `matrix`, `cells`): correct UMIs,
   resolve gene collisions, count one representative per molecule, write a
   10x-style raw matrix, and call cells (OrdMag + EmptyDrops).

| Command | Purpose |
|---|---|
| `morselle extract` | Raw SBX reads → demultiplexed, tagged cDNA FASTQs |
| `morselle illumina` | Paired Illumina 10x FASTQs → the same tagged cDNA FASTQ |
| `morselle minimap` | minimap2, with arguments and streams passed through |
| `morselle annotate` | Gene / transcript / region / collision tags on a BAM |
| `morselle molecules` | UMI correction, collision resolution, representative marking |
| `morselle matrix` | 10x-compatible raw feature-barcode matrix (MEX + H5) |
| `morselle cells` | Cell calling and filtered matrices |
| `morselle audit` | Mapping, annotation, UMI, and molecule metrics for comparisons |
| `morselle probe-demux` | Flex: assign cell barcodes to samples by probe barcode |
| `morselle probeset-align` | Flex: assign GEX reads to genes with the probe-set reference |
| `morselle guides` | Feature Barcoding reads (CRISPR guides, antibody tags) → features |
| `morselle protospacer-call` | CRISPR guide identity per cell |

Run `morselle <COMMAND> --help` for each command's options.

---

## Quick start

### 1. Install

On Ubuntu/Debian or macOS, from a clone of this repository:

```bash
git clone https://github.com/UCSF-Biochem-Center-for-Advanced-Tech/morselle.git
cd morselle
./install.sh --prefix "$HOME/.local"
```

`install.sh` installs the build prerequisites, the Rust toolchain (if absent),
`minimap2`, `samtools`, and `k8`/`paftools.js`, then builds and installs
`morselle`. Omit `--prefix` to install under `/usr/local` with sudo. Check the
installation:

```bash
morselle --version
morselle extract --help
```

See [docs/INSTALL.md](docs/INSTALL.md) for requirements, uninstalling, and
troubleshooting.

### 2. Prepare the reference inputs

You need:

- a 10x cell-barcode whitelist that matches the chemistry (for example
  `737K-arc-v1.txt.gz` for Multiome GEX), from a Cell Ranger installation. It is
  not redistributed here; see [resources/10x_barcodes](resources/10x_barcodes/README.md);
- a minimap2 splice index (`.mmi`), the matching GENCODE GTF, and a junction
  BED. [docs/REFERENCE_PREPARATION.md](docs/REFERENCE_PREPARATION.md) builds
  all three.

### 3. Run the pipeline

`scripts/morselle-pipeline.sh` runs every stage for one sample and can be
restarted. Raw SBX 3′ GEX:

```bash
scripts/morselle-pipeline.sh \
  --input-mode sbx --analysis-mode gex \
  --fastqs raw/device_0/out raw/device_1/out \
  --bclconvert-samplesheet samples.csv --sample M2_GEX \
  --whitelist 737K-arc-v1.txt.gz \
  --mmi GRCh38.splice-sr.mmi \
  --gtf gencode.v44.primary_assembly.basic.annotation.gtf.gz \
  --junctions gencode.v44.primary_assembly.junctions.bed \
  --out-dir results/M2_GEX --threads 16
```

The same library sequenced on Illumina:

```bash
scripts/morselle-pipeline.sh \
  --input-mode illumina --analysis-mode gex \
  --r1 M2_GEX_S8_L003_R1_001.fastq.gz --r2 M2_GEX_S8_L003_R2_001.fastq.gz \
  --sample M2_GEX \
  --whitelist 737K-arc-v1.txt.gz \
  --mmi GRCh38.splice-sr.mmi \
  --gtf gencode.v44.primary_assembly.basic.annotation.gtf.gz \
  --junctions gencode.v44.primary_assembly.junctions.bed \
  --out-dir results/M2_GEX_illumina --threads 16
```

Use `--analysis-mode flex` (with `--probe-set`) or `--analysis-mode tagseq` for
those assays. The pipeline writes a familiar layout:

```text
results/M2_GEX/M2_GEX/outs/
  raw_feature_bc_matrix/          # MEX: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
  raw_feature_bc_matrix.h5
  filtered_feature_bc_matrix/     # called cells
  filtered_feature_bc_matrix.h5
```

---

## Documentation

| Document | Contents |
|---|---|
| [docs/INSTALL.md](docs/INSTALL.md) | Installation, the barcode whitelist, running on a server, troubleshooting |
| [docs/MANUAL.md](docs/MANUAL.md) | End-to-end workflows: SBX and Illumina GEX, 10x Flex, CAT Tag-seq, audits |
| [docs/USAGE.md](docs/USAGE.md) | Command reference: every stage, `extract` options, SAM tags, outputs |
| [docs/REFERENCE_PREPARATION.md](docs/REFERENCE_PREPARATION.md) | Building the minimap2 reference, junction BED, and feature manifest |
| [docs/PROFILING.md](docs/PROFILING.md) | Profiling and benchmarking for developers |
| [CHANGELOG.md](CHANGELOG.md) | Release history |

---

## Development

Building from source needs the Rust toolchain and HDF5. The contributor checks
also need Python 3 with `numpy`.

```bash
cargo build --release     # binaries in target/release/
scripts/check.sh          # what CI runs: build, tests, test inventory, fmt, clippy
```

`tests/golden_outputs.rs` pins the exact output of `molecules` and `matrix` in
every mode, so a refactor that changes counts fails the build. The `examples/`
directory holds diagnostic tools (per-read disposition, loss taxonomy,
barcode-distance and index audits):

```bash
cargo run --release --example disposition -- <fastq> <whitelist> <sheet> <p5-led|p7-led> [max]
```

---

## Contact

MORSELLE is developed at the UCSF Biochemistry Center for Advanced Technology.
Questions, collaborations, and support requests go to Matt Settles at
[matt.settles@ucsf.edu](mailto:matt.settles@ucsf.edu). Please report bugs as
GitHub issues.

## License

BSD 3-Clause. See [LICENSE](LICENSE).
