# v0.6.0 server acceptance: 0.5.6 vs 0.6.0 on real data

For 0.6.1 per-barcode metrics acceptance, compare each new
`*.barcodes.tsv.gz` table between the candidate builds after sorting rows,
and sum every count column against its stage's JSON metric. Test `illumina`,
`extract`, `annotate` and `molecules` on the data and execution modes listed
in `docs/PER_BARCODE_METRICS_PLAN.md`. The 0.6.0 baseline has no barcode
tables; compare its BAM, FASTQ and matrix outputs as before.

Server real-data acceptance was run with this procedure on 2026-10-08 and
passed (results in CHANGELOG 0.6.0). The notes below include what that run
showed needed changing.
This procedure runs the **0.5.6 baseline** and **0.6.0** on identical inputs
from three projects on the analysis server, compares every stage, and
explains every difference before 0.6.0 is accepted.

Purpose of each comparison:

- **Accuracy**: 0.6.0 counts agree with 0.5.6 where no change was intended,
  and with the independent Illumina/Cell Ranger results for the same library.
- **Completeness**: every assay (GEX, Flex, tag-seq) and input mode (SBX,
  Illumina) finishes and produces all outputs.
- **Differences**: every difference matches a changelog entry (table below);
  any other difference is a bug.
- **Bugs**: restart safety, indexed/sharded parity, and multi-sample identity.

## What "0.5.6" means here

Cargo reported 0.5.6 for eight commits. The baseline is **`4d42a18`**. It is
the last 0.5.6 commit before the fractional-default (`d408123`) and
modularization (`c5a64a7`) changes, and it matches the 0.5.6 behaviour that
`tests/golden_outputs.rs` was seeded from (`--tagseq` shared counting, no
fractional mode). Do not use 0.5.7 or 0.5.8: they are marked broken.

The `sbxpress` on the default PATH (`/usr/local/bin`) is **0.3.5**, and
`/mnt/analysis/software_shared/sbxpress/target/release` is a development
build. Neither may be used: `use_build` below checks this before every run.

0.5.6's own `scripts/sbxpress-pipeline.sh` does not parse (CHANGELOG 0.5.6).
Both builds therefore run under the **0.6.0 driver script**, with only
`PATH` switched. Every flag that driver passes by default exists in 0.5.6
(`--compression-level`, `--position-aware-umis` and `--fractional-collisions`
are passed only when requested and are not used for the comparison runs).
Differences can then come only from the binaries.

## Expected differences (0.5.6 → 0.6.0)

Unless this table lists a difference, outputs must be identical. "Identical"
means equal record digests after excluding headers, timing, thread counts,
paths and version strings.

| Stage | Expected | Reason |
|---|---|---|
| extract, illumina, annotate, probeset-align, probe-demux, cells, audit | identical | 0.5.7 moved this code into modules with no intended behaviour change. Golden tests do **not** cover these stages, so this comparison is their main real-data check. |
| molecules, default serial path (driver) | identical BAM records and matrix | 0.5.9 restored the 0.5.6 logic; golden-tested on fixtures |
| molecules `--evidence-workers` > 1 | representative choice may differ on equal-quality ties; matrix identical | 0.6.0 picks the first physical record on a tie |
| molecules, multiplexed Flex (PB) | possible differences where `PB` identity is involved | 0.6.0 uses the same present-`PB` identity for output as for evidence |
| guides, CRISPR molecules/matrix, protospacer-call | identical | integer CRISPR counts; the 0.6.0 protospacer change only affects fractional input |
| molecules `--tagseq` | identical (shared counting in both) | shared is the default again in 0.6.0 |
| MEX `matrix.mtx` line 2, H5 attributes, BAM `@PG` | always differ | they contain the software version; excluded |
| fractional, position-aware, sharded + saturation error | 0.6.0 only | new or changed in 0.6.0; tested on 0.6.0 alone |

## Datasets and representative subsets

| Assay | Project | SBX input | Illumina input | Independent reference |
|---|---|---|---|---|
| GEX (3') | CAT02 | `raw_fastqs_run_5p/device_{0..3}/out` (185 × ~4 GB chunks per device, 14 multiplexed samples) | `illumina_data/20260107_LH00416_0356_A235NMFLT4/CE15204` | Cell Ranger `illumina_data/cellranger_outputs/M2_GEX/outs` |
| Flex v2 multi (GEX + CRISPR), 8 samples + Pos/Neg | CAT06 | `raw_fastqs` (754 × ~2 GB chunks, ~43 M reads each, 4 libraries) | `20260806_LH00416_0423_B23WVVYLT4/RZ16443` | Cell Ranger multi `illumina_data/RZ16443_L1_multi_full/outs` |
| CAT tag-seq, 48 + 48 samples | DZEC04 | none | `DZEC04_set{1,2}_*_L008_R{1,2}_001.fastq.gz` | none (both builds are compared with each other) |

Subset choices:

- **CAT02 GEX: sample `M2_GEX`**, the only sample with a full Cell Ranger run
  (32,357 cells, 2.8 B reads). SBX: two mid-run chunks (`000090`, `000091`)
  from each of the four devices, 8 chunks of ~27.6 M raw reads each
  (~221 M reads, all 14 pooled samples). Illumina: lane 3, depth-matched
  (below).
- **CAT06 Flex multi: library L1 = `TSeq100_L1_GEX` + `TSeq100_L1_CRI`**
  (gene expression and CRISPR Guide Capture from the same cells), multiplexed
  by probe barcode (`[samples]` in `RZ16443_L1_multi_config.csv`). SBX: five consecutive
  mid-run chunks (`000300`–`000304`), as in the earlier 5-shard orchestrator
  validation; ~216 M raw reads across 4 pooled libraries. Illumina: lane 4,
  depth-matched (below). The driver supports
  Illumina Flex as **single-sample only** (no `PB` tag), so compare it
  pooled at gene level. The driver does not count CRISPR, so the L1
  CRISPR library runs as separate steps after it (step 2). L2 is not
  tested; it is an independent replicate library.
- **DZEC04 tag-seq: `set2_post_RT`** (767 M reads, scored with
  `tagseq_set1_barcodes.txt` as in `analysis/README.txt`). Two parts:
  1. **Full depth:** molecules → matrix → cells on the existing
     project-annotated BAM. UMI collisions depend on depth, so the
     collision handling can only be judged at full depth.
  2. **End to end:** the commands of `run_DZEC04_tagseq_pipeline.sh`
     (illumina `--overlap --trim-tso`, non-splice genome index, GENCODE v45,
     molecules, matrix, cells) on a 20 M-pair subset, run once with 0.5.6
     and once with 0.6.0. The pipeline script is not used for DZEC04.

### Depth matching: SBX raw data is pooled, Illumina FASTQs are per sample

SBX device output has **not** been demultiplexed: every chunk holds reads
from every library in the pool, and `sbxpress extract` assigns them to
samples. Illumina FASTQs come from BCL Convert already split by sample. The
same number of input reads therefore gives very different per-sample depth:

| Library | Share of SBX reads | SBX subset | Expected SBX reads for the sample |
|---|---|---|---|
| CAT02 `M2_GEX` | 6.9% of raw reads (measured: 27,461 of 400,000 reads of chunk `000090`, 0.6.0 extract) | 8 chunks, ~221 M raw | ~15 M |
| CAT06 `TSeq100_L1_GEX` | 27.4% of raw reads (8.94 B of 32.6 B, `sbxpress_extract_v2_pinned/metrics_summary.json`) | 5 chunks, ~216 M raw | ~59 M |
| CAT06 `TSeq100_L1_CRI` | 6.8% of raw reads (2.21 B of 32.6 B, same file) | same 5 chunks | ~15 M |

So the Illumina subset is never a fixed size. It is cut to **the number of
read pairs SBX extraction assigned to that sample** (`demux.per_sample`
in the SBX run's `extract/metrics_summary.json`). This matches the two
platforms at the same point: reads assigned to the sample, before
cell-barcode validation, orientation and mapping. Everything after that
point (valid barcodes, usable cDNA, mapped, counted molecules) is a
platform property, to be reported side by side, not normalised away.

What depth matching does and does not cover:

- **0.5.6 vs 0.6.0:** unaffected. Each comparison uses identical inputs.
- **SBX vs Illumina (same build):** like-for-like at matched depth. Compare
  funnels (assigned → valid CB → mapped → molecules), cells called, median
  UMIs/genes per cell, and per-gene concordance.
- **Either platform vs Cell Ranger (full depth, 2.8 B reads for `M2_GEX`):**
  not depth-matched, because the subset holds only ~1% of the reads.
  Compare only CPM-normalised gene-level agreement, on Cell Ranger's
  filtered barcodes. Do not compare totals, cell counts, saturation or
  per-cell UMIs with Cell Ranger.

If SBX assigns too few reads to the sample (under ~10 M), add chunks
(e.g. `000092` from each device) for both builds together, in new output
directories.

## 0. Environment, builds and the fixture-level baseline check

Create one shared environment file and source it in **every** new shell
(tmux, nohup, a second login). It sets the paths and defines `use_build`, which refuses to continue unless the requested
build's `sbxpress` is first on `PATH`.

```bash
VAL=/mnt/analysis/active_projects/sbxpress_validation/v0.5.6_vs_v0.6.0
mkdir -p "$VAL"/{src,install,inputs,runs,compare,logs}
cat > "$VAL/env.sh" <<'EOF'
VAL=/mnt/analysis/active_projects/sbxpress_validation/v0.5.6_vs_v0.6.0
SRC=$VAL/src
REPO=/mnt/analysis/software_shared/sbxpress
ABRF=/mnt/analysis/active_projects/roche_SBX_ABRF
CAT02=$ABRF/CAT02
CAT06=$ABRF/CAT06
DZEC04=/mnt/analysis/active_projects/cat_tagseq/DZEC04
TAGSEQ_WL=/mnt/analysis/active_projects/cat_tagseq/DMTS60/DMTS60_balanced_max_distant_sets1-4_whitelist.txt
export BASE_PATH=${BASE_PATH:-$PATH}
use_build() {
    export PATH="$VAL/install/$1/bin:$BASE_PATH"
    hash -r
    if [[ $(command -v sbxpress) != "$VAL/install/$1/bin/sbxpress" ]]; then
        echo "wrong sbxpress on PATH: $(command -v sbxpress)" >&2
        return 1
    fi
    sbxpress --version
}
python3 -c 'import numpy, h5py, pptx' || echo "python3 lacks numpy/h5py/python-pptx" >&2
EOF
source "$VAL/env.sh"
```

Check out both builds as worktrees of the existing repository (no network
clone), then install each into its own prefix:

```bash
git -C "$REPO" fetch --tags origin
git -C "$REPO" worktree add --detach "$SRC/v0.6.0" v0.6.0
git -C "$REPO" worktree add --detach "$SRC/v0.5.6" 4d42a18
cd "$SRC/v0.6.0"
scripts/check.sh > "$VAL/logs/check-v0.6.0.log" 2>&1
./install.sh --prefix "$VAL/install/v0.6.0" --skip-packages --skip-rust
cd "$SRC/v0.5.6"
./install.sh --prefix "$VAL/install/v0.5.6" --skip-packages --skip-rust
source "$VAL/env.sh"
use_build v0.5.6
use_build v0.6.0
sbxpress molecules --version
sbxpress matrix --version
```

`use_build` must print 0.5.6 and 0.6.0 respectively. Record the commits,
tool versions and machine:

```bash
{ git -C "$SRC/v0.5.6" rev-parse HEAD; git -C "$SRC/v0.6.0" rev-parse HEAD
  samtools --version-only; minimap2 --version; nproc; free -g; } > "$VAL/logs/provenance.txt"
```

Before using real data, run the golden-output test from the 0.6.0 tree
against the 0.5.6 binaries:

```bash
cd "$SRC/v0.6.0"
SBXPRESS_TEST_BIN_DIR="$SRC/v0.5.6/target/release" SBXPRESS_GOLDEN_LEGACY=1 \
  cargo test --release --test golden_outputs > "$VAL/logs/golden-v0.5.6.log" 2>&1
```

Expected result (verified on this server with a scratch 0.5.6 build): the
test reports **FAILED** with exactly one difference, `gex_indexed`, where
only the `bam` digest differs; its `counters` and `matrix` digests are
identical. That is the documented indexed tie-break change. Any other
listed case, or a differing counter or matrix, blocks acceptance.

The `python3` on `PATH` (miniconda) already has numpy, h5py and
python-pptx, which `check.sh` and the comparison scripts need; no venv
is required. `use_build` only adds the sbxpress install in front of `PATH`.

## 1. Build the input subsets (once, shared by both builds)

```bash
source "$VAL/env.sh"
mkdir -p "$VAL/inputs"/{cat02_sbx,cat06_sbx,illumina}
for d in 0 1 2 3; do
  for n in 000090 000091; do
    ln -s "$CAT02/raw_fastqs_run_5p/device_$d/out/PrimaryAnalysis_$n.fastq.gz" \
          "$VAL/inputs/cat02_sbx/device${d}_PrimaryAnalysis_$n.fastq.gz"
  done
done
for n in 000300 000301 000302 000303 000304; do
  ln -s "$CAT06/raw_fastqs/20260901_UCSF_Cat06_843_PT370_C3-$n.fastq.gz" "$VAL/inputs/cat06_sbx/"
done
subset() { pigz -dc "$1" | head -n 80000000 | pigz -p 8 > "$2"; }
subset "$DZEC04/DZEC04_set2_post_RT_S2_L008_R1_001.fastq.gz" "$VAL/inputs/illumina/DZEC04_set2_R1.fastq.gz"
subset "$DZEC04/DZEC04_set2_post_RT_S2_L008_R2_001.fastq.gz" "$VAL/inputs/illumina/DZEC04_set2_R2.fastq.gz"
( cd "$VAL/inputs" && find . \( -type f -o -type l \) | sort | while read -r f; do
    printf '%s\t%s\t%s\n' "$f" "$(readlink -f "$f")" "$(stat -L -c %s "$f")"; done ) > "$VAL/inputs/MANIFEST.tsv"
```

This step creates only the DZEC04 subset (20 M pairs; DZEC04 has no SBX
data to match). The CAT02 and CAT06 Illumina subsets are created in step 2,
after SBX extraction, so their depth can match it. R1 and R2 are cut at the
same record count, so the pairs stay in sync.
Never edit these inputs after the first run: the driver's restart manifest
records each input's size and modification time.

## 2. Run both builds

Run the builds **one after the other** on each dataset so that the
`/usr/bin/time` resource numbers are comparable. Thread counts are fixed and
identical for both builds. Each build/dataset pair gets its own new output
directory.

**Timing.** The server has ~1.1 TB of page cache, and inputs plus
references total ~150 GB. Because 0.5.6 runs first on every dataset, it
would read from NFS while 0.6.0 reads the same files from RAM. Read every
input and reference once before starting, so both builds start with a warm
cache. Repeat this after any long pause or other heavy I/O on the server:

```bash
source "$VAL/env.sh"
cat "$VAL"/inputs/cat02_sbx/* "$VAL"/inputs/cat06_sbx/* "$VAL"/inputs/illumina/*.fastq.gz \
    "$CAT02"/resources/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi \
    "$CAT02"/resources/gencode.v44.primary_assembly.basic.annotation.gtf.gz \
    /mnt/analysis/genomes/GRCh38/GRCh38.primary_assembly.genome.mmi \
    "$DZEC04"/sbxpress_illumina/DZEC04_set2_post_RT/annotation/DZEC04_set2_post_RT.annotated.coordinate.bam > /dev/null
```

Run the timed build pairs with nothing else heavy on the machine. Steps
that are not timed (the DZEC04 mode comparison, step 3 hashing) can run
in parallel with them. Compare per-stage times (`elapsed_secs` in each
metrics JSON, stage log timestamps) as well as the `/usr/bin/time` totals
and peak RSS. Expected wall time for the whole procedure is ~11–15 h run
sequentially, or ~9–11 h with the untimed steps run alongside; the CAT06
probe alignment (16 chunks aligned one after another) is the least certain
part.

Order matters: the SBX runs come first, because the Illumina subsets are
sized from the number of reads SBX extraction assigned to the sample (see
"Depth matching" above).

### CAT02 GEX and CAT06 Flex: SBX (raw, pooled device output)

```bash
source "$VAL/env.sh"
for b in v0.5.6 v0.6.0; do
  use_build "$b" || break
  mkdir -p "$VAL/runs/$b"
  /usr/bin/time -v -o "$VAL/runs/$b/cat02_sbx_gex.time" "$SRC/v0.6.0/scripts/sbxpress-pipeline.sh" \
    --input-mode sbx --analysis-mode gex --sample M2_GEX \
    --fastqs "$VAL/inputs/cat02_sbx" \
    --bclconvert-samplesheet "$CAT02/CE15204_bclconvert_input.csv" \
    --whitelist "$CAT02/resources/737K-arc-v1.txt.gz" \
    --mmi "$CAT02/resources/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi" \
    --gtf "$CAT02/resources/gencode.v44.primary_assembly.basic.annotation.gtf.gz" \
    --junctions "$CAT02/resources/gencode.v44.primary_assembly.junctions.bed" \
    --out-dir "$VAL/runs/$b/cat02_sbx_gex" --threads 16 --extract-threads 64 --map-threads 64 \
    --audit > "$VAL/runs/$b/cat02_sbx_gex.console.log" 2>&1
  /usr/bin/time -v -o "$VAL/runs/$b/cat06_sbx_flex.time" "$SRC/v0.6.0/scripts/sbxpress-pipeline.sh" \
    --input-mode sbx --analysis-mode flex --sample TSeq100_L1_GEX \
    --fastqs "$VAL/inputs/cat06_sbx" \
    --bclconvert-samplesheet "$CAT06/20260806_LH00416_0423_B23WVVYLT4/RZ16443/Reports/SampleSheet.csv" \
    --whitelist "$CAT06/resources/737K-flex-v2.txt.gz" \
    --probe-set "$CAT06/resources/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv" \
    --probe-barcode-whitelist "$CAT06/resources/flex-v2-384.txt" \
    --probe-sample-sheet "$CAT06/illumina_data/RZ16443_L1_multi_config.csv" \
    --out-dir "$VAL/runs/$b/cat06_sbx_flex" --threads 16 --extract-threads 64 --align-threads 32 \
    --audit > "$VAL/runs/$b/cat06_sbx_flex.console.log" 2>&1
done
```

### CAT06 Flex: CRISPR guide capture for L1 (SBX, multiplexed)

CAT06 is a Cell Ranger **multi** run: library L1 is `TSeq100_L1_GEX` plus
`TSeq100_L1_CRI` (CRISPR Guide Capture, `feature_reference_RZ.csv`), and each
sample needs both. The driver counts GEX only, but its SBX extraction above
already wrote the `TSeq100_L1_CRI` reads. These commands follow
`run_full_cat06.sh`: `sbxpress guides` on each extracted chunk, then
`run_flex_multi.sh --crispr-bam`, which reuses the driver's GEX molecule
BAM and probe-demux assignments. Per sample it writes GEX outputs plus CRISPR
matrices and protospacer calls (`flex_samples/<sample>/{gex_outs,crispr_outs}`).
The first chunk runs alone because it creates the guide feature manifest,
which the other chunks then verify.

```bash
source "$VAL/env.sh"
for b in v0.5.6 v0.6.0; do
  use_build "$b" || break
  R=$VAL/runs/$b/cat06_sbx_flex
  C=$R/crispr
  mkdir -p "$C/chunks"
  CRI_FASTQS=("$R"/extract/TSeq100_L1_CRI_S*_R2_*.fastq.gz)
  for i in "${!CRI_FASTQS[@]}"; do
    fq=${CRI_FASTQS[$i]}
    n=$(basename "$fq" .fastq.gz)
    sbxpress guides "$fq" --feature-reference "$CAT06/feature_reference_RZ.csv" \
      --features-out "$C/guide_features.tsv" --metrics "$C/chunks/$n.metrics.json" \
      --output "$C/chunks/$n.bam" 2> "$C/chunks/$n.log" &
    (( i == 0 )) && wait
  done
  wait
  samtools cat -o "$C/TSeq100_L1_CRI.guides.bam" "$C"/chunks/*.bam
  "$SRC/v0.6.0/scripts/run_flex_multi.sh" --crispr-bam "$C/TSeq100_L1_CRI.guides.bam" \
    --guide-features "$C/guide_features.tsv" \
    "$R/TSeq100_L1_GEX/molecules/TSeq100_L1_GEX.molecules.coordinate.bam" \
    "$R/TSeq100_L1_GEX/probe_demux" "$C/flex_samples" > "$C/run_flex_multi.log" 2>&1
done
```

Each sample's `gex_outs` here must be identical to the driver's
`TSeq100_L1_GEX/flex_samples/<sample>/gex_outs`, because both come from
the same molecule BAM and the same assignments.

### Depth-matched Illumina subsets

Take the sample's assigned SBX read count from the extract metrics
(`demux.per_sample[].pairs`). Extraction is expected to be identical in
both builds; if the two counts differ, stop: that difference is a finding
for step 3. Then cut each Illumina lane to the same number of read pairs:

```bash
source "$VAL/env.sh"
sbx_pairs() {  # sbx_pairs RUN SAMPLE -> reads SBX extraction assigned to SAMPLE
  local a b
  a=$(jq --arg s "$2" '.demux.per_sample[] | select(.sample==$s) | .pairs' "$VAL/runs/v0.5.6/$1/extract/metrics_summary.json")
  b=$(jq --arg s "$2" '.demux.per_sample[] | select(.sample==$s) | .pairs' "$VAL/runs/v0.6.0/$1/extract/metrics_summary.json")
  [[ -n $a && $a == "$b" ]] || { echo "extract counts differ or missing: $a vs $b" >&2; return 1; }
  echo "$a"
}
subset_n() { pigz -dc "$1" | head -n $(( 4 * $3 )) | pigz -p 8 > "$2"; }
N02=$(sbx_pairs cat02_sbx_gex M2_GEX)
N06=$(sbx_pairs cat06_sbx_flex TSeq100_L1_GEX)
N06C=$(sbx_pairs cat06_sbx_flex TSeq100_L1_CRI)
echo "CAT02 M2_GEX: $N02; CAT06 TSeq100_L1_GEX: $N06; TSeq100_L1_CRI: $N06C pairs" | tee "$VAL/inputs/illumina/DEPTH.txt"
I02=$CAT02/illumina_data/20260107_LH00416_0356_A235NMFLT4/CE15204
I06=$CAT06/20260806_LH00416_0423_B23WVVYLT4/RZ16443
subset_n "$I02/M2_GEX_S8_L003_R1_001.fastq.gz" "$VAL/inputs/illumina/CAT02_M2_GEX_R1.fastq.gz" "$N02"
subset_n "$I02/M2_GEX_S8_L003_R2_001.fastq.gz" "$VAL/inputs/illumina/CAT02_M2_GEX_R2.fastq.gz" "$N02"
subset_n "$I06/TSeq100_L1_GEX_S1_L004_R1_001.fastq.gz" "$VAL/inputs/illumina/CAT06_L1_GEX_R1.fastq.gz" "$N06"
subset_n "$I06/TSeq100_L1_GEX_S1_L004_R2_001.fastq.gz" "$VAL/inputs/illumina/CAT06_L1_GEX_R2.fastq.gz" "$N06"
subset_n "$I06/TSeq100_L1_CRI_S3_L004_R1_001.fastq.gz" "$VAL/inputs/illumina/CAT06_L1_CRI_R1.fastq.gz" "$N06C"
subset_n "$I06/TSeq100_L1_CRI_S3_L004_R2_001.fastq.gz" "$VAL/inputs/illumina/CAT06_L1_CRI_R2.fastq.gz" "$N06C"
```

One Illumina lane holds far more than either count (lane 3 of `M2_GEX`
alone is ~26 GB of FASTQ), so a single lane is enough. Add these files to
`MANIFEST.tsv` (rerun the manifest command from step 1).

### CAT02 GEX and CAT06 Flex: Illumina (depth-matched)

```bash
source "$VAL/env.sh"
for b in v0.5.6 v0.6.0; do
  use_build "$b" || break
  /usr/bin/time -v -o "$VAL/runs/$b/cat02_illumina_gex.time" "$SRC/v0.6.0/scripts/sbxpress-pipeline.sh" \
    --input-mode illumina --analysis-mode gex --sample M2_GEX \
    --r1 "$VAL/inputs/illumina/CAT02_M2_GEX_R1.fastq.gz" --r2 "$VAL/inputs/illumina/CAT02_M2_GEX_R2.fastq.gz" \
    --whitelist "$CAT02/resources/737K-arc-v1.txt.gz" \
    --mmi "$CAT02/resources/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi" \
    --gtf "$CAT02/resources/gencode.v44.primary_assembly.basic.annotation.gtf.gz" \
    --junctions "$CAT02/resources/gencode.v44.primary_assembly.junctions.bed" \
    --out-dir "$VAL/runs/$b/cat02_illumina_gex" --threads 16 --extract-threads 64 --map-threads 64 \
    --audit > "$VAL/runs/$b/cat02_illumina_gex.console.log" 2>&1
  /usr/bin/time -v -o "$VAL/runs/$b/cat06_illumina_flex.time" "$SRC/v0.6.0/scripts/sbxpress-pipeline.sh" \
    --input-mode illumina --analysis-mode flex --sample TSeq100_L1_GEX \
    --r1 "$VAL/inputs/illumina/CAT06_L1_GEX_R1.fastq.gz" --r2 "$VAL/inputs/illumina/CAT06_L1_GEX_R2.fastq.gz" \
    --whitelist "$CAT06/resources/737K-flex-v2.txt.gz" \
    --probe-set "$CAT06/resources/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv" \
    --out-dir "$VAL/runs/$b/cat06_illumina_flex" --threads 16 --extract-threads 64 --align-threads 32 \
    --audit > "$VAL/runs/$b/cat06_illumina_flex.console.log" 2>&1
done
```

### CAT06 CRISPR, Illumina (pooled, exploratory)

Illumina Flex input has no probe barcode (`PB`), so CRISPR counts cannot be
split by sample and `molecules --flex` cannot be used. Count the
depth-matched `TSeq100_L1_CRI` subset pooled. `sbxpress guides` was built
and tuned on the SBX R2 layout. If its hit rate (`guides.metrics.json`) is
near zero here, record Illumina CRISPR as unsupported; that is not a 0.6.0
regression, as long as 0.5.6 gives the same result. (2026-10-08: 94.6% hit
rate on the CAT06 Illumina subset, so the R2 layout works.)

```bash
source "$VAL/env.sh"
for b in v0.5.6 v0.6.0; do
  use_build "$b" || break
  C=$VAL/runs/$b/cat06_illumina_crispr
  mkdir -p "$C"
  sbxpress illumina --r1 "$VAL/inputs/illumina/CAT06_L1_CRI_R1.fastq.gz" --r2 "$VAL/inputs/illumina/CAT06_L1_CRI_R2.fastq.gz" \
    --whitelist "$CAT06/resources/737K-flex-v2.txt.gz" --output "$C/tagged.R2.fastq.gz" \
    --metrics "$C/illumina.metrics.json" --threads 16 2> "$C/illumina.log"
  sbxpress guides "$C/tagged.R2.fastq.gz" --feature-reference "$CAT06/feature_reference_RZ.csv" \
    --features-out "$C/guide_features.tsv" --metrics "$C/guides.metrics.json" \
    --output "$C/guides.bam" 2> "$C/guides.log"
  sbxpress molecules "$C/guides.bam" --output "$C/molecules.bam" --metrics "$C/molecules.metrics.json" 2> "$C/molecules.log"
  sbxpress matrix "$C/molecules.bam" --features "$C/guide_features.tsv" \
    --cell-barcodes "$CAT06/resources/737K-flex-v2.txt.gz" --output "$C/outs" 2> "$C/matrix.log"
done
```

### DZEC04 tag-seq, part 1: full depth from the project's annotated BAM

The same flags as `run_DZEC04_tagseq_pipeline.sh`. `--evidence-workers 16`
also tests the indexed tie-break (table above).

```bash
source "$VAL/env.sh"
SET2=$DZEC04/sbxpress_illumina/DZEC04_set2_post_RT
for b in v0.5.6 v0.6.0; do
  use_build "$b" || break
  O=$VAL/runs/$b/dzec04_full
  mkdir -p "$O"
  /usr/bin/time -v -o "$O/molecules.time" sbxpress molecules "$SET2/annotation/DZEC04_set2_post_RT.annotated.coordinate.bam" \
    --tagseq --evidence-workers 16 --threads 32 --saturation "$O/saturation.tsv.gz" \
    --output "$O/molecules.bam" --metrics "$O/molecules.metrics.json" 2> "$O/molecules.log"
  sbxpress matrix "$O/molecules.bam" --output "$O/outs" --genome GRCh38 \
    --features "$SET2/reference/genes.mtx" --cell-barcodes "$TAGSEQ_WL" 2> "$O/matrix.log"
  sbxpress cells "$O/outs/raw_feature_bc_matrix.h5" --output "$O/cells" \
    --tagseq-barcodes "$DZEC04/tagseq_set1_barcodes.txt" 2> "$O/cells.log"
done
```

Then use 0.6.0 alone to compare its counting modes at full depth, with the
0.5.6 shared matrix as the reference. No independent count exists for
DZEC04, so this result describes the modes but cannot rank them for accuracy.
Record each mode's expression-dependent bias slope and its total relative
to shared.

```bash
use_build v0.6.0
SBXPRESS=$(command -v sbxpress) bash "$SRC/v0.6.0/scripts/compare_tagseq_modes.sh" \
  --bam "$SET2/annotation/DZEC04_set2_post_RT.annotated.coordinate.bam" \
  --reference "$VAL/runs/v0.5.6/dzec04_full/outs/raw_feature_bc_matrix" \
  --features "$SET2/reference/genes.mtx" --barcodes "$DZEC04/tagseq_set1_barcodes.txt" \
  --out-dir "$VAL/runs/v0.6.0/dzec04_modes" --threads 32 \
  --modes "shared fractional position position_frac" --molecule-args "--evidence-workers 16"
```

`shared` must match the 0.5.6 reference exactly (Pearson r = 1, slope 0,
total ratio 1). The helper reuses any existing `matrix.mtx`, so delete a
mode's subdirectory after an interrupted run.

### DZEC04 tag-seq, part 2: end to end on the subset

```bash
source "$VAL/env.sh"
G=/mnt/analysis/genomes/GRCh38
for b in v0.5.6 v0.6.0; do
  use_build "$b" || break
  O=$VAL/runs/$b/dzec04_subset
  mkdir -p "$O/extract" "$O/reference"
  sbxpress illumina --r1 "$VAL/inputs/illumina/DZEC04_set2_R1.fastq.gz" --r2 "$VAL/inputs/illumina/DZEC04_set2_R2.fastq.gz" \
    --whitelist "$TAGSEQ_WL" --output "$O/extract/set2.tagged.R2.fastq.gz" --failed-output "$O/extract/set2.failed" \
    --metrics "$O/extract/illumina.metrics.json" --threads 8 --split 4 --overlap --trim-tso 2> "$O/illumina.log"
  sbxpress minimap -a -x splice:sr --secondary=no -y -t 64 "$G/GRCh38.primary_assembly.genome.mmi" \
      "$O"/extract/set2.tagged.R2.shard-*.fastq.gz 2> "$O/minimap2.log" \
    | sbxpress annotate - --gtf "$G/gencode.v45.annotation.gtf" --output - --metrics "$O/annotation.metrics.json" \
        --genes-mtx "$O/reference/genes.mtx" --three-prime --unordered --threads 32 2> "$O/annotate.log" \
    | samtools sort -@ 16 -m 2G -T "$O/sort_tmp" -o "$O/annotated.bam" -
  samtools index "$O/annotated.bam"
  sbxpress molecules "$O/annotated.bam" --tagseq --threads 16 \
    --output "$O/molecules.bam" --metrics "$O/molecules.metrics.json" 2> "$O/molecules.log"
  sbxpress matrix "$O/molecules.bam" --output "$O/outs" --genome GRCh38 \
    --features "$O/reference/genes.mtx" --cell-barcodes "$TAGSEQ_WL" 2> "$O/matrix.log"
  sbxpress cells "$O/outs/raw_feature_bc_matrix.h5" --output "$O/cells" \
    --tagseq-barcodes "$DZEC04/tagseq_set1_barcodes.txt" 2> "$O/cells.log"
done
```

After each run, check that every log is free of errors and that
`samtools quickcheck -v` passes on every BAM.

## 3. Compare 0.5.6 with 0.6.0, stage by stage

Create the comparison helpers once. They leave out timing, threads, paths
and version fields. BAMs are first compared in file order, which is
deterministic for molecules output. Only if that differs is the comparison
repeated order-independently (sorted SAM text), because annotation runs
`--unordered` and coordinate ties can come out in either order. Sorting
spills to `$VAL/tmp` on the analysis volume: sorted SAM text for a
full-depth BAM is hundreds of GB, more than the 329 GB free on `/`.

```bash
cat > "$VAL/compare.sh" <<'EOF'
# Drops timing/thread/path keys and the three keys new in 0.6.0 molecules metrics;
# normalises the build name inside path values (genes_mtx, raw_matrix, bam, ...).
jnorm() { jq -S 'walk(if type=="object" then with_entries(select(.key|test("elapsed|secs|per_sec|thread|worker|path|file|input|output|dir|version|command|^fractional_collisions$|^position_aware_umis$|^collision_molecules_fractional$";"i")|not)) elif type=="string" then gsub("/runs/v0\\.5\\.6/|/runs/v0\\.6\\.0/";"/runs/BUILD/") else . end)' "$1"; }
mkdir -p "$VAL/tmp"
fq_digest() { pigz -dc "$@" | paste - - - - | LC_ALL=C sort -T "$VAL/tmp" -S 32G --parallel 16 | md5sum | cut -c1-32; }
bam_digest() { samtools view -@ 8 "$1" | md5sum | cut -c1-32; }
bam_digest_sorted() { samtools view -@ 8 "$1" | LC_ALL=C sort -T "$VAL/tmp" -S 32G --parallel 16 | md5sum | cut -c1-32; }
bam_same() {
  if [[ $(bam_digest "$2") == "$(bam_digest "$3")" ]]; then echo "SAME $1"
  elif [[ $(bam_digest_sorted "$2") == "$(bam_digest_sorted "$3")" ]]; then echo "SAME $1 (record order differs)"
  else echo "DIFF $1"; fi
}
mex_digest() {
  local d=$1 f
  for f in features.tsv barcodes.tsv matrix.mtx; do
    if [[ -e $d/$f.gz ]]; then pigz -dc "$d/$f.gz"; else cat "$d/$f"; fi | grep -v '^%' | md5sum | cut -c1-32
  done | paste -sd,
}
same() { if [[ $2 == "$3" ]]; then echo "SAME $1"; else echo "DIFF $1"; fi; }
# compare_run NAME: every JSON metric, BAM and MEX matrix under runs/v0.5.6/NAME vs runs/v0.6.0/NAME
compare_run() {
  local a=$VAL/runs/v0.5.6/$1 b=$VAL/runs/v0.6.0/$1 f
  echo "== $1"
  while read -r f; do same "$f" "$(jnorm "$a/$f" | md5sum)" "$(jnorm "$b/$f" | md5sum)"
  done < <(cd "$a" && find . -name '*.json' ! -name '.pipeline-state.json' ! -path './tmp/*' | sort)
  while read -r f; do bam_same "$f" "$a/$f" "$b/$f"
  done < <(cd "$a" && find . -name '*.bam' ! -path '*/probeset_chunks/*' ! -path '*/chunks/*' ! -path './tmp/*' | sort)
  while read -r f; do same "$f" "$(mex_digest "$a/$f")" "$(mex_digest "$b/$f")"
  done < <(cd "$a" && find . -name 'matrix.mtx*' -printf '%h\n' | sort -u)
  while read -r f; do same "$f" "$(sort "$a/$f" | md5sum)" "$(sort "$b/$f" | md5sum)"
  done < <(cd "$a" && find . \( -name '*.barcodes.txt' -o -name 'filtered_barcodes.csv' -o -name 'protospacer_calls_*.csv' \) | sort)
}
EOF
source "$VAL/env.sh"
source "$VAL/compare.sh"
for run in cat02_sbx_gex cat02_illumina_gex cat06_sbx_flex cat06_illumina_flex cat06_illumina_crispr dzec04_full dzec04_subset; do
  compare_run "$run"
done | tee "$VAL/compare/stage_digests.txt"
```

Also compare the extracted reads for the sample under test (SBX:
`extract/<SAMPLE>_S*_R2_*.fastq.gz`; Illumina: `extract/<SAMPLE>.tagged.R2.fastq.gz`;
DZEC04: `extract/set2.tagged.R2.shard-*.fastq.gz`), for example:

```bash
same cat02_sbx_extract "$(fq_digest "$VAL"/runs/v0.5.6/cat02_sbx_gex/extract/M2_GEX_S*_R2_*.fastq.gz)" \
                       "$(fq_digest "$VAL"/runs/v0.6.0/cat02_sbx_gex/extract/M2_GEX_S*_R2_*.fastq.gz)"
```

For every `DIFF`:

1. Locate it: the earliest differing stage is the source; later stages
   inherit it. To confirm, run the 0.6.0 stage on the 0.5.6 output of the
   previous stage (e.g. 0.6.0 `sbxpress molecules` on the 0.5.6 annotated BAM)
   in a new directory, and check whether the difference persists.
2. Size it: `jq` diff of the two metrics files, and
   `diff <(samtools view A | sort) <(samtools view B | sort) | head` for BAMs.
3. Classify it: an expected change from the table, or a bug.

### Concordance (how large is any difference, and against what)

`compare_tagseq_modes.py` compares any MEX matrices with a reference (MEX or
10x `.h5`). It reports per-barcode correlation, per-gene log2 ratio,
expression-dependent bias slope and total ratio.

```bash
P=$SRC/v0.6.0/scripts/compare_tagseq_modes.py
R=$VAL/runs
python3 "$P" --reference "$R/v0.5.6/cat02_sbx_gex/M2_GEX/outs/raw_feature_bc_matrix" \
  --mode "v0.6.0=$R/v0.6.0/cat02_sbx_gex/M2_GEX/outs/raw_feature_bc_matrix" \
  --out-dir "$VAL/compare/cat02_sbx_v056_vs_v060"
python3 "$P" --reference "$CAT02/illumina_data/cellranger_outputs/M2_GEX/outs/filtered_feature_bc_matrix.h5" \
  --mode "sbx_v0.5.6=$R/v0.5.6/cat02_sbx_gex/M2_GEX/outs/raw_feature_bc_matrix" \
  --mode "sbx_v0.6.0=$R/v0.6.0/cat02_sbx_gex/M2_GEX/outs/raw_feature_bc_matrix" \
  --mode "illumina_v0.5.6=$R/v0.5.6/cat02_illumina_gex/M2_GEX/outs/raw_feature_bc_matrix" \
  --mode "illumina_v0.6.0=$R/v0.6.0/cat02_illumina_gex/M2_GEX/outs/raw_feature_bc_matrix" \
  --out-dir "$VAL/compare/cat02_vs_cellranger"
for b in v0.5.6 v0.6.0; do
  python3 "$P" --reference "$R/$b/cat02_illumina_gex/M2_GEX/outs/filtered_feature_bc_matrix" \
    --mode "sbx_$b=$R/$b/cat02_sbx_gex/M2_GEX/outs/raw_feature_bc_matrix" \
    --out-dir "$VAL/compare/cat02_sbx_vs_illumina_$b"
done
```

The last comparison is the only depth-matched cross-platform one: SBX and
Illumina were given the same number of `M2_GEX` reads. Alongside it,
tabulate for each platform and build the read funnel from the metrics JSON
(assigned → valid cell barcode → mapped/annotated → counted molecules), the
number of cells called, and median UMIs and genes per cell. The SBX/Illumina
ratios must be the same for 0.5.6 and 0.6.0; a ratio that moves between
builds points to a stage that changed on one platform only.

For CAT06, compare each sample's
`cat06_sbx_flex/TSeq100_L1_GEX/flex_samples/<sample>/gex_outs` against
`RZ16443_L1_multi_full/outs/per_sample_outs/<sample>/sample_filtered_feature_bc_matrix.h5`
for both builds. Cell Ranger's Flex per-sample barcodes are the 16-bp GEM
barcode plus the 10-bp probe barcode (`AAACCAATCAAGCGTGACCCTATATC-1`);
sbxpress keys a sample's cells by the GEM barcode alone. Truncate Cell
Ranger's barcodes to 16 bp and sum a GEM's probe barcodes within the sample
before comparing, or the comparison finds no shared barcodes. Compare the
Illumina Flex run with the pooled
`RZ16443_L1_multi_full/outs/raw_feature_bc_matrix.h5` at gene level.
Check the "genes compared" column: the Cell Ranger and GENCODE v44 feature
IDs must overlap. The subsets are ~1% of the reads Cell Ranger used and
are not depth-matched to it, so judge agreement with Cell Ranger only by
CPM-normalised correlation and bias slope on its filtered barcodes, never
by totals, cell counts or saturation. For CRISPR, compare each sample's
`crispr/flex_samples/<sample>/crispr_outs/protospacer_calls_summary.csv`
with Cell Ranger's `per_sample_outs/<sample>/crispr_analysis/protospacer_calls_summary.csv`:
the same guides should dominate, and the fraction of cells with a call
should be similar in both builds. At ~1% depth, call rates will be lower
than Cell Ranger's. Run the CAT06 SBX-vs-Illumina
comparison at gene level (pooled), since Illumina Flex has no per-sample split. The two
builds should score identically.

## 4. Bug checks on 0.6.0 (and the same checks on 0.5.6 for reference)

### Indexed vs sharded molecules parity

Use the CAT02 SBX annotated BAM and feature manifest from step 2, and run
this for both builds. Matrix output does not depend on input order, so the
sharded BAM is sorted only for `quickcheck`.

```bash
source "$VAL/env.sh"
source "$VAL/compare.sh"
for b in v0.5.6 v0.6.0; do
  use_build "$b" || break
  IN=$VAL/runs/$b/cat02_sbx_gex
  BAM=$IN/M2_GEX/annotation/M2_GEX.annotated.coordinate.bam
  FEATURES=$IN/reference/genes.mtx
  PARITY=$VAL/runs/$b/parity
  mkdir -p "$PARITY"
  sbxpress molecules "$BAM" --evidence-workers 4 --threads 8 \
    --output "$PARITY/indexed.bam" --metrics "$PARITY/indexed.metrics.json"
  sbxpress molecules "$BAM" --shards 8 --shard-concurrency 1 --threads 8 \
    --output "$PARITY/sharded.bam" --metrics "$PARITY/sharded.metrics.json"
  samtools sort -@ 8 -o "$PARITY/sharded.sorted.bam" "$PARITY/sharded.bam"
  samtools quickcheck -v "$PARITY/indexed.bam" "$PARITY/sharded.sorted.bam"
  sbxpress matrix "$PARITY/indexed.bam" --features "$FEATURES" --output "$PARITY/indexed" --no-gzip
  sbxpress matrix "$PARITY/sharded.sorted.bam" --features "$FEATURES" --output "$PARITY/sharded" --no-gzip
  same "parity $b" "$(mex_digest "$PARITY/indexed/raw_feature_bc_matrix")" \
                   "$(mex_digest "$PARITY/sharded/raw_feature_bc_matrix")"
done
```

Integer matrices must be identical. Do not compare fractional matrices
byte for byte: summation order differs between strategies, which changes
the last digit. To test fractional parity (0.6.0 only, with
`--tagseq --fractional-collisions` added to both molecules commands),
require identical entries and a maximum absolute difference below 1e-9:

```bash
paste <(grep -v '^%' "$PARITY/indexed/raw_feature_bc_matrix/matrix.mtx") \
      <(grep -v '^%' "$PARITY/sharded/raw_feature_bc_matrix/matrix.mtx") |
  awk 'NR>1 { if ($1!=$4 || $2!=$5) bad++; d=$3-$6; if (d<0) d=-d; if (d>m) m=d }
       END { print "entry mismatches:", bad+0, "max |diff|:", m+0 }'
```

### Restart and interruption (0.6.0)

Use the CAT02 Illumina inputs, the fastest end-to-end case. Run the
interruption tests in a fresh directory, one per stage. The kill is aimed
at a stage that is still running: [3/5] is molecules, [4/5] is the raw
matrix, [5/5] is cell calling.

```bash
source "$VAL/env.sh"
source "$VAL/compare.sh"
use_build v0.6.0
set -m   # job control: each background job gets its own process group to kill
for stage in 3 4 5; do
  R=$VAL/runs/v0.6.0/restart_stage$stage
  "$SRC/v0.6.0/scripts/sbxpress-pipeline.sh" --input-mode illumina --analysis-mode gex --sample M2_GEX \
    --r1 "$VAL/inputs/illumina/CAT02_M2_GEX_R1.fastq.gz" --r2 "$VAL/inputs/illumina/CAT02_M2_GEX_R2.fastq.gz" \
    --whitelist "$CAT02/resources/737K-arc-v1.txt.gz" \
    --mmi "$CAT02/resources/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi" \
    --gtf "$CAT02/resources/gencode.v44.primary_assembly.basic.annotation.gtf.gz" \
    --junctions "$CAT02/resources/gencode.v44.primary_assembly.junctions.bed" \
    --out-dir "$R" --threads 16 --extract-threads 64 --map-threads 64 > "$R.console.1" 2>&1 &
  PGID=$!
  until grep -q "^\[$stage/5\]" "$R.console.1" 2>/dev/null; do sleep 2; done
  sleep 3   # molecules on this subset can finish within ~20 s
  kill -TERM -- -"$PGID"
  wait "$PGID"
  ls -la "$R/M2_GEX/outs" > "$R.after_kill.txt"
done
```

After each kill: check that `outs/` holds no `.raw.done` / `.cells.done`
marker for the killed stage, and no H5 or MEX file left from it. Then
rerun the same command into the same `R` (console to `$R.console.2`).
Earlier stages must be reported as reused, and the final outputs must be
identical to the uninterrupted `runs/v0.6.0/cat02_illumina_gex` (all
digests except the `restart_stage*` paths).

Then check the configuration guard: rerun one restart directory with an
extra `--no-three-prime`. It must stop with "inputs or analysis settings
changed … use a new output directory". If the kill lands after the stage
has already finished, repeat with a shorter `sleep`.

## Acceptance

- `check.sh` passes; golden legacy test passes (indexed tie-break excepted);
  both installs report the expected version.
- Every run in step 2 completes for both builds, with readable BAM/H5/MEX outputs and clean logs.
- Every `DIFF` in `stage_digests.txt` is explained by the expected-differences
  table. **Any other difference blocks acceptance.** This is especially
  important for extract/illumina/annotate, which no fixture test covers.
- `dzec04_modes` `shared` equals 0.5.6; every mode's bias slope and total are recorded.
- 0.6.0 agrees with Cell Ranger at least as well as 0.5.6 (CAT02 and every CAT06 sample).
- Depth-matched SBX/Illumina ratios (funnel, cells, UMIs/genes per cell)
  are the same for both builds.
- Integer indexed/sharded matrices are identical; fractional matrices are
  within tolerance. Interrupted runs never publish partial matrix or cell
  outputs, and a changed configuration is refused.
- CAT06 L1 produces GEX **and** CRISPR outputs (matrices, protospacer calls)
  for every sample in both builds; CRISPR outputs are identical between builds.
- No sample mixing across CAT06 Flex samples: each sample's barcodes
  and totals match Cell Ranger's sample assignment in proportion.
- Peak RSS and elapsed time per stage are within ~20% of 0.5.6, unless the
  changelog explains the difference (e.g. libdeflate, single-pass tags).

## Evidence to archive

Archive `$VAL/logs` (check, golden, provenance), `$VAL/inputs/MANIFEST.tsv`,
every `*.console.log`, `*.time` and stage log, all metrics JSON,
`compare/stage_digests.txt`, every `compare/*/summary.tsv`, and a short
written explanation for each `DIFF`. Keep the run directories until the
release is accepted.

## After acceptance

Scale to one full library per assay with 0.6.0 (CAT02 `M2_GEX` all devices
and lanes; CAT06 `TSeq100_L1_GEX` all chunks; DZEC04 both sets), and repeat
the resource and Cell Ranger concordance checks. Merge `fix/gemini-review`
into `main`, then replace `/usr/local/bin` (currently 0.3.5). Remove the
worktrees with `git -C "$REPO" worktree remove "$SRC/v0.5.6"` (and `v0.6.0`)
when done.
