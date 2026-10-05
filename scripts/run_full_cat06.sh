#!/usr/bin/env bash
# Full-dataset, full-depth CAT06 Flex + CRISPR run (FLEX_SUPPORT_PLAN.md
# Phase 6's remaining literal bar: all real raw shards, not the 800k-read
# validation sample used through Phase 6's own results section).
#
# CPU budget: 128 cores, allocated in two different ways at two different
# stages, per the shape of the work:
#   - extract:            one process, --threads 128 (rayon pool sized to
#                          128), --split 128 so gzip compression also gets
#                          128-way parallelism (default --split 1 means only
#                          (n_samples+1) mutex-guarded encoders regardless of
#                          --threads -- see PROFILING.md's Phase-1 run, which
#                          hit only 107,862 reads/sec at effectively ~5-way
#                          write parallelism; --split 128 removes that floor).
#   - probeset-align/guides: neither binary is internally parallel (single
#                          reader, no rayon) -- so the 128-way split above is
#                          reused here as 128-way *process* parallelism via
#                          xargs -P128, one chunk per core.
#   - molecules/probe-demux/matrix/cells: single-threaded each; run per
#                          library-pair (L1, L2), negligible against the
#                          128-core budget.
#
# Real inputs (already on disk, read-only, confirmed present before writing
# this script):
#   raw reads:        CAT06/raw_fastqs/*.fastq.gz            (754 shards, ~1.54TB, ~33B reads est.)
#   sample sheet:      CAT06/20260806_.../RZ16443/Reports/SampleSheet.csv (real BCLConvert v2 sheet,
#                       4 libraries: TSeq100_{L1,L2}_{GEX,CRI}, confirmed via distinct i7/i5 pairs)
#   CB whitelist:      cellranger-9.0.1's real 737K-fixed-rna-profiling.txt.gz (the actual Flex
#                       GEM barcode whitelist -- NOT the ad hoc "observed barcodes" whitelist used
#                       for validation earlier in this plan; that one is ground-truth-derived and
#                       circular for a production run)
#   probe barcode wl:  CAT06/resources/flex-v2-384.txt
#   probe-set ref:     CAT06/resources/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv
#   feature ref:       CAT06/feature_reference_RZ.csv
#   probe sample sheet ([samples] section, sample_id -> probe_barcode_ids):
#                       CAT06/illumina_data/RZ16443_L1_multi/outs/config.csv (same file used
#                       throughout this plan)
#
# Known, stated assumption (not verified against any ground truth -- there is
# none for L2): cellranger's own real config.csv only analyzes L1
# (TSeq100_L1_GEX/TSeq100_L1_CRI) -- L2 is a second, independently-indexed
# library that was sequenced but never fed to the real cellranger run this
# plan validates against. L1 and L2 are therefore kept fully separate through
# every stage (never pooled), each producing its own results/L{1,2}/ tree.
# Only L1's results have real ground truth to compare against
# (illumina_data/RZ16443_L1_multi/outs/...); L2's are informational.
#
# Resumable: every stage checks for its expected output before (re)running,
# so a killed/restarted run picks back up rather than redoing finished work.
#
# Usage: scripts/run_full_cat06.sh [output_dir]
# Intended invocation (survives disconnect):
#   setsid nohup scripts/run_full_cat06.sh > full_cat06_driver.log 2>&1 < /dev/null &
#   disown

set -euo pipefail

CAT06=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT06
BIN=/mnt/analysis/software_shared/sbxpress/target/release
OUT=${1:-$CAT06/sbxpress_full_run}
CORES=128

SAMPLE_SHEET="$CAT06/20260806_LH00416_0423_B23WVVYLT4/RZ16443/Reports/SampleSheet.csv"
# 737K-flex-v2.txt.gz, NOT 737K-fixed-rna-profiling.txt.gz (older Fixed RNA
# Profiling kit's whitelist -- confirmed wrong for this dataset: 0/735,036
# real observed CAT06 GEM barcodes matched it beyond coincidence, vs.
# 735,036/735,036 exact matches against 737K-flex-v2.txt.gz). Only shipped
# under cellranger-10.0.0 on this box; 9.0.1/9.0.0/8.0.1 don't have it at all.
WHITELIST=/mnt/analysis/software/cellranger-10.0.0/lib/python/cellranger/barcodes/737K-flex-v2.txt.gz
PB_WHITELIST="$CAT06/resources/flex-v2-384.txt"
PROBE_SET="$CAT06/resources/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv"
FEATURE_REF="$CAT06/feature_reference_RZ.csv"
PROBE_SAMPLE_SHEET="$CAT06/illumina_data/RZ16443_L1_multi/outs/config.csv"

export RAYON_NUM_THREADS=$CORES

mkdir -p "$OUT"/{extracted,aligned,molecules,demux,results,logs}
LOG_DIR="$OUT/logs"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

stage_extract() {
  local marker="$OUT/extracted/.done"
  if [ -f "$marker" ]; then
    log "extract: already done, skipping"
    return
  fi
  log "extract: starting (754 shards, --threads $CORES --split $CORES, gzip level 4 default)"
  "$BIN/sbx10x" \
    --fastqs "$CAT06/raw_fastqs" \
    -o "$OUT/extracted" \
    --bclconvert-samplesheet "$SAMPLE_SHEET" \
    --whitelist "$WHITELIST" \
    --probe-barcode-whitelist "$PB_WHITELIST" \
    --flex \
    --minimap2 \
    --threads "$CORES" \
    --split "$CORES" \
    > "$LOG_DIR/extract.log" 2>&1
  touch "$marker"
  log "extract: done"
}

# One R2 fastq per (library, chunk): TSeq100_L{1,2}_{GEX,CRI}_S{n}_L001_R2_{chunk:03}.fastq.gz
align_one() {
  local mode=$1 fq=$2 out_bam=$3 genes_or_features=$4
  [ -f "$out_bam" ] && return 0
  case "$mode" in
    GEX)
      "$BIN/sbxpress-probeset-align" "$fq" --probe-set "$PROBE_SET" \
        --genes-out "$genes_or_features" -o "$out_bam" \
        --metrics "$out_bam.metrics.json"
      ;;
    CRI)
      "$BIN/sbxpress-guides" "$fq" --feature-reference "$FEATURE_REF" \
        --features-out "$genes_or_features" -o "$out_bam" \
        --metrics "$out_bam.metrics.json"
      ;;
  esac
}
export -f align_one
export BIN PROBE_SET FEATURE_REF

stage_align() {
  local marker="$OUT/aligned/.done"
  if [ -f "$marker" ]; then
    log "align: already done, skipping"
    return
  fi
  mkdir -p "$OUT/aligned"/{TSeq100_L1_GEX,TSeq100_L2_GEX,TSeq100_L1_CRI,TSeq100_L2_CRI}

  local genes_out="$OUT/aligned/genes.tsv"
  local features_out="$OUT/aligned/guide_features.tsv"

  # Warm up the shared, deterministic manifests serially first (both derived
  # purely from the reference files, not from reads -- safe to reuse across
  # every chunk of both libraries in each mode) to avoid a create/verify race
  # once the parallel fan-out starts.
  if [ ! -f "$genes_out" ]; then
    log "align: creating genes manifest (serial warm-up chunk)"
    local first_gex
    first_gex=$(ls "$OUT"/extracted/TSeq100_L1_GEX_S*_R2_001.fastq.gz | head -1)
    align_one GEX "$first_gex" "$OUT/aligned/TSeq100_L1_GEX/$(basename "$first_gex" .fastq.gz).bam" "$genes_out"
  fi
  if [ ! -f "$features_out" ]; then
    log "align: creating features manifest (serial warm-up chunk)"
    local first_cri
    first_cri=$(ls "$OUT"/extracted/TSeq100_L1_CRI_S*_R2_001.fastq.gz | head -1)
    align_one CRI "$first_cri" "$OUT/aligned/TSeq100_L1_CRI/$(basename "$first_cri" .fastq.gz).bam" "$features_out"
  fi

  log "align: fan-out across all chunks, $CORES-way parallel"
  {
    for fq in "$OUT"/extracted/TSeq100_L1_GEX_S*_R2_*.fastq.gz "$OUT"/extracted/TSeq100_L2_GEX_S*_R2_*.fastq.gz; do
      lib=$(basename "$fq" | sed -E 's/_S[0-9]+_.*//')
      echo "GEX|$fq|$OUT/aligned/$lib/$(basename "$fq" .fastq.gz).bam|$genes_out"
    done
    for fq in "$OUT"/extracted/TSeq100_L1_CRI_S*_R2_*.fastq.gz "$OUT"/extracted/TSeq100_L2_CRI_S*_R2_*.fastq.gz; do
      lib=$(basename "$fq" | sed -E 's/_S[0-9]+_.*//')
      echo "CRI|$fq|$OUT/aligned/$lib/$(basename "$fq" .fastq.gz).bam|$features_out"
    done
  } | xargs -P "$CORES" -I{} bash -c '
      IFS="|" read -r mode fq out manifest <<< "{}"
      align_one "$mode" "$fq" "$out" "$manifest"
    ' \
    > "$LOG_DIR/align.log" 2>&1

  log "align: merging per-library chunk BAMs (samtools cat, same deterministic header)"
  for lib in TSeq100_L1_GEX TSeq100_L2_GEX TSeq100_L1_CRI TSeq100_L2_CRI; do
    samtools cat -o "$OUT/aligned/$lib.bam" "$OUT/aligned/$lib"/*.bam
  done
  touch "$marker"
  log "align: done"
}

stage_pipeline() {
  local lib_n=$1  # 1 or 2
  local gex_bam="$OUT/aligned/TSeq100_L${lib_n}_GEX.bam"
  local cri_bam="$OUT/aligned/TSeq100_L${lib_n}_CRI.bam"
  local gex_mol="$OUT/molecules/L${lib_n}_gex.bam"
  local demux_dir="$OUT/demux/L${lib_n}"
  local results_dir="$OUT/results/L${lib_n}"
  local marker="$results_dir/.done"

  if [ -f "$marker" ]; then
    log "pipeline L$lib_n: already done, skipping"
    return
  fi

  if [ ! -f "$gex_mol" ]; then
    log "pipeline L$lib_n: molecules (GEX, normal UMI correction, not --tagseq)"
    # --evidence-workers requires an INDEXED BAM (.bai), which in turn
    # requires a COORDINATE-SORTED BAM -- samtools cat (used to merge the
    # per-chunk align output) only concatenates, it does not sort, so the
    # merged BAM must be sorted first or `samtools index` fails outright
    # ("NO_COOR reads not in a single block", confirmed on the real
    # TSeq100_L1_GEX.bam). Without --evidence-workers, sbx10x-molecules
    # silently falls back to its single-threaded serial path (default
    # --evidence-workers 1), dramatically slower on a multi-billion-read
    # library despite the parallel path (commits 3b7366c/9afc664) being
    # fully present in the binary. ~5,352 real work items at the default
    # --evidence-chunk-bases (53,519 genes / 10Mb chunks) comfortably
    # support high worker counts.
    local gex_bam_sorted="${gex_bam%.bam}.sorted.bam"
    if [ ! -f "$gex_bam_sorted" ]; then
      log "pipeline L$lib_n: sorting aligned BAM by coordinate (required for indexed evidence collection)"
      samtools sort -@ "$CORES" -m 8G -o "$gex_bam_sorted" "$gex_bam"
    fi
    [ -f "$gex_bam_sorted.bai" ] || samtools index -@ "$CORES" "$gex_bam_sorted"
    "$BIN/sbx10x-molecules" "$gex_bam_sorted" -o "$gex_mol" \
      --evidence-workers 64 --threads 16 \
      --metrics "$OUT/molecules/L${lib_n}_gex.metrics.json" \
      > "$LOG_DIR/molecules_L${lib_n}.log" 2>&1
  fi

  mkdir -p "$demux_dir"
  if [ -z "$(ls -A "$demux_dir" 2>/dev/null)" ]; then
    log "pipeline L$lib_n: probe-demux"
    "$BIN/sbxpress-probe-demux" "$gex_mol" --probe-sample-sheet "$PROBE_SAMPLE_SHEET" \
      -o "$demux_dir" > "$LOG_DIR/probe_demux_L${lib_n}.log" 2>&1
  fi

  log "pipeline L$lib_n: per-sample matrix/cells + CRISPR counting/protospacer-call (run_flex_multi.sh)"
  SBXPRESS_LIBEXEC="$BIN" PATH="$BIN:$PATH" \
    /mnt/analysis/software_shared/sbxpress/scripts/run_flex_multi.sh \
    --crispr-bam "$cri_bam" --guide-features "$OUT/aligned/guide_features.tsv" \
    "$gex_mol" "$demux_dir" "$results_dir" \
    > "$LOG_DIR/run_flex_multi_L${lib_n}.log" 2>&1

  touch "$marker"
  log "pipeline L$lib_n: done"
}

log "=== CAT06 full-dataset run starting: $CORES cores, output $OUT ==="
stage_extract
stage_align
stage_pipeline 1
stage_pipeline 2
log "=== CAT06 full-dataset run complete ==="
