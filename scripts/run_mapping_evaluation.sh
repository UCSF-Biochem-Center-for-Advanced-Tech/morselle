#!/usr/bin/env bash
# Run the two-arm CAT02 M2 mapping evaluation described in mapping.md.
#
# Input:  pre-demultiplexed paired Illumina-style R1/R2 FASTQs (sbxpress illumina path)
# Output: two annotated+counted pipeline runs in OUT_DIR/{baseline,multimap}/
#
# Usage:
#   run_mapping_evaluation.sh [--threads N] [--out-dir DIR] [--expect-cells N]
#
# Prerequisites: sbxpress, sbxpress-illumina, samtools, minimap2 on PATH.
# Run tests/validate_cr_bam_to_mtx.sh first to confirm the current binary.
set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
FASTQ_DIR=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/example_fastqs/example_CAT02_fq/M2_GEX
SAMPLE=M2_GEX
WHITELIST=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/737K-arc-v1.txt.gz
MMI=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi
GTF=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/gencode.v44.primary_assembly.basic.annotation.gtf.gz
JUNCTIONS=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/gencode.v44.primary_assembly.junctions.bed
GENES_MTX=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/genes.mtx

OUT_DIR=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/mapping_evaluation
THREADS=16
SORT_MEM=2G
EXPECT_CELLS=

# ── Arg parsing ───────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case $1 in
    --threads)       THREADS=${2:?}; shift 2 ;;
    --out-dir)       OUT_DIR=${2:?}; shift 2 ;;
    --expect-cells)  EXPECT_CELLS=${2:?}; shift 2 ;;
    --help)
      sed -n '1,11p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done

die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"; }
need sbxpress; need samtools

EXTRACT_DIR="$OUT_DIR/extract"
BASELINE_DIR="$OUT_DIR/baseline"
MULTIMAP_DIR="$OUT_DIR/multimap"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$EXTRACT_DIR" "$BASELINE_DIR/outs" "$MULTIMAP_DIR/outs" "$LOG_DIR"

# ── Record provenance ─────────────────────────────────────────────────────────
RUN_INFO="$OUT_DIR/run_info.txt"
{
  echo "run_mapping_evaluation.sh"
  echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "sbxpress: $(sbxpress --version 2>&1 | head -1)"
  echo "sbxpress illumina: (subcommand of sbxpress $(sbxpress --version 2>&1 | head -1))"
  echo "samtools: $(samtools --version | head -1)"
  echo "minimap2: $(sbxpress minimap --version 2>&1 | head -1 || echo unknown)"
  echo ""
  echo "inputs:"
  echo "  fastq_dir:  $FASTQ_DIR"
  echo "  whitelist:  $WHITELIST"
  echo "  mmi:        $MMI"
  echo "  gtf:        $GTF"
  echo "  junctions:  $JUNCTIONS"
  echo "  genes_mtx:  $GENES_MTX"
  echo ""
  echo "checksums:"
  md5sum "$MMI" "$GTF" "$JUNCTIONS" "$WHITELIST" "$GENES_MTX" 2>/dev/null || true
  echo ""
  echo "r2_read_count: $(
    for f in "$FASTQ_DIR"/${SAMPLE}*_R2_*.fastq.gz; do
      zcat "$f"
    done | awk 'NR%4==1' | wc -l
  )"
  echo "threads: $THREADS"
} | tee "$RUN_INFO"

# ── Collect R1 and R2 file lists ──────────────────────────────────────────────
mapfile -t R1_FILES < <(ls "$FASTQ_DIR"/${SAMPLE}*_R1_*.fastq.gz | sort)
mapfile -t R2_FILES < <(ls "$FASTQ_DIR"/${SAMPLE}*_R2_*.fastq.gz | sort)
(( ${#R1_FILES[@]} > 0 )) || die "no R1 FASTQs found in $FASTQ_DIR for $SAMPLE"
(( ${#R1_FILES[@]} == ${#R2_FILES[@]} )) || die "R1/R2 count mismatch: ${#R1_FILES[@]} vs ${#R2_FILES[@]}"
echo ""
echo "R1 files: ${#R1_FILES[@]}"
echo "R2 files: ${#R2_FILES[@]}"

# ── Step 1: sbxpress illumina — extract CB/UMI into R2 FASTQ comment ─────────
# sbxpress illumina outputs one shard per input pair; use metrics JSON as done-check.
ILLUMINA_METRICS="$EXTRACT_DIR/${SAMPLE}.illumina.metrics.json"
TAGGED_R2_BASE="$EXTRACT_DIR/${SAMPLE}_tagged_r2.fastq.gz"
if [[ ! -s "$ILLUMINA_METRICS" ]]; then
  echo ""
  echo "[1/5] converting R1/R2 pairs to minimap2-tagged R2 FASTQ"
  R1_ARGS=(); for f in "${R1_FILES[@]}"; do R1_ARGS+=(--r1 "$f"); done
  R2_ARGS=(); for f in "${R2_FILES[@]}"; do R2_ARGS+=(--r2 "$f"); done
  sbxpress illumina \
    "${R1_ARGS[@]}" "${R2_ARGS[@]}" \
    --whitelist "$WHITELIST" \
    --output "$TAGGED_R2_BASE" \
    --metrics "$ILLUMINA_METRICS" \
    --threads "$THREADS" \
    2>&1 | tee "$LOG_DIR/illumina.log"
else
  echo "[1/5] tagged R2 already present; reusing shards in $EXTRACT_DIR"
fi

# Collect all shards (may be one file or many, depending on thread count)
shopt -s nullglob
mapfile -t TAGGED_R2_FILES < <(ls "$EXTRACT_DIR/${SAMPLE}_tagged_r2"*.fastq.gz | sort)
shopt -u nullglob
(( ${#TAGGED_R2_FILES[@]} > 0 )) || die "no tagged R2 shards found in $EXTRACT_DIR"
echo "  tagged R2 shards: ${#TAGGED_R2_FILES[@]}"

# ── Helper: run annotation, sort, molecules, matrix, and cells for one arm ────
run_arm() {
  local label=$1        # "baseline" or "multimap"
  local arm_dir=$2
  local map_args=("${@:3}")

  local annotated_bam="$arm_dir/${SAMPLE}.annotated.bam"
  local molecule_bam="$arm_dir/${SAMPLE}.molecules.bam"
  local raw_h5="$arm_dir/outs/raw_feature_bc_matrix.h5"
  local filtered_h5="$arm_dir/outs/filtered_feature_bc_matrix.h5"

  mkdir -p "$arm_dir/outs"

  if [[ ! -s "${annotated_bam}.bai" ]]; then
    echo "  mapping, annotation, and coordinate sort"
    if [[ "$label" == "multimap" ]]; then
      local ann_args=(--gtf "$GTF" --output - \
        --metrics "$arm_dir/${SAMPLE}.annotation.metrics.json" \
        --genes-mtx "$GENES_MTX" --multimap-tags --three-prime)
    else
      local ann_args=(--gtf "$GTF" --output - \
        --metrics "$arm_dir/${SAMPLE}.annotation.metrics.json" \
        --genes-mtx "$GENES_MTX" --unordered --three-prime)
    fi
    sbxpress minimap \
      -a -x splice:sr --junc-bed "$JUNCTIONS" "${map_args[@]}" -y \
      -t "$THREADS" "$MMI" "${TAGGED_R2_FILES[@]}" \
      2>"$LOG_DIR/${label}.minimap2.log" \
    | sbxpress annotate - "${ann_args[@]}" \
      2>"$LOG_DIR/${label}.annotate.log" \
    | samtools sort -@ "$THREADS" -m "$SORT_MEM" \
      -T "$arm_dir/tmp" -o "$annotated_bam" -
    samtools index -@ "$THREADS" "$annotated_bam"
    samtools quickcheck -v "$annotated_bam"
  else
    echo "  annotated BAM already present; reusing"
  fi

  if [[ ! -s "${molecule_bam}.bai" ]]; then
    echo "  molecule calling"
    sbxpress molecules "$annotated_bam" \
      --output "$molecule_bam" \
      --metrics "$arm_dir/${SAMPLE}.molecules.metrics.json" \
      2>"$LOG_DIR/${label}.molecules.log"
    samtools index -@ "$THREADS" "$molecule_bam"
    samtools quickcheck -v "$molecule_bam"
  else
    echo "  molecule BAM already present; reusing"
  fi

  if [[ ! -s "$raw_h5" ]]; then
    echo "  matrix construction"
    sbxpress matrix "$molecule_bam" \
      --output "$arm_dir/outs" \
      --features "$GENES_MTX" \
      --cell-barcodes "$WHITELIST" \
      2>"$LOG_DIR/${label}.matrix.log"
  else
    echo "  raw matrix already present; reusing"
  fi

  if [[ ! -s "$filtered_h5" ]]; then
    echo "  cell calling"
    local cell_args=()
    [[ -z "$EXPECT_CELLS" ]] || cell_args+=(--expect-cells "$EXPECT_CELLS")
    sbxpress cells "$raw_h5" --output "$arm_dir/outs" "${cell_args[@]}" \
      2>"$LOG_DIR/${label}.cells.log"
  else
    echo "  filtered matrix already present; reusing"
  fi
}

# ── Step 2a: production baseline (--secondary=no) ─────────────────────────────
echo ""
echo "[2/5] baseline arm (--secondary=no, --unordered)"
run_arm baseline "$BASELINE_DIR" --secondary=no

# ── Step 2b: diagnostic multimap (-N 10, --multimap-tags) ────────────────────
echo ""
echo "[3/5] multimap arm (-N 10, --multimap-tags)"
run_arm multimap "$MULTIMAP_DIR" -N 10

# ── Step 3: brief inline summary of mapping metrics ──────────────────────────
echo ""
echo "[4/5] annotation metrics summary"
for arm in baseline multimap; do
  arm_dir="$OUT_DIR/$arm"
  metrics="$arm_dir/${SAMPLE}.annotation.metrics.json"
  if [[ -s "$metrics" ]]; then
    echo "  $arm:"
    python3 -c "
import json, sys
m = json.load(open('$metrics'))
for k,v in sorted(m.items()):
    print(f'    {k}: {v}')
" 2>/dev/null || cat "$metrics"
  fi
done

# ── Step 4: molecule metrics summary ─────────────────────────────────────────
echo ""
echo "[5/5] molecule metrics summary"
for arm in baseline multimap; do
  arm_dir="$OUT_DIR/$arm"
  metrics="$arm_dir/${SAMPLE}.molecules.metrics.json"
  if [[ -s "$metrics" ]]; then
    echo "  $arm:"
    python3 -c "
import json
m = json.load(open('$metrics'))
for k,v in sorted(m.items()):
    print(f'    {k}: {v}')
" 2>/dev/null || cat "$metrics"
  fi
done

echo ""
echo "Done. Outputs in: $OUT_DIR"
echo "Next: python3 scripts/analyze_mapping.py --bam $MULTIMAP_DIR/${SAMPLE}.annotated.bam --out $OUT_DIR/analysis"
