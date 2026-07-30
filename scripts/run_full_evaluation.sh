#!/usr/bin/env bash
# Two-arm mapping evaluation on the full P5 and P7 CAT02 M2_GEX SBX runs.
# Mirrors run_mapping_evaluation.sh but operates on the pre-extracted paired
# R1/R2 FASTQs in p5_out_7-7/ and p7_out_7-7/. These are large (16 shards
# each, ~8 GB per shard) — expect 12-24 h per orientation on 16 threads.
#
# Usage:
#   run_full_evaluation.sh [--p5-only | --p7-only] [--threads N] [--out-dir DIR]
#
# Restartable: each stage checks for its output before running.
set -euo pipefail

# ── Paths ─────────────────────────────────────────────────────────────────────
P5_FASTQ_DIR=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/p5_out_7-7
P7_FASTQ_DIR=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/p7_out_7-7
SAMPLE=M2_GEX
WHITELIST=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/737K-arc-v1.txt.gz
MMI=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/GRCh38.primary_assembly.gencode-v44.splice-sr.mmi
GTF=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/gencode.v44.primary_assembly.basic.annotation.gtf.gz
JUNCTIONS=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/gencode.v44.primary_assembly.junctions.bed
GENES_MTX=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/resources/genes.mtx

OUT_BASE=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/mapping_evaluation
THREADS=16
SORT_MEM=4G

RUN_P5=1
RUN_P7=1

while [[ $# -gt 0 ]]; do
  case $1 in
    --p5-only)  RUN_P7=0; shift ;;
    --p7-only)  RUN_P5=0; shift ;;
    --threads)  THREADS=${2:?}; shift 2 ;;
    --out-dir)  OUT_BASE=${2:?}; shift 2 ;;
    --help)
      sed -n '1,10p' "$0" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 1 ;;
  esac
done

die()  { echo "ERROR: $*" >&2; exit 1; }
need() { command -v "$1" >/dev/null 2>&1 || die "required: $1"; }
need sbxpress; need samtools

# ── Per-orientation pipeline ───────────────────────────────────────────────────
run_orientation() {
  local orient=$1      # p5 or p7
  local fastq_dir=$2
  local out_dir="$OUT_BASE/$orient"
  local log_dir="$out_dir/logs"
  local extract_dir="$out_dir/extract"

  mkdir -p "$extract_dir" "$out_dir/baseline/outs" "$out_dir/multimap/outs" "$log_dir"

  echo ""
  echo "════════════════════════════════════════════════"
  echo "  $orient  ($fastq_dir)"
  echo "════════════════════════════════════════════════"

  # Provenance
  {
    echo "run_full_evaluation.sh  orient=$orient"
    echo "date: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "sbxpress: $(sbxpress --version 2>&1 | head -1)"
    echo "samtools: $(samtools --version | head -1)"
    echo "whitelist_md5: $(md5sum "$WHITELIST" | cut -d' ' -f1)"
    echo "mmi_md5:       $(md5sum "$MMI"       | cut -d' ' -f1)"
    echo "gtf_md5:       $(md5sum "$GTF"       | cut -d' ' -f1)"
    echo "threads: $THREADS"
  } > "$out_dir/run_info.txt"

  # ── illumina extraction ──────────────────────────────────────────────────────
  local ill_metrics="$extract_dir/${SAMPLE}.illumina.metrics.json"
  local tagged_r2_base="$extract_dir/${SAMPLE}_tagged_r2.fastq.gz"

  if [[ ! -s "$ill_metrics" ]]; then
    echo "[1/4] sbxpress illumina — tagging CB/UMI into R2 FASTQ comment"
    mapfile -t R1_FILES < <(ls "$fastq_dir"/${SAMPLE}*_R1_*.fastq.gz | sort)
    mapfile -t R2_FILES < <(ls "$fastq_dir"/${SAMPLE}*_R2_*.fastq.gz | sort)
    (( ${#R1_FILES[@]} > 0 )) || die "no R1 FASTQs in $fastq_dir for $SAMPLE"
    (( ${#R1_FILES[@]} == ${#R2_FILES[@]} )) || die "R1/R2 count mismatch"
    echo "  ${#R1_FILES[@]} R1/R2 pairs"
    R1_ARGS=(); for f in "${R1_FILES[@]}"; do R1_ARGS+=(--r1 "$f"); done
    R2_ARGS=(); for f in "${R2_FILES[@]}"; do R2_ARGS+=(--r2 "$f"); done
    sbxpress illumina \
      "${R1_ARGS[@]}" "${R2_ARGS[@]}" \
      --whitelist "$WHITELIST" \
      --output "$tagged_r2_base" \
      --metrics "$ill_metrics" \
      --threads "$THREADS" \
      2>&1 | tee "$log_dir/illumina.log"
  else
    echo "[1/4] illumina extraction already done"
  fi

  shopt -s nullglob
  mapfile -t TAGGED_R2_FILES < <(ls "$extract_dir/${SAMPLE}_tagged_r2"*.fastq.gz | sort)
  shopt -u nullglob
  (( ${#TAGGED_R2_FILES[@]} > 0 )) || die "no tagged R2 shards in $extract_dir"
  echo "  tagged R2 shards: ${#TAGGED_R2_FILES[@]}"

  # ── Map, annotate, sort, downstream — for one arm ───────────────────────────
  run_arm() {
    local arm=$1       # baseline or multimap
    local arm_dir="$out_dir/$arm"
    shift
    local map_args=("$@")

    local annotated_bam="$arm_dir/${SAMPLE}.annotated.bam"
    local molecule_bam="$arm_dir/${SAMPLE}.molecules.bam"
    local raw_h5="$arm_dir/outs/raw_feature_bc_matrix.h5"
    local filtered_h5="$arm_dir/outs/filtered_feature_bc_matrix.h5"

    echo ""
    echo "  ── $arm arm ──"

    if [[ ! -s "${annotated_bam}.bai" ]]; then
      echo "  mapping + annotation + sort  ($(date -u +%H:%M:%SZ))"
      if [[ "$arm" == "multimap" ]]; then
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
        2>"$log_dir/${arm}.minimap2.log" \
      | sbxpress annotate - "${ann_args[@]}" \
        2>"$log_dir/${arm}.annotate.log" \
      | samtools sort -@ "$THREADS" -m "$SORT_MEM" \
        -T "$arm_dir/tmp" -o "$annotated_bam" -
      samtools index -@ "$THREADS" "$annotated_bam"
      samtools quickcheck -v "$annotated_bam"
      echo "  annotated BAM done  ($(date -u +%H:%M:%SZ))"
    else
      echo "  annotated BAM already present"
    fi

    if [[ ! -s "${molecule_bam}.bai" ]]; then
      echo "  molecules  ($(date -u +%H:%M:%SZ))"
      sbxpress molecules "$annotated_bam" \
        --output "$molecule_bam" \
        --metrics "$arm_dir/${SAMPLE}.molecules.metrics.json" \
        2>"$log_dir/${arm}.molecules.log"
      samtools index -@ "$THREADS" "$molecule_bam"
    else
      echo "  molecule BAM already present"
    fi

    if [[ ! -s "$raw_h5" ]]; then
      echo "  matrix  ($(date -u +%H:%M:%SZ))"
      sbxpress matrix "$molecule_bam" \
        --output "$arm_dir/outs" \
        --features "$GENES_MTX" \
        --cell-barcodes "$WHITELIST" \
        2>"$log_dir/${arm}.matrix.log"
    else
      echo "  raw matrix already present"
    fi

    if [[ ! -s "$filtered_h5" ]]; then
      echo "  cells  ($(date -u +%H:%M:%SZ))"
      sbxpress cells "$raw_h5" --output "$arm_dir/outs" \
        2>"$log_dir/${arm}.cells.log"
    else
      echo "  filtered matrix already present"
    fi

    echo "  $arm arm complete  ($(date -u +%H:%M:%SZ))"
  }

  echo ""
  echo "[2/4] baseline arm  ($(date -u +%H:%M:%SZ))"
  run_arm baseline --secondary=no

  echo ""
  echo "[3/4] multimap arm  ($(date -u +%H:%M:%SZ))"
  run_arm multimap -N 10

  echo ""
  echo "[4/4] $orient complete  ($(date -u +%H:%M:%SZ))"
  echo "  baseline annotation: $out_dir/baseline/${SAMPLE}.annotation.metrics.json"
  echo "  multimap annotation: $out_dir/multimap/${SAMPLE}.annotation.metrics.json"
  echo ""
  echo "  next:"
  echo "  python3 scripts/analyze_mapping.py \\"
  echo "    --bam $out_dir/multimap/${SAMPLE}.annotated.bam \\"
  echo "    --baseline $out_dir/baseline/${SAMPLE}.annotated.bam \\"
  echo "    --out $out_dir/analysis --threads $THREADS"
}

echo "full P5/P7 mapping evaluation"
echo "started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "threads: $THREADS"

(( RUN_P5 )) && run_orientation p5 "$P5_FASTQ_DIR"
(( RUN_P7 )) && run_orientation p7 "$P7_FASTQ_DIR"

echo ""
echo "all orientations complete: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
