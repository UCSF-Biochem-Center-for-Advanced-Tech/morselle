#!/usr/bin/env bash
# sbxpress-pipeline.sh — Unified multi-mode driver script for sbxpress.
#
# Supports 2 Input Modes:
#   --input-mode {sbx|illumina}  (default: sbx)
#
# Supports 3 Analysis Modes:
#   --analysis-mode {gex|flex|tagseq}  (default: gex)
#
# This script is modular, restartable, and handles extraction/tagging,
# mapping/probe-alignment, molecule deduplication, feature matrix creation,
# cell calling, and optional auditing.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  sbxpress-pipeline.sh \
    --input-mode {sbx|illumina} \
    --analysis-mode {gex|flex|tagseq} \
    --out-dir OUTPUT_DIRECTORY \
    [mode-specific inputs and options]

Modes:
  --input-mode MODE        Data input mode: 'sbx' (default) or 'illumina'.
  --analysis-mode MODE     Analysis mode: 'gex' (default), 'flex', or 'tagseq'.

Input Options:
  --fastqs DIR [...]       (SBX mode) One or more raw SBX device output directories.
  --bclconvert-samplesheet FILE (SBX mode) Sample sheet for demultiplexing.
  --r1 FASTQ [...]         (Illumina mode) Read 1 FASTQ file(s).
  --r2 FASTQ [...]         (Illumina mode) Read 2 FASTQ file(s).
  --sample NAME            Sample ID to process.
  --whitelist FILE         10x cell-barcode whitelist file (.txt or .txt.gz).

Analysis-Specific Options (GEX / Tagseq):
  --mmi FILE               minimap2 splice index (.mmi).
  --gtf FILE               GENCODE-compatible annotation GTF (.gtf or .gtf.gz).
  --junctions FILE         Junction BED file.
  --genome NAME            Genome feature label; default 'GRCh38'.
  --three-prime / --no-three-prime
                           Enable/disable 3' terminal preference tags (default: enabled).

Analysis-Specific Options (Flex):
  --probe-set FILE         10x Flex probe-set reference CSV.
  --probe-barcode-wl FILE  10x Flex probe barcode whitelist.
  --config FILE            Flex sample sheet config.csv (sample_id -> probe_barcode_ids).

Analysis-Specific Options (CAT Tagseq):
  --tagseq-barcodes FILE   Exact barcode whitelist for cell calling (bypasses EmptyDrops).

Performance & Threading Options:
  --threads N              Global default worker threads (default: 16).
  --extract-threads N      Extraction threads (default: 32).
  --map-threads N          minimap2 threads (default: 32).
  --annotate-threads N     Annotation threads (default: 8).
  --sort-threads N         samtools sort threads (default: 8).
  --molecule-threads N     Molecule calling threads (default: 8).
  --align-threads N        Probeset-align threads (default: 16).
  --split N                Extracted shards (default: 16).
  --sort-mem SIZE          samtools sort memory per thread (default: 1G).
  --audit                  Run bam-audit post-processing pass.
  --help                   Show this message.

Output layout:
  OUT_DIR/outs/raw_feature_bc_matrix.h5
  OUT_DIR/outs/filtered_feature_bc_matrix.h5
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
info() { echo "[sbxpress-pipeline] $*"; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"; }
require_path() { [[ -e "$1" ]] || die "required input file/directory not found: $1"; }

# Default parameters
INPUT_MODE="sbx"
ANALYSIS_MODE="gex"
FASTQS=()
R1_FASTQS=()
R2_FASTQS=()
SAMPLE_SHEET=""
SAMPLE=""
WHITELIST=""
OUT_DIR=""

# GEX / Tagseq parameters
MMI=""
GTF=""
JUNCTIONS=""
GENOME="GRCh38"
THREE_PRIME=1

# Flex parameters
PROBE_SET=""
PROBE_BARCODE_WL=""
CONFIG_CSV=""

# Tagseq parameters
TAGSEQ_BARCODES=""

# Threading & Memory
THREADS=16
EXTRACT_THREADS=32
MAP_THREADS=32
ANNOTATE_THREADS=8
SORT_THREADS=8
MOLECULE_THREADS=8
ALIGN_THREADS=16
SPLIT=16
SORT_MEM="1G"
AUDIT=0

# Parse arguments
while (( $# > 0 )); do
    case "$1" in
        --input-mode) INPUT_MODE=${2:?missing value for --input-mode}; shift 2 ;;
        --analysis-mode) ANALYSIS_MODE=${2:?missing value for --analysis-mode}; shift 2 ;;
        --fastqs)
            shift
            while (( $# > 0 )) && [[ $1 != --* ]]; do
                FASTQS+=("$1"); shift
            done
            ;;
        --r1)
            shift
            while (( $# > 0 )) && [[ $1 != --* ]]; do
                R1_FASTQS+=("$1"); shift
            done
            ;;
        --r2)
            shift
            while (( $# > 0 )) && [[ $1 != --* ]]; do
                R2_FASTQS+=("$1"); shift
            done
            ;;
        --bclconvert-samplesheet) SAMPLE_SHEET=${2:?missing value for --bclconvert-samplesheet}; shift 2 ;;
        --sample) SAMPLE=${2:?missing value for --sample}; shift 2 ;;
        --whitelist) WHITELIST=${2:?missing value for --whitelist}; shift 2 ;;
        --out-dir) OUT_DIR=${2:?missing value for --out-dir}; shift 2 ;;
        --mmi) MMI=${2:?missing value for --mmi}; shift 2 ;;
        --gtf) GTF=${2:?missing value for --gtf}; shift 2 ;;
        --junctions) JUNCTIONS=${2:?missing value for --junctions}; shift 2 ;;
        --genome) GENOME=${2:?missing value for --genome}; shift 2 ;;
        --three-prime) THREE_PRIME=1; shift ;;
        --no-three-prime) THREE_PRIME=0; shift ;;
        --probe-set) PROBE_SET=${2:?missing value for --probe-set}; shift 2 ;;
        --probe-barcode-wl) PROBE_BARCODE_WL=${2:?missing value for --probe-barcode-wl}; shift 2 ;;
        --config) CONFIG_CSV=${2:?missing value for --config}; shift 2 ;;
        --tagseq-barcodes) TAGSEQ_BARCODES=${2:?missing value for --tagseq-barcodes}; shift 2 ;;
        --threads) THREADS=${2:?missing value for --threads}; shift 2 ;;
        --extract-threads) EXTRACT_THREADS=${2:?missing value for --extract-threads}; shift 2 ;;
        --map-threads) MAP_THREADS=${2:?missing value for --map-threads}; shift 2 ;;
        --annotate-threads) ANNOTATE_THREADS=${2:?missing value for --annotate-threads}; shift 2 ;;
        --sort-threads) SORT_THREADS=${2:?missing value for --sort-threads}; shift 2 ;;
        --molecule-threads) MOLECULE_THREADS=${2:?missing value for --molecule-threads}; shift 2 ;;
        --align-threads) ALIGN_THREADS=${2:?missing value for --align-threads}; shift 2 ;;
        --split) SPLIT=${2:?missing value for --split}; shift 2 ;;
        --sort-mem) SORT_MEM=${2:?missing value for --sort-mem}; shift 2 ;;
        --audit) AUDIT=1; shift ;;
        --help|-h) usage; exit 0 ;;
        *) die "unknown argument: $1" ;;
    esac
done

# Validate mode selections
[[ "$INPUT_MODE" =~ ^(sbx|illumina)$ ]] || die "invalid --input-mode: '$INPUT_MODE' (must be 'sbx' or 'illumina')"
[[ "$ANALYSIS_MODE" =~ ^(gex|flex|tagseq)$ ]] || die "invalid --analysis-mode: '$ANALYSIS_MODE' (must be 'gex', 'flex', or 'tagseq')"
[[ -n "$OUT_DIR" ]] || die "missing required parameter: --out-dir"
[[ -n "$WHITELIST" ]] || die "missing required parameter: --whitelist"

require_command sbxpress
require_command samtools

# Validate mode-specific prerequisites
if [[ "$INPUT_MODE" == "sbx" ]]; then
    ((${#FASTQS[@]} > 0)) || die "--fastqs is required for SBX input mode"
    [[ -n "$SAMPLE_SHEET" ]] || die "--bclconvert-samplesheet is required for SBX input mode"
    [[ -n "$SAMPLE" ]] || die "--sample is required for SBX input mode"
    require_path "$SAMPLE_SHEET"
elif [[ "$INPUT_MODE" == "illumina" ]]; ((${#R1_FASTQS[@]} > 0 && ${#R2_FASTQS[@]} > 0)) || die "--r1 and --r2 are required for Illumina input mode"; fi

require_path "$WHITELIST"

if [[ "$ANALYSIS_MODE" == "gex" || "$ANALYSIS_MODE" == "tagseq" ]]; then
    [[ -n "$MMI" && -n "$GTF" && -n "$JUNCTIONS" ]] || die "--mmi, --gtf, and --junctions are required for $ANALYSIS_MODE analysis mode"
    require_path "$MMI"
    require_path "$GTF"
    require_path "$JUNCTIONS"
elif [[ "$ANALYSIS_MODE" == "flex" ]]; then
    [[ -n "$PROBE_SET" ]] || die "--probe-set is required for Flex analysis mode"
    require_path "$PROBE_SET"
fi

# Directory structure setup
SAMPLE_NAME="${SAMPLE:-sample}"
EXTRACT_DIR="$OUT_DIR/extract"
ANNOTATION_DIR="$OUT_DIR/$SAMPLE_NAME/annotation"
MOLECULES_DIR="$OUT_DIR/$SAMPLE_NAME/molecules"
OUTS_DIR="$OUT_DIR/$SAMPLE_NAME/outs"
AUDIT_DIR="$OUT_DIR/$SAMPLE_NAME/audit"
GENES_MTX="$OUT_DIR/reference/genes.mtx"

mkdir -p "$EXTRACT_DIR" "$ANNOTATION_DIR" "$MOLECULES_DIR" "$OUTS_DIR" "$OUT_DIR/reference"

info "Starting pipeline run"
info "Input Mode: $INPUT_MODE | Analysis Mode: $ANALYSIS_MODE | Sample: $SAMPLE_NAME"
info "Output Directory: $OUT_DIR"

# Stage 1: Extraction / Tagging
TAGGED_R2="$EXTRACT_DIR/${SAMPLE_NAME}.tagged.R2.fastq.gz"

if [[ ! -s "$TAGGED_R2" ]]; then
    info "Stage 1: Extracting and tagging reads..."
    if [[ "$INPUT_MODE" == "sbx" ]]; then
        if [[ "$ANALYSIS_MODE" == "flex" ]]; then
            sbxpress extract \
                --fastqs "${FASTQS[@]}" \
                --bclconvert-samplesheet "$SAMPLE_SHEET" \
                --whitelist "$WHITELIST" \
                --out-dir "$EXTRACT_DIR" \
                --threads "$EXTRACT_THREADS" \
                --split "$SPLIT" \
                --sample "$SAMPLE_NAME"
        else
            sbxpress extract \
                --fastqs "${FASTQS[@]}" \
                --bclconvert-samplesheet "$SAMPLE_SHEET" \
                --whitelist "$WHITELIST" \
                --out-dir "$EXTRACT_DIR" \
                --threads "$EXTRACT_THREADS" \
                --split "$SPLIT" \
                --sample "$SAMPLE_NAME" \
                --minimap2
        fi
        # Concatenate extracted FASTQ shards for the sample if multiple shards exist
        cat "$EXTRACT_DIR/$SAMPLE_NAME"/*.fastq.gz > "$TAGGED_R2" 2>/dev/null || \
        cp "$EXTRACT_DIR/$SAMPLE_NAME.fastq.gz" "$TAGGED_R2" 2>/dev/null || \
        find "$EXTRACT_DIR" -name "*.fastq.gz" -exec cat {} + > "$TAGGED_R2"
    else
        # Illumina PE Mode
        ILLUMINA_ARGS=()
        for r1 in "${R1_FASTQS[@]}"; do ILLUMINA_ARGS+=(--r1 "$r1"); done
        for r2 in "${R2_FASTQS[@]}"; do ILLUMINA_ARGS+=(--r2 "$r2"); done
        
        sbxpress illumina \
            "${ILLUMINA_ARGS[@]}" \
            --whitelist "$WHITELIST" \
            --output "$TAGGED_R2" \
            --threads "$EXTRACT_THREADS"
    fi
else
    info "Stage 1: Extracted tagged R2 FASTQ already exists ($TAGGED_R2), skipping."
fi

# Stage 2: Alignment & Annotation
ANNOTATED_BAM="$ANNOTATION_DIR/${SAMPLE_NAME}.annotated.bam"

if [[ ! -s "$ANNOTATED_BAM" ]]; then
    info "Stage 2: Running alignment and feature annotation..."
    if [[ "$ANALYSIS_MODE" == "flex" ]]; then
        # Optional Probe Demultiplexing
        if [[ -n "$PROBE_BARCODE_WL" && -n "$CONFIG_CSV" ]]; then
            info "Running probe demultiplexing..."
            sbxpress probe-demux \
                --r2-fastq "$TAGGED_R2" \
                --probe-barcode-wl "$PROBE_BARCODE_WL" \
                --config "$CONFIG_CSV" \
                --output "$OUT_DIR/$SAMPLE_NAME/demux"
        fi

        # Probe-Set Alignment (using multi-threaded sbxpress probeset-align)
        sbxpress probeset-align \
            "$TAGGED_R2" \
            --probe-set "$PROBE_SET" \
            --genes-out "$GENES_MTX" \
            --threads "$ALIGN_THREADS" \
            --output "$ANNOTATED_BAM"
    else
        # 3' GEX / Tagseq genomic splice mapping
        THREE_PRIME_FLAG=""
        [[ "$THREE_PRIME" -eq 1 ]] && THREE_PRIME_FLAG="--three-prime"

        sbxpress minimap -a -x splice:sr \
            --junc-bed "$JUNCTIONS" \
            --secondary=no -y -t "$MAP_THREADS" \
            "$MMI" \
            "$TAGGED_R2" \
        | sbxpress annotate - \
            --gtf "$GTF" \
            --output - \
            --genes-mtx "$GENES_MTX" \
            $THREE_PRIME_FLAG \
            --unordered \
            --threads "$ANNOTATE_THREADS" \
        | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" -T "$ANNOTATION_DIR/tmp" \
            -o "$ANNOTATED_BAM" -

        samtools index -@ "$SORT_THREADS" "$ANNOTATED_BAM"
    fi
else
    info "Stage 2: Annotated BAM already exists ($ANNOTATED_BAM), skipping."
fi

# Stage 3: Molecule Deduplication
MOLECULES_BAM="$MOLECULES_DIR/${SAMPLE_NAME}.molecules.bam"

if [[ ! -s "$MOLECULES_BAM" ]]; then
    info "Stage 3: Calling UMIs and resolving molecules..."
    TAGSEQ_FLAG=""
    [[ "$ANALYSIS_MODE" == "tagseq" ]] && TAGSEQ_FLAG="--tagseq"

    sbxpress molecules "$ANNOTATED_BAM" \
        $TAGSEQ_FLAG \
        --output "$MOLECULES_BAM" \
        --threads "$MOLECULE_THREADS"

    samtools index -@ "$SORT_THREADS" "$MOLECULES_BAM"
else
    info "Stage 3: Molecules BAM already exists ($MOLECULES_BAM), skipping."
fi

# Stage 4: Matrix Generation & Cell Calling
RAW_MATRIX_H5="$OUTS_DIR/raw_feature_bc_matrix.h5"

if [[ ! -s "$RAW_MATRIX_H5" ]]; then
    info "Stage 4a: Generating raw feature-barcode expression matrix..."
    sbxpress matrix "$MOLECULES_BAM" \
        --output "$OUTS_DIR" \
        --genome "$GENOME" \
        --features "$GENES_MTX" \
        --cell-barcodes "$WHITELIST"
fi

info "Stage 4b: Calling cells and writing filtered matrix..."
if [[ -n "$TAGSEQ_BARCODES" ]]; then
    info "Enforcing exact cell barcode whitelist (--tagseq-barcodes $TAGSEQ_BARCODES)..."
    sbxpress cells "$RAW_MATRIX_H5" \
        --tagseq-barcodes "$TAGSEQ_BARCODES" \
        --output "$OUTS_DIR"
else
    sbxpress cells "$RAW_MATRIX_H5" \
        --output "$OUTS_DIR"
fi

# Stage 5: Audit (Optional)
if [[ "$AUDIT" -eq 1 ]]; then
    info "Stage 5: Running audit pass..."
    mkdir -p "$AUDIT_DIR"
    sbxpress audit "$MOLECULES_BAM" \
        --label "${SAMPLE_NAME}_${INPUT_MODE}_${ANALYSIS_MODE}" \
        --cells "$OUTS_DIR/filtered_feature_bc_matrix/barcodes.tsv.gz" \
        --output "$AUDIT_DIR/${SAMPLE_NAME}.audit.json" \
        --threads "$MOLECULE_THREADS"
fi

info "Pipeline complete! Outputs available in: $OUTS_DIR"
