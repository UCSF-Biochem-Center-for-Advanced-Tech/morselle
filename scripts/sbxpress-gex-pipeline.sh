#!/usr/bin/env bash
# Run raw Roche SBX GEX data through extraction, minimap2, annotation, UMI
# correction, matrix creation, and 10x-style cell calling.
#
# This script deliberately keeps extraction and the two coordinate-sorted BAMs:
# they make a failed stage restartable without repeating raw-data processing or
# mapping. The audit pass is optional because it is not needed to create counts.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat <<'EOF'
Usage:
  sbxpress-gex-pipeline.sh \
    --fastqs DEVICE_0_OUT DEVICE_1_OUT [...] \
    --bclconvert-samplesheet SAMPLES.csv --sample SAMPLE_ID \
    --whitelist 10X_BARCODES.txt.gz \
    --mmi REFERENCE.mmi --gtf ANNOTATION.gtf.gz --junctions JUNCTIONS.bed \
    --out-dir OUTPUT_DIRECTORY [options]

Required inputs:
  --fastqs DIR [...]       Raw SBX FASTQ directories, one or more device outputs.
  --bclconvert-samplesheet FILE      BCL Convert-compatible sample sheet.
  --sample NAME            Demultiplexed sample to map and count.
  --whitelist FILE         10x cell-barcode whitelist.
  --mmi FILE               minimap2 splice index.
  --gtf FILE               GENCODE-compatible annotation GTF.
  --junctions FILE         Junction BED used when the minimap index was built.
  --out-dir DIR            Pipeline output directory.

Options:
  --genes-mtx FILE         Canonical feature manifest; default OUT/reference/genes.mtx.
  --genome NAME            Feature genome label; default GRCh38.
  --extract-threads N      Extraction workers; default 32.
  --map-threads N          minimap2 workers; default 32.
  --annotate-threads N     Annotation workers; default 8. Used by the fast
                           unordered mode; multimap-tag mode is streaming.
  --sort-threads N         samtools sort/index workers; default 8.
  --molecule-threads N     Molecule-calling workers; default 8.
  --split N                Extracted FASTQ shards per sample; default 16.
  --sort-mem SIZE          samtools memory per sort worker; default 1G.
  --tmp-dir DIR            Sort temporary directory; default OUT/tmp.
  --expect-cells N         Expected recovered cells; passed to sbxpress cells.
  --no-three-prime         Disable terminal-exon preference tags.
  --multimap-tags          Retain up to ten secondary hits and add NH/HI/mm tags.
                           This uses query-grouped streaming annotation, not
                           unordered parallel annotation.
  --audit                  Also emit comparison audit JSON and molecule sample.
  --audit-label TEXT       Audit label; default SBX_SAMPLE.
  --molecule-sample N      Audit 1-in-N molecule sample; default 64.
  --barcode-metrics        Also write per-barcode .barcodes.tsv.gz tables next to
                           each stage's metrics (adds run time and memory).
  --help                   Show this help.

The output counts are OUT/SAMPLE/outs/raw_feature_bc_matrix.h5 and
OUT/SAMPLE/outs/filtered_feature_bc_matrix.h5, with matching MEX directories.
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"; }
require_path() { [[ -e "$1" ]] || die "required input not found: $1"; }

FASTQS=()
SAMPLE_SHEET=
SAMPLE=
WHITELIST=
MMI=
GTF=
JUNCTIONS=
OUT_DIR=
GENES_MTX=
GENOME=GRCh38
EXTRACT_THREADS=32
MAP_THREADS=32
ANNOTATE_THREADS=8
SORT_THREADS=8
MOLECULE_THREADS=8
SPLIT=16
SORT_MEM=1G
TMP_DIR=
EXPECT_CELLS=
THREE_PRIME=1
MULTIMAP_TAGS=0
AUDIT=0
BARCODE_METRICS=()
AUDIT_LABEL=
MOLECULE_SAMPLE=64

while (( $# > 0 )); do
    case "$1" in
        --fastqs)
            shift
            (( $# > 0 )) || die "missing value for --fastqs"
            while (( $# > 0 )) && [[ $1 != --* ]]; do
                FASTQS+=("$1")
                shift
            done
            ;;
        --bclconvert-samplesheet) SAMPLE_SHEET=${2:?missing value for --bclconvert-samplesheet}; shift 2 ;;
        --sample) SAMPLE=${2:?missing value for --sample}; shift 2 ;;
        --whitelist) WHITELIST=${2:?missing value for --whitelist}; shift 2 ;;
        --mmi) MMI=${2:?missing value for --mmi}; shift 2 ;;
        --gtf) GTF=${2:?missing value for --gtf}; shift 2 ;;
        --junctions) JUNCTIONS=${2:?missing value for --junctions}; shift 2 ;;
        --out-dir) OUT_DIR=${2:?missing value for --out-dir}; shift 2 ;;
        --genes-mtx) GENES_MTX=${2:?missing value for --genes-mtx}; shift 2 ;;
        --genome) GENOME=${2:?missing value for --genome}; shift 2 ;;
        --extract-threads) EXTRACT_THREADS=${2:?missing value for --extract-threads}; shift 2 ;;
        --map-threads) MAP_THREADS=${2:?missing value for --map-threads}; shift 2 ;;
        --annotate-threads) ANNOTATE_THREADS=${2:?missing value for --annotate-threads}; shift 2 ;;
        --sort-threads) SORT_THREADS=${2:?missing value for --sort-threads}; shift 2 ;;
        --molecule-threads) MOLECULE_THREADS=${2:?missing value for --molecule-threads}; shift 2 ;;
        --split) SPLIT=${2:?missing value for --split}; shift 2 ;;
        --sort-mem) SORT_MEM=${2:?missing value for --sort-mem}; shift 2 ;;
        --tmp-dir) TMP_DIR=${2:?missing value for --tmp-dir}; shift 2 ;;
        --expect-cells) EXPECT_CELLS=${2:?missing value for --expect-cells}; shift 2 ;;
        --no-three-prime) THREE_PRIME=0; shift ;;
        --multimap-tags) MULTIMAP_TAGS=1; shift ;;
        --audit) AUDIT=1; shift ;;
        --barcode-metrics) BARCODE_METRICS=(--barcode-metrics); shift ;;
        --audit-label) AUDIT_LABEL=${2:?missing value for --audit-label}; shift 2 ;;
        --molecule-sample) MOLECULE_SAMPLE=${2:?missing value for --molecule-sample}; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) die "unknown option: $1 (use --help)" ;;
    esac
done

(( ${#FASTQS[@]} > 0 )) || die "supply at least one --fastqs directory"
for value in SAMPLE_SHEET SAMPLE WHITELIST MMI GTF JUNCTIONS OUT_DIR; do
    [[ -n ${!value} ]] || die "missing required --${value,,} input"
done
for command in sbxpress samtools; do require_command "$command"; done
for path in "${FASTQS[@]}" "$SAMPLE_SHEET" "$WHITELIST" "$MMI" "$GTF" "$JUNCTIONS"; do
    require_path "$path"
done

EXTRACT_DIR="$OUT_DIR/extract"
SAMPLE_OUT="$OUT_DIR/$SAMPLE"
ANNOTATION_DIR="$SAMPLE_OUT/annotation"
MOLECULES_DIR="$SAMPLE_OUT/molecules"
MATRIX_DIR="$SAMPLE_OUT/outs"
AUDIT_DIR="$SAMPLE_OUT/audit"
LOG_DIR="$SAMPLE_OUT/logs"
GENES_MTX=${GENES_MTX:-"$OUT_DIR/reference/genes.mtx"}
TMP_DIR=${TMP_DIR:-"$OUT_DIR/tmp"}
ANNOTATED_BAM="$ANNOTATION_DIR/$SAMPLE.annotated.coordinate.bam"
MOLECULE_BAM="$MOLECULES_DIR/$SAMPLE.molecules.coordinate.bam"
RAW_H5="$MATRIX_DIR/raw_feature_bc_matrix.h5"
FILTERED_H5="$MATRIX_DIR/filtered_feature_bc_matrix.h5"
FILTERED_BARCODES="$MATRIX_DIR/filtered_feature_bc_matrix/barcodes.tsv.gz"
AUDIT_JSON="$AUDIT_DIR/$SAMPLE.audit.json"
AUDIT_MOLECULES="$AUDIT_DIR/$SAMPLE.molecules.tsv.gz"
AUDIT_LABEL=${AUDIT_LABEL:-"SBX_$SAMPLE"}

source "$SCRIPT_DIR/pipeline_state.sh"
# Mapping preset. minimap2 2.28 accepts `-x splice:sr` but runs the long-read
# `splice` preset (2.29 added a real short-read one), so every validated run
# used `splice` scoring; naming it gives the same mapping on any version
# (2.30's `-x splice` matches 2.28's `-x splice:sr`). -m 20 -s 30 lower the
# minimum chain and alignment scores so short reads are not dropped: on
# CAT02 M2_GEX at matched depth, +16% countable Illumina reads and +5% SBX
# (CHANGELOG 0.6.2). Adding -f 0.00002 gained little for 4.3x the CPU.
MAP_PRESET=(-x splice -m 20 -s 30)
bind_pipeline_state "$EXTRACT_DIR" sbx gex "$WHITELIST" "$SAMPLE_SHEET" "$SPLIT" "${FASTQS[@]}"
bind_pipeline_state "$SAMPLE_OUT" sbx gex "$SAMPLE" "$MMI" "$GTF" "$JUNCTIONS" "$GENOME" \
    "$THREE_PRIME" "$MULTIMAP_TAGS" "$EXPECT_CELLS" "genes-manifest=$GENES_MTX"

mkdir -p "$EXTRACT_DIR" "$ANNOTATION_DIR" "$MOLECULES_DIR" "$MATRIX_DIR" "$LOG_DIR" "$TMP_DIR" "$(dirname "$GENES_MTX")"
(( AUDIT == 0 )) || mkdir -p "$AUDIT_DIR"

echo "sbxpress version: $(sbxpress --version)"
echo "sample:           $SAMPLE"
echo "output:           $SAMPLE_OUT"
echo "map/annotate/sort threads: $MAP_THREADS/$ANNOTATE_THREADS/$SORT_THREADS"
echo "multimap tags:    $([[ $MULTIMAP_TAGS == 1 ]] && echo enabled || echo disabled)"

if [[ ! -s "$EXTRACT_DIR/metrics_summary.json" ]]; then
    echo "[1/5] extracting and demultiplexing raw SBX data"
    sbxpress extract \
        --fastqs "${FASTQS[@]}" \
        --bclconvert-samplesheet "$SAMPLE_SHEET" \
        --whitelist "$WHITELIST" \
        --out-dir "$EXTRACT_DIR" \
        --threads "$EXTRACT_THREADS" \
        --split "$SPLIT" \
        --lenient \
        --minimap2 \
        ${BARCODE_METRICS[@]+"${BARCODE_METRICS[@]}"} \
        2> "$LOG_DIR/extract.log"
else
    echo "[1/5] extraction already complete; reusing $EXTRACT_DIR"
fi

shopt -s nullglob
CDNA_FASTQS=("$EXTRACT_DIR/${SAMPLE}_S"*_R2_*.fastq.gz)
(( ${#CDNA_FASTQS[@]} > 0 )) || die "no extracted R2 FASTQs found for $SAMPLE"

if [[ ! -s "$ANNOTATED_BAM.bai" ]]; then
    echo "[2/5] mapping, annotation, and one coordinate sort"
    ANNOTATE_ARGS=(--gtf "$GTF" --output - --metrics "$ANNOTATION_DIR/$SAMPLE.annotation.metrics.json" --genes-mtx "$GENES_MTX" --threads "$ANNOTATE_THREADS" ${BARCODE_METRICS[@]+"${BARCODE_METRICS[@]}"})
    if (( MULTIMAP_TAGS )); then
        ANNOTATE_ARGS+=(--multimap-tags)
        # Primary + ten secondary alignments gives NH values 1..11, matching
        # the audit's final top-coded NH=11+ ambiguity bucket.
        MAP_MULTIMAP_ARGS=(-N 10)
    else
        ANNOTATE_ARGS+=(--unordered)
        MAP_MULTIMAP_ARGS=(--secondary=no)
    fi
    (( THREE_PRIME == 0 )) || ANNOTATE_ARGS+=(--three-prime)
    sbxpress minimap \
        -a "${MAP_PRESET[@]}" --junc-bed "$JUNCTIONS" "${MAP_MULTIMAP_ARGS[@]}" -y \
        -t "$MAP_THREADS" "$MMI" "${CDNA_FASTQS[@]}" \
        2> "$LOG_DIR/minimap2.log" \
    | sbxpress annotate - "${ANNOTATE_ARGS[@]}" \
        2> "$LOG_DIR/annotate.log" \
    | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" \
        -T "$TMP_DIR/$SAMPLE.annotated" -o "$ANNOTATED_BAM" -
    samtools index -@ "$SORT_THREADS" "$ANNOTATED_BAM"
    samtools quickcheck -v "$ANNOTATED_BAM"
else
    echo "[2/5] annotated BAM already complete; reusing $ANNOTATED_BAM"
fi

if [[ ! -s "$MOLECULE_BAM.bai" ]]; then
    echo "[3/5] correcting UMIs and resolving gene collisions"
    sbxpress molecules "$ANNOTATED_BAM" \
        --output "$MOLECULE_BAM" \
        --metrics "$MOLECULES_DIR/$SAMPLE.molecules.metrics.json" \
        --threads "$MOLECULE_THREADS" \
        ${BARCODE_METRICS[@]+"${BARCODE_METRICS[@]}"} \
        2> "$LOG_DIR/molecules.log"
    samtools index -@ "$SORT_THREADS" "$MOLECULE_BAM"
    samtools quickcheck -v "$MOLECULE_BAM"
else
    echo "[3/5] molecule BAM already complete; reusing $MOLECULE_BAM"
fi

if [[ ! -s "$RAW_H5" || ! -e "$MATRIX_DIR/.raw.done" ]]; then
    RAW_WORK="$MATRIX_DIR/.raw.partial"
    prepare_matrix_stage "$RAW_WORK"
    echo "[4/5] creating complete raw feature-barcode matrix"
    sbxpress matrix "$MOLECULE_BAM" \
        --output "$RAW_WORK" --genome "$GENOME" \
        --features "$GENES_MTX" --cell-barcodes "$WHITELIST" \
        2> "$LOG_DIR/matrix.log"
    publish_matrix_stage "$RAW_WORK" "$MATRIX_DIR" .raw.done
else
    echo "[4/5] raw matrix already complete; reusing $RAW_H5"
fi

if [[ ! -s "$FILTERED_H5" || ! -e "$MATRIX_DIR/.cells.done" ]]; then
    CELL_WORK="$MATRIX_DIR/.cells.partial"
    prepare_matrix_stage "$CELL_WORK"
    echo "[5/5] calling cells and writing filtered counts"
    CELL_ARGS=()
    [[ -z "$EXPECT_CELLS" ]] || CELL_ARGS+=(--expect-cells "$EXPECT_CELLS")
    sbxpress cells "$RAW_H5" --output "$CELL_WORK" ${CELL_ARGS[@]+"${CELL_ARGS[@]}"} \
        2> "$LOG_DIR/cells.log"
    publish_matrix_stage "$CELL_WORK" "$MATRIX_DIR" .cells.done
else
    echo "[5/5] filtered matrix already complete; reusing $FILTERED_H5"
fi

if (( AUDIT == 1 )); then
    if [[ ! -s "$AUDIT_JSON" || ! -s "$AUDIT_MOLECULES" ]]; then
        echo "[audit] generating optional comparison audit"
        sbxpress audit "$MOLECULE_BAM" --label "$AUDIT_LABEL" \
            --output "$AUDIT_JSON" --molecules "$AUDIT_MOLECULES" \
            --molecule-sample "$MOLECULE_SAMPLE" --threads "$SORT_THREADS" \
            --cells "$FILTERED_BARCODES" \
            2> "$LOG_DIR/audit.log"
    else
        echo "[audit] audit already complete; reusing $AUDIT_JSON"
    fi
fi

echo "Pipeline complete"
echo "  raw counts:      $RAW_H5"
echo "  filtered counts: $FILTERED_H5"
echo "  filtered MEX:    $MATRIX_DIR/filtered_feature_bc_matrix"
echo "  molecule BAM:    $MOLECULE_BAM"
