#!/usr/bin/env bash
# Unified, restartable sbxpress pipeline:
#   input:    raw Roche SBX device output (--input-mode sbx, default)
#             or paired Illumina R1/R2 FASTQs (--input-mode illumina)
#   analysis: 3' GEX (default), pooled CAT tagseq, or 10x Flex probe sets.
#
# The GEX path is stage-for-stage the production sbxpress-gex-pipeline.sh
# (same commands, flags, file names, and restart markers). Tagseq differs
# only in `molecules --tagseq` and optional whitelist cell calling. Flex
# follows scripts/run_full_cat06.sh: extract --flex, probeset-align per
# extracted chunk, molecules --flex, then either single-sample counting or
# probe-demux + scripts/run_flex_multi.sh for multiplexed samples.
#
# Every stage checks for its finished output first, so a killed run can be
# restarted with the same command line. Intermediate BAMs are kept for that
# reason. Set SBXPRESS_PIPELINE_DRY_RUN=1 to print commands without running.
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  sbxpress-pipeline.sh [--input-mode sbx|illumina] [--analysis-mode gex|tagseq|flex] \
    --sample NAME --whitelist 10X_BARCODES.txt.gz --out-dir OUTPUT_DIRECTORY \
    <input options> <analysis options> [options]

Input (SBX, default):
  --fastqs DIR [...]       Raw SBX FASTQ directories, one or more device outputs.
  --bclconvert-samplesheet FILE
                           BCL Convert-compatible sample sheet.
Input (Illumina):
  --r1 FASTQ [...]         Read 1 FASTQ(s); one per lane, same order as --r2.
  --r2 FASTQ [...]         Read 2 FASTQ(s).

Common:
  --sample NAME            Sample (SBX: sample-sheet Sample_ID) to map and count.
  --whitelist FILE         10x cell-barcode whitelist.
  --out-dir DIR            Pipeline output directory.

GEX / tagseq:
  --mmi FILE               minimap2 splice index.
  --gtf FILE               GENCODE-compatible annotation GTF.
  --junctions FILE         Junction BED used when the minimap index was built.
  --genes-mtx FILE         Canonical feature manifest; default OUT/reference/genes.mtx
                           (Flex: OUT/reference/probe_set_genes.tsv).
  --no-three-prime         Disable terminal-exon preference tags (default: enabled).
  --multimap-tags          Retain up to ten secondary hits and add NH/HI/mm tags
                           (query-grouped streaming annotation, slower).
Tagseq:
  --tagseq-barcodes FILE   Call exactly these barcodes as cells (skips EmptyDrops).
  --fractional-collisions  Split each cross-gene UMI collision by read support
                           (counts one molecule per collision; real-valued matrix).
                           Default (tagseq): count every gene in a collision in full.
  --position-aware-umis    Deduplicate on (cell, gene, UMI, read 5' start): chance
                           UMI matches stay separate; UMI errors corrected per start.
Flex:
  --probe-set FILE         10x Flex probe-set reference CSV.
  --probe-barcode-whitelist FILE (alias --probe-barcode-wl)
                           Probe-barcode whitelist (e.g. flex-v2-384.txt); SBX only.
  --probe-sample-sheet FILE (alias --config)
                           Multiplexed Flex: config.csv [samples] section. Runs
                           probe-demux + run_flex_multi.sh for per-sample outputs.

Options:
  --genome NAME            Feature genome label; default GRCh38.
  --threads N              Default for every --*-threads not given; default 8.
  --extract-threads N      Extraction / Illumina conversion workers; default 32.
  --map-threads N          minimap2 workers; default 32.
  --annotate-threads N     Annotation workers (unordered mode); default 8.
  --align-threads N        probeset-align workers (Flex); default 16.
  --sort-threads N         samtools sort/index workers; default 8.
  --molecule-threads N     Molecule-calling BAM threads; default 8.
  --split N                Extracted FASTQ shards per sample (SBX); default 16.
  --extract-no-gzip        Write extracted FASTQs uncompressed (SBX). Measured 7.3x
                           faster extraction and 2.2x faster extract+mapping
                           (PROFILING.md) at ~4x the intermediate disk space.
  --molecule-compression-level N
                           BGZF level 0-9 for the molecule BAM (default: htslib's 6;
                           1 was ~24% faster overall for a ~12% larger file).
  --sort-mem SIZE          samtools memory per sort worker; default 1G.
  --tmp-dir DIR            Sort temporary directory; default OUT/tmp.
  --expect-cells N         Expected recovered cells; passed to sbxpress cells.
  --audit                  Also emit comparison audit JSON and molecule sample.
  --audit-label TEXT       Audit label; default SBX_SAMPLE / ILLUMINA_SAMPLE.
  --molecule-sample N      Audit 1-in-N molecule sample; default 64.
  --help                   Show this help.

Outputs: OUT/SAMPLE/outs/{raw,filtered}_feature_bc_matrix.h5 and MEX directories
(multiplexed Flex: OUT/SAMPLE/flex_samples/<sample_id>/gex_outs/...).
EOF
}

die() { echo "ERROR: $*" >&2; exit 1; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"; }
require_path() { [[ -e "$1" ]] || die "required input not found: $1"; }
# run LOG cmd...: execute (or print, in dry-run) with stderr to LOG.
run() {
    local log=$1; shift
    if [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} == 1 ]]; then
        printf '+'; printf ' %q' "$@"; printf '  2> %s\n' "$log"
    else
        "$@" 2> "$log"
    fi
}
# Collect values following a list option until the next --option.
take_list() { LIST=(); while (( $# > 0 )) && [[ $1 != --* ]]; do LIST+=("$1"); shift; done; }

INPUT_MODE=sbx
ANALYSIS_MODE=gex
FASTQS=(); R1_FASTQS=(); R2_FASTQS=()
SAMPLE_SHEET='' SAMPLE='' WHITELIST='' OUT_DIR=''
MMI='' GTF='' JUNCTIONS='' GENES_MTX='' GENOME=GRCh38
PROBE_SET='' PROBE_BARCODE_WL='' PROBE_SAMPLE_SHEET='' TAGSEQ_BARCODES=''
THREADS='' EXTRACT_THREADS='' MAP_THREADS='' ANNOTATE_THREADS='' ALIGN_THREADS=''
SORT_THREADS='' MOLECULE_THREADS=''
SPLIT=16 SORT_MEM=1G TMP_DIR='' EXPECT_CELLS=''
THREE_PRIME=1 MULTIMAP_TAGS=0 FRACTIONAL=0 POSITION_AWARE=0 EXTRACT_NO_GZIP=0
AUDIT=0 AUDIT_LABEL='' MOLECULE_SAMPLE=64 MOLECULE_COMPRESSION=''

while (( $# > 0 )); do
    case "$1" in
        --input-mode) INPUT_MODE=${2:?missing value for --input-mode}; shift 2 ;;
        --analysis-mode) ANALYSIS_MODE=${2:?missing value for --analysis-mode}; shift 2 ;;
        --fastqs|--r1|--r2)
            opt=$1; shift
            (( $# > 0 )) || die "missing value for $opt"
            take_list "$@"; shift ${#LIST[@]}
            case "$opt" in
                --fastqs) FASTQS+=("${LIST[@]}") ;;
                --r1) R1_FASTQS+=("${LIST[@]}") ;;
                --r2) R2_FASTQS+=("${LIST[@]}") ;;
            esac ;;
        --bclconvert-samplesheet) SAMPLE_SHEET=${2:?missing value for $1}; shift 2 ;;
        --sample) SAMPLE=${2:?missing value for $1}; shift 2 ;;
        --whitelist) WHITELIST=${2:?missing value for $1}; shift 2 ;;
        --out-dir) OUT_DIR=${2:?missing value for $1}; shift 2 ;;
        --mmi) MMI=${2:?missing value for $1}; shift 2 ;;
        --gtf) GTF=${2:?missing value for $1}; shift 2 ;;
        --junctions) JUNCTIONS=${2:?missing value for $1}; shift 2 ;;
        --genes-mtx) GENES_MTX=${2:?missing value for $1}; shift 2 ;;
        --genome) GENOME=${2:?missing value for $1}; shift 2 ;;
        --three-prime) THREE_PRIME=1; shift ;;
        --no-three-prime) THREE_PRIME=0; shift ;;
        --multimap-tags) MULTIMAP_TAGS=1; shift ;;
        --probe-set) PROBE_SET=${2:?missing value for $1}; shift 2 ;;
        --probe-barcode-whitelist|--probe-barcode-wl) PROBE_BARCODE_WL=${2:?missing value for $1}; shift 2 ;;
        --probe-sample-sheet|--config) PROBE_SAMPLE_SHEET=${2:?missing value for $1}; shift 2 ;;
        --tagseq-barcodes) TAGSEQ_BARCODES=${2:?missing value for $1}; shift 2 ;;
        --fractional-collisions) FRACTIONAL=1; shift ;;
        # Shared is the default again (fractional is opt-in); accepted for old command lines.
        --no-fractional-collisions) FRACTIONAL=0; shift ;;
        --position-aware-umis) POSITION_AWARE=1; shift ;;
        --threads) THREADS=${2:?missing value for $1}; shift 2 ;;
        --extract-threads) EXTRACT_THREADS=${2:?missing value for $1}; shift 2 ;;
        --map-threads) MAP_THREADS=${2:?missing value for $1}; shift 2 ;;
        --annotate-threads) ANNOTATE_THREADS=${2:?missing value for $1}; shift 2 ;;
        --align-threads) ALIGN_THREADS=${2:?missing value for $1}; shift 2 ;;
        --sort-threads) SORT_THREADS=${2:?missing value for $1}; shift 2 ;;
        --molecule-threads) MOLECULE_THREADS=${2:?missing value for $1}; shift 2 ;;
        --split) SPLIT=${2:?missing value for $1}; shift 2 ;;
        --extract-no-gzip) EXTRACT_NO_GZIP=1; shift ;;
        --molecule-compression-level) MOLECULE_COMPRESSION=${2:?missing value for $1}; shift 2 ;;
        --sort-mem) SORT_MEM=${2:?missing value for $1}; shift 2 ;;
        --tmp-dir) TMP_DIR=${2:?missing value for $1}; shift 2 ;;
        --expect-cells) EXPECT_CELLS=${2:?missing value for $1}; shift 2 ;;
        --audit) AUDIT=1; shift ;;
        --audit-label) AUDIT_LABEL=${2:?missing value for $1}; shift 2 ;;
        --molecule-sample) MOLECULE_SAMPLE=${2:?missing value for $1}; shift 2 ;;
        --help|-h) usage; exit 0 ;;
        *) die "unknown option: $1 (use --help)" ;;
    esac
done

# Thread defaults: explicit per-stage value > --threads > historical default.
EXTRACT_THREADS=${EXTRACT_THREADS:-${THREADS:-32}}
MAP_THREADS=${MAP_THREADS:-${THREADS:-32}}
ANNOTATE_THREADS=${ANNOTATE_THREADS:-${THREADS:-8}}
ALIGN_THREADS=${ALIGN_THREADS:-${THREADS:-16}}
SORT_THREADS=${SORT_THREADS:-${THREADS:-8}}
MOLECULE_THREADS=${MOLECULE_THREADS:-${THREADS:-8}}

[[ $INPUT_MODE =~ ^(sbx|illumina)$ ]] || die "invalid --input-mode '$INPUT_MODE' (sbx or illumina)"
[[ $ANALYSIS_MODE =~ ^(gex|tagseq|flex)$ ]] || die "invalid --analysis-mode '$ANALYSIS_MODE' (gex, tagseq or flex)"
for value in SAMPLE WHITELIST OUT_DIR; do
    [[ -n ${!value} ]] || die "missing required --$(tr '[:upper:]_' '[:lower:]-' <<< "$value") input"
done
INPUTS=("$WHITELIST")
if [[ $INPUT_MODE == sbx ]]; then
    (( ${#FASTQS[@]} > 0 )) || die "supply at least one --fastqs directory"
    [[ -n $SAMPLE_SHEET ]] || die "missing required --bclconvert-samplesheet input"
    INPUTS+=("${FASTQS[@]}" "$SAMPLE_SHEET")
else
    (( ${#R1_FASTQS[@]} > 0 )) || die "supply at least one --r1 FASTQ"
    (( ${#R1_FASTQS[@]} == ${#R2_FASTQS[@]} )) || die "--r1 and --r2 need the same number of files"
    INPUTS+=("${R1_FASTQS[@]}" "${R2_FASTQS[@]}")
fi
if [[ $ANALYSIS_MODE == flex ]]; then
    [[ -n $PROBE_SET ]] || die "--probe-set is required for --analysis-mode flex"
    INPUTS+=("$PROBE_SET")
    if [[ $INPUT_MODE == illumina ]]; then
        # sbxpress illumina does not locate probe barcodes (no PB tag), so
        # multiplexed samples cannot be separated.
        [[ -z $PROBE_SAMPLE_SHEET && -z $PROBE_BARCODE_WL ]] \
            || die "Illumina Flex input is single-sample only: probe-barcode demultiplexing needs --input-mode sbx"
    else
        [[ -n $PROBE_BARCODE_WL ]] || die "--probe-barcode-whitelist is required for SBX Flex extraction"
        INPUTS+=("$PROBE_BARCODE_WL")
    fi
    [[ -z $PROBE_SAMPLE_SHEET ]] || INPUTS+=("$PROBE_SAMPLE_SHEET")
else
    for value in MMI GTF JUNCTIONS; do
        [[ -n ${!value} ]] || die "missing required --$(tr '[:upper:]' '[:lower:]' <<< "$value") input for $ANALYSIS_MODE"
    done
    INPUTS+=("$MMI" "$GTF" "$JUNCTIONS")
fi
[[ -z $TAGSEQ_BARCODES ]] || INPUTS+=("$TAGSEQ_BARCODES")
if [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} != 1 ]]; then
    for command in sbxpress samtools; do require_command "$command"; done
    for path in "${INPUTS[@]}"; do require_path "$path"; done
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXTRACT_DIR="$OUT_DIR/extract"
SAMPLE_OUT="$OUT_DIR/$SAMPLE"
ANNOTATION_DIR="$SAMPLE_OUT/annotation"
MOLECULES_DIR="$SAMPLE_OUT/molecules"
MATRIX_DIR="$SAMPLE_OUT/outs"
AUDIT_DIR="$SAMPLE_OUT/audit"
LOG_DIR="$SAMPLE_OUT/logs"
# Flex gets its own manifest name so a GEX run in the same OUT never
# collides with it (probeset-align verifies, never overwrites, a manifest).
if [[ $ANALYSIS_MODE == flex ]]; then
    GENES_MTX=${GENES_MTX:-"$OUT_DIR/reference/probe_set_genes.tsv"}
else
    GENES_MTX=${GENES_MTX:-"$OUT_DIR/reference/genes.mtx"}
fi
TMP_DIR=${TMP_DIR:-"$OUT_DIR/tmp"}
ANNOTATED_BAM="$ANNOTATION_DIR/$SAMPLE.annotated.coordinate.bam"
MOLECULE_BAM="$MOLECULES_DIR/$SAMPLE.molecules.coordinate.bam"
RAW_H5="$MATRIX_DIR/raw_feature_bc_matrix.h5"
FILTERED_H5="$MATRIX_DIR/filtered_feature_bc_matrix.h5"
FILTERED_BARCODES="$MATRIX_DIR/filtered_feature_bc_matrix/barcodes.tsv.gz"
AUDIT_JSON="$AUDIT_DIR/$SAMPLE.audit.json"
AUDIT_MOLECULES="$AUDIT_DIR/$SAMPLE.molecules.tsv.gz"
[[ $INPUT_MODE == sbx ]] && AUDIT_LABEL=${AUDIT_LABEL:-"SBX_$SAMPLE"} || AUDIT_LABEL=${AUDIT_LABEL:-"ILLUMINA_$SAMPLE"}
MULTIPLEXED_FLEX=0
[[ $ANALYSIS_MODE == flex && -n $PROBE_SAMPLE_SHEET ]] && MULTIPLEXED_FLEX=1

source "$SCRIPT_DIR/pipeline_state.sh"
# Extraction is shared across samples; downstream analysis is sample-specific.
bind_pipeline_state "$EXTRACT_DIR" "$INPUT_MODE" "$ANALYSIS_MODE" "$WHITELIST" "$SAMPLE_SHEET" \
    "$PROBE_BARCODE_WL" "$SPLIT" "$EXTRACT_NO_GZIP" \
    ${FASTQS[@]+"${FASTQS[@]}"} ${R1_FASTQS[@]+"${R1_FASTQS[@]}"} ${R2_FASTQS[@]+"${R2_FASTQS[@]}"}
bind_pipeline_state "$SAMPLE_OUT" "$INPUT_MODE" "$ANALYSIS_MODE" "$SAMPLE" \
    "$MMI" "$GTF" "$JUNCTIONS" "$GENOME" "$PROBE_SET" "$PROBE_SAMPLE_SHEET" \
    "$THREE_PRIME" "$MULTIMAP_TAGS" "$FRACTIONAL" "$POSITION_AWARE" \
    "$EXPECT_CELLS" "$TAGSEQ_BARCODES" "$MOLECULE_COMPRESSION" "genes-manifest=$GENES_MTX"

mkdir -p "$EXTRACT_DIR" "$ANNOTATION_DIR" "$MOLECULES_DIR" "$MATRIX_DIR" "$LOG_DIR" "$TMP_DIR" "$(dirname "$GENES_MTX")"
(( AUDIT == 0 )) || mkdir -p "$AUDIT_DIR"

if [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} != 1 ]]; then
    echo "sbxpress version: $(sbxpress --version)"
fi
echo "sample:           $SAMPLE ($INPUT_MODE input, $ANALYSIS_MODE analysis)"
echo "output:           $SAMPLE_OUT"

# ---- [1/5] extraction / conversion to tagged cDNA FASTQs -------------------
if [[ $INPUT_MODE == sbx ]]; then
    if [[ ! -s "$EXTRACT_DIR/metrics_summary.json" ]]; then
        echo "[1/5] extracting and demultiplexing raw SBX data"
        EXTRACT_ARGS=(--fastqs "${FASTQS[@]}" --bclconvert-samplesheet "$SAMPLE_SHEET"
            --whitelist "$WHITELIST" --out-dir "$EXTRACT_DIR" --threads "$EXTRACT_THREADS"
            --split "$SPLIT" --lenient --minimap2)
        [[ $ANALYSIS_MODE != flex ]] || EXTRACT_ARGS+=(--flex --probe-barcode-whitelist "$PROBE_BARCODE_WL")
        (( EXTRACT_NO_GZIP == 0 )) || EXTRACT_ARGS+=(--no-gzip)
        run "$LOG_DIR/extract.log" sbxpress extract "${EXTRACT_ARGS[@]}"
    else
        echo "[1/5] extraction already complete; reusing $EXTRACT_DIR"
    fi
    shopt -s nullglob
    # .fastq.gz normally; plain .fastq with --extract-no-gzip.
    CDNA_FASTQS=("$EXTRACT_DIR/${SAMPLE}_S"*_R2_*.fastq.gz "$EXTRACT_DIR/${SAMPLE}_S"*_R2_*.fastq)
    shopt -u nullglob
else
    TAGGED_FASTQ="$EXTRACT_DIR/$SAMPLE.tagged.R2.fastq.gz"
    if [[ ! -s "$TAGGED_FASTQ" ]]; then
        echo "[1/5] converting Illumina R1/R2 to tagged cDNA FASTQ"
        ILLUMINA_ARGS=()
        for i in "${!R1_FASTQS[@]}"; do ILLUMINA_ARGS+=(--r1 "${R1_FASTQS[$i]}" --r2 "${R2_FASTQS[$i]}"); done
        # Write to a temporary name so an interrupted run is never mistaken
        # for a finished one by the -s check above.
        run "$LOG_DIR/illumina.log" sbxpress illumina "${ILLUMINA_ARGS[@]}" \
            --whitelist "$WHITELIST" --output "$TAGGED_FASTQ.partial.fastq.gz" \
            --metrics "$EXTRACT_DIR/$SAMPLE.illumina.metrics.json" --threads "$EXTRACT_THREADS"
        [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} == 1 ]] || mv "$TAGGED_FASTQ.partial.fastq.gz" "$TAGGED_FASTQ"
    else
        echo "[1/5] Illumina conversion already complete; reusing $TAGGED_FASTQ"
    fi
    CDNA_FASTQS=("$TAGGED_FASTQ")
fi
if [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} != 1 ]]; then
    (( ${#CDNA_FASTQS[@]} > 0 )) || die "no tagged cDNA FASTQs found for $SAMPLE in $EXTRACT_DIR"
fi

# ---- [2/5] alignment + annotation, one coordinate sort ---------------------
if [[ ! -s "$ANNOTATED_BAM.bai" ]]; then
    if [[ $ANALYSIS_MODE == flex ]]; then
        echo "[2/5] aligning to the Flex probe set"
        ALIGN_DIR="$ANNOTATION_DIR/probeset_chunks"
        mkdir -p "$ALIGN_DIR"
        CHUNK_BAMS=()
        for fq in "${CDNA_FASTQS[@]}"; do
            chunk_name=$(basename "$fq"); chunk_name=${chunk_name%.gz}; chunk_name=${chunk_name%.fastq}
            chunk="$ALIGN_DIR/$chunk_name.bam"
            CHUNK_BAMS+=("$chunk")
            [[ -s "$chunk" ]] && continue
            # probeset-align writes the BAM directly; build under a partial
            # name so a killed chunk is redone rather than reused.
            run "$LOG_DIR/probeset_align.$chunk_name.log" \
                sbxpress probeset-align "$fq" --probe-set "$PROBE_SET" \
                --genes-out "$GENES_MTX" --threads "$ALIGN_THREADS" \
                --metrics "$chunk.metrics.json" --output "$chunk.partial.bam"
            [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} == 1 ]] || mv "$chunk.partial.bam" "$chunk"
        done
        if [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} == 1 ]]; then
            echo "+ samtools cat ${CHUNK_BAMS[*]} | samtools sort ... -o $ANNOTATED_BAM"
        else
            samtools cat -o - "${CHUNK_BAMS[@]}" \
            | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" \
                -T "$TMP_DIR/$SAMPLE.annotated" -o "$ANNOTATED_BAM" -
        fi
    else
        echo "[2/5] mapping, annotation, and one coordinate sort"
        ANNOTATE_ARGS=(--gtf "$GTF" --output - --metrics "$ANNOTATION_DIR/$SAMPLE.annotation.metrics.json"
            --genes-mtx "$GENES_MTX" --threads "$ANNOTATE_THREADS")
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
        if [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} == 1 ]]; then
            echo "+ sbxpress minimap -a -x splice:sr --junc-bed $JUNCTIONS ${MAP_MULTIMAP_ARGS[*]} -y -t $MAP_THREADS $MMI ${CDNA_FASTQS[*]} | sbxpress annotate - ${ANNOTATE_ARGS[*]} | samtools sort ... -o $ANNOTATED_BAM"
        else
            sbxpress minimap \
                -a -x splice:sr --junc-bed "$JUNCTIONS" "${MAP_MULTIMAP_ARGS[@]}" -y \
                -t "$MAP_THREADS" "$MMI" "${CDNA_FASTQS[@]}" \
                2> "$LOG_DIR/minimap2.log" \
            | sbxpress annotate - "${ANNOTATE_ARGS[@]}" \
                2> "$LOG_DIR/annotate.log" \
            | samtools sort -@ "$SORT_THREADS" -m "$SORT_MEM" \
                -T "$TMP_DIR/$SAMPLE.annotated" -o "$ANNOTATED_BAM" -
        fi
    fi
    run "$LOG_DIR/index_annotated.log" samtools index -@ "$SORT_THREADS" "$ANNOTATED_BAM"
    run "$LOG_DIR/quickcheck_annotated.log" samtools quickcheck -v "$ANNOTATED_BAM"
else
    echo "[2/5] annotated BAM already complete; reusing $ANNOTATED_BAM"
fi

# ---- [3/5] UMI correction + gene collision resolution ----------------------
if [[ ! -s "$MOLECULE_BAM.bai" ]]; then
    echo "[3/5] correcting UMIs and resolving gene collisions"
    MOLECULE_ARGS=(--output "$MOLECULE_BAM" --metrics "$MOLECULES_DIR/$SAMPLE.molecules.metrics.json"
        --threads "$MOLECULE_THREADS")
    [[ -z $MOLECULE_COMPRESSION ]] || MOLECULE_ARGS+=(--compression-level "$MOLECULE_COMPRESSION")
    (( POSITION_AWARE == 0 )) || MOLECULE_ARGS+=(--position-aware-umis)
    case "$ANALYSIS_MODE" in
        tagseq)
            MOLECULE_ARGS+=(--tagseq)
            (( FRACTIONAL == 0 )) || MOLECULE_ARGS+=(--fractional-collisions) ;;
        flex)
            # Multiplexed Flex needs PB to keep samples sharing a GEM apart;
            # single-sample Illumina Flex has no PB tag at all.
            [[ $INPUT_MODE != sbx ]] || MOLECULE_ARGS+=(--flex) ;;
    esac
    run "$LOG_DIR/molecules.log" sbxpress molecules "$ANNOTATED_BAM" "${MOLECULE_ARGS[@]}"
    run "$LOG_DIR/index_molecules.log" samtools index -@ "$SORT_THREADS" "$MOLECULE_BAM"
    run "$LOG_DIR/quickcheck_molecules.log" samtools quickcheck -v "$MOLECULE_BAM"
else
    echo "[3/5] molecule BAM already complete; reusing $MOLECULE_BAM"
fi

CELL_ARGS=()
[[ -z $EXPECT_CELLS ]] || CELL_ARGS+=(--expect-cells "$EXPECT_CELLS")
[[ -z $TAGSEQ_BARCODES ]] || CELL_ARGS+=(--tagseq-barcodes "$TAGSEQ_BARCODES")

if (( MULTIPLEXED_FLEX )); then
    # ---- [4/5] multiplexed Flex: per-sample matrices + cells ----------------
    DEMUX_DIR="$SAMPLE_OUT/probe_demux"
    FLEX_OUT="$SAMPLE_OUT/flex_samples"
    if [[ ! -e "$DEMUX_DIR/.done" ]]; then
        echo "[4/5] assigning cell barcodes to Flex samples"
        run "$LOG_DIR/probe_demux.log" sbxpress probe-demux "$MOLECULE_BAM" \
            --probe-sample-sheet "$PROBE_SAMPLE_SHEET" --output "$DEMUX_DIR"
        [[ ${SBXPRESS_PIPELINE_DRY_RUN:-0} == 1 ]] || touch "$DEMUX_DIR/.done"
    fi
    echo "[5/5] per-sample matrices and cell calling (run_flex_multi.sh)"
    run "$LOG_DIR/flex_multi.log" "$SCRIPT_DIR/run_flex_multi.sh" \
        "$MOLECULE_BAM" "$DEMUX_DIR" "$FLEX_OUT" ${CELL_ARGS[@]+"${CELL_ARGS[@]}"}
else
    # ---- [4/5] raw matrix ---------------------------------------------------
    if [[ ! -s "$RAW_H5" || ! -e "$MATRIX_DIR/.raw.done" ]]; then
        RAW_WORK="$MATRIX_DIR/.raw.partial"
        prepare_matrix_stage "$RAW_WORK"
        echo "[4/5] creating complete raw feature-barcode matrix"
        MATRIX_ARGS=(--output "$RAW_WORK" --genome "$GENOME" --cell-barcodes "$WHITELIST"
            --features "$GENES_MTX")
        run "$LOG_DIR/matrix.log" sbxpress matrix "$MOLECULE_BAM" "${MATRIX_ARGS[@]}"
        publish_matrix_stage "$RAW_WORK" "$MATRIX_DIR" .raw.done
    else
        echo "[4/5] raw matrix already complete; reusing $RAW_H5"
    fi
    # ---- [5/5] cell calling -------------------------------------------------
    if [[ ! -s "$FILTERED_H5" || ! -e "$MATRIX_DIR/.cells.done" ]]; then
        CELL_WORK="$MATRIX_DIR/.cells.partial"
        prepare_matrix_stage "$CELL_WORK"
        echo "[5/5] calling cells and writing filtered counts"
        run "$LOG_DIR/cells.log" sbxpress cells "$RAW_H5" --output "$CELL_WORK" ${CELL_ARGS[@]+"${CELL_ARGS[@]}"}
        publish_matrix_stage "$CELL_WORK" "$MATRIX_DIR" .cells.done
    else
        echo "[5/5] filtered matrix already complete; reusing $FILTERED_H5"
    fi
fi

if (( AUDIT == 1 )); then
    if [[ ! -s "$AUDIT_JSON" || ! -s "$AUDIT_MOLECULES" ]]; then
        echo "[audit] generating optional comparison audit"
        AUDIT_ARGS=(--label "$AUDIT_LABEL" --output "$AUDIT_JSON" --molecules "$AUDIT_MOLECULES"
            --molecule-sample "$MOLECULE_SAMPLE" --threads "$SORT_THREADS")
        (( MULTIPLEXED_FLEX )) || AUDIT_ARGS+=(--cells "$FILTERED_BARCODES")
        run "$LOG_DIR/audit.log" sbxpress audit "$MOLECULE_BAM" "${AUDIT_ARGS[@]}"
    else
        echo "[audit] audit already complete; reusing $AUDIT_JSON"
    fi
fi

echo "Pipeline complete"
if (( MULTIPLEXED_FLEX )); then
    echo "  per-sample outputs: $FLEX_OUT/<sample_id>/gex_outs"
else
    echo "  raw counts:      $RAW_H5"
    echo "  filtered counts: $FILTERED_H5"
    echo "  filtered MEX:    $MATRIX_DIR/filtered_feature_bc_matrix"
fi
echo "  molecule BAM:    $MOLECULE_BAM"
