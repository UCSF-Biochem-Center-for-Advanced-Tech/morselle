#!/usr/bin/env bash
# Run an annotated tag-seq BAM through each molecule-counting mode and
# compare every resulting matrix with a reference matrix (e.g. the same
# library sequenced on Illumina and counted independently).
#
# Modes (all with --tagseq):
#   shared        every gene in a UMI collision counted in full (default)
#   fractional    collisions split by read support (--fractional-collisions)
#   position      (cell, gene, UMI, read 5' start) keys (--position-aware-umis)
#   position_frac position-aware + fractional
#
# Usage:
#   scripts/compare_tagseq_modes.sh --bam SAMPLE.annotated.coordinate.bam \
#       --reference ILLUMINA/outs/raw_feature_bc_matrix.h5 --out-dir OUT \
#       [--features genes.mtx] [--barcodes samples.txt] [--threads 8] \
#       [--modes "shared fractional position"] [--molecule-args "--shards 8"]
#
# Each mode's molecules BAM and matrix are kept in OUT/<mode>/ and reused on
# a re-run. The comparison (scripts/compare_tagseq_modes.py, needs numpy;
# h5py for an .h5 reference) writes OUT/summary.tsv and OUT/per_gene.tsv.
# Set MORSELLE=/path/to/morselle to use a specific build.
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
BAM='' REFERENCE='' OUT_DIR='' FEATURES='' BARCODES='' THREADS=8
MODES="shared fractional position" MOLECULE_ARGS=''
while (( $# > 0 )); do
    case "$1" in
        --bam) BAM=${2:?}; shift 2 ;;
        --reference) REFERENCE=${2:?}; shift 2 ;;
        --out-dir) OUT_DIR=${2:?}; shift 2 ;;
        --features) FEATURES=${2:?}; shift 2 ;;
        --barcodes) BARCODES=${2:?}; shift 2 ;;
        --threads) THREADS=${2:?}; shift 2 ;;
        --modes) MODES=${2:?}; shift 2 ;;
        --molecule-args) MOLECULE_ARGS=${2:?}; shift 2 ;;
        --help|-h) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) die "unknown option: $1" ;;
    esac
done
[[ -n $BAM && -n $REFERENCE && -n $OUT_DIR ]] || die "--bam, --reference and --out-dir are required"
[[ -e $BAM && -e $REFERENCE ]] || die "input not found"
MORSELLE=${MORSELLE:-morselle}
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT_DIR"

COMPARE_ARGS=(--reference "$REFERENCE" --out-dir "$OUT_DIR")
[[ -z $BARCODES ]] || COMPARE_ARGS+=(--barcodes "$BARCODES")
for mode in $MODES; do
    case "$mode" in
        shared) flags=(--tagseq) ;;
        fractional) flags=(--tagseq --fractional-collisions) ;;
        position) flags=(--tagseq --position-aware-umis) ;;
        position_frac) flags=(--tagseq --position-aware-umis --fractional-collisions) ;;
        *) die "unknown mode: $mode" ;;
    esac
    dir="$OUT_DIR/$mode"
    mkdir -p "$dir"
    if [[ ! -s "$dir/outs/raw_feature_bc_matrix/matrix.mtx" ]]; then
        echo "[$mode] molecules ${flags[*]}"
        # shellcheck disable=SC2086 # MOLECULE_ARGS is a user-supplied flag list
        "$MORSELLE" molecules "$BAM" "${flags[@]}" $MOLECULE_ARGS \
            --output "$dir/molecules.bam" --threads "$THREADS" 2> "$dir/molecules.log"
        matrix_args=(--output "$dir/outs" --no-gzip)
        [[ -z $FEATURES ]] || matrix_args+=(--features "$FEATURES")
        echo "[$mode] matrix"
        "$MORSELLE" matrix "$dir/molecules.bam" "${matrix_args[@]}" 2> "$dir/matrix.log"
        rm -f "$dir/molecules.bam"  # large; the matrix and metrics are what's compared
    else
        echo "[$mode] already done; reusing $dir"
    fi
    COMPARE_ARGS+=(--mode "$mode=$dir/outs")
done

python3 "$SCRIPT_DIR/compare_tagseq_modes.py" "${COMPARE_ARGS[@]}"
echo "molecule metrics per mode: $OUT_DIR/<mode>/molecules.bam.metrics.json"
