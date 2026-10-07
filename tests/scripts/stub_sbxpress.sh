#!/usr/bin/env bash
# Test stub for `sbxpress` used by pipeline_smoke.sh. For every subcommand it
# (1) checks each --flag the pipeline passes is accepted by the REAL binary
# (REAL_SBXPRESS <sub> --help), and (2) creates the outputs the next stage
# expects. Calls are appended to $STUB_LOG.
set -euo pipefail
sub=${1:?subcommand}; shift
echo "sbxpress $sub $*" >> "$STUB_LOG"
if [[ $sub == --version ]]; then echo "sbxpress stub"; exit 0; fi
if [[ $sub != minimap ]]; then
    help=$("$REAL_SBXPRESS" "$sub" --help 2>&1) || { echo "stub: unknown subcommand $sub" >&2; exit 2; }
    for arg in "$@"; do
        [[ $arg == --* ]] || continue
        flag=${arg%%=*}
        grep -qE -- "(^|[ ,])$flag([ ,=<]|$)" <<< "$help" \
            || { echo "stub: 'sbxpress $sub' does not accept $flag" >&2; exit 2; }
    done
fi
# value OPT...: value following the first matching option
value() { local a=("${ARGS[@]}"); for ((i=0; i<${#a[@]}; i++)); do for o in "$@"; do [[ ${a[$i]} == "$o" ]] && { echo "${a[$((i+1))]}"; return; }; done; done; }
ARGS=("$@")
case "$sub" in
    extract)
        out=$(value --out-dir -o); mkdir -p "$out"; echo '{}' > "$out/metrics_summary.json"
        ext=fastq.gz; for a in "${ARGS[@]}"; do [[ $a == --no-gzip ]] && ext=fastq; done
        for chunk in 001 002; do echo stub > "$out/${STUB_SAMPLE}_S1_L001_R2_${chunk}.$ext"; done ;;
    illumina) echo stub > "$(value --output -o)" ;;
    minimap) echo "@HD	VN:1.6" ;;
    annotate) cat > /dev/null; echo "annotated" ;;
    probeset-align) echo stub > "$(value --output -o)" ;;
    molecules) echo stub > "$(value --output -o)" ;;
    matrix) out=$(value --output -o); mkdir -p "$out"; echo stub > "$out/raw_feature_bc_matrix.h5"
        [[ ${STUB_FAIL_MATRIX:-0} != 1 ]] || exit 9 ;;
    cells)
        list=$(value --tagseq-barcodes || true)
        [[ -z $list || -f $list ]] || { echo "stub: barcode-list argument was split: $list" >&2; exit 2; }
        out=$(value --output -o); mkdir -p "$out/filtered_feature_bc_matrix"
        echo stub > "$out/filtered_feature_bc_matrix.h5"
        printf 'genome,barcode\nG,AAAC-1\n' > "$out/filtered_barcodes.csv" ;;
    probe-demux)
        out=$(value --output -o); mkdir -p "$out"
        echo AAAC-1 > "$out/S1.barcodes.txt"; echo BC001 > "$out/S1.probe_barcodes.txt" ;;
    audit) echo stub > "$(value --output -o)"; m=$(value --molecules); [[ -z $m ]] || echo stub > "$m" ;;
    *) echo "stub: no behavior for $sub" >&2; exit 2 ;;
esac
