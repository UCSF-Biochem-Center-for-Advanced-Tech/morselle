#!/usr/bin/env bash
# Test stub for samtools: creates sort/index outputs, logs calls.
set -euo pipefail
sub=$1; shift
echo "samtools $sub $*" >> "$STUB_LOG"
case "$sub" in
    sort) cat > /dev/null; while (( $# )); do [[ $1 == -o ]] && { echo stub > "$2"; break; }; shift; done ;;
    index) for a in "$@"; do [[ $a == *.bam ]] && echo stub > "$a.bai"; done ;;
    cat) echo catted ;;
    quickcheck) ;;
    *) echo "stub samtools: unsupported $sub" >&2; exit 2 ;;
esac
