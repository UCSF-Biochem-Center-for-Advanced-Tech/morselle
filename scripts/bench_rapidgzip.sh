#!/usr/bin/env bash
# A/B the serial (flate2) and parallel (rapidgzip) gzip input decoders.
#
# Runs `morselle extract` twice over the same inputs — once with the default
# serial reader, once with MORSELLE_RAPIDGZIP=1 — and reports reads/sec for each.
# Both runs use the same binary and the same data, so the difference is the
# decoder and nothing else.
#
#   ./bench_rapidgzip.sh <fastqs> [outdir] [threads] [split] [extra morselle args...]
#
# `fastqs` is a directory of raw SBX FASTQs or a single file. The single-file case
# is the one worth measuring: that is where one serial reader feeds the whole
# worker pool and becomes the ceiling.
#
# Two things to get right before trusting the result:
#   * `split` must be high enough that output compression is not the bottleneck.
#     It defaults to threads/2 here; if the A and B numbers are both far below
#     what the cores should give, raise it and re-run.
#   * Page cache. The first run pays for reading the input from disk and the
#     second does not, which flatters B. This script reads the inputs once up
#     front to warm the cache, and runs A/B/A so drift is visible.
set -euo pipefail

if [[ $# -lt 1 ]]; then
  sed -n '1,21p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
fi

FASTQS=$1
OUT=${2:-rapidgzip_bench_out}
THREADS=${3:-0}
SPLIT=${4:-0}
shift $(( $# < 4 ? $# : 4 ))
EXTRA=("$@")

MORSELLE=${MORSELLE:-morselle}
command -v "$MORSELLE" >/dev/null || { echo "morselle not found (set MORSELLE=/path)" >&2; exit 1; }

if [[ "$THREADS" == "0" ]]; then
  THREADS=$( (nproc 2>/dev/null || sysctl -n hw.ncpu) )
fi
if [[ "$SPLIT" == "0" ]]; then
  SPLIT=$(( THREADS / 2 )); (( SPLIT < 1 )) && SPLIT=1
fi

mkdir -p "$OUT"

echo "inputs : $FASTQS"
echo "threads: $THREADS   split: $SPLIT"
echo

# Warm the page cache so run 1 is not penalised for cold reads.
echo "warming page cache..."
if [[ -d "$FASTQS" ]]; then
  find "$FASTQS" -maxdepth 1 \( -name '*.fastq.gz' -o -name '*.fq.gz' \) -exec cat {} + > /dev/null 2>&1 || true
else
  cat "$FASTQS" > /dev/null 2>&1 || true
fi
echo

# Run one configuration and print its throughput. $1 is a label, $2 is the value
# of MORSELLE_RAPIDGZIP for the run.
run_one() {
  local label=$1 enable=$2 dir="$OUT/$1"
  rm -rf "$dir"
  local log="$OUT/$label.log"
  MORSELLE_RAPIDGZIP="$enable" "$MORSELLE" extract \
    --fastqs "$FASTQS" --out-dir "$dir" \
    --threads "$THREADS" --split "$SPLIT" ${EXTRA[@]+"${EXTRA[@]}"} > "$log" 2>&1

  # Strip ANSI colour before scraping the structured log fields.
  local clean rps decoder
  clean=$(sed 's/\x1b\[[0-9;]*m//g' "$log")
  rps=$(grep -oE 'reads_per_sec=[0-9]+' <<<"$clean" | tail -1 | cut -d= -f2)
  decoder=$(grep -oE 'decoder=.*' <<<"$clean" | tail -1 | cut -d= -f2-)
  : "${rps:=0}"

  printf '%-14s %12s reads/s   %s\n' "$label" "$rps" "${decoder:-unknown}" >&2
  echo "$rps"
}

echo "=== A/B/A ==="
A1=$(run_one flate2_run1  0 | tail -1)
B1=$(run_one rapidgzip    1 | tail -1)
A2=$(run_one flate2_run2  0 | tail -1)
echo

if [[ "$A1" -gt 0 && "$A2" -gt 0 && "$B1" -gt 0 ]]; then
  python3 - "$A1" "$A2" "$B1" <<'PY'
import sys
a1, a2, b = (int(x) for x in sys.argv[1:4])
base = (a1 + a2) / 2
drift = abs(a1 - a2) / base * 100
print(f"flate2 mean : {base:,.0f} reads/s  (run-to-run drift {drift:.1f}%)")
print(f"rapidgzip   : {b:,.0f} reads/s")
print(f"speedup     : {b / base:.2f}x")
if drift > 5:
    print("\nWARNING: the two flate2 runs differ by more than 5%. The machine is")
    print("noisy or thermally throttled; treat the speedup as indicative only.")
if b / base < 1.05:
    print("\nNo meaningful gain. The usual causes, in order of likelihood:")
    print("  1. Input framing. pigz's default dependent-block framing does not")
    print("     parallelize at all. Test it: recompress a slice of the input with")
    print("     bgzip and re-run this script against that copy --")
    print("       zcat IN.fastq.gz | head -4000000 | bgzip -@ 8 > probe.fastq.gz")
    print("     A large speedup on the bgzip copy but none on the original means")
    print("     the framing is the problem, not the decoder.")
    print("  2. Not decode-bound. If files >= threads the per-file path already")
    print("     decodes many files at once; try a single large input.")
    print("  3. Output-bound. Raise --split and re-run.")
PY
fi

echo
echo "logs and outputs under $OUT/"
