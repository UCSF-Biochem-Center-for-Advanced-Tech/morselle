#!/usr/bin/env bash
# Runs scripts/sbxpress-pipeline.sh (every input x analysis mode) and the
# production scripts/sbxpress-gex-pipeline.sh end to end against stub
# sbxpress/samtools commands. Fails if a script has a syntax error, passes a
# flag the real binary rejects, or does not reach "Pipeline complete".
# Usage: pipeline_smoke.sh /path/to/real/sbxpress
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export REAL_SBXPRESS=${1:?path to the real sbxpress binary}
WORK=$(mktemp -d); trap 'rm -rf "$WORK"' EXIT
mkdir -p "$WORK/bin" "$WORK/ref" "$WORK/fastqs"
ln -s "$REPO/tests/scripts/stub_sbxpress.sh" "$WORK/bin/sbxpress"
ln -s "$REPO/tests/scripts/stub_samtools.sh" "$WORK/bin/samtools"
export PATH="$WORK/bin:$PATH" STUB_SAMPLE=S1
for f in wl.txt.gz sheet.csv ref.mmi genes.gtf junc.bed probes.csv pbwl.txt config.csv cells.txt r1.fq.gz r2.fq.gz; do
    : > "$WORK/ref/$f"
done
R="$WORK/ref"
COMMON=(--sample S1 --whitelist "$R/wl.txt.gz")
SBX=(--fastqs "$WORK/fastqs" --bclconvert-samplesheet "$R/sheet.csv")
ILL=(--r1 "$R/r1.fq.gz" --r2 "$R/r2.fq.gz")
GEX=(--mmi "$R/ref.mmi" --gtf "$R/genes.gtf" --junctions "$R/junc.bed")
FLEX=(--probe-set "$R/probes.csv")

check() {
    local name=$1; shift
    export STUB_LOG="$WORK/$name.calls"
    : > "$STUB_LOG"
    if ! "$@" --out-dir "$WORK/$name" > "$WORK/$name.out" 2>&1; then
        echo "FAIL $name"; cat "$WORK/$name.out"; exit 1
    fi
    grep -q "Pipeline complete" "$WORK/$name.out" || { echo "FAIL $name: incomplete"; cat "$WORK/$name.out"; exit 1; }
    echo "ok   $name"
}
expect_call() { grep -qE -- "$2" "$WORK/$1.calls" || { echo "FAIL $1: expected call matching: $2"; cat "$WORK/$1.calls"; exit 1; }; }
reject_call() { ! grep -qE -- "$2" "$WORK/$1.calls" || { echo "FAIL $1: unexpected call matching: $2"; exit 1; }; }

P="$REPO/scripts/sbxpress-pipeline.sh"
check sbx_gex "$P" "${COMMON[@]}" "${SBX[@]}" "${GEX[@]}" --audit --expect-cells 500
expect_call sbx_gex "sbxpress extract .*--lenient --minimap2"
expect_call sbx_gex "sbxpress annotate - .*--unordered.*--three-prime"
expect_call sbx_gex "sbxpress cells .*--expect-cells 500"
reject_call sbx_gex "molecules .*--tagseq"
check sbx_gex_fastbam "$P" "${COMMON[@]}" "${SBX[@]}" "${GEX[@]}" --molecule-compression-level 1
expect_call sbx_gex_fastbam "molecules .*--compression-level 1"
check sbx_gex_multimap "$P" "${COMMON[@]}" "${SBX[@]}" "${GEX[@]}" --multimap-tags --no-three-prime
expect_call sbx_gex_multimap "minimap .*-N 10"
reject_call sbx_gex_multimap "--three-prime"
check sbx_gex_nogzip "$P" "${COMMON[@]}" "${SBX[@]}" "${GEX[@]}" --extract-no-gzip
expect_call sbx_gex_nogzip "extract .*--no-gzip"
expect_call sbx_gex_nogzip "minimap .*S1_S1_L001_R2_001.fastq .*S1_S1_L001_R2_002.fastq$"
check sbx_flex_nogzip "$P" --analysis-mode flex "${COMMON[@]}" "${SBX[@]}" "${FLEX[@]}" --probe-barcode-whitelist "$R/pbwl.txt" --extract-no-gzip
expect_call sbx_flex_nogzip "probeset-align .*R2_002.fastq .*--output .*/R2_002.bam.partial.bam|probeset-align .*S1_S1_L001_R2_002.fastq "
check sbx_tagseq "$P" --analysis-mode tagseq "${COMMON[@]}" "${SBX[@]}" "${GEX[@]}" --tagseq-barcodes "$R/cells.txt"
expect_call sbx_tagseq "sbxpress molecules .*--tagseq"
expect_call sbx_tagseq "sbxpress cells .*--tagseq-barcodes"
check sbx_tagseq_shared "$P" --analysis-mode tagseq "${COMMON[@]}" "${SBX[@]}" "${GEX[@]}" --no-fractional-collisions
expect_call sbx_tagseq_shared "molecules .*--tagseq --no-fractional-collisions"
check sbx_flex "$P" --analysis-mode flex "${COMMON[@]}" "${SBX[@]}" "${FLEX[@]}" --probe-barcode-whitelist "$R/pbwl.txt"
expect_call sbx_flex "extract .*--flex --probe-barcode-whitelist"
expect_call sbx_flex "probeset-align .*S1_S1_L001_R2_002.fastq.gz"
expect_call sbx_flex "molecules .*--flex"
check sbx_flex_multi "$P" --analysis-mode flex "${COMMON[@]}" "${SBX[@]}" "${FLEX[@]}" --probe-barcode-wl "$R/pbwl.txt" --config "$R/config.csv"
expect_call sbx_flex_multi "probe-demux .*--probe-sample-sheet"
expect_call sbx_flex_multi "sbxpress matrix .*--probe-barcodes"
check illumina_gex "$P" --input-mode illumina "${COMMON[@]}" "${ILL[@]}" "${GEX[@]}"
expect_call illumina_gex "sbxpress illumina --r1 .* --r2 "
check illumina_flex "$P" --input-mode illumina --analysis-mode flex "${COMMON[@]}" "${ILL[@]}" "${FLEX[@]}"
reject_call illumina_flex "molecules .*--flex"

# Restartability: a second run reuses every stage and calls no tool.
export STUB_LOG="$WORK/rerun.calls"; : > "$STUB_LOG"
"$P" "${COMMON[@]}" "${SBX[@]}" "${GEX[@]}" --out-dir "$WORK/sbx_gex" > /dev/null 2>&1
if grep -qv -- "--version" "$STUB_LOG"; then echo "FAIL rerun called tools:"; cat "$STUB_LOG"; exit 1; fi
echo "ok   rerun reuses finished stages"

# Bad input combinations fail fast with a message.
if "$P" --input-mode illumina --analysis-mode flex "${COMMON[@]}" "${ILL[@]}" "${FLEX[@]}" \
    --config "$R/config.csv" --out-dir "$WORK/bad" > "$WORK/bad.out" 2>&1; then
    echo "FAIL illumina multiplexed flex should be rejected"; exit 1
fi
echo "ok   rejects multiplexed Illumina Flex"

check gex_script "$REPO/scripts/sbxpress-gex-pipeline.sh" "${COMMON[@]}" "${SBX[@]}" "${GEX[@]}" --audit
echo "all pipeline smoke tests passed"
