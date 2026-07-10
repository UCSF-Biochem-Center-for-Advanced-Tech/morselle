#!/usr/bin/env bash
set -euo pipefail

repo=$(cd "$(dirname "$0")/.." && pwd)
fixtures=${1:-"$repo/../../test_fixtures"}
binary="$repo/target/release/cr-bam-to-mtx"

[[ -d "$fixtures" ]] || {
    echo "fixture directory not found: $fixtures" >&2
    echo "usage: $0 /path/to/test_fixtures" >&2
    exit 2
}

cargo build --release --locked --manifest-path "$repo/Cargo.toml" --bin cr-bam-to-mtx

for fixture in "$fixtures"/*_gene_matrix; do
    expected="$fixture/raw_feature_bc_matrix"
    actual=$(mktemp -d "${TMPDIR:-/tmp}/cr-bam-to-mtx.XXXXXX")
    trap 'rm -rf "$actual"' EXIT
    "$binary" "$fixture/possorted_genome_bam.bam" \
        --features "$expected/features.tsv.gz" \
        --barcodes "$expected/barcodes.tsv.gz" \
        --output "$actual"
    diff -u <(gzip -dc "$expected/features.tsv.gz") <(gzip -dc "$actual/features.tsv.gz")
    diff -u <(gzip -dc "$expected/barcodes.tsv.gz") <(gzip -dc "$actual/barcodes.tsv.gz")
    python3 "$repo/tests/validate_cr_matrix.py" \
        "$expected/matrix.mtx.gz" "$actual/matrix.mtx.gz"
    test -s "$actual/metrics_summary.json"
    rm -rf "$actual"
    trap - EXIT
    echo "PASS $(basename "$fixture")"
done
