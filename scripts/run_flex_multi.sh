#!/usr/bin/env bash
# FLEX_SUPPORT_PLAN.md Phase 3: per-sample GEX matrix + cell-calling.
#
# Zero new Rust here by design -- loops the existing, fully general
# `sbxpress matrix` and `sbxpress cells` commands over each sample's Phase 2
# `probe-demux` barcode list. Needs:
#   - a GEX-library BAM already through `sbxpress annotate` + `sbxpress molecules`
#     (carries CB/PB/xf, same input `probe-demux` itself consumes)
#   - `sbxpress probe-demux`'s output directory (one <sample_id>.barcodes.txt
#     per declared sample)
#
# `--skip-unlisted-barcodes` on the `matrix` call is load-bearing, not
# optional: the GEX BAM is the whole pooled multi-sample library, so every
# per-sample `matrix` invocation will legitimately see plenty of *other*
# samples' barcodes too. Without that flag `matrix` (correctly, for its
# original single-sample use case) hard-errors on the first one it meets.
#
# Usage:
#   scripts/run_flex_multi.sh <gex.bam> <probe_demux_out_dir> <output_dir> [extra `sbxpress cells` args...]
#
# Output, per declared sample:
#   <output_dir>/<sample_id>/gex_raw/{raw_feature_bc_matrix.h5,raw_feature_bc_matrix/}
#   <output_dir>/<sample_id>/gex_outs/{filtered_feature_bc_matrix.h5,filtered_barcodes.csv,...}

set -euo pipefail

usage() {
  echo "usage: $(basename "$0") <gex.bam> <probe_demux_out_dir> <output_dir> [extra 'sbxpress cells' args...]" >&2
  exit 1
}

[ $# -ge 3 ] || usage
gex_bam=$1
demux_dir=$2
out_dir=$3
shift 3
cells_extra_args=("$@")

[ -f "$gex_bam" ] || { echo "error: GEX BAM not found: $gex_bam" >&2; exit 1; }
[ -d "$demux_dir" ] || { echo "error: probe-demux output dir not found: $demux_dir" >&2; exit 1; }

shopt -s nullglob
barcode_files=("$demux_dir"/*.barcodes.txt)
shopt -u nullglob
[ ${#barcode_files[@]} -gt 0 ] || {
  echo "error: no *.barcodes.txt in $demux_dir -- run 'sbxpress probe-demux' first" >&2
  exit 1
}

mkdir -p "$out_dir"

for barcodes in "${barcode_files[@]}"; do
  sample=$(basename "$barcodes" .barcodes.txt)
  gex_raw="$out_dir/$sample/gex_raw"
  gex_outs="$out_dir/$sample/gex_outs"

  if [ ! -s "$barcodes" ]; then
    echo "== $sample: 0 assigned barcodes, skipping (nothing to call) =="
    continue
  fi
  n=$(wc -l < "$barcodes")
  echo "== $sample: $n assigned barcode(s) =="

  echo "-- matrix --"
  sbxpress matrix "$gex_bam" --cell-barcodes "$barcodes" --skip-unlisted-barcodes -o "$gex_raw"

  echo "-- cells --"
  sbxpress cells "$gex_raw/raw_feature_bc_matrix.h5" -o "$gex_outs" "${cells_extra_args[@]}"
done

echo "done: per-sample GEX matrices + cell-calling in $out_dir"
