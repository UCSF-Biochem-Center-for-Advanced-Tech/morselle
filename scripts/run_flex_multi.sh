#!/usr/bin/env bash
# FLEX_SUPPORT_PLAN.md Phase 3 (GEX) + Phase 5 (CRISPR) per-sample counting.
#
# Zero new Rust for the GEX side -- loops the existing, fully general
# `sbxpress matrix` and `sbxpress cells` commands over each sample's Phase 2
# `probe-demux` barcode list. The optional CRISPR side adds Phase 4's
# `sbxpress guides` output and Phase 5's `sbxpress protospacer-call`. Needs:
#   - a GEX-library BAM already through `sbxpress annotate` + `sbxpress molecules`
#     (carries CB/PB/xf, same input `probe-demux` itself consumes)
#   - `sbxpress probe-demux`'s output directory (one <sample_id>.barcodes.txt
#     per declared sample)
#   - optionally, a CRISPR-library guides BAM (`sbxpress guides` output, same
#     pooled multi-sample shape as the GEX BAM) and its `--features-out` TSV,
#     to also run Phase 5's CRISPR counting + protospacer calling per sample
#
# `--skip-unlisted-barcodes` on both `matrix` calls is load-bearing, not
# optional: both BAMs are pooled multi-sample libraries, so every per-sample
# `matrix` invocation legitimately sees other samples' barcodes too. Without
# that flag `matrix` (correctly, for its original single-sample use case)
# hard-errors on the first one it meets.
#
# CRISPR identity piggybacks on GEX-defined cell identity (matching real 10x
# biology, and FLEX_SUPPORT_PLAN.md Phase 5's own design): the CRISPR side
# uses each sample's Phase 3 GEX-*called* cell list (`gex_outs/filtered_barcodes.csv`,
# barcode column only -- not the pre-cell-calling probe-demux list), so a
# sample with no GEX cells called gets no CRISPR counting either.
#
# Usage:
#   scripts/run_flex_multi.sh [--crispr-bam BAM --guide-features TSV] \
#     <gex.bam> <probe_demux_out_dir> <output_dir> [extra `sbxpress cells` args...]
#
# Output, per declared sample:
#   <output_dir>/<sample_id>/gex_raw/{raw_feature_bc_matrix.h5,raw_feature_bc_matrix/}
#   <output_dir>/<sample_id>/gex_outs/{filtered_feature_bc_matrix.h5,filtered_barcodes.csv,...}
#   <output_dir>/<sample_id>/crispr_raw/...        (only with --crispr-bam)
#   <output_dir>/<sample_id>/crispr_outs/{filtered_feature_bc_matrix.h5,
#     protospacer_calls_per_cell.csv,protospacer_calls_summary.csv,
#     protospacer_umi_thresholds.csv,protospacer_umi_thresholds.json}

set -euo pipefail

usage() {
  echo "usage: $(basename "$0") [--crispr-bam BAM --guide-features TSV] <gex.bam> <probe_demux_out_dir> <output_dir> [extra 'sbxpress cells' args...]" >&2
  exit 1
}

crispr_bam=""
guide_features=""
positional=()
while [ $# -gt 0 ]; do
  case "$1" in
    --crispr-bam) crispr_bam=${2:?missing value for --crispr-bam}; shift 2 ;;
    --guide-features) guide_features=${2:?missing value for --guide-features}; shift 2 ;;
    --) shift; positional+=("$@"); break ;;
    *) positional+=("$1"); shift ;;
  esac
done
set -- "${positional[@]}"

[ $# -ge 3 ] || usage
gex_bam=$1
demux_dir=$2
out_dir=$3
shift 3
cells_extra_args=("$@")

[ -f "$gex_bam" ] || { echo "error: GEX BAM not found: $gex_bam" >&2; exit 1; }
[ -d "$demux_dir" ] || { echo "error: probe-demux output dir not found: $demux_dir" >&2; exit 1; }
if [ -n "$crispr_bam" ] || [ -n "$guide_features" ]; then
  [ -n "$crispr_bam" ] && [ -n "$guide_features" ] || {
    echo "error: --crispr-bam and --guide-features must be given together" >&2
    exit 1
  }
  [ -f "$crispr_bam" ] || { echo "error: CRISPR guides BAM not found: $crispr_bam" >&2; exit 1; }
  [ -f "$guide_features" ] || { echo "error: guide features TSV not found: $guide_features" >&2; exit 1; }
fi

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
  sample_dir="$out_dir/$sample"
  gex_raw="$sample_dir/gex_raw"
  gex_outs="$sample_dir/gex_outs"

  if [ ! -s "$barcodes" ]; then
    echo "== $sample: 0 assigned barcodes, skipping (nothing to call) =="
    continue
  fi
  n=$(wc -l < "$barcodes")
  echo "== $sample: $n assigned barcode(s) =="

  echo "-- matrix (GEX) --"
  sbxpress matrix "$gex_bam" --cell-barcodes "$barcodes" --skip-unlisted-barcodes -o "$gex_raw"

  echo "-- cells (GEX) --"
  sbxpress cells "$gex_raw/raw_feature_bc_matrix.h5" -o "$gex_outs" "${cells_extra_args[@]}"

  [ -n "$crispr_bam" ] || continue

  gex_cells="$sample_dir/gex_cells.txt"
  tail -n +2 "$gex_outs/filtered_barcodes.csv" | cut -d, -f2 > "$gex_cells"
  if [ ! -s "$gex_cells" ]; then
    echo "== $sample: 0 GEX cells called, skipping CRISPR counting =="
    continue
  fi

  crispr_molecules="$sample_dir/crispr_molecules.bam"
  crispr_raw="$sample_dir/crispr_raw"
  crispr_outs="$sample_dir/crispr_outs"

  echo "-- molecules (CRISPR, --tagseq) --"
  sbxpress molecules "$crispr_bam" --tagseq -o "$crispr_molecules"

  echo "-- matrix (CRISPR) --"
  sbxpress matrix "$crispr_molecules" --cell-barcodes "$gex_cells" --skip-unlisted-barcodes \
    --features "$guide_features" -o "$crispr_raw"

  echo "-- cells (CRISPR, --tagseq-barcodes) --"
  sbxpress cells "$crispr_raw/raw_feature_bc_matrix.h5" --tagseq-barcodes "$gex_cells" -o "$crispr_outs"

  echo "-- protospacer-call --"
  sbxpress protospacer-call "$crispr_outs/filtered_feature_bc_matrix.h5" -o "$crispr_outs"
done

echo "done: per-sample GEX matrices + cell-calling$([ -n "$crispr_bam" ] && echo " + CRISPR counting/protospacer calling") in $out_dir"
