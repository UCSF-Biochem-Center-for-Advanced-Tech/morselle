#!/usr/bin/env bash
# docs/archive/FLEX_SUPPORT_PLAN.md Phase 3 (GEX) + Phase 5 (CRISPR) per-sample counting.
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
# `--probe-barcodes` on both `matrix` calls is also load-bearing, not
# optional: probe-demux's cell barcode is the composite (GEM, ProbeBarcode)
# pair (real Flex multiplex pooling -- see probe_demux.rs's module doc
# comment), so a GEM legitimately assigned to more than one sample now
# appears in more than one sample's <sample_id>.barcodes.txt. Without this
# flag, each sample's matrix would also pull in every OTHER sample's
# molecules for that shared GEM.
#
# CRISPR identity piggybacks on GEX-defined cell identity (matching real 10x
# biology, and docs/archive/FLEX_SUPPORT_PLAN.md Phase 5's own design): the CRISPR side
# uses each sample's Phase 3 GEX-*called* cell list (`gex_outs/filtered_barcodes.csv`,
# barcode column only -- not the pre-cell-calling probe-demux list), so a
# sample with no GEX cells called gets no CRISPR counting either.
#
# CRISPR `molecules` does NOT use `--tagseq`: that mode exists for a
# genuinely different protocol -- pooled, non-single-cell libraries where one
# barcode is a whole sample at tens-of-millions-of-reads depth (see
# `sbx10x-molecules --help`'s own doc comment on the flag). Real CAT06 CRISPR
# Guide Capture data is ordinary single-cell depth (real per-cell counts:
# 6-7 UMIs/cell against a 3-UMI call floor -- see `protospacer_calling.rs`),
# same order of magnitude as GEX, so `--tagseq`'s skip of cross-guide UMI
# collision adjudication *and* directional UMI error-correction was very
# likely hurting accuracy right at the threshold decision boundary rather
# than helping, not just mismatched by name.
#
# Usage:
#   scripts/run_flex_multi.sh [--crispr-bam BAM --guide-features TSV] \
#     [--jobs N] \
#     <gex.bam> <probe_demux_out_dir> <output_dir> [extra `sbxpress cells` args...]
#
# Output:
#   <output_dir>/crispr_molecules.bam              (only with --crispr-bam; shared
#                                                     across all samples -- see below)
#   <output_dir>/<sample_id>/gex_raw/{raw_feature_bc_matrix.h5,raw_feature_bc_matrix/}
#   <output_dir>/<sample_id>/gex_outs/{filtered_feature_bc_matrix.h5,filtered_barcodes.csv,...}
#   <output_dir>/<sample_id>/crispr_raw/...        (only with --crispr-bam)
#   <output_dir>/<sample_id>/crispr_outs/{filtered_feature_bc_matrix.h5,
#     protospacer_calls_per_cell.csv,protospacer_calls_summary.csv,
#     protospacer_umi_thresholds.csv,protospacer_umi_thresholds.json}
#
# `molecules "$crispr_bam"` runs ONCE, before the per-sample loop, not once
# per sample -- it processes the whole pooled multi-sample guide BAM
# regardless of which sample invokes it (real check on CAT06 CP1 data: the
# output BAM was byte-identical, same size, across 3 different samples' runs
# of the old per-sample-loop version), so repeating it per sample was pure
# waste: tens of real minutes redone 8 extra times for a 9-sample run, hours
# thrown away per library. Matches how the GEX side already builds its
# shared molecules BAM once (by the caller, before this script even starts)
# rather than per sample.
#
# Per-sample work (GEX matrix+cells, then CRISPR matrix+cells+protospacer-call)
# runs in PARALLEL across samples via `--jobs` (default: all samples at
# once), not sequentially. Each sample's own work only reads the shared,
# already-built BAMs (gex.bam / crispr_molecules.bam) and is I/O-bound and
# lightweight on CPU/memory (real measured ~1-2GB RSS per `matrix`-building
# process on CAT06), so concurrent samples are safe on any box with real
# headroom -- this turns an N-samples-times-one-scan sequential wall clock
# into roughly one scan's worth of wall clock. Lower --jobs on a
# memory/IO-constrained box.

set -euo pipefail

usage() {
  echo "usage: $(basename "$0") [--crispr-bam BAM --guide-features TSV] [--jobs N] <gex.bam> <probe_demux_out_dir> <output_dir> [extra 'sbxpress cells' args...]" >&2
  exit 1
}

crispr_bam=""
guide_features=""
jobs=""
positional=()
while [ $# -gt 0 ]; do
  case "$1" in
    --crispr-bam) crispr_bam=${2:?missing value for --crispr-bam}; shift 2 ;;
    --guide-features) guide_features=${2:?missing value for --guide-features}; shift 2 ;;
    --jobs) jobs=${2:?missing value for --jobs}; shift 2 ;;
    --) shift; positional+=("$@"); break ;;
    *) positional+=("$1"); shift ;;
  esac
done
set -- ${positional[@]+"${positional[@]}"}

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

crispr_molecules="$out_dir/crispr_molecules.bam"
if [ -n "$crispr_bam" ] && [ ! -f "$crispr_molecules" ]; then
  echo "-- molecules (CRISPR, whole pooled guide BAM, once) --"
  # --flex: real production bug found 2026-09-21 on CAT06 CP1 -- every
  # sample's real per-sample CRISPR molecule count came back exactly 0
  # despite real, substantial signal existing in the pooled library. Root
  # cause: without --flex, sbx10x-molecules groups UMI evidence by cell
  # barcode alone -- the exact cross-sample contamination bug already found
  # and fixed for the GEX side (RawKey now includes PB) never got applied
  # here, since this script's CRISPR molecules call was never updated to
  # pass the flag. 99.8% of real CBs are shared across pooled Flex samples,
  # so winner-take-all collision resolution without --flex silently hands
  # each shared CB+UMI's "representative" read to an arbitrary sample,
  # collapsing real per-sample signal to zero even though the raw reads and
  # their PB tags are correct.
  sbxpress molecules "$crispr_bam" -o "$crispr_molecules" --flex
fi

[ -n "$jobs" ] || jobs=${#barcode_files[@]}

run_sample() {
  set -euo pipefail
  local barcodes=$1
  local sample sample_dir gex_raw gex_outs probes n
  sample=$(basename "$barcodes" .barcodes.txt)
  sample_dir="$out_dir/$sample"
  gex_raw="$sample_dir/gex_raw"
  gex_outs="$sample_dir/gex_outs"
  probes="$demux_dir/$sample.probe_barcodes.txt"

  if [ ! -s "$barcodes" ]; then
    echo "== $sample: 0 assigned barcodes, skipping (nothing to call) =="
    return 0
  fi
  n=$(wc -l < "$barcodes")
  echo "== $sample: $n assigned barcode(s) =="

  local extra=() argument
  while IFS= read -r -d '' argument; do extra+=("$argument"); done < "$CELLS_EXTRA_ARGS_FILE"

  echo "-- $sample: matrix (GEX) --"
  sbxpress matrix "$gex_bam" --cell-barcodes "$barcodes" --probe-barcodes "$probes" --skip-unlisted-barcodes -o "$gex_raw"

  echo "-- $sample: cells (GEX) --"
  sbxpress cells "$gex_raw/raw_feature_bc_matrix.h5" -o "$gex_outs" ${extra[@]+"${extra[@]}"}

  [ -n "$crispr_bam" ] || { echo "== $sample: done =="; return 0; }

  local gex_cells crispr_raw crispr_outs
  gex_cells="$sample_dir/gex_cells.txt"
  tail -n +2 "$gex_outs/filtered_barcodes.csv" | cut -d, -f2 > "$gex_cells"
  if [ ! -s "$gex_cells" ]; then
    echo "== $sample: 0 GEX cells called, skipping CRISPR counting =="
    return 0
  fi

  crispr_raw="$sample_dir/crispr_raw"
  crispr_outs="$sample_dir/crispr_outs"

  echo "-- $sample: matrix (CRISPR) --"
  sbxpress matrix "$crispr_molecules" --cell-barcodes "$gex_cells" --probe-barcodes "$probes" --skip-unlisted-barcodes \
    --features "$guide_features" -o "$crispr_raw"

  echo "-- $sample: cells (CRISPR, --tagseq-barcodes) --"
  sbxpress cells "$crispr_raw/raw_feature_bc_matrix.h5" --tagseq-barcodes "$gex_cells" -o "$crispr_outs"

  echo "-- $sample: protospacer-call --"
  sbxpress protospacer-call "$crispr_outs/filtered_feature_bc_matrix.h5" -o "$crispr_outs"
  echo "== $sample: done =="
}
export -f run_sample
export out_dir demux_dir gex_bam crispr_bam guide_features crispr_molecules
# NUL-delimited arguments preserve whitespace and quoting across workers.
CELLS_EXTRA_ARGS_FILE=$(mktemp "$out_dir/.cells-args.XXXXXX")
trap 'rm -f "$CELLS_EXTRA_ARGS_FILE"' EXIT
if (( ${#cells_extra_args[@]} > 0 )); then
  printf '%s\0' "${cells_extra_args[@]}" > "$CELLS_EXTRA_ARGS_FILE"
fi
export CELLS_EXTRA_ARGS_FILE

printf '%s\0' "${barcode_files[@]}" | xargs -0 -P "$jobs" -I{} bash -c 'run_sample "$1"' _ {}

echo "done: per-sample GEX matrices + cell-calling$([ -n "$crispr_bam" ] && echo " + CRISPR counting/protospacer calling") in $out_dir"
