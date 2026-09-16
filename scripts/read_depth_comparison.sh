#!/usr/bin/env bash
# Real-vs-real read-depth comparison: 4 checkpoints matching Illumina's own
# reported funnel (CAT06_read_metrics_report.md) against Roche SBX, each
# checkpoint asking "how many raw Roche reads are needed to reach the same
# count Illumina reached at this stage?" and producing its own molecule/
# matrix output folder documenting the real raw-read count used.
#
# Checkpoints (Illumina targets, all real numbers from CAT06_read_metrics_report.md,
# updated 2026-09-16 once the real L2 Cell Ranger run completed):
#   1. Raw reads:                    6,852,776,830   (all 4 libraries)
#   2. Demultiplexed reads:          6,002,061,558   (all 4 libraries)
#   3. Valid cellBC+sampleBC+UMI:    4,831,665,536   (L1_GEX + L2_GEX + L1_CRI
#      -- L2_CRI was never run through Cell Ranger, so the report's own
#      "All processed (3 libraries)" total explicitly excludes it; this is a
#      real, exact Illumina scope, not GEX-only as an earlier plan version
#      assumed before the report was updated)
#   4. Final molecules: computed at run time by summing real
#      filtered_feature_bc_matrix.h5 files, GEX only (L1+L2 -- CRISPR guide
#      UMI counts are a different kind of "molecule" than gene-expression
#      transcript counts, so CP4 stays GEX-scoped even though CP3 now
#      includes L1_CRI). As of 2026-09-16 both real matrices exist
#      (illumina_data/RZ16443_L1_multi/outs/... and
#      RZ16443_L2_multi_complete/outs/...): real combined target is
#      3,704,554,258 (L1 2,115,131,488 + L2 1,589,422,770, both measured
#      directly from the real matrices, not estimated). This script still
#      checks both L2 matrix paths at run time and falls back to L1-only
#      with an explicit scope note if either ever isn't there.
#
# Roche-side yield rates needed to invert "Illumina target -> raw reads
# needed" come from two already-running real full-scale jobs, not new
# calibration compute:
#   - CP2, CP3 numerators: sbxpress_extract_v2/metrics_summary.json
#     (assigned_early, probe_barcode.hit, per-library per_sample_early)
#   - CP4 numerator: sbxpress_full_run/molecules/{L1,L2}_gex.metrics.json
#     (molecules_counted)
#
# CP3's real Illumina target excludes exactly one library (L2_CRI), but
# extract's probe_hit counter is library-agnostic (probe barcodes are read
# on every library together). Approximated by scaling probe_hit by the
# three included libraries' (L1_GEX+L2_GEX+L1_CRI) real share of
# assigned_early -- probe-capture rate has no reason to differ structurally
# by library of the same physical run. Stated here and in the final report
# as an approximation, not represented as exact (confirmed acceptable with
# the user on 2026-09-15 rather than building new per-library
# instrumentation and re-running the multi-hour extract once more).
#
# This script does NOT try to deduplicate processing across checkpoints by
# sharing partial extract output -- each checkpoint runs the real, full,
# independent pipeline (extract -> probeset-align/guides -> molecules ->
# probe-demux -> matrix/cells -> protospacer-call) on its own real shard
# subset. Shard subsets ARE nested (checkpoint N's shards are a strict
# prefix of checkpoint N+1's, in the same fixed shard order) so the four
# runs are directly comparable, but they are NOT deduplicated: re-running a
# smaller checkpoint's shards again inside a larger checkpoint's extract is
# real, bounded, redundant compute, traded deliberately for simplicity and
# per-checkpoint independence/inspectability over a more complex shared-
# extract-output design.
#
# Usage: scripts/read_depth_comparison.sh
# Intended invocation (survives disconnect):
#   setsid nohup scripts/read_depth_comparison.sh > read_depth_comparison_driver.log 2>&1 < /dev/null &
#   disown

set -euo pipefail

CAT06=/mnt/analysis/active_projects/roche_SBX_ABRF/CAT06
BIN=/mnt/analysis/software_shared/sbxpress/target/release
OUT=$CAT06/read_depth_comparison
CORES=128

SAMPLE_SHEET="$CAT06/20260806_LH00416_0423_B23WVVYLT4/RZ16443/Reports/SampleSheet.csv"
WHITELIST=/mnt/analysis/software/cellranger-10.0.0/lib/python/cellranger/barcodes/737K-flex-v2.txt.gz
PB_WHITELIST="$CAT06/resources/flex-v2-384.txt"
PROBE_SET="$CAT06/resources/Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv"
FEATURE_REF="$CAT06/feature_reference_RZ.csv"
PROBE_SAMPLE_SHEET="$CAT06/illumina_data/RZ16443_L1_multi/outs/config.csv"
RAW_FASTQS="$CAT06/raw_fastqs"

EXTRACT_V2_METRICS="$CAT06/sbxpress_extract_v2/metrics_summary.json"
L1_MOLECULES_METRICS="$CAT06/sbxpress_full_run/molecules/L1_gex.metrics.json"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

mkdir -p "$OUT"

# ---------------------------------------------------------------------------
# Step 0: prerequisite check -- refuse to guess, wait for the real numbers.
# ---------------------------------------------------------------------------
check_prereqs() {
  local missing=0
  for f in "$EXTRACT_V2_METRICS" "$L1_MOLECULES_METRICS"; do
    if [ ! -f "$f" ]; then
      echo "not ready: $f does not exist yet" >&2
      missing=1
    fi
  done
  if [ "$missing" = 1 ]; then
    echo "error: prerequisite run(s) not finished. Re-run this script once both exist." >&2
    exit 1
  fi
}

# ---------------------------------------------------------------------------
# Step 1: real yield rates + raw-read targets, computed once prereqs exist.
# ---------------------------------------------------------------------------
compute_targets() {
  python3 - "$EXTRACT_V2_METRICS" "$L1_MOLECULES_METRICS" "$CAT06" <<'PYEOF' > "$OUT/targets.json"
import json, sys, os, subprocess

extract = json.load(open(sys.argv[1]))
molecules_l1 = json.load(open(sys.argv[2]))
cat06 = sys.argv[3]

def matrix_molecule_count(h5_path):
    if not os.path.isfile(h5_path):
        return None
    script = (
        "import h5py,sys;"
        f"f=h5py.File('{h5_path}','r');"
        "print(int(f['matrix']['data'][:].astype('int64').sum()))"
    )
    return int(subprocess.check_output(["python3", "-c", script]).strip())

total_raw = extract["total_reads"]
assigned_early = extract["demux"]["assigned_early"]
probe_hit = extract["probe_barcode"]["hit"]
per_sample_early = {s["sample"]: s["pairs"] for s in extract["demux"]["per_sample_early"]}

# CP3 scope: every library except L2_CRI (Cell Ranger never ran on it) --
# L1_GEX + L2_GEX + L1_CRI. See script header for the approximation this
# feeds (probe_hit is library-agnostic; scaled by this share).
cp3_libs_early = sum(v for k, v in per_sample_early.items() if k != "TSeq100_L2_CRI")
cp3_share = cp3_libs_early / assigned_early
probe_hit_cp3_approx = probe_hit * cp3_share  # approximation -- see script header

l1_gex_early = per_sample_early.get("TSeq100_L1_GEX", 0)
l2_gex_early = per_sample_early.get("TSeq100_L2_GEX", 0)
molecules_counted_l1 = molecules_l1["molecules_counted"]

# CP4: L2 real Illumina matrix + L2 real Roche molecules, only if BOTH exist
# -- never combine a real L1 number with a placeholder/estimated L2 one.
l2_illumina_matrix_candidates = [
    os.path.join(cat06, "illumina_data", "RZ16443_L2_multi_complete", "outs", "filtered_feature_bc_matrix.h5"),
    os.path.join(cat06, "illumina_data", "RZ16443_L2_multi", "outs", "filtered_feature_bc_matrix.h5"),
]
l2_illumina_molecules = None
for cand in l2_illumina_matrix_candidates:
    l2_illumina_molecules = matrix_molecule_count(cand)
    if l2_illumina_molecules is not None:
        break

l2_molecules_metrics_path = os.path.join(cat06, "sbxpress_full_run", "molecules", "L2_gex.metrics.json")
molecules_counted_l2 = None
if os.path.isfile(l2_molecules_metrics_path):
    molecules_counted_l2 = json.load(open(l2_molecules_metrics_path))["molecules_counted"]

cp4_l1_only_illumina = 2115131488
if l2_illumina_molecules is not None and molecules_counted_l2 is not None:
    cp4_scope = "GEX only (L1+L2, real)"
    cp4_illumina_target = cp4_l1_only_illumina + l2_illumina_molecules
    cp4_gex_early = l1_gex_early + l2_gex_early
    cp4_molecules_counted = molecules_counted_l1 + molecules_counted_l2
else:
    cp4_scope = "GEX only (L1 only -- L2 not yet available on " + (
        "the Illumina side" if l2_illumina_molecules is None else "the Roche side"
    ) + ")"
    cp4_illumina_target = cp4_l1_only_illumina
    cp4_gex_early = l1_gex_early
    cp4_molecules_counted = molecules_counted_l1

cp4_gex_share_of_raw = cp4_gex_early / total_raw
r4 = cp4_molecules_counted / cp4_gex_early

r2 = assigned_early / total_raw
r3 = probe_hit_cp3_approx / total_raw

illumina = {
    "cp1_raw": 6852776830,
    "cp2_demuxed": 6002061558,
    "cp3_cb_samplebc_umi": 4831665536,
    "cp3_scope": "L1_GEX + L2_GEX + L1_CRI (excludes L2_CRI, never run through Cell Ranger)",
    "cp4_final_molecules": cp4_illumina_target,
    "cp4_scope": cp4_scope,
}

cp1_raw_needed = illumina["cp1_raw"]
cp2_raw_needed = illumina["cp2_demuxed"] / r2
cp3_raw_needed = illumina["cp3_cb_samplebc_umi"] / r3
# CP4: need enough raw reads so that the GEX-assigned share alone yields
# enough molecules -- raw_needed * cp4_gex_share_of_raw * r4 = target
cp4_raw_needed = illumina["cp4_final_molecules"] / (cp4_gex_share_of_raw * r4)

out = {
    "roche_real": {
        "total_raw": total_raw,
        "assigned_early": assigned_early,
        "probe_hit_all_libraries": probe_hit,
        "cp3_share_of_assigned_early": cp3_share,
        "probe_hit_cp3_approx": probe_hit_cp3_approx,
        "cp4_gex_share_of_raw": cp4_gex_share_of_raw,
        "cp4_molecules_counted": cp4_molecules_counted,
        "molecules_per_gex_assigned_read": r4,
    },
    "illumina_targets": illumina,
    "raw_reads_needed": {
        "cp1": round(cp1_raw_needed),
        "cp2": round(cp2_raw_needed),
        "cp3": round(cp3_raw_needed),
        "cp4": round(cp4_raw_needed),
    },
}
print(json.dumps(out, indent=2))
PYEOF
  cat "$OUT/targets.json"
}

# ---------------------------------------------------------------------------
# Step 2: real per-shard read counts (cached -- computed once, reused every
# time this script runs). Needed to select an exact, real shard subset per
# checkpoint rather than assuming uniform shard size.
# ---------------------------------------------------------------------------
build_shard_counts() {
  local cache="$OUT/shard_read_counts.tsv"
  if [ -s "$cache" ]; then
    log "shard read counts: cached ($cache)"
    return
  fi
  log "shard read counts: computing (754 shards, ${CORES}-way parallel, one real scan)"
  ls "$RAW_FASTQS"/*.fastq.gz | xargs -P "$CORES" -I{} bash -c '
    f="{}"
    n=$(( $(zcat "$f" | wc -l) / 4 ))
    echo -e "$f\t$n"
  ' > "$cache.tmp"
  mv "$cache.tmp" "$cache"
  log "shard read counts: done"
}

# Select the shortest prefix (in a fixed, deterministic shard order) whose
# cumulative real read count reaches >= $1. Writes the shard list to $2 and
# echoes the exact real total.
select_shard_subset() {
  local target=$1 out_list=$2
  python3 - "$OUT/shard_read_counts.tsv" "$target" "$out_list" <<'PYEOF'
import sys
counts_path, target, out_list = sys.argv[1], int(sys.argv[2]), sys.argv[3]
rows = []
with open(counts_path) as f:
    for line in f:
        path, n = line.rstrip("\n").split("\t")
        rows.append((path, int(n)))
rows.sort(key=lambda r: r[0])  # fixed, deterministic order -> real nesting across checkpoints
total = 0
chosen = []
for path, n in rows:
    if total >= target:
        break
    chosen.append(path)
    total += n
with open(out_list, "w") as f:
    f.write("\n".join(chosen) + "\n")
grand_total = sum(n for _, n in rows)
if total < target:
    print(f"SHORTFALL: only {total} of {target} needed reads available across all "
          f"{len(rows)} real shards ({grand_total} total raw reads on disk) -- this "
          f"checkpoint cannot be satisfied with the available Roche dataset.",
          file=sys.stderr)
    print("SHORTFALL", total)
else:
    print("OK", total)
PYEOF
}

# ---------------------------------------------------------------------------
# Step 3: run the real, full pipeline on one checkpoint's shard subset.
# Mirrors scripts/run_full_cat06.sh, generalized to a specific shard list and
# output directory instead of the whole raw_fastqs directory.
# ---------------------------------------------------------------------------
run_checkpoint() {
  local cp=$1 shard_list=$2 raw_used=$3
  local cpdir="$OUT/$cp"
  local logdir="$cpdir/logs"
  mkdir -p "$cpdir"/{extracted,aligned,molecules,demux,results,logs}
  echo "$raw_used" > "$cpdir/raw_reads_used.txt"
  log "$cp: starting on $(wc -l < "$shard_list") shards, $raw_used real raw reads"

  if [ ! -f "$cpdir/extracted/.done" ]; then
    "$BIN/sbx10x" \
      --fastqs $(cat "$shard_list") \
      -o "$cpdir/extracted" \
      --bclconvert-samplesheet "$SAMPLE_SHEET" \
      --whitelist "$WHITELIST" \
      --probe-barcode-whitelist "$PB_WHITELIST" \
      --minimap2 --threads "$CORES" --split "$CORES" \
      > "$logdir/extract.log" 2>&1
    touch "$cpdir/extracted/.done"
  fi
  log "$cp: extract done"

  align_one() {
    local mode=$1 fq=$2 out_bam=$3 manifest=$4
    [ -f "$out_bam" ] && return 0
    case "$mode" in
      GEX) "$BIN/sbxpress-probeset-align" "$fq" --probe-set "$PROBE_SET" --genes-out "$manifest" -o "$out_bam" --metrics "$out_bam.metrics.json" ;;
      CRI) "$BIN/sbxpress-guides" "$fq" --feature-reference "$FEATURE_REF" --features-out "$manifest" -o "$out_bam" --metrics "$out_bam.metrics.json" ;;
    esac
  }
  export -f align_one
  export BIN PROBE_SET FEATURE_REF

  if [ ! -f "$cpdir/aligned/.done" ]; then
    mkdir -p "$cpdir/aligned"/{TSeq100_L1_GEX,TSeq100_L2_GEX,TSeq100_L1_CRI,TSeq100_L2_CRI}
    local genes_out="$cpdir/aligned/genes.tsv" features_out="$cpdir/aligned/guide_features.tsv"
    local first_gex first_cri
    first_gex=$(ls "$cpdir"/extracted/TSeq100_L1_GEX_S*_R2_001.fastq.gz 2>/dev/null | head -1)
    [ -n "$first_gex" ] && [ ! -f "$genes_out" ] && align_one GEX "$first_gex" "$cpdir/aligned/TSeq100_L1_GEX/$(basename "$first_gex" .fastq.gz).bam" "$genes_out"
    first_cri=$(ls "$cpdir"/extracted/TSeq100_L1_CRI_S*_R2_001.fastq.gz 2>/dev/null | head -1)
    [ -n "$first_cri" ] && [ ! -f "$features_out" ] && align_one CRI "$first_cri" "$cpdir/aligned/TSeq100_L1_CRI/$(basename "$first_cri" .fastq.gz).bam" "$features_out"

    {
      for fq in "$cpdir"/extracted/TSeq100_L1_GEX_S*_R2_*.fastq.gz "$cpdir"/extracted/TSeq100_L2_GEX_S*_R2_*.fastq.gz; do
        [ -f "$fq" ] || continue
        lib=$(basename "$fq" | sed -E 's/_S[0-9]+_.*//')
        echo "GEX|$fq|$cpdir/aligned/$lib/$(basename "$fq" .fastq.gz).bam|$genes_out"
      done
      for fq in "$cpdir"/extracted/TSeq100_L1_CRI_S*_R2_*.fastq.gz "$cpdir"/extracted/TSeq100_L2_CRI_S*_R2_*.fastq.gz; do
        [ -f "$fq" ] || continue
        lib=$(basename "$fq" | sed -E 's/_S[0-9]+_.*//')
        echo "CRI|$fq|$cpdir/aligned/$lib/$(basename "$fq" .fastq.gz).bam|$features_out"
      done
    } | xargs -P "$CORES" -I{} bash -c '
        IFS="|" read -r mode fq out manifest <<< "{}"
        align_one "$mode" "$fq" "$out" "$manifest"
      ' > "$logdir/align.log" 2>&1

    for lib in TSeq100_L1_GEX TSeq100_L2_GEX TSeq100_L1_CRI TSeq100_L2_CRI; do
      [ -d "$cpdir/aligned/$lib" ] && ls "$cpdir/aligned/$lib"/*.bam >/dev/null 2>&1 && \
        samtools cat -o "$cpdir/aligned/$lib.bam" "$cpdir/aligned/$lib"/*.bam
    done
    touch "$cpdir/aligned/.done"
  fi
  log "$cp: align done"

  for lib_n in 1 2; do
    local gex_bam="$cpdir/aligned/TSeq100_L${lib_n}_GEX.bam"
    local cri_bam="$cpdir/aligned/TSeq100_L${lib_n}_CRI.bam"
    [ -f "$gex_bam" ] || continue
    local gex_mol="$cpdir/molecules/L${lib_n}_gex.bam"
    local demux_dir="$cpdir/demux/L${lib_n}"
    local results_dir="$cpdir/results/L${lib_n}"
    [ -f "$results_dir/.done" ] && continue

    if [ ! -f "$gex_mol" ]; then
      "$BIN/sbx10x-molecules" "$gex_bam" -o "$gex_mol" --metrics "$cpdir/molecules/L${lib_n}_gex.metrics.json" > "$logdir/molecules_L${lib_n}.log" 2>&1
    fi
    mkdir -p "$demux_dir"
    if [ -z "$(ls -A "$demux_dir" 2>/dev/null)" ]; then
      "$BIN/sbxpress-probe-demux" "$gex_mol" --probe-sample-sheet "$PROBE_SAMPLE_SHEET" -o "$demux_dir" > "$logdir/probe_demux_L${lib_n}.log" 2>&1
    fi
    if [ -f "$cri_bam" ]; then
      SBXPRESS_LIBEXEC="$BIN" PATH="$BIN:$PATH" \
        /mnt/analysis/software_shared/sbxpress/scripts/run_flex_multi.sh \
        --crispr-bam "$cri_bam" --guide-features "$cpdir/aligned/guide_features.tsv" \
        "$gex_mol" "$demux_dir" "$results_dir" > "$logdir/run_flex_multi_L${lib_n}.log" 2>&1
    else
      SBXPRESS_LIBEXEC="$BIN" PATH="$BIN:$PATH" \
        /mnt/analysis/software_shared/sbxpress/scripts/run_flex_multi.sh \
        "$gex_mol" "$demux_dir" "$results_dir" > "$logdir/run_flex_multi_L${lib_n}.log" 2>&1
    fi
    touch "$results_dir/.done"
  done
  log "$cp: done"
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
check_prereqs
compute_targets
build_shard_counts

for cp in cp1 cp2 cp3 cp4; do
  target=$(python3 -c "import json; print(json.load(open('$OUT/targets.json'))['raw_reads_needed']['$cp'])")
  list="$OUT/${cp}_shards.txt"
  actual=$(select_shard_subset "$target" "$list")
  log "$cp: target=$target real_shard_subset_total=$actual ($(wc -l < "$list") shards)"
  run_checkpoint "$cp" "$list" "$actual"
done

log "=== all checkpoints complete -- see $OUT/targets.json and each cpN/raw_reads_used.txt ==="
