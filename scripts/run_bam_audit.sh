#!/usr/bin/env bash
# Gather full-BAM audit data on the server.
#
# Runs `bam-audit` over the Illumina / SBX p5 / SBX p7 Cell Ranger BAMs and writes,
# per BAM, a small JSON summary plus a hash-subsampled molecule dump. Everything it
# produces is small enough to copy back for analysis (see bam_audit_report.py).
#
#   ./run_bam_audit.sh <illumina.bam> <p5.bam> <p7.bam> [outdir] [threads] [sample]
#
# `sample` is the molecule subsample rate (1 in N). The SAME N must be used for all
# three BAMs — the subsample is a hash of (CB,UMI), so identical N selects identical
# molecules in every BAM, which is what makes the dumps joinable. Recommended:
#   full BAMs  -> 256 (Illumina M2 has ~2.8B alignments; keeps dumps ~10-30 MB)
#   small/test -> 1   (dump every molecule)
set -euo pipefail

if [[ $# -lt 3 ]]; then
  sed -n '1,20p' "$0" | sed 's/^# \{0,1\}//'
  exit 1
fi

ILL=$1; P5=$2; P7=$3
OUT=${4:-bam_audit_out}
THREADS=${5:-8}
SAMPLE=${6:-256}

# Prefer an installed bam-audit; otherwise fall back to the release build.
if command -v bam-audit >/dev/null 2>&1; then
  AUDIT=bam-audit
elif [[ -x "$(dirname "$0")/../target/release/bam-audit" ]]; then
  AUDIT="$(cd "$(dirname "$0")/.." && pwd)/target/release/bam-audit"
else
  echo "bam-audit not found. Build it:  cargo build --release --bin bam-audit"
  exit 1
fi
echo "using: $AUDIT"
mkdir -p "$OUT"

run() {
  local bam=$1 label=$2 tag=$3
  [[ -f "$bam" ]] || { echo "missing BAM: $bam"; exit 1; }
  echo "=== $label : $bam ==="
  "$AUDIT" "$bam" \
    --label "$label" \
    --output "$OUT/${tag}.json" \
    --threads "$THREADS" \
    --molecules "$OUT/${tag}.molecules.tsv.gz" \
    --molecule-sample "$SAMPLE" \
    2>&1 | tee "$OUT/${tag}.log"
  echo
}

run "$ILL" "ILLUMINA" "illumina"
run "$P5"  "SBX p5"   "p5"
run "$P7"  "SBX p7"   "p7"

echo "Wrote to $OUT:"
ls -lh "$OUT"
echo
echo "Copy the whole directory back, then:  python3 bam_audit_report.py $OUT"
