#!/usr/bin/env bash
# Map a representative, deterministic sample of DZEC03 cb_unmatched failures.
#
# The failed FASTQ headers contain XF:Z:cb_unmatched. This script takes an equal,
# bounded sample from each failed-output shard, so R1 and R2 receive precisely the
# same set of pairs even though they are streamed separately. It then maps each
# mate independently. That distinction is intentional: R1 asks whether the
# putative barcode-side sequence is genomic; R2 asks whether the insert-side
# sequence is transcript/genome-derived.
#
# Usage:
#   bash scripts/map_dzec03_failed_cb_unmatched.sh
# Optional overrides:
#   THREADS=64 READS_PER_SHARD=50000 OUT_DIR=/path/to/output bash scripts/map_dzec03_failed_cb_unmatched.sh
#
# The default takes 25,000 pairs per shard (100,000 pairs total). It reads only
# enough of each large gzip file to reach that bound, rather than scanning all
# 116 GB of failed FASTQs before mapping.

set -euo pipefail

PROJECT=/mnt/analysis/active_projects/cat_tagseq/DZEC03/sbxpress_illumina
FAILED_DIR="$PROJECT/extract"
MMI=/mnt/analysis/genomes/GRCh38/GRCh38.primary_assembly.genome.mmi
OUT_DIR=${OUT_DIR:-"$PROJECT/failed_cb_unmatched_mapping_sample"}
THREADS=${THREADS:-32}
READS_PER_SHARD=${READS_PER_SHARD:-25000}

for program in gzip awk minimap2 samtools; do
    command -v "$program" >/dev/null || { echo "Missing required program: $program" >&2; exit 1; }
done
[[ -r "$MMI" ]] || { echo "Cannot read minimap2 index: $MMI" >&2; exit 1; }
[[ "$READS_PER_SHARD" =~ ^[1-9][0-9]*$ ]] || { echo "READS_PER_SHARD must be a positive integer" >&2; exit 1; }

mkdir -p "$OUT_DIR"
if compgen -G "$OUT_DIR/failed_cb_unmatched.*" >/dev/null; then
    echo "Refusing to overwrite existing results in $OUT_DIR" >&2
    echo "Choose another OUT_DIR or remove the prior sample deliberately." >&2
    exit 1
fi

tmp_dir=$(mktemp -d "$OUT_DIR/.tmp.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT

sample_mate() {
    local mate=$1
    local out_fastq=$2
    local shard input
    for shard in 001 002 003 004; do
        input="$FAILED_DIR/DZEC03.failed.shard-${shard}_${mate}.fastq.gz"
        [[ -r "$input" ]] || { echo "Cannot read $input" >&2; exit 1; }
        # awk exits once it has the requested number of matching pairs. gzip
        # then normally reports SIGPIPE (141); accept that expected status.
        set +o pipefail
        gzip -dc "$input" | awk -v max_reads="$READS_PER_SHARD" '
            NR % 4 == 1 {
                keep = ($0 ~ /XF:Z:cb_unmatched([[:space:]]|$)/)
            }
            keep {
                print
                if (NR % 4 == 0 && ++kept == max_reads) exit
            }
        '
        status=("${PIPESTATUS[@]}")
        set -o pipefail
        [[ "${status[0]}" == 0 || "${status[0]}" == 141 ]] || exit "${status[0]}"
        [[ "${status[1]}" == 0 ]] || exit "${status[1]}"
    done | gzip -1 > "$out_fastq"
}

echo "Sampling the first $READS_PER_SHARD cb_unmatched pairs from each failed shard..." >&2
sample_mate R1 "$tmp_dir/failed_cb_unmatched.R1.fastq.gz"
sample_mate R2 "$tmp_dir/failed_cb_unmatched.R2.fastq.gz"

r1_records=$(gzip -dc "$tmp_dir/failed_cb_unmatched.R1.fastq.gz" | awk 'END { print NR / 4 }')
r2_records=$(gzip -dc "$tmp_dir/failed_cb_unmatched.R2.fastq.gz" | awk 'END { print NR / 4 }')
[[ "$r1_records" == "$r2_records" ]] || {
    echo "Sample mate count mismatch: R1=$r1_records R2=$r2_records" >&2
    exit 1
}
[[ "$r1_records" -gt 0 ]] || { echo "No cb_unmatched reads selected" >&2; exit 1; }

map_mate() {
    local mate=$1
    local bam="$OUT_DIR/failed_cb_unmatched.${mate}.coordinate.bam"
    echo "Mapping $mate ($r1_records reads)..." >&2
    minimap2 -a -x splice:sr --secondary=no -t "$THREADS" "$MMI" \
        "$tmp_dir/failed_cb_unmatched.${mate}.fastq.gz" \
        | samtools sort -@ "$THREADS" -o "$bam" -
    samtools index -@ "$THREADS" "$bam"
}

map_mate R1
map_mate R2

summary="$OUT_DIR/mapping_summary.tsv"
{
    printf 'mate\treads_sampled\tmapped\tmapped_pct\tmapq_ge_20\tmapq_ge_20_pct\n'
    for mate in R1 R2; do
        bam="$OUT_DIR/failed_cb_unmatched.${mate}.coordinate.bam"
        # Exclude unmapped, secondary, and supplementary records. The latter are
        # alignments rather than extra sampled reads, so they must not inflate
        # the read-level summary.
        total=$r1_records
        mapped=$(samtools view -c -F 2308 "$bam")
        mapq20=$(samtools view -c -F 2308 -q 20 "$bam")
        awk -v mate="$mate" -v total="$total" -v mapped="$mapped" -v mapq20="$mapq20" \
            'BEGIN { printf "%s\t%d\t%d\t%.3f\t%d\t%.3f\n", mate, total, mapped, 100*mapped/total, mapq20, 100*mapq20/total }'
    done
} > "$summary"

echo "Done. Mapping summary: $summary" >&2
if command -v column >/dev/null; then
    column -t -s $'\t' "$summary" >&2
else
    cat "$summary" >&2
fi
