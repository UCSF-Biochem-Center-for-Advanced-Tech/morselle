#!/usr/bin/env bash
# Map every DZEC03 XF:Z:cb_unmatched failure, retaining coordinate-sorted BAMs.
#
# R1 and R2 are mapped independently on purpose. R1 establishes whether the
# sequence replacing the expected CB/UMI is genomic/cDNA; R2 maps the other
# mate. Query names are retained, allowing mate-level joins downstream.
#
# This is a large run: ~827 million pairs. It needs substantial scratch space
# while samtools sorts each BAM, and produces two large BAMs.
#
# Usage:
#   bash scripts/map_dzec03_failed_cb_unmatched_full.sh
#
# Optional overrides:
#   THREADS=64 SORT_MEMORY=2G OUT_DIR=/path/to/output TMP_DIR=/path/to/scratch \
#     bash scripts/map_dzec03_failed_cb_unmatched_full.sh

set -euo pipefail

PROJECT=/mnt/analysis/active_projects/cat_tagseq/DZEC03/sbxpress_illumina
FAILED_DIR="$PROJECT/extract"
MMI=/mnt/analysis/genomes/GRCh38/GRCh38.primary_assembly.genome.mmi
OUT_DIR=${OUT_DIR:-"$PROJECT/failed_cb_unmatched_mapping_full"}
TMP_DIR=${TMP_DIR:-"$OUT_DIR/sort_tmp"}
THREADS=${THREADS:-32}
SORT_MEMORY=${SORT_MEMORY:-2G}

for program in gzip awk minimap2 samtools; do
    command -v "$program" >/dev/null || { echo "Missing required program: $program" >&2; exit 1; }
done
[[ -r "$MMI" ]] || { echo "Cannot read minimap2 index: $MMI" >&2; exit 1; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 1; }

mkdir -p "$OUT_DIR" "$TMP_DIR"
if compgen -G "$OUT_DIR/failed_cb_unmatched.*.coordinate.bam" >/dev/null; then
    echo "Refusing to overwrite existing BAMs in $OUT_DIR" >&2
    echo "Choose another OUT_DIR or remove the prior result deliberately." >&2
    exit 1
fi

stream_mate() {
    local mate=$1
    local shard input
    for shard in 001 002 003 004; do
        input="$FAILED_DIR/DZEC03.failed.shard-${shard}_${mate}.fastq.gz"
        [[ -r "$input" ]] || { echo "Cannot read $input" >&2; exit 1; }
        gzip -dc "$input"
    done | awk '
        NR % 4 == 1 { keep = ($0 ~ /XF:Z:cb_unmatched([[:space:]]|$)/) }
        keep { print }
    '
}

map_mate() {
    local mate=$1
    local bam="$OUT_DIR/failed_cb_unmatched.${mate}.coordinate.bam"
    echo "[$(date -Is)] Mapping all cb_unmatched $mate reads..." >&2
    stream_mate "$mate" \
        | minimap2 -a -x splice:sr --secondary=no -t "$THREADS" "$MMI" - \
        | samtools sort -@ "$THREADS" -m "$SORT_MEMORY" -T "$TMP_DIR/${mate}.sort" -o "$bam" -
    samtools index -@ "$THREADS" "$bam"
    echo "[$(date -Is)] Finished $mate." >&2
}

map_mate R1
map_mate R2

summary="$OUT_DIR/mapping_summary.tsv"
{
    printf 'mate\tprimary_reads\tmapped\tmapped_pct\tmapq_ge_20\tmapq_ge_20_pct\tspliced\tspliced_pct_of_mapped\n'
    for mate in R1 R2; do
        bam="$OUT_DIR/failed_cb_unmatched.${mate}.coordinate.bam"
        total=$(samtools view -c -F 2304 "$bam")
        mapped=$(samtools view -c -F 2308 "$bam")
        mapq20=$(samtools view -c -F 2308 -q 20 "$bam")
        spliced=$(samtools view -F 2308 "$bam" | awk '$6 ~ /[0-9]N/ { count++ } END { print count + 0 }')
        awk -v mate="$mate" -v total="$total" -v mapped="$mapped" -v mapq20="$mapq20" -v spliced="$spliced" \
            'BEGIN { printf "%s\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\n", mate, total, mapped, 100*mapped/total, mapq20, 100*mapq20/total, spliced, 100*spliced/mapped }'
    done
} > "$summary"

echo "[$(date -Is)] Done. Summary: $summary" >&2
if command -v column >/dev/null; then
    column -t -s $'\t' "$summary" >&2
else
    cat "$summary" >&2
fi
