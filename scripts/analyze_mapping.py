#!/usr/bin/env python3
"""
Descriptive mapping-ambiguity analysis for the two-arm CAT02 M2 evaluation.

Three analyses from docs/archive/mapping.md §3:
  1. NH-by-length  — primary-record NH=1 fraction by read-length bin and GX status
  2. secondary classification  — per-secondary-hit class by length bin (requires name-sort)
  3. score gaps  — primary/secondary AS and MAPQ differences by length bin (requires name-sort)

Usage:
  analyze_mapping.py --bam MULTIMAP_ANNOTATED_BAM --out OUTPUT_DIR [options]

Options:
  --bam FILE       Coordinate-sorted annotated multimap BAM (required)
  --baseline FILE  Coordinate-sorted annotated baseline BAM (for side-by-side counts)
  --out DIR        Output directory for TSV files and summary [default: ./mapping_analysis]
  --threads N      samtools threads for name-sort [default: 4]
  --skip-namesort  Skip secondary/score analyses (coordinate-sorted BAM only)
  --help
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

# Length bins proposed from M2 GEX R2 distribution:
#   min=50  p25=279  med=309  p75=347  p90=387  p99=477  max=649
LENGTH_BINS = [0, 250, 300, 350, 400, float("inf")]
LENGTH_LABELS = ["≤250", "251-300", "301-350", "351-400", "≥401"]


def length_bin(n: int) -> str:
    for lo, hi, label in zip(LENGTH_BINS, LENGTH_BINS[1:], LENGTH_LABELS):
        if lo < n <= hi:
            return label
    return LENGTH_LABELS[-1]


def clipped_bases(read: pysam.AlignedSegment) -> int:
    """Return total soft+hard clipped bases."""
    if read.cigartuples is None:
        return 0
    return sum(l for op, l in read.cigartuples if op in (4, 5))


def aligned_bases(read: pysam.AlignedSegment) -> int:
    """Query bases consumed by alignment (M, I, S, =, X operations on query)."""
    if read.cigartuples is None:
        return 0
    return sum(l for op, l in read.cigartuples if op in (0, 1, 4, 7, 8))


# ── Analysis 1: primary record NH-by-length ───────────────────────────────────

def analyze_primaries(bam_path: str) -> pd.DataFrame:
    """
    For every primary record, record: length bin, NH, GX assigned, mm tag.
    Returns a DataFrame with one row per primary read.
    """
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            nh = read.get_tag("NH") if read.has_tag("NH") else 1
            gx = read.get_tag("GX") if read.has_tag("GX") else ""
            mm = read.get_tag("mm") if read.has_tag("mm") else ""
            rlen = read.query_length or 0
            mapq = read.mapping_quality
            clip = clipped_bases(read)
            rows.append({
                "length_bin": length_bin(rlen),
                "read_length": rlen,
                "NH": nh,
                "nh1": int(nh == 1),
                "gx_assigned": int(bool(gx) and gx != "-"),
                "mm_rescue": int(bool(mm) and mm not in ("", ".")),
                "mapq": mapq,
                "clip_frac": clip / max(rlen, 1),
            })
    return pd.DataFrame(rows)


def summarize_primaries(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate primary records by length bin."""
    cat = pd.CategoricalDtype(categories=LENGTH_LABELS, ordered=True)
    df["length_bin"] = df["length_bin"].astype(cat)
    grp = df.groupby("length_bin", observed=True)
    summary = grp.agg(
        reads=("nh1", "count"),
        nh1=("nh1", "sum"),
        gx_assigned=("gx_assigned", "sum"),
        mm_rescue=("mm_rescue", "sum"),
        mean_mapq=("mapq", "mean"),
        mean_clip_frac=("clip_frac", "mean"),
        mean_length=("read_length", "mean"),
    ).reset_index()
    summary["nh1_frac"] = summary["nh1"] / summary["reads"]
    summary["gx_frac"] = summary["gx_assigned"] / summary["reads"]
    summary["mm_frac"] = summary["mm_rescue"] / summary["reads"]
    return summary


def nh_distribution(df: pd.DataFrame) -> pd.DataFrame:
    """NH value distribution by length bin (top-coded at 11+)."""
    cat = pd.CategoricalDtype(categories=LENGTH_LABELS, ordered=True)
    df = df.copy()
    df["length_bin"] = df["length_bin"].astype(cat)
    df["NH_bin"] = df["NH"].clip(upper=11).astype(str)
    df.loc[df["NH"] >= 11, "NH_bin"] = "11+"
    pivot = (
        df.groupby(["length_bin", "NH_bin"], observed=True)
        .size()
        .rename("count")
        .reset_index()
    )
    total = df.groupby("length_bin", observed=True).size().rename("total").reset_index()
    pivot = pivot.merge(total, on="length_bin")
    pivot["fraction"] = pivot["count"] / pivot["total"]
    return pivot


# ── Analysis 2: secondary hit classification (requires name-sorted BAM) ───────

def classify_secondaries(bam_path: str) -> pd.DataFrame:
    """
    For every NH>1 query, classify each secondary alignment relative to the
    primary gene assignment:
      - same_gene: GX matches primary GX
      - diff_gene: GX is a different gene
      - unannotated: no GX on secondary
      - intergenic: mapped but RE tag indicates intergenic/other

    Processes a name-sorted BAM where all records for a query are consecutive.
    """
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        query_group: list[pysam.AlignedSegment] = []
        current_qname = None

        def process_group(group):
            if not group:
                return
            primaries = [r for r in group if not r.is_secondary and not r.is_supplementary and not r.is_unmapped]
            secondaries = [r for r in group if r.is_secondary and not r.is_unmapped]
            if not primaries or not secondaries:
                return
            primary = primaries[0]
            nh = primary.get_tag("NH") if primary.has_tag("NH") else 1
            if nh <= 1:
                return
            p_gx = primary.get_tag("GX") if primary.has_tag("GX") else ""
            p_len = primary.query_length or 0
            p_as = primary.get_tag("AS") if primary.has_tag("AS") else None
            p_mapq = primary.mapping_quality

            for sec in secondaries:
                s_gx = sec.get_tag("GX") if sec.has_tag("GX") else ""
                s_as = sec.get_tag("AS") if sec.has_tag("AS") else None
                s_mapq = sec.mapping_quality
                s_chrom = sec.reference_name or ""

                if not s_gx or s_gx == "-":
                    if s_chrom.startswith("chr") and "random" not in s_chrom and "Un_" not in s_chrom:
                        sec_class = "unannotated_mapped"
                    else:
                        sec_class = "decoy_other"
                elif s_gx == p_gx:
                    sec_class = "same_gene"
                else:
                    sec_class = "diff_gene"

                as_gap = (p_as - s_as) if (p_as is not None and s_as is not None) else None
                rows.append({
                    "query": primary.query_name,
                    "length_bin": length_bin(p_len),
                    "read_length": p_len,
                    "NH": nh,
                    "primary_gx": p_gx,
                    "secondary_gx": s_gx,
                    "sec_class": sec_class,
                    "primary_AS": p_as,
                    "secondary_AS": s_as,
                    "as_gap": as_gap,
                    "primary_mapq": p_mapq,
                    "secondary_mapq": s_mapq,
                    "primary_chrom": primary.reference_name or "",
                    "secondary_chrom": s_chrom,
                })

        for read in bam.fetch(until_eof=True):
            if read.query_name != current_qname:
                process_group(query_group)
                query_group = []
                current_qname = read.query_name
            query_group.append(read)
        process_group(query_group)

    return pd.DataFrame(rows)


def summarize_secondaries(df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate secondary classifications by length bin."""
    if df.empty:
        return df
    cat = pd.CategoricalDtype(categories=LENGTH_LABELS, ordered=True)
    df = df.copy()
    df["length_bin"] = df["length_bin"].astype(cat)
    pivot = (
        df.groupby(["length_bin", "sec_class"], observed=True)
        .size()
        .rename("count")
        .reset_index()
    )
    total = df.groupby("length_bin", observed=True).size().rename("total").reset_index()
    pivot = pivot.merge(total, on="length_bin")
    pivot["fraction"] = pivot["count"] / pivot["total"]
    return pivot


def summarize_score_gaps(df: pd.DataFrame) -> pd.DataFrame:
    """AS gap and MAPQ difference statistics by length bin and secondary class."""
    if df.empty or "as_gap" not in df.columns:
        return df
    cat = pd.CategoricalDtype(categories=LENGTH_LABELS, ordered=True)
    df = df.copy()
    df["length_bin"] = df["length_bin"].astype(cat)
    valid = df.dropna(subset=["as_gap"])
    if valid.empty:
        print("  note: no AS tags found in BAM; score gap analysis skipped")
        return pd.DataFrame()
    grp = valid.groupby(["length_bin", "sec_class"], observed=True)["as_gap"]
    summary = grp.agg(["count", "mean", "median",
                        lambda x: np.percentile(x, 25),
                        lambda x: np.percentile(x, 75),
                        lambda x: np.percentile(x, 5)]).reset_index()
    summary.columns = ["length_bin", "sec_class", "n", "mean_gap",
                        "median_gap", "p25_gap", "p75_gap", "p5_gap"]
    return summary


# ── name-sort helper ──────────────────────────────────────────────────────────

def namesort_bam(bam_path: str, threads: int, tmp_dir: str) -> str:
    """Return path to a name-sorted copy of bam_path in tmp_dir."""
    out = os.path.join(tmp_dir, "namesorted.bam")
    cmd = ["samtools", "sort", "-n", f"-@{threads}", "-o", out, bam_path]
    print(f"  name-sorting BAM: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    return out


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam", required=True, help="Coordinate-sorted multimap annotated BAM")
    ap.add_argument("--baseline", default=None, help="Coordinate-sorted baseline annotated BAM (optional)")
    ap.add_argument("--out", default="mapping_analysis", help="Output directory [mapping_analysis]")
    ap.add_argument("--threads", type=int, default=4, help="samtools threads for name-sort [4]")
    ap.add_argument("--skip-namesort", action="store_true",
                    help="Skip secondary/score analyses (analysis 1 only)")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)

    # ── Analysis 1: primary records ───────────────────────────────────────────
    print(f"[1/3] analyzing primary records: {args.bam}")
    df_prim = analyze_primaries(args.bam)
    print(f"  {len(df_prim):,} primary mapped records")

    prim_summary = summarize_primaries(df_prim)
    out_prim = os.path.join(args.out, "nh_by_length.tsv")
    prim_summary.to_csv(out_prim, sep="\t", index=False, float_format="%.4f")
    print(f"  wrote {out_prim}")

    nh_dist = nh_distribution(df_prim)
    out_nhdist = os.path.join(args.out, "nh_distribution.tsv")
    nh_dist.to_csv(out_nhdist, sep="\t", index=False, float_format="%.4f")
    print(f"  wrote {out_nhdist}")

    if args.baseline:
        print(f"  analyzing baseline arm: {args.baseline}")
        df_base = analyze_primaries(args.baseline)
        base_summary = summarize_primaries(df_base)
        out_base = os.path.join(args.out, "nh_by_length_baseline.tsv")
        base_summary.to_csv(out_base, sep="\t", index=False, float_format="%.4f")
        print(f"  wrote {out_base}")

    # ── Analysis 2+3: secondary records (need name-sort) ─────────────────────
    if args.skip_namesort:
        print("[2/3] skipping secondary analysis (--skip-namesort)")
        df_sec = pd.DataFrame()
    else:
        with tempfile.TemporaryDirectory(prefix="analyze_mapping_") as tmp:
            ns_bam = namesort_bam(args.bam, args.threads, tmp)
            print(f"[2/3] classifying secondary hits")
            df_sec = classify_secondaries(ns_bam)

        if df_sec.empty:
            print("  no NH>1 primary records found (or no secondaries)")
        else:
            print(f"  {len(df_sec):,} secondary records classified")
            sec_summary = summarize_secondaries(df_sec)
            out_sec = os.path.join(args.out, "secondary_classes.tsv")
            sec_summary.to_csv(out_sec, sep="\t", index=False, float_format="%.4f")
            print(f"  wrote {out_sec}")

            print(f"[3/3] computing alignment score gaps")
            gap_summary = summarize_score_gaps(df_sec)
            if not gap_summary.empty:
                out_gaps = os.path.join(args.out, "score_gaps.tsv")
                gap_summary.to_csv(out_gaps, sep="\t", index=False, float_format="%.4f")
                print(f"  wrote {out_gaps}")

    # ── Text summary ──────────────────────────────────────────────────────────
    summary_path = os.path.join(args.out, "summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"mapping ambiguity analysis\n")
        f.write(f"bam: {args.bam}\n\n")

        f.write("NH=1 fraction by read-length bin (multimap arm)\n")
        f.write(f"{'bin':<12} {'reads':>10} {'NH=1':>10} {'NH=1%':>8} {'GX%':>8} {'mm%':>8}\n")
        for _, row in prim_summary.iterrows():
            f.write(f"{row['length_bin']:<12} {int(row['reads']):>10,} {int(row['nh1']):>10,} "
                    f"{row['nh1_frac']:>7.1%} {row['gx_frac']:>7.1%} {row['mm_frac']:>7.1%}\n")

        if not df_sec.empty:
            f.write("\nsecondary hit class breakdown by length bin\n")
            sec_summary = summarize_secondaries(df_sec)
            for _, row in sec_summary.iterrows():
                f.write(f"  {row['length_bin']:<12} {row['sec_class']:<22} "
                        f"{int(row['count']):>8,} / {int(row['total']):>8,} ({row['fraction']:.1%})\n")

    with open(summary_path) as f:
        print("\n" + f.read())
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    main()
