#!/usr/bin/env python3
"""
Fast GX-degree-by-length analysis for any morselle annotated BAM.

Works on coordinate-sorted BAMs with or without NH tags. Uses the GX tag
(semicolon-separated gene IDs) and RE tag (E/I/N) that are present in all
morselle annotate outputs. No name-sort required.

Output: TSV tables and a plain-text summary to stdout and --out directory.

Usage:
  analyze_gx_by_length.py --bam BAM [--bam BAM ...] --labels LABEL [LABEL ...] \
                           [--out DIR] [--bins BIN_EDGES] [--sample N]
"""

import argparse
import os
import sys
from collections import defaultdict

import pandas as pd
import pysam

DEFAULT_BINS   = [0, 250, 300, 350, 400, float("inf")]
DEFAULT_LABELS = ["≤250", "251-300", "301-350", "351-400", "≥401"]


def bin_label(n: int, bins=DEFAULT_BINS, labels=DEFAULT_LABELS) -> str:
    for lo, hi, lab in zip(bins, bins[1:], labels):
        if lo < n <= hi:
            return lab
    return labels[-1]


def clipped_bases(read: pysam.AlignedSegment) -> int:
    if read.cigartuples is None:
        return 0
    return sum(l for op, l in read.cigartuples if op in (4, 5))


def analyze_bam(bam_path: str, sample_every: int = 1) -> pd.DataFrame:
    """
    Scan primary records and return one row per (length_bin, gx_degree, re_class).
    sample_every=N samples 1-in-N reads by position in file for speed on huge BAMs.
    """
    rows = []
    n = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            n += 1
            if sample_every > 1 and n % sample_every != 0:
                continue

            rlen = read.query_length or 0
            clip  = clipped_bases(read)
            mapq  = read.mapping_quality

            gx_raw = read.get_tag("GX") if read.has_tag("GX") else ""
            if not gx_raw or gx_raw == "-":
                gx_degree = 0
            else:
                gx_degree = len(gx_raw.split(";"))

            nh = read.get_tag("NH") if read.has_tag("NH") else 1

            re_class = read.get_tag("RE") if read.has_tag("RE") else "?"

            rows.append({
                "length_bin":  bin_label(rlen),
                "read_length": rlen,
                "gx_degree":   min(gx_degree, 5),  # top-code at 5+
                "gx_degree_raw": gx_degree,
                "re_class":    re_class,
                "nh1":         int(nh == 1),
                "nh":          nh,
                "mapq":        mapq,
                "clip_frac":   clip / max(rlen, 1),
            })
    return pd.DataFrame(rows)


def length_summary(df: pd.DataFrame) -> pd.DataFrame:
    cat = pd.CategoricalDtype(categories=DEFAULT_LABELS, ordered=True)
    df = df.copy()
    df["length_bin"] = df["length_bin"].astype(cat)
    grp = df.groupby("length_bin", observed=True)
    s = grp.agg(
        reads          = ("gx_degree", "count"),
        gx0            = ("gx_degree", lambda x: (x == 0).sum()),
        gx1            = ("gx_degree", lambda x: (x == 1).sum()),
        gx2plus        = ("gx_degree", lambda x: (x >= 2).sum()),
        mean_length    = ("read_length", "mean"),
        mean_mapq      = ("mapq", "mean"),
        mean_clip_frac = ("clip_frac", "mean"),
        nh1_sum        = ("nh1", "sum"),
    ).reset_index()
    s["gx0_frac"]   = s["gx0"]   / s["reads"]
    s["gx1_frac"]   = s["gx1"]   / s["reads"]
    s["gx2p_frac"]  = s["gx2plus"] / s["reads"]
    s["nh1_frac"]   = s["nh1_sum"] / s["reads"]
    return s


def re_summary(df: pd.DataFrame) -> pd.DataFrame:
    """RE class (E/I/N) breakdown by length bin."""
    cat = pd.CategoricalDtype(categories=DEFAULT_LABELS, ordered=True)
    df = df.copy()
    df["length_bin"] = df["length_bin"].astype(cat)
    pivot = (
        df.groupby(["length_bin", "re_class"], observed=True)
        .size().rename("count").reset_index()
    )
    total = df.groupby("length_bin", observed=True).size().rename("total").reset_index()
    pivot = pivot.merge(total, on="length_bin")
    pivot["fraction"] = pivot["count"] / pivot["total"]
    return pivot


def gx_degree_dist(df: pd.DataFrame) -> pd.DataFrame:
    """GX degree distribution (0,1,2,3,4,5+) by length bin."""
    cat = pd.CategoricalDtype(categories=DEFAULT_LABELS, ordered=True)
    df = df.copy()
    df["length_bin"] = df["length_bin"].astype(cat)
    df["deg_bin"] = df["gx_degree"].apply(lambda x: str(min(x, 5)) if x < 5 else "5+")
    pivot = (
        df.groupby(["length_bin", "deg_bin"], observed=True)
        .size().rename("count").reset_index()
    )
    total = df.groupby("length_bin", observed=True).size().rename("total").reset_index()
    pivot = pivot.merge(total, on="length_bin")
    pivot["fraction"] = pivot["count"] / pivot["total"]
    return pivot


def print_summary(label: str, ls: pd.DataFrame, f=None):
    out = f or sys.stdout
    print(f"\n{'='*60}", file=out)
    print(f"  {label}", file=out)
    print(f"{'='*60}", file=out)
    print(f"{'bin':<12} {'reads':>12} {'GX=0%':>8} {'GX=1%':>8} {'GX≥2%':>8} "
          f"{'MAPQ':>7} {'clip%':>7} {'NH=1%':>8}", file=out)
    print("-" * 75, file=out)
    for _, r in ls.iterrows():
        print(f"{r['length_bin']:<12} {int(r['reads']):>12,} "
              f"{r['gx0_frac']:>7.1%} {r['gx1_frac']:>8.1%} {r['gx2p_frac']:>8.1%} "
              f"{r['mean_mapq']:>7.1f} {r['mean_clip_frac']:>7.1%} "
              f"{r['nh1_frac']:>8.1%}", file=out)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bam",    required=True, nargs="+", help="Annotated BAM(s)")
    ap.add_argument("--labels", required=True, nargs="+", help="Label for each BAM")
    ap.add_argument("--out",    default="gx_analysis", help="Output directory")
    ap.add_argument("--sample", type=int, default=10,
                    help="Sample 1-in-N primary reads for speed on large BAMs [10]")
    args = ap.parse_args()

    if len(args.bam) != len(args.labels):
        ap.error("--bam and --labels must have the same number of entries")

    os.makedirs(args.out, exist_ok=True)

    all_summaries = []
    summary_path  = os.path.join(args.out, "summary.txt")

    with open(summary_path, "w") as sf:
        for bam_path, label in zip(args.bam, args.labels):
            print(f"scanning {label}: {bam_path}  (1-in-{args.sample} sample)")
            df = analyze_bam(bam_path, sample_every=args.sample)
            print(f"  sampled {len(df):,} primary records")

            ls = length_summary(df)
            ls.insert(0, "label", label)
            all_summaries.append(ls)

            ls.to_csv(os.path.join(args.out, f"{label}_length_summary.tsv"),
                      sep="\t", index=False, float_format="%.4f")

            re = re_summary(df)
            re.insert(0, "label", label)
            re.to_csv(os.path.join(args.out, f"{label}_re_by_length.tsv"),
                      sep="\t", index=False, float_format="%.4f")

            gd = gx_degree_dist(df)
            gd.insert(0, "label", label)
            gd.to_csv(os.path.join(args.out, f"{label}_gx_degree.tsv"),
                      sep="\t", index=False, float_format="%.4f")

            print_summary(label, ls)
            print_summary(label, ls, f=sf)

    # Combined table across all BAMs
    if len(all_summaries) > 1:
        combined = pd.concat(all_summaries, ignore_index=True)
        combined.to_csv(os.path.join(args.out, "combined_length_summary.tsv"),
                        sep="\t", index=False, float_format="%.4f")
        print(f"\nwrote combined_length_summary.tsv")

    print(f"\nwrote summary to {summary_path}")


if __name__ == "__main__":
    main()
