#!/usr/bin/env python3
"""
Compare downstream outputs of the baseline and multimap evaluation arms (mapping.md §4).

Reads metrics JSONs and matrix outputs from both arms and the three existing
CellRanger/sbxpress comparators (Illumina, P5, P7 M2_GEX) and prints a
side-by-side concordance table.

Usage:
  compare_downstream.py --eval-dir DIR [--comparators DIR [DIR ...]] [--out FILE]

  --eval-dir DIR    mapping_evaluation root (contains baseline/ and multimap/)
  --comparators DIR existing cellranger/sbxpress outputs to include as reference
                    (each must have outs/filtered_feature_bc_matrix.h5 and a
                     *molecules.metrics.json or equivalent)
  --out FILE        output TSV [default: stdout]
"""

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd


def load_metrics(metrics_path: str) -> dict:
    if not os.path.exists(metrics_path):
        return {}
    with open(metrics_path) as f:
        return json.load(f)


def load_cellranger_summary(outs_dir: str) -> dict:
    """Read CellRanger metrics_summary.csv if present."""
    csv_path = os.path.join(outs_dir, "metrics_summary.csv")
    if not os.path.exists(csv_path):
        return {}
    try:
        df = pd.read_csv(csv_path)
        # CellRanger uses comma-formatted numbers and % signs
        result = {}
        for col in df.columns:
            val = str(df[col].iloc[0]).replace(",", "").replace("%", "").strip()
            try:
                result[col] = float(val)
            except ValueError:
                result[col] = df[col].iloc[0]
        return result
    except Exception:
        return {}


def h5_cell_count(h5_path: str) -> int | None:
    if not os.path.exists(h5_path):
        return None
    try:
        import scipy.io
        import h5py
        with h5py.File(h5_path, "r") as f:
            barcodes = f["matrix/barcodes"][:]
        return len(barcodes)
    except Exception:
        pass
    # fallback: count barcodes.tsv.gz next to the h5
    mex_dir = os.path.join(os.path.dirname(h5_path), "filtered_feature_bc_matrix")
    bc_path = os.path.join(mex_dir, "barcodes.tsv.gz")
    if os.path.exists(bc_path):
        import gzip
        with gzip.open(bc_path, "rt") as fh:
            return sum(1 for _ in fh)
    return None


def h5_total_counts(h5_path: str) -> int | None:
    if not os.path.exists(h5_path):
        return None
    # Try MEX directory (always present alongside h5)
    mex_dir = os.path.join(os.path.dirname(h5_path), "filtered_feature_bc_matrix")
    mtx_path = os.path.join(mex_dir, "matrix.mtx.gz")
    if os.path.exists(mtx_path):
        import gzip
        total = 0
        with gzip.open(mtx_path, "rt") as fh:
            for line in fh:
                if line.startswith("%"):
                    continue
                parts = line.split()
                if len(parts) == 3:
                    try:
                        total += int(parts[2])
                    except ValueError:
                        pass
        return total
    return None


def describe_arm(label: str, arm_dir: str) -> dict:
    """Collect metrics for one evaluation arm."""
    row = {"label": label}

    ann_metrics = load_metrics(os.path.join(arm_dir, "M2_GEX.annotation.metrics.json"))
    mol_metrics = load_metrics(os.path.join(arm_dir, "M2_GEX.molecules.metrics.json"))

    # Annotation metrics (actual keys from sbxpress annotate JSON)
    row["primary_reads"]  = ann_metrics.get("primary")
    row["gene_assigned"]  = ann_metrics.get("gene_assigned")
    row["gx_degree_1"]    = ann_metrics.get("gene_degree_1")
    row["gx_degree_2plus"] = (
        (ann_metrics.get("gene_degree_2") or 0) +
        (ann_metrics.get("gene_degree_3plus") or 0)
    ) or None

    # Molecule metrics (actual keys from sbxpress molecules JSON)
    row["molecules"]      = mol_metrics.get("molecules_counted")
    row["eligible_reads"] = mol_metrics.get("eligible_records")

    filt_h5 = os.path.join(arm_dir, "outs", "filtered_feature_bc_matrix.h5")
    row["cells"]          = h5_cell_count(filt_h5)
    row["filtered_umis"]  = h5_total_counts(filt_h5)

    return row


def describe_comparator(label: str, cellranger_outs: str) -> dict:
    """Collect metrics for an existing CellRanger or sbxpress comparator output."""
    row: dict = {"label": label}
    # CellRanger metrics_summary.csv column names vary by version; read what we can
    # and fall back to counting from MEX files.
    cr = load_cellranger_summary(cellranger_outs)
    row["primary_reads"]  = None  # not in CR summary
    row["gene_assigned"]  = None  # not in CR summary
    row["gx_degree_1"]    = None  # not in CR summary
    row["gx_degree_2plus"] = None
    row["molecules"]      = cr.get("Median UMI Counts per Cell",
                            cr.get("Median UMI counts per cell", None))
    row["eligible_reads"] = None
    row["cells"]          = cr.get("Estimated Number of Cells",
                            cr.get("Cells", None))

    filt_h5 = os.path.join(cellranger_outs, "filtered_feature_bc_matrix.h5")
    if row["cells"] is None:
        row["cells"] = h5_cell_count(filt_h5)
    row["filtered_umis"]  = h5_total_counts(filt_h5)

    return row


def fmt(val) -> str:
    if val is None:
        return "—"
    if isinstance(val, float):
        return f"{val:,.1f}"
    if isinstance(val, int):
        return f"{val:,}"
    return str(val)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--eval-dir", required=True,
                    help="mapping_evaluation root directory")
    ap.add_argument("--comparators", nargs="*", default=None,
                    help="Paths to CellRanger/sbxpress outs/ directories to include")
    ap.add_argument("--out", default=None, help="Output TSV file [stdout]")
    args = ap.parse_args()

    rows = []

    # Evaluation arms
    baseline_dir = os.path.join(args.eval_dir, "baseline")
    multimap_dir = os.path.join(args.eval_dir, "multimap")
    if os.path.isdir(baseline_dir):
        rows.append(describe_arm("sbxpress-baseline", baseline_dir))
    if os.path.isdir(multimap_dir):
        rows.append(describe_arm("sbxpress-multimap", multimap_dir))

    # Existing comparators
    default_comparators = [
        ("illumina-CR",
         "/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/illumina_data/cellranger_outputs/M2_GEX/outs"),
        ("sbxpress-p5",
         "/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/p5_out_7-7/cellranger_outputs/M2_GEX/outs"),
        ("sbxpress-p7",
         "/mnt/analysis/active_projects/roche_SBX_ABRF/CAT02/p7_out_7-7/cellranger_outputs/M2_GEX/outs"),
    ]
    comparator_paths = [(os.path.basename(p), p) for p in (args.comparators or [])]
    if not comparator_paths:
        comparator_paths = default_comparators

    for label, outs_dir in comparator_paths:
        if os.path.isdir(outs_dir):
            rows.append(describe_comparator(label, outs_dir))
        else:
            print(f"  skipping {label}: {outs_dir} not found", file=sys.stderr)

    if not rows:
        print("no data found", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(rows)
    cols_ordered = ["label", "primary_reads", "gene_assigned", "gx_degree_1",
                    "gx_degree_2plus", "molecules", "eligible_reads", "cells", "filtered_umis"]
    cols_present = [c for c in cols_ordered if c in df.columns]
    df = df[cols_present]

    # Print formatted table to stdout always
    print(f"\n{'label':<22} " + "  ".join(f"{c:<14}" for c in cols_present[1:]))
    print("-" * (22 + 16 * len(cols_present)))
    for _, row in df.iterrows():
        vals = "  ".join(f"{fmt(row.get(c)):>14}" for c in cols_present[1:])
        print(f"{row['label']:<22} {vals}")

    if args.out:
        df.to_csv(args.out, sep="\t", index=False)
        print(f"\nwrote {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
