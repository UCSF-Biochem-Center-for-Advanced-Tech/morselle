#!/usr/bin/env python3
"""Compare morselle molecule-counting modes against a reference matrix.

Reads one raw matrix per mode (MEX directories written by `morselle matrix`)
and a reference matrix (MEX directory or 10x .h5, e.g. the same library
sequenced on Illumina and counted independently), restricts all of them to
the shared barcodes and genes, and reports per mode how well it agrees with
the reference.

Why these metrics: SBX and the reference will differ in depth, so raw totals
are not comparable. What separates the modes is *expression-dependent* bias:
UMI-collision undercounting (fractional splitting, or same-gene molecules
merged on one UMI) hits highly expressed genes hardest, so after library-size
normalisation their log2(mode / reference) ratio falls as expression rises.
A mode with no collision bias has a flat ratio (slope ~ 0).

Output (tab-separated, also printed):
  summary.tsv   one row per mode: barcodes/genes compared, median per-barcode
                Pearson r of log1p CPM, median Spearman rho, median |log2 CPM
                ratio| for genes with >= --min-count reference counts, the
                slope of log2 CPM ratio vs log2 reference CPM, and the mode's
                total relative to the reference.
  per_gene.tsv  per mode and gene: summed counts, CPM, log2 ratio.

Requires numpy (and h5py only for an .h5 reference).
"""
import argparse
import gzip
import os
import sys

import numpy as np


def open_text(path):
    if os.path.exists(path):
        return open(path, "rt")
    if os.path.exists(path + ".gz"):
        return gzip.open(path + ".gz", "rt")
    raise FileNotFoundError(path)


def strip_suffix(barcode):
    return barcode.split("-")[0]


def strip_version(gene_id):
    return gene_id.split(".")[0] if gene_id.startswith("ENS") else gene_id


def read_mex(directory):
    with open_text(os.path.join(directory, "features.tsv")) as f:
        genes = [strip_version(line.split("\t")[0]) for line in f if line.strip()]
    with open_text(os.path.join(directory, "barcodes.tsv")) as f:
        barcodes = [strip_suffix(line.strip()) for line in f if line.strip()]
    counts = {}
    with open_text(os.path.join(directory, "matrix.mtx")) as f:
        header_seen = False
        for line in f:
            if line.startswith("%"):
                continue
            if not header_seen:
                header_seen = True
                continue
            g, b, v = line.split()
            key = (barcodes[int(b) - 1], genes[int(g) - 1])
            counts[key] = counts.get(key, 0.0) + float(v)
    return counts


def read_h5(path):
    import h5py  # only needed for .h5 references

    with h5py.File(path, "r") as f:
        m = f["matrix"]
        genes = [strip_version(x.decode()) for x in m["features"]["id"][:]]
        barcodes = [strip_suffix(x.decode()) for x in m["barcodes"][:]]
        data, indices, indptr = m["data"][:], m["indices"][:], m["indptr"][:]
    counts = {}
    for col, barcode in enumerate(barcodes):
        for pos in range(indptr[col], indptr[col + 1]):
            key = (barcode, genes[indices[pos]])
            counts[key] = counts.get(key, 0.0) + float(data[pos])
    return counts


def read_matrix(path):
    if os.path.isdir(path):
        inner = os.path.join(path, "raw_feature_bc_matrix")
        return read_mex(inner if os.path.isdir(inner) else path)
    return read_h5(path)


def dense(counts, barcodes, genes):
    bi = {b: i for i, b in enumerate(barcodes)}
    gi = {g: i for i, g in enumerate(genes)}
    out = np.zeros((len(barcodes), len(genes)))
    for (b, g), v in counts.items():
        if b in bi and g in gi:
            out[bi[b], gi[g]] += v
    return out


def rankdata(x):
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x))
    ranks[order] = np.arange(len(x))
    # average ties
    sorted_x = x[order]
    start = 0
    for i in range(1, len(x) + 1):
        if i == len(x) or sorted_x[i] != sorted_x[start]:
            ranks[order[start:i]] = (start + i - 1) / 2.0
            start = i
    return ranks


def corr(a, b):
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--reference", required=True, help="reference matrix: MEX dir, outs dir, or .h5")
    ap.add_argument("--mode", action="append", required=True, metavar="NAME=MATRIX",
                    help="mode name and its matrix (MEX/outs dir or .h5); repeat per mode")
    ap.add_argument("--barcodes", help="restrict to these barcodes (one per line, -1 suffix optional)")
    ap.add_argument("--min-count", type=float, default=10.0,
                    help="reference count floor for the |log2 ratio| and slope metrics (default 10)")
    ap.add_argument("--out-dir", default=".")
    args = ap.parse_args()

    reference = read_matrix(args.reference)
    modes = {}
    for spec in args.mode:
        name, _, path = spec.partition("=")
        if not path:
            sys.exit(f"--mode expects NAME=MATRIX, got {spec!r}")
        modes[name] = read_matrix(path)

    ref_barcodes = {b for b, _ in reference}
    barcodes = set.intersection(ref_barcodes, *({b for b, _ in m} for m in modes.values()))
    if args.barcodes:
        with open_text(args.barcodes) as f:
            wanted = {strip_suffix(line.strip()) for line in f if line.strip()}
        barcodes &= wanted
    genes = set.intersection({g for _, g in reference}, *({g for _, g in m} for m in modes.values()))
    barcodes, genes = sorted(barcodes), sorted(genes)
    if not barcodes or not genes:
        sys.exit("no shared barcodes/genes between the reference and every mode "
                 "(check barcode suffixes and gene id versions)")

    ref = dense(reference, barcodes, genes)
    ref_cpm = ref / ref.sum(axis=1, keepdims=True).clip(min=1) * 1e6
    os.makedirs(args.out_dir, exist_ok=True)
    rows = []
    with open(os.path.join(args.out_dir, "per_gene.tsv"), "w") as per_gene:
        per_gene.write("mode\tgene\tmode_count\treference_count\tmode_cpm\treference_cpm\tlog2_cpm_ratio\n")
        ref_gene = ref.sum(axis=0)
        ref_gene_cpm = ref_gene / ref_gene.sum() * 1e6
        for name, counts in modes.items():
            m = dense(counts, barcodes, genes)
            cpm = m / m.sum(axis=1, keepdims=True).clip(min=1) * 1e6
            pearson = [corr(np.log1p(cpm[i]), np.log1p(ref_cpm[i])) for i in range(len(barcodes))]
            spearman = [corr(rankdata(cpm[i]), rankdata(ref_cpm[i])) for i in range(len(barcodes))]
            mode_gene = m.sum(axis=0)
            mode_gene_cpm = mode_gene / max(mode_gene.sum(), 1) * 1e6
            keep = ref_gene >= args.min_count
            ratio = np.log2((mode_gene_cpm[keep] + 1) / (ref_gene_cpm[keep] + 1))
            level = np.log2(ref_gene_cpm[keep] + 1)
            slope = float(np.polyfit(level, ratio, 1)[0]) if keep.sum() >= 3 else float("nan")
            rows.append({
                "mode": name,
                "barcodes": len(barcodes),
                "genes": len(genes),
                "median_pearson_log1p_cpm": float(np.nanmedian(pearson)),
                "median_spearman_cpm": float(np.nanmedian(spearman)),
                "median_abs_log2_cpm_ratio": float(np.median(np.abs(ratio))) if keep.any() else float("nan"),
                "log2_ratio_vs_expression_slope": slope,
                "total_vs_reference": float(m.sum() / max(ref.sum(), 1)),
            })
            for j, gene in enumerate(genes):
                r = np.log2((mode_gene_cpm[j] + 1) / (ref_gene_cpm[j] + 1))
                per_gene.write(f"{name}\t{gene}\t{mode_gene[j]:.4f}\t{ref_gene[j]:.4f}\t"
                               f"{mode_gene_cpm[j]:.3f}\t{ref_gene_cpm[j]:.3f}\t{r:.4f}\n")

    columns = list(rows[0].keys())
    with open(os.path.join(args.out_dir, "summary.tsv"), "w") as out:
        out.write("\t".join(columns) + "\n")
        for row in rows:
            out.write("\t".join(f"{row[c]:.4f}" if isinstance(row[c], float) else str(row[c])
                                for c in columns) + "\n")
    print("\t".join(columns))
    for row in rows:
        print("\t".join(f"{row[c]:.4f}" if isinstance(row[c], float) else str(row[c]) for c in columns))


if __name__ == "__main__":
    main()
