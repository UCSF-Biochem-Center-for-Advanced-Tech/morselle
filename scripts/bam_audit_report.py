#!/usr/bin/env python3
"""Analyze the output of run_bam_audit.sh — run this locally on the returned files.

Reads the per-BAM JSON summaries and the hash-subsampled molecule dumps and prints:

  1. R2 CONTAMINATION AUDIT  — residual adapter / polyT / polyA / TSO, as % of mapped
     reads. Residual contamination cannot align, so STAR soft-clips it; that is the
     signal. Illumina is the trimming baseline.
  2. MAPPING & GENE ASSIGNMENT — uniqueness, region (exonic/intronic/intergenic),
     gene-call rate, multi-gene ambiguity, mismatches. Then the same by read-length
     bin, which isolates the long-read effect.
  3. MOLECULE CONCORDANCE — for the SAME molecule (CB+UMI) seen in two platforms, do
     they call the same gene? Uses Cell Ranger's counted representative (xf bit 8),
     one gene call per molecule.

     WARNING (learned the hard way): do NOT union GX across all reads of a molecule.
     Platforms differ in reads/molecule (Illumina ~3.7 vs SBX ~1.7), so the deeper
     platform accumulates more gene labels and you get a large phantom "disagreement".

Usage:  python3 bam_audit_report.py <bam_audit_out_dir>
"""
import gzip
import json
import math
import sys
from collections import Counter
from pathlib import Path

ORDER = [("illumina", "ILLUMINA"), ("p5", "SBX p5"), ("p7", "SBX p7")]


def load_json(d, tag):
    p = d / f"{tag}.json"
    return json.loads(p.read_text()) if p.exists() else None


REGION_RANK = {"E": 3, "N": 2, "I": 1, ".": 0, "": 0}


def load_molecules(d, tag):
    """(CB,UMI) -> dict of per-molecule mapping evidence (all reads aggregated):
       counted, gene (rep, single), n_reads, min_nh, max_mapq, region, max_readlen."""
    p = d / f"{tag}.molecules.tsv.gz"
    if not p.exists():
        return {}
    mols = {}
    with gzip.open(p, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            mols[(f[0], f[1])] = {
                "counted": f[2] == "1",
                "gene": frozenset(g for g in f[3].split(";") if g and g != "-"),
                "n_reads": int(f[4]),
                "min_nh": int(f[5]),
                "max_mapq": int(f[6]),
                "region": f[7],
                "max_readlen": int(f[8]),
            }
    return mols


def pct(x, d):
    return 100.0 * x / d if d else 0.0


def pearson(a, b):
    """Pearson correlation over the union of keys, treating absence as zero."""
    keys = a.keys() | b.keys()
    if len(keys) < 2:
        return float("nan")
    av = sum(a.get(k, 0) for k in keys) / len(keys)
    bv = sum(b.get(k, 0) for k in keys) / len(keys)
    cov = sum((a.get(k, 0) - av) * (b.get(k, 0) - bv) for k in keys)
    va = sum((a.get(k, 0) - av) ** 2 for k in keys)
    vb = sum((b.get(k, 0) - bv) ** 2 for k in keys)
    return cov / math.sqrt(va * vb) if va and vb else float("nan")


def contamination(js):
    print("\n" + "=" * 78)
    print("1. R2 CONTAMINATION AUDIT   (% of mapped reads; Illumina = trimming baseline)")
    print("=" * 78)
    hdr = f"{'':<26}" + "".join(f"{lab:>16}" for _, lab in ORDER if js.get(_))
    print(hdr)
    rows = [
        ("adapter (both ends)", lambda a: pct(a["clip5"]["adapter"] + a["clip3"]["adapter"], a["primary"] - a["unmapped"])),
        ("polyT (5' of R2)", lambda a: pct(a["clip5"]["poly_t"], a["primary"] - a["unmapped"])),
        ("polyA (3' of R2)", lambda a: pct(a["clip3"]["poly_a"], a["primary"] - a["unmapped"])),
        ("TSO (5' of R2)", lambda a: pct(a["clip5"]["tso"] + a["clip3"]["tso"], a["primary"] - a["unmapped"])),
        ("soft-clipped bases", lambda a: pct(a["softclip_bases"], a["read_bases"])),
        ("heavily clipped (>50%)", lambda a: pct(a["heavily_clipped"], a["primary"] - a["unmapped"])),
    ]
    for name, fn in rows:
        line = f"{name:<26}"
        for tag, _ in ORDER:
            a = js.get(tag)
            if a:
                line += f"{fn(a):>15.3f}%"
        print(line)
    print("\n  Residual adapter lands at R2's 5' end and residual polyT at its 3' end as")
    print("  polyA, because sbx10x emits R2 reverse-complemented. Both motifs are tested")
    print("  at both ends with indel-robust 8-mer matching.")


def mapping(js):
    print("\n" + "=" * 78)
    print("2. MAPPING & GENE ASSIGNMENT")
    print("=" * 78)
    print(f"{'':<26}" + "".join(f"{lab:>16}" for tag, lab in ORDER if js.get(tag)))
    rows = [
        ("mean read length", lambda a: a["read_bases"] / max(a["primary"] - a["unmapped"], 1), "{:>15.1f} "),
        ("unmapped", lambda a: pct(a["unmapped"], a["primary"]), "{:>15.2f}%"),
        ("unique (NH==1)", lambda a: pct(a["unique_nh1"], a["primary"] - a["unmapped"]), "{:>15.2f}%"),
        ("MAPQ 255", lambda a: pct(a["mapq255"], a["primary"] - a["unmapped"]), "{:>15.2f}%"),
        ("exonic", lambda a: pct(a["exonic"], a["primary"] - a["unmapped"]), "{:>15.2f}%"),
        ("intronic", lambda a: pct(a["intronic"], a["primary"] - a["unmapped"]), "{:>15.2f}%"),
        ("intergenic", lambda a: pct(a["intergenic"], a["primary"] - a["unmapped"]), "{:>15.2f}%"),
        ("gene called (GX)", lambda a: pct(a["gene_called"], a["primary"] - a["unmapped"]), "{:>15.2f}%"),
        ("gene collision |GX|>1", lambda a: pct(a["multi_gene"], a["gene_called"]), "{:>15.3f}%"),
        ("gene degree exactly 2", lambda a: pct(a.get("gene_degree", [0]*5)[2], a["gene_called"]), "{:>15.3f}%"),
        ("gene degree 3+", lambda a: pct(sum(a.get("gene_degree", [0]*5)[3:]), a["gene_called"]), "{:>15.3f}%"),
        ("transcript collision >1", lambda a: pct(sum(a.get("transcript_degree", [0]*5)[2:]), sum(a.get("transcript_degree", [0]*5)[1:])), "{:>15.3f}%"),
        ("mismatches /100bp", lambda a: pct(a["mismatches"], a["aligned_bases"]), "{:>15.3f} "),
        ("reads / molecule", lambda a: a["reads_per_molecule"], "{:>15.2f} "),
        ("UMIs / mapped read", lambda a: a["molecules"] / max(a["primary"] - a["unmapped"], 1), "{:>15.4f} "),
    ]
    for name, fn, fmt in rows:
        line = f"{name:<26}"
        for tag, _ in ORDER:
            a = js.get(tag)
            if a:
                line += fmt.format(fn(a))
        print(line)

    if any(a.get("terminal_collision_evaluated", 0) for a in js.values()):
        print("\n  --- terminal-exon preference on gene collisions (non-destructive) ---")
        print(f"  {'':<26}" + "".join(f"{lab:>16}" for tag, lab in ORDER if js.get(tag)))
        for name, key in [
            ("collisions evaluated", "terminal_collision_evaluated"),
            ("resolved to one gene", "terminal_collision_resolved"),
            ("still multiple genes", "terminal_collision_unresolved"),
            ("no terminal support", "terminal_collision_unsupported"),
        ]:
            print(f"  {name:<24}" + "".join(f"{js[tag].get(key, 0):>16,}" for tag, _ in ORDER if js.get(tag)))

    for tag, lab in ORDER:
        a = js.get(tag)
        if not a:
            continue
        print(f"\n  --- {lab}: by read length (isolates the long-read effect) ---")
        print(f"  {'len':>11} {'reads':>13} {'uniq%':>8} {'exonic%':>9} {'geneGX%':>9} {'multi%':>8} {'nM/100bp':>9}")
        for b in a["bins"]:
            n = b["reads"]
            if n < 500:
                continue
            print(f"  {b['len_lo']:>4}-{b['len_hi']:<6} {n:>13,} "
                  f"{pct(b['unique_nh1'], n):>7.2f}% {pct(b['exonic'], n):>8.2f}% "
                  f"{pct(b['gene_called'], n):>8.2f}% {pct(b['multi_gene'], n):>7.3f}% "
                  f"{pct(b['mismatches'], b['aligned_bases']):>8.3f}")


def concordance(mols, sample_rate):
    print("\n" + "=" * 78)
    print("3. MOLECULE CONCORDANCE   (same CB+UMI, counted representative: same gene?)")
    print("=" * 78)
    pairs = [("illumina", "p5"), ("illumina", "p7"), ("p5", "p7")]
    names = dict(ORDER)
    for x, y in pairs:
        a, b = mols.get(x), mols.get(y)
        if not a or not b:
            continue
        shared = [k for k in a.keys() & b.keys() if a[k]["counted"] and b[k]["counted"]]
        ac = {k for k, value in a.items() if value["counted"]}
        bc = {k for k, value in b.items() if value["counted"]}
        n = len(shared)
        agree = sum(1 for k in shared if a[k]["gene"] == b[k]["gene"])
        print(f"\n  {names[x]} vs {names[y]}{'  (long-vs-long control)' if x == 'p5' else ''}")
        print(f"    counted molecule keys (sampled; ~estimated at 1/{sample_rate}):")
        print(f"      shared       {len(ac & bc):>9,}  ~{len(ac & bc) * sample_rate:>12,}")
        print(f"      {names[x]:<10} only {len(ac - bc):>9,}  ~{len(ac - bc) * sample_rate:>12,}")
        print(f"      {names[y]:<10} only {len(bc - ac):>9,}  ~{len(bc - ac) * sample_rate:>12,}")
        print(f"    both counted a molecule with a gene: {n:,}")
        if n:
            print(f"    agree            {agree:>9,}   {100.0*agree/n:6.2f}%")
            print(f"    different gene   {n-agree:>9,}   {100.0*(n-agree)/n:6.2f}%")
        cells_a = Counter(k[0] for k in ac)
        cells_b = Counter(k[0] for k in bc)
        genes_a = Counter(g for k in ac for g in a[k]["gene"])
        genes_b = Counter(g for k in bc for g in b[k]["gene"])
        print(f"    per-cell molecule Pearson r: {pearson(cells_a, cells_b):.4f}")
        print(f"    per-gene molecule Pearson r: {pearson(genes_a, genes_b):.4f}")


def advantage(mols):
    """Section 4: WHERE do longer reads map better? For molecules shared by two
    platforms, adjudicate uniqueness / region / countability. The reference (short)
    platform is Illumina; 'LONG' is the SBX run. A positive asymmetry (LONG-only >>
    SHORT-only) is the long-read advantage."""
    print("\n" + "=" * 78)
    print("4. LONG-READ ADVANTAGE   (same molecule, which platform mapped it better?)")
    print("=" * 78)
    names = dict(ORDER)
    for short, long_ in [("illumina", "p5"), ("illumina", "p7")]:
        s, l = mols.get(short), mols.get(long_)
        if not s or not l:
            continue
        shared = list(s.keys() & l.keys())
        n = len(shared)
        # uniqueness: min_nh == 1 means placed uniquely somewhere
        l_uniq_s_multi = sum(1 for k in shared if l[k]["min_nh"] == 1 and s[k]["min_nh"] > 1)
        s_uniq_l_multi = sum(1 for k in shared if s[k]["min_nh"] == 1 and l[k]["min_nh"] > 1)
        # in a gene (exonic OR intronic) vs intergenic. NOT exonic-vs-intronic: a long
        # read spanning an exon-intron junction is intronic yet still the same gene, so
        # ranking exonic above intronic would wrongly penalize length.
        ingene = lambda r: r in ("E", "N")
        l_ingene_s_inter = sum(1 for k in shared if ingene(l[k]["region"]) and not ingene(s[k]["region"]))
        s_ingene_l_inter = sum(1 for k in shared if ingene(s[k]["region"]) and not ingene(l[k]["region"]))
        # confidently counted (got a gene) on one platform only
        l_only_counted = sum(1 for k in shared if l[k]["counted"] and not s[k]["counted"])
        s_only_counted = sum(1 for k in shared if s[k]["counted"] and not l[k]["counted"])
        # among molecules BOTH counted but calling different genes, who has better evidence
        disc = [k for k in shared if l[k]["counted"] and s[k]["counted"] and l[k]["gene"] != s[k]["gene"]]
        l_disc_better = sum(1 for k in disc if (REGION_RANK[l[k]["region"]], -l[k]["min_nh"]) > (REGION_RANK[s[k]["region"]], -s[k]["min_nh"]))
        s_disc_better = sum(1 for k in disc if (REGION_RANK[s[k]["region"]], -s[k]["min_nh"]) > (REGION_RANK[l[k]["region"]], -l[k]["min_nh"]))

        SH, LO = names[short], names[long_]
        print(f"\n  {LO} (long) vs {SH} (short) — {n:,} shared molecules")
        print(f"  {'axis':<24}{'LONG better':>14}{'SHORT better':>14}{'asymmetry':>12}")
        def row(axis, lo, sh):
            asy = f"{lo/max(sh,1):.1f}x" if sh else ("inf" if lo else "-")
            print(f"  {axis:<24}{lo:>9,} {100.0*lo/max(n,1):>4.1f}%{sh:>9,} {100.0*sh/max(n,1):>4.1f}%{asy:>12}")
        row("uniquely mapped", l_uniq_s_multi, s_uniq_l_multi)
        row("in a gene (vs interg.)", l_ingene_s_inter, s_ingene_l_inter)
        row("counted (got a gene)", l_only_counted, s_only_counted)
        if disc:
            print(f"  among {len(disc):,} both-counted DISAGREEMENTS: "
                  f"{LO} stronger evidence {l_disc_better:,} vs {SH} {s_disc_better:,}")
    print("\n  'uniquely mapped' = one platform placed the molecule at a single locus (NH==1)")
    print("  while the other multi-mapped it (this is conservative for SBX, which has")
    print("  FEWER reads/molecule and thus fewer chances at a unique hit, yet still wins).")
    print("  'in a gene' = exonic or intronic vs intergenic. 'counted' = Cell Ranger")
    print("  confidently assigned a gene on one platform only. A large LONG-better")
    print("  asymmetry is the long-read mapping advantage.")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    d = Path(sys.argv[1])
    js = {tag: load_json(d, tag) for tag, _ in ORDER}
    js = {k: v for k, v in js.items() if v}
    if not js:
        print(f"no *.json found in {d}")
        sys.exit(1)
    samples = {v["molecule_sample"] for v in js.values()}
    if len(samples) > 1:
        print(f"WARNING: molecule subsample rates differ {samples} — dumps are NOT joinable.")
    mols = {tag: load_molecules(d, tag) for tag, _ in ORDER}
    contamination(js)
    mapping(js)
    concordance(mols, next(iter(samples)))
    advantage(mols)


if __name__ == "__main__":
    main()
