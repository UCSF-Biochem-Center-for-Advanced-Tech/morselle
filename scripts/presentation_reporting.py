"""Pipeline metric discovery, insert reporting, and sample-specific slides."""

import json
import os
from pathlib import Path

# ── Metrics Ingestion & Export Helpers ────────────────────────────────────────

def load_json(path_str):
    if not path_str or not os.path.exists(path_str):
        return None
    with open(path_str, "r") as f:
        return json.load(f)

def export_insert_sizes(illumina_metrics, tsv_path):
    if not illumina_metrics or "insert_length_histogram" not in illumina_metrics:
        print(f"Warning: No insert_length_histogram found in metrics for {tsv_path}")
        return
    histo = illumina_metrics["insert_length_histogram"]
    with open(tsv_path, "w") as f:
        f.write("insert_length\tcount\n")
        for length in sorted(int(k) for k in histo.keys()):
            f.write(f"{length}\t{histo[str(length)]}\n")
    print(f"Exported insert-size distribution to {tsv_path} ({len(histo)} length bins)")

def slide_insert_size_distribution(prs, illumina_metrics):
    from presentation_layout import blank_slide, bg, header_bar, bullet_block, simple_table, box, Inches, UCSF_TEAL
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "Insert-size distribution", "Illumina PE overlap geometry & cDNA fragment recovery")

    if not illumina_metrics or not illumina_metrics.get("overlap_enabled"):
        bullet_block(s, [
            (0, "• PE overlap mode was not active or metrics JSON not supplied."),
            (0, "• Enable --overlap on sbxpress illumina to reconstruct PE fragment inserts."),
        ], font_size=16)
        return

    histo = illumina_metrics.get("insert_length_histogram", {})
    sins = illumina_metrics.get("overlap_short_sins", 0)
    mins = illumina_metrics.get("overlap_medium_mins", 0)
    lins = illumina_metrics.get("overlap_long_lins", 0)
    total_merged = sins + mins
    total = total_merged + lins

    headers = ["Category", "Insert class", "Read count", "Percentage"]
    rows = [
        ["Short insert", "sins (R1+R2 read-through)", f"{sins:,}", f"{100.0 * sins / max(1, total):.1f}%"],
        ["Medium insert", "mins (R1+R2 overlap)", f"{mins:,}", f"{100.0 * mins / max(1, total):.1f}%"],
        ["Long insert", "lins (unmerged mates)", f"{lins:,}", f"{100.0 * lins / max(1, total):.1f}%"],
        ["Total merged", "sins + mins", f"{total_merged:,}", f"{100.0 * total_merged / max(1, total):.1f}%"],
    ]
    simple_table(s, headers, rows, top=Inches(1.4), col_widths=[Inches(2.2), Inches(3.8), Inches(2.5), Inches(2.2)], font_size=13)

    if histo:
        lengths = [int(k) for k in histo.keys()]
        min_l, max_l = min(lengths), max(lengths)
        mean_l = sum(int(k) * v for k, v in histo.items()) / max(1, sum(histo.values()))
        box(s, Inches(0.4), Inches(3.8), Inches(12.5), Inches(0.4),
            f"Observed cDNA insert length span: {min_l}–{max_l} bp (mean {mean_l:.1f} bp across {total_merged:,} merged fragments)",
            font_size=14, bold=True, color=UCSF_TEAL)

    bullets = [
        (0, "Key Insights:"),
        (1, "• Overlap geometry removes adapters without sequence-matching false positives."),
        (1, "• Discordant base rate across overlap region: {:.4f}%.".format(
            100.0 * illumina_metrics.get("discordant_bases", 0) / max(1, illumina_metrics.get("merged_bases", 1))
        )),
        (1, "• Merged reads carry YC duplex-consensus tags for downstream tracking."),
    ]
    bullet_block(s, bullets, top=Inches(4.4), font_size=14)



def load_metric(explicit, report_dir, kind):
    """Find real pipeline filenames; ambiguous sample directories fail clearly."""
    if explicit:
        if not Path(explicit).is_file():
            raise FileNotFoundError(f"metrics not found: {explicit}")
        return load_json(explicit)
    if not report_dir:
        return None
    root = Path(report_dir)
    patterns = {
        "illumina": ["illumina_metrics.json", "*.illumina.metrics.json"],
        "molecules": ["molecules_metrics.json", "*.molecules.metrics.json", "*.bam.metrics.json"],
        "cell": ["cell_calling_metrics.json"],
    }
    paths = sorted({p for pattern in patterns[kind] for p in root.rglob(pattern) if p.is_file()})
    if kind == "illumina" and not paths:
        candidate = root.parent / "extract" / f"{root.name}.illumina.metrics.json"
        if candidate.is_file():
            paths = [candidate]
    if len(paths) > 1:
        raise ValueError(f"multiple {kind} metrics in {root}; provide an explicit metrics path")
    return load_json(paths[0]) if paths else None


def metric_rows(molecules, cells):
    """Values shown in the current-run slide, also usable without rendering."""
    rows = []
    for data, fields in [
        (molecules, [("Eligible reads", "eligible_records"), ("Counted molecules", "molecules_counted"),
                     ("Representatives", "representatives_tagged"), ("UMI nodes corrected", "umi_nodes_collapsed"),
                     ("Collision groups", "collision_groups")]),
        (cells, [("Called cells/barcodes", "filtered_cells"), ("Filtered molecules", "filtered_molecules")]),
    ]:
        if data is not None:
            rows.extend([label, f"{data[key]:,}"] for label, key in fields if key in data)
    return rows


def slide_run_metrics(prs, molecules, cells):
    from presentation_layout import blank_slide, bg, header_bar, bullet_block, simple_table, box, Inches
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "Current run metrics", "Values from the supplied pipeline outputs")
    rows = metric_rows(molecules, cells)
    if rows:
        simple_table(s, ["Metric", "Value"], rows,
                     col_widths=[Inches(7), Inches(5.5)], font_size=16)
    else:
        bullet_block(s, [(0, "No molecule or cell-calling metrics supplied.")])
    if molecules:
        mode = "fractional" if molecules.get("fractional_collisions") else "shared" if molecules.get("tagseq") else "winner-take-all"
        box(s, Inches(0.4), Inches(5.0), Inches(12.5), Inches(0.5), f"Collision policy: {mode}", font_size=16)
