#!/usr/bin/env python3
"""
Build sbxpress mapping-analysis presentation as a .pptx file.
Covers slides from docs/archive/PRESENTATION_PLAN.md that are marked presentable now.

Usage: python3 scripts/build_presentation.py [--out FILE]
"""
import argparse
from presentation_layout import *
from presentation_history import *
from presentation_reporting import *

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="Build sbxpress presentation slide deck from pipeline metrics.")
    ap.add_argument("--out", default="sbxpress_mapping_analysis.pptx", help="Output PPTX presentation file.")
    ap.add_argument("--report-dir", help="Directory containing pipeline outputs and metrics JSON files.")
    ap.add_argument("--illumina-metrics", help="Path to sbxpress illumina metrics JSON.")
    ap.add_argument("--molecules-metrics", help="Path to sbxpress molecules metrics JSON.")
    ap.add_argument("--cell-metrics", help="Path to sbxpress cells metrics JSON.")
    ap.add_argument("--export-insert-sizes", help="Export insert-size distribution histogram as TSV to this path.")
    args = ap.parse_args()

    illumina_metrics = load_metric(args.illumina_metrics, args.report_dir, "illumina")
    molecules_metrics = load_metric(args.molecules_metrics, args.report_dir, "molecules")
    cell_metrics = load_metric(args.cell_metrics, args.report_dir, "cell")

    if args.export_insert_sizes:
        export_insert_sizes(illumina_metrics, args.export_insert_sizes)

    prs = new_prs()

    slide_title(prs)
    slide_run_metrics(prs, molecules_metrics, cell_metrics)
    historical_start = len(prs.slides)
    slide_core_message(prs)
    slide_pipeline(prs)
    slide_insert_size_distribution(prs, illumina_metrics)
    slide_m1_paradox(prs)
    slide_m2_orientation(prs)
    slide_m3_secondary(prs)
    slide_m4_churn(prs)
    slide_validation_status(prs)
    slide_next_steps(prs)

    # The remaining mapping slides are historical examples, never inferred
    # from this sample's metrics. Label every one explicitly.
    for i, slide in enumerate(prs.slides):
        if i < historical_start or i == historical_start + 2:
            continue
        box(slide, Inches(9.8), Inches(7.1), Inches(3.2), Inches(0.25),
            "Historical mapping example", font_size=9, color=RED_WARN)
    prs.save(args.out)
    print(f"saved {args.out}  ({len(prs.slides)} slides)")


if __name__ == "__main__":
    main()

