#!/usr/bin/env python3
"""
Build sbxpress mapping-analysis presentation as a .pptx file.
Covers slides from PRESENTATION_PLAN.md that are marked presentable now.

Usage: python3 scripts/build_presentation.py [--out FILE]
"""
import argparse
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Inches, Pt

# ── Palette ───────────────────────────────────────────────────────────────────
UCSF_NAVY   = RGBColor(0x05, 0x2F, 0x5F)
UCSF_TEAL   = RGBColor(0x00, 0x82, 0x8B)
UCSF_GOLD   = RGBColor(0xB3, 0x9B, 0x13)
WHITE       = RGBColor(0xFF, 0xFF, 0xFF)
LIGHT_GREY  = RGBColor(0xF2, 0xF2, 0xF2)
MID_GREY    = RGBColor(0xCC, 0xCC, 0xCC)
DARK_GREY   = RGBColor(0x44, 0x44, 0x44)
RED_WARN    = RGBColor(0xC0, 0x39, 0x2B)

W = Inches(13.33)   # widescreen 16:9
H = Inches(7.5)

# ── Helpers ───────────────────────────────────────────────────────────────────

def new_prs() -> Presentation:
    prs = Presentation()
    prs.slide_width  = W
    prs.slide_height = H
    return prs

def blank_slide(prs):
    layout = prs.slide_layouts[6]   # completely blank
    return prs.slides.add_slide(layout)

def bg(slide, color=WHITE):
    fill = slide.background.fill
    fill.solid()
    fill.fore_color.rgb = color

def box(slide, left, top, width, height, text="", font_size=18,
        bold=False, color=DARK_GREY, bg_color=None, align=PP_ALIGN.LEFT,
        wrap=True):
    txb = slide.shapes.add_textbox(left, top, width, height)
    if bg_color:
        txb.fill.solid()
        txb.fill.fore_color.rgb = bg_color
    tf = txb.text_frame
    tf.word_wrap = wrap
    p = tf.paragraphs[0]
    p.alignment = align
    run = p.add_run()
    run.text = text
    run.font.size  = Pt(font_size)
    run.font.bold  = bold
    run.font.color.rgb = color
    return txb

def header_bar(slide, title, subtitle=""):
    # Navy bar across top
    bar = slide.shapes.add_shape(
        1,  # MSO_SHAPE_TYPE.RECTANGLE
        Inches(0), Inches(0), W, Inches(1.1)
    )
    bar.fill.solid()
    bar.fill.fore_color.rgb = UCSF_NAVY
    bar.line.fill.background()

    box(slide, Inches(0.3), Inches(0.08), Inches(12), Inches(0.55),
        text=title, font_size=28, bold=True, color=WHITE, align=PP_ALIGN.LEFT)
    if subtitle:
        box(slide, Inches(0.3), Inches(0.62), Inches(12), Inches(0.38),
            text=subtitle, font_size=14, color=UCSF_GOLD, align=PP_ALIGN.LEFT)

def framing_box(slide, text, top=Inches(6.9)):
    """Small italic caveat box at the bottom."""
    b = slide.shapes.add_shape(1, Inches(0.3), top, Inches(12.7), Inches(0.5))
    b.fill.solid()
    b.fill.fore_color.rgb = LIGHT_GREY
    b.line.color.rgb = MID_GREY
    tf = b.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.alignment = PP_ALIGN.LEFT
    run = p.add_run()
    run.text = "⚠ " + text
    run.font.size = Pt(9)
    run.font.italic = True
    run.font.color.rgb = DARK_GREY

def simple_table(slide, headers, rows,
                 left=Inches(0.4), top=Inches(1.3),
                 col_widths=None, font_size=14):
    n_cols = len(headers)
    n_rows = len(rows)
    if col_widths is None:
        col_widths = [Inches(12.5 / n_cols)] * n_cols
    row_h = Inches(0.38)
    tbl = slide.shapes.add_table(
        n_rows + 1, n_cols,
        left, top,
        sum(col_widths), row_h * (n_rows + 1)
    ).table
    tbl.columns[0].width  # touch to init

    for ci, w in enumerate(col_widths):
        tbl.columns[ci].width = w

    def cell_fmt(cell, text, hdr=False):
        cell.text = text
        p = cell.text_frame.paragraphs[0]
        p.alignment = PP_ALIGN.CENTER
        run = p.runs[0] if p.runs else p.add_run()
        run.font.size = Pt(font_size)
        run.font.bold = hdr
        run.font.color.rgb = WHITE if hdr else DARK_GREY
        if hdr:
            cell.fill.solid()
            cell.fill.fore_color.rgb = UCSF_NAVY
        else:
            cell.fill.solid()
            cell.fill.fore_color.rgb = WHITE

    for ci, h in enumerate(headers):
        cell_fmt(tbl.cell(0, ci), h, hdr=True)
    for ri, row in enumerate(rows):
        for ci, val in enumerate(row):
            cell_fmt(tbl.cell(ri + 1, ci), str(val))
    return tbl

def bullet_block(slide, items, left=Inches(0.5), top=Inches(1.3),
                 width=Inches(12.3), font_size=16, color=DARK_GREY):
    txb = slide.shapes.add_textbox(left, top, width, Inches(5.5))
    tf = txb.text_frame
    tf.word_wrap = True
    first = True
    for item in items:
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.space_before = Pt(4)
        indent, text = (item[0], item[1]) if isinstance(item, tuple) else (0, item)
        p.level = indent
        run = p.add_run()
        run.text = text
        run.font.size = Pt(font_size - indent * 2)
        run.font.color.rgb = color
        if text.startswith("▶") or text.startswith("•"):
            run.font.bold = False

# ── Slide builders ────────────────────────────────────────────────────────────

def slide_title(prs):
    s = blank_slide(prs)
    bg(s, UCSF_NAVY)
    # Main title
    box(s, Inches(0.6), Inches(1.8), Inches(12), Inches(1.4),
        "sbxpress: transparent minimap2-based\nGEX processing for Roche SBX reads",
        font_size=36, bold=True, color=WHITE, align=PP_ALIGN.LEFT)
    box(s, Inches(0.6), Inches(3.4), Inches(10), Inches(0.6),
        "Mapping ambiguity, multimap behaviour, and P5/P7 orientation — preliminary results",
        font_size=18, color=UCSF_GOLD, align=PP_ALIGN.LEFT)
    box(s, Inches(0.6), Inches(6.5), Inches(10), Inches(0.4),
        "CAT02 M2 GEX  ·  2026-07-30  ·  In progress — full-scale P5/P7 run pending",
        font_size=12, color=MID_GREY, align=PP_ALIGN.LEFT)


def slide_core_message(prs):
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "Core message", "What sbxpress is and isn't")

    bullets = [
        (0, "sbxpress replaces the opaque short-read mapping/assignment path with explicit,"),
        (0, "inspectable minimap2 → annotation → molecule → matrix → cell stages."),
        (0, ""),
        (0, "Two objectives in strict order:"),
        (1, "① Trusted 10x-style baseline  — reproduce relevant decisions; validate against CellRanger."),
        (1, "② SBX-specific optimisation  — only after baseline behaviour is validated."),
        (0, ""),
        (0, "What sbxpress does not claim yet:"),
        (1, "• Drop-in replacement for every reference workflow decision."),
        (1, "• Full validation of molecule calling, cell calling, or matrix construction."),
        (1, "• Demonstrated improvement from longer SBX cDNA reads (that experiment is in progress)."),
    ]
    bullet_block(s, bullets, font_size=17)


def slide_pipeline(prs):
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "Pipeline overview", "Two input paths, one shared tagged-cDNA handoff")

    box(s, Inches(0.4), Inches(1.25), Inches(12.5), Inches(0.5),
        "Raw SBX reads  →  sbxpress extract  ┐",
        font_size=15, color=DARK_GREY)
    box(s, Inches(0.4), Inches(1.65), Inches(12.5), Inches(0.5),
        "Illumina R1/R2  →  sbxpress illumina  ┘  →  tagged cDNA FASTQ  →  minimap2  →  annotation  →  molecules  →  matrix  →  cells",
        font_size=15, color=DARK_GREY)

    # Two mapping modes
    headers = ["Mode", "minimap2 flag", "Annotation flag", "Valid for"]
    rows = [
        ["Production baseline", "--secondary=no", "--unordered", "Counts, molecules, matrix, cell calling"],
        ["Multimap diagnostic", "-N 10", "--multimap-tags", "NH/HI/mm accounting; secondary-hit classification"],
    ]
    simple_table(s, headers, rows,
                 top=Inches(2.6),
                 col_widths=[Inches(2.4), Inches(2.2), Inches(2.2), Inches(5.5)],
                 font_size=13)

    box(s, Inches(0.4), Inches(4.0), Inches(12.5), Inches(0.4),
        "NH=1 means one reported alignment under the chosen mapper policy — not one biologically unique gene.",
        font_size=13, bold=True, color=UCSF_TEAL)

    box(s, Inches(0.4), Inches(4.5), Inches(12.5), Inches(0.5),
        "Tags written by sbxpress annotate:", font_size=13, bold=True, color=DARK_GREY)
    box(s, Inches(0.4), Inches(4.9), Inches(12.5), Inches(1.4),
        "GX  compatible gene IDs  ·  GN  gene names  ·  TX  transcript IDs  "
        "·  RE  region class (E/N/I = exonic/intronic/intergenic)  "
        "·  NH / HI  reported hit count / index  ·  mm  single-gene multimap rescue flag",
        font_size=13, color=DARK_GREY)


def slide_m1_paradox(prs):
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "M1 — The central paradox",
               "Full-scale P5/P7 production BAMs  ·  ~61M sampled primary records  ·  1-in-10 systematic sample")

    box(s, Inches(0.4), Inches(1.2), Inches(12.5), Inches(0.4),
        "Longer reads align more confidently (MAPQ ↑) but land in annotated gene territory less often (GX=1% ↓).",
        font_size=15, bold=True, color=UCSF_TEAL)

    headers = ["Length bin", "P5  GX=1%", "P5  MAPQ", "P5  clip%", "P7  GX=1%", "P7  MAPQ", "P7  clip%"]
    rows = [
        ["≤250 bp",  "60.0%", "52.1", "3.1%",  "62.3%", "53.3", "3.3%"],
        ["251–300 bp","55.5%", "54.1", "4.0%",  "57.6%", "54.7", "3.4%"],
        ["301–350 bp","53.4%", "55.0", "5.2%",  "55.7%", "55.6", "3.9%"],
        ["351–400 bp","52.5%", "55.4", "7.4%",  "55.3%", "56.4", "4.3%"],
        ["≥401 bp",  "51.6%", "54.8", "15.6%", "55.7%", "56.8", "5.2%"],
    ]
    simple_table(s, headers, rows,
                 top=Inches(1.75),
                 col_widths=[Inches(1.6), Inches(1.6), Inches(1.5), Inches(1.5),
                              Inches(1.6), Inches(1.5), Inches(1.5)],
                 font_size=13)

    box(s, Inches(0.4), Inches(4.55), Inches(12.5), Inches(0.5),
        "P5 clip% diverges sharply at ≥401 bp (15.6% vs 5.2% for P7) — see slide M2.",
        font_size=14, color=UCSF_GOLD, bold=True)

    box(s, Inches(0.4), Inches(5.1), Inches(12.5), Inches(0.8),
        "GX=0 (no gene assigned) rises 39%→48% P5, 37%→44% P7 across the length range. "
        "Multi-gene (GX≥2) stays flat at ~0.6% — the lost single-gene reads go to "
        "unannotated space, not to a competing gene.",
        font_size=13, color=DARK_GREY)

    framing_box(s, "Baseline (--secondary=no) BAMs; NH=1 is trivially 100% here and is not shown. "
                "GX=1% decline is a separate signal from the NH=1 question and is not yet decomposed "
                "by secondary-hit class at this scale.")


def slide_m2_orientation(prs):
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "M2 — P5 vs P7 orientation effect at long lengths",
               "Same data source as M1")

    box(s, Inches(0.4), Inches(1.2), Inches(12.5), Inches(0.5),
        "The orientation effect is invisible in short reads and only surfaces in the longest-read tail.",
        font_size=15, bold=True, color=UCSF_TEAL)

    headers = ["Length bin", "P5 soft-clip%", "P7 soft-clip%", "Δ (P5−P7)"]
    rows = [
        ["≤250 bp",   "3.1%",  "3.3%", "−0.2 pp"],
        ["251–300 bp","4.0%",  "3.4%", "+0.6 pp"],
        ["301–350 bp","5.2%",  "3.9%", "+1.3 pp"],
        ["351–400 bp","7.4%",  "4.3%", "+3.1 pp"],
        ["≥401 bp",   "15.6%", "5.2%", "+10.4 pp ◀"],
    ]
    simple_table(s, headers, rows,
                 top=Inches(1.8),
                 col_widths=[Inches(2.2), Inches(2.8), Inches(2.8), Inches(3.5)],
                 font_size=14)

    bullets = [
        (0, "Working hypothesis: SBX signal quality degrades over read position."),
        (0, "  P5 reads CB+UMI first → cDNA is decoded second → longer cDNA reads accumulate more degradation → excess clipping."),
        (0, "  P7 reads cDNA first → CB/UMI degrades instead → lower parse rate, but surviving cDNA reads are cleaner."),
        (0, ""),
        (0, "Same mechanism explains P7's better barcode accuracy (76.7% vs 68.4% whitelist hit rate in the cross-platform audit)."),
        (0, "The length-dependent clip divergence is new evidence — not previously characterised at full scale."),
    ]
    bullet_block(s, bullets, top=Inches(4.05), font_size=14)

    framing_box(s, "Confirm before this slide ships: P5/P7 here are the same molecule in two orientations, "
                "not 10x 5'/3' capture-chemistry variants (different exon coverage would create a similar pattern for unrelated reasons).")


def slide_m3_secondary(prs):
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "M3 — Where secondary hits actually go",
               "M2 compact cohort  ·  sbxpress illumina path  ·  ~9.4M reads  ·  250–650 bp range")

    box(s, Inches(0.4), Inches(1.2), Inches(12.5), Inches(0.45),
        "When a longer read picks up a secondary hit, it almost never lands on a competing annotated gene.",
        font_size=15, bold=True, color=UCSF_TEAL)

    headers = ["Secondary class", "Typical fraction", "What it means"]
    rows = [
        ["unannotated_mapped", "55–58%", "Genomic hit, no annotated gene — pseudogene, repeat, annotation gap"],
        ["decoy_other",        "19–27%", "Random/decoy chromosome — AS gap ≈ 0 (tied placement to low-complexity sequence)"],
        ["diff_gene",          "20–22%", "Hit on a different gene — paralog or shared exon"],
        ["same_gene",          "1–2%",   "Secondary confirms primary gene — mm rescue candidate"],
    ]
    simple_table(s, headers, rows,
                 top=Inches(1.75),
                 col_widths=[Inches(2.8), Inches(2.0), Inches(7.5)],
                 font_size=13)

    box(s, Inches(0.4), Inches(4.05), Inches(12.5), Inches(0.35),
        "AS score gap (primary − secondary) by class:", font_size=13, bold=True, color=DARK_GREY)

    bullets = [
        (0, "diff_gene:  median gap grows with read length (6 bp at ≤250 → 9 bp at ≥401).  Longer reads are better differentiated from alternative gene placements, not worse."),
        (0, "decoy_other:  median gap ≈ 0 across all bins.  These are tied placements to repeat/low-complexity sequence — a filter candidate, not a length effect."),
        (0, "same_gene:  median gap 0–3 bp (nearly tied, expected for the same locus).  mm rescue pool is real but tiny (~1–2% of secondaries)."),
    ]
    bullet_block(s, bullets, top=Inches(4.45), font_size=13)

    framing_box(s, "M2 compact cohort only (sbxpress illumina path). NH=1 fraction is flat at 91–92% across "
                "length bins in this cohort. Full-scale P5/P7 diagnostic run in progress — will confirm "
                "whether secondary-hit composition shifts when NH=1 does decline with length.")


def slide_m4_churn(prs):
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "M4 — Multimap mode changes almost nothing downstream",
               "Preliminary  ·  ~9.5M read subsets  ·  full-scale P5/P7 run in progress")

    box(s, Inches(0.4), Inches(1.2), Inches(6.0), Inches(0.35),
        "Baseline vs multimap arm (M2 compact):", font_size=13, bold=True, color=DARK_GREY)

    headers = ["Metric", "Baseline", "Multimap", "Δ"]
    rows = [
        ["Molecules counted", "5,798,640", "5,798,793", "+153  (0.003%)"],
        ["Cells called",      "26,303",    "26,305",    "+2"],
        ["Filtered UMI total","7,676,750", "7,677,264", "+514  (0.007%)"],
    ]
    simple_table(s, headers, rows,
                 top=Inches(1.6),
                 col_widths=[Inches(2.8), Inches(2.3), Inches(2.3), Inches(3.0)],
                 font_size=13)

    box(s, Inches(0.4), Inches(3.2), Inches(12.5), Inches(0.35),
        "Read-name-matched GX churn (all 3 datasets):", font_size=13, bold=True, color=DARK_GREY)

    headers2 = ["Dataset", "Reads", "Lost single-gene", "Lost → 0 / → 2+ genes", "Gained single-gene", "mm-tagged gains", "Net Δ"]
    rows2 = [
        ["M2 compact (sbxpress illumina)", "9,419,913", "501  (0.009%)", "475 / 26", "652  (0.018%)", "115 / 652 (17.6%)", "+151"],
        ["SBX P5 prelim (shard prefix)",   "9,497,915", "379  (0.007%)", "351 / 28", "624  (0.014%)", "89 / 624 (14.3%)",  "+232"],
        ["SBX P7 prelim (shard prefix)",   "9,513,521", "284  (0.005%)", "263 / 21", "591  (0.015%)", "87 / 591 (14.7%)",  "+303"],
    ]
    simple_table(s, headers2, rows2,
                 top=Inches(3.6),
                 col_widths=[Inches(2.6), Inches(1.4), Inches(1.7), Inches(1.85), Inches(1.7), Inches(1.85), Inches(1.2)],
                 font_size=11)

    box(s, Inches(0.4), Inches(5.85), Inches(12.5), Inches(0.7),
        "Key finding: ~93–95% of 'losses' drop to zero compatible genes (not multi-gene ambiguity). "
        "Only 14–18% of 'gains' carry the mm tag (genuine multimap rescue). "
        "The remaining churn traces to minimap2 choosing a slightly different primary alignment "
        "between --secondary=no and -N 10 runs — expected mapper behaviour.",
        font_size=13, color=DARK_GREY)

    framing_box(s, "Preliminary shard-prefix subsamples, not random draws or full-scale runs. "
                "Full-scale P5/P7 run in progress — confirms whether rates hold at 1.2B reads per orientation. "
                "These results support keeping --secondary=no as the production default.")


def slide_validation_status(prs):
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "Validation status", "What is done, what is pending")

    headers = ["Layer", "Target", "Status"]
    rows = [
        ["BAM → matrix",        "Feature order, barcode order, sparse entries, MEX/HDF5 layout",          "✓ Exact match validated"],
        ["Mapping / annotation","Tag assignment, multimap behaviour, gene/transcript evidence",             "⏳ Partial — compact cohort complete; full-scale in progress"],
        ["Molecule calling",    "CB+UB+GX groups, UMI correction, collision resolution, xf selection",    "⏳ Pending"],
        ["Cell calling",        "Called barcode set and filtered per-cell matrix values",                  "⏳ Pending"],
    ]
    simple_table(s, headers, rows,
                 top=Inches(1.3),
                 col_widths=[Inches(2.2), Inches(7.0), Inches(3.6)],
                 font_size=13)

    box(s, Inches(0.4), Inches(3.8), Inches(12.5), Inches(0.4),
        "Acceptance criteria:", font_size=13, bold=True, color=DARK_GREY)
    bullets = [
        (0, "• Exact equality where semantics are deterministic (e.g. BAM→matrix reconstruction)."),
        (0, "• Documented, explainable deviation everywhere else — not just 'close enough'."),
        (0, "• No SBX-specific pipeline change ships until the baseline validation is complete."),
    ]
    bullet_block(s, bullets, top=Inches(4.2), font_size=14)


def slide_next_steps(prs):
    s = blank_slide(prs)
    bg(s)
    header_bar(s, "Next steps", "In order")

    bullets = [
        (0, "① Full-scale P5/P7 diagnostic run  (in progress — launched 2026-07-30)"),
        (1, "NH=1-by-length at 1.2B reads per orientation."),
        (1, "Baseline vs multimap downstream impact at production scale."),
        (1, "Secondary-hit composition on real long SBX reads (not compact fixture)."),
        (0, ""),
        (0, "② Complete baseline validation"),
        (1, "Molecule calling: UMI correction, collision resolution, xf representative selection."),
        (1, "Cell calling: called barcode set and filtered matrix values vs CellRanger."),
        (0, ""),
        (0, "③ Fragment-length optimisation  (after ① and ②)"),
        (1, "Extend length bins below 250 bp; compute unmapped-rate-by-length from pre-alignment distribution."),
        (1, "Get reads-per-run vs target-molecule-length from device throughput data."),
        (1, "Only then claim a fragment-length sweet spot."),
        (0, ""),
        (0, "④ SBX-specific A/B experiments  (after ① ② ③)"),
        (1, "Gene-aware multimap rescue  ·  length-aware trimming/filtering  ·  revised collision model."),
    ]
    bullet_block(s, bullets, top=Inches(1.2), font_size=14)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="sbxpress_mapping_analysis.pptx")
    args = ap.parse_args()

    prs = new_prs()

    slide_title(prs)
    slide_core_message(prs)
    slide_pipeline(prs)
    slide_m1_paradox(prs)
    slide_m2_orientation(prs)
    slide_m3_secondary(prs)
    slide_m4_churn(prs)
    slide_validation_status(prs)
    slide_next_steps(prs)

    prs.save(args.out)
    print(f"saved {args.out}  ({len(prs.slides)} slides)")


if __name__ == "__main__":
    main()
