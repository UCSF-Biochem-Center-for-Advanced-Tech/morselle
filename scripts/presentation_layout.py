"""Shared slide layout and palette for the sbxpress reports."""
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

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

