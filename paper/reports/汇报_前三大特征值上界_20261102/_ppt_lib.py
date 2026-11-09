# -*- coding: utf-8 -*-
"""自编译 PPT：不依赖任何模板 skill，直接用 python-pptx 版式引擎生成。"""
import os
from pptx import Presentation
from pptx.util import Inches, Pt, Emu
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import qn
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(HERE, "assets")

INK   = RGBColor.from_string("1A1A1A")
BLUE  = RGBColor.from_string("1F3864")
BLUE2 = RGBColor.from_string("2F5597")
BLUEL = RGBColor.from_string("E8EDF6")
GREEN = RGBColor.from_string("2E7D32")
GREENL= RGBColor.from_string("E7F2E8")
AMBER = RGBColor.from_string("B26A00")
AMBERL= RGBColor.from_string("FBF0DF")
RED   = RGBColor.from_string("B3261E")
REDL  = RGBColor.from_string("FBE9E7")
GREY  = RGBColor.from_string("5B6470")
GREYL = RGBColor.from_string("F2F4F7")
LINE  = RGBColor.from_string("D3D8DF")
WHITE = RGBColor.from_string("FFFFFF")
FONT  = "Microsoft YaHei"

prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)
BLANK = prs.slide_layouts[6]
SW, SH = 13.333, 7.5

def _setfont(run, size, color=INK, bold=False, italic=False, name=FONT):
    f = run.font
    f.size = Pt(size); f.bold = bold; f.italic = italic
    f.color.rgb = color; f.name = name
    rPr = run._r.get_or_add_rPr()
    for tag in ("a:latin", "a:ea", "a:cs"):
        el = rPr.find(qn(tag))
        if el is None:
            el = rPr.makeelement(qn(tag), {}); rPr.append(el)
        el.set("typeface", name)

def tb(slide, x, y, w, h, anchor=MSO_ANCHOR.TOP):
    box = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = box.text_frame
    tf.word_wrap = True
    tf.vertical_anchor = anchor
    tf.margin_left = 0; tf.margin_right = 0; tf.margin_top = 0; tf.margin_bottom = 0
    return box, tf

def para(tf, first=False):
    if first and not tf.paragraphs[0].runs:
        return tf.paragraphs[0]
    return tf.add_paragraph()

def line(tf, text, size=15, color=INK, bold=False, italic=False, space_after=6,
         space_before=0, align=PP_ALIGN.LEFT, first=False, indent=0.0):
    p = para(tf, first)
    p.alignment = align
    p.space_after = Pt(space_after); p.space_before = Pt(space_before)
    if indent:
        p.level = 0
        pPr = p._pPr if p._pPr is not None else p.get_or_add_pPr()
        pPr.set("marL", str(int(indent * 914400)))
        pPr.set("indent", str(int(-0.18 * 914400)))
    r = p.add_run(); r.text = text
    _setfont(r, size, color, bold, italic)
    return p

def rich(tf, parts, size=15, space_after=6, first=False, align=PP_ALIGN.LEFT):
    """parts: list of (text, color, bold)"""
    p = para(tf, first); p.alignment = align; p.space_after = Pt(space_after)
    for text, color, bold in parts:
        r = p.add_run(); r.text = text
        _setfont(r, size, color, bold)
    return p

def rect(slide, x, y, w, h, fill=None, linec=None, lw=1.0, shape=MSO_SHAPE.RECTANGLE, radius=None):
    sp = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    if fill is None:
        sp.fill.background()
    else:
        sp.fill.solid(); sp.fill.fore_color.rgb = fill
    if linec is None:
        sp.line.fill.background()
    else:
        sp.line.color.rgb = linec; sp.line.width = Pt(lw)
    sp.shadow.inherit = False
    if radius is not None and shape == MSO_SHAPE.ROUNDED_RECTANGLE:
        try:
            sp.adjustments[0] = radius
        except Exception:
            pass
    sp.text_frame.word_wrap = True
    return sp

def chip(slide, x, y, text, color, fill, w=None, h=0.32, size=11.5):
    w = w or (0.30 + 0.155 * len(text))
    sp = rect(slide, x, y, w, h, fill=fill, shape=MSO_SHAPE.ROUNDED_RECTANGLE, radius=0.5)
    tf = sp.text_frame; tf.margin_left = Inches(0.08); tf.margin_right = Inches(0.08)
    tf.margin_top = 0; tf.margin_bottom = 0
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = text; _setfont(r, size, color, True)
    return x + w

PAGE = [0]
def base(title, section=None, subtitle=None):
    PAGE[0] += 1
    s = prs.slides.add_slide(BLANK)
    rect(s, 0, 0, SW, SH, fill=WHITE)
    rect(s, 0, 0, 0.16, SH, fill=BLUE)
    rect(s, 0.55, 0.42, 12.25, 0.035, fill=LINE)
    y = 0.62
    if section:
        _, tf = tb(s, 0.58, y - 0.06, 8.0, 0.34)
        line(tf, section, 12.5, BLUE2, True, first=True, space_after=0)
        y += 0.34
    _, tf = tb(s, 0.58, y, 12.2, 0.62)
    line(tf, title, 27, BLUE, True, first=True, space_after=0)
    y += 0.66
    if subtitle:
        _, tf = tb(s, 0.58, y, 12.2, 0.34)
        line(tf, subtitle, 13, GREY, False, True, first=True, space_after=0)
        y += 0.36
    rect(s, 0.58, 7.02, 12.2, 0.012, fill=LINE)
    _, tf = tb(s, 0.58, 7.08, 8.0, 0.3)
    line(tf, "谱极值图论 · 硕士学位论文汇报 · 2026-11-09", 10, GREY, first=True, space_after=0)
    _, tf = tb(s, 10.6, 7.08, 2.2, 0.3)
    line(tf, "%d / 26" % PAGE[0], 10, GREY, first=True, space_after=0, align=PP_ALIGN.RIGHT)
    return s, y + 0.16

def pic(slide, name, x, y, w=None, h=None):
    path = os.path.join(A, name)
    iw, ih = Image.open(path).size
    ar = ih / iw
    if w is not None and h is None:
        h = w * ar
    elif h is not None and w is None:
        w = h / ar
    slide.shapes.add_picture(path, Inches(x), Inches(y), Inches(w), Inches(h))
    return w, h

def table(slide, x, y, w, rows, widths, size=11.5, header=True, rowh=0.34):
    n_r, n_c = len(rows), len(rows[0])
    gf = slide.shapes.add_table(n_r, n_c, Inches(x), Inches(y), Inches(w), Inches(rowh * n_r))
    t = gf.table
    t.first_row = header
    for j, cw in enumerate(widths):
        t.columns[j].width = Inches(cw)
    for i, row in enumerate(rows):
        t.rows[i].height = Inches(rowh)
        for j, val in enumerate(row):
            c = t.cell(i, j)
            c.margin_left = Inches(0.08); c.margin_right = Inches(0.06)
            c.margin_top = Inches(0.02); c.margin_bottom = Inches(0.02)
            c.vertical_anchor = MSO_ANCHOR.MIDDLE
            c.fill.solid()
            if i == 0 and header:
                c.fill.fore_color.rgb = BLUE
            else:
                c.fill.fore_color.rgb = WHITE if i % 2 else GREYL
            tf = c.text_frame; tf.word_wrap = True
            p = tf.paragraphs[0]
            r = p.add_run(); r.text = str(val)
            _setfont(r, size, WHITE if (i == 0 and header) else INK, bold=(i == 0 and header))
    return gf
