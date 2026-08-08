#!/usr/bin/env python3
"""Typeset PAPER.md into the two-column manuscript PDF.

Regenerates hf-crosstalk-manuscript.pdf from PAPER.md and results/figures/.
Committed because the manuscript is a deliverable: an earlier version of this
code lived only in a notebook kernel and was lost when the session reset, which
left the PDF the one artifact in the package that could not be rebuilt.

Usage:
    python build_manuscript.py [PAPER.md] [out.pdf]

Requires: reportlab, pypdf, pillow. Fonts: Times New Roman and Arial (macOS
paths below; override with --fontdir or set FONTDIR).
"""
import os
import re
import sys

from PIL import Image as PILImage
from reportlab.lib.enums import TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (BaseDocTemplate, Frame, Image, KeepTogether,
                                NextPageTemplate, PageBreak, PageTemplate,
                                Paragraph, Spacer)

PW, PH = A4
TM, BM, LM, RM = 17 * mm, 17 * mm, 17 * mm, 17 * mm
GUT = 6 * mm
COLW = (PW - LM - RM - GUT) / 2
FULLW = PW - LM - RM
FONTDIRS = ["/System/Library/Fonts/Supplemental/", "/Library/Fonts/",
            "/usr/share/fonts/truetype/msttcorefonts/"]

# Figure placement offsets, in paragraphs after the figure's first citation.
# Tuned by the search in tune_offsets(); see README for how to re-tune.
OFFSETS = {2: 3, 3: -4, 6: 2, 7: 2, 5: -2}
# Figures given a dedicated page rather than an inline strip (tall multi-panel).
TALL_FRACTION = 0.62


def _font(name):
    for d in FONTDIRS:
        p = os.path.join(d, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError("font not found: %s (looked in %s)" % (name, FONTDIRS))


def register_fonts():
    for alias, fname in [("TNR", "Times New Roman.ttf"),
                         ("TNR-B", "Times New Roman Bold.ttf"),
                         ("TNR-I", "Times New Roman Italic.ttf"),
                         ("ARI", "Arial.ttf"), ("ARI-B", "Arial Bold.ttf")]:
        pdfmetrics.registerFont(TTFont(alias, _font(fname)))
    pdfmetrics.registerFontFamily("TNR", normal="TNR", bold="TNR-B", italic="TNR-I")


def styles():
    base = dict(fontName="TNR", fontSize=9.2, leading=11.4, alignment=TA_JUSTIFY,
                spaceShrinkage=0.15, hyphenationLang="en_US",
                embeddedHyphenation=1, hyphenationMinWordLength=5,
                allowWidows=0, allowOrphans=0)
    S = {
        "title": ParagraphStyle("title", fontName="TNR-B", fontSize=15.5, leading=18.5,
                                alignment=TA_LEFT, spaceAfter=7),
        "h2": ParagraphStyle("h2", fontName="TNR-B", fontSize=10.4, leading=12.4,
                             spaceBefore=7, spaceAfter=2.5, alignment=TA_LEFT),
        "h3": ParagraphStyle("h3", fontName="TNR-B", fontSize=9.5, leading=11.5,
                             spaceBefore=5.5, spaceAfter=2, alignment=TA_LEFT),
        "p": ParagraphStyle("p", **base),
        "li": ParagraphStyle("li", leftIndent=9, bulletIndent=2, **base),
        "cap": ParagraphStyle("cap", fontName="TNR", fontSize=7.8, leading=9.4,
                              alignment=TA_JUSTIFY, spaceBefore=3),
        "ref": ParagraphStyle("ref", fontName="TNR", fontSize=7.9, leading=9.5,
                              leftIndent=9, firstLineIndent=-9, alignment=TA_LEFT),
        "abs": ParagraphStyle("abs", fontName="TNR", fontSize=8.6, leading=10.6,
                              alignment=TA_JUSTIFY, spaceShrinkage=0.15,
                              allowWidows=0, allowOrphans=0),
    }
    for k in ("p", "li", "ref", "cap", "abs"):
        S[k].bulletFontName = "TNR"
    return S


SUP = str.maketrans("0123456789+-=()n", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079"
                                        "\u207a\u207b\u207c\u207d\u207e\u207f")


def inline(t):
    """Markdown inline spans -> reportlab markup."""
    t = t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t, flags=re.S)
    t = re.sub(r"(?<!\*)\*([^*\n]+?)\*(?!\*)", r"<i>\1</i>", t)
    t = re.sub(r"`([^`]+?)`", r'<font name="ARI" size="8.2">\1</font>', t)
    # non-breaking hyphen in compounds reportlab breaks while dropping the hyphen
    for w in ("fibroblast-derived", "cardiac-plasma"):
        t = t.replace(w, w.replace("-", "\u2011"))
    return t


def md_blocks(text):
    """Parse the manuscript into (kind, value) blocks."""
    out, lines, i = [], text.split("\n"), 0
    while i < len(lines):
        s = lines[i].strip()
        if not s:
            i += 1
            continue
        if s.startswith("# "):
            out.append(("title", s[2:]))
        elif s.startswith("### "):
            out.append(("h3", s[4:]))
        elif s.startswith("## "):
            out.append(("h2", s[3:]))
        elif s.startswith("!["):
            m = re.match(r"!\[Figure (\d)\]\(([^)]*)\)", s)
            if m:
                out.append(("fig", (int(m.group(1)), m.group(2))))
        elif s.startswith(("- ", "* ")):
            out.append(("li", s[2:]))
        elif s.startswith("|") and i + 1 < len(lines) and set(lines[i + 1].strip()) <= set("|-: "):
            rows, j = [], i
            while j < len(lines) and lines[j].strip().startswith("|"):
                rows.append([c.strip() for c in lines[j].strip().strip("|").split("|")])
                j += 1
            out.append(("table", [r for r in rows if not set("".join(r)) <= set("-: ")]))
            i = j
            continue
        else:
            para = [s]
            j = i + 1
            while j < len(lines) and lines[j].strip() and not lines[j].lstrip().startswith(
                    ("#", "!", "|", "- ", "* ")):
                para.append(lines[j].strip())
                j += 1
            out.append(("p", " ".join(para)))
            i = j
            continue
        i += 1
    return out


class SpanDoc(BaseDocTemplate):
    """Two-column body with a full-width banner frame for title and abstract."""

    def __init__(self, fn, banner_h, **kw):
        BaseDocTemplate.__init__(self, fn, pagesize=A4, leftMargin=LM,
                                 rightMargin=RM, topMargin=TM, bottomMargin=BM, **kw)
        body_top = PH - TM - banner_h
        banner = Frame(LM, body_top, FULLW, banner_h, 0, 0, 0, 0, id="banner")
        c1 = Frame(LM, BM, COLW, body_top - BM, 0, 0, 0, 0, id="c1")
        c2 = Frame(LM + COLW + GUT, BM, COLW, body_top - BM, 0, 0, 0, 0, id="c2")
        f1 = Frame(LM, BM, COLW, PH - TM - BM, 0, 0, 0, 0, id="f1")
        f2 = Frame(LM + COLW + GUT, BM, COLW, PH - TM - BM, 0, 0, 0, 0, id="f2")
        full = Frame(LM, BM, FULLW, PH - TM - BM, 0, 0, 0, 0, id="full")
        self.addPageTemplates([
            PageTemplate("first", [banner, c1, c2], onPage=self._num),
            PageTemplate("rest", [f1, f2], onPage=self._num),
            PageTemplate("plate", [full], onPage=self._num),
        ])

    @staticmethod
    def _num(canvas, doc):
        canvas.saveState()
        canvas.setFont("TNR", 8)
        canvas.drawCentredString(PW / 2, BM * 0.45, str(canvas.getPageNumber()))
        canvas.restoreState()

    def handle_pageBreak(self, slow=None):
        # A PageBreak in a two-column template abandons the remaining column.
        # Suppress it when the current page has no content yet.
        f = self.frame
        if f and f.id in ("c1", "f1") and not f._atTop:
            return BaseDocTemplate.handle_pageBreak(self, slow)
        if f and f._atTop and f.id in ("c1", "f1"):
            return
        return BaseDocTemplate.handle_pageBreak(self, slow)


def figure_flowables(num, path, caption, S, width):
    w, h = PILImage.open(path).size
    img = Image(path, width=width, height=width * h / w)
    return [img, Spacer(1, 3), Paragraph(inline(caption), S["cap"])]


def build(md_text, out, figdir="results/figures", offsets=None, art2file=None):
    """Render md_text to `out`. Returns the banner height used."""
    S = styles()
    offsets = dict(OFFSETS if offsets is None else offsets)
    blocks = md_blocks(md_text)

    caps = {int(re.match(r"\*\*Figure (\d)", v.strip()).group(1)): v
            for k, v in blocks if k == "p" and v.lstrip().startswith("**Figure")}
    figs = {}
    for k, v in blocks:
        if k == "fig":
            n, ref = v
            fn = (art2file or {}).get(n) or ref
            figs[n] = os.path.join(figdir, os.path.basename(fn))

    # Banner = everything before section 1 (title, authors, highlights, abstract,
    # keywords, abbreviations) rendered full width. Measured, not assumed.
    BANNER_HEADS = {"highlights", "abstract", "keywords", "abbreviations"}
    banner = []
    for k, v in blocks:
        if k == "title":
            banner.append(Paragraph(inline(v), S["title"]))
            continue
        if k == "h2":
            if v.strip().lower() in BANNER_HEADS:
                banner.append(Paragraph(v.strip().upper(), S["h3"]))
                continue
            break
        if banner and k in ("p", "li") and not v.lstrip().startswith("**Figure"):
            banner.append(Paragraph(inline(v), S["abs"] if k == "p" else S["li"]))
    bh = sum(f.wrap(FULLW, 10 ** 6)[1] + getattr(f, "spaceAfter", 0) for f in banner) + 8
    # Cap the banner so the body still gets usable column space on page 1.
    bh = min(bh, (PH - TM - BM) * 0.86)

    story, in_abs, pending, para_i = [], False, {}, {}
    for k, v in blocks:
        if k == "title":
            story.append(Paragraph(inline(v), S["title"]))
            continue
        if k == "h2" and v.strip().lower() in BANNER_HEADS:
            story += [Paragraph(v.strip().upper(), S["h3"])]
            in_abs = True
            continue
        if k == "h2" and in_abs:
            in_abs = False
            story.append(NextPageTemplate("rest"))
        if k == "fig":
            continue
        if k == "p" and v.lstrip().startswith("**Figure"):
            continue
        if k == "title":
            continue
        sty = S["abs"] if in_abs else S[{"h2": "h2", "h3": "h3", "li": "li"}.get(k, "p")]
        if k == "table":
            for row in v:
                story.append(Paragraph(inline(" | ".join(row)), S["cap"]))
            continue
        if k == "h2" and v.startswith("References"):
            story.append(Paragraph(inline(v), S["h2"]))
            sty = S["ref"]
            continue
        story.append(Paragraph(inline(v), sty))
        if k == "p":
            para_i[len(story) - 1] = True
            # release any figure whose citation has been passed by its offset
            for n in sorted(list(pending)):
                cited, off = pending[n]
                if len([i for i in para_i if i > cited]) >= off:
                    story += figure_flowables(n, figs[n], caps[n], S, COLW)
                    del pending[n]
            for n in sorted(figs):
                if n not in caps or n in pending:
                    continue
                if re.search(r"Fig\. %d\b" % n, v):
                    pending[n] = (len(story) - 1, max(0, offsets.get(n, 0)))
    for n in sorted(pending):
        story += figure_flowables(n, figs[n], caps[n], S, COLW)

    doc = SpanDoc(out, bh)
    doc.build(story)
    return bh


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "PAPER.md"
    out = sys.argv[2] if len(sys.argv) > 2 else "hf-crosstalk-manuscript.pdf"
    register_fonts()
    md_text = open(src).read()
    # PAPER.md embeds figures by artifact id so the web renderer resolves them;
    # this manifest maps figure number -> filename on disk for the PDF build.
    a2f = {}
    for cand in ("figure_manifest.json", "results/figure_manifest.json"):
        if os.path.exists(cand):
            import json
            a2f = {int(k): v for k, v in json.load(open(cand)).items()}
            break
    figdir = "results/figures" if os.path.isdir("results/figures") else "figs"
    bh = build(md_text, out, figdir=figdir, art2file=a2f)
    from pypdf import PdfReader
    r = PdfReader(out)
    print("%s -> %s | %d pages | banner %.0f mm" % (src, out, len(r.pages), bh / mm))


if __name__ == "__main__":
    main()
