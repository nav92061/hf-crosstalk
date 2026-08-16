#!/usr/bin/env python3
"""Regenerate every manuscript figure from the released result tables.

Writes a vector PDF and a 300 dpi PNG for each figure. The PDF is the
submission format: Type-42 fonts are embedded (Elsevier requires embedded
fonts in EPS/PDF artwork), so the figure is resolution-independent at any
print size rather than upsampled from a raster.

    python figures/build_figures.py            # -> results/figures/
    python figures/build_figures.py /some/dir  # -> /some/dir

Inputs are results/tables/*.csv from this repository. No network access and
no language-model call is required.
"""
import os
import sys
import glob
import runpy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TABLES = os.path.join(ROOT, "results", "tables")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "results", "figures")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Vector output with embedded fonts, not paths: keeps text selectable and
# searchable for the typesetter and satisfies the embedded-font requirement.
matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})


def apply_figure_style(sizes=(8, 7, 6)):
    """House style: no top/right spines, grey grid, sans-serif, tight labels."""
    t, l, s = sizes
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": l,
        "axes.titlesize": t, "axes.labelsize": l,
        "xtick.labelsize": s, "ytick.labelsize": s, "legend.fontsize": s,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "figure.dpi": 110,
        "axes.titlelocation": "left", "axes.titlepad": 8.0,
    })


def panel_letter(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top", ha="left")


def artifact_path(_vid):
    raise RuntimeError(
        "artifact_path() is a platform helper and is unavailable outside the "
        "original session; the scripts read results/tables/ instead."
    )


def main():
    os.makedirs(OUT, exist_ok=True)
    scripts = sorted(p for p in glob.glob(os.path.join(HERE, "fig*.py"))
                     if not p.endswith("build_figures.py"))
    ns = dict(TABLES=TABLES, OUT=OUT, apply_figure_style=apply_figure_style,
              panel_letter=panel_letter, artifact_path=artifact_path)
    made = []
    for s in scripts:
        g = dict(ns)
        g["__name__"] = "__main__"
        g["__file__"] = s
        try:
            exec(compile(open(s).read(), s, "exec"), g)
            plt.close("all")
            made.append(os.path.basename(s))
        except Exception as exc:
            print("FAILED %s: %s: %s" % (os.path.basename(s), type(exc).__name__, exc))
    pdfs = sorted(os.path.basename(p) for p in glob.glob(os.path.join(OUT, "*.pdf")))
    print("scripts run: %d of %d" % (len(made), len(scripts)))
    print("vector PDFs written: %d" % len(pdfs))
    for p in pdfs:
        print("   ", p)


if __name__ == "__main__":
    main()
