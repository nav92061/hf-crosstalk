import os, re, json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# skill:figure-style kernel.py (auto-injected on skill load)
META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")

    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


apply_figure_style(sizes=(8, 7, 6))

FOCAL = "#1a5fb4"
GREY = "#9a9996"
ALARM = "#c01c28"
ACC2 = "#e5a50a"

lig = pd.read_csv(os.path.join(TABLES, "ligand_tissue_specificity.csv"))
rtn = pd.read_csv(os.path.join(TABLES, "receptor_tumor_vs_normal.csv"))

fig = plt.figure(figsize=(7.4, 3.6))
gsy = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.34)
axA = fig.add_subplot(gsy[0])
axB = fig.add_subplot(gsy[1])

# Panel A: ligand cardiac specificity — heart rank vs enrichment ratio
L = lig.dropna(subset=["ratio_best_heart_chamber", "best_heart_rank"]).copy()
iscard = L["classification"].str.startswith("cardiac-enriched")
lowab = L["classification"].str.contains("low-abundance")
axA.scatter(L.loc[~iscard, "best_heart_rank"], L.loc[~iscard, "ratio_best_heart_chamber"].clip(lower=0.05),
            s=18, c=GREY, alpha=.7, linewidth=0, label="not cardiac-enriched")
axA.scatter(L.loc[iscard & ~lowab, "best_heart_rank"], L.loc[iscard & ~lowab, "ratio_best_heart_chamber"],
            s=52, c=FOCAL, zorder=5, edgecolor="white", linewidth=.6, label="cardiac-enriched, abundant")
axA.scatter(L.loc[lowab, "best_heart_rank"], L.loc[lowab, "ratio_best_heart_chamber"],
            s=42, facecolor="none", edgecolor=ACC2, linewidth=1.3, zorder=5, label="cardiac-enriched but <1 TPM")
axA.set_yscale("log")
axA.set_xscale("log")
axA.axhline(5, color="0.45", lw=.8, ls="--")
axA.text(52, 5.6, "5× enrichment", fontsize=5.5, color="0.45", ha="right")
for g in ["NPPA", "NPPB", "GDF6"]:
    r_ = L[L["gene"] == g]
    if len(r_):
        dy = 7 if g in ("NPPA",) else (-9 if g in ("NPPB",) else 4)
        axA.annotate(g, (r_["best_heart_rank"].iloc[0], r_["ratio_best_heart_chamber"].iloc[0]),
                     xytext=(5, dy), textcoords="offset points", fontsize=5.4, style="italic",
                     color=FOCAL if g in ("NPPA", "NPPB") else (ACC2 if g == "GDF6" else "0.4"))
axA.set_xlabel("Heart's rank among 52 GTEx tissues (log)")
axA.set_ylabel("Heart / median-tissue expression ratio (log)")
axA.set_title("Under a rank-based definition only\nNPPA and NPPB are cardiac-enriched", loc="left")
axA.legend(loc="lower left", frameon=False, fontsize=5.6)

offs = {"WNT9A": (6, 7), "CCL5": (6, -9), "MDK": (-30, -2), "POSTN": (6, 3), "PENK": (-26, 4), "GDF15": (6, -3)}
for g, (dx, dy) in offs.items():
    r_ = L[L["gene"] == g]
    if len(r_):
        axA.annotate(g, (r_["best_heart_rank"].iloc[0], r_["ratio_best_heart_chamber"].iloc[0]),
                     xytext=(dx, dy), textcoords="offset points", fontsize=5.4, style="italic", color="0.4")

# Panel B: NPR1 tumor vs adjacent normal
N = rtn[(rtn["receptor"] == "NPR1") & rtn["log2FC_bg_corrected_within_tcga"].notna()].copy()
N = N.sort_values("log2FC_bg_corrected_within_tcga")
yy = np.arange(len(N))
sig = N["fdr_within_tcga"] < 0.05
axB.barh(yy, N["log2FC_bg_corrected_within_tcga"], color=[ALARM if s else GREY for s in sig],
         height=.68, edgecolor="white", linewidth=.5)
axB.axvline(0, color="0.3", lw=.9)
axB.set_yticks(yy)
axB.set_yticklabels(N["tumor_type"], fontsize=6)
axB.set_xlabel("log₂FC tumor vs adjacent normal\n(background-corrected, unit-matched)", fontsize=6.5)
axB.set_title("NPR1 is depleted in tumor in all 15\ntypes with adjacent normals", loc="left")
axB.text(0.04, 0.97, f"all 15 types FDR < 0.05\nmedian log₂FC = {N['log2FC_bg_corrected_within_tcga'].median():.2f}",
         transform=axB.transAxes, fontsize=6, color=ALARM, va="top", ha="left")
axB.margins(y=0.02)

panel_letter(axA, "a")
panel_letter(axB, "b")
fig.savefig(os.path.join(OUT, "fig5_specificity_checks.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig5_specificity_checks.png"), dpi=300, bbox_inches="tight")