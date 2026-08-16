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

fa = pd.read_csv(os.path.join(TABLES, "floor_audit.csv"))
dec = pd.read_csv(os.path.join(TABLES, "twofactor_decomposition.csv"))

fig = plt.figure(figsize=(7.2, 3.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1], wspace=0.34)
axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1])

# Panel A: floor audit — log scale, points not bars (per §3.3)
f = fa.sort_values("max_capacity").reset_index(drop=True)
floorv = float(f["floor"].iloc[0])
xx = np.arange(len(f))
c = [ALARM if not p else FOCAL for p in f["passes"]]
axA.scatter(xx, f["max_capacity"].clip(lower=0.2), s=7, c=c, linewidth=0)
axA.axhline(floorv, color="0.3", lw=1, ls="--")
axA.set_yscale("log")
axA.annotate(f"floor = {floorv:.0f}\n(0.5% of GAPDH)", (len(f)*0.42, floorv), xytext=(0,-26),
             textcoords="offset points", fontsize=6, color="0.3", ha="left")
axA.set_xlabel("Ligand–receptor interaction (ranked)")
axA.set_ylabel("Max receptor capacity across\n30 tumor types (RSEM norm. count)")
axA.set_title("28 of 188 interactions have receptors that are\nuntranscribed in every tumor type", loc="left")
for nm,yv in [("PENK→OPRD1",3.10),("PENK→OPRM1",34.4),("PENK→OPRK1",None)]:
    row=f[f["interaction"].str.contains(nm.split("→")[1], na=False)]
    if len(row):
        i=row.index[0]
        axA.annotate(nm, (i, row["max_capacity"].iloc[0]), xytext=(6,-2), textcoords="offset points",
                     fontsize=5.5, color=ALARM)
axA.set_ylim(0.15, 2e5)
from matplotlib.lines import Line2D
axA.legend(handles=[Line2D([],[],marker="o",ls="",color=ALARM,ms=4,label="below floor (dropped)"),
                    Line2D([],[],marker="o",ls="",color=FOCAL,ms=4,label="scored")],
           loc="lower right", frameon=False, fontsize=6)

# Panel B: two-factor decomposition
dd = dec.sort_values("availability_weighted_capacity", ascending=False)
top = dd.head(12)
yb = np.arange(len(top))[::-1]
axB.scatter(top["mean_capacity_z"], yb, s=30, color=GREY, label="unweighted capacity", zorder=3)
axB.scatter(top["availability_weighted_capacity"], yb, s=30, color=FOCAL, label="availability-weighted", zorder=3)
for i,(a,b) in enumerate(zip(top["mean_capacity_z"], top["availability_weighted_capacity"])):
    axB.plot([a,b],[yb[i],yb[i]], color="0.8", lw=.9, zorder=1)
axB.axvline(0, color="0.35", lw=.8, ls=":")
axB.set_yticks(yb); axB.set_yticklabels(top["tumor_type"], fontsize=6)
axB.set_xlabel("Receiver capacity (z)")
axB.set_title("Ligand availability is constant across tumors;\nonly receiver capacity varies", loc="left")
axB.legend(loc="lower right", frameon=False, fontsize=6)
axB.margins(y=0.02)
panel_letter(axA,"a"); panel_letter(axB,"b")
fig.savefig(os.path.join(OUT, "fig3_floor_and_decomposition.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig3_floor_and_decomposition.png"), dpi=300, bbox_inches="tight")
print("saved")