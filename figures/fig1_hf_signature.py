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
import json
from scipy import stats as st

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

de1 = pd.read_csv(os.path.join(TABLES, "hf_de_GSE116250.csv"))
de2 = pd.read_csv(os.path.join(TABLES, "hf_de_GSE141910.csv"))

m = de1.merge(de2, on="gene", suffixes=("_d", "_v"))

rp = json.load(open("/Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/46dfe399-9e00-4bd1-b3b2-19547e9394eb/v8d4419d5_hf_replication_stats.json"))

allsh = m.dropna(subset=["log2FC_d", "log2FC_v"])
rho_all = st.spearmanr(allsh["log2FC_d"], allsh["log2FC_v"]).statistic

both = m[(m["fdr_d"] < 0.05) & (m["fdr_v"] < 0.05)]
rho_both = st.spearmanr(both["log2FC_d"], both["log2FC_v"]).statistic
r_both = st.pearsonr(both["log2FC_d"], both["log2FC_v"]).statistic
agree = np.sign(both["log2FC_d"]) == np.sign(both["log2FC_v"])

fig = plt.figure(figsize=(7.2, 3.4))
gs_ = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1], wspace=0.36)
axA = fig.add_subplot(gs_[0])
axB = fig.add_subplot(gs_[1])
axC = fig.add_subplot(gs_[2])

# A: discovery volcano
v = de1.dropna(subset=["log2FC", "fdr"]).copy()
v["nl"] = -np.log10(v["fdr"].clip(lower=1e-300))
sig = v["fdr"] < 0.05
axA.scatter(v.loc[~sig, "log2FC"], v.loc[~sig, "nl"], s=2, c="0.82", linewidth=0)
axA.scatter(v.loc[sig, "log2FC"], v.loc[sig, "nl"], s=2, c=FOCAL, linewidth=0, alpha=.55)
offs = {"NPPA": (6, -6), "NPPB": (6, 4), "PENK": (-26, -4), "WNT9A": (4, 3), "MYH6": (-24, 2), "POSTN": (5, -2)}
for g, (dx, dy) in offs.items():
    r_ = v[v["gene_symbol"] == g]
    if len(r_):
        axA.scatter(r_["log2FC"], r_["nl"], s=16, c=ALARM, zorder=5, linewidth=0)
        axA.annotate(g, (r_["log2FC"].iloc[0], r_["nl"].iloc[0]), xytext=(dx, dy),
                     textcoords="offset points", fontsize=5.5, style="italic", color=ALARM)
axA.axhline(-np.log10(0.05), color="0.4", lw=.7, ls=":")
axA.set_xlabel("log₂ fold change (failing vs non-failing)")
axA.set_ylabel("−log₁₀ FDR")
axA.set_title("Discovery: GSE116250\n(50 failing, 14 non-failing)", loc="left")

# B: effect-size concordance
axB.scatter(both.loc[agree, "log2FC_d"], both.loc[agree, "log2FC_v"], s=2.5, c=FOCAL, linewidth=0, alpha=.5)
axB.scatter(both.loc[~agree, "log2FC_d"], both.loc[~agree, "log2FC_v"], s=2.5, c=ALARM, linewidth=0, alpha=.5)
lim = [both[["log2FC_d", "log2FC_v"]].min().min() - .3, both[["log2FC_d", "log2FC_v"]].max().max() + .3]
axB.plot(lim, lim, color="0.4", lw=.7, ls="--")
axB.axhline(0, color="0.8", lw=.5)
axB.axvline(0, color="0.8", lw=.5)
axB.set_xlim(lim); axB.set_ylim(lim)
axB.set_xlabel("log₂FC discovery")
axB.set_ylabel("log₂FC validation")
axB.set_title(f"Replication: ρ = {rho_both:.2f} in genes\nsignificant in both cohorts", loc="left")
axB.text(.03, .95, f"88.1% same direction\nn = {len(both):,} sig in both\n(ρ = {rho_all:.2f} over all {len(allsh):,})",
         transform=axB.transAxes, fontsize=6, va="top", color="0.35")

# C: concordance strengthens with effect size + permutation control
bins = [(0, .5), (.5, 1), (1, 1.5), (1.5, 10)]
xs = []; ys = []
for lo, hi in bins:
    s = both[(both["log2FC_d"].abs() >= lo) & (both["log2FC_d"].abs() < hi)]
    if len(s) > 5:
        xs.append(f"{lo}–{hi if hi < 10 else '∞'}")
        ys.append((np.sign(s["log2FC_d"]) == np.sign(s["log2FC_v"])).mean() * 100)
axC.bar(np.arange(len(xs)), ys, color=FOCAL, width=.6, edgecolor="white", linewidth=.6)
axC.axhline(50, color=ALARM, lw=.9, ls="--")
axC.annotate("chance", (len(xs) - .4, 51), fontsize=5.5, color=ALARM, ha="right")
for i, yv in enumerate(ys):
    axC.annotate(f"{yv:.0f}", (i, yv), xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6)
axC.set_xticks(np.arange(len(xs)))
axC.set_xticklabels(xs, fontsize=6)
axC.set_xlabel("|log₂FC| in discovery")
axC.set_ylabel("% same direction in validation")
axC.set_ylim(0, 108)
axC.set_title("Concordance rises with effect size;\npermuted labels give 0 hits", loc="left")

for a, l in [(axA, "a"), (axB, "b"), (axC, "c")]:
    panel_letter(a, l)

fig.savefig(os.path.join(OUT, "fig1_hf_signature.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig1_hf_signature.png"), dpi=300, bbox_inches="tight")