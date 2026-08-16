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

dep = pd.read_csv(os.path.join(TABLES, "depmap_lineage_dependency.csv"))
cor = pd.read_csv(os.path.join(TABLES, "depmap_compatibility_correlation.csv"))
pan = pd.read_csv(os.path.join(TABLES, "cell_line_panel.csv"))

fig = plt.figure(figsize=(7.6, 3.6))
gsx = fig.add_gridspec(1, 3, width_ratios=[1.1, 1, 1], wspace=0.48)
axA = fig.add_subplot(gsx[0])
axB = fig.add_subplot(gsx[1])
axC = fig.add_subplot(gsx[2])
rng = np.random.default_rng(3)
groups = [("Essential\ncontrols", ["RPL3", "PLK1", "POLR2B", "EEF2", "POLR2D", "RPL7"], "0.35"),
          ("Dependency\nreceptors", ["NCL", "ITGAV", "EGFR", "ITGB5", "ITGB1", "ERBB2", "FGFR1", "IGF1R"], ACC2),
          ("HF-axis\nreceptors", ["NPR1", "FZD1", "FZD2", "FZD5", "FZD6", "FZD7", "FZD10", "LRP5", "LRP6"], FOCAL)]
for i, (lab, gl, col) in enumerate(groups):
    v = dep[dep["gene"].isin(gl)]["mean_gene_effect"].dropna().values
    axA.scatter(np.full(len(v), i) + rng.normal(0, .075, len(v)), v, s=7, color=col, alpha=.5, linewidth=0)
    axA.hlines(np.median(v), i - .3, i + .3, color=col, lw=2.2, zorder=5)
    axA.annotate(f"{np.median(v):.2f}", (i, np.median(v)), xytext=(0, 5), textcoords="offset points",
                 ha="center", fontsize=5.5, color=col, weight="bold")
axA.axhline(-0.5, color=ALARM, lw=1, ls="--")
axA.text(-0.5, -0.44, "dependency threshold", fontsize=5.3, color=ALARM, ha="left", va="bottom")
axA.axhline(0, color="0.75", lw=.6)
axA.set_xticks(range(3))
axA.set_xticklabels([g[0] for g in groups], fontsize=6)
axA.set_ylabel("Mean Chronos gene effect\n(gene × lineage)")
axA.set_title("No HF-axis receptor is a\ndependency in any lineage", loc="left")
axA.set_xlim(-.55, 2.55)

axB.scatter(cor["spearman_rho"], -np.log10(cor["p_value"]), s=14, c=FOCAL, alpha=.65, linewidth=0)
axB.axhline(-np.log10(0.05), color="0.4", lw=.8, ls=":")
axB.text(-0.55, -np.log10(0.05) + .04, "nominal p = 0.05", fontsize=5.5, color="0.4")
axB.set_xlabel("Spearman ρ:\ngene effect vs compatibility", fontsize=6.5)
axB.set_ylabel("−log₁₀ p")
axB.set_title(f"0 of {len(cor)} receptors survive\nFDR correction", loc="left")
axB.axvline(0, color="0.8", lw=.6)

mix = pan[~pan["arm"].str.contains("subtype_matched")].sort_values("receptor_expr_log2tpm1", ascending=False)
for i, (a, lab, col) in enumerate([("responder_high_receptor", "high", FOCAL), ("control_low_receptor", "low", GREY)]):
    v = mix[mix["arm"] == a].reset_index(drop=True)
    for j, (cl, yv) in enumerate(zip(v["cell_line"], v["receptor_expr_log2tpm1"].values)):
        axC.scatter([i], [yv], s=34, color=col, zorder=3, edgecolor="white", linewidth=.5)
        axC.annotate(cl, (i, yv), xytext=(9, (j % 2) * 8 - 4), textcoords="offset points",
                     fontsize=5, va="center", color="0.35")
axC.set_xticks([0, 1])
axC.set_xticklabels(["high FZD7\n(responder)", "low FZD7\n(control)"], fontsize=6)
axC.set_ylabel("FZD7 expression, log₂(TPM+1)")
axC.set_title("Uterine panel: 48.7× separation,\nhousekeeping-independent", loc="left")
axC.set_xlim(-.4, 1.95)
axC.set_ylim(-0.6, 6.6)
for a, l in [(axA, "a"), (axB, "b"), (axC, "c")]:
    panel_letter(a, l)

for t in list(axA.texts):
    if t.get_text() in ("-2.39", "-0.32", "0.01", "0.00", "-0.00"):
        t.remove()
meds = {}
for i, (lab, gl, col) in enumerate(groups):
    v = dep[dep["gene"].isin(gl)]["mean_gene_effect"].dropna().values
    meds[lab] = float(np.median(v))
    axA.annotate(f"median {np.median(v):.2f}", (i + 0.34, np.median(v)), xytext=(0, 0),
                 textcoords="offset points", ha="left", va="center", fontsize=5.3, color=col)
axA.set_xlim(-.55, 3.05)
for t in list(axB.texts):
    t.remove()
axB.text(0.02, 0.62, "nominal p = 0.05", transform=axB.transAxes, fontsize=5.5, color="0.4")

panel_letter(axB, "b")
axA.set_xticklabels(["Essential\ncontrols", "Dependency\nreceptors", "HF-axis\nreceptors"], fontsize=5.6)
axA.tick_params(axis="x", pad=2)

fig.savefig(os.path.join(OUT, "fig4_depmap_validation.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig4_depmap_validation.png"), dpi=300, bbox_inches="tight")