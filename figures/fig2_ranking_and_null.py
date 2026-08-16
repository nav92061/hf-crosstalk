import os, re, json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

rank = pd.read_csv(os.path.join(TABLES, "pancancer_susceptibility_ranking.csv"))

def apply_figure_style(sizes=(8, 7, 6)):
    mpl.rcParams.update({
        'font.size': sizes[1],
        'axes.titlesize': sizes[0],
        'axes.labelsize': sizes[1],
        'xtick.labelsize': sizes[2],
        'ytick.labelsize': sizes[2],
        'legend.fontsize': sizes[2],
        'axes.spines.top': False,
        'axes.spines.right': False,
    })

def panel_letter(ax, letter):
    ax.text(-0.08, 1.04, letter, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='bottom', ha='right')

FOC, GREY = "#4C6EA8", "#9A9A9A"

r = rank.sort_values("score", ascending=False).reset_index(drop=True)
sig = r.null_p_primary < 0.05
y = np.arange(len(r))[::-1]

vals = [0.323, 0.363, 0.107]
cols = ["#B5433A", "#B5433A", FOC]

apply_figure_style(sizes=(8, 7, 6))
fig = plt.figure(figsize=(7.4, 5.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.62, 1], wspace=0.34)
axA = fig.add_subplot(gs[0]); axB = fig.add_subplot(gs[1])
for k, (_, row) in enumerate(r.iterrows()):
    c = FOC if row.null_p_primary < 0.05 else GREY
    axA.plot([row.ci_low, row.ci_high], [y[k], y[k]], color=c, lw=1.4, solid_capstyle="butt", zorder=2)
    axA.plot(row.score, y[k], "o", ms=4.6, color=c, zorder=3)
axA.axvline(0, color="#555555", lw=0.7, ls=":", zorder=1)
axA.set_yticks(y); axA.set_yticklabels(r.tumor_type)
axA.set_xlabel("Compatibility score (availability-weighted receptor z)")
axA.set_title("Four of 30 tumor types exceed an\nexpression-matched random null", loc="left")
for k in [0, 1]:
    axA.annotate("%.2f" % r.score.iloc[k], (r.ci_high.iloc[k], y[k]), xytext=(4, 0),
                 textcoords="offset points", va="center",
                 color=FOC if sig.iloc[k] else GREY, fontsize=7)
axA.legend(handles=[Line2D([],[],color=FOC,marker="o",ms=4.6,lw=1.4,label="null p < 0.05 (receptor-specific)"),
                    Line2D([],[],color=GREY,marker="o",ms=4.6,lw=1.4,label="not distinguishable from random"),
                    Line2D([],[],color=GREY,lw=1.4,label="bootstrap 95% CI")],
           loc="lower right", frameon=False, fontsize=6.5, borderaxespad=0.2)
axA.margins(y=0.012); axA.set_xlim(-0.72, 1.02)
axB.bar(range(3), vals, color=cols, width=0.58)
for i, v in enumerate(vals):
    axB.annotate("%.3f" % v, (i, v), xytext=(0, 3), textcoords="offset points", ha="center", fontsize=7)
axB.axhline(0.5, color="#555555", lw=0.7, ls="--")
axB.annotate("confounded", (2.40, 0.508), ha="right", fontsize=6.5, color="#555555")
axB.annotate("lower = more receptor-specific", (-0.44, 0.578), fontsize=6.5, color="#555555")
axB.set_xticks(range(3))
axB.set_xticklabels(["Unmatched\npool", "Matched,\nno floor", "Matched,\nfloor"])
axB.set_ylabel("Mean Spearman ρ vs observed ranking")
axB.set_title("Effect of the expression floor on\nreceptor specificity", loc="left")
axB.set_ylim(0, 0.63)
for ax, l in [(axA,"a"), (axB,"b")]: panel_letter(ax, l)
fig.savefig(os.path.join(OUT, "fig2_ranking_and_null.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig2_ranking_and_null.png"), dpi=300, bbox_inches="tight")