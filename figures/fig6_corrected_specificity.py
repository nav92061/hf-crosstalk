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
import shutil
import os

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

def A(fn):
    return pd.read_csv(os.path.join(TABLES, fn))

net = pd.read_csv(os.path.join(TABLES, "crosstalk_network.csv"))
ci = A(os.path.join(TABLES, "complex_independence.csv"))
fs = A(os.path.join(TABLES, "full_secretome_screen.csv"))

tb = lambda s: s.astype(str).str.upper().isin(["TRUE","1"])
S = fs[fs.gene.isin(sorted(set(net.ligand_subunit_gene.dropna())))].copy()
vmin = S[["log2FC_disc","log2FC_val"]].min(axis=1)
m_sig = (S.fdr_disc<0.05)&(S.fdr_val<0.05)
m_ind = (vmin>=1.0)&m_sig
m_ab  = tb(S.abundant)
m_base = S.best_heart_rank <= 26

steps = [len(S), int(m_ind.sum()), int((m_ind&m_ab).sum()), int((m_ind&m_ab&m_base).sum())]

fails = pd.DataFrame({"induction": (~m_ind).values, "abundance": (~m_ab).values,
                      "baseline": (~m_base).values}, index=S.index)
combo = fails.apply(lambda rr: tuple(sorted([k for k in fails.columns if rr[k]])), axis=1)
cc = combo.value_counts()

apply_figure_style(sizes=(8,7,6))
RED, GY = "#B5433A", "#9A9A9A"
fig = plt.figure(figsize=(7.4, 5.9))
gs = fig.add_gridspec(2, 2, hspace=0.60, wspace=0.36)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1])
axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])
slabs = ["all 51", "≥2-fold\ninduction", "+ abundant", "+ heart rank\n≤ 26 of 52"]
axA.bar(range(4), steps, color=[GY,GY,GY,RED], width=0.62)
for i,vv in enumerate(steps):
    axA.annotate(str(vv), (i,vv), xytext=(0,3), textcoords="offset points", ha="center", fontsize=7)
axA.set_xticks(range(4)); axA.set_xticklabels(slabs); axA.set_ylabel("Ligands")
axA.set_title("Sequential attrition across the three\npre-specified criteria", loc="left")
axA.set_ylim(0, 58)
LAB = {("abundance","baseline","induction"):"all three", ("abundance","induction"):"induction\n+ abundance",
       ("induction",):"induction only", ("abundance",):"abundance only", ("baseline",):"baseline only",
       ("baseline","induction"):"induction\n+ baseline", ("abundance","baseline"):"abundance\n+ baseline",
       ():"passes all"}
order = sorted(cc.index, key=lambda k: -cc[k])
axB.barh(range(len(order)), [cc[k] for k in order],
         color=[RED if k==() else GY for k in order], height=0.60)
for i,k in enumerate(order):
    axB.annotate(str(cc[k]), (cc[k], i), xytext=(3,0), textcoords="offset points", va="center", fontsize=7)
axB.set_yticks(range(len(order))); axB.set_yticklabels([LAB[k] for k in order]); axB.invert_yaxis()
axB.set_xlabel("Ligands"); axB.set_xlim(0, max(cc)*1.20)
axB.set_title("Criteria failed, by combination", loc="left")
sel = (m_ind & m_ab).values
axC.scatter(S.best_heart_rank[~sel], vmin[~sel], s=13, color=GY, zorder=2, label="did not clear")
axC.scatter(S.best_heart_rank[sel], vmin[sel], s=34, color=RED, zorder=3, label="cleared induction and abundance")
for g, dxy in [("NPPA",(6,2)),("NPPB",(6,-9)),("PTN",(6,2)),("MDK",(6,2))]:
    rw = S[S.gene==g].iloc[0]
    axC.annotate(g, (rw.best_heart_rank, vmin[S.gene==g].iloc[0]), xytext=dxy,
                 textcoords="offset points", fontsize=7, color=RED, style="italic")
axC.axhline(1.0, color="#555555", lw=0.7, ls="--")
axC.axvline(26, color="#555555", lw=0.7, ls=":")
axC.annotate("rank 26 cutoff", (26.8, 0.13), fontsize=6.3, color="#555555")
axC.set_xlabel("Best heart rank among 52 GTEx tissues")
axC.set_ylabel("Induction (min log$_2$FC)")
axC.set_title("Induction magnitude against\nbaseline cardiac rank", loc="left")
axC.legend(frameon=False, fontsize=6.5, loc="upper right")
ptn = ci[(ci.ligand=="PTN") & (ci.n_types_enriched_min_rule>0)].sort_values("n_types_enriched_min_rule", ascending=False)
yy = np.arange(len(ptn))
axD.hlines(yy, 0, ptn.n_types_enriched_min_rule, color=GY, lw=0.8)
axD.plot(ptn.n_types_enriched_min_rule, yy, "o", ms=5.4, color=RED)
for i,(_,rw) in enumerate(ptn.iterrows()):
    axD.annotate("%d" % rw.n_types_enriched_min_rule, (rw.n_types_enriched_min_rule, i),
                 xytext=(5,0), textcoords="offset points", va="center", fontsize=7)
axD.set_yticks(yy); axD.set_yticklabels(["PTN → %s" % x for x in ptn.complex]); axD.invert_yaxis()
axD.set_xlabel("Tumor types with receptor enriched")
axD.set_title("Surviving pleiotrophin axes\n(complex-level, minimum-subunit rule)", loc="left")
axD.set_xlim(0, ptn.n_types_enriched_min_rule.max()*1.25)
for ax,l in [(axA,"a"),(axB,"b"),(axC,"c"),(axD,"d")]: panel_letter(ax, l)
fig.savefig(os.path.join(OUT, "fig6_corrected_specificity.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig6_corrected_specificity.png"), dpi=300, bbox_inches="tight")