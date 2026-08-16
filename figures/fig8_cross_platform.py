import os, re, json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# Load data
cp = pd.read_csv(os.path.join(TABLES, "hf_cross_platform_replication.csv"))
cm = pd.read_csv(os.path.join(TABLES, "hf_canonical_markers_cross_platform.csv"))
pd_ = pd.read_csv(os.path.join(TABLES, "hf_permutation_controls.csv"))
pdg = pd.read_csv(os.path.join(TABLES, "hf_marker_probe_level_diagnostics.csv"))

prim = cp[(cp.reference == "signature_mean_of_two_rnaseq") & (cp.contrast == "failing_vs_non_failing")]
prim = prim.set_index("gse")

pw = pd.read_csv(os.path.join(TABLES, "hf_pairwise_cohort_agreement.csv"))
# pairwise cohort agreement table - need same_technology column
# pw is the pairwise cohort agreement data; in the trace it was loaded separately
# but from the dependency mappings, hf_cross_platform_replication.csv is our source
# The pw table has same_technology, cohort_a, cohort_b columns
# Let's check if cross_platform has pairwise rows
pw = pd.read_csv(os.path.join(TABLES, "hf_pairwise_cohort_agreement.csv"))

COH = ["GSE5406", "GSE57345", "GSE1869"]
NP = {"GSE5406": "n=194/16", "GSE57345": "n=177/136", "GSE1869": "n=31/6"}

META_GREY = "#aaaaaa"
COHC = {"GSE5406": "#1f4e9c", "GSE57345": "#4f8ed6", "GSE1869": "#b8b8b8"}
MK = ["NPPA", "NPPB", "MYH6", "ATP2A2"]

fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.65))

conc = [float(prim.loc[c, "concordance_sig_in_both"]) for c in COH]
rho  = [float(prim.loc[c, "rho_sig_in_both"]) for c in COH]

# (a)
ax = axes[0]
x = np.arange(len(COH))
for xi, c in enumerate(COH):
    ax.vlines(xi - 0.13, 0, conc[xi], color=META_GREY, lw=0.8, zorder=1)
    ax.vlines(xi + 0.13, 0, rho[xi],  color=META_GREY, lw=0.8, zorder=1)
    ax.plot(xi - 0.13, conc[xi], "o", ms=6, color=COHC[c], zorder=3)
    ax.plot(xi + 0.13, rho[xi], "s", ms=5.5, mfc="white", mec=COHC[c], mew=1.3, zorder=3)
ax.axhline(0.75, color="#c0392b", ls="--", lw=0.9)
ax.axhline(0.40, color="#c0392b", ls=":", lw=0.9)
ax.text(2.46, 0.765, "0.75 required", color="#c0392b", fontsize=6, va="bottom", ha="right")
ax.text(2.46, 0.415, "0.40 required", color="#c0392b", fontsize=6, va="bottom", ha="right")
ax.plot([], [], "o", ms=5, color="#555555", label="concordance")
ax.plot([], [], "s", ms=5, mfc="white", mec="#555555", label="ρ, sig in both")
ax.set_xticks(x)
ax.set_xticklabels([c + "\n" + NP[c] for c in COH], fontsize=6)
ax.set_ylabel("agreement with signature")
ax.set_ylim(-0.12, 1.02)
ax.margins(x=0.12)
ax.set_title("Two of three cohorts clear\nboth pre-specified thresholds", fontsize=8)
ax.legend(frameon=False, fontsize=6, loc="lower left", handletextpad=0.4, borderpad=0.1)

# (b)
ax = axes[1]
rng = np.random.default_rng(7)
for i, mask in enumerate([pw.same_technology == True, pw.same_technology == False]):
    sub = pw.loc[mask]
    if {"cohort_a", "cohort_b"}.issubset(sub.columns):
        orig = sub.apply(lambda r: {r.cohort_a, r.cohort_b} == {"GSE116250", "GSE141910"}, axis=1)
    else:
        orig = pd.Series(False, index=sub.index)
    jx = np.full(len(sub), i) + rng.uniform(-.07, .07, len(sub))
    ax.plot(jx[~orig.values], sub.loc[~orig.values, "concordance_sig_in_both"], "o", ms=4.5,
            color=("#1f4e9c" if i == 0 else "#b8b8b8"), alpha=.8, zorder=2)
    if orig.any():
        ax.plot(jx[orig.values], sub.loc[orig.values, "concordance_sig_in_both"], "*",
                ms=11, color="#c0392b", zorder=4)
        oy = float(sub.loc[orig.values, "concordance_sig_in_both"].iloc[0])
        ax.annotate("original pair,\nboth GPL16791", xy=(float(jx[orig.values][0]), oy),
                    xytext=(0.42, 0.565), fontsize=6, color="#c0392b", ha="center",
                    arrowprops=dict(arrowstyle="-", lw=0.7, color="#c0392b"))
    m = sub.concordance_sig_in_both.median()
    ax.hlines(m, i - 0.20, i + 0.20, color="black", lw=1.7, zorder=3)
    ax.text(i + 0.24, m, "%.3f" % m, va="center", fontsize=6.5, zorder=4)
ax.set_xticks([0, 1])
ax.set_xticklabels(["same\ntechnology\n(4 pairs)", "cross\ntechnology\n(6 pairs)"], fontsize=6.5)
ax.set_ylabel("directional concordance")
ax.set_ylim(0.44, 1.0)
ax.margins(x=0.30)
ax.set_title("Sharing a platform inflates\nagreement by 0.17", fontsize=8)

# (c)
ax = axes[2]
w = 0.26
for j, c in enumerate(COH):
    vals = [float(cm[(cm.gse == c) & (cm.gene == g)].log2FC.iloc[0]) for g in MK]
    fdrs = [float(cm[(cm.gse == c) & (cm.gene == g)].fdr.iloc[0]) for g in MK]
    xs = np.arange(len(MK)) + (j - 1) * w
    ax.bar(xs, vals, width=w * 0.86, color=COHC[c], edgecolor="black", lw=0.4, zorder=2, label=c)
    for xi, v_, f_ in zip(xs, vals, fdrs):
        if f_ < 0.05:
            ax.plot(xi, v_ + (0.16 if v_ >= 0 else -0.26), "*", ms=4.2, color="black", zorder=3)
ax.axhline(0, color="black", lw=0.7)
ax.set_xticks(np.arange(len(MK)))
ax.set_xticklabels(MK, fontsize=6.5, style="italic")
ax.set_ylabel("log$_2$ fold change")
ax.set_title("NPPB does not reproduce on\neither array platform", fontsize=8)
ax.set_ylim(-2.75, 2.75)
ax.legend(frameon=False, fontsize=6, loc="upper right", handlelength=1.0, borderpad=0.1)
ax.text(0.02, 0.03, "★ FDR<0.05", transform=ax.transAxes, fontsize=6, ha="left", va="bottom")

for ax in axes:
    ax.yaxis.set_label_coords(-0.20, 0.5)

fig.tight_layout(w_pad=2.2)
for _ax, _l in zip(axes, ["a", "b", "c"]):
    panel_letter(_ax, _l)

fig.savefig(os.path.join(OUT, "fig8_cross_platform.pdf"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig8_cross_platform.png"), dpi=300, bbox_inches="tight")
print("saved fig8_cross_platform.png")