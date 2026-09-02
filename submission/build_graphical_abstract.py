import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

try:
    apply_figure_style(sizes=(9,8,7))
except NameError:
    mpl.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False})

CL = pd.read_csv("endpoint1_clean_stratum.csv").set_index("structure")
obs = pd.read_csv("endpoint1_pairs.csv")

obs_df = pd.read_csv("endpoint2_discordance_distribution.csv")
nul_df = pd.read_csv("endpoint2_null_discordance.csv")

obs = obs_df.pct_discordant.dropna()
nul = nul_df.null_pct_discordant

fig = plt.figure(figsize=(13/2.54, 5.2/2.54), dpi=300)
gs = fig.add_gridspec(1, 3, left=.085, right=.985, top=.72, bottom=.30, wspace=.78)
BLUE, RED, GREY = "#2166AC", "#B2182B", "#7F7F7F"

# a — IMAGING: reconstruction effect sizes (the paper's primary result)
ax = fig.add_subplot(gs[0])
order = ["heart","aorta","superior_vena_cava","inferior_vena_cava","pulmonary_vein"]
lbl = {"heart":"heart","aorta":"aorta","superior_vena_cava":"SVC",
       "inferior_vena_cava":"IVC","pulmonary_vein":"pulm.\nvein"}
y = np.arange(len(order))
med = [CL.loc[s,"median_pct"] for s in order]
lo  = [float(CL.loc[s,"med_ci"].split("-")[0]) for s in order]
hi  = [float(CL.loc[s,"med_ci"].split("-")[1]) for s in order]
ax.errorbar(med, y, xerr=[np.array(med)-np.array(lo), np.array(hi)-np.array(med)],
            fmt="o", ms=4, color=BLUE, ecolor=BLUE, elinewidth=1.1, capsize=2, lw=0)
ax.axvline(10, color=RED, ls="--", lw=1)
ax.set_xlim(0, 13.5)
ax.text(10.4, 0.15, "10% limit", color=RED, fontsize=5.2, va="center", ha="left")
ax.set_yticks(y); ax.set_yticklabels([lbl[s] for s in order], fontsize=6)
ax.invert_yaxis()
ax.set_xlabel("volume change (%)", fontsize=6.5); ax.tick_params(labelsize=6)
ax.set_title("Same scan, two reconstructions\n45 pairs, 31 patients", fontsize=6.8, loc="left")

# b — IMAGING: the association endpoint is a null
ax = fig.add_subplot(gs[1])
bp = ax.boxplot([obs, nul], widths=.55, patch_artist=True, showfliers=False,
                medianprops=dict(lw=1.4), tick_labels=["observed","null"])
for p,c in zip(bp["boxes"], ["#9AB6D4", "#D9D9D9"]): p.set_facecolor(c); p.set_edgecolor("none")
for m,c in zip(bp["medians"], [BLUE, GREY]): m.set_color(c)
ax.set_ylabel("associations discordant (%)", fontsize=6.5)
ax.tick_params(labelsize=6); ax.set_ylim(0, 60)
ax.set_title("Association endpoint\nnull result, $p$ = 0.33", fontsize=6.8, loc="left")

# c — TRANSCRIPTOMIC: the same guardrail withdrawing a finding
ax = fig.add_subplot(gs[2])
steps = ["start","floor","source","adversarial"]; vals = [188,160,7,0]
ax.step(range(len(steps)), vals, where="post", color=BLUE, lw=1.6)
ax.fill_between(range(len(steps)), vals, step="post", color=BLUE, alpha=.15)
ax.annotate("7 axes withdrawn\nby pre-set criteria", xy=(2.55,7), xytext=(0.55,105),
            fontsize=5.4, color=RED, ha="left",
            arrowprops=dict(arrowstyle="->", lw=.7, color=RED))
ax.set_xticks(range(4))
ax.set_xticklabels(["start","floor","source","adv."], fontsize=5.4, rotation=0)
ax.set_ylabel("ligand-receptor axes", fontsize=6.5); ax.tick_params(labelsize=6)
ax.set_xlim(-0.4, 3.4)
ax.set_title("Same rule, second modality\n30 tumor types", fontsize=6.8, loc="left")

fig.text(.085, .90, "Executable assertions that halt, not caveats that do not",
         fontsize=8.4, fontweight="bold")
fig.text(.085, .055,
 "An assumption we wrote as prose, that reconstructions of one acquisition share coverage, was false; only measurement caught it.",
 fontsize=5.2, color="#333")
fig.savefig("Graphical Abstract.png", dpi=300, facecolor="white")
fig.savefig("Graphical Abstract.pdf", facecolor="white")