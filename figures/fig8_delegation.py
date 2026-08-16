import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

apply_figure_style(sizes=(8, 7, 6))
META_GREY = "#9A9A9A"
COL = {"completed": "#4C6EA8", "failed": "#B5433A", "cancelled": "#D9A03C"}

dr = pd.read_csv(os.path.join(TABLES, "delegation_record.csv"))
T = dr.sort_values(["wave", "duration_min"]).reset_index(drop=True)

fig = plt.figure(figsize=(7.4, 3.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.5, 1], wspace=0.30)
axA = fig.add_subplot(gs[0])
axB = fig.add_subplot(gs[1])

# (a) per-track wall-clock duration, grouped by dispatch wave
y = np.arange(len(T))
for k, r in T.iterrows():
    axA.plot([T.duration_min.min() * 0.55, r.duration_min], [k, k],
             color=META_GREY, lw=0.7, zorder=1)
    axA.scatter(r.duration_min, k, s=30, color=COL[r.status], zorder=3)
axA.set_yticks(y)
axA.set_yticklabels(T.track, fontsize=6)
axA.set_xscale("log")
axA.set_xlim(T.duration_min.min() * 0.55, T.duration_min.max() * 2.4)
axA.set_xticks([10, 100, 1000])
axA.set_xticklabels(["10 min", "100 min", "1000 min"])
axA.set_xlabel("Wall-clock duration per track (log)", labelpad=5)
axA.set_title("Wall-clock duration of 23 sub-agent tracks\nacross 5 dispatch waves", loc="left")
axA.set_ylim(-0.8, len(T) - 0.2)
for w in sorted(T.wave.unique())[:-1]:
    cut = T[T.wave == w].index.max() + 0.5
    axA.axhline(cut, color=META_GREY, lw=0.5, ls=":")
for w in sorted(T.wave.unique()):
    idx = T[T.wave == w].index
    axA.text(1.015, (idx.min() + idx.max()) / 2, "wave %d" % w,
             transform=axA.get_yaxis_transform(), va="center", ha="left",
             fontsize=6, color=META_GREY)
axA.legend(handles=[Line2D([], [], marker="o", ls="", color=COL[s], markersize=5, label=s)
                    for s in ["completed", "failed", "cancelled"]],
           loc="center right", frameon=False, fontsize=6,
           bbox_to_anchor=(0.98, 0.30))
panel_letter(axA, "a")

# (b) token consumption, input against output
tot_in = T.input_tokens.sum() / 1e6
tot_out = T.output_tokens.sum() / 1e6
bars = axB.bar(["input\n(context)", "output\n(generated)"], [tot_in, tot_out],
               color=["#4C6EA8", "#8A8A8A"], width=0.55)
for b, v in zip(bars, [tot_in, tot_out]):
    axB.text(b.get_x() + b.get_width() / 2, v + tot_in * 0.02,
             "%.0fM" % v if v >= 10 else "%.1fM" % v,
             ha="center", va="bottom", fontsize=7, fontweight="bold")
axB.set_ylabel("Tokens (millions)")
axB.set_ylim(0, tot_in * 1.16)
axB.set_title("Token consumption by category\n(%.0f:1 input-to-output ratio)" % (tot_in / tot_out), loc="left")
axB.text(0.97, 0.62, "$%.0f total\n%d tracks\nmedian %.0f min"
         % (T.cost_usd.sum(), len(T), T.duration_min.median()),
         transform=axB.transAxes, ha="right", va="center", fontsize=7, color="#3A3A3A")
panel_letter(axB, "b")

fig.savefig(os.path.join(OUT, "fig8_delegation.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(OUT, "fig8_delegation.png"), dpi=300, bbox_inches="tight")
print("ratio %.0f:1 | total $%.0f | median %.0f min | in %.1fM out %.1fM"
      % (tot_in / tot_out, T.cost_usd.sum(), T.duration_min.median(), tot_in, tot_out))
