"""Pre-specified summary for PRESPEC_external_accuracy.md.
Median (IQR) Dice and HD95; Bland-Altman mean bias and 95% LoA for volume (TS - reference).
Failures/exclusions listed with reason, not dropped silently."""
import sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

raw = pd.read_csv(sys.argv[1])
raw.to_csv("external_accuracy_cases.csv", index=False)
ok = raw[raw.status == "ok"].copy()

def miqr(x):
    q = np.percentile(x, [25, 50, 75]); return q[1], q[0], q[2]

rows = []
for name, col in [("Dice", "dice"), ("HD95_mm", "hd95_mm"), ("HD95_TS_to_ref_mm", "hd95_ts_to_ref_mm"),
                  ("HD95_ref_to_TS_mm", "hd95_ref_to_ts_mm"), ("Volume_ref_mL", "vol_ref_ml"), ("Volume_TS_mL", "vol_ts_ml")]:
    m, lo, hi = miqr(ok[col].astype(float))
    rows.append(dict(metric=name, n=len(ok), median=m, q1=lo, q3=hi, mean=np.nan, sd=np.nan, loa_lower=np.nan, loa_upper=np.nan))
for name, col in [("Volume_diff_mL (TS-ref)", "vol_diff_ml"), ("Volume_diff_pct (TS-ref)/ref", "vol_diff_pct")]:
    d = ok[col].astype(float); mu, sd = d.mean(), d.std(ddof=1)
    m, lo, hi = miqr(d)
    rows.append(dict(metric=name, n=len(ok), median=m, q1=lo, q3=hi, mean=mu, sd=sd,
                     loa_lower=mu - 1.96 * sd, loa_upper=mu + 1.96 * sd))
summ = pd.DataFrame(rows)
summ["n_total_cases"] = len(raw)
summ["n_not_ok"] = int((raw.status != "ok").sum())
summ.to_csv("external_accuracy_summary.csv", index=False)
print(summ.round(3).to_string(index=False))
print(raw[raw.status != "ok"][["PatientID", "status", "reason"]].to_string(index=False))

# Bland-Altman (mL)
mean_v = (ok.vol_ref_ml + ok.vol_ts_ml) / 2; diff = ok.vol_diff_ml
mu, sd = diff.mean(), diff.std(ddof=1)
fig, ax = plt.subplots(figsize=(5.5, 4.2))
ax.scatter(mean_v, diff, s=22, color="#2b6a99", edgecolor="white", linewidth=0.5)
for y, ls, lab in [(mu, "-", f"Bias {mu:.0f} mL"), (mu + 1.96 * sd, "--", f"+1.96 SD {mu+1.96*sd:.0f} mL"),
                   (mu - 1.96 * sd, "--", f"−1.96 SD {mu-1.96*sd:.0f} mL")]:
    ax.axhline(y, color="#444444", ls=ls, lw=1)
    ax.text(ax.get_xlim()[1] if False else mean_v.max(), y, " " + lab, va="bottom", ha="right", fontsize=8)
ax.axhline(0, color="#999999", lw=0.6, ls=":")
ax.set_xlabel("Mean of TotalSegmentator and reference heart volume (mL)")
ax.set_ylabel("TotalSegmentator − reference (mL)")
ax.set_title(f"LCTSC heart volume agreement (n = {len(ok)})", fontsize=10)
for s in ("top", "right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig("fig_external_accuracy_bland_altman.png", dpi=200)
