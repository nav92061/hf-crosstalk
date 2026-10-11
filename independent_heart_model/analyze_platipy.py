"""Independent-model whole-heart analysis per PRESPEC_independent_heart_model.md.
Run from the repository root: python independent_heart_model/analyze_platipy.py
Writes independent_heart_model/independent_heart_model_results_replayed.csv."""
import os, re, sys, shutil, subprocess, numpy as np, pandas as pd
from scipy.stats import spearmanr
S5 = ["heart", "aorta", "pulmonary_vein", "superior_vena_cava", "inferior_vena_cava"]
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
IMG = os.path.join(ROOT, "imaging"); INP = HERE
VCSV = os.path.join(HERE, "volumes_platipy_heart.csv")

def replay(vol_df, tag, script, drop_pairs=None):
    d = os.path.join(HERE, "_replay", tag); shutil.rmtree(d, ignore_errors=True); os.makedirs(d)
    pp = pd.read_csv(f"{IMG}/recon_pairs_confirmed.csv")
    if drop_pairs is not None:
        key = set(map(frozenset, drop_pairs)); pp = pp[[frozenset((a, b)) not in key for a, b in zip(pp.ua, pp.ub)]]
    pp.to_csv(d + "/recon_pairs_confirmed.csv", index=False)
    shutil.copy(f"{IMG}/series_geometry_ipp.csv", d); shutil.copy(f"{IMG}/colab_series_manifest.csv", d)
    shutil.copy(f"{IMG if script == 'endpoint1_tables.py' else HERE}/{script}", d + "/endpoint1_tables.py")
    vol_df.to_csv(d + "/volumes_full.csv", index=False)
    r = subprocess.run([sys.executable, "endpoint1_tables.py"], cwd=d, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr[-1500:]
    return pd.read_csv(d + "/endpoint1_clean_stratum.csv").set_index("structure"), pd.read_csv(d + "/endpoint1_pairs.csv")

def hi(x): return float(str(x).split("-")[-1])

ts = pd.read_csv(f"{IMG}/volumes_full.csv")
# ---- acceptance checks
uids87 = pd.read_csv(os.path.join(HERE, "series.csv")).series_uid
v = pd.read_csv(VCSV).drop_duplicates("series_uid", keep="last")
chk = pd.DataFrame({"series_uid": uids87}).merge(v, on="series_uid", how="left")
chk["missing_row"] = chk.error.isna() & chk.heart_vol_mm3.isna()
chk["run_error"] = chk.error.notna()
ok_run = ~chk.run_error & ~chk.missing_row
chk["extent_fail"] = ok_run & ~chk.extent_check_pass.fillna(False).astype(bool)
chk["label_fail"] = ok_run & ((chk.n_labels.fillna(0) != 1) | chk.heart_empty.fillna(True).astype(bool))
chk["failed"] = chk.missing_row | chk.run_error | chk.extent_fail | chk.label_fail
checks_summary = {"series_expected": len(uids87), "rows_returned": int((~chk.missing_row).sum()),
                  "run_errors": int(chk.run_error.sum()), "extent_fail": int(chk.extent_fail.sum()),
                  "label_or_empty_fail": int(chk.label_fail.sum()), "accepted": int((~chk.failed).sum()),
                  "heart_edge_contact": int(chk.loc[~chk.failed, "heart_edge"].astype(bool).sum())}

# ---- volume table: TotalSegmentator table with heart replaced by PlatiPy (other structures unused)
pl = chk.set_index("series_uid")
vv = ts.copy().set_index("series_uid")
sel = vv.index.isin(pl.index)
vv.loc[sel, "heart_vol_mm3"] = pl.loc[vv.index[sel], "heart_vol_mm3"]
bad = pl.index[pl.failed]
vv.loc[vv.index.isin(bad), "failed"] = True
vv.loc[vv.index.isin(bad), "error"] = "excluded by acceptance check"
vv = vv.reset_index()

# FOV-matched rule (same as manuscript: ReconstructionDiameter differs by >0.5 mm)
hdr = pd.read_csv(f"{IMG}/headers_merged.csv", low_memory=False).set_index("series_uid")
fov = pd.to_numeric(hdr.ReconstructionDiameter, errors="coerce")
allp = pd.read_csv(f"{IMG}/recon_pairs_confirmed.csv")
fov_drop = [tuple(x) for x in allp[(allp.ua.map(fov) - allp.ub.map(fov)).abs() > 0.5][["ua", "ub"]].values]

res = {}
for mdl, vol in [("TotalSegmentator", ts), ("PlatiPy", vv)]:
    for tag, script, drop in [("45-pair, image count x thickness (primary)", "endpoint1_tables.py", None),
                              ("47-pair, slice positions", "endpoint1_tables_ipp.py", None),
                              ("FOV-matched, image count x thickness", "endpoint1_tables.py", fov_drop)]:
        res[(mdl, tag)] = replay(vol, re.sub(r"\W+", "_", f"{mdl}_{tag}")[:40], script, drop)
rows = []
for (mdl, tag), (tb, _) in res.items():
    r = tb.loc["heart"]
    rows.append(dict(model=mdl, analysis=tag, n_pairs=int(r.n_pairs), n_patients=int(r.n_patients),
                     median_pct=r.median_pct, median_ci=r.med_ci, mm_mean=r.mm_mean, mm_ci=r.mm_ci,
                     both_below_10=(hi(r.med_ci) < 10) and (hi(r.mm_ci) < 10)))
heart_tab = pd.DataFrame(rows)

# ---- secondary: between-model agreement (accepted series)
acc = pl[~pl.failed]
a = acc.heart_vol_mm3; b = ts.set_index("series_uid").loc[a.index, "heart_vol_mm3"]
okv = (a > 0) & (b > 0); sd = 100 * (a[okv] - b[okv]) / ((a[okv] + b[okv]) / 2)
agree = {"n_series": int(okv.sum()), "median_signed_pct_platipy_minus_ts": round(float(sd.median()), 2),
         "iqr": f"{sd.quantile(.25):.1f} to {sd.quantile(.75):.1f}", "spearman_rho": round(float(spearmanr(a[okv], b[okv])[0]), 3)}

# ---- secondary: direction (sharper kernel larger), same ordering as Table 4
GE = {"SOFT": 1, "STANDARD": 2, "DETAIL": 3, "LUNG": 4, "BONEPLUS": 5, "BONE": 5}
PH = {"A": 1, "B": 2, "C": 3, "D": 4, "YA": 5, "YB": 6}
def sharp(k):
    k = str(k).strip()
    if k in GE: return GE[k]
    m = re.match(r"^B(\d+)", k)
    if m: return int(m.group(1)) / 10
    return PH.get(k)
def direction(pp):
    st = pp[(pp.z_pctdiff <= 5) & pp.kernel_arm & (~pp.thickness_arm)].copy()
    st["sa"] = st.k_a.map(sharp); st["sb"] = st.k_b.map(sharp); st = st.dropna(subset=["sa", "sb"]); st = st[st.sa != st.sb]
    va = st.ua.map(VOL); vb = st.ub.map(VOL)
    flip = st.sa > st.sb
    smooth = np.where(flip, vb, va); sharp_ = np.where(flip, va, vb)
    d = pd.Series(100 * (sharp_ - smooth) / ((sharp_ + smooth) / 2), index=st.index)
    d[(smooth == 0) | (sharp_ == 0) | pd.isna(smooth) | pd.isna(sharp_)] = np.nan
    pat = d.groupby(st.patient).median().dropna()
    return {"n_patients": len(pat), "n_sharper_larger": int((pat > 0).sum()), "median_sharp_minus_smooth_pct": round(float(pat.median()), 2)}
dirn = {}
for mdl, vol in [("TotalSegmentator", ts), ("PlatiPy", vv)]:
    VOL = vol.set_index("series_uid").heart_vol_mm3.where(~vol.set_index("series_uid").get("failed", pd.Series(False, index=vol.series_uid)).fillna(False).astype(bool))
    dirn[mdl] = direction(res[(mdl, "45-pair, image count x thickness (primary)")][1])

heart_tab.to_csv(os.path.join(HERE, "independent_heart_model_results_replayed.csv"), index=False)
print(checks_summary); print(heart_tab.to_string(index=False)); print("agreement", agree); print("direction", dirn)
