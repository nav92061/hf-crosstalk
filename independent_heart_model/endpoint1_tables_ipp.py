import hashlib
import os
import re
import secrets

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


import pandas as pd
import numpy as np

# Load manifest
man = pd.read_csv("colab_series_manifest.csv")

# Load the full volumes and pairs data
full = pd.read_csv("volumes_full.csv")
pp = pd.read_csv("recon_pairs_confirmed.csv")
sa = pp[pp.same_acquisition].copy()
sa["kernel_arm"] = sa.k_a.astype(str) != sa.k_b.astype(str)
sa["thickness_arm"] = sa.th_a != sa.th_b

STR = ["heart", "aorta", "pulmonary_vein", "superior_vena_cava", "inferior_vena_cava"]

vol = full[~full.failed].set_index("series_uid") if "failed" in full.columns else full.set_index("series_uid")

# Build failed flag if needed
if "failed" not in full.columns:
    full["failed"] = full.error.notna() & full.error.astype(str).str.strip().ne("")
    vol = full[~full.failed].set_index("series_uid")

good = set(vol.index)
sa["complete"] = sa.ua.isin(good) & sa.ub.isin(good)

# Build paired table
rows = []
for _, r in sa[sa.complete].iterrows():
    a, b = vol.loc[r.ua], vol.loc[r.ub]
    d = {"patient": r.patient, "ua": r.ua, "ub": r.ub,
         "kernel_arm": r.kernel_arm, "thickness_arm": r.thickness_arm,
         "k_a": r.k_a, "k_b": r.k_b, "th_a": r.th_a, "th_b": r.th_b}
    for s in STR:
        va, vb = a[f"{s}_vol_mm3"], b[f"{s}_vol_mm3"]
        d[f"{s}_a"], d[f"{s}_b"] = va, vb
        d[f"{s}_pctdiff"] = np.nan if (va == 0 or vb == 0) else 100 * abs(va - vb) / ((va + vb) / 2)
        d[f"{s}_either_zero"] = (va == 0) or (vb == 0)
    rows.append(d)
P = pd.DataFrame(rows)

# Compute z-coverage difference
ni = man.set_index("series_uid").n_images
th = man.set_index("series_uid").slice_thickness_mm
geo = pd.read_csv("series_geometry_ipp.csv").set_index("series_uid").extent_ipp
P["zmm_a"] = P.ua.map(geo)
P["zmm_b"] = P.ub.map(geo)
P["z_pctdiff"] = 100 * abs(P.zmm_a - P.zmm_b) / ((P.zmm_a + P.zmm_b) / 2)
P.to_csv("endpoint1_pairs.csv", index=False)

import warnings
import statsmodels.formula.api as smf

# CLEANEST stratum: same acquisition, z-matched, kernel differs, thickness IDENTICAL
C = P[(P.z_pctdiff <= 5) & P.kernel_arm & (~P.thickness_arm)]
print(f"CLEAN STRATUM: {len(C)} pairs, {C.patient.nunique()} patients")
print("  (same acquisition, z-coverage within 5%, kernel differs, slice thickness identical)")
rng = np.random.default_rng(0)


def bootci(x, B=10000):
    x = np.asarray(x, float)
    bs = np.median(rng.choice(x, size=(B, len(x)), replace=True), axis=1)
    return np.percentile(bs, [2.5, 97.5])


rows = []
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for s in STR:
        pp_ = C.groupby("patient")[f"{s}_pctdiff"].median().dropna()
        d = C[["patient", f"{s}_pctdiff"]].dropna().rename(columns={f"{s}_pctdiff": "y"})
        mm = smf.mixedlm("y ~ 1", d, groups=d.patient).fit(reml=True)
        ic = mm.conf_int().loc["Intercept"]
        lo, hi = bootci(pp_)
        rows.append(dict(structure=s, n_patients=len(pp_), n_pairs=len(d),
                         median_pct=round(pp_.median(), 2), med_ci=f"{lo:.2f}-{hi:.2f}",
                         max_pct=round(pp_.max(), 2),
                         mm_mean=round(float(mm.params["Intercept"]), 2),
                         mm_ci=f"{float(ic[0]):.2f}-{float(ic[1]):.2f}",
                         both_exclude_10=bool(hi < 10 and float(ic[1]) < 10)))
CL = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print(CL.to_string(index=False))
CL.to_csv("endpoint1_clean_stratum.csv", index=False)
print("\nall five structures: both aggregations exclude the 10% threshold ->",
      bool(CL.both_exclude_10.all()))