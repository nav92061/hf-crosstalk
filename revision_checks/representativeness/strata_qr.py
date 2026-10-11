# Analysis B (PRESPEC_representativeness.md): principal no-cutoff model (M4) within strata.
import pandas as pd, numpy as np, statsmodels.formula.api as smf, warnings, json, sys
warnings.filterwarnings("ignore")
P = pd.read_csv(sys.argv[1]); COV = pd.read_csv("pair_covariates_186.csv")
assert (P.ua.values == COV.ua.values).all() and (P.ub.values == COV.ub.values).all()
for c in ["thick_cat", "collection"]: P[c] = COV[c].values
STRUCTS = ["heart","aorta","pulmonary_vein","superior_vena_cava","inferior_vena_cava"]
P['th_same'] = np.isclose(P.th_a.astype(float), P.th_b.astype(float))
P['absz'] = P.z_pctdiff.abs()
def build(df, s):
    return pd.DataFrame({'patient': df.patient.values, 'y': df[f"{s}_pctdiff"].abs().values,
                         'absz': df.absz.values}).dropna(subset=['y','absz']).reset_index(drop=True)
def qr_cluster_boot(d, formula="y ~ absz", nboot=2000, seed=2026, q=0.5):
    base = smf.quantreg(formula, d).fit(q=q); terms = list(base.params.index)
    rng = np.random.default_rng(seed); pats = d.patient.unique()
    idx_by_pat = {pt: d.index[d.patient == pt].values for pt in pats}
    draws = []; nfail = 0
    for _ in range(nboot):
        sel = rng.choice(pats, len(pats), replace=True)
        bd = d.loc[np.concatenate([idx_by_pat[pt] for pt in sel])]
        try:
            f = smf.quantreg(formula, bd).fit(q=q); draws.append([f.params.get(t, np.nan) for t in terms])
        except Exception:
            draws.append([np.nan]*len(terms)); nfail += 1
    D = np.array(draws, float)
    return base, terms, np.nanpercentile(D, 2.5, axis=0), np.nanpercentile(D, 97.5, axis=0), nfail
rows = []
def fit(df, s, dataset, stratum_var, stratum):
    d = build(df, s); npair, npat = len(d), d.patient.nunique()
    r = dict(structure=s, dataset=dataset, stratum_var=stratum_var, stratum=stratum, n_pairs=npair, n_patients=npat,
             n_primary_pairs=int(df.primary.sum()) if 'primary' in df else np.nan)
    if npair < 8 or npat < 5:
        r.update(fitted=False, note="below pre-specified minimum (>=8 pairs from >=5 patients)"); rows.append(r); return
    base, terms, lo, hi, nf = qr_cluster_boot(d)
    for t, l, h in zip(terms, lo, hi):
        k = 'intercept' if t == 'Intercept' else t
        r[f'{k}_est'] = float(base.params[t]); r[f'{k}_ci_lo'] = float(l); r[f'{k}_ci_hi'] = float(h)
    r.update(fitted=True, n_boot_failed=nf, absz_range=f"{d.absz.min():.1f}-{d.absz.max():.1f}",
             frac_absz_le5=round(float((d.absz <= 5).mean()), 3), note="")
    rows.append(r); print(dataset, stratum_var, stratum, s, npair, npat, round(r['intercept_est'], 3), flush=True)
COV_prim = COV.primary.values
P['primary'] = COV_prim
# 1) reproduce published principal model (all pairs) and kernel-only thickness-identical fit
for s in STRUCTS:
    fit(P, s, "all_186_pairs", "all", "all")
K = P[P.kernel_arm & P.th_same].copy()
for s in STRUCTS:
    fit(K, s, "kernel_only_thickness_identical", "all", "all")
for var, levels in [("thick_cat", ["<=1.25","1.25-3",">3"]), ("collection", ["LUAD","LUSC","ESCA"])]:
    for lev in levels:
        sub = K[K[var] == lev]
        for s in STRUCTS:
            fit(sub, s, "kernel_only_thickness_identical", var, lev)
pd.DataFrame(rows).to_csv("strata_qr_raw.csv", index=False)
