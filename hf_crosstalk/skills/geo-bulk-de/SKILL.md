---
name: geo-bulk-de
description: Download a GEO bulk RNA-seq series (supplementary matrices, per-sample count CSVs, or series matrix metadata), parse disease-group labels out of the series matrix characteristics, and run two-group differential expression with empirical-Bayes moderated t statistics implemented in pure Python (a limma-trend / limma-voom equivalent, not a call into limma). Also computes cross-cohort replication statistics — overlap hypergeometric p, directional concordance, log2FC Spearman rho. Use for GEO series accessions (GSE…), MAGNet/LV heart cohorts, or any bulk RNA-seq case-vs-control contrast where group labels must be derived from GEO metadata rather than hardcoded.
---

# geo-bulk-de

Differential expression from a GEO series, end to end, with no R dependency and
no hardcoded sample labels.

## Purpose

Three jobs, in order:

1. **Fetch** — build GEO FTP paths from a GSE accession, download supplementary
   expression matrices (or a `_RAW.tar` of per-sample tables) and the series
   matrix.
2. **Label** — parse `!Sample_characteristics_ch1` out of the series matrix into
   a per-sample DataFrame. Group assignment comes from this table; you never
   hand-write which GSM is a case.
3. **Test** — moderated-t differential expression appropriate to the data scale
   (raw counts vs RPKM/FPKM/TPM), plus cross-cohort replication statistics.

Everything lives in `kernel.py`: pure stdlib + pandas/numpy/scipy/statsmodels,
functions taking paths/DataFrames and returning DataFrames/dicts. No agent, no
platform SDK, no network beyond `ftp.ncbi.nlm.nih.gov`.

## Statistics: what is actually implemented

This module **reimplements limma's empirical-Bayes moderation in Python**. It
does not call limma, and results will not be bit-identical to limma. The
estimators are the published ones (Smyth 2004, *Stat Appl Genet Mol Biol* 3:3;
Law et al. 2014, *Genome Biol* 15:R29), implemented as follows.

Per gene, ordinary least squares on the design `X = [intercept, group_indicator,
covariates…]` gives contrast estimate `beta`, residual variance `s2`, and
residual df `df = n - p`. Then:

**Prior estimation** (`fit_f_dist`, moment matching on `log s2`):

```
e     = log(s2) - digamma(df/2) + log(df/2)
evar  = var(e) - mean(trigamma(df/2))
d0    = 2 * trigamma^-1(evar)                        if evar > 0
s0^2  = exp(mean(e) + digamma(d0/2) - log(d0/2))
d0    = inf ,  s0^2 = exp(mean(e))                   if evar <= 0
```

`trigamma_inverse` is a Newton iteration on `polygamma(1, ·)`.

**Posterior variance and moderated t** (`moderated_t`):

```
s~^2 = (d0 * s0^2 + df * s2) / (d0 + df)
t    = beta / (s~ * sqrt(c' (X'X)^-1 c))
df_total = df + d0            (two-sided t; normal limit when d0 = inf)
```

**The trend** (`limma-trend`): with `trend=True`, `s0^2` is gene-specific.
`log(s2)` is lowess-smoothed against mean expression, and `fit_f_dist` is
applied to `s2 / trend` so that `s0^2_g = scale * trend(mean_expr_g)`. This is
limma's `eBayes(trend=TRUE)` behaviour.

**voom weights** (`_voom_weights`, used by `design="voom_like"`): lowess of
`sqrt(residual sd)` on mean log2 count, evaluated at each fitted value to give
per-observation precision weights `w = 1/sqrt_sd^4`, then weighted least squares
with a non-trended prior. This is the voom construction, but the lowess span and
the fitted-value interpolation differ in detail from limma's `voom()`.

**FDR** is Benjamini–Hochberg (`bh_fdr`).

### Which `design` to use

| data you have | `design` | transform applied |
|---|---|---|
| raw counts | `"voom_like"` | log2-CPM + voom precision weights, non-trended prior |
| raw counts | `"limma_trend_counts"` | log2-CPM, trended prior |
| RPKM / FPKM / TPM | `"limma_trend"` | `log2(x + 1)`, trended prior |
| already log-scale | `"limma_trend"` + `already_logged=True` | none |

RPKM/FPKM input is a real limitation, not a preference: gene-length-normalised
values cannot be given a count-based mean-variance model, library-size
normalisation has already happened upstream in a way you cannot inspect, and
composition bias (a few very high genes absorbing signal) is baked in. Say so
when you report results from `"limma_trend"` on RPKM.

## Function signatures

Fetching:

```python
gse_ftp_dir(gse) -> str
list_geo_supplementary(gse, timeout=300) -> list[str]
fetch_geo_supplementary(gse, dest_dir, pattern=None, timeout=300,
                        overwrite=False) -> list[str]     # local paths
fetch_series_matrix(gse, dest_dir, timeout=300, overwrite=False) -> str
extract_tar(tar_path, dest_dir, pattern=None) -> list[str]
```

Metadata and loading:

```python
parse_series_matrix_metadata(path) -> DataFrame   # indexed by GSM
load_expression(path, sep="\t", gene_col=0, drop_cols=(), aggregate="sum")
    -> DataFrame                                  # genes x samples
load_per_sample_tables(paths, value_col=None, gene_col=0, sep=",",
                       sample_name_regex=r"(GSM\d+)") -> DataFrame
strip_ensembl_version(index) -> list[str]
```

Testing:

```python
filter_by_expression(expr, min_value=1.0, min_fraction=0.2, groups=None) -> DataFrame
cpm(counts, log=False, prior_count=0.5, lib_sizes=None) -> DataFrame
check_de_inputs(expr, groups, case, control, min_per_group=2) -> Series
fit_f_dist(s2, df_resid) -> {"d0", "s0_2"}
moderated_t(beta, s2, df_resid, stdev_unscaled, mean_expr, trend=True,
            trend_frac=0.4) -> dict
bh_fdr(p) -> ndarray

run_de(expr, groups, case, control, design="limma_trend", covariates=None,
       min_per_group=2, filter_min_value=None, filter_min_fraction=0.2,
       trend_frac=0.4, already_logged=False) -> DataFrame
```

`run_de` returns one row per tested gene, sorted by `p_value`, columns:
`gene, log2FC, t, p_value, fdr, mean_expr, n_case, n_ctrl`, with
`.attrs = {"design", "d0", "df_resid", "n_genes_tested"}`. `log2FC` is
case-minus-control on the log2 scale used by the chosen design. Samples whose
group label is neither `case` nor `control` are dropped, so a three-arm series is
handled by calling `run_de` once per contrast.

Replication:

```python
replication_stats(de_a, de_b, fdr=0.05, gene_col="gene") -> dict
replicated_signature(de_a, de_b, fdr=0.05, gene_col="gene",
                     suffixes=("discovery","validation")) -> DataFrame
```

`replication_stats` keys: `n_genes_a, n_genes_b, n_shared_genes, n_sig_a,
n_sig_b, n_replicated, n_replicated_same_sign, expected_overlap, hypergeom_p,
directional_concordance, spearman_rho, spearman_p, pearson_r_sig,
fdr_threshold`. `n_sig_a`/`n_sig_b` are counted **among shared genes** so the
hypergeometric test is internally consistent; `hypergeom_p` is
`P(overlap >= observed)` under independence with the shared-gene set as the
universe.

`replicated_signature` returns `gene, log2FC_<a>, fdr_<a>, log2FC_<b>,
fdr_<b>, direction` where `direction` is `"up"`/`"down"`, sorted by the worse of
the two FDRs.

## Worked example

GSE116250 (LV RNA-seq, RPKM; non-failing vs dilated vs ischemic cardiomyopathy),
pooling both cardiomyopathies into one failing arm:

```python
import pandas as pd
import kernel as k

paths = k.fetch_geo_supplementary("GSE116250", "data/GSE116250", pattern=r"rpkm")
sm    = k.fetch_series_matrix("GSE116250", "data/GSE116250")
meta  = k.parse_series_matrix_metadata(sm)
print(meta.columns.tolist())          # ALWAYS inspect before assuming a column name
# -> ['gsm','title','source_name_ch1','organism_ch1','library_strategy',
#     'tissue','disease','age','sex']            (the label column here is 'disease')
print(meta["disease"].value_counts())
# -> dilated cardiomyopathy 37 / non-failing 14 / ischemic cardiomyopathy 13

# this file carries a gene-symbol column alongside the ENSG ids -> drop it
expr = k.load_expression(paths[0], gene_col="Gene", drop_cols=["Common_name"])
expr.index = k.strip_ensembl_version(expr.index)

# QC: this submitter's RPKM file encodes a sentinel 999999999 in some cells.
# Drop affected genes rather than imputing them.
expr = expr.loc[~(expr == 999999999.0).any(axis=1)]

# meta is indexed by GSM but the RPKM columns are sample titles -> map explicitly
title2gsm = dict(zip(meta["title"], meta["gsm"]))
assert set(expr.columns) <= set(title2gsm)
DMAP = {"non-failing": "non_failing",
        "dilated cardiomyopathy": "DCM",
        "ischemic cardiomyopathy": "ICM"}
arm = pd.Series({c: DMAP[meta.loc[title2gsm[c], "disease"].strip().lower()]
                 for c in expr.columns})
groups = arm.replace({"DCM": "failing", "ICM": "failing"})

de = k.run_de(expr, groups, case="failing", control="non_failing",
              design="limma_trend", filter_min_value=1.0)
# secondary contrasts: pass `arm` instead — samples labelled neither
# case nor control are dropped automatically
de_dcm = k.run_de(expr, arm, case="DCM", control="non_failing",
                  design="limma_trend", filter_min_value=1.0)
print(de.head(), de.attrs)
de.to_csv("hf_de_GSE116250.csv", index=False)
```

Validation cohort with one per-sample CSV inside a `_RAW.tar`:

```python
tar   = k.fetch_geo_supplementary("GSE141910", "data/GSE141910", pattern=r"RAW\.tar$")[0]
files = k.extract_tar(tar, "data/GSE141910/csv", pattern=r"\.csv(\.gz)?$")
mat   = k.load_per_sample_tables(files)      # columns named from the GSM prefix
mat.index = k.strip_ensembl_version(mat.index)
```

**Never trust a `_RAW.tar` to contain raw counts — check the scale.** A GEO
series named "RAW" and described only as "STAR alignment" may still ship
normalised values, and picking the wrong `design` silently invalidates the
variance model:

```python
v = mat.iloc[:, 0]
print(bool((v == v.round()).all()),        # integers?
      float(v.min()), float(v.max()),      # any exact zeros? plausible ceiling?
      int((v == 0).sum()),
      float(mat.sum().std() / mat.sum().mean()))   # column-sum CV
```

Non-integer values, no zeros, a max around 20, and a column-sum CV near 0
(GSE141910: `False 1.82 21.03 0 0.0066`, verified by running the block above) mean the data are already log2-scale
and library-normalised — a variance-stabilising transform, not counts. Use the
already-logged path and say so when reporting:

```python
de2 = k.run_de(mat, groups2, "failing", "non_failing",
               design="limma_trend", already_logged=True, filter_min_value=None)
```

Integer values with exact zeros and a wide column-sum CV would instead mean real
counts — then use `design="voom_like", filter_min_value=1.0`.

Harmonise identifiers **before** comparing (`replication_stats` raises rather
than silently returning zero overlap), then:

```python
stats = k.replication_stats(de, de2, fdr=0.05)
sig   = k.replicated_signature(de, de2, fdr=0.05)
```

Read `n_shared_genes` against `n_sig_a`/`n_sig_b` before quoting `hypergeom_p`:
when a large fraction of shared genes is significant in both cohorts the
hypergeometric test saturates (a tiny p-value at ~1x fold enrichment), and
`directional_concordance` plus `spearman_rho` are the metrics that carry
information.

## Preconditions (enforced — these raise, they are not advice)

`check_de_inputs`, called by every `run_de`, raises `ValueError` when:

- the expression matrix is empty in either dimension;
- sample identifiers (columns) are not unique;
- any gene row is entirely NA;
- any sample column is non-numeric after coercion;
- the group vector is not aligned to `expr.columns` (any column without a
  label — the error names the offenders rather than silently dropping them);
- `case` or `control` is absent from the group vector (the error lists the labels
  that *are* present);
- either group has fewer than `min_per_group` (default 2) samples.

`run_de` additionally raises on an unknown `design`, a rank-deficient design
matrix, NA covariates after alignment, and an expression filter that removes
every gene. `cpm` raises on a non-positive library size. `parse_series_matrix_metadata`
raises on a missing `!Sample_geo_accession` line, ragged characteristics rows, or
duplicate GSM ids. `replication_stats` raises on missing columns, duplicate gene
ids, or zero shared genes (the "you forgot to harmonize identifiers" case).

## Postconditions

- `run_de` output has no duplicate genes, `fdr >= p_value` elementwise, and rows
  sorted ascending by `p_value`.
- `log2FC > 0` means higher in `case`.
- `parse_series_matrix_metadata` returns exactly one row per GSM in series order.
- `load_expression` / `load_per_sample_tables` return an all-numeric frame with
  unique gene index and no all-NA rows.
- `replicated_signature` rows satisfy `fdr < threshold` in both cohorts and
  matching `log2FC` signs.

## Sanity checks worth running

The module was validated on synthetic negative-binomial-like counts (4000 genes,
15 vs 15, 200 genes spiked 3x): all 200 spikes recovered at FDR<0.05 under both
counts designs, and under a permuted null the p-value distribution was uniform
(4.8% below 0.05) with zero genes at FDR<0.05. Re-run that pattern after any
change to `fit_f_dist` or `moderated_t` — a broken prior shows up as either a
degenerate `d0` (0 or inf) or an anti-conservative null.

## Standalone use (no agent)

`kernel.py` imports only `gzip`, `io`, `os`, `re`, `tarfile`,
`urllib.request`, `numpy`, `pandas`, `scipy`, `statsmodels`. There are no
references to any agent/platform API, so:

```bash
pip install pandas numpy scipy statsmodels
cd skills/geo-bulk-de
python - <<'PY'
import kernel as k
sm = k.fetch_series_matrix("GSE116250", "/tmp/gse")
print(k.parse_series_matrix_metadata(sm).head())
PY
```

Or copy `kernel.py` next to your script and `import kernel`. The only network
host required is `ftp.ncbi.nlm.nih.gov` over HTTPS; all functions that touch the
network take an explicit `dest_dir` and cache by filename, so a second run with
the same `dest_dir` is offline unless `overwrite=True`.
