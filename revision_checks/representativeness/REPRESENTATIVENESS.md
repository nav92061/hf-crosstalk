# Representativeness of the retained primary stratum: results

Analysis run 2026-10-10 (UTC), after and according to `PRESPEC_representativeness.md` (written 2026-10-10 21:01 UTC). The manuscript was not edited.

## 0. Sanity check: reproduction of published Table 2
`endpoint1_tables.py` was run unchanged except that the absolute path of `colab_series_manifest.csv` was replaced by a local copy (inputs: `recon_pairs_confirmed.csv`, `volumes_full.csv`, `colab_series_manifest.csv`). The output `endpoint1_clean_stratum.csv` is **byte-identical** to published Table 2 (45 pairs, 31 patients), and the regenerated `endpoint1_pairs.csv` is exactly equal to the supplied 186-pair file.

## 1. Analysis A: raking of the 45 retained pairs to the 186-pair margins

### Methods as run
- **Population**: the 45 primary-stratum pairs (kernel differs, thickness identical, z-coverage difference ≤5%) from 31 patients.
- **Target**: the 186 segmented pairs (`endpoint1_pairs.csv`). Covariates were derived per pair from `headers_merged.csv` (joined on series UID) and the series manifest (`pair_covariates_186.csv`).
  - *Slice-thickness category* (≤1.25, >1.25–3, >3 mm) from **series A** thickness (`th_a`, identical to the manifest `slice_thickness_mm`). See Deviation 1: 102/186 pairs have different A and B thicknesses.
  - *Collection*: `collection` (TCGA-LUAD/LUSC/ESCA). No pair fell in "other". Identical for both series of every pair.
  - *Vendor*: `vendor` (GE/Siemens/Philips). No "other". Identical within every pair.
  - *Contrast*: DICOM **`ContrastBolusAgent` (0018,0010)**. A non-empty value was coded "agent recorded", a literal NONE/NO/N/A "none recorded", and empty/missing "unknown" (its own level, per plan). This coding coincides exactly with the derived `contrast` flag in `headers_merged.csv`, and the status was concordant within every pair (6 pairs differ in the agent string, but both members record an agent). No pair had an explicit "none" value, so the levels present were agent recorded (37/186) and unknown (149/186).
- **Raking**: iterative proportional fitting over the four margins (unit starting weights, up to 1,000 cycles, converged when the maximum weight change < 1e-10). Weights were normalised to mean 1 and then capped at 5× the mean (single pass, no re-raking).
- **Outcome**: per patient, the median PD of its pairs, with patient weight = sum of its pair weights; then the weighted median across patients (smallest value with cumulative weight ≥ half, midpoint on an exact half). 95% CI = 2.5/97.5 percentiles of 10,000 patient-level bootstrap resamples (NumPy `default_rng(0)`), with raking and trimming re-run in every resample; duplicated patients are treated as distinct clusters.
- **Effective sample size**: Kish, n_eff = (Σw)²/Σw².

### Margins
| Variable | Level | 186 pairs n (share) | Retained n (share) | Share after raking | Share after 5× trim |
|---|---|---|---|---|---|
| thick_cat | <=1.25 | 28 (0.150) | 1 (0.022) | 0.150 | 0.170 |
| thick_cat | 1.25-3 | 90 (0.484) | 11 (0.244) | 0.484 | 0.412 |
| thick_cat | >3 | 68 (0.366) | 33 (0.733) | 0.366 | 0.417 |
| collection | LUAD | 37 (0.199) | 15 (0.333) | 0.261 | 0.339 |
| collection | LUSC | 76 (0.409) | 29 (0.644) | 0.535 | 0.491 |
| collection | ESCA | 73 (0.393) | 1 (0.022) | 0.204 | 0.170 |
| vendor | GE | 142 (0.763) | 30 (0.667) | 0.763 | 0.780 |
| vendor | Siemens | 38 (0.204) | 11 (0.244) | 0.204 | 0.170 |
| vendor | Philips | 6 (0.032) | 4 (0.089) | 0.032 | 0.049 |
| contrast_status | agent_recorded | 37 (0.199) | 16 (0.356) | 0.199 | 0.244 |
| contrast_status | unknown | 149 (0.801) | 29 (0.644) | 0.801 | 0.756 |

Every target level is represented in the retained stratum at the point estimate, so **no level is formally uncovered**. Coverage of three levels is nominal. ESCA (39% of the target) and ≤1.25 mm (15%) are each represented by one retained pair (TCGA-VR-A8ER and TCGA-50-6595, respectively), and Philips (3%) by 4 pairs from 2 patients. The only ESCA pair is a Siemens pair and the Siemens target is 20.4%, so the four margins cannot all be met. Raking converged to a fixed point with the collection margin unattained (ESCA 0.204 vs target 0.393; maximum margin gap 0.188, rising to 0.222 after trimming). Four pairs exceeded 5× the mean and were trimmed. Ten pairs (8 patients, mostly Siemens non-ESCA) received weights < 0.001.

**Effective sample size**: 8.23 pairs (of 45) and 4.90 patients (of 31). In 6,511 of 10,000 bootstrap resamples, at least one target level had no retained pair (ESCA 3,568; ≤1.25 mm 3,515; Philips 1,281; 1.25–3 mm 1). In those resamples, the target was renormalised over the levels present (Deviation 5). Raking converged in 9,932 of 10,000 resamples.

### Results (PD, %)
| Structure | Published Table 2 median (95% CI) | Plan estimator, uniform pair weights | **Raked weighted median (95% CI)** | Upper limit < 10%? |
|---|---|---|---|---|
| heart | 0.63 (0.50-0.97) | 0.82 | **1.28** (0.55–1.39) | yes |
| aorta | 0.98 (0.69-1.17) | 0.98 | **1.16** (0.47–1.56) | yes |
| pulmonary_vein | 4.96 (3.20-5.63) | 5.03 | **5.04** (2.81–8.64) | yes |
| superior_vena_cava | 2.12 (1.90-3.46) | 2.01 | **0.81** (0.60–2.82) | yes |
| inferior_vena_cava | 1.52 (1.15-2.22) | 1.69 | **3.42** (1.10–7.56) | yes |

The "uniform pair weights" column applies the plan's estimator with all pair weights = 1. It differs from Table 2 because the plan weights each patient by its number of pairs, whereas Table 2 weights patients equally. The difference between this column and the raked column is the effect of the reweighting itself.

**Pre-specified criterion (every weighted upper limit < 10%): met.** The largest upper limits were pulmonary vein 8.64% and IVC 7.56%. Some weighted point estimates moved appreciably (IVC 1.52 → 3.42%; SVC 2.12 → 0.81%; heart 0.63 → 1.28%). The heart point estimate lies near the upper end of its percentile interval. With an effective sample of about five patients, the weighted estimates are dominated by a few patients: one LUSC 1.25–3 mm patient with two pairs, the single ESCA pair, and the single ≤1.25 mm pair. The analysis shows that reweighting to the full-cohort mix does not push any structure past 10%. It does not provide precise representative estimates, and the ESCA margin could not be reached.

## 2. Analysis B: principal no-cutoff model within strata

### Methods as run
The principal no-cutoff model (`M4_median_quantreg_absz`: median regression of |PD| on |ZD|, `statsmodels` QuantReg q = 0.5, patient-cluster bootstrap, 2,000 resamples, `default_rng(2026)` re-seeded for every fit) was re-implemented verbatim from the lineage code of `coverage_continuous_model.csv` (`strata_qr.py`). Population: all segmented kernel-only, thickness-identical pairs (84 pairs, 54 patients). As in the published model, pairs with an either-member-zero PD were dropped per structure. Strata were slice-thickness category (unambiguous here because thickness is identical within each pair) and collection. A stratum was fitted only if it had ≥8 pairs from ≥5 patients.

**Reproduction check**: the all-186-pair fit and the 84-pair kernel-only fit reproduce every published M4 estimate and CI bound (intercept and slope, all five structures) with a maximum absolute difference of 0.0 and identical n.

### Intercepts (PD at ZD = 0, %; 95% CI)
| Stratum | Pairs/patients (max over structures) | Primary-stratum pairs in stratum | Heart | Aorta | Pulmonary vein | SVC | IVC |
|---|---|---|---|---|---|---|---|
| All 186 segmented pairs (published M4 reproduction) | 185/58 | 45 | 0.98 (0.55 to 1.19) | 0.97 (0.49 to 1.34) | 4.94 (3.89 to 5.83) | 2.01 (1.33 to 3.05) | 2.29 (1.37 to 4.45) |
| Kernel-only, thickness-identical (all) | 84/54 | 45 | 0.82 (0.50 to 1.01) | 1.02 (0.57 to 1.38) | 4.15 (2.98 to 5.11) | 1.99 (1.32 to 3.39) | 2.29 (1.67 to 5.44) |
| thick_cat=<=1.25 | 6/5 | 1 | not fitted | not fitted | not fitted | not fitted | not fitted |
| thick_cat=1.25-3 | 34/21 | 11 | 0.67 (0.03 to 0.97) | 0.54 (0.22 to 3.37) | 3.40 (2.28 to 6.77) | 1.78 (0.59 to 2.82) | 3.79 (1.48 to 21.95) |
| thick_cat=>3 | 44/33 | 33 | 0.88 (0.58 to 1.39) | 1.16 (0.68 to 1.38) | 4.28 (2.46 to 5.63) | 2.03 (1.51 to 3.87) | 2.29 (1.10 to 4.52) |
| collection=LUAD | 20/16 | 15 | 0.82 (0.49 to 1.59) | 0.67 (0.38 to 1.38) | 4.28 (3.51 to 6.49) | 1.83 (1.08 to 4.33) | 1.67 (1.00 to 4.45) |
| collection=LUSC | 46/25 | 29 | 0.92 (0.47 to 1.01) | 0.84 (0.50 to 1.26) | 3.40 (2.31 to 5.03) | 2.01 (0.84 to 3.41) | 2.29 (1.47 to 6.55) |
| collection=ESCA | 18/13 | 1 | 0.60 (-38.41 to 28.13) | 4.48 (-35.71 to 73.37) | 6.50 (2.12 to 206.86) | 2.03 (1.13 to 4.39) | 24.66 (1.52 to 49.57) |

- ≤1.25 mm stratum: 6 pairs from 5 patients, below the pre-specified minimum, so it was not fitted.
- In the >3 mm, 1.25–3 mm, LUAD and LUSC strata, the intercepts lie in the same range as the primary result (heart 0.67–0.92%, aorta 0.54–1.16%, pulmonary vein 3.4–4.3%, SVC 1.8–2.0%, IVC 1.7–3.8%). Upper CI limits are below 10% except for the IVC in the 1.25–3 mm stratum (upper limit 21.9%).
- ESCA (18 pairs, 13 patients): only 1 pair per structure has |ZD| ≤ 5%, so the intercept is an extrapolation. The intercepts were heart 0.60, aorta 4.48, pulmonary vein 6.50, SVC 2.03 and IVC 24.7%, with very wide CIs (aorta −35.7 to 73.4; pulmonary vein 2.1 to 206.9; IVC 1.5 to 49.6). The ESCA data cannot confirm that the primary result holds in oesophageal-cancer scans. For the IVC, the ESCA point estimate exceeds 10%.
- Analysis B is descriptive. No criterion was pre-specified.

## 3. Deviations from the plan and implementation choices
1. **Thickness category for the target**: the plan does not say which member's thickness defines the category, and 102/186 pairs differ (exactly the thickness-arm pairs). Series A thickness was used, following the `slice_thickness_a_mm` convention of `retained_vs_excluded.csv`. *Exploratory* re-run with series B thickness (`representativeness_exploratory_seriesB_thickness.csv`): weighted medians 1.39/1.26/5.04/0.81/3.42% for heart/aorta/PV/SVC/IVC, all upper limits ≤ 8.64%, n_eff 9.4 pairs and 4.3 patients. The criterion conclusion is unchanged. This run is labelled exploratory and was not pre-specified.
2. **Contrast field**: the plan said "as recorded in headers_merged.csv". `ContrastBolusAgent` was used (§1). The "none recorded" level was empty in the target, so contrast had effectively two levels.
3. **Infeasible margins**: raking could not match the collection margin because the single ESCA pair is also a Siemens pair. The IPF fixed point was used unchanged; the plan did not specify a procedure for infeasibility.
4. **Trimming**: the plan says "trimmed at 5× the mean" without saying whether to renormalise or re-rake afterwards. A single cap at 5× the mean of the mean-1 raked weights was applied, with no re-raking.
5. **Uncovered levels within bootstrap resamples**: the plan defines uncovered levels for the observed sample only (there were none). In resamples lacking a level, that level's target share was dropped and the remaining shares of that variable renormalised. This was an implementation choice.
6. **Weighted-median convention**: lower weighted median with midpoint on an exact half; not specified in the plan.
7. **Table 2 reproduction**: only the manifest path was changed, as instructed.

No other subgroups, thresholds or exclusions were added.

## Files
- `representativeness_reweighted.csv`: Analysis A results per structure (weighted medians, CIs, ESS, trimming, bootstrap diagnostics).
- `representativeness_strata.csv`: Analysis B intercepts and slopes per stratum and structure, including the reproduction rows.
- `representativeness_margins.csv`: target vs retained vs weighted margins.
- `representativeness_pair_weights.csv`: raw and trimmed weight for each of the 45 pairs.
- `pair_covariates_186.csv`: derived covariates for all 186 pairs.
- `strata_qr.py`: Analysis B script.
- `representativeness_exploratory_seriesB_thickness.csv`: exploratory sensitivity analysis (Deviation 1).
