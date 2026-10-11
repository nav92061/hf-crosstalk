# Representativeness of the retained stratum: pre-specified plan

Written 2026-10-10 21:01 UTC, before any computation in this analysis.

## Question
Does the primary result (45 retained kernel-only pairs) depend on the retained stratum's mix of slice thickness, cancer collection, vendor and contrast status, which differ from the 186 segmented pairs (retained_vs_excluded.csv)?

## Analysis A (reweighting)
- Population: the 45 primary-stratum pairs (endpoint1_clean_stratum definition, image count x thickness extent).
- Target: the distribution among all 186 segmented pairs of: slice-thickness category (<=1.25, 1.25-3, >3 mm), collection (LUAD, LUSC, ESCA, other), vendor (GE, Siemens, Philips, other) and contrast status (as recorded in headers_merged.csv; unknown as its own level).
- Method: raking (iterative proportional fitting) of pair weights to the target margins, weights trimmed at 5x the mean. Levels absent from the retained stratum cannot be reweighted and are reported as uncovered with their share of the target.
- Outcome: weighted patient-level median PD per structure (patient weight = sum of its pair weights), 95% CI from 10,000 patient-level bootstrap resamples with weights recomputed in each resample, seed 0; effective sample size reported.
- Criterion: the primary result is considered robust to the retained mix if every structure's weighted upper limit is below 10%.

## Analysis B (stratification without a coverage cutoff)
- Population: all segmented kernel-only, thickness-identical pairs.
- Method: the principal no-cutoff model (median regression of PD on |ZD| with patient-cluster bootstrap, 2,000 resamples, seed 2026), fitted separately within slice-thickness category and within collection, for strata with >=8 pairs from >=5 patients.
- Outcome: intercept (PD at ZD = 0) and 95% CI per stratum and structure. Descriptive; no criterion.

## Reporting commitment
The result is reported in the manuscript or Supplement whatever it shows. No additional subgroups, thresholds or exclusions will be added after results are seen; any further analysis is labelled exploratory.
