# External accuracy check of TotalSegmentator whole-heart segmentation: pre-specified plan

Written 2026-10-10 21:01 UTC, before any segmentation in this analysis.

## Question
The study measures reproducibility, not accuracy. Does TotalSegmentator 2.18.0, run as in the main study, segment the whole heart accurately against independent manual contours on comparable CT (lung-cancer radiotherapy planning CT)?

## Data
- TCIA LCTSC (Lung CT Segmentation Challenge 2017; 60 thoracic CT scans with expert heart contours following RTOG 1106 guidelines). Eligibility of the dataset is confirmed before use: public download, licence allowing this use, and not listed among TotalSegmentator's training sources. If any of these fails, the analysis is not run and the reason is reported.

## Method
- TotalSegmentator 2.18.0, task total, default (not fast) mode, same version and settings as the main study, on GPU. One-case pilot first.
- Reference: RTSTRUCT contour named Heart (case-insensitive match), rasterised to the CT grid.
- Metrics per case: Dice, 95th-percentile Hausdorff distance (mm), volume difference (TotalSegmentator minus reference, mL and %).
- Summary: median and IQR of Dice and HD95; Bland-Altman mean bias and 95% limits of agreement for volume.
- Cases where TotalSegmentator fails, the heart contour is absent, or the CT and RTSTRUCT do not align are listed with the reason, not dropped silently.

## Interpretation (fixed in advance)
- Radiotherapy heart contours (RTOG) start at the pulmonary trunk bifurcation and include pericardium; TotalSegmentator's heart class has a different superior boundary. A volume bias is therefore expected and is reported, not corrected.
- Median Dice >= 0.90 will be described as consistent with accurate whole-heart segmentation on this CT type; below 0.90 will be reported as a limitation of accuracy. The check does not establish accuracy on the study's own scans, for great vessels, or for chambers.

## Reporting commitment
The result is reported in the manuscript or Supplement whatever it shows. No additional subgroups, thresholds or exclusions will be added after results are seen; any further analysis is labelled exploratory.
