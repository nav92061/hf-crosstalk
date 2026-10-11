# External accuracy check of TotalSegmentator whole-heart segmentation (TCIA LCTSC)

Plan: PRESPEC_external_accuracy.md (written 2026-10-10 21:01 UTC, before any segmentation). Analysis run 2026-10-10 22:42 to 2026-10-11 01:27 UTC. Report written 2026-10-11 01:28 UTC.

## Result
- TotalSegmentator ran on 60 cases. 59 have metrics; 1 failed (listed below).
- **Dice: median 0.756 (IQR 0.734–0.791).** The range was 0.683–0.837; no case reached 0.90.
- **HD95: median 23.7 mm (IQR 22.3–25.9).** By direction:
  - TotalSegmentator surface to reference: median 29.6 mm.
  - Reference surface to TotalSegmentator: median 12.4 mm.
- **Volume, TotalSegmentator − reference (Bland–Altman):** mean bias -174.5 mL, 95% limits of agreement -308.6 to -40.4 mL.
  - As a percentage of the reference volume: bias -23.0%, limits -36.7% to -9.2%.
  - Median volumes: reference 732 mL, TotalSegmentator 573 mL.
- **Interpretation under the pre-specified rule:** the median Dice is below 0.90, so this is a **limitation of accuracy**. Against RTOG-1106 radiotherapy heart contours on lung-cancer planning CT, TotalSegmentator's whole-heart mask overlaps only moderately and is systematically about 23% smaller. As the plan expected, part of this reflects different definitions of the heart, and the bias is reported, not corrected. The check does not establish accuracy on the study's own scans, for great vessels or for chambers. It does not bear on reproducibility across kernels, which the main study measures as within-scanner differences that cancel any constant definitional bias. That last point is a reasoning remark, not a result of this analysis.

## Dataset eligibility
- **Availability:** the TCIA (NBIA API v1) query for Collection=LCTSC returned 60 CT and 60 RTSTRUCT series from 60 patients, one pair per study. Download needed no credentials.
- **Licence:** CC BY 3.0 Unported on all series. Collection DOI 10.7937/K9/TCIA.2017.3r3fvz08; cite it as such.
- **Not in TotalSegmentator's training data:** LCTSC is not listed in Wasserthal et al. 2023 (University Hospital Basel PACS sample), in the Zenodo datasets v1 (6802614) or v2.0.1 (10047292), or in the GitHub README or improvements_in_v2.md.
- **Residual uncertainty:** about 420 of the v2 `total` training images are not public and are described only as including "other institutions". Nothing indicates that LCTSC is among them.
- Evidence is in lctsc_eligibility_evidence.json.

## Methods as run
- **Software:** TotalSegmentator 2.18.0, task `total`, default mode (not `--fast`), multilabel output with heart = label 51. Command: `TotalSegmentator -i <CT DICOM dir> -o seg.nii.gz --ml --task total`. Environment: torch 2.14.1+cu130, nnunetv2 2.8.1, NVIDIA A10 (Modal).
- **Image and weights:** Modal image im-yv1flLIoZ7jwNBBrTZTMDi; weights were downloaded when the image was built. Package versions are in pip_freeze.txt.
- **Jobs:** the pilot (job a9e087f5, LCTSC-Test-S1-101) was approved and ran in 151 s. The other 59 cases ran as job 491d8fd2 (9366 s) and resumed from the pilot's row.
- **Reference:** the RTSTRUCT ROI named "Heart" (case-insensitive exact match) was rasterised onto the CT grid with the even-odd polygon rule.
- **Alignment checks:** FrameOfReferenceUID matched in 60/60 cases; the referenced CT series matched in 60/60; the largest contour off-plane distance was 0.0 mm.
- **Grid:** TotalSegmentator labels were resampled to the native CT grid with nearest-neighbour interpolation. CT slice spacings were 3.0 mm (n=27), 2.5 mm (24), 2.0 mm (7) and 1.25 mm (1).
- **Metrics:**
  - Dice.
  - HD95: the 95th percentile of pooled bidirectional surface distances, using 6-connected surface voxels, in mm.
  - Volume difference, TotalSegmentator − reference, in mL and as a % of the reference.
  - Bland–Altman: mean ± 1.96 SD (sample SD).

## Cases without metrics (not dropped silently)
| PatientID | status | reason |
|---|---|---|
| LCTSC-Test-S2-202 | failed | TotalSegmentator's internal DICOM-to-NIfTI conversion (dicom2nifti) raised ConversionError MISSING_DICOM_FILES, meaning it detected a gap in the slice sequence. The reference Heart ROI was present (37 contour slices, aligned). |

## Exploratory (not pre-specified)
These are descriptive diagnostics of where the disagreement lies, computed from the superior and inferior z-extents already recorded per case (external_accuracy_exploratory_extent.csv). They are not part of the pre-specified analysis.
- **Superior boundary:** the TotalSegmentator heart mask extends a median of 35.0 mm (IQR 32.2–39.5, range 24–57) further superiorly than the RTOG contour.
- **Inferior boundary:** the two agree, median difference 0.0 mm (IQR 0.0–3.0).
- **Asymmetric HD95:** TotalSegmentator surface points lie far from the reference, mostly at the extra superior extent, while reference points lie close to TotalSegmentator.
- Taken together with the 23% smaller volume, the superior extent beyond the RTOG boundary and the asymmetric HD95 point to systematic differences in how the boundary is defined, not misplacement. The excess volume of the RTOG contour presumably lies outside the cranio-caudal extent difference, for example pericardium or pericardial fat included in the RTOG heart. This was not tested.
- The inferior-boundary agreement near 0 mm also argues against a registration or orientation error.

## Deviations from the plan
- **Pilot approval:** the brief allowed up to 15 minutes for approval of the pilot submission. It waited about 31 minutes, after which I stopped the local submission call and wrote an interim "blocked" report. The pilot had in fact been approved late and ran to completion (job a9e087f5). The interim report's statements that no job was created and nothing was billed were wrong, and this report supersedes them. The full run was then submitted and approved normally.
- **Analysis:** no deviations. No subgroups, thresholds or exclusions were added. Implementation choices fixed before any results: pooled bidirectional HD95; labels resampled to the CT grid; the reference excluded if FrameOfReferenceUID mismatched or contours lay more than half a slice off-plane (neither occurred).
