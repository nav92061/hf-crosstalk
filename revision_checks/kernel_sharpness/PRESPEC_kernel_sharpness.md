# Objective kernel sharpness: pre-specified plan

Written 2026-10-10 21:01 UTC, before any pixel data are downloaded for this analysis.

## Question
The signed-difference analyses order each pair from smoother to sharper kernel using vendor kernel families, which is approximate across vendors. Does the direction result hold when sharpness is measured from the images?

## Measurement
- Series: the 87 series of the 47-pair slice-position stratum (both members of every pair in the 45- and 47-pair strata).
- Slices: five slices at 40%, 45%, 50%, 55% and 60% of the stack (mid-thorax), identical relative positions in both members of a pair.
- Body mask: HU > -500, largest connected component, eroded by 10 pixels.
- Sharpness index (primary): fraction of 2D power-spectrum energy above 0.4 cycles/mm (frequencies expressed in cycles/mm so that pixel-size differences do not drive the index), computed within a 128 x 128 mm central body ROI, median over the five slices.
- Secondary: noise, measured as the SD of the residual after a 5 x 5 median filter within the body mask (HU -100 to 200), median over slices.
- Images are deleted after measurement.

## Analysis
- For each pair, the member with the higher sharpness index is defined as sharper. Report agreement with the vendor-family ordering (fraction of pairs; pairs with indices within 2% of each other reported as ties and excluded from the agreement fraction).
- Recompute, for all five structures in the 45-pair primary stratum, the patient-level count of sharper-larger and the median signed difference, using TotalSegmentator volumes, with a two-sided sign test.
- Criterion: the manuscript's direction claim is retained for a structure if, under the objective ordering, the sharper member is larger in at least two-thirds of patients with sign-test p < 0.05. Otherwise the claim is withdrawn or qualified for that structure.

## Reporting commitment
The result is reported in the manuscript or Supplement whatever it shows. No additional subgroups, thresholds or exclusions will be added after results are seen; any further analysis is labelled exploratory.
