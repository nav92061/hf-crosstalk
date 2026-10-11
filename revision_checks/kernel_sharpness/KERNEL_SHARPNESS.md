# Objective kernel sharpness and direction of kernel-related volume differences

Executed 2026-10-10 22:48 UTC following `PRESPEC_kernel_sharpness.md` (plan written 2026-10-10 21:01 UTC, before any pixel data were downloaded).

## Sanity check
`endpoint1_tables.py` was run unchanged, except that the absolute manifest path was replaced with the local copy. Its output matched the published Table 2 (`endpoint1_clean_stratum.csv`) cell for cell: 45 pairs, 31 patients. The regenerated `endpoint1_pairs.csv` was numerically identical to the supplied 186-pair input.

## Measurement (step 1)
- Series: 87, the union of both members of every pair in the 45-pair primary stratum and the 47-pair slice-position stratum. The 45-pair set is a subset of the 47-pair set, which adds 2 pairs. All 87 were measured, 5 slices each (435 slices). 60 series were GE, 22 Siemens and 5 Philips.
- Source: public TCIA NBIA API (getSOPInstanceUIDs, getSingleImage). Instances were sorted by their projection of ImagePositionPatient onto the slice normal (cross product of the ImageOrientationPatient row and column vectors). The slices at rank round(p*(N-1)) for p = 0.40, 0.45, 0.50, 0.55 and 0.60 were taken, and RescaleSlope/Intercept were applied.
- Body mask: HU > -500, largest connected component, eroded by 10 pixels (disk).
- Primary index: fraction of 2D power-spectrum energy at radial spatial frequency > 0.4 cycles/mm. Frequencies were computed from PixelSpacing. The DC term was excluded, and a 128 x 128 mm ROI was centred on the eroded-mask centroid. Median over the 5 slices.
- Secondary index: SD of (image minus 5 x 5 median-filtered image) within the eroded body mask restricted to HU -100 to 200. Median over the 5 slices.
- Images were held in memory only, never written to disk, and discarded after each slice was measured.

Per-kernel summary (all 87 series):

| Vendor | Kernel | n series | Median HF index | Median noise SD (HU) |
|---|---|---|---|---|
| GE | SOFT | 12 | 0.00132 | 15.6 |
| GE | STANDARD | 17 | 0.00155 | 16.6 |
| GE | DETAIL | 1 | 0.00141 | 14.3 |
| GE | LUNG | 29 | 0.03873 | 46.5 |
| GE | PET AC | 1 | 0.00001 | 3.2 |
| PHILIPS | B | 1 | 0.00035 | 11.1 |
| PHILIPS | C | 1 | 0.00077 | 12.8 |
| PHILIPS | D | 1 | 0.00173 | 18.1 |
| PHILIPS | YA | 1 | 0.03526 | 76.1 |
| PHILIPS | YB | 1 | 0.01543 | 48.4 |
| SIEMENS | B20f | 4 | 0.00024 | 12.1 |
| SIEMENS | B30s | 1 | 0.00590 | 16.9 |
| SIEMENS | B31f | 2 | 0.00099 | 13.4 |
| SIEMENS | B40s | 3 | 0.00091 | 18.6 |
| SIEMENS | B41s | 1 | 0.00009 | 6.1 |
| SIEMENS | B60f | 4 | 0.02585 | 63.1 |
| SIEMENS | B70s | 5 | 0.00904 | 32.2 |
| SIEMENS | B80f | 2 | 0.05416 | 70.6 |

## Agreement with vendor-family ordering
Vendor ordering: GE SOFT<STANDARD<DETAIL<LUNG<BONE/BONEPLUS; Philips A<B<C<D<YA<YB; Siemens by the number in the B/I kernel name.

| Set | Pairs | Vendor-defined | Ties (<=2%) | Agree / non-tied | Fraction |
|---|---|---|---|---|---|
| 45-pair primary | 45 | 45 | 0 | 45/45 | 1.000 |
| 47-pair slice-position | 47 | 46 | 0 | 46/46 | 1.000 |
| union 47 pairs/87 series | 47 | 46 | 0 | 46/46 | 1.000 |

The measured index ordered every vendor-defined pair the same way as the vendor family. There were no ties: within-pair relative differences in the index ranged from 77% to 197% (median 185%) in the primary stratum. The secondary noise index also agreed in all 46 vendor-defined pairs. The one pair the vendor scheme cannot order, GE STANDARD vs PET AC (TCGA-50-8460), is outside the 45-pair stratum. The index ranks PET AC as the smoother member.

## Direction analysis (step 2): 45-pair primary stratum, TotalSegmentator volumes, objective ordering
Signed difference = 100*(V_sharper - V_smoother)/mean. Patient-level median over the patient's pairs; patient bootstrap, 10,000 resamples, seed 0; two-sided exact sign test on patient medians.
Criterion (pre-specified): the claim is retained if the sharper member is larger in >= 2/3 of patients and sign-test p < 0.05.

| Structure | Patients | Sharper larger | Sharper smaller | Fraction | Median signed PD % (95% CI) | Sign-test p | Criterion |
|---|---|---|---|---|---|---|---|
| heart | 31 | 27 | 4 | 0.871 | 0.63 (0.49 to 0.94) | 3.4e-05 | met: retained |
| aorta | 31 | 26 | 5 | 0.839 | 0.97 (0.67 to 1.16) | 0.00019 | met: retained |
| pulmonary_vein | 31 | 28 | 3 | 0.903 | 4.15 (2.32 to 5.31) | 4.6e-06 | met: retained |
| superior_vena_cava | 31 | 27 | 4 | 0.871 | 2.01 (1.78 to 3.39) | 3.4e-05 | met: retained |
| inferior_vena_cava | 31 | 15 | 16 | 0.484 | -0.06 (-1.69 to 0.70) | 1 | not met: withdraw/qualify |

Because objective and vendor orderings agree in all 45 primary pairs, the vendor-family results are numerically identical. They also match `agreement_primary_stratum.csv` (frac_patients_sharper_larger and patient_median_sharp_minus_smooth). Full rows for both orderings are in `kernel_sharpness_direction.csv`.

## Conclusion against the pre-specified criterion
- Retained under objective ordering: heart, aorta, pulmonary vein, superior vena cava.
- Not met: inferior vena cava, with 15/31 patients sharper-larger, median signed difference -0.06% and p = 1.0. Any direction statement for the IVC must be withdrawn or qualified as showing no consistent direction. This matches the existing vendor-ordering result; the objective ordering does not change it.

## Deviations and implementation details (not departures from the plan's decision rules)
1. The plan does not specify these details, so they were chosen before results were seen:
   - a Hann window was applied to the ROI before the FFT to limit edge leakage;
   - the ROI was centred on the eroded body-mask centroid and clipped to the image;
   - the ROI is a square in the image and is not restricted to mask pixels, so it can include lung (minimum in-body fraction 0.17);
   - the 2% tie rule was applied as |a-b|/mean(a,b).
2. The task text asks for the "equivalent axis" for the 2 series labelled coronal in `series_geometry_ipp.csv` (TCGA pair B41s/B70s). Sorting used the IOP-derived slice normal, which is (0,0,1) for both series. The measured slices are therefore axial images ordered along z, and the method is identical for all series.
3. Three series (one pair in the primary stratum, TCGA-50-5066 STANDARD/LUNG, plus the PET AC series) have 1.37 mm pixels, so Nyquist is 0.366 cycles/mm, below the 0.4 cutoff. For these the index reflects only the diagonal spectral corners. It is small in absolute terms but still separated the pair members about 10-fold, in the vendor-predicted direction. Four primary pairs differ slightly in pixel spacing; the cycles/mm formulation handles this as planned.
4. Operational: the first two background download runs hung after 49 series (the process stalled with no API error). They were stopped and the resumable script was rerun from series 50 in the analysis kernel with the same code; measurement method and outputs were unaffected.
5. No subgroups, thresholds or exclusions were added. No exploratory analyses were run.
