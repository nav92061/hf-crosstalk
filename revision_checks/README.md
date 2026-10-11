# Revision checks added after the R1 response (each with a dated plan written before computation)

| Folder | Plan | Question | Result file |
|---|---|---|---|
| `representativeness/` | `PRESPEC_representativeness.md` | Does the 45-pair result depend on the retained mix? Raking to the 186-pair margins; no-cutoff model within strata | `representativeness_reweighted.csv`, `representativeness_strata.csv` |
| `kernel_sharpness/` | `PRESPEC_kernel_sharpness.md` | Does an image-measured sharpness index reproduce the vendor kernel ordering? | `kernel_sharpness_series.csv`, `kernel_sharpness_direction.csv` |
| `external_accuracy/` | `PRESPEC_external_accuracy.md` | TotalSegmentator 2.18.0 whole heart against expert RTOG contours on TCIA LCTSC (60 CT) | `external_accuracy_cases.csv`, `external_accuracy_summary.csv` |

Each folder's results memo (`REPRESENTATIVENESS.md`, `KERNEL_SHARPNESS.md`, `EXTERNAL_ACCURACY.md`) lists deviations from its plan.
The sharpness and accuracy scripts download images from TCIA; the accuracy job needs an NVIDIA GPU.
