# Independent-model whole-heart check: pre-specified plan

Written 2026-10-10 12:00 EDT before any PlatiPy output exists. Saved as an artifact; its timestamp precedes the first job submission.

## Why this model
The PlatiPy open whole-heart model (Finnegan et al., Phys Eng Sci Med 2023; Zenodo record 6585664) is an nnU-Net 3d_lowres network trained on the TCIA LCTSC and NSCLC-Radiomics collections, whose heart labels are clinical radiotherapy contours, not TotalSegmentator output. None of the study's TCGA collections are in its training data.

## Model and inputs
- platipy 0.7.2, nnunet 1.7.1; model Task400_OPEN_HEART_1FOLD (3d_lowres), weights zip MD5 feb63660bce90c2b6ecf561c8463b6a9 (verified at build; job halts otherwise).
- Inference with PlatiPy's own whole-heart settings (HYBRID_SETTINGS_DEFAULTS["nnunet_settings"]), unmodified. Only the whole heart is segmented; the atlas substructure step is not run.
- Same 87 series as the primary stratum under either z-extent definition (vista_job/series.csv), downloaded from TCIA one at a time.
- Volume = voxel count x voxel volume on the native image grid. No post-processing beyond the model's own.

## Acceptance checks (applied in code before any volume is used)
1. Weights MD5 matches; otherwise the run halts.
2. Every series yields a row: a volume or a recorded error. Nothing is dropped silently.
3. Converted-image z-extent agrees with the slice-position extent within one slice spacing; failing series are excluded and listed.
4. The output contains exactly one label (whole heart) and is non-empty.
5. Whether the heart touches the first or last slice is recorded.

## Analysis (fixed now)
- The released primary-endpoint script (imaging/endpoint1_tables.py) is applied unchanged, with PlatiPy heart volumes substituted for TotalSegmentator heart volumes. Only the heart rows are reported; other structures are not evaluated.
- Primary: 45-pair stratum (original z-extent definition). Also reported: 47-pair stratum (slice-position extent) and the 41-pair field-of-view-matched stratum already in the manuscript.
- Reproduction criterion: the whole-heart result is reproduced if both the patient-level median bootstrap interval and the mixed-model interval lie below 10%, as in the primary analysis.
- Secondary: number of patients in whom the sharper kernel gave the larger heart volume; between-model agreement per series (median signed percentage difference PlatiPy minus TotalSegmentator, Spearman rho).
- No other subgroups, thresholds or outlier rules will be added after results are seen.

## Reporting commitment
The result is reported in the manuscript or Supplement whatever it shows, including failure of the reproduction criterion. Absolute heart volumes are expected to differ between models because radiotherapy contouring guidelines define the heart differently; this does not affect within-pair differences, which use the same model on both members of a pair.
