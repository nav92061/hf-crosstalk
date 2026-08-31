# Draft manuscript sections, imaging arm

Insert as new subsections at the end of Methods, Results and Discussion. Numbering assumes
the cuts in `MANUSCRIPT_RESTRUCTURE.md` have been made (3.6, 3.8, 3.10 removed), so the new
Results subsections are 3.9–3.11 in the revised paper; they are labelled 2.11 / 3.11–3.13 /
4.7 here to match the planning documents.

Every number below is in a saved artifact. Nothing is estimated.

---

## 2.11 Imaging cohort, reconstruction pairing, and segmentation

Computed tomography series were drawn from four TCGA collections in The Cancer Imaging
Archive (TCGA-LUSC, TCGA-LUAD, TCGA-ESCA, TCGA-SARC). DICOM headers for one representative
slice per series were retrieved through the archive's `getSOPInstanceUIDs` and
`getSingleImage` endpoints, which return a complete header for approximately 0.5 MB per
series and require no volume transfer; 573 series were probed with no fetch errors.

The analysis unit is a *reconstruction pair*: two series belonging to one study that share
`FrameOfReferenceUID` and `AcquisitionTime`, and therefore derive from a single irradiation
event, but differ in `ConvolutionKernel` or `SliceThickness`. Shared frame of reference and
acquisition time were verified by a second header pull for the 234 series involved rather
than inferred from study membership. This yielded **193 confirmed pairs across 217 distinct
series in 59 patients**. Twenty-one distinct reconstruction kernels appear in the cohort
(GE `STANDARD`, `LUNG`, `SOFT`, `BONEPLUS`, `DETAIL`; Siemens `B10f`–`B80f`; Philips `YA`,
`YB`), and slice thickness spans 0.625–5.0 mm.

Segmentation used TotalSegmentator 2.18.0 (nnU-Net 2.8.1, PyTorch 2.10.0+cu128, Tesla T4),
open `total` task, with the six model checkpoints hashed and recorded
(`model_provenance.json`). The five pre-specified structures are **heart, aorta, pulmonary
vein, superior vena cava and inferior vena cava**. The pre-specification originally named
the four cardiac chambers and myocardium; inspection of the distributed package showed
those come from task `heartchambers_highres`, which the package marks `commercial: True`
and gates behind a licence. The structure list was re-specified on the open model before
any segmentation was run, and the amendment is dated in `PRESPECIFICATION_IMAGING.md`.

Volumes were taken from the tool's own `statistics.json` (`-s -sx`) rather than recomputed
downstream. `--stats_include_incomplete` was set: by default the tool reports volume 0.0
for any structure touching the image boundary, and thoracic staging CT routinely clips the
great vessels. Volumes for truncated structures are therefore **partial-field measurements,
not anatomical totals**, and are not reported as absolute anatomy. Fast mode was not used;
it resamples to 3.0 mm, which would attenuate the thickness differences under study.

---

## 3.11 A pre-specified eligibility check halts the cross-study endpoint

The pre-specified primary endpoint compared chamber volumes between two *separate* studies
of the same patient, requiring acquisition discordance, an interval of 90 days or less, and
concordant contrast status. Applying those criteria to the header table before any
segmentation ran reduced 35 patients with repeat imaging to 31 with acquisition
discordance, 17 within 90 days, and **5 with verifiable contrast concordance**.

The binding constraint is the contrast field. `ContrastBolusAgent` is populated in 21 of 81
primary series; recovering contrast status from every series in each study, not only the
segmentable one, gained a single study. Of 74 within-patient study pairs, only 20 have
contrast status known on both sides, and 36% of those are contrast-discordant, so the
pairs with unknown status cannot be assumed concordant.

At n = 5 the endpoint cannot reach its own 10% kill threshold under any plausible
within-patient variance. The endpoint was withdrawn before any segmentation was run.
This is reported rather than replaced silently: the eligibility criteria were fixed in
advance, they were evaluated on metadata that cost no compute, and they failed.

---

## 3.12 Reconstruction kernel changes AI-derived volumes by under 5%

A design available in the same data replaces the withdrawn endpoint. Within a single study,
the same acquisition is frequently reconstructed more than once. Comparing those
reconstructions removes contrast, interval and anatomical change by construction: the
patient, the irradiation event and the moment are identical, and only the reconstruction
differs.

One pre-specification assumption proved false and is the reason this section reports a
restricted stratum. We had asserted that pairs sharing an acquisition therefore share
z-coverage. Measured, within-pair z-coverage differs by a median of 15% and by up to 158%.
Reconstructions of one event are routinely exported over different ranges: a chest–abdomen
series and a heart-focused reformat share raw data but not extent. Volume difference correlates with
coverage difference (Spearman ρ = 0.70 for aorta, 0.68 for inferior vena cava, 0.50 for
heart), so a substantial part of the apparent reconstruction effect is field-of-view
truncation. Analysis is therefore restricted to pairs matched on z-coverage within 5% with
slice thickness held identical, leaving **45 pairs in 31 patients**. This filter is
post-hoc and declared as such; it was not pre-specified, it was forced by a pre-specified
assumption failing.

With the patient as the unit of inference and the pair as the unit of measurement, both
pre-specified aggregations agree (Table X; Fig. Y):

| Structure | Median per patient | 95% CI | Mixed-model mean | 95% CI |
|---|---|---|---|---|
| Heart | 0.63% | 0.50–0.97 | 0.85% | 0.64–1.07 |
| Aorta | 0.98% | 0.69–1.17 | 1.04% | 0.80–1.29 |
| Pulmonary vein | 4.96% | 3.20–5.63 | 5.41% | 3.93–6.89 |
| Superior vena cava | 2.12% | 1.90–3.46 | 2.53% | 1.86–3.20 |
| Inferior vena cava | 1.52% | 1.15–2.22 | 2.73% | 1.54–3.91 |

Every interval excludes the pre-specified 10% threshold on both aggregations. Whole-heart
volume is robust to the reconstruction kernel; the thin-walled pulmonary vein is roughly
eight times more sensitive than the heart, the expected direction for a structure whose
apparent volume depends on partial-volume effects at its boundary.

Segmentation failed on 4 of 217 series (1.8%), all with error `MISSING_DICOM_FILES`. These
are incomplete instance sequences in the archive, and the failure arises during format
conversion, before the model is invoked. All four are GE
(4/150) against 0/58 Siemens and 0/9 Philips. Six further series required relaxed DICOM
validation and then succeeded. This is an archive-integrity rate, not a model-reliability
rate, and is reported as such.

---

## 3.13 Association stability under reconstruction change: a null result

The remaining pre-specified endpoint asked whether the screen's own receptor panel produces
stable imaging associations. All 100 floor-passing receptors were scored per patient from
TCGA primary-tumour expression (obligate complexes scored as the minimum of their
subunits), and each receptor–structure Spearman correlation was computed twice from the
same patients and the same acquisitions, differing only in reconstruction.

Correlations are highly reproducible: median |Δρ| = 0.014 across 500 associations.

The endpoint was framed in terms of significance verdicts, and that framing does not
survive scrutiny. Of associations reaching nominal p < 0.05 under at least one
reconstruction, a median of **28% (95% range 5–47%)** are significant under one
reconstruction and not the other, integrating over the two arbitrary analyst choices
available, which pair is used for patients contributing several, and which member of a
pair is labelled first. A directional statement is not available: the pair labelling
carries no kernel ordering, and swapping the labels inverts any "proportion lost" figure.

Calibration against a permutation null settles the interpretation. Permuting patient labels
on the expression matrix destroys any real receptor–volume association while preserving the
pair structure, the reconstruction difference, n, and the p-value density near the
threshold. The null discordance is **22% (95% range 8–41%)**, and the observed value does
not exceed it (p = 0.33). Discordance of this magnitude is a property of dichotomising at
p < 0.05 with n = 31, arising equally for real and spurious associations, and is not a
specific instability of the segmentation output. Separately, **no association survives
Benjamini–Hochberg correction across the 500 tests** (minimum q = 0.98), so there are no
findings here whose stability could be at issue.

Endpoint 2 therefore reports a null. NPR1, the receptor carried forward from the
transcriptomic screen, shows no sign change and no threshold crossing on any structure
(maximum |Δρ| = 0.046); its heart association (ρ ≈ 0.32, p ≈ 0.08) is non-significant under
both reconstructions.

---

## 4.7 What the imaging arm establishes

The arm supports one measurement claim and two methodological ones.

The measurement claim is narrow and quantified: for a widely used open segmentation model,
whole-heart volume derived from thoracic CT is reproducible to about 1% across
reconstruction kernels when field of view is matched, while thin vascular structures vary
several times more. Anyone using automated chamber or vessel volumetry across
heterogeneously reconstructed archives now has an effect size for the reconstruction
component of their measurement error, on this model and this population.

The first methodological claim concerns what prose cannot enforce. Our own pre-specification
asserted that same-acquisition pairs share z-coverage. It was written in advance, in good
faith, and it was false, with ρ up to 0.70 between coverage difference and apparent volume
difference, an unguarded analysis would have attributed field-of-view truncation to the
reconstruction kernel. Measurement caught it; the pre-specification did not, because an
assumption stated as a sentence cannot fail. Had the same requirement been written as an
assertion over the header table, it would have halted at the audit stage rather than after
segmentation. This is the paper's thesis applied to the paper's own authors, and it is the
strongest evidence in it that the distinction between executable and asserted requirements
is not stylistic.

The second concerns thresholds. Endpoint 2 was designed to ask whether reconstruction
choice destabilises radiogenomic associations, and the honest answer is that at n = 31 the
question cannot be separated from the fragility of the p < 0.05 cut itself: the observed
discordance matches a permutation null in which no association is real. A study reporting
only dichotomised outcomes on a single reconstruction would have reported associations that
another equally valid reconstruction of the same scan does not support, but so would a
study perturbed in any other way by a comparable amount.

Scope is bounded and does not improve with more analysis. Fifty-nine patients, one
pre-trained model, one archive, and cancer patients scanned for cancer staging rather than
patients with heart failure. The imaging population and the transcriptomic ligand source
(explanted failing hearts, GSE116250 and GSE141910) have no biological overlap, so the arm
tests the durability of the screen's output under measurement variation and does not bear
on the heart-failure hypothesis. Any sentence implying otherwise is a defect.

The arm also produced the project's third pre-specified withdrawal, after the seven
pleiotrophin axes and the cross-study imaging endpoint. Two of the three were withdrawn
before compute was spent, on metadata alone.
