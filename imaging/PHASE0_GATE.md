# Phase 0 — Header audit and decision gate

**Verdict: PROCEED to Phase 1.** Both gate conditions are met. But the audit found three
things that change the Phase 2 design, and one of them is decisive on power grounds.

Everything below is computed from DICOM headers actually retrieved, not from series-level
metadata. Route: TCIA `getSOPInstanceUIDs` → `getSingleImage` for the **middle slice** of
each series, header read with `stop_before_pixels=True`. 573 series probed, **0 fetch
errors**, ~0.5 MB per series. No image volume was downloaded — the plan's instruction to
check this before pulling volumes was correct and the answer is that headers are
effectively free.

## Cohort reconstruction

The 83-patient cohort reproduces exactly from `cardiac_ct_feasibility.json`:
LUSC 35, LUAD 29, ESCA 16, SARC 3, spanning 573 expression-matched CT series with
`BodyPartExamined` in LUNG / CHEST / ESOPHAGUS / CHESTABDPELVIS.

**81 of 83 patients have at least one segmentable series.** Series-level exclusions
(non-exclusive categories): 145 localizer/scout, 61 secondary capture, 212 under 40
slices, 6 with no `SliceThickness` → **304 seg-usable series in 81 patients**.

The two dropped patients (TCGA-38-4625, TCGA-38-4628) have only 8 mm axial series with
<40 slices plus localizers. 8 mm through-plane resolution cannot support cardiac chamber
segmentation; they are excluded before Phase 1 rather than allowed to fail silently in it.

## Header coverage

Two populations matter and they differ, so both are given: the 81 primary series that
enter the analysis, and all 573 probed series.

| Field | 81 primary series | All 573 probed |
|---|---|---|
| Manufacturer, ManufacturerModelName | 100% | 100% |
| StudyDate | 100% | 100% |
| KVP | 100% | 99.0% |
| SliceThickness, PixelSpacing | 100% | 99.0% |
| ConvolutionKernel | 100% | 95.8% |
| XRayTubeCurrent | 100% | 92.0% |
| SoftwareVersions | 100% | 94.6% |
| **ContrastBolusAgent** | **25.9% (21/81)** | **16.8% (96/573)** |

Every acquisition parameter the design needs is present in **all 81** patients selected
for analysis; the sub-100% figures apply to the full probed set, which includes localizers
and secondary captures that are excluded anyway. `ContrastBolusAgent` is the sole field
that is sparse in the analysis population.

## Gate 1 — vendor strata: PASS

Patient-level, using one primary series per patient (thinnest slice, then most images):

- **GE 55**, **Siemens 23**, Philips 3

Two strata clear n ≥ 15. The kill rule ("if one vendor holds 70+ of 83, stop") is **not
tripped**: the largest vendor holds 55/81 = 68%.

## Gate 2 — continuous-parameter spread: PASS

| Parameter | Levels | Range |
|---|---|---|
| Slice thickness | 8 | 0.625 – 5.0 mm (8×) |
| Tube current | 64 | 40 – 713 mA |
| In-plane spacing | 55 | 0.43 – 1.37 mm |
| kVp | 3 | 120 / 130 / 140 |

Slice thickness is the strongest axis and, importantly, **varies within both vendors**
(GE 0.625–5.0, Siemens 1.0–5.0), so it is not merely a vendor proxy.

## What the audit found that the plan did not anticipate

### 1. kVp is partially confounded with vendor — estimable, but on thin strata

| | 120 kVp | 130 kVp | 140 kVp |
|---|---|---|---|
| GE | 38 | 0 | 17 |
| Siemens | 9 | 14 | 0 |
| Philips | 2 | 0 | 1 |

**Corrected from an earlier draft of this document.** I first wrote that the planned model
`~ vendor + slice thickness + contrast + kVp` was "not identifiable." That was wrong, and
checking it rather than asserting it is what caught it. Because GE spans 120/140 *and*
Siemens spans 120/130, within-vendor kVp contrasts exist, the design matrix is **full rank
(5 of 5)**, and all VIFs are **below 2.3** — no serious inflation. kVp is estimable.

What is true is narrower: **130 kVp occurs only in Siemens and 140 kVp never in Siemens**,
so those two cells are vendor-confounded, and the kVp effect rests on small within-vendor
strata (Siemens 9 vs 14; GE 38 vs 17). An unqualified "140 kVp effect" would be carried
almost entirely by GE patients.

**Revised recommendation:** kVp may stay in the model, but report it as **within-vendor
contrasts** (GE 120 vs 140; Siemens 120 vs 130) rather than a single pooled coefficient,
and state the per-cell n alongside. This is a weaker claim than my first draft made and a
more defensible one.

Worth noting for the paper: `np.linalg.lstsq` returns finite coefficients for this design
with no warning, and would do so even in the genuinely rank-deficient case. Standard
tooling does not halt on structural confounding; a guardrail has to be written for it.

Vendor × slice-thickness is *associated* but not collinear (χ² p = 0.0010, Cramér's V =
0.412), and thickness retains within-vendor spread, so both can stay in the model.

### 2. Contrast cannot be modelled as a binary

`ContrastBolusAgent` is populated in only **21 of 81** primary series. Worse, in 17
further series the tag is absent while the series description explicitly implies contrast
("with contrast", "arterial", "portal"). Treating absence as "no contrast" would
mislabel at least 17 patients.

After reconciling the tag against series/protocol/study descriptions:
**contrast 38, non-contrast 4, unknown 39.**

**Correction:** contrast enters as a **three-level factor including `unknown`**, never a
binary. Only 4 patients are confirmed non-contrast, so the contrast-vs-non contrast is
uninformative at this n and should not be reported as an estimated effect. This is worth
stating in the paper: a widely used DICOM field is unusable for retrospective
stratification in 48% of a real archive cohort.

### 3. Test-retest is not a conditional side-analysis — it is the stronger design

The plan treated Phase 3 as optional ("if any patient has repeat scans... even 5–10
pairs"). The header audit finds far more:

- **35 of 81 patients have ≥2 seg-usable studies** (52 within-patient study pairs)
- **31** were scanned on a **different scanner model** between studies
- **27** at a **different slice thickness** (median 2.0× change, up to 5.0×)
- 7 on a **different vendor**
- median interval 46 days (IQR 22–130)

**This inverts the design.** Minimum detectable effect at 80% power, Bonferroni-corrected
across the 5 pre-specified structures:

| Design | n | Min. detectable effect |
|---|---|---|
| Between-patient, GE vs Siemens | 55 vs 23 | Cohen d = 0.849 |
| Within-patient, thickness-discordant | 27 | paired d = 0.658 |
| Within-patient, model-discordant | 31 | paired d = 0.614 |

The paired design detects a smaller effect with a third of the patients, and the
comparison understates its advantage: pairing removes between-patient anatomical
variance, which is the dominant variance component for chamber volume.

**The between-patient design cannot answer the question at this n.** Translating d = 0.849
through a published between-subject CoV for cardiac chamber volume on non-gated CT
(20–25%) puts its detection floor at **17–21% of mean volume**. A between-patient null
would therefore be consistent with a 15% acquisition-driven volume shift — which is
exactly the magnitude a radiology reader cares about. Reported as a null, it would be
uninterpretable; reported as a finding, it would be unsupportable.

The paired comparison has no such problem, because the same heart is measured twice.

## Recommended design change

Make the **within-patient paired comparison the primary endpoint** and the
between-patient regression a secondary, explicitly under-powered, descriptive analysis.

Primary: in the 27 thickness-discordant and 31 model-discordant patients, does the
segmentation return a different volume for the *same heart* scanned under different
acquisition parameters? That is a direct measurement of segmentation instability with the
patient as its own control, it needs no body-surface-area normalization (the plan's
Phase 2 step 1 becomes unnecessary for the primary analysis), and it is the quantity a
Meta-Radiology reader would cite.

Confounder to pre-register: interval between studies. Real anatomical change over a
median 46 days is a competing explanation. 25 of 35 repeat patients are within 90 days;
pre-specify the ≤90-day subset as primary and use interval as a covariate to show the
effect does not scale with time — if it does, it is disease progression, not the model.

## Correction to a saved artifact

The first version of `dicom_header_audit.csv` was built with a per-column
`groupby().first()`, which skips missing values **independently per column**. For the 7
patients whose primary series has no `ContrastBolusAgent`, it silently filled the field
from a *different* series of the same patient, so that CSV showed 28 populated values
where the primary series carry 21. Only `ContrastBolusAgent` was affected — it is the one
sparse column; vendor, kVp, tube current, thickness, spacing, kernel, software version and
series UID were identical under a whole-row selection (0 differing rows of 81).

The table is regenerated by selecting a whole row per patient. The 21/81 figure in the
gate table, `phase0_gate.json`, and the prose was computed pre-groupby and was correct;
the artifact was the defective copy. Coverage is 21/81 = 25.9%.

Worth noting for Phase 1: in all 7 cases a sibling series *does* record contrast, so the
patient likely received it — the series that would actually be segmented just doesn't say
so. That strengthens rather than weakens the three-level-factor decision.

## Deliverables

- `dicom_header_audit.csv` — 81 patients × 28 fields, one primary series each
- `dicom_headers_raw.csv` — all 573 series, full tag set
- `phase0_gate.json` — machine-readable gate verdicts
- `phase3_testretest.csv`, `thickness_discordant_pairs.csv` — paired-design cohort
- `fig_phase0_gate.png` — the six gate panels
- `PRESPECIFICATION_IMAGING.md` — kill criteria, fixed before Phase 1

## What has not been done

No segmentation has been run. No model has been downloaded. TotalSegmentator weights are
not pinned yet — the plan's requirement to log model version and weights hash applies to
Phase 1 and is carried into the pre-specification, not satisfied here.
