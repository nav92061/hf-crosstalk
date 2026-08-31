# Pre-specification: segmentation-stability arm — VERSION 2

**Supersedes v1 (2026-08-29).** v1 named a cross-study paired endpoint as primary. The
intersection guardrail halted it at n=5 before any segmentation ran; that halt is reported,
not hidden. This version re-issues the arm on a within-study reconstruction design and adds
a second endpoint that uses the paper's own receptor panel.

**Fixed before any segmentation is run.** No segmentation output of any kind exists. No
volume has been downloaded. TotalSegmentator is not yet pinned.

## 0. What changed from v1 and why

| | v1 | v2 |
|---|---|---|
| Primary | cross-study paired, n=5 | within-study reconstruction, n=58 kernel / n=28 thickness |
| Contrast | 3-level factor | irrelevant — identical by construction |
| Interval | ≤90 days, covariate | zero by construction |
| Unit of inference | **undeclared** | patient, explicitly |
| Headline statistic | p-values, Bonferroni ×5 | % of mean with CI; p demoted |
| Uses the paper's findings | no | yes — endpoint 2 |

The v1 endpoint is retained in the manuscript **as a reported halt**, with its funnel
(35 → 31 → 17 → 5). It is not re-run under relaxed criteria.

## 1. Cohort — closed

59 patients with ≥1 confirmed same-acquisition reconstruction pair, drawn from the 81
seg-usable patients of the Phase 0 audit. Kernel arm 58, thickness arm 28 (overlapping).
193 confirmed pairs, 217 distinct series, 20.2 GiB.

Same-acquisition confirmation: identical `FrameOfReferenceUID` **and** identical
`AcquisitionTime` or `AcquisitionNumber`. 234/234 tags fetched, 0 errors. Pairs failing
this test (53 of 246) are excluded and not reconsidered.

Cohort closed at this file's timestamp. No patient may be added or removed on the basis of
a segmentation or expression result.

## 2. Unit of inference — declared

**The patient is the unit of inference. The pair is the unit of measurement.**

This is not a formality here. In the kernel arm, 151 pairs come from 58 patients: median 2
pairs per patient, maximum 10, and 24 patients contribute exactly one. Analysing pairs as
independent would weight one patient ten times another and inflate significance through
pseudo-replication.

Pre-specified aggregation: **one value per patient per structure**, taken as the median
across that patient's pairs. Confirmatory analysis: linear mixed model with patient as a
random intercept, fitted on all pairs.

**Both are reported.** Where the median-aggregate and the mixed model disagree, the
disagreement is the result and is reported as such — not resolved by choosing the more
favourable one. The spread across aggregation rules is a reportable quantity.

## 3. Structures — five, fixed

**Amended 2026-08-29, before any segmentation ran.** v2 named LV, RV, LA, RA and
myocardium. Inspection of the TotalSegmentator 2.18.0 package shows those five come from
task `heartchambers_highres` (task_id 301), which the package marks `commercial: True` and
gates behind a licence. The open `total` task (92 structures) contains whole-heart but no
chambers. Rather than take a licence dependency, the structure list is re-specified on the
open model:

**heart, aorta, pulmonary vein, superior vena cava, inferior vena cava.**

Still exactly five, so the multiplicity correction is unchanged. Two consequences worth
stating rather than hiding: the readout is whole-heart volume, not per-chamber, which a
cardiac reader may find less familiar; and the great vessels are thin structures whose
apparent volume is *more* sensitive to slice thickness than the heart is, which sharpens
the reconstruction question rather than weakening it.

This amendment is recorded because it was made before segmentation. Any change to this
list after results exist would be a different matter and is not permitted.

All five reported regardless of outcome. Per-series outputs recorded: volume (mm3),
n_voxels, mean intensity, centroid, bounding box, runtime, exception status — taken from
TotalSegmentator's own `statistics.json` rather than recomputed downstream.

**Border-truncated structures.** TotalSegmentator reports volume 0.0 for any structure
touching the image boundary unless `--stats_include_incomplete` is set. Thoracic staging
CT routinely clips the great vessels, and four of the five structures here are vessels, so
the flag is set and volumes are computed for truncated structures. Consequence, declared
before results: these are **partial-field volumes, not anatomical totals**, and are not
reported as absolute anatomy. The paired comparison remains valid because both members of
a pair are reconstructions of the same acquisition and therefore share z-coverage; a
per-structure count of border-touching cases is reported alongside the effect sizes.

## 4. Endpoint 1 (primary) — reconstruction sensitivity

**Does the segmentation return a different volume for the same acquisition reconstructed
differently?**

- Kernel arm: n=58 patients, 151 pairs
- Thickness arm: n=28 patients, 108 pairs
- Reported separately. No pooling.

Estimand: **within-patient percentage difference in structure volume**, reported as % of
that structure's mean with a 95% CI. Thin-slice / sharper-kernel series is the reference.

### Statistical framing — deliberately not p-value-led

A null of "reconstruction has no effect" is false a priori: the images differ, the model is
deterministic, so the masks differ. With enough pairs it will be rejected regardless, and
rejection carries no information.

**The 10% kill criterion is the test.** Lead with magnitude and CI. Bonferroni-corrected
p-values across the five structures are reported in a footnote for completeness, not as the
inferential basis. A reviewer asking "what did you think you were testing" should find the
answer in the first sentence of the results.

## 5. Endpoint 2 (co-primary) — radiogenomic association stability

**Declared before segmentation, not added afterward.** This is what connects the arm to the
paper rather than leaving it a standalone QA note.

For each patient, correlate chamber volume against tumour expression of the receptors from
the existing screen — NPR1 and the 100 receptor units in the 160 floor-passing axes of
`floor_audit.csv` — then **recompute the same correlations from the other reconstruction of
the same acquisition** and ask how many associations survive.

Feasibility verified before writing this: **all 59 recon-arm patients have primary-tumour
RNA-seq** in `tcga_expression_subset.csv.gz` (1436 genes × 9538 aliquots); NPR1 and the
PTN/MDK receptor set are all present.

Pre-specified metrics, computed per structure:
1. **Sign-flip rate** — fraction of receptor associations reversing direction between reconstructions
2. **Significance-crossing rate** — fraction crossing p<0.05 in either direction
3. **Rank instability** — Spearman correlation between the two reconstructions' receptor rankings

At n=58 the critical |r| is 0.258 and roughly 25 of 500 tests are expected significant by
chance alone. **Abundant spurious associations help rather than hurt**: the claim is
instability, not discovery. No association found here is offered as a biological finding,
in either direction — a point the manuscript must state explicitly so a reader cannot mine
the table.

**Kill criterion for endpoint 2:** if fewer than 10% of receptor associations change sign or
significance across reconstructions, the arm reports that radiogenomic associations are
stable to reconstruction, and says so.

## 6. Endpoint 3 — segmentation failure rate

Failure/exception rate per series by vendor and thickness, reported whether or not
endpoints 1 and 2 return anything. Pre-specified now so it cannot serve as a post-hoc
salvage. Model version and weights hash logged; an unpinned run is discarded and re-run.

## 7. Kill criterion — endpoint 1

> **If no reconstruction parameter shifts any of the five structures by more than 10% of
> that structure's mean volume, with the 95% CI excluding 10%, the arm reports a null.**

Threshold set below the design's detection floor (kernel MDE d=0.45, thickness d=0.65) and
above typical inter-observer variability for CT chamber volumetry.

## 8. Scope — narrowed from v1, and stated plainly

This design measures **reconstruction sensitivity, not acquisition sensitivity.** Both
series come from one irradiation event, so it cannot speak to scanner-to-scanner, vendor, or
protocol variation. v1 claimed the broader thing; v2 does not, and the manuscript title and
abstract must match v2.

Unchanged from v1: n is small, one pretrained model, one archive. These are cancer patients
imaged for staging, not heart-failure patients. Endpoint 2 uses the paper's receptor panel
but **cannot recover the heart-failure ligand side** — the ligand source is explanted
failing hearts and no imaging cohort here can supply it. Relevance to the paper improves;
it is not solved, and the manuscript should not imply it is.

## 9. Withdrawal precedent

The v1 endpoint was withdrawn against a pre-specified criterion, before compute was spent.
That is the second such withdrawal in this project after the seven pleiotrophin axes. If
endpoint 1 or 2 triggers its kill criterion, the same applies.
