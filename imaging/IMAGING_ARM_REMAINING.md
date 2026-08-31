# What imaging work is left

The pre-specified cross-study endpoint is dead (n=5). A different design in the same data
is alive and cleaner. This document is what remains, in order.

## The pivot

**Dead:** compare a patient's chamber volumes across two *separate studies* on different
scanners. Requires acquisition discordance + ≤90 days + verified contrast concordance
simultaneously → **5 patients, MDE d=1.53 ≈ 15% of mean**, above the 10% kill threshold.
The endpoint cannot reach its own criterion.

**Alive:** compare *two reconstructions of the same acquisition*. Many studies here
reconstructed one irradiation event at multiple slice thicknesses and kernels.

| Design | n | pairs | MDE (d) | Confounds |
|---|---|---|---|---|
| **Recon, kernel differs** | **58** | 151 | **0.45** | none — identical raw data |
| **Recon, thickness differs** | **28** | 108 | **0.65** | none — identical raw data |
| Cross-study, contrast assumed | 17 | 24 | 0.83 | ~36% silently contrast-discordant |
| Cross-study, contrast verified | 5 | 7 | 1.53 | contrast, interval, anatomy |

Confirmed by identical `FrameOfReferenceUID` plus identical `AcquisitionTime` or
`AcquisitionNumber` — 234/234 tags fetched, 0 errors, 193 of 246 pairs confirmed.

Why it is cleaner: contrast, cardiac phase, anatomy, patient position and elapsed time are
identical **by construction**. The confound that killed the cross-study endpoint cannot
arise. No BSA normalization needed. Interval is zero, so disease progression is excluded.

**Scope cost, which must be stated:** this measures *reconstruction* sensitivity, not
scanner-to-scanner sensitivity. It cannot speak to vendor or protocol variation. The claim
narrows from "acquisition shift" to "reconstruction shift" and the manuscript must say so.

## Work remaining

### 1. Re-write the pre-specification — DONE
`PRESPECIFICATION_IMAGING.md` **v2** is issued and supersedes v1. Reconstruction design
primary, kernel (n=58) and thickness (n=28) reported separately, v1 endpoint demoted to a
reported halt. Three changes beyond the pivot, all fixed before segmentation:
- **Unit of inference declared**: patient, not pair. 151 pairs come from 58 patients
  (median 2, max 10, 24 with exactly one), so treating pairs as independent is
  pseudo-replication. Median-per-patient aggregation plus a patient-random-effect mixed
  model, **both reported**, with disagreement between them reported as a result.
- **Magnitude-led, not p-led**: "reconstruction has no effect" is false a priori, so
  rejection carries no information. The 10% criterion is the test; Bonferroni goes in a
  footnote.
- **Endpoint 2 added** (below).

### 2. Pin TotalSegmentator
Weights hash, model version, container digest. The plan requires it and Phase 0 did not do
it. An unpinned run is discarded.

### 3. Download volumes — the first real data transfer

**Corrected.** An earlier version of this document said "~118 series," which was 234/2 —
wrong, because series are shared across multiple pairs and cannot be halved. Counted
properly from `recon_pairs_confirmed.csv`:

| Scope | pairs | distinct series | total size |
|---|---|---|---|
| All confirmed | 193 | **217** | **20.2 GiB** |
| Kernel arm only | 151 | 210 | 19.0 GiB |
| Thickness arm only | 108 | 121 | 12.9 GiB |

The two arms overlap by 114 series, so running both costs only 7 more series than kernel
alone — the arm choice is not a download-cost decision.

Median series is 70 MiB. Against the 2.8 GiB free-disk ceiling from prior feasibility work,
**staged transfer is mandatory, not optional**: download one series → segment → write masks
and volumes → delete the volume → next. Peak disk stays under ~200 MiB. At 217 series this
is the longest-running step in the arm and should be checkpointed so it can resume.

### 4. Run segmentation, log failures
Per-series: 5 structure volumes, surface areas, centroids, runtime, exception status.
Failure rate by vendor and thickness is a pre-specified endpoint in its own right.

### 5. Endpoint 1 — reconstruction sensitivity
Per-patient median of paired log-ratios, kernel and thickness arms separately, effect
sizes as % of mean with CIs, all five structures regardless of outcome. Mixed model as
confirmatory. Bonferroni footnoted.

### 5b. Endpoint 2 — radiogenomic association stability
Correlate chamber volume against tumour expression of the screen's own receptor panel
(NPR1 + the 100 receptor units in the 160 floor-passing axes), then recompute from the
other reconstruction of the same acquisition and measure sign-flip rate, significance-
crossing rate, and rank instability.

**Feasibility verified:** all 59 recon-arm patients have primary-tumour RNA-seq
(1436 genes × 9538 aliquots); NPR1 and the PTN/MDK receptor set all present. At n=58 the
critical |r| is 0.258 and ~25 of 500 tests are expected significant by chance — abundant
spurious hits are the point, since the claim is instability rather than discovery.

This is the step that stops the arm being a standalone QA note. It costs a few lines on
top of segmentation output already being computed, and it must be declared now: added
after segmentation it is a post-hoc rescue.

### 6. The guardrail instance
The intersection check that halted the cross-study endpoint is the paper's demonstrable
safeguard — it fired before compute was spent, on the authors' own data. Write it up as a
result, with the funnel figure.

## Two things deliberately not done

**Dropping expression-matching.** 40 more patients are available (LUAD 38, LUSC 2) at zero
data cost, since the imaging arm never uses the RNA-seq. Projected gain to the strict
cross-study intersection: **~2 patients**. The binding constraint is the contrast field,
not cohort size, so this does not rescue anything and is not worth the cohort-definition
churn. It would roughly double the reconstruction cohort, which is worth revisiting only
if the n=28 thickness arm proves underpowered in practice.

**Widening beyond the four TCGA collections.** Same reasoning; defer unless the
reconstruction design fails.

## The relevance limit that remains

Endpoint 2 uses the paper's receptor panel and its expression matrix, which closes most of
the drift. It does **not** close all of it: the cohort is still LUSC 26 / LUAD 17 / ESCA 15
/ SARC 1, selected by where the heart falls in the field of view, and no imaging cohort
here can supply the heart-failure ligand side, whose source is explanted failing hearts.
The honest claim to the editor is that the arm now tests the stability of the paper's own
receptor panel under reconstruction change — not that it validates the biology.

## Deliverables

- `phase1_intersection_gate.json` — the halt, the cross-tab, the rescue design
- `recon_pairs_confirmed.csv` — 246 pairs with same-acquisition confirmation flags
- `fig_design_pivot.png` — funnel collapse, contrast cross-tab, design comparison
