# PRESPECIFICATION — adversarial stress test of the PTN axes

Written and saved **before any analysis cell in this phase was executed**. Thresholds,
directions, null models and seeds below are fixed. A result that misses a threshold is recorded
as a FAIL and reported as such.

## What is being tested and why

The corrected specificity screen (§2.10) returned 7 surviving axes, all with PTN as ligand,
against 0 under the original criterion. Two things could make that a false positive:

1. **PTN was admitted by a relaxation I introduced.** Criterion A3 was loosened from "heart ranks
   first among 52 GTEx tissues" to "heart ranks in the upper half". PTN sits at rank 23 of 52.
   A hit that appears only because a criterion was relaxed is exactly what a relaxation
   manufactures.
2. **The receptor side may not discriminate.** PTN's receptors are syndecans, nucleolin and
   integrin alpha-V — broadly expressed adhesion and proteoglycan genes. If most receptors pass
   the tumour-enrichment test regardless of ligand, "7 axes" describes receptor promiscuity, not PTN.

**Base rate measured before setting thresholds** (from `receptor_tumor_vs_normal.csv`, the 81
floor-passing receptor subunits): 60 of 81 receptors (**74.1%**) are enriched versus matched
adjacent normal in at least one tumour type. Pass rates: 55.6% at >=2 types, 23.5% at >=4,
11.1% at >=6. PTN's seven receptors pass in a mean of 3.86 types against an all-receptor mean of
2.52 (medians 4 and 2). PTN's receptors are therefore somewhat above average but **well inside**
the distribution — this is the concern the tracks below must resolve, and the numbers are recorded
here so the thresholds cannot be chosen to flatter the result.

**Prior stated in advance:** the expected outcome is that PTN fails Track 1 or Track 2. A failure
is a good result — it converts the parent conclusion from "we used the wrong filter" to "we used
the wrong filter, corrected it, and the conclusion held."

---

## Global prohibitions

- **P1.** No cross-platform absolute comparison (TCGA RSEM / GTEx TPM / RPKM / normalised counts).
  Cross-dataset comparisons are ordinal only. Carried over from the parent pre-specification.
- **P2.** No threshold adjustment after seeing a result. No near-miss promotions.
- **P3.** No substituting a different cohort, endpoint or statistic because the pre-specified one
  gave an unwelcome answer. A cohort that cannot be obtained is reported as not obtained.
- **P4.** Literature records must be retrieved and verified. No citation written from memory.

---

## Overall retention rule, fixed now

**PTN is retained as a candidate axis only if it survives all three tracks.** Failing any one
track is sufficient to withdraw the seven-axis claim as a positive finding. Partial survival is
reported as a failure with the surviving components named.

---

## TRACK 1 — Does the receptor side discriminate at all?

### T1.1 Receptor-side base-rate audit
Compute, for every one of the 81 floor-passing receptors, the number of tumour types in which it
satisfies the unchanged receptor criterion, independent of any ligand.

- **KILL CRITERION:** if **>= 60%** of floor-passing receptors pass in at least one tumour type,
  the seven-axis count is declared **UNINFORMATIVE** as evidence for PTN specifically. The
  measured base rate is already 74.1%, so on current evidence this criterion is expected to
  trigger; the audit exists to confirm it under the exact per-axis logic used in §2.10 and to
  quantify how far PTN's receptors sit above the base rate.
- Report the full distribution, not just the summary.

### T1.2 Matched-random null on the corrected axis set
Substitute PTN's seven receptors with expression-decile-matched random receptor sets of the same
size, 1,000 permutations, **seed 20260729**. Statistic: total number of (receptor, tumour type)
enriched-and-floor-passing pairs across the set, observed value 27.

- **PASS:** observed exceeds the 95th percentile of the null.
- **FAIL:** otherwise. A FAIL means PTN's receptor set is unremarkable for its expression class.

### T1.3 Receptor family independence
Test whether SDC1, SDC2, SDC3, SDC4, NCL, ITGAV and PTPRZ1 behave as independent receptors or as
one correlated block, using Spearman correlation of per-tumour-type expression across the 30 types
and hierarchical clustering.

- **Reported, not a kill criterion**, with one exception: if the seven receptors form a single
  block at **mean pairwise rho >= 0.7**, the seven axes are reported as **one observation, not
  seven**, and the count is presented that way in the manuscript.

**Track 1 verdict:** PASS requires T1.2 to pass AND T1.1 not to trigger the uninformative
criterion.

---

## TRACK 2 — Is PTN a cardiac signal or a generic injury signal?

### T2.1 Induction outside the heart
Process public bulk cohorts for non-cardiac organ failure or fibrosis (candidate contexts:
idiopathic pulmonary fibrosis, cirrhotic liver, chronic kidney disease / renal fibrosis) using the
`geo-bulk-de` skill, with each cohort's normalisation kept separate (P1). At least **two**
non-cardiac cohorts must be successfully processed for the track to return a verdict; if fewer are
obtainable, the track is reported as INCONCLUSIVE and, by the overall retention rule, PTN is not
retained on unverified generality.

- **KILL CRITERION:** PTN is a generic injury signal — and the heart is not a distinctive source —
  if PTN is significantly induced (FDR < 0.05, same direction) in **>= 2** non-cardiac
  failure/fibrosis cohorts with a median effect size **>= 50%** of the cardiac effect
  (cardiac reference: min replicated log2FC 1.5023).
- **PASS:** PTN induction is absent or substantially weaker outside the heart.

### T2.2 Cell of origin in failing myocardium
Identify a public human heart-failure single-nucleus or single-cell RNA-seq dataset and determine
which cell population expresses PTN and whether that population's PTN rises in failure. Marker-level
or pseudobulk analysis is acceptable on an 8 GiB machine; state which was used.

- **KILL CRITERION:** if PTN is expressed predominantly by fibroblasts and its induction tracks a
  general fibrotic programme rather than a cardiomyocyte-specific one, PTN is reported as a
  **fibroblast wound-healing signal**, which does not make the failing heart a distinctive
  systemic source. Combined with a T2.1 failure this is decisive; alone it is a strong caveat that
  must appear in the abstract.
- Cardiomyocyte-dominant expression would be the result that supports a cardiac-source reading.
- If no suitable dataset is obtainable, report as such — do not substitute mouse data for a human
  claim without labelling it.

**Track 2 verdict:** PASS requires T2.1 to pass; T2.2 is reported alongside and can independently
force the fibroblast caveat.

---

## TRACK 3 — Has plasma PTN been measured in heart failure?

Search the literature connectors for circulating pleiotrophin in heart failure or cardiac injury,
and separately for PTN as a secreted factor in cancer. Every record retrieved and verified (P4),
with identifiers recorded and each classified as supporting, contradicting, or irrelevant to a
cardiac plasma source.

- This track **cannot kill PTN on its own** and is not part of the retention conjunction in the
  sense of requiring a positive: an absence of measurement is an informative gap, not evidence
  against.
- **However:** if published data show plasma PTN is **not** elevated in heart failure, that is
  a direct contradiction of the endocrine premise and **is** a kill criterion.
- If plasma PTN has never been measured in heart failure, record that explicitly as the gap that
  the whole axis rests on.

---

## Seeds

| Purpose | Seed |
|---------|------|
| T1.2 matched-random receptor null | 20260729 |
| Any bootstrap in this phase | 20260729 |
| Inherited, unchanged | parent permutation 20260727, corrected-phase 20260728 |

## Multiplicity

Pre-specified primary tests: **3** — T1.2 null, T2.1 cross-tissue induction, T3 contradiction
check. T1.1, T1.3 and T2.2 are characterisations with stated interpretive consequences rather than
significance tests.

## Recording

Verdicts append to `validation_summary.csv` continuing the existing 41 rows, same columns. A FAIL
is recorded with equal prominence to a PASS. Figures are built **only if PTN survives all three
tracks**; if it does not, the manuscript reports the withdrawal in text.
