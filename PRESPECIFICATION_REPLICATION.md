# Pre-specification: protein-level falsification and cross-platform replication

**Status: fixed before any of the analyses below were run.** Written after a
bias audit (`bias_audit_protein.csv`) and before any tumour type other than
ccRCC was touched. Anything already observed is declared in §0.

The purpose of this document is to make the two tracks falsifiable in advance,
so that a negative cannot be softened afterwards and a positive cannot be
assembled from whichever comparison happened to clear significance.

---

## 0. What has already been seen (declared, not blinded)

Verifying the machine-level controls required computing on real data. The
following were therefore observed *before* these criteria were fixed, and
ccRCC is **not a blind test set**:

| Quantity | Observed value |
|---|---|
| ccRCC NPR1 complete-case log2-ratio difference (tumour − normal) | **−0.153** |
| ccRCC NPR1 missingness, tumour vs normal | 16.4% vs 19.0% (Fisher p=0.704) |
| ccRCC GAPDH log2-ratio, tumour vs normal | 0.295 vs −0.427 |
| ccRCC plex × arm uniformity | χ²=15.69, dof=22, p=0.831 |
| ccRCC tumour/normal pairs sharing a plex | 84 of 84 |

Consequences, binding:

1. **Thresholds below were chosen without reference to −0.153.** They are
   inherited from the transcript-level analysis or set on statistical grounds.
2. **The primary test set is the seven tumour types not yet examined**: OV,
   BRCA, PAAD, UCEC, LUAD, HNSC, GBM. ccRCC is reported alongside them but
   flagged as unblinded, and the headline verdict is computed both with and
   without it.
3. **GAPDH is not invariant between arms** (0.295 vs −0.427). The transcript
   floor (0.005×GAPDH) is therefore **not portable to protein**. The protein
   inclusion criterion is peptide evidence, §2.1 — not a housekeeping ratio.

---

## 1. The claim under test

The paper's central negative, at transcript level, in two independent forms:

- **Unit-matched adjacent normal**: NPR1 is depleted in tumour in **all 15**
  tumour types with ≥10 adjacent normals, median background-corrected
  log₂FC **−1.91**, all FDR<0.05.
- **GTEx rank-based**: NPR1 is enriched in **0 of 30** types, depleted in 12,
  comparable in 18.

Both are mRNA. The protein-level test asks whether that depletion is a
property of the transcript or of the receptor.

---

## 2. Track A — NPR1 protein in tumour vs matched normal (CPTAC)

### 2.1 Inclusion: what makes a study testable

A (study, gene) cell is **testable** only if all hold:

- the gene appears in the study's `quantDataMatrix` gene index;
- **≥ 10 tumour and ≥ 10 normal aliquots** carry a non-missing value;
- **≥ 60%** of aliquots in each arm are non-missing.

A cell failing any of these is reported **UNTESTABLE**. It is never reported as
a negative, and never contributes to a count of types showing depletion. Genes
absent from a study matrix (in ccRCC: NPR2, NPR3, SDC3, PTPRZ1) are untestable
by definition.

This criterion exists because absence of detection and absence of protein are
indistinguishable in a TMT matrix, and conflating them would manufacture the
paper's own conclusion.

### 2.2 Primary test: paired within plex

Because all 84 ccRCC pairs share a plex with their partner, the primary test is
a **paired test on within-plex tumour−normal differences**, one pair per case:

- Wilcoxon signed-rank on the paired differences (does not assume normality);
- BH-FDR across tumour types within the gene;
- effect size reported as the median paired difference in log₂ ratio units.

Pairing within plex is what removes TMT batch. It is preferred over any
statistical batch correction because the shared reference channel cancels
arithmetically.

**If a study has no within-plex pairs**, the fallback is an unpaired
Mann-Whitney with plex as a blocking factor, and the result is flagged as
`batch_corrected` rather than `batch_free`.

**Single-arm plexes.** Not every plex is mixed: in ccRCC, 21 of 23 contain both
arms, but plexes 19 and 20 hold 9 tumours and 0 normals. Such plexes contribute
no pairs and drop out of the paired test automatically. They must **not** enter
an unpaired analysis without plex as a blocking factor, because their tumours
would otherwise supply an arm difference that is partly a plex difference. Any
unpaired sensitivity analysis (§2.5.4) reports how many aliquots sit in
single-arm plexes.

### 2.3 Direction and the three possible verdicts

Fixed now, with no discretion later:

| Outcome | Criterion | Meaning |
|---|---|---|
| **REPLICATES** | median paired difference < 0 and FDR<0.05 in **≥ 5 of 7** blind types | The transcript negative holds at protein level. Paper's conclusion survives an independent molecule. |
| **REVERSES** | median paired difference > 0 and FDR<0.05 in **≥ 3 of 7** blind types | **The paper's central negative is wrong.** NPR1 protein is gained where its mRNA is lost; the screen must be reopened. |
| **INCONCLUSIVE** | anything else, including mixed directions or too few testable types | Reported as inconclusive. Not spun as support for either side. |

The asymmetry (5 for replication, 3 for reversal) is deliberate: it is easier
to overturn the paper than to confirm it. A conclusion should be cheap to
falsify and expensive to defend.

### 2.3b Which types permit a direct transcript-vs-protein contrast

The 15 types with a unit-matched transcript baseline and the 8 CPTAC proteome
types overlap in **five**: BRCA, HNSC, KIRC, LUAD, UCEC. Of the seven blind
types, **four** have a transcript baseline (BRCA, UCEC, LUAD, HNSC).

**GBM, OV and PAAD have no unit-matched transcript baseline** — TCGA has too few
adjacent normals in those types. A protein result there stands on its own and is
**not** evidence of agreement or contradiction with the transcript analysis. It
still counts toward the §2.3 tallies, because the question "is NPR1 protein
depleted in tumour" is well posed regardless of whether mRNA was measurable in
the same type; but the manuscript must not describe those three as confirming or
refuting the transcript result.

Where both exist, report the transcript log₂FC and the protein paired difference
side by side, and state the sign agreement per type.

### 2.4 Mandatory base rate — the guard that killed PTN

Before NPR1's result is interpreted, compute over **all genes** quantified in
each study the fraction depleted in tumour at FDR<0.05 by the same test.

- If that base rate **exceeds 40%**, then "NPR1 is depleted" is a statement
  about the assay or about tumour tissue generally, **not about NPR1**, and the
  NPR1 result is reported **UNINFORMATIVE** regardless of its own p-value.
- NPR1's percentile within the genome-wide effect distribution is reported
  alongside its FDR, always.

This is the same guard that withdrew the seven PTN axes, where 74% of receptors
passed the receiver criterion with no ligand involved.

### 2.5 Sensitivity analyses, fixed in advance

Each is reported whether or not it agrees with the primary:

1. `unshared_log2_ratio` (unique peptides only) in place of `log2_ratio`.
2. Dropout-aware: treat missing as left-censored at the per-aliquot minimum
   observed value, and re-test. Complete-case remains primary.
3. ccRCC excluded (it is unblinded).
4. Unpaired test ignoring plex, to show what batch handling contributes.

### 2.6 Track A kill criterion

If fewer than **4 of the 7 blind types** are testable under §2.1, Track A is
**abandoned as underpowered** and reported as such. No verdict is issued on the
strength of two or three cohorts.

---

## 3. Track B — PTN receptor panel at protein level

Panel: SDC1, SDC2, SDC3, SDC4, NCL, ITGAV, PTPRZ1 (the seven axes surviving
the corrected-specificity phase), plus NPR2 and NPR3.

- Same testability, pairing, and FDR rules as Track A.
- **Null**: abundance-matched random receptors, matched on decile of mean
  protein abundance within the study, 300 draws, seed 20260727. An unmatched
  null is not acceptable: in the transcript phase an unmatched (naive) pool gave
  ρ=0.323 against ρ=0.363 for the decile-matched pool at the same floor setting,
  so an unmatched pool understates the null and would wrongly credit a score with
  specificity. (The separate 0.363 → 0.107 contrast is the effect of the
  expression floor, not of matching, and is not the justification for this rule.)
- **Kill criterion**: if the panel's count of enriched-in-tumour receptors does
  not exceed the **95th percentile** of the matched null, the panel is
  **withdrawn at protein level**, exactly as it was at transcript level.

PTN itself: reported if quantified, but CPTAC measures the tumour, not plasma.
Absence of PTN protein in a tumour matrix is **uninformative** about a
circulating cardiac ligand and will not be presented as evidence either way.

---

## 4. Track C — cross-platform HF signature replication

Cohorts: **GSE5406** (n=210, Affymetrix GPL96), **GSE1869** (n=37, GPL96),
**GSE57345** (n=319). The existing discovery and validation pair
(GSE116250, GSE141910) **both ran on GPL16791** — a shared-platform weakness
found during this audit and not previously reported.

### 4.1 Criteria

| Quantity | Threshold to call the signature platform-independent |
|---|---|
| Directional concordance, genes significant in both | **≥ 75%** |
| log₂FC Spearman ρ, genes significant in both | **≥ 0.40** |
| NPPA direction | up in failing, FDR<0.05, in ≥2 of 3 cohorts |
| NPPB direction | up in failing, FDR<0.05, in ≥2 of 3 cohorts |
| Label-permutation control | **0** significant genes at FDR<0.05, every cohort |

The concordance threshold (75%) is below the 0.881 already observed between the
two RNA-seq cohorts, because cross-technology comparison is expected to be
noisier. It is set on that reasoning, not on any microarray result.

### 4.2 Two correlations, never conflated

A prior audit finding: ρ across all shared genes (0.541) and ρ among genes
significant in both (0.759) are different quantities. **Both are reported
separately for every cohort pair.** Quoting one as the other is a defect.

### 4.3 Track C kill criterion

If concordance falls below 75% or ρ below 0.40 in **2 or more** of the three
cohorts, the signature is declared **platform-dependent**, and every downstream
claim resting on the 5,036-gene signature is re-scoped in the manuscript.

---

## 5. Analysis hygiene, binding on all tracks

- **float64 throughout.** Verified: row-order invariance 0.000e+00,
  sample-order 8.9e-16, repeat-matmul bit-identical, float32 drift 1.3e-07 log₂
  units against effect sizes ~1e-1.
- **Seeds pinned**: 20260727 for every resampling procedure.
- **No optional-stopping.** All tumour types and all cohorts specified above are
  analysed and reported. Dropping one requires a §2.1 or §2.6 reason recorded in
  `validation_summary.csv`.
- **Verdict vocabulary**: PASS / FAIL / CORRECTED / UNTESTABLE / UNINFORMATIVE /
  INCONCLUSIVE. A partial pass is a FAIL; `conjunction_verdict()` enforces this
  in code.
- **Figures only for surviving results.** If nothing survives, no figure is
  built, and the absence is stated.

---

## 6. What would change the paper's conclusion

Stated plainly so it cannot be renegotiated after the fact:

- **NPR1 protein enriched in tumour in ≥3 of 7 blind types (FDR<0.05)** →
  the central negative is refuted; the hypothesis reopens at protein level.
- **Base rate >40%** → the transcript-level NPR1 result is reinterpreted as an
  assay/tissue property rather than a receptor-specific one, which weakens the
  paper's mechanism even though it does not support the original hypothesis.
- **Signature platform-dependent in ≥2 cohorts** → the ligand list itself is
  partly a platform artefact, and §3 of the paper is re-scoped.
- **PTN panel exceeding its matched null at protein level** → the withdrawal is
  reversed and the axes return as candidates.

Any other outcome leaves the paper's conclusions standing, and that will be
stated without inflation.
