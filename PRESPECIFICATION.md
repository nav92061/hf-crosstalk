# PRESPECIFICATION — corrected specificity screen and NPR1 exploratory follow-up

Written and saved **before any analysis cell in this phase was executed**. Thresholds,
directions, null models and seeds below are fixed. If a result misses a threshold, it is
recorded as a FAIL and reported as such; thresholds are not adjusted after seeing data.

Provenance of the design: the original screen returned zero surviving axes. Three design
choices were identified that could each have manufactured that negative. One of the three
(organ-mass plasma weighting) was **dropped before pre-specification** because it requires
translation efficiency, secretion rate, clearance half-life and volume of distribution,
none of which are present in the available data or reliably in the literature for these
ligands. It is not tested here and its absence remains a stated limitation.

---

## Global prohibitions

- **P1. No cross-platform absolute comparison.** TCGA RSEM, GTEx TPM, GSE116250 RPKM and
  GSE141910 normalised counts are never placed on a common numeric scale. Comparisons that
  cross datasets are ordinal (rank/percentile) only. This rule already governs the parent
  manuscript; the corrected specificity definition must not violate it.
- **P2. No threshold adjustment post hoc.** No changing a cutoff, direction, covariate set,
  outcome variable or correction method after seeing a result.
- **P3. No outcome substitution.** Overall survival is the pre-specified survival endpoint.
  Progression-free interval is not a fallback if OS fails.
- **P4. Near-misses are fails.** A result at p = 0.06 against a p < 0.05 criterion is a FAIL.

---

## TRACK A — Failure-induced cardiac distinctiveness (corrects the primary analysis)

**Flaw being corrected.** The original ligand filter asked whether a ligand is enriched in
GTEx left ventricle, which is *non-failing* myocardium. A ligand that is unremarkable in the
healthy heart but strongly induced in failure, while unchanged elsewhere, is a distinctive
failing-cardiac source and the original filter discarded it. This is limitation #4 of the
parent manuscript.

**Design constraint.** Because of P1, the corrected definition may not estimate absolute
failing-heart abundance and compare it against GTEx TPM in other tissues. All three criteria
below are computed **within a single dataset each**.

### Criteria (a ligand must satisfy A1 AND A2 AND A3)

| ID | Criterion | Source (single dataset) | Threshold |
|----|-----------|------------------------|-----------|
| A1 | Failure induction, replicated | `hf_de_GSE116250.csv`, `hf_de_GSE141910.csv` | `log2FC > 0` and `FDR < 0.05` in **both** cohorts, and `min(log2FC across cohorts) >= 1.0` (≥2-fold) |
| A2 | Abundance guard, housekeeping-relative | each heart cohort separately | `2^mean_expr(ligand) / 2^mean_expr(GAPDH) >= 0.005` in **both** cohorts |
| A3 | Heart is not an also-ran at baseline | GTEx only (`ligand_tissue_specificity.csv`) | `best_heart_rank <= 26` of 52 tissues (upper half) |

A2 mirrors the housekeeping-relative floor (0.005 × GAPDH) that the parent manuscript already
validates for TCGA, applied inside the heart cohorts rather than across platforms. A3 is
deliberately weaker than the original criterion (which effectively required heart rank 1);
relaxing it is the correction under test.

### Receptor side — unchanged

The receptor criteria are **not** modified: tumour-enrichment versus matched adjacent normal
and the 0.005 × GAPDH TCGA expression floor are applied exactly as in the parent analysis,
from `receptor_tumor_vs_normal.csv` and `crosstalk_network.csv`.

### Pass criterion

- **PASS** if `n_surviving_axes >= 1` (the parent result is 0).
- **FAIL** if `n_surviving_axes == 0`.
- **CRITERION UNINFORMATIVE** if `n_surviving_axes > 20`. A relaxation that admits more than
  20 of the 160 floor-passing axes is not discriminating between ligands and will be reported
  as an uninformative criterion, **not** as a positive result. This guard is set in advance
  precisely so that a permissive fix cannot be read as a discovery.

### Reported alongside, whatever the outcome

Count under the original criterion (expected 0) for contrast; the identity of every ligand
newly admitted by A3; and per-criterion attrition (how many ligands each of A1, A2, A3 removes).

---

## TRACK B — NPR1 as a lost antiproliferative receptor (EXPLORATORY)

**Status label, to travel with every result from this track.** This hypothesis was generated
*by* the refutation, not predicted by the original screen. The original hypothesis was that a
cardiac ligand engages a tumour receptor; this track tests close to the opposite — that tumours
benefit from shedding that receptor. It is reported in a separate, clearly labelled exploratory
section. Presenting it as a confirmed prediction of the original design would be HARKing.

### The conjunction rule (fixed here, before any test)

NPR1 has three sub-analyses. **All three must pass and agree in direction.** Passing one or
two is recorded as a FAIL for the track, reported as a partial result with the failing
conjunct named. Rationale: three independent chances to find a hit is how a negative becomes
a false positive; requiring simultaneous agreement is a materially harder test to pass by chance.

Direction consistency required across all three: low NPR1 ⇒ worse survival, higher
proliferation, and genomic loss rather than transcriptional state alone.

### B1 — Survival

- Data: TCGA-CDR (GDC `1b5f413e-a8d1-4d10-92eb-7c4ae739ed81`), endpoint `OS` / `OS.time`.
- Cohort: the 15 tumour types with confirmed NPR1 depletion in `receptor_tumor_vs_normal.csv`.
- Model: Cox proportional hazards, NPR1 expression z-scored **within tumour type** (avoids P1),
  stratified by tumour type, covariates age at diagnosis and sex.
- Pre-specified direction: **HR < 1** (higher NPR1 → better survival).
- **PASS**: pooled stratified model `p < 0.05` two-sided **and** `HR < 1`.
- Also reported: per-type Cox with BH-FDR across the 15 types; count of types with HR < 1.

### B2 — Proliferation coupling

- Data: `tcga_tme_scores.csv` column `Module11_Prolif_score`; NPR1 from `tcga_expression_subset.csv.gz`.
- Test: Spearman correlation, NPR1 versus proliferation score, computed **within each tumour type**.
- Pre-specified direction: **negative** rho.
- **PASS**: median within-type rho `< 0` **and** at least 8 of 15 types individually
  significant at BH-FDR `< 0.05` with negative rho.

### B3 — Genomic loss versus transcriptional silencing

- Data: TCGA ABSOLUTE allele-specific copy number, `TCGA_mastercalls.abs_segtabs.fixed.txt`
  (GDC `0f4f5701-7b61-41ae-bda9-2805d1ca9781`), confirmed reachable and hg19-based.
- NPR1 locus, GRCh37, from Ensembl GRCh37 REST: **chr1:153,651,113–153,666,468** (ENSG00000169418).
- Per sample: `Modal_Total_CN` of the segment(s) overlapping the locus; a sample is called
  *lost* if NPR1 copy number is below that sample's modal ploidy.
- Null: 1,000 random loci matched on chromosome arm (1q) and segment length distribution,
  **seed 20260728**.
- **PASS**: NPR1 loss frequency exceeds the **95th percentile** of the arm-matched null.
- Methylation is a **conditional secondary** read, not part of the conjunction. The merged
  PanCanAtlas HumanMethylation27/450 beta matrix (GDC `d82e2c44-89eb-43d9-b6d3-712732bf6a53`)
  is reachable, but a probe-to-gene mapping is not currently in hand. If a mapping cannot be
  obtained from an allowlisted source, B3 is scored on copy number alone and that limitation
  is recorded. Absence of methylation data is **not** grounds to relax B3.

### Track B pass criterion

**PASS** only if B1 AND B2 AND B3 all pass with directions as specified. Any other outcome is
a FAIL, reported with the specific conjunct(s) that failed.

---

## Seeds

| Purpose | Seed |
|---------|------|
| B3 arm-matched locus null | 20260728 |
| Any bootstrap in this phase | 20260728 |
| Inherited from parent analysis (unchanged) | permutation 20260727, bootstrap 20260727, DE permutation 42 |

New seed 20260728 is deliberately distinct from the parent analysis seeds so that no result
here can be a re-draw of a previously seen null.

---

## Multiplicity accounting

Pre-specified primary tests in this phase: **4** — Track A axis count (1), B1 pooled Cox (1),
B2 median-rho conjunction (1), B3 null exceedance (1). Per-type analyses inside B1 and B2 are
BH-corrected within their own family and are secondary. The Track B conjunction requirement is
itself the multiplicity control for that track.

## Recording

Verdicts are appended to `validation_summary.csv` as new rows continuing the existing 31, with
the same columns (`n`, `check`, `expected`, `observed`, `verdict`, `bearing`). A FAIL is
recorded with equal prominence to a PASS.
