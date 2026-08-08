---
name: crosstalk-score
description: Score how compatible each receiver context (tumor type, tissue, cell line) is with a set of ligands released by a source tissue, decomposed into ligand availability and receiver receptor capacity. Enforces an absolute expression floor and a matched-random specificity null as code-level validators. Use for reverse cardio-oncology, organ-crosstalk, secretome-to-receptor, or any ligand-receptor compatibility ranking across many receiver contexts.
---

# crosstalk-score

Ranks receiver contexts by their compatibility with a source tissue's secreted
ligands, and — more importantly — refuses to produce a ranking that a
matched-random receptor set could have produced.

## Why the validators exist

A ligand-receptor compatibility score is usually built by z-standardizing
receptor expression across receivers and averaging. That design has a specific
failure mode: **z-scoring a gene that is not transcribed anywhere still yields
large z-values**, because the standardization only sees relative variation. A
receptor sitting at noise level will show apparent "high compatibility" in
whichever receivers have marginally more noise.

This is not hypothetical. On real TCGA data (30 solid tumor types, 188
candidate interactions from a human heart-failure secretome):

| Configuration | Matched-random null ρ (mean) | Fraction of random sets with ρ>0.5 |
|---|---|---|
| No expression floor | 0.363 | 21.5% |
| **Floor at 0.005 × GAPDH** | **0.107** | **1.7%** |
| Unmatched (naive) random pool | 0.323 | 17.0% |

Two lessons are encoded in code as a result:

1. **The floor is what makes the score receptor-specific.** Removing
   below-floor receptors cut the null correlation from 0.363 to 0.107. The
   confound people attribute to "a global expression axis" is substantially
   *caused by* scoring untranscribed receptors.
2. **The null pool must be expression-matched.** A naive random pool gives
   ρ=0.323 and would wrongly flag a specific score as confounded.

`require_expression_floor()` **raises** by default and
`matched_random_null()` matches on expression decile. Guardrails live in
assertions, not in prose that a reader may skip.

## Functions

```python
housekeeping_floor(expr, reference="GAPDH", fraction=0.005) -> float
```
Unit-agnostic absolute floor. Works for RSEM normalized counts, TPM or FPKM
because it is defined relative to a housekeeping gene's median rather than as a
hardcoded number. **Never hardcode "2 TPM"** — that is meaningless in RSEM
normalized counts, where GAPDH's median is ~69,000.

```python
receptor_capacity(expr, interactions, rule="min") -> DataFrame  # interactions x receivers
```
Aggregates subunits into per-interaction receptor level. `rule="min"` is the
default and encodes obligate-subunit semantics: a heteromeric receptor cannot
signal if any required subunit is absent. `"geomean"` is the sensitivity
alternative.

```python
require_expression_floor(capacity, floor, raise_on_violation=True) -> (passing, audit)
```
**MANDATORY VALIDATOR — raises `AssertionError`.** Returns an audit table so
dropped interactions are disclosed rather than silently removed.

```python
ligand_availability(de_table, ligand_genes, tier_map=..., tier_weights=...) -> Series
matched_random_null(expr, interactions, availability, background_genes, n_iter=1000, seed=0) -> dict
crosstalk_score(expr, interactions, availability, rule="min", floor=None, adjust="none") -> dict
bootstrap_ci(expr, interactions, availability, n_boot=1000) -> DataFrame
leave_one_out(expr, interactions, availability, by="ligand_subunit_gene") -> DataFrame
decompose(availability, capacity, interactions) -> DataFrame
double_center(mat) -> DataFrame
```

`decompose()` separates the two factors. Availability is **constant across
receivers** (it is a property of the source), so only capacity can explain
differences between receivers — worth stating explicitly, because a two-factor
narrative can otherwise imply both factors vary per receiver when they do not.

### Source specificity by induction, and conjunction scoring

```python
induction_specificity(disease_de, baseline_rank, housekeeping_gene="GAPDH",
                      min_log2fc=1.0, max_fdr=0.05, min_housekeeping_ratio=0.005,
                      max_baseline_rank_frac=0.5, log_scale=True) -> DataFrame
```
Scores source-tissue ligands on **disease induction** instead of healthy-tissue
specificity. Three criteria, each computed within a single dataset so no
cross-platform absolute comparison is implied: A1 replicated induction in every
cohort supplied, A2 abundance as a fraction of a housekeeping gene inside the
source cohorts, A3 the source tissue not excluded at baseline.

Why it exists: a screen scoring cardiac specificity in **healthy** myocardium
returned zero surviving axes. Re-scoring by failure induction recovered seven.
A ligand unremarkable at baseline but strongly induced in disease is a plausible
disease-source ligand, and a healthy-baseline filter discards it.

Read `binding_constraint` before believing `pass_all`. If one criterion removes
nearly everything (in our run A1 removed 43 of 51), the screen is testing that
criterion, not the biology. A3 is the deliberate relaxation — keep it visible so
a reader can see whether a hit depended on it.

A ligand absent from `baseline_rank` **fails A3**; it is never silently skipped.

```python
conjunction_verdict(results, require_all=True) -> dict
```
Scores several sub-analyses as a conjunction rather than a best-of. Where one
hypothesis implies several independent predictions, testing each and reporting
whichever clears significance inflates the false-positive rate by the number of
predictions. Pass `direction_ok` explicitly per conjunct: a result can be
significant in the **wrong** direction, which is evidence against, not for.

This must be fixed **before** the tests run. In our use, one of three conjuncts
passed and the other two pointed the opposite way; reporting the single hit
would have been a false positive.

```python
receptor_pass_rate(receptor_table, receptor_col="receptor", pass_col="recpass",
                   min_types=1) -> dict
```
Base rate at which receptors satisfy the receiver criterion with **no ligand
involved**. Ask this before reading an axis count as evidence for a particular
ligand. We measured 74% — 60 of 81 floor-passing receptors passed with some
tumour type regardless of ligand, which made a seven-axis result uninformative
on its own.

## Preconditions

- `expr` genes x receivers, non-negative, absolute units (not pre-z-scored).
- `interactions` long-form: one row per (interaction, receptor subunit).
- The floor reference gene must be present in `expr` — extract housekeeping
  genes alongside your receptors.
- Availability index must match `ligand_subunit_gene` values, or all weights
  are zero and `crosstalk_score` raises.

## Postconditions

- No scored interaction is below the floor in every receiver.
- `floor_audit` accounts for every candidate interaction (passed or dropped).
- Report `matched_random_null` alongside any ranking. A high null ρ means the
  ranking is not receptor-specific and must not be presented as a biological
  result.

## Worked example

```python
import kernel as ck

floor = ck.housekeeping_floor(expr_by_type, "GAPDH", 0.005)
avail = ck.ligand_availability(de_table, ligands, tier_map=tiers)

# Fails loudly if any receptor is untranscribed everywhere:
try:
    res = ck.crosstalk_score(expr_by_type, interactions, avail, floor=floor)
except AssertionError as e:
    print("below-floor receptors:", e)
    res = ck.crosstalk_score(expr_by_type, interactions, avail,
                             floor=floor, enforce_floor=False)

null = ck.matched_random_null(expr_by_type, interactions, avail,
                              background_genes, n_iter=300, seed=42,
                              floor=floor, enforce_floor=False)
assert null["rho_mean"] < 0.3, "ranking is not receptor-specific"

boot = ck.bootstrap_ci(expr_by_type, interactions, avail, n_boot=400)
dec  = ck.decompose(avail, res["capacity"], interactions)
```

## Interpretation limits

- Receptor expression is **necessary, not sufficient**. The score measures
  compatibility, never ligand exposure, receptor occupancy, or activation.
- Bulk receiver expression mixes malignant, stromal and immune compartments.
- Passing the specificity null means the ranking is receptor-specific. It does
  **not** establish that signalling occurs.
- Compare the primary and `adjust="double_center"` rankings. If they disagree
  (we observed ρ=0.318), the ranking is correction-sensitive and mid-table
  positions should not be over-read.
- **A relaxed criterion needs a base-rate guard.** Loosening a filter to recover
  hits will recover them. Run `receptor_pass_rate()` first and state the number:
  if most receptors pass ligand-free, an axis count says nothing about the
  ligand. Our recovered axes ran on receptors 74% of the universe could match,
  and the recovered set scored *below* the median of expression-matched random
  receptor sets (27 pairs vs null median 30, p=0.688).
- **Induction is not a systemic source.** `induction_specificity()` admits
  ligands on fold-change, which does not establish the source tissue dominates
  the circulating pool. Converting tissue mRNA to a plasma contribution needs
  translation efficiency, secretion rate, clearance and volume of distribution —
  none of which are in expression data. Do not substitute mass-weighted mRNA for
  a measured concentration.
- **Check the cell of origin before calling a tissue a source.** Our recovered
  ligand was 71% fibroblast-derived in single-nucleus data, showed no per-cell
  induction in disease, and was induced in lung and kidney fibrosis at half the
  cardiac effect — a generic wound-healing transcript, not an organ-specific
  signal. Bulk induction that does not reproduce per-cell is a warning, not a
  finding.

## Standalone use (no agent)

`kernel.py` is pure `numpy`/`pandas`/`scipy` — no platform APIs, no network, no
`host.*`. Use it directly:

```bash
pip install numpy pandas scipy statsmodels
python -c "
import kernel as ck, pandas as pd
expr = pd.read_csv('expr_by_receiver.csv', index_col=0)
inter = pd.read_csv('interactions_long.csv')
de = pd.read_csv('source_de.csv')
avail = ck.ligand_availability(de, inter['ligand_subunit_gene'].unique())
print(ck.crosstalk_score(expr, inter, avail, enforce_floor=False)['score'])
"
```

Portable to any agent harness or plain script. The validators hold regardless
of which model (or no model) is driving.
