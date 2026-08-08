# Numeric audit of PAPER.md

Audited version: `cd9cc5a5-cf76-41fb-b388-9b06786b505b` (26,066 bytes, 451 lines).
Every numeric and structural claim was checked against the saved artifacts, recomputing from
source tables wherever the artifact stored a summary rather than the underlying data.

**61 claims checked: 49 verified exactly, 12 findings.** Three are outright numerical errors.
No finding changes any conclusion of the study.

| Severity | Count | Meaning |
|---|---|---|
| ERROR | 3 | Number in the paper does not match the data |
| INCONSISTENCY | 2 | Number is real but attached to the wrong denominator or population |
| HAZARD | 1 | Index points a reader at a table that cannot support the claim |
| STALE | 3 | True when written, superseded by later phases |
| AMBIGUOUS | 1 | Correct but requires the reader to reconstruct it |
| COSMETIC | 1 | Presentation only |
| INCOMPLETE | 1 | Checked and found sound; recorded for completeness |

---

## The three errors

### E1. 710 TCGA adjacent-normal samples
**Location:** §2.7 line 263; abstract
**Paper says:** 710 TCGA adjacent-normal samples
**Data says:** 680 — sum of per-type counts in receptor_tumor_vs_normal.csv (BRCA 114, KIRC 72, THCA 59, LUAD 58, PRAD 52, LUSC 51, LIHC 50, HNSC 43, COAD 40, KIRP 32, STAD 32, KICH 25, UCEC 22, BLCA 19, ESCA 11)
**Why it happened:** 710 is n_mapped_to_tumor_type from specificity_summary.json — samples mapped to a tumour type BEFORE the >=10-normals-per-type filter. 30 of them fall in types that were then dropped. The figure actually used is 680.
**Fix:** Change to 680, or say '710 mapped, 680 in the 15 analysed types'

### E2. NPR1 fails the expression floor in 26 of 30 types
**Location:** §2.7 line 272; abstract; §2.9 summary
**Paper says:** NPR1 fails the expression floor in 26 of 30 types
**Data says:** 24 of 30 — NPR1 passes the floor in KIRC, KIRP, MESO, OV, PCPG, SARC (passes_tcga_floor_in_type in receptor_tumor_vs_normal.csv)
**Why it happened:** 26 traces to specificity_summary.json prose written by the sub-agent; the column it summarises says 24. Never independently recomputed in the original audit.
**Fix:** Change to 24 of 30. Conclusion unaffected — NPR1 still fails in the large majority.

### E3. RPL3 -2.53, PLK1 -2.69, POLR2B -2.42, EEF2 -2.35
**Location:** §2.6 line 224
**Paper says:** RPL3 -2.53, PLK1 -2.69, POLR2B -2.42, EEF2 -2.35
**Data says:** From depmap_lineage_dependency.csv (mean over 25 lineages): RPL3 -2.548, PLK1 -2.674, POLR2B -2.395, EEF2 -2.366
**Why it happened:** Paper quotes depmap_summary.json's positive_control_gene_effects (-2.533, -2.693, -2.418, -2.350), which are pan-CELL-LINE means. The per-lineage table gives different values. Both defensible, but the paper cites neither source explicitly and rounds inconsistently.
**Fix:** State which aggregation is quoted; -2.53/-2.69 are pan-cell-line, -2.55/-2.67 are pan-lineage.

## Numbers that are real but misattached

### I4. §2.7 line 275 vs §2.8 tables
**Claim:** zero of 160 floor-passing interactions survive
**Verified:** axis_survival.csv has 146 rows, all survives_both=False. 160 is the count of floor-passing INTERACTIONS (complex-level, e.g. BMP4_BMPR1A_ACVR2A); 146 is ligand->receptor-SUBUNIT pairs (BMP4->ACVR2A).
**Problem:** Two different units of analysis given one number. Both counts are individually correct; the sentence attaches the survival test to the wrong denominator.
**Fix:** 'zero of 146 ligand-receptor-subunit axes derived from the 160 floor-passing interactions'

### I5. §2.8 line 293
**Claim:** FZD7 vs housekeeping mean rho=0.131, p=0.447 — reported in the sentence describing the 4v4 panel
**Verified:** 0.131/0.447 is receptor_vs_housekeeping_pool_rho over all DepMap Uterus lines. Recomputed on the 16 panel rows: rho=0.437, p=0.091; on the 8 primary-panel rows: rho=0.214, p=0.610.
**Problem:** A pool-level statistic placed in a sentence about the panel. The claim it supports (receptor level is not a library-size artefact) holds at every level, but the number quoted is not the panel's.
**Fix:** Say 'across the 36 screened Uterus lines' rather than implying the panel.

### H10. TABLE_INDEX.csv rows 17, 18, 20
**Claim:** version_id pointers for receptor_tumor_vs_normal.csv, ligand_tissue_specificity.csv, tumor_normal_tissue_map.csv
**Verified:** Each filename exists as TWO separate artifacts. TABLE_INDEX points at b9491a67 (13 cols) but the paper's numbers come from 9c318e4e (25 cols); it points at c4b0258c (51 rows) but the '2 of 42' claim needs a82f4063 (42 rows).
**Problem:** A reader following TABLE_INDEX gets a table lacking the columns behind the claims — n_tcga_adjacent_normals, log2FC_bg_corrected_within_tcga, fdr_within_tcga, passes_tcga_floor_in_type, consensus_direction are all absent from the indexed version.
**Fix:** Repoint rows 17/18/20 at the versions the results use, or list both with roles.

### A12. §2.7 line 252-256
**Claim:** only 2 of 42 ... The remaining 39 peak elsewhere
**Verified:** 2 + GDF6 + 39 = 42. Arithmetic is exact.
**Problem:** 'Remaining 39' silently excludes GDF6, which was described in the preceding clause. Correct but requires the reader to do the subtraction.
**Fix:** 'The other 39' -> 'The remaining 39 (excluding GDF6)'.

## Superseded by later phases

### S7. Abstract; §3 items 3,5; §2.7 heading
**Claim:** Zero axes satisfied both specificity criteria / specificity checks refute the axis-level hypothesis
**Current state:** Superseded by the corrected-specificity phase: 7 axes survive when cardiac source is scored by failure induction (axis_survival_corrected.csv), all PTN. Those 7 were then withdrawn by adversarial testing (validation rows 42-50).
**Note:** Not wrong as written — it is true under the healthy-heart filter — but the paper does not yet contain §2.10/§2.11 or the PTN verdict, so a reader sees a conclusion the project has since tested twice.
**Fix:** Merge REVISED_SECTIONS.md and PTN_VERDICT.md into PAPER.md; apply the abstract wording in PTN_VERDICT.md.

### S8. §6 opening; §7
**Claim:** Six errors were found and corrected ... rows 26-31 of validation_summary.csv / 31 controls
**Current state:** validation_summary.csv now has 50 rows (19 PASS, 13 FAIL, 7 CORRECTED, 4 NULL RESULT, 2 UNINFORMATIVE, 2 LARGELY FAIL, 1 each ENFORCED/FLAG/LARGELY NULL). Rows 26-31 are still the six audit corrections.
**Note:** Row references remain valid; the totals in §7 do not.
**Fix:** Update §7 counts to 50 rows, and add the errors found in THIS audit as new rows.

### S9. §7; TABLE_INDEX.csv
**Claim:** Twenty-two numbered tables (TABLE_INDEX.csv)
**Current state:** TABLE_INDEX.csv has 22 rows and does not list any of the 11 tables produced after it: ligand_specificity_failure_based, axis_survival_corrected, npr1_survival/proliferation/copynumber, receptor_base_rate_audit, ptn_receptor_null, ptn_receptor_correlation, ptn_cross_tissue_induction, ptn_cell_of_origin, ptn_literature.
**Note:** Index frozen at the first phase.
**Fix:** Extend TABLE_INDEX to 33 rows; FIGURE_INDEX to 7 (fig6, fig7 missing).

## Presentation

### C11. lines 169 and 201
Figure 3 is embedded at line 169, Figure 2 at line 201. Each embed correctly matches its own caption and FIGURE_INDEX entry; only the document order is inverted. Figures appear 1, 3, 2, 4, 5. There are no in-text 'Figure N' cross-references, so nothing resolves incorrectly — but numbering should follow appearance.
**Fix:** Swap the two blocks, or renumber.

---

## What verified exactly

All 49 are listed in `PAPER_AUDIT_VERIFIED.csv`. The load-bearing ones:

- **Replication:** 5,036 replicated genes, 18,061 shared, 5,717 significant in both, concordance
  0.8809, Spearman 0.5408 all shared and 0.7593 among sig-in-both, hypergeometric 2.15×10⁻³⁵,
  fold enrichment 1.073, permutation control 0 genes with rho −0.1344. Concordance rises 0.846
  (|log₂FC|<0.5) to 0.9525 (>1.0), matching the paper's "85% to 95%". Canonical markers exact.
- **Floor:** GAPDH median 69,134.5, floor 345.67, 28 interactions dropped, 160 retained, 81
  receptors. All four PENK→opioid axes below floor at 4.48×10⁻⁵, 4.97×10⁻⁴, 1.37×10⁻³ and
  2.97×10⁻³ of GAPDH.
- **Ranking:** UCS 0.6785 (CI 0.3997–0.9142, rank SD 0.549), exactly 4 of 30 beating the null
  (UCS 0.0033, BRCA 0.0100, KIRC 0.0133, PAAD 0.0332), OV second by score at p=0.0598, ACC first
  under double-centring with UCS sixth.
- **Nulls and stability:** 0.363 no-floor, 0.107 floor, 0.323 naive pool, 0.034 double-centred;
  21.5% and 1.7% of draws above 0.5; p95 0.4263; all eight stability controls at their stated values.
- **Receptors:** two methods agree at ρ=0.545, p=8.91×10⁻⁹⁴, n=1,200 — recomputed from raw columns.
  NPR1 depleted in 15 of 15 unit-matched types, median −1.9135, maximum FDR 2.04×10⁻³.
- **DepMap:** 81 receptors, 22 mapped lineages (22 of 31 flagged `included_in_correlation`), zero
  at FDR<0.05, minimum FDR 0.784. FZD7 panel separation 48.72×, gene effects +0.0583 to +0.3291,
  housekeeping arms Mann-Whitney p=1.000.
- **Figure integrity:** all five embeds resolve to the artifacts their captions describe and match
  FIGURE_INDEX.

## Method

Where an artifact stored only a summary, the underlying data was recomputed rather than trusted:
concordance and both correlations from the two DE tables; the receptor cross-method ρ from
`delta_percentile_bg_corrected` against `log2FC_bg_corrected_within_tcga`; the opioid ratios from
`tcga_expr_by_tumor_type.csv` against the GAPDH median; the panel fold separation from
`2**log2tpm1 - 1`. Errors 1, 2 and 3 were all found this way — each is a case where the paper
quoted a sub-agent's prose summary instead of the column the summary describes.
