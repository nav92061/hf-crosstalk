# Supplementary material

*Companion to* "A reusable agent-skill pipeline for organ-crosstalk screens, and what
it found looking for heart-failure signals in 30 tumor types"

## S1. Exploratory: NPR1 is not a lost tumor suppressor on the evidence tested

This section is reported as supplementary because it tests a hypothesis generated
by the refutation in Section 3.8 of the main text, not one predicted by the
original design. The original hypothesis was that a cardiac ligand engages a
tumor receptor; the hypothesis here is close to its opposite, that tumors benefit
from shedding that receptor. Presenting it in the main results would misrepresent
the order in which it was formed. The conjunction rule below was fixed before any
of the three predictions was tested, and the analysis is retained in full because
it failed.

### 3.10 Exploratory: NPR1 is not a lost tumor suppressor on the evidence tested

**Status.** This section tests a hypothesis generated *by* the refutation in Section 3.8, not one predicted
by the original design. The original hypothesis was that a cardiac ligand engages a tumor
receptor; the hypothesis here is close to its opposite, that tumors benefit from *shedding* that
receptor. It is reported separately for that reason. Presenting it as a prediction of the screen
would misrepresent the order in which it was formed.

Section 3.7 established that NPR1, the sole receptor for the only two cardiac-enriched ligands, is
depleted in tumor relative to matched adjacent normal in all 15 testable types. Natriuretic
peptide signaling through NPR1 is antiproliferative, which admits a reading in which NPR1 loss
is advantageous to the tumor: the failing heart would not be engaging tumor receptors but
altering the selective landscape in favor of tumors that have already lost them.

Three predictions follow, and a conjunction rule requiring all three to hold in the specified
direction was fixed before any of them was tested. The rule exists because three independent
opportunities to find a hit is how a negative result becomes a false positive.

**B1, survival: FAIL.** In a Cox model over 6,764 patients and 1,741 events, NPR1 z-scored within
tumor type and stratified by type, the hazard ratio is 1.0285 (95% CI 0.9809–1.0784, p = 0.2446).
The pre-specified direction was HR < 1; the point estimate is on the wrong side of 1 and the result
is not significant. Per-type, 9 of 15 have HR < 1, but only KIRP reaches FDR < 0.05 in the
predicted direction while STAD and UCEC are significant in the *opposite* direction.

**B2, proliferation coupling: PASS.** Within-type Spearman correlation between NPR1 and the
proliferation signature is negative in 13 of 15 types at FDR < 0.05 and positive in none, with a
median rho of −0.2032 (range −0.578 in STAD to +0.007 in UCEC; 6,596 samples matched). This is the
predicted direction and the pre-specified threshold of at least 8 of 15 is met comfortably.

**B3, genomic loss: FAIL, and informatively so.** Against a null of 1,000 loci matched on
chromosome arm 1q and on overlapping-segment length, NPR1 falls below its ploidy in 3.10% of
6,648 samples. The pass threshold was the null's 95th percentile, 4.93%. NPR1 does not merely miss
that threshold. At the 0.8th percentile of the null, below the null median of 4.09%, NPR1 is
deleted less often than comparable loci. A sensitivity analysis using ABSOLUTE table ploidy
instead of segment-derived ploidy fails likewise (4.45% against 6.83%).

**Conjunction verdict: FAIL, 1 of 3.** Two conjuncts point away from the hypothesis and the one
that passes is the weakest evidentially: an inverse correlation between a receptor and a
proliferation score is compatible with the hypothesis but also with proliferating tumor cells
simply expressing less of a differentiation-associated receptor, which is not a selective claim.
B3 in particular argues against selection: if losing NPR1 conferred advantage, the locus should
be deleted more often than matched loci, and it is deleted less. Promoter methylation was
pre-specified as a conditional secondary read available only if all three conjuncts passed; that
precondition never held, so it was not assessed, and its absence is not treated as a reason to
soften B3.

The finding of Section 3.7 therefore stands as reported: NPR1 depletion is real and it defeats the
natriuretic peptide axis, but on the evidence tested it is not a selected loss of a tumor
suppressor. It remains most parsimoniously explained as tumor tissue expressing less of a
receptor characteriztic of the differentiated normal tissue (Fig. 7).

![Figure 7](art_cf1f415e-3cf5-4885-8bed-b117db44868e)

**Figure 7. Exploratory NPR1 analysis; the conjunction fails.** (a) Per-type Cox hazard ratios per
standard deviation of NPR1, dark points at FDR < 0.05; the pooled estimate is HR 1.03, p = 0.24
against a pre-specified HR < 1. (b) Within-type Spearman correlation with the proliferation
signature; green points are significant at FDR < 0.05, negative in 13 of 15 types. (c) NPR1 locus
deletion frequency against the arm-1q matched null; NPR1 is below the null median, not above the
95th percentile. All three panels report an exploratory hypothesis generated by the refutation.



## S2. Implementation validation against canonical limma

![Figure S2](art_5f24b105-b818-4c4e-8f44-f2a0d615a82b)

**Figure S2. Implementation validation against canonical limma.** (a, b) Per-gene log₂ fold-change agreement between the implementation used here and limma 3.66.0 with `eBayes(trend=TRUE)`, for the discovery and validation cohorts; the dashed line is y = x. Pearson and Spearman correlations are 1.000000 in both cohorts, with maximum absolute differences of 4×10⁻¹⁴ and 8×10⁻¹⁴. (c, d) Overlap of the FDR<0.05 gene sets, Jaccard 0.99913 and 0.99978. The validation criterion (Pearson ≥ 0.99 and set overlap ≥ 0.95 in both cohorts) was fixed before the comparison ran.

`limma_validation.csv` reports all 83 comparison rows, including the pre-fixed
criterion, both cohorts, three limma variants and the one blocked comparison.
`limma_de_GSE116250.csv` and `limma_de_GSE141910.csv` carry the canonical limma
per-gene results for independent checking.

The implementation was validated against canonical limma 3.66.0 with the criterion
fixed beforehand: agreement required a log₂ fold-change Pearson correlation of at
least 0.99 and at least 95% overlap of the FDR<0.05 gene sets, in both cohorts.
Both cohorts were re-analyzed from source with `lmFit` and `eBayes(trend=TRUE)`.
Agreement is exact to floating-point precision: log₂ fold-change Pearson and
Spearman correlations are 1.000000 in both cohorts, with maximum absolute
differences of 4×10⁻¹⁴ and 8×10⁻¹⁴; moderated *t* correlations are 0.9999991 and
0.9999999. The FDR<0.05 sets agree at Jaccard 0.9991 and 0.9998, sharing 9,226 of
9,230 and 13,612 of 13,614 genes. The empirical-Bayes hyperparameters were
estimated independently and match: prior degrees of freedom 3.041 against limma's
3.037 in the discovery cohort, and 4.818 against 4.825 in the validation cohort.
Applying the replication rule to the limma output recovers 5,036 genes sharing
5,033 with the signature reported here. A `trend=FALSE` variant, whose global prior
differs materially (prior degrees of freedom 2.80), also agrees at Pearson
1.000000. One comparison could not be made: GSE141910 has no count matrix deposited
in GEO, only per-sample values that are already log₂-scale and library-normalized,
so edgeR TMM normalization with voom precision weights is not applicable and
`eBayes(trend=TRUE)` is the correct canonical analogue rather than a convenient one.

![Figure 1](art_5f24b105-b818-4c4e-8f44-f2a0d615a82b)

**Figure S1. Implementation validation against canonical limma.** (a, b) Per-gene log₂ fold-change agreement between the implementation used here and limma 3.66.0 with `eBayes(trend=TRUE)`, for the discovery and validation cohorts; the dashed line is y = x. Pearson and Spearman correlations are 1.000000 in both cohorts, with maximum absolute differences of 4×10⁻¹⁴ and 8×10⁻¹⁴. (c, d) Overlap of the FDR<0.05 gene sets, Jaccard 0.99913 and 0.99978. The validation criterion (Pearson ≥ 0.99 and set overlap ≥ 0.95 in both cohorts) was fixed before the comparison ran.

## S3. Reviewer-requested sensitivity analyses

The six sweeps summarized in the main text are reported in full in the accompanying
tables: `floor_sensitivity.csv` (expression floor across six fractions and two
alternative definitions), `specificity_definition_sensitivity.csv` (six cardiac
specificity definitions), `permutation_resolution.csv` (300 to 20,000 draws per
tumor type with Monte Carlo standard errors), `null_pool_width.csv` (2 to 40
matching bins), `complex_independence.csv` (per-complex versus per-subunit scoring
for all 171 receptor complexes) and `etiology_stratified.csv` (per-etiology
differential expression for the 42 secreted ligands).

## S4. Analysis timeline

`ANALYSIS_TIMELINE.md` gives artifact-store version timestamps for every
pre-specification and every result it governs, the full revision history of each
pre-specification, and an explicit statement that these are provenance records
rather than independent preregistration.
