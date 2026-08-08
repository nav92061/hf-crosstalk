# PTN verdict — the recovered axes do not survive adversarial testing

Replaces the provisional §2.10 wording in [REVISED_SECTIONS.md]({{artifact:3ea2b852-d958-493b-b50a-71feb4daa1d8}}).
Kill criteria fixed in advance: [PRESPECIFICATION_PTN.md]({{artifact:1babf914-e798-402b-aab2-c132187d47a3}}).
Verdicts: rows 42–50 of [validation_summary.csv]({{artifact:b0f219b6-40a7-40a3-9a66-52ba22053274}}).

**No figures accompany this section.** The pre-specification stated that figures would be built only
if PTN survived all three tracks. It did not.

---

## 2.10 (revised) Correcting the specificity filter recovers axes, and the recovered axes do not survive

§2.8 reported zero surviving axes under a ligand filter that scored cardiac distinctiveness in
GTEx left ventricle — non-failing myocardium. That filter is wrong for a hypothesis about the
failing heart, and correcting it changes the count: re-scoring cardiac source by **failure
induction**, with the receptor criteria unchanged, leaves **7 of 146 axes surviving instead of 0**,
all with pleiotrophin (PTN) as ligand.

That is a real methodological finding about the filter. It is not a biological finding about PTN,
and the distinction is the substance of this section. PTN entered only because criterion A3 was
relaxed from "heart ranks first of 52 GTEx tissues" to "heart ranks in the upper half" — PTN sits
at rank 23. A hit that appears when a criterion is loosened is precisely what loosening a criterion
produces, so three kill criteria were fixed before any of the tests below ran, with PTN retained
only if it survived all three. **It failed all three.**

### The receptor test is not discriminating (Track 1)

Before asking whether PTN's receptors are special, we asked how special any receptor is. Of the 81
floor-passing receptors, **60 (74.1%) are enriched versus matched adjacent normal in at least one
tumour type** with no ligand involved. The pre-specified threshold for declaring an axis count
uninformative was 60%. Seven receptors pass in more tumour types than any PTN receptor does
(PLXNA3 in 13, PLXNA1 in 11, TNFRSF12A in 10, ERBB3 in 9, ITGA6, F2R and ITGB4 in 8 each). PTN's
receptors average 3.86 types against an all-receptor average of 2.52 — above the mean, but at the
76th percentile of a distribution most members of which already pass.

The formal test is worse than merely unremarkable. Substituting PTN's seven receptors with
expression-decile-matched random receptor sets of the same size (1,000 permutations, seed 20260729)
gives a null median of 30 enriched pairs against **27 observed** — PTN's set sits at the 31st
percentile of the null, empirical one-sided p = 0.688. Three null configurations agree. **PTN's
receptor set performs slightly worse than random receptor sets matched for expression level.**

One check passed: the seven receptors are not one co-expressed family counted seven times
(mean pairwise Spearman rho −0.091 across the 30 tumour types, maximum 0.692). The axes are
genuinely seven receptors. Seven independent receptors that random draws beat.

### PTN is an injury gene, not a cardiac gene (Track 2)

PTN was tested in five non-cardiac organ-failure and fibrosis cohorts, each processed with the
same moderated-t pipeline and each cohort's normalisation kept separate. COL1A1 rose in all five
(log₂FC 1.17–2.52), confirming every cohort captured a fibrotic signal. PTN was significantly
induced in two: lung IPF (GSE134692, +0.646, FDR 4.6×10⁻³) and diabetic nephropathy (GSE142025,
+0.871, FDR 2.5×10⁻⁵). The median of 0.7589 is **50.5% of the cardiac effect** (1.5023) against a
kill threshold of 50%. The criterion triggers by 1.0% — marginal, pre-specified, and reported as
triggered without adjustment.

The marginality does not rescue PTN, because the direction is inconsistent in a way the threshold
does not capture: in the **largest** cohort tested, lung IPF with 103 cases and 103 controls
(GSE150910), PTN is significantly **down** (−0.556, FDR 2.5×10⁻⁴). The two lung cohorts disagree in
sign. A ligand whose fibrosis response reverses between two cohorts of the same disease is not a
robust injury signal in either direction, let alone a cardiac-specific one.

Cell of origin is more decisive than the induction test. In two independent human heart
single-nucleus atlases, PTN is **fibroblast-derived**: 70.7% of PTN counts in fibroblasts in the
Reichart DCM/ACM atlas and 61.6% in the Kuppe myocardial-infarction atlas, against 6.2% and 16.2%
in cardiomyocytes, with fibroblast:cardiomyocyte mean expression ratios of 35.3× and 13.5×. The failing
heart's PTN is a wound-healing transcript from its stroma, not a cardiomyocyte product.

And the induction does not survive the change of resolution. Per-fibroblast PTN does **not** rise
in disease in either atlas (DCM +0.028, FDR 0.87; ACM −0.311, FDR 0.54; MI −0.430, FDR 0.30), nor
does tissue-level pseudobulk (DCM +0.041, p = 0.72; MI −0.010, p = 0.68). The same MI fibroblast
comparison detects POSTN at +1.77 and COL1A1 at +1.03, so the test had power. **The bulk cardiac
induction of log₂FC 1.50 that admitted PTN in the first place does not reproduce as a per-cell
induction in either single-nucleus dataset.** That discrepancy is unresolved — it may reflect
aetiology or compositional differences between the bulk and single-nucleus cohorts — and it
undermines the premise independently of any kill criterion.

### The measurement the premise requires has never been made (Track 3)

Across 52 PubMed queries over a 991-record pleiotrophin corpus, plus a full-text scan of 478
cardiac plasma-proteomics papers, **609 records were retrieved and verified, and none measures
circulating PTN in a heart-failure cohort against controls**. Four measure circulating PTN in some
human cardiac context, and two of the four point away from the premise:

- **PMID 39523983** (Circ Heart Fail 2024): plasma PTN by SomaScan is **lower** in hypertrophic
  cardiomyopathy than in three comparator cardiomyopathies, in both training and test sets, and
  after adjustment for 19 clinical parameters.
- **PMID 21685905** (Nat Biotechnol 2011): coronary-sinus sampling during planned myocardial
  infarction detected PTN release — and the authors **explicitly eliminated it** as not specific to
  myocardial injury, because catheterisation alone moved it.

The two supporting records are outside heart failure: anthracycline cardiomyopathy in childhood
cancer survivors (PMID 35861824; only 8 of 28 cases had heart failure, in a population where serum
PTN is independently raised by cancer), and coronary collateral grade (PMID 28885394).

PTN's secreted role in cancer is by contrast well established. The half of the axis that is
supported is the tumour half; the cardiac-plasma half is unmeasured, and where adjacent
measurements exist they run the wrong way. The one study designed to detect cardiac PTN release
into blood found a signal and rejected it as an artefact of the procedure.

### Verdict

**The seven PTN axes are withdrawn as a positive finding.** All three pre-specified tracks failed:
the receptor test is passed by three-quarters of candidate receptors and PTN's set underperforms
matched-random draws; PTN is a fibroblast injury transcript induced in lung and kidney fibrosis at
half the cardiac effect, with no per-cell cardiac induction; and the plasma measurement the
endocrine premise depends on has never been made, with the nearest published evidence contradicting
it.

The paper's conclusion is therefore unchanged from §2.8, but it is now a stronger claim than it was.
The original negative could have been an artifact of scoring cardiac specificity in healthy tissue.
That objection has been tested directly: correcting the filter does recover axes, and the recovered
axes do not survive. A negative that has been attacked on its most obvious weakness and held is
worth more than one that has not.

Tables: [receptor_base_rate_audit.csv]({{artifact:be8dbf3d-c405-436b-b242-63307d4c42ae}}), [ptn_receptor_null.csv]({{artifact:a5cb65bd-bc49-4bc0-9d71-6f9b8f79a1f0}}),
[ptn_receptor_correlation.csv]({{artifact:a0730b02-43e8-4f8c-b49a-266dc836c6bc}}), [ptn_cross_tissue_induction.csv]({{artifact:9e390931-ffea-4a7c-8e6f-1f2a2cd91c9a}}),
[ptn_cell_of_origin.csv]({{artifact:097cd740-0258-42f7-aa58-cab111d76ad4}}), [ptn_literature.csv]({{artifact:7aaa38a5-5d5e-4f6c-9509-530d3a4d1a46}}),
[ptn_literature_summary.md]({{artifact:f8cbeafd-be14-4261-a6d9-4e4c32d0fea0}}).

---

## Consequential edits

**Abstract.** Replace the provisional wording proposed in the previous revision with:

> *Zero axes satisfied both specificity criteria when cardiac distinctiveness was scored in healthy
> myocardium. Because that filter cannot see a ligand induced only in failure, we re-scored cardiac
> source by failure induction: seven axes survive, all on pleiotrophin. Pre-specified adversarial
> testing withdraws them — 74% of candidate receptors pass the receptor test ligand-free and PTN's
> set underperforms expression-matched random draws (p = 0.69); PTN is a fibroblast transcript
> (71% of counts) induced in lung and kidney fibrosis at half the cardiac effect, with no per-cell
> cardiac induction; and circulating PTN has never been measured in heart failure, while plasma PTN
> is lower in hypertrophic cardiomyopathy and was rejected as procedurally non-specific in
> coronary-sinus sampling. No axis survives at transcript level.*

**§3, item 10.** *Zero axes pass when cardiac source is defined in healthy myocardium. Seven pass
when it is defined by failure induction, and all seven are withdrawn by pre-specified adversarial
testing (§2.10). No axis survives.*

**§3, new not-supported item.** *That pleiotrophin is a cardiac-specific ligand. It is a
fibroblast wound-healing transcript, induced in non-cardiac fibrosis at half the cardiac effect,
with no per-cell induction in failing myocardium.*

**§4, new limitation.** *The seven-axis result was recovered by relaxing one criterion and withdrawn
by testing the base rate of the criterion it was recovered against. Any screen that loosens a filter
to recover hits should report the base rate at which the unmodified half of the test passes; ours
was 74%.*

**§4, new limitation.** *Bulk cardiac PTN induction (log₂FC 1.50, replicated across two cohorts) does
not reproduce as per-cell induction in two single-nucleus atlases. Bulk differential expression in a
tissue whose composition changes with disease can reflect cell-proportion shifts rather than
per-cell regulation; this screen's ligand-availability term is computed from bulk and inherits that
ambiguity for every ligand, not only PTN.*

**§5, Methods.** The `crosstalk-score` skill now carries `induction_specificity()`,
`conjunction_verdict()` and `receptor_pass_rate()`, implementing the corrected criteria, the
conjunction rule and the base-rate guard as assertions with regression tests reproducing the numbers
in §2.10–§2.11.
