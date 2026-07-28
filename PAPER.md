# A ligand–receptor compatibility screen finds no transcriptome-level evidence that the failing heart meets a tumour-specific receptor

**Status:** working manuscript. Every number below is traceable to a saved artifact
produced by the pipeline described in Methods. Controls that failed are reported
alongside those that passed.

---

## Abstract

Heart failure (HF) has been proposed to promote cancer through systemically
released factors, a hypothesis supported by mouse experiments but unresolved in
humans. Because published work has examined one tumour model at a time, it is
unknown whether some cancer types are intrinsically more able to receive such
signals than others. We formalised that question as a two-factor model —
**ligand availability** from the failing heart × **receiver capacity** of each
tumour type — and tested it across 30 solid tumour types.

We derived a human HF signature from two independent left-ventricular RNA-seq
cohorts (GSE116250 discovery, n=64; GSE141910 validation, n=360), retaining
5,036 genes significant at FDR<0.05 with concordant direction in both
(directional concordance 0.881; effect-size Spearman ρ=0.541 across all 18,061
shared genes). Filtering to Human Protein Atlas secreted proteins and mapping to
CellChatDB yielded 188 candidate secreted-signalling interactions from 51
HF-upregulated ligands. Receptor expression was measured in 9,538 TCGA primary
solid tumours.

The methodological result is that a mandatory absolute-expression floor is what
makes such a score interpretable. Without it, expression-matched random receptor
sets reproduced our ranking at ρ=0.363 (21.5% of draws exceeding ρ=0.5); with a
floor at 0.005×GAPDH, this fell to ρ=0.107 (1.7%). The confound usually
attributed to a "global expression axis" is therefore substantially *caused by*
scoring receptors that are not transcribed. The floor removed 28 of 188
interactions, including all four proenkephalin→opioid-receptor axes.

The biological result is negative and specific. Uterine carcinosarcoma ranked
highest (score 0.679, bootstrap 95% CI 0.400–0.914), but the ranking reordered
substantially under double-centring (ρ=0.318), and receiver capacity did not
track tumour microenvironment state (1 of 29 features at FDR<0.05). Of 42
HF-upregulated secreted ligands, only NPPA and NPPB are cardiac-enriched and
abundant in GTEx; NPR1, their sole floor-passing receptor, is **depleted** in
tumour relative to adjacent normal in all 15 unit-matched tumour types (median
log₂FC −1.91, all FDR<0.05) and fails the expression floor in 26 of 30 types. No
receptor showed compatibility-correlated CRISPR dependency (0 of 81, FDR<0.05).
**Zero axes satisfied both specificity criteria.**

We conclude that transcriptome-level data do not support a model in which a
heart-specific ligand meets a tumour-specific receptor. The hypothesis survives
only in the weaker form that HF may raise systemic levels of broadly secreted
factors. We release the six analysis skills as reusable, platform-independent
modules, with the floor and the specificity null enforced as code-level
assertions rather than documentation.

---

## 1. Introduction

Reverse cardio-oncology asks whether cardiovascular disease promotes cancer
rather than the reverse. The foundational experiment is mechanistic: heterotopic
transplantation of a failing heart into tumour-bearing mice accelerated tumour
growth, implying a blood-borne signal rather than a haemodynamic one
(Meijers 2018, PMID: 29459363). Myocardial-infarction models likewise accelerate
breast cancer growth through altered myeloid output (Koelwyn 2020, PMID:
32661390; Avraham 2020, PMID: 32475164).

Human evidence is far weaker and openly conflicting. Several cohorts report
excess cancer incidence after HF (Hasin 2013, PMID: 23810869; Banke 2016, PMID:
26751260; Roderburg 2021, PMID: 34180146), but the largest study with
adjudicated cancer outcomes and 19.9-year median follow-up found no association
(adjusted HR 1.05, 95% CI 0.86–1.29; Selvaraj 2018, PMID: 29622155), and
introducing a two-year diagnostic lag attenuated a national cohort's hazard
ratio from 1.64 to 1.09 (Kwak 2020, PMID: 32863081). Conventional
cardiovascular risk factors predict incident cancer as well as incident HF
(PMID: 41014601), and clonal haematopoiesis and low-grade inflammation are
shared substrates for both diseases (Jaiswal 2017, PMID: 28636844; Fuster 2017,
PMID: 28104796). The direction of effect is not even settled: tumour-bearing
mice subjected to pressure overload developed *less* cardiac hypertrophy and
fibrosis (Awwad 2023, PMID: 37508517). Our literature base comprises 107
verified PubMed records, 21 of which we classified as contradicting a causal
HF→cancer link.

This work does not attempt to settle causality; no available dataset links HF
status to cancer incidence with the molecular resolution required. We instead
address a narrower, tractable question that the field has not systematically
asked: **if the failing heart releases signals, which tumour types are
transcriptionally equipped to receive them, and does any candidate axis have
both a distinctively cardiac source and a tumour-specific receptor?**

Framing it as compatibility rather than causation makes the analysis falsifiable
in a useful way, because it can fail in three distinct places — the ligand may
not be cardiac-specific, the receptor may not be expressed, or the receptor may
not be tumour-enriched — and we test all three.

---

## 2. Results

### 2.1 A replicated human heart-failure signature

Differential expression in GSE116250 (50 failing hearts — 37 dilated, 13
ischaemic cardiomyopathy — vs 14 non-failing) recovered the expected cardiac
stress programme: NPPA log₂FC +3.91, NPPB +3.97, with MYH6 (−0.61) and ATP2A2
(−0.44) down, as expected in failing myocardium. GSE116250 distributes RPKM
rather than counts, so a limma-trend-equivalent moderated *t* on log₂(RPKM+1)
was used; 65 genes carrying a sentinel value of 999999999 were detected and
removed before analysis.

Validation in the larger MAGNet cohort (GSE141910; 194 cardiomyopathy vs 166
non-failing donors) replicated the signature: of 18,061 shared genes, 5,717 were
significant in both cohorts and 5,036 shared direction — a directional
concordance of 0.881, rising monotonically with discovery effect size (85% at
|log₂FC|<0.5 to 95% above 1.0). Effect sizes correlated at Spearman ρ=0.541
across all shared genes and ρ=0.759 among those significant in both.

Two cautions matter. First, the hypergeometric overlap test returns
p=2.1×10⁻³⁵ but only 1.07× enrichment: with 8,052 and 11,949 significant genes
respectively, the test is near-saturated and is *not* evidence of replication.
Directional concordance and effect-size correlation are the informative metrics.
Second, we ran a label-permutation control through the identical pipeline: it
returned **0 significant genes** and ρ=−0.134 against the real validation
result, confirming the pipeline does not manufacture concordance.

![Figure 1](art_32f64f02-cc63-4a6a-865e-9e0c0eb0c07b)

**Figure 1.** Human HF signature. (a) Discovery volcano, GSE116250. (b)
Cross-cohort effect-size concordance among genes significant in both; red points
disagree in direction. (c) Concordance increases with discovery effect size.

### 2.2 From signature to candidate crosstalk network

Of 3,084 replicated HF-upregulated genes, 295 are annotated as secreted by the
Human Protein Atlas and 72 are CellChatDB ligands. Restricting to
`Secreted Signaling` — the only annotation class mechanistically plausible at a
distance — gave **188 interactions from 51 ligands across 38 pathways**, with
receptor complexes expanded to 103 constituent subunits. Ligand availability was
computed as the mean replicated effect size scaled by secretion tier; the highest
were NPPA (3.05), NPPB (3.03), PENK (2.85) and WNT9A (2.32).

### 2.3 The expression floor is the load-bearing methodological choice

TCGA PanCanAtlas values are RSEM normalized counts, not TPM: GAPDH's median
across the 30 per-tumour-type medians is 69,134, so any floor expressed in
transcript-per-million terms is meaningless here. We therefore defined the floor
relative to a housekeeping gene — 0.005×GAPDH = 345.7 in these units.

Twenty-eight of 188 interactions fall below it in *every* one of the 30 tumour
types. These include **all four** proenkephalin→opioid-receptor axes (OPRD1,
OPRK1, OPRL1, OPRM1): OPRD1 reaches at most 4×10⁻⁵ of GAPDH, OPRM1 5×10⁻⁴ and
OPRK1 1.4×10⁻³. PENK is the third most available HF ligand, so an unfiltered
analysis would nominate the opioid axis prominently; the entire pathway is
non-viable on receptor-expression grounds.

The floor is not hygiene, it is what makes the score mean anything. Under an
expression-matched random receptor null (300 permutations, seed 20260727):

| Null configuration | Mean ρ vs observed | Draws with ρ>0.5 |
|---|---|---|
| No floor applied | 0.363 | 21.5% |
| **Floor applied (primary)** | **0.107** | **1.7%** |
| Floor applied, double-centred | 0.034 | 1.0% |
| Unmatched (naive) random pool | 0.323 | 17.0% |

Two conclusions follow. The "global expression axis" confound that afflicts
co-expression rankings is substantially *caused by* scoring untranscribed
receptors — remove them and the ranking becomes receptor-specific. And the null
pool must be expression-matched: a naive pool gives ρ=0.323 and would wrongly
condemn a specific score.

![Figure 3](art_62ff7487-a8e3-476a-83df-92f5b5f22f80)

**Figure 3.** (a) Floor audit; 28 interactions (red) never clear the floor in any
tumour type. (b) Two-factor decomposition. Availability is a property of the
source and is constant across receivers, so only capacity can explain
differences between tumour types.

### 2.4 Pan-cancer ranking, and how much of it to believe

Scoring the 160 floor-passing interactions ranked uterine carcinosarcoma first
(0.679, bootstrap 95% CI 0.400–0.914, rank SD 0.55 over 400 resamples), followed
by ovarian serous cystadenocarcinoma (0.432), clear-cell renal carcinoma (0.377),
thyroid (0.285), prostate (0.272) and pancreatic adenocarcinoma (0.244). Only
**4 of 30** types beat the matched-random null at p<0.05 — UCS (p=0.003), KIRC
(p=0.013), BRCA (p=0.010) and PAAD (p=0.033) — and notably ovarian carcinoma,
ranked second by score, is **not** among them (p=0.060). Rank order and
receptor-specificity are therefore not the same thing, and a high score alone is
insufficient grounds to nominate a tumour type.

Stability controls pass: leave-one-ligand-out preserves the ranking at minimum
ρ=0.907 (worst case WNT9A) and leave-one-pathway-out at ρ=0.842 (WNT, the
largest single contributor at 80 of the floor-passing interactions). Switching
the receptor-complex rule from minimum-of-subunits to geometric mean gives
ρ=0.908, and removing availability weighting altogether gives ρ=0.922 — so no
single ligand, pathway, or parameter choice drives the result.

One control does not pass. Double-centring, which removes additive row and
column effects, reorders the ranking substantially: ρ=0.318 against the primary,
promoting adrenocortical carcinoma to first and demoting UCS to sixth. Ovarian
carcinoma is the only type in the top two under both variants. **The ranking is
correction-sensitive, and mid-table positions should not be interpreted.**

![Figure 2](art_b2c7e5cf-eb22-48ce-a162-0ceeb6fc09e1)

**Figure 2.** (a) Ranking with bootstrap CIs; the 4 blue types beat the
matched-random null at p<0.05. (b) The floor is what makes the score
receptor-specific.

### 2.5 The "exploit" half of the hypothesis is unsupported

Our hypothesis held that susceptibility depends on the ability to *receive and
exploit* cardiac signals. Availability is constant across tumour types by
construction, so only receiver capacity can explain differences. Testing capacity
against 29 Thorsson immune-landscape features (substituting for TIMER2, which was
unreachable), **only one survives FDR correction** — B-cell receptor score
(ρ=0.571, FDR=0.028). Leukocyte fraction, stromal fraction, TGF-β score,
proliferation and wound-healing signatures are all null. Receptor compatibility
and microenvironment state are largely independent axes; we found no evidence for
the exploitation component as formulated.

### 2.6 Receptor expression does not correspond to lineage dependency

Using DepMap 24Q4 Public (1,178 cell lines, resolved via figshare article
27993248 because depmap.org serves a bot interstitial), we correlated per-lineage
CRISPR gene effect against compatibility for 81 floor-passing receptors across 22
mapped lineages. Positive controls confirm correct parsing: RPL3 −2.53, PLK1
−2.69, POLR2B −2.42, EEF2 −2.35, with non-essential MYT1 −0.02 and HBB +0.16.

**No receptor is significant after FDR correction (0 of 81; 3 nominal.)** More
pointedly, the receptors on the leading axes are not dependencies anywhere: NPR1
has pan-lineage mean gene effect −0.163, and every WNT9A receptor subunit lies
between LRP5 (−0.163) and FZD6 (+0.146), with no lineage cell crossing −0.5. The
eight receptors that *are* dependencies (NCL, ITGAV, EGFR, ITGB5, ITGB1, ERBB2,
FGFR1, IGF1R) are pan-essential or canonical adhesion/RTK genes whose
essentiality is lineage-driven and unrelated to HF ligand availability.

This is a meaningful distinction rather than a disappointment: receptor
expression licenses a cell to *receive* a signal, but its loss need not impair
proliferation in vitro. The nominated axes are candidate signalling conduits, not
drug targets, and should be tested by ligand stimulation with a pathway reporter
rather than by viability knockout.

![Figure 4](art_19780104-a0e1-43b7-8cf4-58699de2941e)

**Figure 4.** (a) HF-axis receptors cluster at zero gene effect while essential
controls sit near −2.4. (b) No receptor survives FDR. (c) Proposed uterine panel.

### 2.7 Specificity checks refute the axis-level hypothesis

Two questions determine whether any axis is meaningful: is the ligand
distinctively cardiac, and is the receptor tumour-enriched?

**Ligands.** Across 52 GTEx tissues, only **2 of 42** floor-passing
HF-upregulated secreted ligands are cardiac-enriched and abundant: NPPA (heart
rank 1, 35.8 TPM in left ventricle) and NPPB (rank 1, 26.8 TPM); GDF6 is
relatively cardiac-enriched but at 0.18 TPM is too scarce to be a credible
circulating signal. The remaining 39 peak elsewhere — GDF15 in kidney cortex,
CCL5 in whole blood, WNT9A in sigmoid colon (rank 29 of 52), MDK in ovary. The
failing heart is a distinctive source only for the natriuretic peptides.

**Receptors.** Comparison used two independent strategies, because TCGA
(RSEM normalized counts) and GTEx (TPM) cannot be placed on a common scale: a
background-corrected percentile-rank comparison against GTEx, and a genuinely
unit-matched log₂ fold change against 710 TCGA adjacent-normal samples streamed
from the same source matrix (15 tumour types with ≥10 normals). The two agree at
ρ=0.545 (p=8.9×10⁻⁹⁴, n=1,200 pairs), validating the rank proxy without equating
units. Of 78 receptors compared, 11 are enriched under both methods, 25 are
comparable to normal tissue, and 15 are depleted.

**The decisive result:** NPR1 — the receptor for the only two genuinely
cardiac-enriched ligands — is significantly **depleted** in tumour relative to
adjacent normal in **all 15** unit-matched tumour types (median
background-corrected log₂FC −1.91, all FDR<0.05), and independently fails the
expression floor in 26 of 30 types. Applying both criteria strictly (ligand
cardiac-enriched and abundant; receptor tumour-enriched and floor-passing in a
top-6 tumour type), **zero of 160 floor-passing interactions survive.**

The axis with a distinctive cardiac source has the least tumour-enriched
receptor, and the axes with tumour-enriched receptors are driven by ligands the
heart does not distinctively secrete.

![Figure 5](art_caffe368-443d-42c9-8c25-5ff4d6d020cc)

**Figure 5.** (a) Only NPPA and NPPB are cardiac-enriched among 42 HF ligands.
(b) NPR1 is depleted in tumour versus adjacent normal in every testable type.

### 2.8 A falsifiable experiment, despite the negative result

The screen's most testable axis is WNT9A→FZD7, chosen from the data: WNT is the
largest pathway contributor, and FZD7 has the widest dynamic range of any WNT9A
receptor in DepMap Uterus (the lineage mapping to UCS and UCEC). We specify a
4-versus-4 panel — responders JHUEM3, SNU685, EMTOKA, SNU1077; controls HEC108,
HOUAI, HEC151, HEC251 — with 48.7-fold mean linear TPM separation in FZD7.
Receptor level is not a library-size artefact (FZD7 vs housekeeping mean
ρ=0.131, p=0.447; arms indistinguishable in housekeeping expression, p=1.000),
responders have lower baseline WNT target activity so induction has headroom, and
FZD7 knockout does not kill these lines (gene effect +0.058 to +0.329), making a
knockdown arm interpretable as a signalling-competence manipulation. A
subtype-matched endometrial-carcinoma stratum is included. Note that WNT9A is
*not* cardiac-enriched, so this tests receptor-mediated signalling competence,
not a heart-specific mechanism.

---

## 3. What this study supports, and what it does not

**Supported by the data reported here**

1. A human HF transcriptomic signature replicates across two independent LV cohorts (concordance 0.881; permutation control 0 hits).
2. An absolute expression floor is necessary for a receptor-compatibility ranking to be receptor-specific, and quantitatively so (null ρ 0.363 → 0.107).
3. Specificity nulls must be expression-matched; a naive pool misjudges specificity (ρ=0.323).
4. Among floor-passing interactions, gynaecologic and renal tumours rank highest, with UCS robust to resampling.
5. Only NPPA and NPPB are distinctively and abundantly cardiac among the candidate ligands.
6. NPR1 is depleted in tumour relative to adjacent normal in all 15 testable types.
7. No nominated receptor shows compatibility-correlated CRISPR dependency.

**Not supported, and not claimed**

1. That heart failure causes cancer, in humans or in these data. No dataset here links HF status to cancer incidence.
2. That the ranking identifies clinically susceptible cancer types. It ranks transcriptional compatibility, and reorders under double-centring (ρ=0.318).
3. That any specific axis mediates HF-driven tumour growth. Zero axes pass both specificity checks.
4. That receptor expression implies signalling. We measured neither ligand exposure, receptor occupancy, nor pathway activation.
5. That tumour susceptibility depends on microenvironment state, as our own hypothesis proposed (1 of 29 features at FDR<0.05).
6. That these receptors are drug targets. They are not dependencies.
7. Anything about protein-level or circulating ligand abundance. This is transcript-level throughout, and plasma concentration is what an endocrine mechanism would require.

---

## 4. Limitations

**Transcript-only.** Every measurement is mRNA. An endocrine mechanism requires
circulating protein, and cardiac transcript induction need not produce elevated
plasma ligand. Cell-specific knockout of CCN2 showed that a secreted cardiac
factor can act autocrine on neighbouring fibroblasts rather than distally
(PMID: 30040954).

**Bulk tissue on both sides.** TCGA and GTEx values mix malignant, stromal and
immune compartments, so receptor signal cannot be assigned to tumour cells.
Compartment markers were extracted but a deconvolution was not performed.

**Method heterogeneity across cohorts.** GSE116250 provides RPKM and GSE141910
log-scale normalized counts, so cross-cohort effect sizes are not strictly
comparable; we rely on rank concordance rather than pooled magnitudes. The
empirical-Bayes moderated *t* is our own implementation of a limma-trend
equivalent, not a call into limma.

**GTEx heart is non-failing.** Tissue specificity in healthy myocardium bounds
baseline distinctiveness, not HF inducibility. A ligand could be
non-cardiac-enriched at baseline yet become a major cardiac output in failure.

**Unit mismatch.** TCGA RSEM and GTEx TPM were never placed on a common scale;
cross-dataset direction claims are ordinal only. The unit-matched
adjacent-normal comparison is the stronger evidence and covers 15 of 30 types.

**Mapping judgement calls.** DepMap-lineage→TCGA and TCGA→GTEx mappings are
disclosed with per-row confidence (Tables 15, 20). UCS is the top-ranked type yet
has a poor normal-tissue counterpart, only 57 tumours, and shares a single DepMap
lineage with UCEC — its specificity result is the weakest in the table.
Adrenocortical carcinoma, first under double-centring, has one screened DepMap
line and is untestable there.

**TIMER2 not used.** The Thorsson immune-landscape table substituted for TIMER2.
To be precise about why: `compbio.cn` is reachable, but TIMER2 is served as an
interactive Shiny application with no bulk-download endpoint we could locate
(`/timer2/data/` and `/timer2/estimation.php` both return 404), so its
deconvolution estimates are not retrievable programmatically at pan-cancer scale.
Thorsson is the published source for the same immune and stromal quantities over
the same TCGA samples, but it is not an identical computation, and we did not
cross-validate the two.

**Power.** Dependency correlations use 22 lineage-level aggregates; only a large
effect would be detectable. The null result bounds effect size, it does not
exclude a small one.

---

## 5. Methods

Six reusable skills, each a `SKILL.md` plus a `kernel.py` of pure
stdlib/numpy/pandas/scipy code containing no platform API calls, so they run
under any agent harness or as plain scripts.

**`geo-bulk-de`** — fetches GEO supplementary matrices, parses sample
characteristics from the series matrix (group labels are parsed, never
hardcoded), and runs differential expression with an explicit empirical-Bayes
variance-shrinkage implementation. Preconditions raise on duplicate sample IDs,
groups smaller than two, and all-NA genes.

**`hpa-secretome`** — Human Protein Atlas download API; classifies genes into
core (secreted and plasma-detectable), extended, and none, with a rule-based
structural-ECM flag and an explicit exclusion log. Raises if the response is not
a parseable TSV with >10,000 rows, guarding against saving an error page.

**`cellchat-lr`** — parses `CellChatDB.human.rda` in Python and expands receptor
complexes to subunits, so a heteromeric receptor is scored on its constituents.

**`tcga-pancan`** — streams the ~1.9 GB PanCanAtlas matrix row-wise, never
loading it whole (8 GiB machine), maps aliquot barcodes to tumour type, and
restricts to primary solid tumours (`01`), excluding LAML, DLBC and THYM.

**`crosstalk-score`** — the two-factor model. `require_expression_floor()` raises
by default; `matched_random_null()` matches substitutes on expression decile;
`decompose()` separates availability from capacity; `bootstrap_ci()` and
`leave_one_out()` provide stability. Guardrails are assertions, not prose, so
they hold regardless of which model or script drives the skill.

**`depmap-lineage`** — resolves the DepMap release through the figshare API
(depmap.org returns a bot-verification page that a naive download would save as
if it were a CSV; `fetch_depmap_file` asserts the payload is not HTML), reads
only requested gene columns via `usecols`, and requires the lineage mapping to be
passed explicitly rather than inferred.

Data: GSE116250, GSE141910, HPA, CellChatDB, GDC PanCanAtlas
(EBPlusPlusAdjustPANCAN, Thorsson immune landscape), DepMap 24Q4 Public via
figshare, GTEx v8 via an expression connector. Seeds recorded (permutation
20260727, bootstrap 20260727, DE permutation 42). Full provenance and re-run
instructions in `REPRODUCTION.md`; all controls in `validation_summary.csv`.

---

## 6. Corrections made during internal audit

Every numeric claim in this manuscript was machine-checked against the saved
artifacts, and an independent review pass re-checked the prose against the tool
outputs behind it. Six errors were found and corrected; they are logged as rows
26–31 of `validation_summary.csv` rather than silently fixed:

1. The floor drops **four** PENK→opioid-receptor axes (OPRD1, OPRK1, OPRL1, OPRM1), not three. The whole OPIOID pathway is non-viable.
2. **Four** of 30 tumour types beat the specificity null at p<0.05, not six. Ovarian carcinoma ranks second by score but does not pass (p=0.060) — rank and receptor-specificity are not the same thing.
3. GAPDH's median across per-tumour-type medians is **69,134**, not ≈64,000. The floor itself (345.7) was always computed rather than typed, so no result changes.
4. `compbio.cn` **is** grantable; an earlier note claimed TIMER2 was "not grantable" having probed a different host and never requested the real one. TIMER2 is nonetheless unusable here for a different reason (Shiny app, no bulk endpoint).
5. A documented scale-check example in `geo-bulk-de` reported output values that a real run does not produce (min/max 1.00/22.39 versus the actual 1.82/21.03). The qualitative conclusion it supports — that this "RAW" series ships log₂-scale values, not counts — is unaffected and was independently confirmed.
6. The claim that the skill catalog held "only Anthropic built-ins" was overstated; one unrelated personal skill pre-existed. None of the six analysis skills did.

None of these alter the study's conclusions, and one (item 2) makes the primary
ranking claim weaker than first written.

## 7. Data and code availability

Twenty-two numbered tables (`TABLE_INDEX.csv`), five figures
(`FIGURE_INDEX.csv`), 31 controls and audit corrections with verdicts
(`validation_summary.csv`), 107 verified citations
(`literature_citations.csv`), and claims the literature does *not* support
(`literature_gaps.json`). The six skills are published and loadable; each
documents standalone use without an agent.
