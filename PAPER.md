# An Agent-Skill Pipeline Finds No Match Between Heart-Failure Signals and Tumor Receptors in 30 Cancer Types

Nicolas Avelar¹, Ruidong Zhang², Yi Pan²

¹ North Oconee High School, Bogart, GA, USA
² University of Georgia, Athens, GA, USA

## Highlights

- A hierarchical multi-agent pipeline (1 orchestrator, 23 sub-agents, 6 deterministic analytical “skills”) screened receptor compatibility between heart-failure signals and 30 TCGA tumor types.

- No axis combined a distinctively cardiac ligand with a tumor-enriched receptor: NPR1, the sole receptor for the two cardiac-enriched ligands NPPA/NPPB, is depleted in tumor relative to matched adjacent normal in all 15 testable types.

- An expression floor was load-bearing: without it, expression-matched random receptor sets reproduced the tumor-type ranking (ρ = 0.363, 21.5% of draws > 0.5); with it, ρ fell to 0.107 (1.7%).

- Re-scoring cardiac specificity by failure-induction recovered seven candidate axes (all pleiotrophin); three kill criteria fixed in advance withdrew all seven.

- Two of the three corrections above were discovered as pipeline errors and are now encoded as code-level assertions in reusable skills, not left as prose recommendations a future study could silently omit.

## Abstract

**Background.** Heart failure has been proposed to accelerate tumor growth through factors released into the circulation. Testing that proposition across many tumor types is primarily a data-integration problem, involving five public resources, three unit systems, and a methodological choice at every join. It is also the kind of problem an autonomous agent can now execute end to end, with the attendant risk that the agent’s methodological errors leave no trace in its output.

**Methods.** The screen was executed by a hierarchical multi-agent system comprising one orchestrating agent and 23 isolated sub-agents dispatched in five waves over 91 hours. Each sub-agent received a task statement rather than the orchestrator’s conversation and returned artifacts. The procedures they called are packaged as six reusable skills, each a documentation file an agent reads to decide whether the procedure fits its task together with a Python module that executes it: 4,858 lines, 114 functions, and 141 guardrail checks written as assertions that halt execution rather than as advice in a methods section. The modules contain no language-model call, no API key and no platform function, so identical inputs return identical numbers irrespective of which agent invoked them. A heart-failure signature was derived from two RNA-seq cohorts (GSE116250, n=64; GSE141910, n=360), re-tested in three microarray cohorts, intersected with the secretome and a ligand–receptor database, and scored against receptor expression in 9,538 TCGA primary solid tumors across 30 types.

**Results.** No axis has both a distinctively cardiac ligand and a receptor enriched in tumor. NPR1, the sole floor-passing receptor of the only two cardiac-enriched ligands NPPA and NPPB, is depleted in tumor relative to matched adjacent normal in all 15 testable types (median logFC −1.91, all FDR<0.05). Two guardrails proved load-bearing and both were written after the pipeline made the corresponding error. Scoring without an absolute expression floor allowed expression-matched random receptor sets to reproduce the tumor-type ranking at ρ=0.363, with 21.5% of draws exceeding ρ=0.5; a floor at 0.005×GAPDH reduced these to ρ=0.107 and 1.7%. An unmatched random pool gave ρ=0.323, establishing that the null must be matched as well as floored. Re-scoring the source criterion against failure induction subsequently recovered seven candidate axes involving pleiotrophin, and kill criteria fixed before that analysis ran withdrew all seven, because 45 to 77% of candidate receptors satisfied the receptor criterion with no ligand involved. Reliability of the executing system is reported as measured quantities: 19 of 23 tracks completed normally, and all four that did not had already persisted their results, so no analysis was lost.

**Conclusions.** The negative is reported as a cardio-oncology finding and applies to receptor-mediated signaling by ligands induced in the failing human ventricle, not to secreted mediators in general. Two requirements of the endocrine premise are asserted rather than measured here: myocardial transcript induction does not establish that a ligand is elevated in the circulation, and bulk tumor tissue does not identify which cell expresses a receptor. Receptor enrichment over adjacent normal is a specificity criterion rather than a necessary condition for response, and is available for 15 of the 30 types. The transferable contribution is the form the corrections took. A requirement written as an assertion inside a reusable skill cannot be silently omitted by the next study, whereas the same requirement written as a methods paragraph routinely is.

## Keywords

Cardio-oncology; heart failure; tumor microenvironment; ligand-receptor signaling; TCGA; secretome; natriuretic peptide; pleiotrophin; multi-agent systems; computational reproducibility

## Abbreviations

**HF** — Heart failure

**DCM / ICM / HCM** — Dilated / ischemic / hypertrophic cardiomyopathy

**TCGA** — The Cancer Genome Atlas

**GTEx** — Genotype-Tissue Expression project

**CPTAC** — Clinical Proteomic Tumor Analysis Consortium

**PDC** — Proteomic Data Commons

**TMT** — Tandem mass tag (quantitative proteomics)

**HPA** — Human Protein Atlas

**GEO** — Gene Expression Omnibus

**DepMap** — Cancer Dependency Map

**RSEM** — RNA-Seq by Expectation-Maximization (transcript quantification)

**RPKM / TPM / CPM** — Reads/Transcripts/Counts per kilobase million

**FDR** — False discovery rate (Benjamini-Hochberg)

**logFC** — Log fold change

**LR** — Ligand-receptor

**NPPA / NPPB** — Natriuretic peptide A- / B-type

**NPR1** — Natriuretic peptide receptor 1

**PTN** — Pleiotrophin

**AGT / AGTR1** — Angiotensinogen / angiotensin II receptor type 1

**CCN2 / CTGF** — Cellular communication network factor 2 / connective tissue growth factor

**MAGNet** — Myocardial Applied Genomics Network

**IPF** — Idiopathic pulmonary fibrosis

**NAFLD / NASH** — Non-alcoholic fatty liver disease / steatohepatitis

**CI** — Confidence interval

**ρ** — Spearman rank correlation coefficient

## 1. Introduction

Cardio-oncology has largely studied how cancer treatment damages the heart. The converse direction, whether established heart disease influences the course of cancer, has received attention more recently, and the experimental evidence for it is specific. Transplanting a failing heart into tumor-bearing mice accelerated tumor growth even though the recipient’s own circulation was otherwise unchanged, which points to a factor carried in the blood rather than to altered perfusion or oxygen delivery (Meijers 2018 [1]). Myocardial infarction accelerates breast cancer growth in mice through changes in the output of bone-marrow-derived myeloid cells (Koelwyn 2020 [2]; Avraham 2020 [3]).

Human evidence is weaker and openly conflicting. Several cohorts report excess cancer incidence following a heart-failure diagnosis (Hasin 2013 [4]; Banke 2016 [5]; Roderburg 2021 [6]). The largest study with adjudicated cancer outcomes and a 19.9-year median follow-up found no association (adjusted HR 1.05, 95% CI 0.86–1.29; Selvaraj 2018 [7]). In a national cohort, allowing a two-year lag between diagnoses, which removes cancers likely to have been present but undetected when heart failure was diagnosed, attenuated the hazard ratio from 1.64 to 1.09 (Kwak 2021 [8]). Conventional cardiovascular risk factors predict incident cancer about as well as they predict incident heart failure (van den Berg 2026 [9]), and clonal hematopoiesis of indeterminate potential, in which an expanded blood-cell clone carries a somatic mutation without overt disease, is a shared substrate for both (Jaiswal 2017 [10]; Fuster 2017 [11]). Even the direction of effect is unsettled: tumor-bearing mice subjected to pressure overload developed less cardiac hypertrophy and fibrosis than tumor-free controls (Awwad 2023 [12]). Of 107 verified records assembled for this work, 21 report findings inconsistent with a causal link from heart failure to cancer.

We do not attempt to settle causality, because no available dataset links heart-failure status to cancer incidence with molecular resolution. We address a narrower question that follows from the mouse experiments and can be tested with existing human data. If the failing heart releases signaling proteins, a tumor can respond to them only if it expresses the corresponding receptors. Two conditions must therefore hold simultaneously for any specific heart-to-tumor signal: the ligand must be released distinctively by the heart rather than by tissues generally, and its receptor must be present in tumor tissue, ideally at higher levels than in the matched normal tissue from which the tumor arose.

Treating the problem as receptor compatibility rather than causation makes it falsifiable, because it can fail in three separable ways. The ligand may not be cardiac-specific, in which case the heart is not a distinctive source. The receptor may not be transcribed in tumor tissue at all, in which case no signal can be received. Or the receptor may be present but no more abundant than in normal tissue, in which case there is no reason to expect a tumor-selective effect. We tested all three across 30 solid tumor types and report where each candidate axis failed.

There is a second reason to report this screen. It was executed by an autonomous agent, and the work involved is the kind now routinely delegated: five public resources in three unit systems, a join at every boundary, and a methodological choice at each join that determines the answer while leaving no trace in the output when it goes wrong. Two of those choices went wrong here, in ways that produced confident and entirely spurious rankings. Both were caught, and the question this paper takes seriously is what to do with a correction once it has been made. Recorded as a sentence in a methods section, it binds nobody, and the next study inherits the text rather than the constraint. We instead packaged the analysis as six reusable skills in which each correction is an assertion that halts execution, and we report both what that changed and what it cost.

## 2. Materials and methods

All analysis is transcript-level and uses public data. Parameter values given here are the values used, not library defaults.

Each step is implemented as a skill: a directory holding a documentation file and a Python module. The documentation is addressed to an agent deciding whether the procedure fits the task in front of it and carries a structured header naming what the module does and when to invoke it. The module is what executes, and is where the methodological requirements live as assertions. Six were written for this analysis, and each subsection below names the one that performs its step. The unit is deliberately larger than a function and smaller than a pipeline, corresponding to the amount of work that recurs across studies of similar shape. Section 2.10 gives the measured properties of the whole set.

### 2.1 Data sources

Retrieval is the geo-bulk-de and tcga-pancan skills. Heart-failure transcriptomes were GSE116250 (discovery; left-ventricular RNA-seq, submitter-supplied RPKM [13]; 37 dilated cardiomyopathy, 13 ischemic cardiomyopathy, 14 non-failing donors, with the primary contrast pooling the 50 failing hearts against the 14 non-failing) and GSE141910 (validation; MAGNet left ventricle, per-sample log-scale normalized values [14]; 166 dilated cardiomyopathy, 28 hypertrophic cardiomyopathy and 166 non-failing donors, with 6 peripartum cardiomyopathy samples excluded, giving 194 versus 166). Group labels were parsed from the series-matrix characteristics fields rather than hardcoded, since hardcoding is the mechanism by which a mislabeled cohort becomes a silent result.

Cross-platform replication (Section 3.10) used three further cohorts, all Affymetrix microarray and none used elsewhere: GSE5406 (194 failing, 16 non-failing; GPL96), GSE57345 (177 failing, 136 non-failing; GPL11532) and GSE1869 (31 failing, 6 non-failing; GPL96) [15][16]. Probe-level values were collapsed to genes by maximum mean intensity with a deterministic tie-break. Probes annotated to several genes were dropped rather than assigned to the first listed gene, and log scale was detected from the values rather than assumed.

Protein abundance (Section 3.11) came from eight CPTAC tandem-mass-tag proteome studies in the Proteomic Data Commons [17][18][19]: ccRCC (PDC000127), ovarian (PDC000110), breast (PDC000120), pancreatic (PDC000270), endometrial (PDC000125), lung adenocarcinoma (PDC000153), head and neck (PDC000221) and glioblastoma (PDC000204). Log ratios were taken against the within-plex reference, with aliquots labeled tumor or adjacent normal from the study biospecimen records and assigned to TMT plexes from the experimental design. The GPL9052 sub-series of GSE57345 is a stub carrying no value rows and 6 samples, and was not analyzed.

Tumor transcriptomes came from the GDC PanCanAtlas EBPlusPlusAdjustPANCAN gene-expression matrix. Values are RSEM normalized_count, batch-corrected and upper-quartile scaled; they are neither TPM nor FPKM [20][21][22][23], a distinction that becomes consequential in Section 3.3. Aliquot barcodes were parsed to sample type and study, and analysis retained primary solid tumors only (sample-type code 01), excluding LAML, DLBC and THYM, giving 9,538 tumors across 30 types. Adjacent normal tissue (sample-type code 11) was streamed from the same matrix: 710 samples mapped to a tumor type, of which 680 fall in the 15 types retaining at least 10 normals.

Supporting resources were the Human Protein Atlas protein-class and secretome annotations [24][25] (search_download.php, columns g,eg,pc,secl,scml,rnats,blconcms); CellChatDB human (CellChatDB.human.rda) [26]; GTEx v8 median TPM per tissue across 52 tissues [27]; the Thorsson immune-landscape feature tables, leukocyte fraction, ABSOLUTE purity and TCGA-CDR [28] clinical endpoints from GDC [29]; and DepMap 24Q4 Public (figshare article 27993248). Because depmap.org serves a bot-verification interstitial [30][31][32], the DepMap release is resolved through the figshare API and the payload is asserted not to be HTML, since storing an error page as data yields plausible but meaningless output. Non-cardiac fibrosis cohorts for Section 3.9 were GSE150910 and GSE134692 (lung, idiopathic pulmonary fibrosis), GSE126848 and GSE130970 (liver, NASH/NAFLD) and GSE142025 (kidney, diabetic nephropathy) [33][34][35][36]. Human heart single-nucleus atlases (Reichart dilated and arrhythmogenic cardiomyopathy; Kuppe myocardial infarction) were obtained via CELLxGENE Census [37][38][39].

### 2.2 Differential expression and cross-cohort replication

This section is the geo-bulk-de skill. Each cohort was analyzed separately, on its own scale, with an empirical-Bayes moderated t-test equivalent to limma-trend implemented directly in NumPy and SciPy [40][41][42]. Genes were retained when expression exceeded 1.0 (RPKM or CPM as appropriate) in at least 20% of the smaller group’s samples, and 65 genes carrying a sentinel value of 999999999.0 were detected and dropped from GSE116250. The prior degrees of freedom fitted to 3.04 in discovery and 4.82 in validation, and the variance trend used a 0.4 lowess fraction. Multiplicity was controlled by Benjamini–Hochberg FDR at 0.05 [43]. Ensembl identifiers were version-stripped before intersecting the cohorts, yielding 18,061 shared genes.

A gene entered the replicated signature when it reached FDR<0.05 in both cohorts with the same direction of effect. Replication was quantified by directional concordance, log fold-change Spearman correlation, and a hypergeometric overlap test. The overlap test is near-saturated here, since 44.6% and 66.2% of shared genes are significant in the two cohorts, so concordance and effect-size correlation are the informative metrics and the hypergeometric p is reported alongside its fold enrichment. As a negative control, the discovery group labels were permuted (seed 42) and the identical pipeline re-run against the real validation result.

The implementation was validated against canonical limma 3.66.0 with the criterion fixed beforehand: a log fold-change Pearson correlation of at least 0.99 and at least 95% overlap of the FDR<0.05 gene sets, in both cohorts. Agreement is exact to floating-point precision. Pearson and Spearman correlations on log fold change are 1.000000 in both cohorts, with maximum absolute differences of 4×10⁻¹⁴ and 8×10⁻¹⁴, and the FDR<0.05 sets agree at Jaccard 0.9991 and 0.9998. The empirical-Bayes hyperparameters were estimated independently and match (prior degrees of freedom 3.041 against 3.037, and 4.818 against 4.825). Applying the replication rule to the limma output recovers 5,036 genes sharing 5,033 with the signature reported here. Full comparison tables, the least-agreeing variant and the one comparison that could not be made are given in the supplement, and agreement is shown in Supplementary Figure S2.

### 2.3 Secretome annotation and ligand–receptor mapping

These steps are the hpa-secretome and cellchat-lr skills. Up-replicated genes were annotated against HPA protein classes and tiered as core (predicted secreted and plasma-detectable), extended (secreted only) or none. A rule-based regex flag excluded structural extracellular-matrix families, which are secreted but are not plausible long-range endocrine signals: collagens, laminins, fibrillins, fibulins, elastin and microfibril proteins, fibronectin, nidogens, tenascins, structural proteoglycans and cartilage matrix. Exclusions were logged individually rather than aggregated into a count. The HPA response was required to parse as TSV with more than 10,000 rows, which guards against storing an error page as data.

Surviving ligands were mapped to CellChatDB, with analysis restricted to the Secreted Signaling annotation class. Heteromeric receptors were expanded to their constituent subunits, so that a receptor complex is scored on its obligate parts rather than on a single representative gene.

### 2.4 Two-factor compatibility score

This section and Sections 2.5 to 2.8 are the crosstalk-score skill, where the methodological risk concentrates and where 13 of the 141 guardrail checks sit. For each ligand–receptor interaction and each tumor type:

**Receptor capacity.** Subunit expression is aggregated per interaction using a minimum rule, which encodes obligate-subunit semantics: a complex cannot signal above the level of its scarcest component. A geometric-mean rule was run as a sensitivity analysis.

**Expression floor (mandatory).** Aggregated capacity must exceed an absolute floor set at 0.005 × the median GAPDH level across per-tumor-type medians. Expressing the floor relative to a housekeeping gene makes it unit-agnostic and therefore portable across platforms. In this matrix the median GAPDH value is 69,134.5, giving a floor of 345.7 RSEM normalized_count. The validator raises by default and names the offending receptor; interactions failing it in every tumor type are dropped before scoring rather than down-weighted.

**Standardisation and weighting.** Floor-passing capacities are z-standardized per receptor across tumor types, then averaged across interactions weighted by ligand availability, defined as the mean replicated log fold change of the ligand across the two heart cohorts.

**Uncertainty and stability.** Bootstrap 95% confidence intervals over 400 resamples (seed 20260727); leave-one-ligand-out and leave-one-pathway-out re-ranking; the geometric-mean complex rule; an unweighted variant; and a double-centered variant removing both row and column means.

### 2.5 Specificity null

The load-bearing validator asks whether a ranking is receptor-specific or merely tracks a global per-tumor expression axis. Each real receptor subunit is replaced by a random gene matched on mean-expression decile, the score is recomputed, and the null ranking is correlated against the observed one, over 300 iterations at seed 20260727. The background pool comprises 200 genes drawn 20 per decile from equal-count deciles of transcriptome-wide median expression across the 9,538 retained tumors, excluding the target panel and its legacy aliases. Four configurations were run: with the floor (primary), without it, from a uniform rather than decile-matched pool, and double-centered. A per-tumor-type empirical p is the fraction of null draws whose score equals or exceeds the observed one.

### 2.6 Ligand and receptor specificity criteria

**Ligand (source) criterion.** A ligand is cardiac-distinctive if the heart ranks first among the 52 GTEx tissues by median TPM and reaches an abundance plausible for a circulating signal.

**Receptor (receiver) criterion.** Because TCGA RSEM and GTEx TPM cannot be placed on a common scale, two independent strategies were used and required to agree. The first is a background-corrected percentile-rank comparison, in which each gene is converted to a within-dataset percentile among the same 1,379 genes shared by both matrices, separately per tumor type and per GTEx tissue, and the percentile delta is placed against the distribution of the same delta for all 1,379 background genes in that pair; enriched means top decile and positive, depleted bottom decile and negative. The second is a genuinely unit-matched log fold change against TCGA adjacent normal from the same source matrix, with BH-FDR within tumor type. No absolute value was subtracted, divided or plotted across the two platforms, and cross-dataset direction claims are ordinal only.

An axis survives only if the ligand criterion and the receptor criterion both hold.

### 2.7 Induction-based source criteria

Because the criterion in Section 2.6 scores the healthy heart, it is structurally unable to detect a ligand induced only in failure, which is a material limitation for a hypothesis about failing myocardium. A corrected definition was pre-specified before any computation, with three criteria, each computed within a single dataset so that no cross-platform absolute comparison is implied. A1 requires replicated induction of at least 1.0 log unit at FDR<0.05 in both heart cohorts. A2 requires abundance of at least 0.005 of GAPDH inside each heart cohort, the same housekeeping-relative device used for the tumor floor. A3 requires the heart to rank in the upper half of the 52 GTEx tissues rather than first. The receptor side was left unchanged. A pass required at least one surviving axis, with a guard declaring the criterion uninformative above twenty.

### 2.8 Adversarial tests of recovered axes

Kill criteria were fixed before any computation, with the receptor base rate measured before thresholds were chosen so that they could not be set to flatter the result. Retention required all three tracks to survive.

**T1, receiver discrimination.** The ligand-free base rate at which the 81 floor-passing receptors satisfy the receiver criterion in at least one tumor type, declared uninformative at 60% or above; a matched-random null substituting the candidate receptor set with expression-decile-matched sets of equal size, 1,000 permutations at seed 20260729, requiring the observed count to exceed the null’s 95th percentile; and pairwise co-expression among the candidate receptors, where a mean Spearman ρ of 0.7 or above would indicate that one observation was being reported several times.

**T2, tissue generality.** The same moderated-t pipeline applied to five non-cardiac fibrosis and organ-failure cohorts, each cohort’s normalization kept separate and COL1A1 required as a positive control. The criterion fails if the ligand is significantly induced in at least two cohorts at a median of at least 50% of the cardiac effect. Cell of origin was resolved in two single-nucleus atlases by per-cell-type expression share and by donor-level pseudobulk tests within the dominant cell type.

**T3, literature.** A PubMed search for circulating measurements of the ligand in heart failure. Published evidence contradicting a cardiac plasma source is a kill criterion.

**Conjunction rule.** Where a hypothesis generated three predictions (Section 3.10), all three were required to hold in the pre-specified direction. A single significant result with the others opposing is recorded as a failure, and any secondary read conditional on the conjunction was not performed when the conjunction failed. Survival analysis used Cox models [44] with the receptor z-scored within tumor type and stratified by type, and copy-number loss was tested against 1,000 loci matched on chromosome arm and overlapping-segment length.

### 2.9 Dependency screen

This section is the depmap-lineage skill. Per-lineage CRISPR gene effect was averaged across cell lines within each DepMap OncotreeLineage, then correlated (Spearman, BH-FDR) against the compatibility score of the TCGA types to which that lineage maps. The mapping was passed explicitly rather than inferred, and at least 5 screened lines per lineage were required, retaining 22 of 31 lineages. A gene effect below −0.5 was treated as a dependency, and essential and non-essential genes were carried as parsing controls.

### 2.10 The six skills, measured

**Table 1.** The six analytical skills, their function, and their load-bearing guardrails. Line, function, and guardrail counts are totals across all six modules as measured by static scan (Section 2.10); the distribution across modules is given in the text.

- Skill (module)

- Function

- Representative guardrails

- geo-bulk-de

- GEO retrieval; empirical-Bayes moderated-t differential expression (limma-trend equivalent)

- Detects and drops sentinel-value placeholders; detects log scale rather than assuming it; parses group labels from series-matrix fields rather than hardcoding them

- hpa-secretome

- HPA protein-class annotation and secretion tiering (core / extended / none)

- Rejects a response that is not a >10,000-row TSV (guards against silently storing an HTML error page as data)

- cellchat-lr

- Ligand-receptor mapping via CellChatDB; receptor-complex subunit expansion

- Restricts mapping to the Secreted Signaling annotation class; scores obligate complexes on their scarcest subunit

- tcga-pancan

- Memory-bounded streaming of the PanCanAtlas expression matrix; sample-type parsing

- Excludes non-solid-tumor types (LAML, DLBC, THYM); separates primary tumor from adjacent-normal by aliquot barcode

- crosstalk-score

- Two-factor (availability × capacity) compatibility scoring

- **Mandatory expression floor** (0.005×GAPDH) before any receptor is scored; **expression-matched null construction** (not a naive random pool); **conjunction rule** that raises rather than records a partial pass. Carries 13 of the 141 total guardrail checks.

- depmap-lineage

- CRISPR dependency screen via DepMap, aggregated per lineage

- Resolves the DepMap release via the figshare API and asserts the payload is not HTML; requires an explicit (not inferred) lineage-to-tumor-type mapping and a minimum of 5 screened lines per lineage

The six modules total 4,858 lines of Python containing 114 functions, with 1,522 lines of accompanying documentation. Their third-party dependencies are NumPy, pandas, SciPy, statsmodels and requests [45][46], with rdata required for one optional file reader; everything else is standard library. They contain no language-model call, no API key and no platform function, and a scan of all six executable modules for those returns zero in every category. That property is what makes the rest of this paper checkable, since every figure and table below is the output of code no model participated in executing, and the modules return the same numbers whichever agent invoked them or whether one did at all. The claim covers the module rather than the whole skill, since the documentation is addressed to an agent and is not executable, and no result here depends on how an agent reads it.

Across the six there are 141 executable guardrail checks: 78 assertions and 63 explicit raises. The distribution is informative rather than uniform. The scoring module, where the methodological risk concentrates, carries 13; the retrieval module, which handles heterogeneous public data and must reject malformed input in many ways, carries 47. Three are load-bearing for the results below and are named where they act. The expression floor refuses to score a receptor whose abundance falls below an absolute threshold, and raises an error naming the gene (Section 3.3). The null constructor refuses to draw a background set not matched on expression level (Section 3.3). The conjunction rule refuses to record a partial conjunction as a pass, so an axis satisfying two of three criteria cannot be counted as surviving (Section 3.9).

The distinction this draws is between a check a reader must choose to honor and a check a caller cannot avoid. A reader of a methods section can omit the floor without producing anything that looks wrong. A caller of this module cannot, because omitting it produces an exception rather than a number.

Random seeds were fixed as follows: 20260727 for the specificity permutation and the bootstrap, 42 for the differential-expression label permutation, 20260728 for the induction-criteria phase, and 20260729 for the adversarial receptor null. All computation used double precision.

### 2.11 The multi-agent system that ran the analysis

The analysis was executed by a hierarchical multi-agent system rather than by a single agent working serially. Because the call for reliability evidence in this area is specific, we report the topology and its failures as measured quantities rather than as description.

One orchestrating agent held the research question and the artifact store, and independent sub-agents were dispatched to tracks that could proceed without each other’s intermediate state. Twenty-three sub-agent frames ran in five waves over 91.4 hours of wall-clock time (Fig. 1). Wave 1 split the four data-acquisition tracks that share no inputs (heart-failure signature, secretome and ligand–receptor mapping, tumor expression axis, literature corpus) together with the dependency screen and normal-tissue comparison. Wave 2 ran five follow-up tests of the negative. Wave 3 ran the two orthogonal-molecule tracks, proteome and platform replication. Wave 4 ran six verification tracks, including two independent audits of the manuscript’s own claims. Wave 5 ran the implementation validation, the sensitivity sweeps and the citation expansion.

Each sub-agent received a task statement and a context summary rather than the orchestrator’s conversation, and returned results as artifacts plus a structured output object. That isolation is the reason the split is worth making, since a sub-agent cannot inherit an assumption the orchestrator formed earlier, and two of the wave-4 audit tracks were dispatched specifically to check claims the orchestrator had already written down.

![Figure 1](art_9a2d30bd-e8dc-4c0b-995a-bf881ba0d38f)

**Figure 1.** The multi-agent system that executed the analysis. (a) Wall-clock duration of each of the 23 sub-agent tracks, grouped into the five dispatch waves and colored by terminal status. The two longest tracks (Proteome, Platforms) include time parked on user approvals for network access. All four non-completing tracks had already written their results to the artifact store, and every one of those artifacts is cited in this paper. (b) Token consumption across the fleet. Cost is context rather than generation, which is why re-running the analysis from the six skills requires no model at all. Per-track figures are given in delegation_record.csv.

Nineteen of twenty-three tracks reached completion normally. Three terminated in a failed state and one was cancelled after 199 minutes when it exceeded its memory budget on an 8 GiB machine, and was re-dispatched with a streaming implementation. The instructive part is what happened to the work. All four non-completing tracks had already written their results to the artifact store, and every one of those artifacts is cited in this paper. Failure occurred at the reporting boundary rather than in the analysis, because results were persisted as artifacts as they were produced rather than returned only in a final message. A system that returns work only on successful completion would have lost four tracks and roughly 27 hours of compute.

The sub-agent fleet consumed $164 at 2026 API prices, with an input:output token ratio of 106:1, so essentially all of the cost is context rather than generation. Median track duration was 44 minutes. This is the figure that makes the division of Section 4.6 practical rather than rhetorical, since the deterministic modules those agents called cost nothing per invocation, and re-running the entire computational analysis from the six skills requires no model at all.

### 2.12 Statistical analysis conventions

This section collects, in one place, the statistical conventions applied throughout Sections 2 and 3; each is also introduced at the point it is first used.

**Software.** All statistical procedures were implemented directly in NumPy and SciPy, with pandas for data handling, statsmodels where noted, and requests for retrieval [45][46]; rdata was used for one optional file reader. No procedure calls a language model, an external API for computation, or a platform-specific function (Section 2.10).

**Differential expression.** Each cohort was analyzed on its own scale with an empirical-Bayes moderated t-test equivalent to limma-trend [40-42], validated against canonical limma 3.66.0 with a criterion fixed in advance (Section 2.2). Multiplicity was controlled by Benjamini-Hochberg FDR at 0.05 within each cohort or tumor type [43].

**Correlation and replication.** Cross-cohort and cross-platform replication is reported as directional concordance and Spearman/Pearson correlation of log fold change among genes significant in both cohorts, not as a hypergeometric overlap test alone; the latter is reported alongside its fold enrichment because it saturates when a large fraction of genes are significant in both cohorts (Sections 2.2, 3.1, 3.10).

**Resampling and permutation.** Bootstrap confidence intervals (400 resamples) and permutation nulls (300-20,000 draws depending on the test) used fixed random seeds, documented per procedure in Section 2.10 (20260727 for the specificity permutation and bootstrap; 42 for the differential-expression label permutation; 20260728 for the induction-criteria phase; 20260729 for the adversarial receptor null). Permutation p-values are one-sided empirical estimates (fraction of null draws meeting or exceeding the observed statistic) unless stated otherwise, and resolution limits relative to the number of draws are reported explicitly where a p-value approaches 1/N (Section 3.3).

**Survival and genomic association.** Cox proportional-hazards models [44] were stratified by tumor type with the receptor z-scored within type; copy-number associations were tested against a matched null of 1,000 loci sampled on chromosome arm and segment length (Section 2.8).

**Batch and design checks.** Proteomic contrasts were tested for arm-batch confounding by chi-square test before any effect estimate was interpreted, and a study was excluded from testing rather than analyzed through a confound (Section 3.11).

**Pre-specification and alpha level.** All thresholds used to declare a result significant, an axis surviving, or a test uninformative were fixed before the analysis phase in which they were applied (Sections 2.5-2.8, 3.9), and are reported as fixed even where the observed value fell close to the threshold (e.g., Section 3.9, Track 2). Unless otherwise noted, statistical tests are two-sided at α = 0.05 after multiplicity correction; permutation and bootstrap procedures are described above.

## 3. Results

### 3.1 A replicated human heart-failure signature

Differential expression in GSE116250 (50 failing hearts, of which 37 dilated and 13 ischemic cardiomyopathy, against 14 non-failing) recovered the expected cardiac stress programme: NPPA logFC +3.91 and NPPB +3.97, with MYH6 (−0.61) and ATP2A2 (−0.44) down, as expected in failing myocardium. GSE116250 distributes RPKM rather than counts, so a limma-trend-equivalent moderated t on log(RPKM+1) was used, and 65 genes carrying a sentinel value of 999999999 were detected and removed before analysis.

Validation in the larger MAGNet cohort (GSE141910; 194 cardiomyopathy against 166 non-failing donors) replicated the signature. Of 18,061 shared genes, 5,717 were significant in both cohorts and 5,036 shared direction, a directional concordance of 0.881 that rises monotonically with discovery effect size, from 85% at |logFC|<0.5 to 95% above 1.0. Effect sizes correlated at Spearman ρ=0.541 across all shared genes and ρ=0.759 among those significant in both.

Two cautions matter. The hypergeometric overlap test returns p=2.1×10⁻³ but only 1.07× enrichment: with 8,052 and 11,949 significant genes respectively, the test is near-saturated and is not evidence of replication despite the p-value, so directional concordance and effect-size correlation are the informative metrics. Separately, a label-permutation control run through the identical pipeline returned 0 significant genes and ρ=−0.134 against the real validation result, confirming that the pipeline does not manufacture concordance (Fig. 2).

![Figure 2](art_32f64f02-cc63-4a6a-865e-9e0c0eb0c07b)

**Figure 2.** Replication of the heart-failure expression signature. (a) Discovery volcano, GSE116250. (b) Cross-cohort effect-size concordance among genes significant in both; red points disagree in direction. (c) Concordance increases with discovery effect size.

Pooling cardiomyopathy etiologies does not mask a secreted signal. Re-deriving differential expression within each etiology separately against the same controls, 12 of the 42 upregulated secreted ligands are induced in one etiology and not the other (POSTN, INHBA, HBEGF, CCL3, NTF3, LGALS9, CCL5, CCL4, BMP4, WNT10B, CX3CL1 and MSTN), but none is cardiac-enriched at GTEx rank 1 and only five are abundant, so none would have passed the criteria NPPA and NPPB passed. Label-permutation controls returned no significant genes in all four stratified contrasts. The ischemic and hypertrophic strata are small, at 13 and 28 samples, so a null there is weak evidence of absence rather than evidence of a null.

### 3.2 From signature to candidate crosstalk network

Of 3,084 replicated HF-upregulated genes, 295 are annotated as secreted by the Human Protein Atlas and 70 are CellChatDB ligands of any annotation class, 50 of them within Secreted Signaling. Restricting to Secreted Signaling, the only annotation class mechanistically plausible at a distance, gave 188 interactions from 51 ligands across 38 pathways, with receptor complexes expanded to 103 constituent subunits. Ligand availability was computed as the mean replicated effect size scaled by secretion tier; the highest were NPPA (3.05), NPPB (3.03), PENK (2.85) and WNT9A (2.32).

### 3.3 The expression floor is the load-bearing methodological choice

TCGA PanCanAtlas values are RSEM normalized counts rather than TPM. The median GAPDH value across the 30 per-tumor-type medians is 69,134, so any floor expressed in transcript-per-million terms is meaningless in this matrix, and a familiar TPM threshold imported from elsewhere would return a number that carried no information. The floor was therefore defined relative to a housekeeping gene: 0.005×GAPDH = 345.7 in these units.

Twenty-eight of 188 interactions fall below it in every one of the 30 tumor types. These include all four proenkephalin-to-opioid-receptor axes (OPRD1, OPRK1, OPRL1, OPRM1), where OPRD1 reaches at most 4×10⁻⁴ of GAPDH, OPRM1 5×10⁻⁴ and OPRK1 1.4×10⁻³. PENK is the third most available HF ligand, so an unfiltered analysis would have nominated the opioid axis prominently, whereas the entire pathway is non-viable on receptor-expression grounds.

The floor is not hygiene; it is what makes the score interpretable. Under an expression-matched random receptor null (300 permutations, seed 20260727), the mean ρ against the observed ranking and the fraction of draws exceeding ρ=0.5 were as follows.

No floor applied: ρ = 0.363, with 21.5% of draws exceeding 0.5.

Floor applied (primary): ρ = 0.107, with 1.7%.

Floor applied, double-centered: ρ = 0.034, with 1.0%.

Unmatched (naive) random pool: ρ = 0.323, with 17.0%.

Two conclusions follow. The global expression axis confound that afflicts co-expression rankings is substantially caused by scoring untranscribed receptors; remove them and the ranking becomes receptor-specific. Separately, the null pool must itself be expression-matched, since a naive pool gives ρ=0.323 and would wrongly condemn a specific score (Fig. 3).

![Figure 3](art_62ff7487-a8e3-476a-83df-92f5b5f22f80)

**Figure 3.** The expression floor and the two-factor decomposition. (a) Floor audit; 28 interactions (red) never clear the floor in any tumor type. (b) Two-factor decomposition. Availability is a property of the source and is constant across receivers, so only capacity can explain differences between tumor types.

The floor value is not tuned, but neither does it sit on a broad plateau. A pre-registered plateau test asked whether, across a fourfold window from 0.0025 to 0.01 of GAPDH, the null correlation stayed below 0.20, the count of types beating the null stayed within one of four, and pass-set membership changed by at most one member. It fails two of the three. The pass set is identical at 0.0025 and 0.005 (UCS, BRCA, KIRC, PAAD), so the published value is not perched on a knife edge. However, the null correlation falls monotonically across the whole sweep with no flat region (ρ = 0.295, 0.201, 0.113, −0.021, −0.161 and −0.184 at fractions 0.001, 0.0025, 0.005, 0.01, 0.025 and 0.05), and a twofold move to 0.01 changes the pass set to BRCA, KIRC, OV, THCA and UCS, gaining two types and losing one. The defensible reading is that a floor is necessary and any value in this range makes the score receptor-specific, while the identity of the four passing types is sensitive to the exact fraction. Combined with the background-pool result below, the count is better reported as four of 30 with a range of three to seven across defensible choices.

Three hundred permutations resolve the reported p-values. Because the smallest reported value, 0.003, is 1/300 and therefore at the resolution limit, the principal null was re-run at 1,000, 5,000 and 20,000 draws. The set passing p<0.05 is exactly UCS, BRCA, KIRC and PAAD at every depth, and no type crosses the threshold as depth increases. At 20,000 draws UCS sits at 0.0022 with room below it, so the published value is a real quantity rather than a floor artifact, and ovarian carcinoma reaches 0.055 without ever passing. Monte Carlo standard error for UCS falls from 0.0047 to 0.00033 over that range, so deeper permutation would not change any conclusion.

The decile-matched null is insensitive to bin width and sensitive to pool composition. Varying the matching from 2 to 40 bins moves the null correlation only between 0.113 and 0.143, always far below the 0.30 threshold at which the score would be judged non-specific, and the passing count stays within one of four. At 40 bins the correlation begins to rise as bins become too thin to draw from, which places the sound operating range at 10 to 20 deciles. Pool composition matters more than width: a stratified background pool gives seven passing types rather than four, which is the upper end of the range quoted above.

### 3.4 Pan-cancer ranking, and how much of it to believe

Scoring the 160 floor-passing interactions ranked uterine carcinosarcoma first (0.679, bootstrap 95% CI 0.400–0.914, rank SD 0.55 over 400 resamples), followed by ovarian serous cystadenocarcinoma (0.432), clear-cell renal carcinoma (0.377), thyroid (0.285), prostate (0.272) and pancreatic adenocarcinoma (0.244). Only 4 of 30 types beat the matched-random null at p<0.05: UCS (p=0.003), KIRC (p=0.013), BRCA (p=0.010) and PAAD (p=0.033). Ovarian carcinoma, ranked second by score, is not among them (p=0.060). Rank order and receptor-specificity are therefore not the same thing, and a high score alone is insufficient grounds to nominate a tumor type.

Stability controls pass. Leave-one-ligand-out preserves the ranking at minimum ρ=0.907 (worst case WNT9A) and leave-one-pathway-out at ρ=0.842 (WNT, the largest single contributor at 40 of the 160 floor-passing interactions, or 80 rows once receptor complexes are expanded to subunits). Switching the receptor-complex rule from minimum-of-subunits to geometric mean gives ρ=0.908, and removing availability weighting altogether gives ρ=0.922, so no single ligand, pathway or parameter choice drives the result.

One control does not pass. Double-centering, which removes additive row and column effects, reorders the ranking substantially, at ρ=0.318 against the primary, promoting adrenocortical carcinoma to first and demoting UCS to sixth. Ovarian carcinoma is the only type in the top two under both variants. The ranking is correction-sensitive, and mid-table positions should not be interpreted (Fig. 4).

![Figure 4](art_b2c7e5cf-eb22-48ce-a162-0ceeb6fc09e1)

**Figure 4.** Pan-cancer compatibility ranking against an expression-matched null. (a) Ranking with bootstrap CIs; the 4 blue types beat the matched-random null at p<0.05. (b) The floor is what makes the score receptor-specific.

### 3.5 The exploitation component of the hypothesis is unsupported

Our hypothesis held that susceptibility depends on the ability to receive and then exploit cardiac signals. Availability is constant across tumor types by construction, so only receiver capacity can explain differences. Testing capacity against 29 Thorsson immune-landscape features, substituting for TIMER2, which was unreachable, only one survives FDR correction: B-cell receptor score (ρ=0.571, FDR=0.028). Leukocyte fraction, stromal fraction, TGF-β score, proliferation and wound-healing signatures are all null. Receptor compatibility and microenvironment state are therefore largely independent axes, and we found no evidence for the exploitation component as formulated.

### 3.6 Receptor expression does not correspond to lineage dependency

Using DepMap 24Q4 Public (1,178 cell lines, resolved via figshare article 27993248 because depmap.org serves a bot interstitial), we correlated per-lineage CRISPR gene effect against compatibility for 81 floor-passing receptors across 22 mapped lineages. Positive controls confirm correct parsing: as pan-cell-line means, RPL3 −2.53, PLK1 −2.69, POLR2B −2.42 and EEF2 −2.35, with non-essential MYT1 −0.02 and HBB +0.16. Aggregated per lineage instead, the same genes give −2.55, −2.67, −2.40 and −2.37. The two aggregations are reported separately because only the per-lineage values enter the correlation below.

No receptor is significant after FDR correction (0 of 81, 3 nominal). More pointedly, the receptors on the leading axes are not dependencies anywhere. NPR1 has a pan-lineage mean gene effect of −0.163, and every WNT9A receptor subunit lies between LRP5 (−0.163) and FZD6 (+0.146), with no lineage cell crossing −0.5. The eight receptors that are dependencies (NCL, ITGAV, EGFR, ITGB5, ITGB1, ERBB2, FGFR1, IGF1R) are pan-essential or canonical adhesion and RTK genes whose essentiality is lineage-driven and unrelated to HF ligand availability.

This is a meaningful distinction rather than a disappointment. Receptor expression licenses a cell to receive a signal, but its loss need not impair proliferation in vitro. The nominated axes are candidate signaling conduits rather than drug targets, and should be tested by ligand stimulation with a pathway reporter rather than by viability knockout (Fig. 5).

![Figure 5](art_19780104-a0e1-43b7-8cf4-58699de2941e)

**Figure 5.** Receptor expression and CRISPR dependency. (a) HF-axis receptors cluster at zero gene effect while essential controls sit near −2.4. (b) No receptor survives FDR. (c) Proposed uterine panel.

### 3.7 Specificity checks under a healthy-heart ligand filter

Two questions determine whether any axis is meaningful: whether the ligand is distinctively cardiac, and whether the receptor is tumor-enriched.

Across 52 GTEx tissues, only 2 of 42 floor-passing HF-upregulated secreted ligands are both cardiac-enriched and abundant: NPPA (heart rank 1, 35.8 TPM in left ventricle) and NPPB (rank 1, 26.8 TPM). GDF6 is relatively cardiac-enriched but at 0.18 TPM is too scarce to be a credible circulating signal. The other 39 peak elsewhere, with GDF15 in kidney cortex, CCL5 in whole blood, WNT9A in sigmoid colon (rank 29 of 52) and MDK in ovary. The failing heart is a distinctive source only for the natriuretic peptides.

The negative depends on the specificity metric rather than on the rank depth. Because heart ranking first is one of several defensible definitions, the criterion was re-scored under six of them, with GTEx v8 medians re-fetched and the 52-tissue set reproduced exactly. The three rank-based definitions are interchangeable: heart rank 1, top 3 and top 5 all give NPPA and NPPB as the only ligands that are both specific and abundant, and all give zero surviving axes. A fold-enrichment criterion is not interchangeable with them. At a fivefold threshold AGT joins as a third specific and abundant ligand (heart median 138.2 TPM against an all-tissue median 19.1, fold 7.2), though it fails every rank cutoff at best rank 11, and the axis count stays at zero because the AGT receptor AGTR1 is enriched in no tumor type. The two permissive definitions break the result: requiring only that heart sit above the tissue median admits 19 further ligands, 14 of them abundant, and returns 29 axes with a tumor-enriched receptor, while a twofold criterion returns 14. The negative is therefore a claim about strictly cardiac-restricted ligands. Under a cardiac-expressed-above-the-tissue-median reading the screen returns candidates, and which reading is correct is a question about what an endocrine source must look like rather than a question the data settle.

For receptors, comparison used two independent strategies, because TCGA (RSEM normalized counts) and GTEx (TPM) cannot be placed on a common scale: a background-corrected percentile-rank comparison against GTEx, and a genuinely unit-matched log fold change against TCGA adjacent-normal samples streamed from the same source matrix. Of those, 710 mapped to a tumor type, of which 680 fall in the 15 types retaining at least 10 normals and were analyzed. The two agree at ρ=0.545 (p=8.9×10⁻⁹⁴, n=1,200 pairs), validating the rank proxy without equating units. Of 78 receptors compared, 11 are enriched under both methods, 25 are comparable to normal tissue, and 15 are depleted.

The decisive result concerns NPR1. The receptor for the only two genuinely cardiac-enriched ligands is significantly depleted in tumor relative to adjacent normal in all 15 unit-matched tumor types (median background-corrected logFC −1.91, all FDR<0.05), and independently fails the expression floor in 24 of 30 types, passing only in KIRC, KIRP, MESO, OV, PCPG and SARC. Applying both criteria strictly, requiring the ligand to be cardiac-enriched and abundant and the receptor to be tumor-enriched and floor-passing in a top-6 tumor type, zero of the 146 ligand–receptor-subunit axes derived from the 160 floor-passing interactions survive.

The pattern is therefore inverted at both ends. The axis with a distinctive cardiac source has the least tumor-enriched receptor, and the axes with tumor-enriched receptors are driven by ligands the heart does not distinctively secrete (Fig. 6).

![Figure 6](art_caffe368-443d-42c9-8c25-5ff4d6d020cc)

**Figure 6.** Ligand cardiac specificity and receptor tumor enrichment. (a) Only NPPA and NPPB are cardiac-enriched among 42 HF ligands. (b) NPR1 is depleted in tumor versus adjacent normal in every testable type.

### 3.8 A falsifiable experiment, despite the negative result

A negative screen should still nominate something testable. The most tractable axis in the data is WNT9A to FZD7, since WNT is the largest pathway contributor and FZD7 has the widest dynamic range of any WNT9A receptor in DepMap Uterus. A 4-versus-4 cell-line panel with 48.7-fold mean separation in FZD7 is specified in the reproduction package, together with the controls that make it interpretable: receptor level is not a library-size artifact (ρ=0.131, p=0.447 against housekeeping mean across all 36 screened lines), and FZD7 knockout does not kill the lines, so a knockdown arm reads as a signaling-competence manipulation rather than a viability effect. WNT9A is not cardiac-enriched, so this tests whether the compatibility score predicts receptor-mediated signaling competence at all. It does not test a heart-specific mechanism, and it is offered as a check on the scoring skill rather than as a rescue of the hypothesis.

### 3.9 Correcting the specificity filter recovers axes, and the recovered axes do not survive

Both the corrected specificity criteria and the three kill criteria applied to the recovered axes were fixed before any computation in this phase, with the receptor base rate measured before thresholds were chosen so that they could not be set to flatter the result (Sections 2.7 and 2.8).

Section 3.7 reported zero surviving axes under a ligand filter that scored cardiac distinctiveness in GTEx left ventricle, which is non-failing myocardium. That filter is wrong for a hypothesis about the failing heart, and correcting it changes the count: re-scoring cardiac source by failure induction, with the receptor criteria unchanged, leaves 7 of 146 axes surviving instead of 0, all with pleiotrophin (PTN) as ligand.

That is a methodological finding about the filter rather than a biological finding about PTN, and the distinction is the substance of this section. PTN entered only because criterion A3 was relaxed from heart ranking first of 52 GTEx tissues to heart ranking in the upper half, and PTN sits at rank 23. A hit that appears when a criterion is loosened is precisely what loosening a criterion produces, so three kill criteria were fixed before any of the tests below ran, with PTN retained only if it survived all three. It failed all three.

**Table 2.** Summary of the three pre-specified kill criteria applied to the seven pleiotrophin (PTN) axes recovered after correcting the ligand-specificity filter (Section 2.8, 3.9). All three criteria, and the receptor base rate against which Track 1 was judged, were fixed before any of the tests below were run. Retention required all three to survive.

- Track

- Test

- Pre-specified threshold

- Observed result

- Verdict

- T1 – Receiver discrimination

- Ligand-free base rate at which candidate receptors satisfy the receiver criterion; matched-random null (1,000 permutations); pairwise receptor co-expression

- Base rate ≥60% declared uninformative; observed count must exceed the null’s 95th percentile

- 74.1% of all floor-passing receptors pass with no ligand involved (45-77% across aggregation rules); the 7 PTN receptors sit at the 31st percentile of a matched-random null (one-sided p = 0.688)

- **Fail**

- T2 – Tissue generality

- Moderated-t differential expression in 5 non-cardiac fibrosis/organ-failure cohorts; cell-of-origin in 2 single-nucleus cardiac atlases

- Fails if significantly induced in ≥2 cohorts at a median ≥50% of the cardiac effect

- Significant in 2 of 5 cohorts at 50.5% of the cardiac effect (triggers by 1.0%); direction reverses in the largest lung cohort; PTN is 62-71% fibroblast-derived with no significant per-cell induction in either atlas

- **Fail**

- T3 – Literature

- PubMed and full-text search for circulating PTN measured in a heart-failure cohort against controls

- Published evidence against a cardiac plasma source is a kill criterion

- No heart-failure cohort measurement exists in 609 verified records; the one study designed to detect cardiac PTN release (coronary-sinus sampling) explicitly rejected it as non-specific to myocardial injury

- **Fail**

- **Conjunction rule**

- All three tracks required to pass

- —

- 0 of 3 tracks pass

- **All 7 axes withdrawn**

#### The receptor test is not discriminating (Track 1)

Before asking whether the PTN receptors are special, we asked how special any receptor is. Of the 81 floor-passing receptors, 60 (74.1%) are enriched versus matched adjacent normal in at least one tumor type with no ligand involved. That figure is rule-dependent: scoring receptor complexes whole rather than per subunit gives 46 of 103 (44.7%) under the obligate-subunit rule and 79 of 103 (76.7%) under a geometric mean, so the base rate should be read as 45 to 77%. The conclusion it was raised to support, that most receptors pass with no ligand involved and that a six or seven axis count is therefore uninformative on its own, holds under every rule and is strongest under the geometric mean. The pre-specified threshold for declaring an axis count uninformative was 60%.

Seven receptors pass in more tumor types than any PTN receptor does: PLXNA3 in 13, PLXNA1 in 11, TNFRSF12A in 10, ERBB3 in 9, and ITGA6, F2R and ITGB4 in 8 each. The PTN receptors average 3.86 types against an all-receptor average of 2.52, placing them above the mean but at the 76th percentile of a distribution most members of which already pass. Exceeding the average within a population in which nearly every member passes is not evidence of specificity.

The formal test is worse than merely unremarkable. Substituting the seven PTN receptors with expression-decile-matched random receptor sets of the same size (1,000 permutations, seed 20260729) gives a null median of 30 enriched pairs against 27 observed. The PTN set sits at the 31st percentile of the null, with empirical one-sided p = 0.688, and three null configurations agree. The PTN receptor set performs slightly worse than random receptor sets matched for expression level.

One check passed. The seven receptors are not a single co-expressed family counted seven times, with mean pairwise Spearman ρ of −0.091 across the 30 tumor types and a maximum of 0.692. The axes are genuinely seven independent receptors, and random draws beat them.

#### PTN is an injury gene rather than a cardiac gene (Track 2)

PTN was tested in five non-cardiac organ-failure and fibrosis cohorts, each processed with the same moderated-t pipeline and each cohort’s normalization kept separate. COL1A1 rose in all five (logFC 1.17–2.52), confirming that every cohort captured a fibrotic signal. PTN was significantly induced in two: lung IPF (GSE134692, +0.646, FDR 4.6×10⁻³) and diabetic nephropathy (GSE142025, +0.871, FDR 2.5×10⁻⁵). The median of 0.7589 is 50.5% of the cardiac effect (1.5023) against a kill threshold of 50%. The criterion triggers by 1.0%, and is reported as triggered without adjustment, since revising a pre-specified threshold after observing which side of it a value falls would defeat the purpose of pre-specifying it.

The marginality does not rescue PTN, because the direction is inconsistent in a way the threshold does not capture. In the largest cohort tested, lung IPF with 103 cases and 103 controls (GSE150910), PTN is significantly down [47] (−0.556, FDR 2.5×10⁻⁵). The two lung cohorts disagree in sign. A ligand whose fibrosis response reverses between two cohorts of the same disease is not a robust injury signal in either direction, let alone a cardiac-specific one.

Cell of origin is more decisive than the induction test. In two independent human heart single-nucleus atlases, PTN is fibroblast-derived: 70.7% of PTN counts fall in fibroblasts in the Reichart DCM/ACM atlas and 61.6% in the Kuppe myocardial-infarction atlas, against 6.2% and 16.2% in cardiomyocytes, with fibroblast-to-cardiomyocyte mean expression ratios of 35.3× and 13.5×. The PTN of the failing heart is a wound-healing transcript from its stroma rather than a cardiomyocyte product.

The induction also does not survive the change of resolution. Per-fibroblast PTN does not rise in disease in either atlas (DCM +0.028, FDR 0.87; ACM −0.311, FDR 0.54; MI −0.430, FDR 0.30), nor does tissue-level pseudobulk (DCM +0.041, p = 0.72; MI −0.010, p = 0.68). The same MI fibroblast comparison detects POSTN at +1.77 and COL1A1 at +1.03, so the test had power. The bulk cardiac induction of logFC 1.50 that admitted PTN in the first place does not reproduce as a per-cell induction in either single-nucleus dataset. That discrepancy is unresolved. It may reflect aetiology or compositional differences between the bulk and single-nucleus cohorts, and it undermines the premise independently of any kill criterion.

#### The measurement the premise requires has never been made (Track 3)

Across 52 PubMed queries over a 991-record pleiotrophin corpus, together with a full-text scan of 478 cardiac plasma-proteomics papers, 609 records were retrieved and verified, and none measures circulating PTN in a heart-failure cohort against controls. Four measure circulating PTN in some human cardiac context, and two of those four point away from the premise. In the first (Circ Heart Fail 2024 [48]), plasma PTN measured by SomaScan is lower in hypertrophic cardiomyopathy than in three comparator cardiomyopathies, in both training and test sets, and after adjustment for 19 clinical parameters. In the second (Nat Biotechnol 2011 [49]), coronary-sinus sampling during planned myocardial infarction detected PTN release, and the authors explicitly eliminated it as not specific to myocardial injury, because catheterisation alone moved it.

The two supporting records fall outside heart failure: anthracycline cardiomyopathy in childhood cancer survivors [50], in which only 8 of 28 cases had heart failure, in a population where serum PTN is independently raised by cancer, and coronary collateral grade [51].

The secreted role of PTN in cancer is by contrast well established. The half of the axis that is supported is the tumor half. The cardiac-plasma half is unmeasured, and where adjacent measurements exist they run in the opposite direction. The one study designed to detect cardiac PTN release into blood found a signal and rejected it as an artifact of the procedure.

#### Verdict

The seven PTN axes are withdrawn as a positive finding. All three pre-specified tracks failed. The receptor test is passed by three-quarters of candidate receptors and the PTN set underperforms matched-random draws; PTN is a fibroblast injury transcript induced in lung and kidney fibrosis at half the cardiac effect, with no per-cell cardiac induction; and the plasma measurement on which the endocrine premise depends has never been made, with the nearest published evidence contradicting it.

The conclusion of the paper is therefore unchanged from Section 3.7, but it is now a stronger claim than it was. The original negative could have been an artifact of scoring cardiac specificity in healthy tissue. That objection has been tested directly: correcting the filter does recover axes, and the recovered axes do not survive. A negative that has been attacked on its most obvious weakness and held is worth more than one that has not (Fig. 7).

![Figure 7](art_df78cf23-161c-4d79-89c9-984651f1e760)

**Figure 7.** Cardiac specificity re-scored against failure induction. (a) Sequential attrition across the three pre-specified criteria; 3 of 51 ligands clear all three. (b) Which criteria each ligand fails, showing failure induction as the binding constraint. (c) Induction magnitude against baseline heart rank; PTN clears the induction and abundance criteria but sits 23rd of 52 tissues at baseline. (d) The seven surviving axes, all with PTN as ligand, by number of tumor types in which the receptor is enriched versus matched adjacent normal. Ligand criteria are computed within single datasets only; the receptor side is unchanged from Section 3.7. n = 51 ligands, 146 axes.

### 3.10 Cross-platform replication of the signature

The discovery and validation cohorts of Section 3.1 both ran on GPL16791, so the replication reported there lies within a single sequencing technology. To test whether the signature is platform-independent, it was re-derived in three cohorts not used previously, all Affymetrix microarray: GSE5406 (194 failing, 16 non-failing, GPL96), GSE57345 (177 failing, 136 non-failing, GPL11532) and GSE1869 (31 failing, 6 non-failing, GPL96). Thresholds were fixed before these analyses ran: directional concordance of at least 0.75 and log fold-change Spearman ρ of at least 0.40 among genes significant in both, with NPPA and NPPB required to rise in at least two of the three cohorts. Label permutation returned 0 significant genes in all 6 contrasts.

The two adequately powered cohorts replicate the signature. GSE5406 gives directional concordance 0.846 and ρ=0.713 among genes significant in both, and GSE57345 gives 0.819 and ρ=0.774. Correlation across all shared genes is lower in both (0.347 and 0.479), as expected, and is reported separately throughout because the two quantities are not interchangeable. GSE1869 falls below both thresholds (0.602, ρ=0.187), which with 6 control samples is uninformative rather than contradictory. One cohort below threshold does not reach the pre-specified kill criterion of two.

Two findings qualify Section 3.1. The shared platform inflated the apparent agreement: across all 10 cohort pairs, the 4 same-technology pairs have median concordance 0.918 and median ρ=0.722 among genes significant in both, against 0.745 and 0.545 for the 6 cross-technology pairs. The original pair sat at 0.881 and 0.759, inside the same-technology range. The concordance figure quoted in Section 3.1 is therefore about 0.17 units optimistic as an estimate of platform-independent agreement, and should be read as a within-technology number.

Second, and specific to this hypothesis, NPPB does not reproduce on either array platform, at −0.142 (FDR 0.90) in GSE5406 and +0.303 (FDR 0.27) in GSE57345, against +3.97 in RNA-seq discovery. This is not a probe-collapsing artifact, because NPPB is measured by a single probe on each platform (206801_at on GPL96, 7912520 on GPL11532), so no collapsing choice arises for it. Cross-platform support for the cardiac-ligand arm therefore rests mainly on NPPA, which behaves correctly in both (+1.856 and +1.733, FDR<10⁻⁸), as do MYH6 (−1.579 in GSE57345) and ATP2A2. Because the pre-specified criterion required both natriuretic peptides to rise in at least two of three cohorts and NPPB rises significantly in none, the criterion is not met, and the verdict on platform-independence is inconclusive rather than passed (Fig. 8).

![Figure 8](art_268435da-0544-4f30-be6e-abd2a231ba77)

**Figure 8.** Cross-platform replication of the heart-failure signature. (a) Agreement between each microarray cohort and the RNA-seq signature, against thresholds fixed before the analysis ran (dashed 0.75 for directional concordance, dotted 0.40 for log fold-change Spearman ρ among genes significant in both). GSE5406 and GSE57345 clear both; GSE1869, with 6 control samples, clears neither. One cohort below threshold does not reach the pre-specified kill criterion of two. (b) Directional concordance for all 10 pairwise cohort comparisons, split by whether the pair shares a sequencing technology. Bars are group medians. The original discovery/validation pair (star) sits inside the same-technology group because both cohorts ran on GPL16791. (c) The four canonical heart-failure markers in each cohort; stars mark FDR<0.05. NPPB is not significant on either array platform despite +3.97 in RNA-seq discovery, and is measured by a single probe on each platform, so this is not a probe-collapsing artifact. n = 194/16 (GSE5406), 177/136 (GSE57345), 31/6 (GSE1869) failing versus non-failing. Unit of replication is the donor heart.

The consequence for the argument is limited but should be stated. NPPB was one of the two ligands carrying the cardiac-specific arm of the hypothesis, and one of them is not measurable as induced on array platforms. This does not rescue the hypothesis, because the axis failed on the receptor side, where NPR1 is depleted in tumor in all 15 testable types, and NPPA, which shares that receptor, does replicate. It does mean that the ligand evidence rests on one peptide rather than two when platform is varied, and that a plasma measurement would be the way to settle whether NPPB is induced at all.

### 3.11 Protein-level test of the receptor result

Every measurement above is transcript abundance, so the central negative was tested in a different molecule. CPTAC tandem-mass-tag proteomics provides tumor and adjacent-normal protein abundance for eight of the tumor types ranked here. Criteria were fixed before any type other than ccRCC was examined: a gene is testable in a study only if it appears in that study’s quantified gene index and at least 10 aliquots per arm and 60% of each arm carry a value; the primary test is a paired comparison within TMT plex; NPR1 protein enriched in tumor in at least 3 of the 7 previously unexamined types would refute the negative, depletion in at least 5 would replicate it, and fewer than 4 testable types would abandon the test as underpowered. ccRCC was examined during the design of these controls and is reported as unblinded throughout.

Batch structure was checked in every study rather than assumed, and one study fails outright. In the breast cohort all 18 adjacent normals sit in two normal-only plexes, giving 0 mixed-arm plexes of 17 and 0 within-plex pairs (χ²=143.0, p=1.9×10⁻²²), so arm is perfectly confounded with batch, no batch-free contrast exists, and the study is untestable regardless of coverage. The glioblastoma cohort has only 10 normals, none paired within plex. Six studies have balanced designs, ccRCC among them (21 of 23 plexes mixed, χ²=15.7, p=0.83).

The test could not be run. NPR1 is testable in only 3 of the 8 studies, and in 2 of the 7 blind types, against a pre-specified minimum of 4. It is absent from the breast and glioblastoma gene indices entirely, and falls below the coverage floor in pancreatic (38%/33% of aliquots), lung (47%/50%) and head-and-neck (29%/39%). Absence of detection in a mass-spectrometry matrix is not absence of protein, so these are recorded as untestable rather than as null results. By the pre-specified rule, Track A is abandoned as underpowered, and no verdict on protein-level NPR1 is issued.

What the three testable types show is consistent with the transcript result but cannot carry it. NPR1 protein is lower in tumor in all three: ovarian −0.697 (13 pairs, FDR 0.017), endometrial −0.929 (20 pairs, FDR 0.0014) and ccRCC −0.114 (68 pairs, unblinded). All four pre-specified sensitivity analyses agree in direction, and none restores power. The mandatory base-rate check clears: across studies, 35.7% of quantified proteins are depleted in tumor, below the 40% bar above which the result would have been uninformative, and NPR1 sits at the 6.9th percentile of the genome-wide effect distribution in ovarian and the 7.0th in endometrial cancer. Three types agreeing is suggestive. It is not the five the criteria required, and reporting it as replication would be exactly the after-the-fact reasoning the pre-specification exists to prevent.

The pleiotrophin receptor panel was tested on the same footing and is withdrawn at protein level. Its count of receptors enriched in tumor fails to exceed the 95th percentile of an abundance-matched random null in all 6 studies where the null is defined; the closest is endometrial cancer at 4 observed against a null 95th percentile of 4 (empirical p=0.12), and in two studies the observed count falls below the null mean. Coverage again dominates, since 29 of 72 panel cells are untestable, 23 of them because the gene is absent from that study’s matrix. Pleiotrophin itself is quantified in all eight matrices and is mostly lower in tumor, which is the opposite of what the recovered axes required.

The conclusion is therefore unchanged but its scope is now explicit. The transcript-level negative stands, the protein-level test that could have overturned it is not answerable with current public proteomics at the coverage NPR1 achieves, and the measurement that would decide the question remains a plasma one.

### 3.12 A named mediator from the review literature runs in the opposite direction

A 2024 review of this field [55] nominates SerpinA3 as the principal secreted mediator linking heart failure to tumor growth, on the basis of the mouse study that established the phenomenon. SERPINA3 is a strong test case for the screen because it is annotated as secreted to blood at core tier, so it should have been available to the pipeline. It is in the replicated signature, but with the opposite sign to the one the hypothesis requires: logFC −2.64 in the discovery cohort (FDR 1.2×10⁻⁸) and −2.12 in validation (FDR 7.4×10⁻²¹), a replicated decrease in failing human left ventricle. The upregulation filter therefore excluded it, and it is also absent from CellChatDB as a ligand, so it would not have produced an interaction even had the direction agreed.

The discrepancy is most likely a species difference rather than a contradiction of the mouse work, since the induction was demonstrated in murine myocardium and in mouse plasma, whereas these are human explanted ventricles from patients with established end-stage disease, many on mechanical support. It nonetheless bears directly on the scope of the screen. The one mediator the review literature names most confidently is one this design could not have recovered, for two independent reasons, and that is worth stating rather than leaving to a reader who knows the field.

The same review [55] states that angiotensin II acting through AT1R promotes breast cancer proliferation. That axis is present in this screen and is informative about the specificity question in Section 3.7. AGT is upregulated in both cohorts, though weakly (logFC +0.22 and +0.13, FDR 0.042 and 0.047), and it is the third ligand admitted when cardiac specificity is defined by fold enrichment rather than rank. Its receptor is the reason the axis does not survive: AGTR1 is enriched in no tumor type against matched adjacent normal, and in the validation cohort it is itself strongly reduced in failing myocardium (logFC −1.16, FDR 4.7×10⁻²¹). A mechanism can be real in a treated mouse xenograft and still fail a compatibility screen built on human tumor transcript abundance, because the two tests ask different questions. This is the same distinction on which the negative result rests, visible here in a case the literature considers settled.
### 3.13 The closest human study nominates different ligands, and they fail different gates

A 2025 European Heart Journal report is the closest human work to this screen and reaches
a partially different answer, so its nominated ligands were run through the gates used
here rather than discussed in the abstract [56]. It used single-cell data from ischemic
and failing human hearts with NicheNet, implicated cardiac mesenchymal stromal cells as
the source, nominated POSTN, NGF and PDGF-family ligands, and reported increased
lung-cancer-cell proliferation in response to conditioned medium.

Those four ligand families fail at three different points, and none of the failures is a
disagreement about the underlying measurement. NGF is in the replicated signature but
runs the wrong way for the hypothesis, falling in failing ventricle in both cohorts
(log2FC −0.73 and −0.57, FDR 2.9×10⁻⁶ and 1.6×10⁻¹³), so it cannot be a heart-failure-
induced secreted signal in bulk myocardium whatever its behavior in a stromal
subpopulation. PDGFA, PDGFB and PDGFC are not in the replicated signature at all. POSTN
and PDGFD are replicated and induced (POSTN +1.25 and +0.84; PDGFD +0.58 and +0.34), pass
the receptor expression floor, and were scored — POSTN through ITGAV_ITGB3 and
ITGAV_ITGB5, PDGFD through PDGFRB — but both fail cardiac specificity: heart ranks 12th
and 14th of 52 GTEx tissues, with aorta the top tissue for each, and POSTN's heart
abundance is below the cross-tissue median (ratio 0.43). By comparison NPPA and NPPB rank
first with ratios of 12 and 88.

The disagreement is therefore about cell resolution and about what the source criterion
requires, not about whether these transcripts change. A single-cell design can attribute
a ligand to a stromal subpopulation that bulk ventricle averages away, which is a real
advantage and the same one exploited here for pleiotrophin (Section 3.9). What that
design does not do is ask whether the heart is a *distinctive* source: a fibroblast-
derived ligand that many organs' fibroblasts also make is not an endocrine cardiac
signal, and aorta out-expressing heart for both surviving candidates is exactly that
situation. Their conditioned-medium experiment is functional evidence this screen has no
counterpart to, and it establishes that failing-heart stromal secretions can affect
tumor-cell proliferation. It does not establish that the responsible ligand is
cardiac-specific or that it reaches a tumor through the circulation. The two results are
compatible: a paracrine, non-tissue-specific mechanism can be real while no
transcriptomically supported cardiac-endocrine axis exists.

## 4. Discussion

### 4.1 Where the hypothesis fails, and why the location matters

The two-factor framing was chosen because it can fail in identifiable places rather than merely returning a null. It failed in two of them, and the location matters more than the verdict. The receiver side is not the constraint: receptor transcripts for HF-upregulated ligands are present, often abundantly, across tumor types, and a compatibility ranking built on them is stable under resampling and under leaving out any single ligand or pathway. The constraint is specificity at both ends of the axis. Only the natriuretic peptides are distinctively and abundantly cardiac among the candidate ligands, and their sole receptor is the one receptor in the panel that is consistently depleted in tumor relative to matched adjacent normal. An axis needs both ends, and no axis has both.

That the stability of the score coexists with the failure of the axes is the substantive point. A compatibility ranking can be internally robust and still carry no specific signal, which is why the ranking is reported with its specificity null rather than on its own.

### 4.2 Three requirements that outlive this dataset, and the form they were written in

Three findings from Sections 3.3 and 3.9 are not about the heart. Each is a property of compatibility scoring in general, and each is now an assertion rather than a recommendation.

The first is that any score standardizing expression before checking that the gene is expressed at all inherits an artifact. The confound usually attributed to a global expression axis is substantially produced by scoring untranscribed receptors, since standardizing a gene sitting at instrument noise still yields large values, and those values are what a random draw reproduces.

The second is that the null used to detect this must itself be matched on expression. The naive pool and the floor address different failures, and neither substitutes for the other.

The third is a base rate. A screen that relaxes any criterion should report the rate at which the criterion it did not relax is passed, because without that number an axis count cannot be read as evidence about a particular ligand.

The sequence in Section 3.9, comprising a negative result, a filter identified as capable of having caused it, a correction that recovers hits, and a base-rate check that withdraws them, is generic enough that others will repeat it. That is the argument for writing the three requirements where they bind. Stated as prose here, they would be inherited by the next study as text, and a study that omits the floor would produce a result indistinguishable in form from one that honors it. Stated as assertions in the scoring skill, omitting the floor produces an exception naming the offending gene.

The claim is not that assertions are better than prose. Prose can express nuance that code cannot, and an assertion encodes a threshold rather than the reasoning that justified it, which may be wrong for another data type. The claim is that the two fail differently, and that for requirements of this kind, where the wrong answer resembles a right answer, the failure mode of the assertion is the one to prefer. An untranscribed receptor scored as though it were expressed does not announce itself; a missing exception does.

### 4.3 Bulk expression as a source measure

Two independent observations undercut bulk differential expression as a proxy for tissue output. The recovered ligand localises overwhelmingly to fibroblasts rather than cardiomyocytes, and its cardiac induction does not reproduce per-cell in either single-nucleus atlas, while established fibrotic markers do rise in the same comparison. Bulk differential expression in a tissue whose cellular composition changes with disease can reflect proportion shifts rather than regulation. The ligand-availability term of this screen is computed from bulk for every ligand and inherits that ambiguity generally, not only for the ligand that happened to be examined at single-cell resolution.

### 4.4 What would decide the question

The measurement the endocrine premise requires has never been made for the recovered ligand: circulating concentration in a heart-failure cohort against controls. The nearest published measurements run against it, including one study that detected the ligand during a planned myocardial infarction and rejected it as not specific to myocardial injury because catheterisation alone moved it. No reanalysis of expression data can substitute. Converting tissue mRNA into a circulating contribution requires translation efficiency, secretion rate, clearance and volume of distribution, none of which are available for these ligands, and that calculation was excluded rather than attempted with fabricated parameters.

The decisive experiments are therefore measurement rather than analysis: plasma proteomics across a heart-failure cohort with cancer-free controls, and, for any axis surviving such a screen, ligand stimulation with a pathway reporter rather than viability knockout. Receptor expression licenses a cell to receive a signal, and its loss need not impair proliferation in culture, which is why the dependency screen was negative for reasons unrelated to the hypothesis.

### 4.5 Limitations

**Transcript-only.** Every measurement is mRNA. An endocrine mechanism requires circulating protein, and cardiac transcript induction need not produce elevated plasma ligand. Cell-specific knockout of CCN2 showed that a secreted cardiac factor can act autocrine on neighbouring fibroblasts rather than distally [52].

**Platform dependence of the signature.** The discovery and validation cohorts both ran on GPL16791. Same-technology cohort pairs agree substantially better than cross-technology pairs (median concordance 0.918 against 0.745), so the 0.881 reported in Section 3.1 is a within-technology figure. The signature does reproduce on two independent array platforms, but NPPB does not, so the cardiac-ligand arm rests on NPPA alone once platform is varied (Section 3.10).

**Protein-level coverage, not protein-level agreement.** The proteomic test of NPR1 was abandoned as underpowered rather than passed or failed (Section 3.11). NPR1 is absent from two of eight CPTAC gene indices and below a 60% per-arm coverage floor in three more, and absence of detection in a mass-spectrometry matrix is not evidence of absent protein. The breast study is untestable for a different reason: all its adjacent normals sit in normal-only plexes, so batch is perfectly confounded with the comparison.

**Ambiguous probes were dropped, not resolved.** On the array platforms, 1,223 of 22,283 probes on GPL96 (5.5%) and 2,354 of 33,297 on GPL11532 (7.1%) map to more than one gene symbol. These were excluded rather than assigned to their first listed gene, which loses a small fraction of measurable genes (2.6% of signature symbols on GPL96, 4.9% on GPL11532) but avoids attributing the signal of a probe to an arbitrary one of its targets. No canonical marker is affected.

**Bulk tissue on both sides.** TCGA and GTEx values mix malignant, stromal and immune compartments, so receptor signal cannot be assigned to tumor cells. Compartment markers were extracted but a deconvolution was not performed.

**Method heterogeneity across cohorts.** GSE116250 provides RPKM and GSE141910 log-scale normalized counts, so cross-cohort effect sizes are not strictly comparable, and we rely on rank concordance rather than pooled magnitudes. The empirical-Bayes moderated t is our own implementation of a limma-trend equivalent rather than a call into limma [42].

**GTEx heart is non-failing.** Tissue specificity in healthy myocardium bounds baseline distinctiveness rather than HF inducibility. A ligand could be non-cardiac-enriched at baseline yet become a major cardiac output in failure. That possibility is no longer hypothetical, since Section 3.9 tests it directly and recovers PTN, which then fails on independent grounds.

**Unit mismatch.** TCGA RSEM and GTEx TPM were never placed on a common scale, and cross-dataset direction claims are ordinal only. The unit-matched adjacent-normal comparison is the stronger evidence and covers 15 of 30 types.

**Mapping judgment calls.** The DepMap-lineage to TCGA and TCGA to GTEx mappings are disclosed with per-row confidence (Tables 15, 20). UCS is the top-ranked type yet has a poor normal-tissue counterpart, only 57 tumors, and shares a single DepMap lineage with UCEC, so its specificity result is the weakest in the table. Adrenocortical carcinoma, first under double-centering, has one screened DepMap line and is untestable there.

**TIMER2 not used.** The Thorsson immune-landscape table substituted for TIMER2. To be precise about why, compbio.cn is reachable, but TIMER2 is served as an interactive Shiny application with no bulk-download endpoint we could locate (/timer2/data/ and /timer2/estimation.php both return 404), so its deconvolution estimates are not retrievable programmatically at pan-cancer scale. Thorsson is the published source for the same immune and stromal quantities over the same TCGA samples, but it is not an identical computation, and we did not cross-validate the two.

**Induction is not a systemic source.** The corrected specificity criterion in Section 3.9 admits ligands on induction magnitude without establishing that the failing heart dominates the circulating pool. Converting tissue mRNA into a circulating contribution requires parameters unavailable for these ligands, so that analysis was pre-emptively excluded rather than attempted with fabricated values.

**A relaxed criterion needs a base-rate guard.** The seven-axis result was recovered by relaxing one criterion and withdrawn by measuring the base rate of the criterion it was recovered against. Any screen that loosens a filter to recover hits should report the rate at which the unmodified half of the test passes; ours was 74%.

**Bulk induction may not be per-cell induction.** Cardiac PTN induction (logFC 1.50, replicated across two cohorts) does not reproduce as per-cell induction in two single-nucleus atlases. Bulk differential expression in a tissue whose composition changes with disease can reflect cell-proportion shifts rather than per-cell regulation, and the ligand-availability term of this screen is computed from bulk for every ligand, not only PTN.

**Power.** Dependency correlations use 22 lineage-level aggregates, so only a large effect would be detectable. The null result bounds effect size rather than excluding a small one.

### 4.6 What the model decided, and what the code enforces

Describing this pipeline as agentic invites a fair question: how much of the result depends on the language model that assembled it, and would a different or cheaper model have produced something else. The answer follows from where the work sits. None of the 4,858 lines contains a call to a language model, an API key or a platform function, so given the same inputs they return the same numbers whichever agent invoked them, and every figure and table here is output of that layer.

The contribution of the model was orchestration and judgment: deciding that ligand distinctiveness and receptor abundance are the two factors worth scoring, recognising that untranscribed receptors were inflating the null and that a floor was therefore required, choosing an expression-matched rather than a naive background, and requiring a base rate before the hits recovered by a relaxed filter could be believed. Those are the decisions this paper argues about. Once made, each was written into a module as an assertion, at which point it ceased to be a judgment and became a property of the software.

That division has a practical consequence for cost. The expensive capability is needed at the design stage, where the number of decisions is small, whereas the deterministic layer, where essentially all the compute goes, needs no model at all. A cheaper or smaller model driving these skills inherits the same guardrails, because a below-floor receptor raises an error regardless of which agent supplied the arguments. What a weaker model risks is not a corrupted number but a worse experiment: a poorly chosen comparison, or a filter relaxed without checking its base rate. Running the modules is cheap and safe; deciding what to run is where model quality matters, and it is also where a mistake is hardest to detect afterwards.

### 4.7 Reuse beyond this question

The same six skills apply wherever the secreted output of a source tissue is scored against receiver contexts, including any organ-crosstalk screen, any secretome-to-receptor ranking, and any compatibility score across many tissues or cell lines. Nothing in the scoring module depends on the source being heart or the receivers being tumors, since it consumes a feature matrix and a vector of source-derived weights.

What does not transfer automatically is the threshold. The floor is calibrated as a fraction of GAPDH in RSEM units, and a study on a different platform must recalibrate it rather than inherit the number, since the assertion enforces that a floor exists and was applied rather than that 0.005×GAPDH is correct everywhere. That is the honest limit of the mechanism: it makes a requirement unavoidable, not universally parameterised.

### 4.8 What the multi-agent split bought, and what it cost

Splitting the work bought three things, and the third changed a result. These are observations from one execution, not a controlled comparison: no single-agent baseline, no scripted-workflow baseline and no context-isolation ablation were run, so what follows describes what this configuration did and not what an alternative would have failed to do.

The first is throughput on independent tracks. The four wave-1 acquisition tracks share no inputs and ran concurrently; serially they would have consumed the sum of their durations rather than the maximum.

The second is context isolation as a correctness property rather than an efficiency one. Each sub-agent received a task statement and a summary, never the conversation of the orchestrator, so it could not inherit an assumption formed earlier. This is what made the wave-4 audits usable, since tracks were dispatched specifically to re-derive claims the orchestrator had already written into the manuscript, working from the result tables rather than from its prose, and they found errors. The adjacent-normal sample count was reported as 710 where the table sums to 680, because 710 counts samples before the minimum-normals filter. The number of tumor types in which NPR1 fails the expression floor was reported as 26 of 30 where the column gives 24. One survival claim mixed units, stating zero of 160 when the table holds 146 subunit-level axes and 160 is the complex-level count. A review track separately found a pre-specification justifying its matched-null requirement by citing the expression-floor contrast (ρ=0.363 against 0.107) rather than the matching contrast (ρ=0.323), which is a reasoning error rather than a typographical one. None of these changes a conclusion, and that is the point: they are the class of error that survives self-review, because an agent auditing its own prose in its own context has already been persuaded by it.

The third is that isolation is what allowed the kill criteria to bite. The adversarial phase of Section 3.9 was dispatched with its criteria in the task statement and no access to the reasoning that had produced the seven candidate axes, and it withdrew all seven. We do not claim that a single-context agent could not have done this, but the design removes the mechanism by which it would most plausibly have failed to.

The architecture claim this paper can support is therefore narrow. The audits found errors and the pre-specified criteria withdrew a positive result, which shows the guardrails did work; it does not show that 23 agents were required to produce it, and a comparative evaluation against a single-agent or scripted pipeline on the same task is the obvious next experiment. The costs are real. The fleet consumed $164, essentially all of it context rather than generation at an input:output ratio of 106:1, which is the characteristic cost profile of this architecture, since every sub-agent pays to be told what it needs to know. Four of 23 tracks did not complete normally, and one of those consumed 199 minutes before exceeding its memory budget. Handoff is the fragile boundary, because the work of a track survives a failure only if it was persisted as it went, which is a design requirement rather than a platform guarantee. Isolation also cuts both ways, since a sub-agent lacking context can make a locally sensible choice that is wrong for the study, which is why the output of every track was checked against the artifact it claimed to produce rather than accepted from its summary.

## 5. Conclusions

Across 30 solid tumor types, no ligand–receptor axis has both a distinctively cardiac source and a receptor enriched in tumor tissue, where enrichment is measured against unit-matched adjacent normal in the 15 types that have it. Enrichment is a specificity criterion: a tumor could respond through a receptor expressed adequately but not enriched, so this is evidence against a tumor-selective cardiac axis rather than against receptor-mediated response in general. The conclusion survives correction of the design choice most capable of having produced it. Scoring cardiac specificity in healthy myocardium cannot detect a ligand induced only in disease, and re-scoring by failure induction recovered seven candidate axes, all involving pleiotrophin. Criteria fixed in advance withdrew all seven, because the receptor criterion was satisfied by 60 of 81 candidate receptors with no ligand involved.

Two further tests were pre-specified and run after the negative was reached. The signature reproduces on independent microarray platforms, but NPPB does not, so the cardiac-ligand arm rests on NPPA alone once platform is varied. A protein-level test of the receptor result was attempted in CPTAC proteomics and abandoned under its own kill criterion, because NPR1 reaches the required coverage in too few tumor types. The three types where it is testable all show lower NPR1 protein in tumor, which is directional agreement with the transcript result but not independent confirmation of it: three types cannot validate a conclusion drawn across 15, and the analysis was abandoned precisely because it is underpowered. The pleiotrophin panel is withdrawn at protein level as it was at transcript level.

The hypothesis is not refuted in general. It remains tenable in the weaker form that heart failure raises circulating concentrations of factors secreted by many tissues, which would require neither a cardiac-specific ligand nor a tumor-specific receptor. Distinguishing that from the specific form tested here requires plasma measurements in heart-failure patients, and for pleiotrophin, the strongest candidate this screen produced, such measurements do not exist. Tumor proteomics does not substitute, since it measures the receiver rather than the circulating signal. The screen also cannot recover a mediator that acts without a curated cognate receptor, is regulated in the opposite direction in human end-stage disease than in the mouse model, or acts through protease activity rather than receptor binding, which is the case of SerpinA3 (Section 3.12). The negative applies to receptor-mediated signaling by ligands induced in the failing human ventricle rather than to secreted mediators in general.

What we would hand to the next group is not the negative but the instrument. The three requirements of Section 4.2 were each discovered as an error this pipeline made, and each now exists as an assertion inside a skill rather than as a paragraph here. That is the difference the paper argues for: a study built on the scoring module cannot omit the expression floor or substitute an unmatched null without raising an error that names the offending gene, whereas a study inheriting the same three findings as prose can omit all of them and produce output that looks correct.

The untested premise is reuse. We can demonstrate the mechanism within one study, since the kill criteria that withdrew the only positive result of this study were themselves evaluated against a floored, expression-matched background because the module enforced it. Demonstrating that the constraint binds a different group would require a second analysis built on the same skills, and that has not happened yet.

## 6. Future directions: integrating biomedical imaging

Every measurement in this screen is molecular: transcript abundance, protein abundance, or a dependency score. Imaging was not used, and its absence is a gap rather than a neutral omission, because several of the open questions in Sections 4.4 and 4.5 are more naturally imaging questions than sequencing questions.

**Cardiac phenotype as a continuous exposure.** Heart-failure status was treated here as a binary contrast (failing versus non-failing donor hearts), which discards the variation in severity, chamber geometry, and fibrotic burden that cardiac MRI and echocardiography quantify routinely in living patients. Late gadolinium enhancement and native T1 mapping index diffuse and focal myocardial fibrosis directly, and strain imaging indexes subclinical dysfunction before an ejection-fraction threshold is crossed. If a circulating ligand’s release scales with fibrotic burden rather than with a diagnostic label, an imaging-derived severity score is a better exposure variable than the case-control contrast used throughout this paper, and cohorts that pair cardiac imaging with plasma proteomics, rather than tissue transcriptomes, would let the compatibility score in Section 2.4 be regressed against a continuous phenotype instead of a group mean.

**Spatial and histological readouts of the tumor side.** The receptor side of the screen used bulk RSEM values that mix malignant, stromal, and immune compartments (Section 4.5, Bulk tissue on both sides), the same limitation that Section 3.9’s single-nucleus analysis exists to address on the cardiac side. TCGA carries diagnostic whole-slide histopathology images for the great majority of the tumors scored here, and spatial transcriptomics and imaging mass cytometry, of the kind already used for the Reichart and Kuppe cardiac atlases in this study (Section 2.1), exist for a growing number of tumor types. Neither was used to ask where, within a tumor, a floor-passing receptor is expressed. A receptor enriched in bulk tumor tissue but confined to stroma or infiltrating leukocytes, by analogy to the fibroblast-restricted pleiotrophin result in Section 3.9, would not support a tumor-cell-intrinsic signaling claim regardless of its bulk enrichment, and digital pathology or spatial platforms are the direct way to test that rather than infer it from cell-of-origin proxies.

**Imaging as a circulating-factor endpoint.** Section 4.4 identifies the decisive missing measurement as circulating ligand concentration in a heart-failure cohort, and that is a plasma-proteomics problem rather than an imaging one. But imaging enters downstream of it: were a candidate axis to survive plasma measurement, whether a tumor actually responds would be testable non-invasively. Longitudinal tumor growth-rate imaging (serial cross-sectional imaging in patients or PDX models), or PET tracers sensitive to a specific receptor’s signaling output rather than to generic glucose uptake, would let a nominated axis be tested for tumor-level consequence without a biopsy at every timepoint. Section 3.6 already makes the analogous point for genetic dependency, that receptor expression licenses a response without guaranteeing one; an imaging-based functional readout is the same argument applied to whole-tumor behavior instead of cell viability.

**Which heart failure?** This screen treated heart failure as one exposure, and every
limitation of that choice is a question a subsequent study can answer with existing
cohort designs rather than new biology.

**Severity is not binary, and the design assumed it was.** Cases were failing versus
non-failing explanted ventricles, which collapses a continuum onto a label. If a ligand's
release scales with fibrotic burden or wall stress rather than with a diagnostic
category, the case-control contrast used here is the wrong estimator and would attenuate
a real signal toward the null. The natriuretic peptides make the point concrete: NPPA
and NPPB are clinically graded biomarkers whose plasma concentration tracks filling
pressure across orders of magnitude, so a design that scores them as present-or-absent in
tissue discards precisely the axis on which they carry information. The imaging-derived
severity indices in the preceding subsection are one route to a continuous exposure;
NYHA class, natriuretic-peptide concentration, and ejection fraction are others already
recorded in most cardiology cohorts.

**Etiology was tested and did not change the conclusion, but only for three etiologies.**
Stratifying the 42 upregulated secreted ligands by cardiomyopathy subtype found 12 with
etiology-specific induction in at least one cohort, none of which is cardiac-enriched at
rank 1, and two of which clear the abundance bar set by NPPA and NPPB — so pooling did
not conceal a candidate axis (Section 3.1). That test covers dilated, ischemic, and
hypertrophic cardiomyopathy, with strata of 37, 13, and 28 patients, and a null in the
two small strata is low power rather than demonstrated absence. It says nothing about
the etiologies absent from these cohorts entirely. Pressure overload from aortic
stenosis, volume overload from mitral regurgitation, and inflammatory myocarditis impose
mechanically and immunologically distinct stresses on the ventricle, and there is no
reason to assume their secreted programs are the union of what dilated and ischemic
disease produce. Myocarditis in particular is an inflammatory state, and the one
mediator this screen was structurally unable to recover — SerpinA3, an acute-phase
protease inhibitor (Section 3.12) — is exactly the class of molecule an inflammatory
etiology would be expected to elevate.

**Stage is confounded with tissue availability.** Explanted ventricles come from
transplant recipients, so every case here is end-stage disease under maximal therapy.
Early-stage heart failure, which is where the epidemiological association with cancer
incidence would have to originate if the relationship is causal, is not represented in
any tissue cohort and cannot be, since those patients are not undergoing
cardiectomy. Endomyocardial biopsy series and circulating markers are the only accessible
windows onto that stage, and both trade tissue depth for accessibility. This is a
structural limit on the entire tissue-transcriptomic approach to this question, not a
gap in the present study.

**Treatment is an unmodeled exposure that acts on the same genes.** Patients supplying
these ventricles were on guideline-directed therapy, and the major drug classes
demonstrably alter the transcripts this screen scores: neurohormonal blockade changes
natriuretic-peptide signaling, and mechanical unloading by a left ventricular assist
device produces partial reverse remodeling with measurable transcriptional recovery.
What is measured here is therefore the treated failing ventricle, not the failing
ventricle. Whether the secreted program of interest is a property of the disease or a
property of the disease under therapy is an answerable question — paired pre- and
post-LVAD samples exist — and it should be answered before a circulating-ligand
hypothesis is tested prospectively in a treated population.

Each of these is a stratification the compatibility score can accept without
modification, because the scoring module consumes a feature matrix and a vector of
source-derived weights and is indifferent to how the source group was defined. The
guardrails travel with it: a severity-stratified or etiology-stratified rerun inherits
the expression floor and the matched null whether or not whoever runs it has read this
paper.

**What this implies for how medicine gets done.** Each stratification above shares a
shape, and it is the shape most translational questions now have: the biological
hypothesis is cheap to state, the data to test it already exist in five separate public
resources, and the work is integration under methodological choices that determine the
answer. That work is now delegable. What this study argues is that delegating it safely
depends less on the capability of the agent than on where the methodological knowledge
is written down. This screen made two errors that produced confident and entirely
spurious rankings, and neither would have been visible in the output; both were caught,
and both now halt execution rather than appear as advice a subsequent study can skip.

That has a specific implication for clinical translation. The pipeline nominated seven
candidate axes and then destroyed them, at no cost beyond compute, because a criterion
fixed in advance was enforced by code rather than by the analyst's memory. In a
translational chain the alternative is expensive: a nominated target that survives to a
validation cohort, an animal model, or a trial before the base rate that should have
withdrawn it is computed. Screens of this kind are cheap to run and cheap to believe,
and the discipline that makes them worth running is not the model doing the analysis but
the encoded refusal to score what should not be scored.

The pattern generalizes past cardio-oncology to any question with the same structure:
one organ's secreted output against many receiver contexts, a source phenotype
measurable at scale, and public reference data on both sides. Kidney disease and cancer
risk, hepatic injury and extrahepatic malignancy, and the systemic effects of chronic
inflammation all fit it. For those, what transfers from this work is not the negative
result but the requirement that a compatibility score check absolute abundance before
standardizing it, judge specificity against an expression-matched background, and
measure the base rate of any criterion it relaxes — three requirements that cost nothing
to honor and, on this dataset, changed the answer.

**Extending the skill architecture rather than the hypothesis.** None of these directions requires abandoning the two-factor framing of Section 1 or the guardrail architecture of Section 2.10; they require a seventh skill. An imaging-integration module would need the same properties enforced elsewhere in this pipeline: a deterministic, model-free feature extraction step (for example, radiomic or histological feature computation) whose output is checked against a floor and a matched null exactly as receptor expression is checked in Section 2.4, so that an imaging feature correlated with the compatibility score is held to the same specificity standard that excluded 87.5% of the pan-cancer ranking in Section 3.4 and withdrew all seven pleiotrophin axes in Section 3.9. The risk an imaging module would most plausibly reproduce is the one this paper spent the most space on: a feature that tracks a global axis, such as tumor size or slide staining intensity, rather than the receptor-specific signal it is meant to measure, unless the same floor-and-matched-null discipline is written into it as an assertion rather than left as a recommendation.

We have not built this module, and the discussion above is a specification of what a next study would need rather than a result. It is included because Section 5 hands the next group an instrument, not a conclusion, and imaging integration is where that instrument’s guardrails have not yet been tested.

## 7. Data and code availability

All data are public: GSE116250 and GSE141910 (Gene Expression Omnibus) [53][54], TCGA PanCanAtlas expression and the Thorsson immune landscape, GTEx v8 median tissue expression, Human Protein Atlas protein-class annotations, CellChatDB, and DepMap 24Q4. Accession numbers, endpoints and retrieval dates are given in Section 2.1.

Analysis code is available at https://github.com/nav92061/hf-crosstalk and is released as six skills, each comprising a documentation file plus a Python module: geo-bulk-de (GEO retrieval and moderated-t differential expression), hpa-secretome (protein-class annotation and secretion tiering), cellchat-lr (ligand–receptor parsing with receptor-complex expansion), tcga-pancan (memory-bounded streaming of the PanCanAtlas matrix), crosstalk-score (two-factor compatibility scoring with the floor, the matched-random null and the conjunction rule), and depmap-lineage (dependency screening with release resolution). The expression floor, the matched-random null and the conjunction rule are assertions that halt execution on violation rather than documented recommendations.

The skills carry no dependency on the environment used here: no language-model call, no API key, no platform function. They can be driven by another agent framework, or by a plain Python script with no agent at all, and the guardrails hold in every case because they are assertions in the code rather than instructions in a prompt. Continuous integration fails the build if a platform call appears in any module. The repository additionally carries the result tables, figures, pre-specifications, control log, and build_manuscript.py, which regenerates the typeset manuscript from source.

What is not released is the agent layer itself. The orchestrator and the 23 sub-agents were runtime configurations of a general-purpose agent on a hosted platform rather than software we wrote, and they cannot be packaged. What is reproducible is the delegation record, comprising the tracks, their waves, their status, duration and token cost, in delegation_record.csv, together with the skills the agents called, which is where the guardrails live. A reader wanting to reproduce the analysis runs the skills directly; a reader wanting to reproduce the system would need to reconstruct the delegation from that record on whatever agent framework they use. We state this explicitly because a paper about agents that ships no agent should say so.

Sixty-eight numbered result tables, 8 figures, the pre-specifications for all three analysis phases, and a log of 186 controls and verdicts, including those that failed, accompany the code.

## Author contributions

**Nicolas Avelar:** Conceptualization, Methodology, Software, Formal analysis,
Investigation, Data curation, Visualization, Writing – original draft, Writing – review
& editing, Project administration, Funding acquisition. **Ruidong Zhang:** Methodology,
Software (original pipeline on which the present analysis was developed), Writing –
review & editing. **Yi Pan:** Conceptualization, Methodology, Writing – review &
editing.

The analysis was executed by a hierarchical multi-agent system (one orchestrating agent
and 23 sub-agent instances dispatched across five waves) operating on six deterministic,
model-free analytical modules, under human design and audit (Sections 2.11, 4.6 and
4.8). N.A. specified the guardrails, audited every reported quantity against its source
table, and is responsible for the final content. All authors reviewed the completed
manuscript and take responsibility for the integrity of the analysis, including the
portion executed by the agent pipeline described in Section 2.11.

## Funding

This research received no external funding. All compute and API costs, including the
$164 consumed by the sub-agent fleet at 2026 pricing (Section 2.11), were funded by the
first author.

## Conflicts of interest

The authors declare no conflicts of interest.

## Ethics statement

This study exclusively re-analyzed de-identified, publicly available human genomic, transcriptomic, and proteomic datasets (Gene Expression Omnibus, TCGA PanCanAtlas, CPTAC/Proteomic Data Commons, GTEx, DepMap, CELLxGENE Census; Section 2.1). No new human subjects data were collected and no additional ethics approval was sought beyond that obtained by the original data-generating studies, consistent with their terms of use.

## Reporting guidelines

As a computational reanalysis rather than a clinical or animal study, no single existing reporting checklist (e.g., STROBE, ARRIVE) applies in full. The manuscript instead follows, and reports adherence to, the pre-specification, multiplicity-correction, matched-null, and negative-control conventions detailed in Section 2.12 and cross-referenced at each point of use (Sections 2.5-2.8, 3.9). Pre-registration of thresholds within phases is stated explicitly wherever a criterion is applied.

## Acknowledgments

The authors thank the patients and donor families whose tissue made the primary datasets
possible, and the consortia that curate and release them: the Gene Expression Omnibus,
TCGA and the PanCanAtlas, GTEx, CPTAC and the Proteomic Data Commons, DepMap, the Human
Protein Atlas, and the maintainers of CellChatDB. This work is entirely a reanalysis of
their data.

## References

1. Meijers WC, Maglione M, Bakker SJL, et al. Heart Failure Stimulates Tumor Growth by Circulating Factors. Circulation 2018;138(7):678-691. doi:10.1161/CIRCULATIONAHA.117.030816.

2. Koelwyn GJ, Newman AAC, Afonso MS, et al. Myocardial infarction accelerates breast cancer via innate immune reprogramming. Nat Med 2020;26(9):1452-1458. doi:10.1038/s41591-020-0964-7.

3. Avraham S, Abu-Sharki S, Shofti R, et al. Early Cardiac Remodeling Promotes Tumor Growth and Metastasis. Circulation 2020;142(7):670-683. doi:10.1161/CIRCULATIONAHA.120.046471.

4. Hasin T, Gerber Y, McNallan SM, et al. Patients with heart failure have an increased risk of incident cancer. J Am Coll Cardiol 2013;62(10):881-6. doi:10.1016/j.jacc.2013.04.088.

5. Banke A, Schou M, Videbaek L, et al. Incidence of cancer in patients with chronic heart failure: a long-term follow-up study. Eur J Heart Fail 2016;18(3):260-6. doi:10.1002/ejhf.472.

6. Roderburg C, Loosen SH, Jahn JK, et al. Heart failure is associated with an increased incidence of cancer diagnoses. ESC Heart Fail 2021;8(5):3628-3633. doi:10.1002/ehf2.13421.

7. Selvaraj S, Bhatt DL, Claggett B, et al. Lack of Association Between Heart Failure and Incident Cancer. J Am Coll Cardiol 2018;71(14):1501-1510. doi:10.1016/j.jacc.2018.01.069.

8. Kwak S, Kwon S, Lee SY, et al. Differential risk of incident cancer in patients with heart failure: A nationwide population-based cohort study. J Cardiol 2021;77(3):231-238. doi:10.1016/j.jjcc.2020.07.026.

9. van den Berg PF, Yousif LI, Koop Y, et al. Framingham risk score associates with incident cancer and heart failure. Eur J Prev Cardiol 2026;33(4):490-497. doi:10.1093/eurjpc/zwaf618.

10. Jaiswal S, Natarajan P, Silver AJ, et al. Clonal Hematopoiesis and Risk of Atherosclerotic Cardiovascular Disease. N Engl J Med 2017;377(2):111-121. doi:10.1056/NEJMoa1701719.

11. Fuster JJ, MacLauchlan S, Zuriaga MA, et al. Clonal hematopoiesis associated with TET2 deficiency accelerates atherosclerosis development in mice. Science 2017;355(6327):842-847. doi:10.1126/science.aag1381.

12. Awwad L, Shofti R, Haas T, et al. Tumor Growth Ameliorates Cardiac Dysfunction. Cells 2023;12(14). doi:10.3390/cells12141853.

13. Sweet ME, Cocciolo A, Slavov D, et al. Transcriptome analysis of human heart failure reveals dysregulated cell adhesion in dilated cardiomyopathy and activated immune pathways in ischemic heart failure. BMC Genomics 2018;19(1):812. doi:10.1186/s12864-018-5213-9. PMID:

14. Liu Y, Morley M, Brandimarto J, et al. RNA-Seq identifies novel myocardial gene expression signatures of heart failure. Genomics 2014;105(2):83-9. doi:10.1016/j.ygeno.2014.12.002. PMID: 25528681.

15. Hannenhalli S, Putt ME, Gilmore JM, et al. Transcriptional genomics associates FOX transcription factors with human heart failure. Circulation 2006;114(12):1269-76. doi:10.1161/CIRCULATIONAHA.106.632430. PMID: 16952980.

16. Kittleson MM, Ye SQ, Irizarry RA, et al. Identification of a gene expression profile that differentiates between ischemic and nonischemic cardiomyopathy. Circulation 2004;110(22):3444-51. doi:10.1161/01.CIR.0000148178.19465.11. PMID: 15557369.

17. Ellis MJ, Gillette M, Carr SA, et al. Connecting genomic alterations to cancer biology with proteomics: the NCI Clinical Proteomic Tumor Analysis Consortium. Cancer Discov 2013;3(10):1108-12. doi:10.1158/2159-8290.CD-13-0219. PMID: 24124232.

18. Lindgren CM, Adams DW, Kimball B, et al. Simplified and Unified Access to Cancer Proteogenomic Data. J Proteome Res 2021;20(4):1902-1910. doi:10.1021/acs.jproteome.0c00919. PMID: 33560848.

19. Edwards NJ, Oberti M, Thangudu RR, et al. The CPTAC Data Portal: A Resource for Cancer Proteomics Research. J Proteome Res 2015;14(6):2707-13. doi:10.1021/pr501254j. PMID: 25873244.

20. Weinstein JN, Collisson EA, Mills GB, et al. The Cancer Genome Atlas Pan-Cancer analysis project. Nat Genet 2013;45(10):1113-20. doi:10.1038/ng.2764. PMID: 24071849.

21. Hoadley KA, Yau C, Hinoue T, et al. Cell-of-Origin Patterns Dominate the Molecular Classification of 10,000 Tumors from 33 Types of Cancer. Cell 2018;173(2):291-304.e6. doi:10.1016/j.cell.2018.03.022. PMID: 29625048.

22. Grossman RL, Heath AP, Ferretti V, et al. Toward a Shared Vision for Cancer Genomic Data. N Engl J Med 2016;375(12):1109-12. doi:10.1056/NEJMp1607591. PMID: 27653561.

23. Li B, Dewey CN. RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC Bioinformatics 2011;12:323. doi:10.1186/1471-2105-12-323. PMID: 21816040.

24. Uhlen M, Fagerberg L, Hallstrom BM, et al. Proteomics. Tissue-based map of the human proteome. Science 2015;347(6220):1260419. doi:10.1126/science.1260419. PMID: 25613900.

25. Uhlen M, Karlsson MJ, Hober A, et al. The human secretome. Sci Signal 2019;12(609). doi:10.1126/scisignal.aaz0274. PMID: 31772123.

26. Jin S, Guerrero-Juarez CF, Zhang L, et al. Inference and analysis of cell-cell communication using CellChat. Nat Commun 2021;12(1):1088. doi:10.1038/s41467-021-21246-9. PMID: 33597522.

27. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 2020;369(6509):1318-1330. doi:10.1126/science.aaz1776. PMID: 32913098.

28. Carter SL, Cibulskis K, Helman E, et al. Absolute quantification of somatic DNA alterations in human cancer. Nat Biotechnol 2012;30(5):413-21. doi:10.1038/nbt.2203. PMID: 22544022.

29. Thorsson V, Gibbs DL, Brown SD, et al. The Immune Landscape of Cancer. Immunity 2018;48(4):812-830.e14. doi:10.1016/j.immuni.2018.03.023. PMID: 29628290.

30. Tsherniak A, Vazquez F, Montgomery PG, et al. Defining a Cancer Dependency Map. Cell 2017;170(3):564-576.e16. doi:10.1016/j.cell.2017.06.010. PMID: 28753430.

31. Dempster JM, Boyle I, Vazquez F, et al. Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. Genome Biol 2021;22(1):343. doi:10.1186/s13059-021-02540-7. PMID: 34930405.

32. Barretina J, Caponigro G, Stransky N, et al. The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 2012;483(7391):603-7. doi:10.1038/nature11003. PMID:

33. Sivakumar P, Thompson JR, Ammar R, et al. RNA sequencing of transplant-stage idiopathic pulmonary fibrosis lung reveals unique pathway regulation. ERJ Open Res 2019;5(3). doi:10.1183/23120541.00117-2019. PMID: 31423451.

34. Suppli MP, Rigbolt KTG, Veidal SS, et al. Hepatic transcriptome signatures in patients with varying degrees of nonalcoholic fatty liver disease compared with healthy normal-weight individuals. Am J Physiol Gastrointest Liver Physiol 2019;316(4):G462-G472. doi:10.1152/ajpgi.00358.2018. PMID: 30653341.

35. Hoang SA, Oseini A, Feaver RE, et al. Gene Expression Predicts Histological Severity and Reveals Distinct Molecular Profiles of Nonalcoholic Fatty Liver Disease. Sci Rep 2019;9(1):12541. doi:10.1038/s41598-019-48746-5. PMID: 31467298.

36. Fan Y, Yi Z, D'agati VD, et al. Comparison of Kidney Transcriptomic Profiles of Early and Advanced Diabetic Nephropathy Reveals Potential New Mechanisms for Disease Progression. Diabetes 2019;68(12):2301-2314. doi:10.2337/db19-0204. PMID: 31578193.

37. Reichart D, Lindberg EL, Maatz H, et al. Pathogenic variants damage cell composition and single cell transcription in cardiomyopathies. Science 2022;377(6606):eabo1984. doi:10.1126/science.abo1984. PMID:

38. Kuppe C, Ramirez Flores RO, Li Z, et al. Spatial multi-omic map of human myocardial infarction. Nature 2022;608(7924):766-777. doi:10.1038/s41586-022-05060-x. PMID: 35948637.

39. Abdulla S, Aevermann B, Assis P, et al. CZ CELLxGENE Discover: a single-cell data platform for scalable exploration, analysis and modeling of aggregated data. Nucleic Acids Res 2025;53(D1):D886-D900. doi:10.1093/nar/gkae1142. PMID: 39607691.

40. Ritchie ME, Phipson B, Wu D, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res 2015;43(7):e47. doi:10.1093/nar/gkv007. PMID: 25605792.

41. Law CW, Chen Y, Shi W, et al. voom: Precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biol 2014;15(2):R29. doi:10.1186/gb-2014-15-2-r29. PMID: 24485249.

42. Smyth GK. Linear models and empirical bayes methods for assessing differential expression in microarray experiments. Stat Appl Genet Mol Biol 2004;3:Article3. doi:10.2202/1544-6115.1027. PMID: 16646809.

43. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. Journal of the Royal Statistical Society Series B 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x. Not indexed in PubMed; record retrieved from the Crossref REST API.

44. Cox DR. Regression Models and Life-Tables. Journal of the Royal Statistical Society Series B 1972;34(2):187-202. doi:10.1111/j.2517-6161.1972.tb00899.x. Not indexed in PubMed; record retrieved from the Crossref REST API.

45. Harris CR, Millman KJ, Van Der Walt SJ, et al. Array programming with NumPy. Nature 2020;585(7825):357-362. doi:10.1038/s41586-020-2649-2. PMID: 32939066.

46. Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 2020;17(3):261-272. doi:10.1038/s41592-019-0686-2. PMID: 32015543.

47. Furusawa H, Cardwell JH, Okamoto T, et al. Chronic Hypersensitivity Pneumonitis, an Interstitial Lung Disease with Distinct Molecular Signatures. Am J Respir Crit Care Med 2020;202(10):1430-1444. doi:10.1164/rccm.202001-0134OC. PMID: 32602730.

48. Akita K, Maurer MS, Tower-Rader A, et al. Comprehensive Proteomics Profiling Identifies Circulating Biomarkers to Distinguish Hypertrophic Cardiomyopathy From Other Cardiomyopathies With Left Ventricular Hypertrophy. Circ Heart Fail 2025;18(1):e012434. doi:10.1161/CIRCHEARTFAILURE.124.012434.

49. Addona TA, Shi X, Keshishian H, et al. A pipeline that integrates the discovery and verification of plasma protein biomarkers reveals candidate markers for cardiovascular disease. Nat Biotechnol 2011;29(7):635-43. doi:10.1038/nbt.1899.

50. Leerink JM, Feijen EAM, Moerland PD, et al. Candidate Plasma Biomarkers to Detect Anthracycline-Related Cardiomyopathy in Childhood Cancer Survivors: A Case Control Study in the Dutch Childhood Cancer Survivor Study. J Am Heart Assoc 2022;11(14):e025935. doi:10.1161/JAHA.121.025935.

51. Türker Duyuler P, Duyuler S, Gök M, et al. Pleiotrophin levels are associated with improved coronary collateral circulation. Coron Artery Dis 2018;29(1):68-73. doi:10.1097/MCA.0000000000000556.

52. Dorn LE, Petrosino JM, Wright P, et al. CTGF/CCN2 is an autocrine regulator of cardiac fibrosis. J Mol Cell Cardiol 2018;121:205-211. doi:10.1016/j.yjmcc.2018.07.130.

53. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional genomics data sets--update. Nucleic Acids Res 2012;41(Database issue):D991-5. doi:10.1093/nar/gks1193. PMID: 23193258.

54. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Res 2002;30(1):207-10. doi:10.1093/nar/30.1.207. PMID: 11752295.

55. Seuthe K, Picard FSR, Winkels H, Pfister R. Cancer Development and Progression in Patients with Heart Failure. Curr Heart Fail Rep 2024;21(6):515-529. doi:10.1007/s11897-024-00680-y. PMID: 39340596.

56. Caller T, Weiss L, Sharon E, et al. Human heart-tumor interaction in ischemic / failing heart is mediated by cardiac mesenchymal stromal cells. Eur Heart J 2025. doi:10.1093/eurheartj/ehaf784.4672.
