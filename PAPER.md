# Reproducibility of automated cardiac and vascular volumetry across chest CT reconstruction kernels: a same-acquisition paired study with executable guardrails

Nicolas Avelar{^a,b}, Ruidong Zhang{^b}, Yi Pan{^b}

{^a} North Oconee High School, Bogart, GA, USA

{^b} School of Computing, University of Georgia, 200 D.W. Brooks Drive, Athens, GA 30602, USA

*Corresponding author:* Nicolas Avelar, School of Computing, University of Georgia, 200 D.W. Brooks Drive, Athens, GA 30602, USA. E-mail: navelar773@gmail.com

## Highlights

- Same-acquisition CT pairs isolate reconstruction effects from anatomy: 193 pairs
- AI heart volume shifts 0.6% across kernels; pulmonary vein 5.0%, an 8-fold gap
- Same-acquisition CT reconstructions differ in z-coverage: 59% of pairs by over 5%
- That false assumption was in our own prespecification, where it could not fail
- Segmentation failed on 4 of 217 series, all archive faults, not model faults

## Abstract

**Background.** Reliability of automated radiological analysis is limited less by model accuracy than by having no mechanism to refuse. A pre-trained segmentation network returns a measurement for any image, including acquisitions outside its validated range, and when analysis is delegated to agents nothing is obliged to notice.

**Methods.** In The Cancer Imaging Archive we identified 193 reconstruction pairs in 59 patients: two CT series from one irradiation event, matched on DICOM frame of reference and acquisition time, differing only in convolution kernel or slice thickness. TotalSegmentator 2.18.0 measured heart, aorta, pulmonary vein and venae cavae across 217 series. Three endpoints and a 10% threshold were fixed in advance as assertions.

**Results.** An eligibility check halted one endpoint at five patients, on metadata. Reconstructions of one acquisition need not share z-coverage: median 15%, and 59% of pairs exceeded 5% despite matched frame of reference; coverage predicts volume difference (rho 0.70), confounding any study pairing series on those keys. In 45 coverage-matched pairs of identical thickness, kernel change altered heart volume by 0.63% (95% CI 0.50-0.97) and pulmonary vein by 4.96% (3.20-5.63); all intervals excluded 10%. Segmentation failed on 4 of 217 series. Radiogenomic association discordance, volumes against tumor RNA-seq in the same patients, did not exceed a permutation null (p=0.33).

**Conclusions.** Whole-heart volumetry on chest CT is reproducible across kernels to about 1% when field of view matches; thin vessels vary several times more. The failed assumption was prose, where it could not halt; measurement caught it.

## Keywords

Computed tomography; Image reconstruction; Automated segmentation; Measurement reproducibility; Radiogenomics; Autonomous agent

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

**NPR1 / NPR2 / NPR3** — Natriuretic peptide receptor 1 / 2 / 3

**PTN** — Pleiotrophin

**AGT / AGTR1** — Angiotensinogen / angiotensin II receptor type 1

**CCN2 / CTGF** — Cellular communication network factor 2 / connective tissue growth factor

**MAGNet** — Myocardial Applied Genomics Network

**IPF** — Idiopathic pulmonary fibrosis

**NAFLD / NASH** — Non-alcoholic fatty liver disease / steatohepatitis

**HR** — Hazard ratio
**GDC** — Genomic Data Commons
**TCIA** — The Cancer Imaging Archive
**kVp** — Peak kilovoltage (X-ray tube potential)
**CI** — Confidence interval

**ρ** — Spearman rank correlation coefficient

## 1. Introduction

Automated organ segmentation is now applied to imaging that was acquired for some other purpose, on scanners and protocols the model never saw during training. A pre-trained network will return a volume for any input it is given. It has no way to signal that the acquisition falls outside the range where its output was validated, and a measurement returned outside that range is indistinguishable, in the output file, from one returned inside it. This is the failure mode that matters as radiological analysis is delegated to autonomous systems: not that a model is wrong, but that nothing in the pipeline is obliged to notice.

The usual remedy is a caveat in a methods section, which binds nobody and which the next study inherits as text rather than as constraint. We test a different one. A methodological requirement can be written as an assertion inside the analysis code, where violating it halts execution rather than producing a footnote. This paper asks what that changes in practice, and it reports the answer from two analyses that share the mechanism and share nothing else.

The imaging analysis uses a design that isolates the model from the patient. Within a single CT study, one irradiation event is frequently reconstructed more than once, at a different convolution kernel or slice thickness. Comparing two such reconstructions holds the patient, the anatomy, the contrast bolus and the moment fixed by construction, so any difference in a measured volume is attributable to the reconstruction and to the model reading it. From four public collections in The Cancer Imaging Archive we assembled 193 such pairs in 59 patients, confirming shared acquisition from DICOM frame of reference and acquisition time rather than inferring it from study membership, and measured five thoracic structures with TotalSegmentator. Three endpoints were fixed in advance: how far volumes move under reconstruction change, whether downstream associations survive that movement, and how often segmentation fails outright.

Two of the three endpoints returned nulls, and both were reached by criteria set before the data could speak. The result we did not anticipate came from our own pre-specification. We had written, as a sentence, that two reconstructions of one acquisition necessarily share z-coverage. They do not: within-pair coverage differs by a median of 15%, and coverage difference predicts apparent volume difference at Spearman rho 0.70 for the aorta. Left unmeasured, that assumption would have converted field-of-view truncation into a reconstruction effect and reported a larger, cleaner-looking number. A sentence cannot fail. The same requirement, written as an assertion over the DICOM header table, would have halted the analysis at the audit stage, before any GPU time was spent. Section 4.8 relates this design to the interpretability, reliability, privacy and ethical concerns that delegated radiological analysis raises.

The second analysis is a transcriptomic screen, included because the mechanism under test should not be specific to imaging. It asks whether the failing human heart releases a signaling protein for which tumors carry a selectively enriched receptor, a question with a specific experimental basis and openly conflicting human evidence{^50}, and one whose answer turns on methodological choices at every join between public datasets. The two cohorts have no biological overlap and are not meant to. Cancer patients scanned for staging cannot speak to heart failure, and we make no claim that they do.

The basis for that second question is as follows. Cardio-oncology has largely studied how cancer treatment damages the heart. The converse direction, whether established heart disease influences the course of cancer, has received attention more recently, and the experimental evidence for it is specific. Transplanting a failing heart into tumor-bearing mice accelerated tumor growth even though the recipient's own circulation was otherwise unchanged, which points to a factor carried in the blood rather than to altered perfusion or oxygen delivery (Meijers 2018{^1}). Myocardial infarction accelerates breast cancer growth in mice through changes in the output of bone-marrow-derived myeloid cells (Koelwyn 2020{^2}; Avraham 2020{^3}).

Human evidence is weaker and openly conflicting. Several cohorts report excess cancer incidence following a heart-failure diagnosis (Hasin 2013{^4}; Banke 2016{^5}; Roderburg 2021{^6}). The largest study with adjudicated cancer outcomes and a 19.9-year median follow-up found no association (adjusted HR 1.05, 95% CI 0.86–1.29; Selvaraj 2018{^7}). In a national cohort, allowing a two-year lag between diagnoses, which removes cancers likely to have been present but undetected when heart failure was diagnosed, attenuated the hazard ratio from 1.64 to 1.09 (Kwak 2021{^8}). Conventional cardiovascular risk factors predict incident cancer about as well as they predict incident heart failure (van den Berg 2026{^9}), and clonal hematopoiesis of indeterminate potential, in which an expanded blood-cell clone carries a somatic mutation without overt disease, is a shared substrate for both (Jaiswal 2017{^10}; Fuster 2017{^11}). Even the direction of effect is unsettled: tumor-bearing mice subjected to pressure overload developed less cardiac hypertrophy and fibrosis than tumor-free controls (Awwad 2023{^12}). Of 107 verified records assembled for this work, 21 report findings inconsistent with a causal link from heart failure to cancer.

We do not attempt to settle causality, because no available dataset links heart-failure status to cancer incidence with molecular resolution. We address a narrower question that follows from the mouse experiments and can be tested with existing human data. If the failing heart releases signaling proteins, a tumor can respond to them only if it expresses the corresponding receptors. Two conditions must therefore hold simultaneously for any specific heart-to-tumor signal: the ligand must be released distinctively by the heart rather than by tissues generally, and its receptor must be present in tumor tissue, ideally at higher levels than in the matched normal tissue from which the tumor arose.

Treating the problem as receptor compatibility rather than causation makes it falsifiable, because it can fail in three separable ways. The ligand may not be cardiac-specific, in which case the heart is not a distinctive source. The receptor may not be transcribed in tumor tissue at all, in which case no signal can be received. Or the receptor may be present but no more abundant than in normal tissue, in which case there is no reason to expect a tumor-selective effect. We tested all three across 30 solid tumor types and report where each candidate axis failed.

The screen was executed by an autonomous agent, and the work involved is the kind now routinely delegated: five public resources in three unit systems, a join at every boundary, and a methodological choice at each join that determines the answer while leaving no trace in the output when it goes wrong. Two of those choices went wrong here, in ways that produced confident and entirely spurious rankings. Both were caught, and the question this paper takes seriously is what to do with a correction once it has been made. Recorded as a sentence in a methods section, it binds nobody, and the next study inherits the text rather than the constraint. We instead packaged the analysis as six reusable skills in which each correction is an assertion that halts execution, and we report both what that changed and what it cost. The second imaging endpoint is radiogenomic in design, pairing image-derived volumes with primary-tumor RNA-seq in the same patients, and it is reported as a null.

## 2. Material and methods

The imaging methods are given first (Section 2.1). Sections 2.2 to 2.10 describe the transcriptomic screen, and Sections 2.11 and 2.12 the reusable modules and the agent system that executed both arms.

All analysis is transcript-level and uses public data. Parameter values given here are the values used, not library defaults.

Each step is implemented as a skill: a directory holding a documentation file and a Python module. The documentation is addressed to an agent deciding whether the procedure fits the task in front of it and carries a structured header naming what the module does and when to invoke it. The module is what executes, and is where the methodological requirements live as assertions. Six were written for this analysis, and each subsection below names the one that performs its step. The unit is deliberately larger than a function and smaller than a pipeline, corresponding to the amount of work that recurs across studies of similar shape. Section 2.11 gives the measured properties of the whole set.

### 2.1 Imaging cohort, reconstruction pairing, and segmentation

Computed tomography series were drawn from four TCGA collections in The Cancer Imaging
Archive (TCGA-LUSC, TCGA-LUAD, TCGA-ESCA, TCGA-SARC). DICOM headers for one representative
slice per series were retrieved through the archive's `getSOPInstanceUIDs` and
`getSingleImage` endpoints, which return a complete header for approximately 0.5 MB per
series and require no volume transfer; 573 series were probed with no fetch errors.

The analysis unit is a *reconstruction pair*: two series belonging to one study that share
`FrameOfReferenceUID` and `AcquisitionTime`, and therefore derive from a single irradiation
event, but differ in `ConvolutionKernel` or `SliceThickness`. Shared frame of reference and
acquisition time were verified by a second header pull for the 234 series involved rather
than inferred from study membership. This yielded **193 confirmed pairs across 217 distinct
series in 59 patients**. Twenty-one distinct reconstruction kernels appear in the cohort
(GE `STANDARD`, `LUNG`, `SOFT`, `BONEPLUS`, `DETAIL`; Siemens `B10f`–`B80f`; Philips `YA`,
`YB`), and slice thickness spans 0.625–5.0 mm.

Segmentation used TotalSegmentator 2.18.0 (nnU-Net 2.8.1, PyTorch 2.10.0+cu128, Tesla T4),
open `total` task, with the six model checkpoints hashed and recorded
(`model_provenance.json`). The five pre-specified structures are **heart, aorta, pulmonary
vein, superior vena cava and inferior vena cava**. The pre-specification originally named
the four cardiac chambers and myocardium; inspection of the distributed package showed
those come from task `heartchambers_highres`, which the package marks `commercial: True`
and gates behind a licence. The structure list was re-specified on the open model before
any segmentation was run, and the amendment is dated in `PRESPECIFICATION_IMAGING.md`.

Volumes were taken from the tool's own `statistics.json` (`-s -sx`) rather than recomputed
downstream. `--stats_include_incomplete` was set: by default the tool reports volume 0.0
for any structure touching the image boundary, and thoracic staging CT routinely clips the
great vessels. Volumes for truncated structures are therefore **partial-field measurements,
not anatomical totals**, and are not reported as absolute anatomy. Fast mode was not used;
it resamples to 3.0 mm, which would attenuate the thickness differences under study.

### 2.2 Data sources

Retrieval is the geo-bulk-de and tcga-pancan skills. Heart-failure transcriptomes were GSE116250 (discovery; left-ventricular RNA-seq, submitter-supplied RPKM{^13}; 37 dilated cardiomyopathy, 13 ischemic cardiomyopathy, 14 non-failing donors, with the primary contrast pooling the 50 failing hearts against the 14 non-failing) and GSE141910 (validation; MAGNet left ventricle, per-sample log-scale normalized values{^14}; 166 dilated cardiomyopathy, 28 hypertrophic cardiomyopathy and 166 non-failing donors, with 6 peripartum cardiomyopathy samples excluded, giving 194 versus 166). Group labels were parsed from the series-matrix characteristics fields rather than hardcoded, since hardcoding is the mechanism by which a mislabeled cohort becomes a silent result.

Cross-platform replication used three further cohorts, all Affymetrix microarray and none used elsewhere: GSE5406 (194 failing, 16 non-failing; GPL96), GSE57345 (177 failing, 136 non-failing; GPL11532) and GSE1869 (31 failing, 6 non-failing; GPL96).{^15,16} Probe-level values were collapsed to genes by maximum mean intensity with a deterministic tie-break. Probes annotated to several genes were dropped rather than assigned to the first listed gene, and log scale was detected from the values rather than assumed.

Protein abundance for the abandoned protein-level test came from eight CPTAC tandem-mass-tag proteome studies in the Proteomic Data Commons{^17,18,19}: ccRCC (PDC000127), ovarian (PDC000110), breast (PDC000120), pancreatic (PDC000270), endometrial (PDC000125), lung adenocarcinoma (PDC000153), head and neck (PDC000221) and glioblastoma (PDC000204). Log ratios were taken against the within-plex reference, with aliquots labeled tumor or adjacent normal from the study biospecimen records and assigned to TMT plexes from the experimental design. The GPL9052 sub-series of GSE57345 is a stub carrying no value rows and 6 samples, and was not analyzed.

Tumor transcriptomes came from the GDC PanCanAtlas EBPlusPlusAdjustPANCAN gene-expression matrix. Values are RSEM normalized_count, batch-corrected and upper-quartile scaled; they are neither TPM nor FPKM{^20,21,22,23} a distinction that becomes consequential in Section 3.6. Aliquot barcodes were parsed to sample type and study, and analysis retained primary solid tumors only (sample-type code 01), excluding LAML, DLBC and THYM, giving 9,538 tumors across 30 types. Adjacent normal tissue (sample-type code 11) was streamed from the same matrix: 710 samples mapped to a tumor type, of which 680 fall in the 15 types retaining at least 10 normals.

Supporting resources were the Human Protein Atlas protein-class and secretome annotations{^24,25} (search_download.php, columns g,eg,pc,secl,scml,rnats,blconcms); CellChatDB human (CellChatDB.human.rda){^26}; GTEx v8 median TPM per tissue across 52 tissues{^27}; the Thorsson immune-landscape feature tables, leukocyte fraction, ABSOLUTE purity and TCGA-CDR{^28} clinical endpoints from GDC{^29}; and DepMap 24Q4 Public (figshare article 27993248). Because depmap.org serves a bot-verification interstitial{^30,31,32} the DepMap release is resolved through the figshare API and the payload is asserted not to be HTML, since storing an error page as data yields plausible but meaningless output. Non-cardiac fibrosis cohorts for Section 3.10 were GSE150910 and GSE134692 (lung, idiopathic pulmonary fibrosis), GSE126848 and GSE130970 (liver, NASH/NAFLD) and GSE142025 (kidney, diabetic nephropathy).{^33,34,35,36} Human heart single-nucleus atlases (Reichart dilated and arrhythmogenic cardiomyopathy; Kuppe myocardial infarction) were obtained via CELLxGENE Census.{^37,38,39}

### 2.3 Differential expression and cross-cohort replication

This section is the geo-bulk-de skill. Each cohort was analyzed separately, on its own scale, with an empirical-Bayes moderated t-test equivalent to limma-trend implemented directly in NumPy and SciPy.{^40,41,42} Genes were retained when expression exceeded 1.0 (RPKM or CPM as appropriate) in at least 20% of the smaller group's samples, and 65 genes carrying a sentinel value of 999999999.0 were detected and dropped from GSE116250. The prior degrees of freedom fitted to 3.04 in discovery and 4.82 in validation, and the variance trend used a 0.4 lowess fraction. Multiplicity was controlled by Benjamini–Hochberg FDR at 0.05.{^43} Ensembl identifiers were version-stripped before intersecting the cohorts, yielding 18,061 shared genes.

A gene entered the replicated signature when it reached FDR < 0.05 in both cohorts with the same direction of effect. Replication was quantified by directional concordance, log fold-change Spearman correlation, and a hypergeometric overlap test. The overlap test is near-saturated here, since 44.6% and 66.2% of shared genes are significant in the two cohorts, so concordance and effect-size correlation are the informative metrics and the hypergeometric p is reported alongside its fold enrichment. As a negative control, the discovery group labels were permuted (seed 42) and the identical pipeline re-run against the real validation result.

The implementation was validated against canonical limma 3.66.0 with the criterion fixed beforehand: a log fold-change Pearson correlation of at least 0.99 and at least 95% overlap of the FDR < 0.05 gene sets, in both cohorts. Agreement is exact to floating-point precision. Pearson and Spearman correlations on log fold change are 1.000000 in both cohorts, with maximum absolute differences of 4×10{^-14} and 8×10{^-14}, and the FDR < 0.05 sets agree at Jaccard 0.9991 and 0.9998. The empirical-Bayes hyperparameters were estimated independently and match (prior degrees of freedom 3.041 against 3.037, and 4.818 against 4.825). Applying the replication rule to the limma output recovers 5,036 genes sharing 5,033 with the signature reported here. Full comparison tables, the least-agreeing variant and the one comparison that could not be made are given in the supplement, and agreement is shown in Supplementary Figure S2.

### 2.4 Secretome annotation and ligand–receptor mapping

These steps are the hpa-secretome and cellchat-lr skills. Up-replicated genes were annotated against HPA protein classes and tiered as core (predicted secreted and plasma-detectable), extended (secreted only) or none. A rule-based regex flag excluded structural extracellular-matrix families, which are secreted but are not plausible long-range endocrine signals: collagens, laminins, fibrillins, fibulins, elastin and microfibril proteins, fibronectin, nidogens, tenascins, structural proteoglycans and cartilage matrix. Exclusions were logged individually rather than aggregated into a count. The HPA response was required to parse as TSV with more than 10,000 rows, which guards against storing an error page as data.

Surviving ligands were mapped to CellChatDB, with analysis restricted to the Secreted Signaling annotation class. Heteromeric receptors were expanded to their constituent subunits, so that a receptor complex is scored on its obligate parts rather than on a single representative gene.

### 2.5 Two-factor compatibility score

This section and Sections 2.6 to 2.8 are the crosstalk-score skill, where the methodological risk concentrates and where 13 of the 141 guardrail checks sit. For each ligand–receptor interaction and each tumor type:

**Receptor capacity.** Subunit expression is aggregated per interaction using a minimum rule, which encodes obligate-subunit semantics: a complex cannot signal above the level of its scarcest component. A geometric-mean rule was run as a sensitivity analysis.

**Expression floor (mandatory).** Aggregated capacity must exceed an absolute floor set at 0.005 × the median GAPDH level across per-tumor-type medians. Expressing the floor relative to a housekeeping gene makes it unit-agnostic and therefore portable across platforms. In this matrix the median GAPDH value is 69,134.5, giving a floor of 345.67 RSEM normalized_count. The validator raises by default and names the offending receptor; interactions failing it in every tumor type are dropped before scoring rather than down-weighted.

**Standardization and weighting.** Floor-passing capacities are z-standardized per receptor across tumor types, then averaged across interactions weighted by ligand availability, defined as the mean replicated log fold change of the ligand across the two heart cohorts.

**Uncertainty and stability.** Bootstrap 95% confidence intervals over 400 resamples (seed 20260727); leave-one-ligand-out and leave-one-pathway-out re-ranking; the geometric-mean complex rule; an unweighted variant; and a double-centered variant removing both row and column means.

### 2.6 Specificity null

The load-bearing validator asks whether a ranking is receptor-specific or merely tracks a global per-tumor expression axis. Each real receptor subunit is replaced by a random gene matched on mean-expression decile, the score is recomputed, and the null ranking is correlated against the observed one, over 300 iterations at seed 20260727. The background pool comprises 200 genes drawn 20 per decile from equal-count deciles of transcriptome-wide median expression across the 9,538 retained tumors, excluding the target panel and its legacy aliases. Four configurations were run: with the floor (primary), without it, from a uniform rather than decile-matched pool, and double-centered. A per-tumor-type empirical p is the fraction of null draws whose score equals or exceeds the observed one.

### 2.7 Ligand and receptor specificity criteria

**Ligand (source) criterion.** A ligand is cardiac-distinctive if the heart ranks first among the 52 GTEx tissues by median TPM and reaches an abundance plausible for a circulating signal.

**Receptor (receiver) criterion.** Because TCGA RSEM and GTEx TPM cannot be placed on a common scale, two independent strategies were used and required to agree. The first is a background-corrected percentile-rank comparison, in which each gene is converted to a within-dataset percentile among the same 1,379 genes shared by both matrices, separately per tumor type and per GTEx tissue, and the percentile delta is placed against the distribution of the same delta for all 1,379 background genes in that pair; enriched means top decile and positive, depleted bottom decile and negative. The second is a genuinely unit-matched log fold change against TCGA adjacent normal from the same source matrix, with BH-FDR within tumor type. No absolute value was subtracted, divided or plotted across the two platforms, and cross-dataset direction claims are ordinal only.

An axis survives only if the ligand criterion and the receptor criterion both hold.

### 2.8 Induction-based source criteria

Because the criterion in Section 2.7 scores the healthy heart, it is structurally unable to detect a ligand induced only in failure, which is a material limitation for a hypothesis about failing myocardium. A corrected definition was pre-specified before any computation, with three criteria, each computed within a single dataset so that no cross-platform absolute comparison is implied. A1 requires replicated induction of at least 1.0 log unit at FDR < 0.05 in both heart cohorts. A2 requires abundance of at least 0.005 of GAPDH inside each heart cohort, the same housekeeping-relative device used for the tumor floor. A3 requires the heart to rank in the upper half of the 52 GTEx tissues rather than first. The receptor side was left unchanged. A pass required at least one surviving axis, with a guard declaring the criterion uninformative above twenty.

### 2.9 Adversarial tests of recovered axes

Kill criteria were fixed before any computation, with the receptor base rate measured before thresholds were chosen so that they could not be set to flatter the result. Retention required all three tracks to survive.

**T1, receiver discrimination.** The ligand-free base rate at which the 81 floor-passing receptors satisfy the receiver criterion in at least one tumor type, declared uninformative at 60% or above; a matched-random null substituting the candidate receptor set with expression-decile-matched sets of equal size, 1,000 permutations at seed 20260729, requiring the observed count to exceed the null's 95th percentile; and pairwise co-expression among the candidate receptors, where a mean Spearman ρ of 0.7 or above would indicate that one observation was being reported several times.

**T2, tissue generality.** The same moderated-t pipeline applied to five non-cardiac fibrosis and organ-failure cohorts, each cohort's normalization kept separate and COL1A1 required as a positive control. The criterion fails if the ligand is significantly induced in at least two cohorts at a median of at least 50% of the cardiac effect. Cell of origin was resolved in two single-nucleus atlases by per-cell-type expression share and by donor-level pseudobulk tests within the dominant cell type.

**T3, literature.** A PubMed search for circulating measurements of the ligand in heart failure. Published evidence contradicting a cardiac plasma source is a kill criterion.

**Conjunction rule.** Where a hypothesis generated three predictions, all three were required to hold in the pre-specified direction. A single significant result with the others opposing is recorded as a failure, and any secondary read conditional on the conjunction was not performed when the conjunction failed. Survival analysis used Cox models{^44} with the receptor z-scored within tumor type and stratified by type, and copy-number loss was tested against 1,000 loci matched on chromosome arm and overlapping-segment length.

### 2.10 Dependency screen

This section is the depmap-lineage skill. Per-lineage CRISPR gene effect was averaged across cell lines within each DepMap OncotreeLineage, then correlated (Spearman, BH-FDR) against the compatibility score of the TCGA types to which that lineage maps. The mapping was passed explicitly rather than inferred, and at least 5 screened lines per lineage were required, retaining 22 of 31 lineages. A gene effect below −0.5 was treated as a dependency, and essential and non-essential genes were carried as parsing controls.

### 2.11 The six skills, measured

**Table 1.** The six analytical skills, their function, and their load-bearing guardrails. Line, function, and guardrail counts are totals across all six modules as measured by static scan (Section 2.11); the distribution across modules is given in the text.

| Skill (module) | Function | Representative guardrails |
|---|---|---|
| geo-bulk-de | GEO retrieval; empirical-Bayes moderated-t differential expression (limma-trend equivalent) | Detects and drops sentinel-value placeholders; detects log scale rather than assuming it; parses group labels from series-matrix fields rather than hardcoding them |
| hpa-secretome | HPA protein-class annotation and secretion tiering (core / extended / none) | Rejects a response that is not a >10,000-row TSV (guards against silently storing an HTML error page as data) |
| cellchat-lr | Ligand-receptor mapping via CellChatDB; receptor-complex subunit expansion | Restricts mapping to the Secreted Signaling annotation class; scores obligate complexes on their scarcest subunit |
| tcga-pancan | Memory-bounded streaming of the PanCanAtlas expression matrix; sample-type parsing | Excludes non-solid-tumor types (LAML, DLBC, THYM); separates primary tumor from adjacent-normal by aliquot barcode |
| crosstalk-score | Two-factor (availability × capacity) compatibility scoring | **Mandatory expression floor** (0.005×GAPDH) before any receptor is scored; **expression-matched null construction** (not a naive random pool); **conjunction rule** that raises rather than records a partial pass. Carries 13 of the 141 total guardrail checks. |
| depmap-lineage | CRISPR dependency screen via DepMap, aggregated per lineage | Resolves the DepMap release via the figshare API and asserts the payload is not HTML; requires an explicit (not inferred) lineage-to-tumor-type mapping and a minimum of 5 screened lines per lineage |
The six modules total 4,858 lines of Python containing 114 functions, with 1,522 lines of accompanying documentation. Their third-party dependencies are NumPy, pandas, SciPy, statsmodels and requests{^45,46} with rdata required for one optional file reader; everything else is standard library. They contain no language-model call, no API key and no platform function, and a scan of all six executable modules for those returns zero in every category. That property is what makes the rest of this paper checkable, since every figure and table below is the output of code no model participated in executing, and the modules return the same numbers whichever agent invoked them or whether one did at all. The claim covers the module rather than the whole skill, since the documentation is addressed to an agent and is not executable, and no result here depends on how an agent reads it.

Across the six there are 141 executable guardrail checks: 78 assertions and 63 explicit raises. The distribution is informative rather than uniform. The scoring module, where the methodological risk concentrates, carries 13; the retrieval module, which handles heterogeneous public data and must reject malformed input in many ways, carries 47. Three are load-bearing for the results below and are named where they act. The expression floor refuses to score a receptor whose abundance falls below an absolute threshold, and raises an error naming the gene (Section 3.6). The null constructor refuses to draw a background set not matched on expression level (Section 3.6). The conjunction rule refuses to record a partial conjunction as a pass, so an axis satisfying two of three criteria cannot be counted as surviving (Section 3.10).

The distinction this draws is between a check a reader must choose to honor and a check a caller cannot avoid. A reader of a methods section can omit the floor without producing anything that looks wrong. A caller of this module cannot, because omitting it produces an exception rather than a number.

Random seeds were fixed as follows: 20260727 for the specificity permutation and the bootstrap, 42 for the differential-expression label permutation, 20260728 for the induction-criteria phase, and 20260729 for the adversarial receptor null. All computation used double precision.

### 2.12 The multi-agent system that ran the analysis

The analysis was executed by a hierarchical multi-agent system rather than by a single agent working serially. Because the call for reliability evidence in this area is specific, we report the topology and its failures as measured quantities rather than as description.

One orchestrating agent held the research question and the artifact store, and independent sub-agents were dispatched to tracks that could proceed without each other's intermediate state. Twenty-three sub-agent frames ran in five waves over 91.4 hours of wall-clock time (Fig. 1). Wave 1 split the four data-acquisition tracks that share no inputs (heart-failure signature, secretome and ligand–receptor mapping, tumor expression axis, literature corpus) together with the dependency screen and normal-tissue comparison. Wave 2 ran five follow-up tests of the negative. Wave 3 ran the two orthogonal-molecule tracks, proteome and platform replication. Wave 4 ran six verification tracks, including two independent audits of the manuscript's own claims. Wave 5 ran the implementation validation, the sensitivity sweeps and the citation expansion.

Each sub-agent received a task statement and a context summary rather than the orchestrator's conversation, and returned results as artifacts plus a structured output object. This isolation is the principal justification for partitioning the work, since a sub-agent cannot inherit an assumption the orchestrator formed earlier, and two of the wave-4 audit tracks were dispatched specifically to check claims the orchestrator had already written down.

![Figure 1](art_9a2d30bd-e8dc-4c0b-995a-bf881ba0d38f)

**Figure 1.** The multi-agent system that executed the analysis. (a) Wall-clock duration of each of the 23 sub-agent tracks, grouped into the five dispatch waves and colored by terminal status. The two longest tracks (Proteome, Platforms) include time parked on user approvals for network access. All four non-completing tracks had already written their results to the artifact store, and every one of those artifacts is cited in this paper. (b) Token consumption across the fleet. Cost is context rather than generation, which is why re-running the analysis from the six skills requires no model at all. Per-track figures are given in delegation_record.csv.

Nineteen of twenty-three tracks reached completion normally. Three terminated in a failed state and one was cancelled after 199 minutes when it exceeded its memory budget on an 8 GiB machine, and was re-dispatched with a streaming implementation. The instructive part is what happened to the work. All four non-completing tracks had already written their results to the artifact store, and every one of those artifacts is cited in this paper. Failure occurred at the reporting boundary rather than in the analysis, because results were persisted as artifacts as they were produced rather than returned only in a final message. A system that returns work only on successful completion would have lost four tracks and roughly 27 hours of compute.

The sub-agent fleet consumed $164 at 2026 API prices, with an input:output token ratio of 106:1, so essentially all of the cost is context rather than generation. Median track duration was 44 minutes. This is the figure that makes the division of Section 4.7 practical rather than rhetorical, since the deterministic modules those agents called cost nothing per invocation, and re-running the entire computational analysis from the six skills requires no model at all.

## 3. Results

The imaging analysis is reported first (Sections 3.1 to 3.3), followed by the transcriptomic screen that tests the same guardrail mechanism on a different modality (Sections 3.4 to 3.10). The two share no biological premise; readers interested only in the imaging result need Sections 3.1 to 3.3 and 4.1.

### 3.1 A pre-specified eligibility check halts the cross-study endpoint

The pre-specified primary endpoint compared chamber volumes between two *separate* studies
of the same patient, requiring acquisition discordance, an interval of 90 days or less, and
concordant contrast status. Applying those criteria to the header table before any
segmentation ran reduced 35 patients with repeat imaging to 31 with acquisition
discordance, 17 within 90 days, and **5 with verifiable contrast concordance**.

The binding constraint is the contrast field. `ContrastBolusAgent` is populated in 21 of 81
primary series; recovering contrast status from every series in each study, not only the
segmentable one, gained a single study. Of 74 within-patient study pairs, only 20 have
contrast status known on both sides, and 36% of those are contrast-discordant, so the
pairs with unknown status cannot be assumed concordant.

At n = 5 the endpoint cannot reach its own 10% kill threshold under any plausible
within-patient variance. The endpoint was withdrawn before any segmentation was run.
This is reported rather than replaced silently: the eligibility criteria were fixed in
advance, they were evaluated on metadata that cost no compute, and they failed.

### 3.2 Reconstruction kernel changes AI-derived volumes by under 5%

A design available in the same data replaces the withdrawn endpoint. Within a single study,
the same acquisition is frequently reconstructed more than once. Comparing those
reconstructions removes contrast, interval and anatomical change by construction: the
patient, the irradiation event and the moment are identical, and only the reconstruction
differs.

Reconstructions of one acquisition do not share reconstructed extent, and this is the
first result the design returns. Within-pair z-coverage differs by a median of 15% and by
up to 158%, and 110 of 186 pairs (59%) differ by more than 5%. The cause is ordinary
radiological practice rather than archive corruption: a chest–abdomen series and a
heart-focused reformat are built from one irradiation event and exported as two series of
the same study, sharing raw data but not extent.

The consequence for quantitative imaging is direct. Volume difference correlates with
coverage difference (Spearman ρ = 0.70 for aorta, 0.68 for inferior vena cava, 0.50 for
heart), so comparing two reconstructions of one event measures field-of-view truncation
together with kernel and attributes the sum to the kernel. **Matching on
`FrameOfReferenceUID` and `AcquisitionTime` does not remove this**: every pair here
satisfies both, which is how the pairs were built, and 59% still differ in coverage by
more than 5%. A shared frame of reference constrains the coordinate system, not the
reconstructed range. Any retrospective study that pairs series by study, by acquisition
time or by frame of reference, and then attributes the residual difference to
reconstruction, inherits this confound. It is invisible in the series description and
requires reading `ImagePositionPatient` across the stack to detect.

Our own pre-specification asserted the opposite, which is why the effect was measured
rather than assumed (Section 4.1). Analysis is therefore restricted to pairs matched on
z-coverage within 5% with slice thickness held identical, leaving **45 pairs in 31
patients**. This filter is post-hoc and declared as such; it was not pre-specified, it was
forced by a pre-specified assumption failing.

With the patient as the unit of inference and the pair as the unit of measurement, both
pre-specified aggregations agree (Table 2; Fig. 2).

Every interval excludes the pre-specified 10% threshold on both aggregations. Whole-heart
volume is robust to the reconstruction kernel; the thin-walled pulmonary vein is roughly
eight times more sensitive than the heart, the expected direction for a structure whose
apparent volume depends on partial-volume effects at its boundary.

What these magnitudes mean for practice depends on the comparison being made. Reported
interscan repeatability for non-gated CT chamber volumetry is of the order of 5 to 10% of the
mean, so a reconstruction-attributable component near 1% for whole heart is small relative to
the measurement's own noise floor and would not be expected to alter a threshold-based reading.
The vascular structures are the caution. At about 5% for pulmonary vein, reconstruction alone
approaches the lower end of published repeatability, which matters for any longitudinal
comparison assembled from archives where kernel is not held fixed, and for multicentre work
where kernel varies by site rather than by protocol. The practical recommendation from these
data is narrow and specific: when automated vessel volumetry is compared across time points or
sites, record `ConvolutionKernel` and `SliceThickness` and either match on them or report them
as covariates, and match reconstructed z-coverage explicitly rather than assuming that series
from one study share it.

Segmentation failed on 4 of 217 series (1.8%), all with error `MISSING_DICOM_FILES`. These
are incomplete instance sequences in the archive, and the failure arises during format
conversion, before the model is invoked. All four are GE
(4/150) against 0/58 Siemens and 0/9 Philips. Six further series required relaxed DICOM
validation and then succeeded. This is an archive-integrity rate, not a model-reliability
rate, and is reported as such.

**Table 2.** Reconstruction reproducibility of five thoracic structures in 45 coverage-matched same-acquisition pairs (31 patients), kernel differing and slice thickness identical. Both pre-specified aggregations are given with 95% confidence intervals; the pre-specified kill threshold was 10%.

| Structure | Median per patient | 95% CI | Mixed-model mean | 95% CI |
|---|---|---|---|---|
| Heart | 0.63% | 0.50-0.97 | 0.85% | 0.64-1.07 |
| Aorta | 0.98% | 0.69-1.17 | 1.04% | 0.80-1.29 |
| Pulmonary vein | 4.96% | 3.20-5.63 | 5.41% | 3.93-6.89 |
| Superior vena cava | 2.12% | 1.90-3.46 | 2.53% | 1.86-3.20 |
| Inferior vena cava | 1.52% | 1.15-2.22 | 2.73% | 1.54-3.91 |

![Figure 2](art_e2c0dcf3-d415-4488-a4ad-3c8a3aeb30cd)

**Figure 2.** Reconstruction reproducibility of AI-derived thoracic volumes. (a) Volume difference within a pair against z-coverage difference, for all 186 segmented pairs; the dashed line marks the 5% coverage-matching criterion applied for the primary analysis. Coverage difference accounts for a substantial share of the apparent reconstruction effect (Spearman rho 0.70 for aorta). (b) Patient-level volume differences in the 45 coverage-matched pairs (31 patients) where the convolution kernel differs and slice thickness is identical; the dashed line is the pre-specified 10% threshold. (c) The two pre-specified aggregations, median per patient and a patient-random-intercept mixed model, with 95% confidence intervals. Every interval excludes 10% under both. Tables: endpoint1_clean_stratum.csv, endpoint1_mixedmodel.csv.

### 3.3 Association stability under reconstruction change: a null result

The remaining pre-specified endpoint asked whether the screen's own receptor panel produces
stable imaging associations. Primary-tumor RNA-seq is available for 58 of the 59 patients in
the reconstruction cohort, and for 31 of the 31 in the coverage-matched stratum that carries
this endpoint. All 100 floor-passing receptors were scored per patient (obligate complexes
scored as the minimum of their subunits), and each receptor–structure Spearman correlation was
computed twice from the same patients and the same acquisitions, differing only in
reconstruction.

Correlations are highly reproducible: median |Δρ| = 0.014 across 500 associations (Fig. 3).

The endpoint was framed in terms of significance verdicts, and that framing does not
survive scrutiny. Of associations reaching nominal p < 0.05 under at least one
reconstruction, a median of **28% (95% range 5–47%)** are significant under one
reconstruction and not the other, integrating over the two arbitrary analyst choices
available, which pair is used for patients contributing several, and which member of a
pair is labeled first. A directional statement is not available: the pair labeling
carries no kernel ordering, and swapping the labels inverts any "proportion lost" figure.

Calibration against a permutation null settles the interpretation. Permuting patient labels
on the expression matrix destroys any real receptor–volume association while preserving the
pair structure, the reconstruction difference, n, and the p-value density near the
threshold. The null discordance is **22% (95% range 8–41%)**, and the observed value does
not exceed it (p = 0.33). Discordance of this magnitude is a property of dichotomizing at
p < 0.05 with n = 31, arising equally for real and spurious associations, and is not a
specific instability of the segmentation output. Separately, **no association survives
Benjamini–Hochberg correction across the 500 tests** (minimum q = 0.98), so there are no
findings here whose stability could be at issue.

Endpoint 2 therefore reports a null. NPR1, the receptor carried forward from the
transcriptomic screen, shows no sign change and no threshold crossing on any structure
(maximum |Δρ| = 0.046); its heart association (ρ ≈ 0.32, p ≈ 0.08) is non-significant under
both reconstructions.

![Figure 3](art_d1c73ebc-ef4d-4e5a-91d0-3a0c9d1f659e)

**Figure 3.** Association stability under reconstruction change is a null result. (a) Spearman correlations between receptor score and structure volume, computed twice from the same patients and the same acquisitions and differing only in reconstruction (500 associations; median absolute change 0.014). (b) Percentage of nominally significant associations discordant between reconstructions, integrated over the two arbitrary analyst choices (which pair is used for patients contributing several, and which member is labeled first) across 200 realizations, against a permutation null in which patient labels on the expression matrix are shuffled (60 realizations). Observed 28% does not exceed the null 22% (p = 0.33). (c) Distribution of nominal p values across all 500 tests; none survives Benjamini-Hochberg correction (minimum q = 0.98). Tables: endpoint2_associations.csv, endpoint2_verdict.json.

### 3.4 A replicated human heart-failure signature

Differential expression in GSE116250 (50 failing hearts, of which 37 dilated and 13 ischemic cardiomyopathy, against 14 non-failing) recovered the expected cardiac stress programme: NPPA logFC +3.91 and NPPB +3.97, with MYH6 (−0.61) and ATP2A2 (−0.44) down, as expected in failing myocardium. GSE116250 distributes RPKM rather than counts, so a limma-trend-equivalent moderated t on log(RPKM+1) was used, and 65 genes carrying a sentinel value of 999999999 were detected and removed before analysis.

Validation in the larger MAGNet cohort (GSE141910; 194 cardiomyopathy against 166 non-failing donors) replicated the signature. Of 18,061 shared genes, 5,717 were significant in both cohorts and 5,036 shared direction, a directional concordance of 0.881 that rises monotonically with discovery effect size, from 85% at |logFC| < 0.5 to 95% above 1.0. Effect sizes correlated at Spearman ρ=0.541 across all shared genes and ρ=0.759 among those significant in both.

Two cautions matter. The hypergeometric overlap test returns p=2.1×10{^-3} but only 1.07× enrichment: with 8,052 and 11,949 significant genes respectively, the test is near-saturated and is not evidence of replication despite the p-value, so directional concordance and effect-size correlation are the informative metrics. Separately, a label-permutation control run through the identical pipeline returned 0 significant genes and ρ=−0.134 against the real validation result, confirming that the pipeline does not manufacture concordance (Fig. 4).

![Figure 4](art_32f64f02-cc63-4a6a-865e-9e0c0eb0c07b)

**Figure 4.** Replication of the heart-failure expression signature. (a) Discovery volcano, GSE116250 (50 failing, 14 non-failing left ventricle); empirical-Bayes moderated t, Benjamini-Hochberg FDR. (b) Cross-cohort effect-size concordance among the 5,717 genes significant in both cohorts (Spearman rho=0.759; rho=0.541 across all 18,061 shared genes); red points disagree in direction. (c) Directional concordance rises with discovery effect size, from 85% below |log2FC| 0.5 to 95% above 1.0. Validation cohort GSE141910 (194 cardiomyopathy, 166 non-failing). Label permutation returned 0 significant genes.

Pooling cardiomyopathy etiologies does not mask a secreted signal. Re-deriving differential expression within each etiology separately against the same controls, 12 of the 42 upregulated secreted ligands are induced in one etiology and not the other (POSTN, INHBA, HBEGF, CCL3, NTF3, LGALS9, CCL5, CCL4, BMP4, WNT10B, CX3CL1 and MSTN), but none is cardiac-enriched at GTEx rank 1 and only five are abundant, so none would have passed the criteria NPPA and NPPB passed. Label-permutation controls returned no significant genes in all four stratified contrasts. The ischemic and hypertrophic strata are small, at 13 and 28 samples, so a null there is weak evidence of absence rather than evidence of a null.

The signature was also re-derived in three Affymetrix cohorts not used for discovery (GSE5406, GSE57345, GSE1869) against thresholds fixed in advance. The two adequately powered cohorts replicate it and the third, at 31 against 6 samples, does not; label permutation returned no significant genes in any contrast. Sharing a sequencing platform inflates apparent agreement by 0.17, so the replication of Section 3.4 is an upper bound. Full per-cohort values are in the repository.
### 3.5 From signature to candidate crosstalk network

Of 3,084 replicated HF-upregulated genes, 295 are annotated as secreted by the Human Protein Atlas and 70 are CellChatDB ligands of any annotation class, 50 of them within Secreted Signaling. Restricting to Secreted Signaling, the only annotation class mechanistically plausible at a distance, gave 188 interactions from 51 ligands across 38 pathways, with receptor complexes expanded to 103 constituent subunits. Ligand availability was computed as the mean replicated effect size scaled by secretion tier; the highest were NPPA (3.05), NPPB (3.03), PENK (2.85) and WNT9A (2.32). None of these quantities decided anything on its own: the criteria that act on them were fixed before this step ran, which is what makes the attrition reported in Sections 3.9 and 3.10 a test rather than a selection.

### 3.6 The absolute expression floor as a determinant of the result

TCGA PanCanAtlas values are RSEM normalized counts rather than TPM. The median GAPDH value across the 30 per-tumor-type medians is 69,134.5, so any floor expressed in transcript-per-million terms is meaningless in this matrix, and a familiar TPM threshold imported from elsewhere would return a number that carried no information. The floor was therefore defined relative to a housekeeping gene: 0.005×GAPDH = 345.67 in these units.

Twenty-eight of 188 interactions fall below it in every one of the 30 tumor types. These include all four proenkephalin-to-opioid-receptor axes (OPRD1, OPRK1, OPRL1, OPRM1), where OPRD1 reaches at most 4×10{^-4} of GAPDH, OPRM1 5×10{^-4} and OPRK1 1.4×10{^-3}. PENK is the third most available HF ligand, so an unfiltered analysis would have nominated the opioid axis prominently, whereas the entire pathway is non-viable on receptor-expression grounds.

The floor is not hygiene; it is what makes the score interpretable. Under an expression-matched random receptor null (300 permutations, seed 20260727), the mean ρ against the observed ranking and the fraction of draws exceeding ρ=0.5 were as follows.

No floor applied: ρ = 0.363, with 21.5% of draws exceeding 0.5.

Floor applied (primary): ρ = 0.107, with 1.7%.

Floor applied, double-centered: ρ = 0.034, with 1.0%.

Unmatched (naive) random pool: ρ = 0.323, with 17.0%.

The global expression-axis confound affecting co-expression rankings is therefore attributable in large part to the scoring of untranscribed receptors; their removal renders the ranking receptor-specific. The null pool must also be expression-matched, since an unmatched pool yields ρ=0.323 and would incorrectly reject a receptor-specific score (Fig. 5).

![Figure 5](art_62ff7487-a8e3-476a-83df-92f5b5f22f80)

**Figure 5.** The expression floor and the two-factor decomposition. (a) Floor audit across all 188 candidate interactions; 28 (red) never clear the floor of 0.005xGAPDH (345.67 RSEM normalized counts) in any of the 30 tumor types, including all four PENK-opioid axes at 4x10-5 to 3x10-3 of GAPDH. (b) Two-factor decomposition of the compatibility score into ligand availability and receiver capacity. Availability is a property of the source and is constant across receivers, so only capacity can explain differences between tumor types.

The floor value is not tuned, but neither does it sit on a broad plateau. A pre-registered plateau test asked whether, across a fourfold window from 0.0025 to 0.01 of GAPDH, the null correlation stayed below 0.20, the count of types beating the null stayed within one of four, and pass-set membership changed by at most one member. It fails two of the three. The pass set is identical at 0.0025 and 0.005 (UCS, BRCA, KIRC, PAAD), so the published value is not perched on a knife edge. However, the null correlation falls monotonically across the whole sweep with no flat region (ρ = 0.295, 0.201, 0.113, −0.021, −0.161 and −0.184 at fractions 0.001, 0.0025, 0.005, 0.01, 0.025 and 0.05), and a twofold move to 0.01 changes the pass set to BRCA, KIRC, OV, THCA and UCS, gaining two types and losing one. The defensible reading is that a floor is necessary and any value in this range makes the score receptor-specific, while the identity of the four passing types is sensitive to the exact fraction. Combined with the background-pool result below, the count is better reported as four of 30 with a range of three to seven across defensible choices.

Both the permutation resolution and the sensitivity of the decile-matched null to bin width and pool composition are reported in the reproduction package.

### 3.7 Pan-cancer compatibility ranking and its statistical resolution

Scoring the 160 floor-passing interactions ranked uterine carcinosarcoma first (0.679, bootstrap 95% CI 0.400–0.914, rank SD 0.55 over 400 resamples), followed by ovarian serous cystadenocarcinoma (0.432), clear-cell renal carcinoma (0.377), thyroid (0.285), prostate (0.272) and pancreatic adenocarcinoma (0.244). Only 4 of 30 types beat the matched-random null at p < 0.05: UCS (p=0.003), KIRC (p=0.013), BRCA (p=0.010) and PAAD (p=0.033). Ovarian carcinoma, ranked second by score, is not among them (p=0.060). Rank order and receptor-specificity are therefore not the same thing, and a high score alone is insufficient grounds to nominate a tumor type.

Stability controls pass. Leave-one-ligand-out preserves the ranking at minimum ρ=0.907 (worst case WNT9A) and leave-one-pathway-out at ρ=0.842 (WNT, the largest single contributor at 40 of the 160 floor-passing interactions, or 80 rows once receptor complexes are expanded to subunits). Switching the receptor-complex rule from minimum-of-subunits to geometric mean gives ρ=0.908, and removing availability weighting altogether gives ρ=0.922, so no single ligand, pathway or parameter choice drives the result.

One control does not pass. Double-centering, which removes additive row and column effects, reorders the ranking substantially, at ρ=0.318 against the primary, promoting adrenocortical carcinoma to first and demoting UCS to sixth. Ovarian carcinoma is the only type in the top two under both variants. The ranking is correction-sensitive, and mid-table positions should not be interpreted (Fig. 6).

![Figure 6](art_b2c7e5cf-eb22-48ce-a162-0ceeb6fc09e1)

**Figure 6.** Pan-cancer compatibility ranking against an expression-matched null. (a) Compatibility score for each of 30 tumor types with 95% bootstrap confidence intervals (400 resamples, seed 20260727); the 4 types shown in blue exceed the expression-decile-matched random null at p < 0.05 (UCS 0.003, BRCA 0.010, KIRC 0.013, PAAD 0.033). Ovarian carcinoma ranks second by score but does not pass (p=0.060). Individual positions are not interpretable: the ranking reorders substantially under double-centring (ρ=0.318 against the primary ranking). (b) Mean null correlation under three background constructions, showing that the absolute expression floor rather than the expression matching is what renders the score receptor-specific.

### 3.8 The exploitation component of the hypothesis is unsupported

Our hypothesis held that susceptibility depends on the ability to receive and then exploit cardiac signals. Availability is constant across tumor types by construction, so only receiver capacity can explain differences. Testing capacity against 29 Thorsson immune-landscape features, substituting for TIMER2, which was unreachable, only one survives FDR correction: B-cell receptor score (ρ=0.571, FDR=0.028). Leukocyte fraction, stromal fraction, TGF-β score, proliferation and wound-healing signatures are all null. Receptor compatibility and microenvironment state are therefore largely independent axes, and we found no evidence for the exploitation component as formulated.

Receptor expression also failed to correspond to CRISPR lineage dependency in DepMap 24Q4 (0 of 81 floor-passing receptors significant after FDR correction; NPR1 pan-lineage gene effect -0.163), which is reported in full in the reproduction package.

The other two natriuretic peptide receptors were tested and do not change the result. Within CellChatDB, NPR2 is paired only with NPPC and NPR3 only with OSTN, neither of which is induced in the failing ventricle, so NPR1 is the only receptor reachable from NPPA or NPPB in the curated set. That is a property of the database rather than of receptor pharmacology, since NPR3 binds both peptides as a clearance receptor in vivo; a curation gap of this kind cannot be detected from inside the screen. Independently of curation, all three receptors are depleted rather than enriched in tumor across the unit-matched types, so the exploitation component fails on the data as well as on the database.
### 3.9 Specificity checks under a healthy-heart ligand filter

Two questions determine whether any axis is meaningful: whether the ligand is distinctively cardiac, and whether the receptor is tumor-enriched.

Across 52 GTEx tissues, only 2 of 42 floor-passing HF-upregulated secreted ligands are both cardiac-enriched and abundant: NPPA (heart rank 1, 35.8 TPM in left ventricle) and NPPB (rank 1, 26.8 TPM). GDF6 is relatively cardiac-enriched but at 0.18 TPM is too scarce to be a credible circulating signal. The other 39 peak elsewhere, with GDF15 in kidney cortex, CCL5 in whole blood, WNT9A in sigmoid colon (rank 29 of 52) and MDK in ovary. The failing heart is a distinctive source only for the natriuretic peptides.

The negative depends on the specificity metric rather than on the rank depth. Because heart ranking first is one of several defensible definitions, the criterion was re-scored under six of them, with GTEx v8 medians re-fetched and the 52-tissue set reproduced exactly. The three rank-based definitions are interchangeable: heart rank 1, top 3 and top 5 all give NPPA and NPPB as the only ligands that are both specific and abundant, and all give zero surviving axes. A fold-enrichment criterion is not interchangeable with them. At a fivefold threshold AGT joins as a third specific and abundant ligand (heart median 138.2 TPM against an all-tissue median 19.1, fold 7.2), though it fails every rank cutoff at best rank 11, and the axis count stays at zero because the AGT receptor AGTR1 is enriched in no tumor type. The two permissive definitions break the result: requiring only that heart sit above the tissue median admits 19 further ligands, 14 of them abundant, and returns 29 axes with a tumor-enriched receptor, while a twofold criterion returns 14. The negative is therefore a claim about strictly cardiac-restricted ligands. Under a cardiac-expressed-above-the-tissue-median reading the screen returns candidates, and which reading is correct is a question about what an endocrine source must look like rather than a question the data settle.

For receptors, comparison used two independent strategies, because TCGA (RSEM normalized counts) and GTEx (TPM) cannot be placed on a common scale: a background-corrected percentile-rank comparison against GTEx, and a genuinely unit-matched log fold change against TCGA adjacent-normal samples streamed from the same source matrix. Of those, 710 mapped to a tumor type, of which 680 fall in the 15 types retaining at least 10 normals and were analyzed. The two agree at ρ=0.545 (p=8.9×10{^-94}, n=1,200 pairs), validating the rank proxy without equating units. Of 78 receptors compared, 11 are enriched under both methods, 25 are comparable to normal tissue, and 15 are depleted.

The decisive result concerns NPR1. The receptor for the only two genuinely cardiac-enriched ligands is significantly depleted in tumor relative to adjacent normal in all 15 unit-matched tumor types (median background-corrected logFC −1.91, all FDR < 0.05), and independently fails the expression floor in 24 of 30 types, passing only in KIRC, KIRP, MESO, OV, PCPG and SARC. Applying both criteria strictly, requiring the ligand to be cardiac-enriched and abundant and the receptor to be tumor-enriched and floor-passing in a top-6 tumor type, zero of the 146 ligand–receptor-subunit axes derived from the 160 floor-passing interactions survive.

The pattern is therefore inverted at both ends. The axis with a distinctive cardiac source has the least tumor-enriched receptor, and the axes with tumor-enriched receptors are driven by ligands the heart does not distinctively secrete (Fig. 7).

![Figure 7](art_caffe368-443d-42c9-8c25-5ff4d6d020cc)

**Figure 7.** Ligand cardiac specificity and receptor tumor enrichment. (a) Under a rank-based definition (heart rank 1 of 52 GTEx tissues with left-ventricle median at least 5 TPM) only NPPA and NPPB are cardiac-enriched among the 42 HF-upregulated secreted ligands; the criterion is definition-dependent, and a fold-enrichment definition also admits AGT (Section 3.9). (b) NPR1, the sole floor-passing receptor for those two ligands, is depleted in tumor relative to unit-matched adjacent normal in all 15 types with at least 10 adjacent normals (median log2FC -1.91, all FDR < 0.05, 680 normal samples).

### 3.10 Correction of the specificity filter and adversarial testing of the recovered axes

Both the corrected specificity criteria and the three kill criteria applied to the recovered axes were fixed before any computation in this phase, with the receptor base rate measured before thresholds were chosen so that they could not be set to flatter the result (Sections 2.8 and 2.9).

Section 3.9 reported zero surviving axes under a ligand filter that scored cardiac distinctiveness in GTEx left ventricle, which is non-failing myocardium. That filter is wrong for a hypothesis about the failing heart, and correcting it changes the count: re-scoring cardiac source by failure induction, with the receptor criteria unchanged, leaves 7 of 146 axes surviving instead of 0, all with pleiotrophin (PTN) as ligand.

That is a methodological finding about the filter rather than a biological finding about PTN, and the distinction is the substance of this section. PTN entered only because criterion A3 was relaxed from heart ranking first of 52 GTEx tissues to heart ranking in the upper half, and PTN sits at rank 23. A hit that appears when a criterion is loosened is precisely what loosening a criterion produces, so three kill criteria were fixed before any of the tests below ran, with PTN retained only if it survived all three. It failed all three.

**Table 3.** Summary of the three pre-specified kill criteria applied to the seven pleiotrophin (PTN) axes recovered after correcting the ligand-specificity filter (Sections 2.9 and 3.10). All three criteria, and the receptor base rate against which Track 1 was judged, were fixed before any of the tests below were run. Retention required all three to survive.

| Track | Test | Pre-specified threshold | Observed result | Verdict |
|---|---|---|---|---|
| T1 – Receiver discrimination | Ligand-free base rate at which candidate receptors satisfy the receiver criterion; matched-random null (1,000 permutations); pairwise receptor co-expression | Base rate ≥60% declared uninformative; observed count must exceed the null's 95th percentile | 74.1% of all floor-passing receptors pass with no ligand involved (45-77% across aggregation rules); the 7 PTN receptors sit at the 31st percentile of a matched-random null (one-sided p = 0.688) | **Fail** |
| T2 – Tissue generality | Moderated-t differential expression in 5 non-cardiac fibrosis/organ-failure cohorts; cell-of-origin in 2 single-nucleus cardiac atlases | Fails if significantly induced in ≥2 cohorts at a median ≥50% of the cardiac effect | Significant in 2 of 5 cohorts at 50.5% of the cardiac effect (triggers by 1.0%); direction reverses in the largest lung cohort; PTN is 62-71% fibroblast-derived with no significant per-cell induction in either atlas | **Fail** |
| T3 – Literature | PubMed and full-text search for circulating PTN measured in a heart-failure cohort against controls | Published evidence against a cardiac plasma source is a kill criterion | No heart-failure cohort measurement exists in 609 verified records; the one study designed to detect cardiac PTN release (coronary-sinus sampling) explicitly rejected it as non-specific to myocardial injury | **Fail** |
| **Conjunction rule** | All three tracks required to pass | — | 0 of 3 tracks pass | **All 7 axes withdrawn** |
#### Track 1: receiver discrimination

Before asking whether the PTN receptors are special, we asked how special any receptor is. Of the 81 floor-passing receptors, 60 (74.1%) are enriched versus matched adjacent normal in at least one tumor type with no ligand involved. That figure is rule-dependent: scoring receptor complexes whole rather than per subunit gives 46 of 103 (44.7%) under the obligate-subunit rule and 79 of 103 (76.7%) under a geometric mean, so the base rate should be read as 45 to 77%. The conclusion it was raised to support, that most receptors pass with no ligand involved and that a six or seven axis count is therefore uninformative on its own, holds under every rule and is strongest under the geometric mean. The pre-specified threshold for declaring an axis count uninformative was 60%.

Seven receptors pass in more tumor types than any PTN receptor does: PLXNA3 in 13, PLXNA1 in 11, TNFRSF12A in 10, ERBB3 in 9, and ITGA6, F2R and ITGB4 in 8 each. The PTN receptors average 3.86 types against an all-receptor average of 2.52, placing them above the mean but at the 76th percentile of a distribution most members of which already pass. Exceeding the average within a population in which nearly every member passes is not evidence of specificity.

The formal test is worse than merely unremarkable. Substituting the seven PTN receptors with expression-decile-matched random receptor sets of the same size (1,000 permutations, seed 20260729) gives a null median of 30 enriched pairs against 27 observed. The PTN set sits at the 31st percentile of the null, with empirical one-sided p = 0.688, and three null configurations agree. The PTN receptor set performs slightly worse than random receptor sets matched for expression level.

One check passed. The seven receptors are not a single co-expressed family counted seven times, with mean pairwise Spearman ρ of −0.091 across the 30 tumor types and a maximum of 0.692. The axes are genuinely seven independent receptors, and random draws beat them.

#### Track 2: tissue generality

PTN was tested in five non-cardiac organ-failure and fibrosis cohorts, each run through the same moderated-t pipeline with normalization kept separate. Induction is not cardiac-restricted: PTN rises in non-cardiac fibrosis, and the direction is inconsistent across the largest cohorts in a way the pre-specified threshold does not capture. Cell of origin is more decisive. In two independent human heart single-nucleus atlases PTN is fibroblast-derived rather than cardiomyocyte-derived, and per-fibroblast PTN does not rise in disease in either atlas, so the bulk induction reflects a shift in fibroblast fraction rather than increased output per cell. Track 2 fails.

#### Track 3: literature evidence for a circulating source

A structured literature search (52 PubMed queries over a 991-record pleiotrophin corpus, plus a full-text scan of 478 cardiac plasma-proteomics papers) found no report of elevated circulating PTN in heart failure. The two supporting records fall outside heart failure. The secreted role of PTN in cancer is by contrast well established, so the half of the axis that is supported is the tumor half. Track 3 fails.

#### Verdict

The seven PTN axes are withdrawn as a positive finding. All three pre-specified tracks failed. The receptor test is passed by three-quarters of candidate receptors and the PTN set underperforms matched-random draws; PTN is a fibroblast injury transcript induced in lung and kidney fibrosis at half the cardiac effect, with no per-cell cardiac induction; and the plasma measurement on which the endocrine premise depends has never been made, with the nearest published evidence contradicting it.

The conclusion of the paper is therefore unchanged from Section 3.9, but it is now a stronger claim than it was. The original negative could have been an artifact of scoring cardiac specificity in healthy tissue. That objection has been tested directly: correcting the filter does recover axes, and the recovered axes do not survive. A negative that has been attacked on its most obvious weakness and held is worth more than one that has not (Fig. 8).

![Figure 8](art_df78cf23-161c-4d79-89c9-984651f1e760)

**Figure 8.** Cardiac specificity re-scored against failure induction. (a) Sequential attrition across the three pre-specified criteria; 3 of 51 ligands clear all three. The third criterion is A3, requiring the heart to rank in the upper half of the 52 GTEx tissues. (b) Combinations of criteria failed across the same 51 ligands; the 3 ligands failing none are NPPA, NPPB and PTN. (c) Induction magnitude, taken as the minimum log2 fold change across the two cohorts, against baseline cardiac rank; the four ligands clearing both induction and abundance are shown in red, of which MDK is excluded at rank 36 while PTN is retained at rank 23. (d) The surviving pleiotrophin axes scored at complex level under the minimum-subunit rule, by number of tumor types in which the receptor is enriched relative to matched adjacent normal. Ligand criteria are computed within single datasets only; the receptor side is unchanged from Section 3.9. n = 51 ligands, 146 axes.

## 4. Discussion

This paper reports one methodological result and tests it twice. The imaging arm is the test a radiology reader can check directly: an off-the-shelf segmentation model applied to scans acquired for cancer staging, with the reproducibility of its cardiac and vascular volumes measured against reconstruction settings that vary for reasons unrelated to the heart. The transcriptomic arm is the same mechanism applied where we could not see the answer in advance, and it is included because a requirement that halts on one dataset proves little. What both arms test is whether a methodological requirement written as code behaves differently from the same requirement written as prose. Sections 4.1 and 4.8 are the ones that bear on imaging practice; Sections 4.2 to 4.4 concern the biological question the second arm was applied to and can be read or skipped independently.

### 4.1 What the imaging arm establishes

The arm supports two measurement claims and two methodological ones.

The first measurement claim is narrow and quantified: for a widely used open segmentation
model, whole-heart volume derived from thoracic CT is reproducible to about 1% across
reconstruction kernels when field of view is matched, while thin vascular structures vary
several times more. Anyone using automated chamber or vessel volumetry across
heterogeneously reconstructed archives now has an effect size for the reconstruction
component of their measurement error, on this model and this population.

The second is a property of the archives themselves and applies to any retrospective
quantitative imaging study, not only to those using a segmentation model. Reconstructions
of a single acquisition differ in z-coverage by a median of 15%, and 59% of pairs differ by
more than 5%, even though all of them share `FrameOfReferenceUID` and `AcquisitionTime`.
Coverage difference correlates with apparent volume difference at ρ up to 0.70. Pairing
series by study, acquisition time or frame of reference is therefore not sufficient to
isolate a reconstruction effect: the natural matching keys constrain the coordinate system
and the irradiation event, but not the exported range, and the resulting truncation is
attributed to whichever variable the study set out to test. The check is inexpensive —
compare `ImagePositionPatient` extent across the two stacks before comparing anything
derived from them — and it is not something a series description reveals. We report it
because a study that omits it will not see it fail; it will see a plausible effect size for
the wrong variable.

The first methodological claim concerns what prose cannot enforce, and the coverage result
is where we learned it. Our own pre-specification asserted that same-acquisition pairs
share z-coverage. It was written in advance, in good faith, and it was false; an unguarded
analysis would have attributed field-of-view truncation to the reconstruction kernel.
Measurement caught it; the pre-specification did not, because an assumption stated as a
sentence cannot fail. Had the same requirement been written as an assertion over the header
table, it would have halted at the audit stage rather than after segmentation. This is the
paper's thesis applied to the paper's own authors, and it is the strongest evidence in it
that the distinction between executable and asserted requirements is not stylistic.

The second concerns thresholds. Endpoint 2 was designed to ask whether reconstruction
choice destabilises radiogenomic associations, and the honest answer is that at n = 31 the
question cannot be separated from the fragility of the p < 0.05 cut itself: the observed
discordance matches a permutation null in which no association is real. A study reporting
only dichotomized outcomes on a single reconstruction would have reported associations that
another equally valid reconstruction of the same scan does not support, but so would a
study perturbed in any other way by a comparable amount.

Scope is bounded and does not improve with more analysis. Fifty-nine patients, one
pre-trained model, one archive, and cancer patients scanned for cancer staging rather than
patients with heart failure. The imaging population and the transcriptomic ligand source
(explanted failing hearts, GSE116250 and GSE141910) have no biological overlap, so the arm
tests the durability of the screen's output under measurement variation and does not bear
on the heart-failure hypothesis. Any sentence implying otherwise is a defect.

The arm also produced the project's third pre-specified withdrawal, after the seven
pleiotrophin axes and the cross-study imaging endpoint. Two of the three were withdrawn
before compute was spent, on metadata alone.

### 4.2 Location of the failure within the two-factor model
The two-factor framing was chosen because it can fail in identifiable places rather than merely returning a null, and it failed in two of them. The receiver side is not the constraint: receptor transcripts for HF-upregulated ligands are present across tumor types, and a compatibility ranking built on them is stable under resampling and under leaving out any single ligand or pathway. The constraint is specificity at both ends. Only the natriuretic peptides are distinctively and abundantly cardiac among the candidate ligands, and their sole receptor is the one consistently depleted in tumor relative to matched adjacent normal. An axis needs both ends, and no axis has both. The transferable observation is not about this pathway but about where the framing put the failure: because the two factors were declared separately in advance, the null could be localized to one of them instead of reported as an undifferentiated negative. A screen that scores a single composite cannot do that, and the same applies to an imaging pipeline that reports one number where two conditions were required.
### 4.3 Generalizable methodological requirements and their implementation

Three findings from Sections 3.6 and 3.8 are not about the heart. Each is a property of compatibility scoring in general, and each is now an assertion rather than a recommendation.

The first is that any score standardizing expression before checking that the gene is expressed at all inherits an artifact. The confound usually attributed to a global expression axis is substantially produced by scoring untranscribed receptors, since standardizing a gene sitting at instrument noise still yields large values, and those values are what a random draw reproduces.

The second is that the null used to detect this must itself be matched on expression. The naive pool and the floor address different failures, and neither substitutes for the other.

The third is a base rate. A screen that relaxes any criterion should report the rate at which the criterion it did not relax is passed, because without that number an axis count cannot be read as evidence about a particular ligand.

The sequence in Section 3.10, comprising a negative result, a filter identified as capable of having caused it, a correction that recovers hits, and a base-rate check that withdraws them, is generic enough that others will repeat it. This is the rationale for encoding the three requirements at the point of enforcement. Stated as prose here, they would be inherited by the next study as text, and a study that omits the floor would produce a result indistinguishable in form from one that honors it. Stated as assertions in the scoring skill, omitting the floor produces an exception naming the offending gene.

The claim is not that assertions are better than prose. Prose can express nuance that code cannot, and an assertion encodes a threshold rather than the reasoning that justified it, which may be wrong for another data type. The claim is that the two fail differently, and that for requirements of this kind, where the wrong answer resembles a right answer, the failure mode of the assertion is the one to prefer. An untranscribed receptor scored as though it were expressed does not announce itself; a missing exception does.

### 4.4 Measurements required to resolve the question
The measurement the endocrine premise requires has never been made for the recovered ligand: circulating concentration in a heart-failure cohort against controls. Converting tissue mRNA into a circulating contribution requires translation efficiency, secretion rate, clearance and volume of distribution, none of which are available for these ligands, so that calculation was excluded rather than attempted with fabricated parameters. The decisive experiments are measurement rather than analysis.
### 4.5 Limitations

*Imaging arm.*

**Imaging cohort and model.** Fifty-nine patients, one pre-trained model, one archive, and cancer patients scanned for staging rather than patients with heart failure. The reconstruction result is a property of TotalSegmentator 2.18.0 on thoracic CT from four TCGA collections, and does not transfer to other models, other body regions, or gated cardiac acquisitions without retesting.

**Structures are not the pre-specified ones.** The four cardiac chambers originally pre-specified sit behind a commercial licence in the segmentation package. Whole heart plus four great vessels were substituted on the open model before segmentation ran. Chamber-level reproducibility is therefore untested here.

**Partial-field volumes.** Thoracic staging CT clips the great vessels, so vessel volumes are partial-field measurements rather than anatomical totals. Paired comparison remains interpretable only where both members of a pair share coverage, which is why the primary analysis is restricted to coverage-matched pairs.

**The coverage filter is post-hoc.** Coverage matching was not pre-specified. It was forced by a pre-specified assumption proving false, and it costs 186 pairs down to 45 and 58 patients down to 31. Readers who prefer the unfiltered estimate will find it in `endpoint1_results.csv`, where the effect is larger and confounded.

**Threshold fragility, not measurement fragility.** The association endpoint cannot separate reconstruction sensitivity from the behavior of a p < 0.05 cut at n=31. It reports a null against a permutation calibration, and no association survives multiplicity correction, so nothing here supports a claim that reconstruction choice threatens radiogenomic findings.

*Transcriptomic arm.*

**Transcript-only.** Every measurement is mRNA. An endocrine mechanism requires circulating protein, and cardiac transcript induction need not produce elevated plasma ligand. Cell-specific knockout of CCN2 showed that a secreted cardiac factor can act autocrine on neighbouring fibroblasts rather than distally.{^47}

**Platform dependence of the signature.** Discovery and validation shared a sequencing platform, which inflates apparent agreement by 0.17. NPPA replicates on independent array platforms and NPPB does not, so the cardiac-ligand arm rests on NPPA alone once platform is varied.

**Protein-level coverage.** The proteomic test of the receptor was abandoned as underpowered rather than passed or failed, so the screen's central claim is untested at the protein level in the direction that matters.

**Bulk tissue on both sides.** TCGA and GTEx values mix malignant, stromal and immune compartments, so receptor signal cannot be assigned to tumor cells. Compartment markers were extracted but a deconvolution was not performed.

**Units are not comparable across sources.** TCGA RSEM, GTEx TPM, GSE116250 RPKM and GSE141910 normalized counts were never placed on a common scale. Cross-dataset direction claims are therefore ordinal, and the unit-matched tumor-versus-adjacent-normal analysis is the one that supports magnitude statements.

**GTEx heart is non-failing.** Tissue specificity in healthy myocardium bounds baseline distinctiveness rather than HF inducibility. A ligand could be non-cardiac-enriched at baseline yet become a major cardiac output in failure. That possibility is no longer hypothetical, since Section 3.10 tests it directly and recovers PTN, which then fails on independent grounds.

**Induction as distinct from systemic release.** The corrected specificity criterion in Section 3.10 admits ligands on induction magnitude without establishing that the failing heart dominates the circulating pool. Converting tissue mRNA into a circulating contribution requires parameters unavailable for these ligands, so that analysis was pre-emptively excluded rather than attempted with fabricated values.

**Base-rate guard for relaxed criteria.** The seven-axis result was recovered by relaxing one criterion and withdrawn by measuring the base rate of the criterion it was recovered against. Any screen that loosens a filter to recover hits should report the rate at which the unmodified half of the test passes; ours was 74%.

**Bulk versus per-cell induction.** Cardiac PTN induction replicates across two bulk cohorts but does not reproduce as per-cell induction in single-cell data, so the bulk signal may reflect a change in cell composition rather than in per-cell output. That distinction is not resolvable from bulk transcriptomes.

**Power.** Dependency correlations use 22 lineage-level aggregates, so only a large effect would be detectable. The null result bounds effect size rather than excluding a small one.

### 4.6 Reuse beyond this question

The same six skills apply wherever the secreted output of a source tissue is scored against receiver contexts, including any organ-crosstalk screen, any secretome-to-receptor ranking, and any compatibility score across many tissues or cell lines. Nothing in the scoring module depends on the source being heart or the receivers being tumors, since it consumes a feature matrix and a vector of source-derived weights.

What does not transfer automatically is the threshold. The floor is calibrated as a fraction of GAPDH in RSEM units, and a study on a different platform must recalibrate it rather than inherit the number, since the assertion enforces that a floor exists and was applied rather than that 0.005×GAPDH is correct everywhere. The limit of the mechanism is therefore that it renders a requirement unavoidable without making its parameterization universal.

### 4.7 Effects and costs of the multi-agent execution model

Three effects of partitioning the work were observed, the third of which altered a reported result. These are observations from a single execution rather than a controlled comparison: no single-agent baseline, no scripted-workflow baseline and no context-isolation ablation were run, so what follows describes what this configuration did and not what an alternative would have failed to do.

Where the work sits determines how much of it depends on the language model. None of the
4,858 lines contains a call to a language model, an API key or a platform function, so
given the same inputs they return the same numbers whichever agent invoked them, and every
figure and table here is output of that layer. The model contributed orchestration and
judgment: choosing ligand distinctiveness and receptor abundance as the two factors,
recognizing that an unfloored score was measuring the expression axis, and deciding that a
recovered axis had to be tested against a ligand-free base rate. Those choices are
reported so they can be disputed, and a different model making different ones would
produce a different paper from the same modules.

The first effect was increased throughput on independent tracks. The four wave-1 acquisition tracks share no inputs and ran concurrently; serially they would have consumed the sum of their durations rather than the maximum.

The third is that isolation is what allowed the kill criteria to bite. The adversarial phase of Section 3.10 was dispatched with its criteria in the task statement and no access to the reasoning that had produced the seven candidate axes, and it withdrew all seven. We do not claim that a single-context agent could not have done this, but the design removes the mechanism by which it would most plausibly have failed to.

The second effect is context isolation as a correctness property rather than an efficiency one. Each sub-agent received a written brief rather than the orchestrator's conversation and returned artifacts, so a sub-agent could not inherit an unstated assumption from an earlier track, and the audits that caught both methodological errors were run by agents that had not seen the reasoning which produced them.

The architecture claim this paper can support is therefore narrow. The audits found errors and the pre-specified criteria withdrew a positive finding, both under a division of labour in which the auditing agent could not see the audited agent's context. We cannot say the same result would have been unreachable by a single agent with the same guardrails, because that comparison was not run. What the delegation record supports is a cost and a mechanism, not a superiority claim.

### 4.8 Interpretability, reliability, privacy and ethics in delegated analysis

Four concerns recur in discussions of autonomous agents in radiology, and this study bears on
each of them differently.

**Interpretability.** The concern is usually directed at the model. Here the opaque component
was the analysis around it. TotalSegmentator's internal representation is not interpretable and
we did not attempt to make it so; what we made interpretable is the chain of decisions that
turns its output into a claim. Each methodological requirement is a named assertion with a
threshold fixed in advance, so a reader can see which condition bound which result, and the
z-coverage failure is legible precisely because the requirement was written down before the data
were seen. Interpretability of the pipeline is separable from interpretability of the network,
and it is the former that determines whether a reported number can be audited.

**Reliability.** We report three reliability measures rather than asserting robustness: a
segmentation failure rate of 4 of 217 series, all traceable to archive integrity rather than to
the model; a reconstruction-attributable measurement variation of 0.63% for whole heart and
4.96% for pulmonary vein; and an association-level result that proved indistinguishable from a
permutation null. The last is the informative one for delegated analysis. An agent that had
reported the nominal discordance without calibrating it would have produced a plausible,
publishable, and unsupported claim, and nothing in the output would have indicated the problem.

**Privacy.** Both arms use only de-identified, publicly released data under existing data-use
terms, so no protected health information entered the pipeline. That is a property of this study
design and not a solution to the general problem. Agents that fetch, cache and transform imaging
data create copies at every step, and the pre-specification records which series were retrieved
and what was retained; a pipeline operating on institutional data would need that record to be an
enforced access boundary rather than a log.

**Ethical standards.** The specific hazard when analysis is delegated is that a wrong result
arrives with the same confidence as a right one. Our answer is procedural: criteria that can
withdraw a finding, fixed before the finding exists. It was exercised three times, twice at no
computational cost, and once against a result we had already written up. We also declare where
the agent system was used as a method rather than as a writing aid, and which author is
accountable for the reported content. Neither the assertions nor the declaration removes the need
for human judgement; they make the point at which judgement was exercised visible to a reviewer.

## 5. Conclusions

On routine chest CT, automated whole-heart volumetry is reproducible across reconstruction kernels to a patient-level median of 0.63% when reconstructed field of view is matched, and every confidence interval across five thoracic structures excludes the 10% threshold fixed in advance. Thin vascular structures are roughly eight times more sensitive than the heart, with pulmonary vein at 4.96%, which is the magnitude that matters for longitudinal or multicentre work assembled from archives where kernel is not held constant. Segmentation failed on 4 of 217 series, and all four failures were archive integrity faults arising before the model was invoked rather than model failures.

That measurement came out of a design that tests the guardrail mechanism rather than the anatomy. Two of its three endpoints report nulls, and both nulls were reached by criteria fixed before the data could speak: an eligibility check halted the cross-study endpoint at five evaluable patients on metadata alone, and the association endpoint proved indistinguishable from a permutation null once the arbitrary analyst choices were integrated over rather than fixed at one convenient value. The endpoint that returned a measurement returned a reassuring one: automated whole-heart volumetry is reproducible across reconstruction kernels to about 1% when field of view is matched.

The most instructive result was not planned. A requirement we had written into our own pre-specification, that reconstructions of one acquisition share z-coverage, was false, and coverage difference correlates with apparent volume difference at rho 0.70 for the aorta. Written as prose, it could not fail; the analysis would have attributed field-of-view truncation to the reconstruction kernel and reported a larger effect with a clean audit trail. Written as an assertion over the header table, it would have halted before any GPU time was spent. That is the distinction this paper argues for, and we found it by violating it.

Across 30 solid tumor types, no ligand–receptor axis has both a distinctively cardiac source and a receptor enriched in tumor tissue, where enrichment is measured against unit-matched adjacent normal in the 15 types that have it. Enrichment is a specificity criterion: a tumor could respond through a receptor expressed adequately but not enriched, so this is evidence against a tumor-selective cardiac axis rather than against receptor-mediated response in general. The conclusion survives correction of the design choice most capable of having produced it. Scoring cardiac specificity in healthy myocardium cannot detect a ligand induced only in disease, and re-scoring by failure induction recovered seven candidate axes, all involving pleiotrophin. Criteria fixed in advance withdrew all seven, because the receptor criterion was satisfied by 60 of 81 candidate receptors with no ligand involved.

Two further tests were pre-specified and run after the negative was reached. The signature reproduces on independent microarray platforms, but NPPB does not, so the cardiac-ligand arm rests on NPPA alone once platform is varied. A protein-level test of the receptor result was attempted in CPTAC proteomics and abandoned under its own kill criterion, because NPR1 reaches the required coverage in too few tumor types. The three types where it is testable all show lower NPR1 protein in tumor, which is directional agreement with the transcript result but not independent confirmation of it: three types cannot validate a conclusion drawn across 15, and the analysis was abandoned precisely because it is underpowered. The pleiotrophin panel is withdrawn at protein level as it was at transcript level.

The hypothesis is not refuted in general. It remains tenable in the weaker form that heart failure raises circulating concentrations of factors secreted by many tissues, which would require neither a cardiac-specific ligand nor a tumor-specific receptor. Distinguishing that from the specific form tested here requires plasma measurements in heart-failure patients, and for pleiotrophin, the strongest candidate this screen produced, such measurements do not exist. Tumor proteomics does not substitute, since it measures the receiver rather than the circulating signal. The screen also cannot recover a mediator that acts without a curated cognate receptor, is regulated in the opposite direction in human end-stage disease than in the mouse model, or acts through protease activity rather than receptor binding, which is the case of SerpinA3, a mediator nominated by other groups.{^51} The negative applies to receptor-mediated signaling by ligands induced in the failing human ventricle rather than to secreted mediators in general.

The transferable output of this study is the instrument rather than the negative result itself. The three requirements of Section 4.3 were each discovered as an error this pipeline made, and each now exists as an assertion inside a skill rather than as a paragraph here. The distinction is as follows: a study built on the scoring module cannot omit the expression floor or substitute an unmatched null without raising an error that names the offending gene, whereas a study inheriting the same three findings as prose can omit all of them and produce output that looks correct.

What remains untested is reuse by others. Both arms here were built by the same group on the same skills, so we can show that a constraint written as code binds the study that wrote it. Whether it binds a group that inherits the module, and would rather it did not, is the claim that matters and the one we cannot make from two arms of our own.

## 6. Data and code availability

All data are public: GSE116250 and GSE141910 (Gene Expression Omnibus){^48,49} TCGA PanCanAtlas expression and the Thorsson immune landscape, GTEx v8 median tissue expression, Human Protein Atlas protein-class annotations, CellChatDB, and DepMap 24Q4. Accession numbers, endpoints and retrieval dates are given in Section 2.2. Imaging series are public collections in The Cancer Imaging Archive; every SeriesInstanceUID analyzed is listed in `colab_series_manifest.csv`, so the exact scans are retrievable. Per-series segmentation output (`volumes.csv`), the pinned model record with checkpoint hashes (`model_provenance.json`), the dated pre-specification and its amendments, and every endpoint table are in the same repository.

Analysis code is available at https://github.com/nav92061/hf-crosstalk and is archived at
Zenodo (doi:10.5281/zenodo.21968364), which resolves to the most recent release. It is
released as the six skills itemized in Table 1, each comprising a documentation file plus a
Python module. The expression floor, the matched-random null and the conjunction rule are assertions that halt execution on violation rather than documented recommendations.

The skills carry no dependency on the environment used here: no language-model call, no API key, no platform function. They can be driven by another agent framework, or by a plain Python script with no agent at all, and the guardrails hold in every case because they are assertions in the code rather than instructions in a prompt. Continuous integration fails the build if a platform call appears in any module. The repository also carries the result tables, figures, pre-specifications, control log, and build_manuscript.py, which regenerates the typeset manuscript from source. Figures are released at 300 dpi, and every quantity plotted in the imaging figures is traceable to a named table in the repository, so the reported values are checkable against the released data without network access or a language-model call.

What is not released is the agent layer itself. The orchestrator and the 23 sub-agents were runtime configurations of a general-purpose agent on a hosted platform rather than software we wrote, and they cannot be packaged. What is reproducible is the delegation record, comprising the tracks, their waves, their status, duration and token cost, in delegation_record.csv, together with the skills the agents called, which is where the guardrails live. A reader wanting to reproduce the analysis runs the skills directly; a reader wanting to reproduce the system would need to reconstruct the delegation from that record on whatever agent framework they use. We state this explicitly because a paper about agents that ships no agent should say so.

Sixty-eight numbered result tables, 8 figures, the pre-specifications for all three analysis phases, and a log of 264 controls and verdicts, including those that failed, accompany the code.

## Author contributions

**Nicolas Avelar:** Conceptualization, Methodology, Software, Formal analysis,
Investigation, Data curation, Visualization, Writing – original draft, Writing – review
& editing, Project administration, Funding acquisition. **Ruidong Zhang:** Methodology,
Software (original transcriptomic pipeline on which that arm was developed), Writing –
review & editing. **Yi Pan:** Conceptualization, Methodology, Writing – review &
editing.

The analysis was executed by a hierarchical multi-agent system (one orchestrating agent
and 23 sub-agent instances dispatched across five waves) operating on six deterministic,
model-free analytical modules, under human design and audit (Sections 2.12 and 4.7). N.A. specified the guardrails, audited every reported quantity against its source
table, and is responsible for the final content. All authors reviewed the completed
manuscript and take responsibility for the integrity of the analysis, including the
portion executed by the agent pipeline described in Section 2.12.

## Funding

This research received no external funding. All compute and API costs, including the
$164 consumed by the sub-agent fleet at 2026 pricing (Section 2.12), were funded by the
first author.

## Conflicts of interest

The authors declare no conflicts of interest. No employment, consultancy, stock ownership, honoraria, paid expert testimony, patent application, grant or other funding relates to this work, and no author holds an editorial or advisory role with this journal.

## Ethics statement

This study exclusively re-analyzed de-identified, publicly available human genomic, transcriptomic, and proteomic datasets (Gene Expression Omnibus, TCGA PanCanAtlas, CPTAC/Proteomic Data Commons, GTEx, DepMap, CELLxGENE Census; Section 2.2) together with de-identified public CT series from The Cancer Imaging Archive (TCGA-LUSC, TCGA-LUAD, TCGA-ESCA, TCGA-SARC; Section 2.1). No new human subjects data were collected and no additional ethics approval was sought beyond that obtained by the original data-generating studies, consistent with their terms of use.

## Acknowledgments

The authors thank the patients and donor families whose tissue made the primary datasets
possible, and the consortia that curate and release them: the Gene Expression Omnibus,
TCGA and the PanCanAtlas, GTEx, CPTAC and the Proteomic Data Commons, DepMap, the Human
Protein Atlas, and the maintainers of CellChatDB. This work is entirely a reanalysis of
their data.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During the preparation of this work, the authors used Claude (Anthropic) for drafting and
revising the manuscript text. The authors reviewed and edited the output as
needed and take full responsibility for the content of the published article.

The same system also executed the computational analysis. That use is part of the
research process rather than of manuscript preparation, and it is therefore reported as
a method of the study in Sections 2.12 and 4.7, which give the delegation record,
per-track costs and failure modes. Every numerical claim in the manuscript was verified
against its source table and every correction was logged in `validation_summary.csv`
(Section 6). No figure or image was generated or altered by generative AI; all figures
are renderings of the analysis outputs produced by the analysis code.

## References

1. Meijers WC, Maglione M, Bakker SJL, et al. Heart Failure Stimulates Tumor Growth by Circulating Factors. Circulation 2018;138(7):678-691. doi:10.1161/CIRCULATIONAHA.117.030816.

2. Koelwyn GJ, Newman AAC, Afonso MS, et al. Myocardial infarction accelerates breast cancer via innate immune reprogramming. Nat Med 2020;26(9):1452-1458. doi:10.1038/s41591-020-0964-7.

3. Avraham S, Abu-Sharki S, Shofti R, et al. Early Cardiac Remodeling Promotes Tumor Growth and Metastasis. Circulation 2020;142(7):670-683. doi:10.1161/CIRCULATIONAHA.120.046471.

4. Hasin T, Gerber Y, McNallan SM, et al. Patients with heart failure have an increased risk of incident cancer. J Am Coll Cardiol 2013;62(10):881-6. doi:10.1016/j.jacc.2013.04.088.

5. Banke A, Schou M, Videbaek L, et al. Incidence of cancer in patients with chronic heart failure: a long-term follow-up study. Eur J Heart Fail 2016;18(3):260-6. doi:10.1002/ejhf.472.

6. Roderburg C, Loosen SH, Jahn JK, et al. Heart failure is associated with an increased incidence of cancer diagnoses. ESC Heart Fail 2021;8(5):3628-3633. doi:10.1002/ehf2.13421.

7. Selvaraj S, Bhatt DL, Claggett B, et al. Lack of Association Between Heart Failure and Incident Cancer. J Am Coll Cardiol 2018;71(14):1501-1510. doi:10.1016/j.jacc.2018.01.069.

8. Kwak S, Kwon S, Lee SY, et al. Differential risk of incident cancer in patients with heart failure: A nationwide population-based cohort study. J Cardiol 2021;77(3):231-238. doi:10.1016/j.jjcc.2020.07.026.

9. van den Berg PF, Yousif LI, Koop Y, et al. Framingham risk score associates with incident cancer and heart failure. Eur J Prev Cardiol 2026;33(4):490-497. doi:10.1093/eurjpc/zwaf618. Published online 27 September 2025.

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

43. Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. J R Stat Soc Series B 1995;57(1):289-300. doi:10.1111/j.2517-6161.1995.tb02031.x. Not indexed in PubMed; record retrieved from the Crossref REST API.

44. Cox DR. Regression Models and Life-Tables. J R Stat Soc Series B 1972;34(2):187-202. doi:10.1111/j.2517-6161.1972.tb00899.x. Not indexed in PubMed; record retrieved from the Crossref REST API.

45. Harris CR, Millman KJ, Van Der Walt SJ, et al. Array programming with NumPy. Nature 2020;585(7825):357-362. doi:10.1038/s41586-020-2649-2. PMID: 32939066.

46. Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 2020;17(3):261-272. doi:10.1038/s41592-019-0686-2. PMID: 32015543.

47. Dorn LE, Petrosino JM, Wright P, et al. CTGF/CCN2 is an autocrine regulator of cardiac fibrosis. J Mol Cell Cardiol 2018;121:205-211. doi:10.1016/j.yjmcc.2018.07.130.

48. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional genomics data sets--update. Nucleic Acids Res 2012;41(Database issue):D991-5. doi:10.1093/nar/gks1193. PMID: 23193258.

49. Edgar R, Domrachev M, Lash AE. Gene Expression Omnibus: NCBI gene expression and hybridization array data repository. Nucleic Acids Res 2002;30(1):207-10. doi:10.1093/nar/30.1.207. PMID: 11752295.

50. Seuthe K, Picard FSR, Winkels H, Pfister R. Cancer Development and Progression in Patients with Heart Failure. Curr Heart Fail Rep 2024;21(6):515-529. doi:10.1007/s11897-024-00680-y. PMID: 39340596.

51. Caller T, Weiss L, Sharon E, et al. Human heart-tumor interaction in ischemic / failing heart is mediated by cardiac mesenchymal stromal cells [abstract]. Eur Heart J 2025;46(Suppl 1):ehaf784.4672. doi:10.1093/eurheartj/ehaf784.4672.
