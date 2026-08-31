# Prior art: has this methodology been used elsewhere?

**Short answer: yes, every individual component has decades of precedent. What appears
novel is the specific combination and the venue.** A reviewer who knows any one of these
literatures will recognise the mechanism, so the paper should cite them rather than
present the idea as new.

26 references verified against OpenAlex (DOIs in `prior_art_table.csv`). Below, what each
field already established and how the imaging arm differs.

## 1. Halting on a violated precondition — software engineering, since 1992

Meyer's *design by contract* (10.1109/2.161279, 2135 citations) is the direct ancestor:
routines declare preconditions, and a violated precondition is a **failure**, not a
warning. This is exactly the "assert rather than caveat" mechanism, thirty years old.

In scientific computing specifically, Wilson et al. *Best Practices for Scientific
Computing* (10.1371/journal.pbio.1001745) and Taschuk & Wilson *Ten simple rules for making
research software more robust* (10.1371/journal.pcbi.1005412) both recommend
fail-fast assertions in research code.

**Nothing about "requirements as executable assertions" is new.** The paper must not claim
it is.

## 2. Pipelines that refuse on bad input — production ML, since ~2017

- Breck et al., *The ML test score* (10.1109/bigdata.2017.8258038) — a rubric where
  input-validation tests are scored production-readiness requirements
- Breck et al., *Data Validation for Machine Learning* — Google's TFX validation layer,
  which **blocks training** on schema violation
- Schelter et al., *Automating large-scale data quality verification*
  (10.14778/3229863.3229867) — Amazon Deequ: declarative "unit tests for data" that fail a
  pipeline
- Amershi et al. (10.1109/icse-seip.2019.00042), Mitchell et al. *Model Cards*
  (10.1145/3287560.3287596), Gebru et al. *Datasheets* (10.1145/3458723)

Model Cards are the closest conceptual sibling and the most important contrast: they
document intended use and out-of-scope applications **in prose**. The paper's claim should
be positioned precisely here — not "we invented scope documentation" but *"a Model Card
that the code enforces instead of a Model Card the reader is trusted to have read."*

## 3. Refusing to answer outside the validated range — three independent traditions

This is the single most important precedent, because it is the imaging arm's exact
mechanism and it already has a name in three fields:

- **Reject option** — Chow 1970 (10.1109/tit.1970.1054406) formalised abstention 55 years
  ago; surveyed by Hendrickx et al. 2024 (10.1007/s10994-024-06534-x); SelectiveNet
  (Geifman & El-Yaniv, arXiv:1705.08500)
- **Applicability domain** — regulatory QSAR. A model's prediction is only accepted for
  compounds inside a defined descriptor space; outside it, the model **returns nothing**.
  This is an OECD validation principle, i.e. legally binding in some contexts
  (Sahigara et al., 10.3390/molecules17054791)
- **Conformal prediction** — distribution-free validity guarantees
  (Angelopoulos & Bates, 10.1561/2200000101)

**The QSAR applicability domain is the strongest precedent for the paper's thesis** and
should be cited prominently: a whole regulated field already requires that a model declare
its domain and refuse outside it. The paper's argument is that radiology AI has no
equivalent convention.

## 4. Knowing when a segmentation is wrong — imaging, since 2017

- Valindria et al., *Reverse Classification Accuracy* (10.1109/tmi.2017.2665165) —
  predicts per-case segmentation quality **without ground truth**
- Robinson et al. (10.1186/s12968-019-0523-x) — automated QC on UK Biobank cardiac MR
- Mårtensson et al. (10.1016/j.media.2020.101714) — *The reliability of a deep learning
  model in clinical out-of-distribution MRI data*, the closest imaging analogue
- Wasserthal et al., TotalSegmentator itself (10.1148/ryai.230024)

This literature is well developed and largely **per-case, post-hoc**: score the output,
then decide. The paper's variant is **pre-hoc and cohort-level** — check acquisition
parameters against the validated envelope *before* segmenting, and refuse. Complementary
rather than competing, and the distinction is worth one sentence.

## 5. That acquisition changes measurements — radiomics, extensively established

- Zwanenburg et al., IBSI (10.1148/radiol.2020191145, 4059 citations)
- Traverso et al. (10.1016/j.ijrobp.2018.05.053), Berenguer et al.
  (10.1148/radiol.2018172361), Shafiq-ul-Hassan et al. (10.1002/mp.12123) — slice
  thickness, kernel and voxel size all shift extracted features

**The empirical finding is not novel and must not be sold as one.** Radiomics settled this
for hand-crafted features. What is less documented is the same question for a *deep-learning
anatomical segmentation* tool's volume outputs, measured within-patient. That is a
narrower gap and the honest framing.

## 6. Dataset shift in clinical AI — the motivating problem

Finlayson et al., *The Clinician and Dataset Shift in AI* (NEJM, 10.1056/nejmc2104626) and
Kelly et al. (10.1186/s12916-019-1426-2) both argue deployed models fail silently under
shift. This is the "why it matters" citation.

## 7. Pre-specification — established, and the paper's own second mechanism

Nosek et al., *The preregistration revolution* (10.1073/pnas.1708274114); Munafò et al.,
*A manifesto for reproducible science* (10.1038/s41562-016-0021). `PRESPECIFICATION_IMAGING.md`
is a registered-report-style commitment. The novelty is not pre-specification but that the
kill criterion is **machine-checked** rather than honour-system.

## 8. Nearest overall match found

**Schintke et al. 2024, *Validity constraints for data analysis workflows*, Future
Generation Computer Systems (10.1016/j.future.2024.03.037).** This is the closest single
paper: it proposes declaring "validity constraints" alongside a workflow so that hidden
assumptions become explicit and violations cause **graceful termination with a meaningful
error** rather than silent wrong results. Implemented as a proof of concept in Nextflow.

The overlap with the paper's thesis is substantial and it **must be cited**. Differences:
it is a workflow-language contribution aimed at portability and cluster execution, with no
scientific claim attached; the paper applies the pattern to *methodological* validity
(cohort scope, confounding structure, effect thresholds) in a domain-science manuscript.

## What is actually novel — the defensible claim

1. **Not** the assert-instead-of-caveat mechanism (Meyer 1992; Schintke 2024).
2. **Not** refuse-outside-validated-range (Chow 1970; QSAR applicability domain).
3. **Not** acquisition-dependence of imaging measurements (IBSI and successors).
4. **Plausibly novel:** applying the pattern to *methodological* preconditions of a
   biological claim — cohort scope, confounding structure, pre-specified kill thresholds —
   rather than to data schemas or software portability.
5. **Plausibly novel:** demonstrating the guardrails firing to **withdraw the authors' own
   finding** (the seven PTN axes), and reporting that as the result. Pre-registration
   literature advocates this; instances where the machinery actually executed a withdrawal
   are rarer.
6. **Novel for this venue:** carrying one mechanism across a transcriptomic and an imaging
   analysis in a single paper.

## Consequence for the manuscript

The Introduction as currently written is exposed on this point. It should concede items 1-3
explicitly and claim only 4-6. A reviewer from software engineering, cheminformatics, or
radiomics will otherwise recognise their own field's decades-old idea presented as new, and
that impression is harder to recover from than a narrow claim is to defend.

One consequence for the imaging arm specifically: since applicability-domain checking is a
*regulatory requirement* in QSAR, the sharpest framing available is not "we built a
safeguard" but **"radiology AI lacks a convention that another regulated modelling field
has had for twenty years, and here is what it costs."**
