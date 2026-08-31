# Restructure: cut the loose biology, fix the two-papers problem

Net effect: **8,986 → ~9,900 words (+10%), Results stays at 10 subsections, figures stay
at 7.** The imaging arm is absorbed rather than appended.

## The two-papers problem is not solved by trimming

Cutting words does not fix it. The problem is that a reader meets a transcriptomic screen
for ten sections and then meets a CT study. The fix is a **shared spine**: both halves are
instances of one claim — that a methodological precondition written as executable code
catches errors that prose caveats do not.

Under that spine the paper has one argument and two demonstrations:

- **Demonstration 1 (transcriptomic).** A pre-specified expression floor and an
  expression-matched null withdrew seven pleiotrophin axes the screen had recovered.
- **Demonstration 2 (imaging).** A pre-specified intersection check halted an endpoint
  before compute was spent, and a reconstruction swap tests whether the screen's own
  receptor panel produces stable imaging associations.

This requires rewriting the Introduction and Abstract — roughly 250 words, replaced not
added. Without it the imaging is an appendix and the relevance objection survives.

## Cuts — 1,060 words, 3 Results subsections, 2 figures

### 3.10 Exploratory: NPR1 as lost tumour suppressor — 604 words, Fig 7
Its own text concedes this tests a hypothesis "close to the opposite" of the design,
generated post-hoc by the refutation. It is a **third direction** in a paper about to carry
two modalities, and a radiology reviewer cannot assess it.
**Keep:** one sentence in 3.9 recording that the reverse-direction hypothesis was tested
and failed its pre-specified conjunction rule. Full analysis → supplement.

### 3.6 Receptor expression vs lineage dependency — 222 words, Fig 4
Introduces a **fourth data modality** (DepMap CRISPR) to support a null the paper already
establishes twice. Nothing in the guardrail argument depends on it.
**Keep:** one clause in 3.5 — receptor compatibility does not track CRISPR dependency
(0 of 81 after FDR). Full analysis → supplement.

### 2.9 Dependency screen — 83 words
Methods for 3.6. Goes with it.

### 3.8 Falsifiable experiment, WNT9A→FZD7 panel — 151 words
Proposes a wet-lab panel the paper does not run, on an axis its own text concedes is **not
cardiac-enriched**. Forward-looking, not load-bearing.
**Keep:** one sentence in 4.5 "What would decide the question."

### Not cut, deliberately
3.3 (expression floor), 3.7 (specificity checks), 3.9 (filter correction and withdrawal)
are the guardrail thesis in transcriptomic form. They are the half of the paper the imaging
arm rhymes with, and cutting them would leave the spine with one leg.

## Additions — 1,975 words

| Section | Content | Words | Status |
|---|---|---|---|
| 2.11 | Imaging cohort, pair confirmation, TotalSegmentator pin, structure amendment | 400 | ready |
| 3.11 | Precondition check halts the cross-study endpoint (35→31→17→5 funnel) | 300 | **ready** |
| 3.12 | Endpoint 1 — reconstruction changes AI volumes | 400 | needs Colab |
| 3.13 | Endpoint 2 — association stability under reconstruction swap | 400 | needs Colab |
| 4.7 | Scope, prior art, what it does not show | 300 | ready |
| 4.6 | 3–4 imaging bullets | 175 | after results |

Figures: old Fig 4 and Fig 7 retire; new **Fig 6** (imaging design: funnel, contrast
cross-tab, design comparison — already rendered) and **Fig 7** (endpoints 1 and 2).

## Result

| | before | after |
|---|---|---|
| Words | 8,986 | ~9,900 |
| Results subsections | 10 | 10 |
| Figures | 7 | 7 |
| Data modalities | 4 (HF arrays, TCGA, GTEx, DepMap) | 3 (HF arrays, TCGA/GTEx, CT) |

The paper gets **shorter in structure while gaining a modality**. Modality count drops
because DepMap leaves with 3.6.

## Order of work

1. Rewrite Introduction + Abstract to the shared spine — **after** Colab returns, since
   the framing depends on whether associations prove stable.
2. Write 2.11, 3.11, 4.7 now; all are fully determined.
3. Make the cuts and move the material to supplement.
4. Fill 3.12 and 3.13 from the Colab output.
