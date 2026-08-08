# hf-crosstalk

**A multi-agent bioinformatics pipeline with executable methodological guardrails.** Six skills — each a
documentation file an agent reads to decide whether the procedure fits its task, plus a
Python module that executes it — applied to the question of whether the failing heart
releases a signal any of 30 solid tumor types is equipped to receive.

**What ran this, and what ships.** The analysis was executed by a hierarchical
multi-agent system: one orchestrating agent and 23 isolated sub-agents in five waves
over 91 hours, $164 in API cost at a 106:1 input-to-output token ratio. 19 of 23 tracks
completed normally; all four that did not had already persisted their results, so no
analysis was lost. That record is in `results/tables/delegation_record.csv`. The agent
layer itself is **not** shipped — it was a runtime configuration of a hosted
general-purpose agent, not software written here. What ships is the six skills the
agents called, which is where the guardrails live, plus `build_manuscript.py` to
regenerate the PDF.

**The architecture.** 4,858 lines of Python, 114 functions, 141 guardrail checks written
as assertions that halt execution rather than as advice in a methods section. No
language-model call, no API key, no platform function in any module: the same inputs
return the same numbers whichever agent invoked them, or none. Every figure and table in
the paper is output of that layer.

**Two guardrails exist because the pipeline made the error first.** Scoring without an
absolute expression floor let expression-matched random receptor sets reproduce the
tumor-type ranking at rho=0.363, with 21.5% of draws exceeding 0.5; a floor at
0.005xGAPDH cut this to 0.107 and 1.7%. An unmatched random pool gave 0.323, so the null
must be matched as well as floored. Standardizing an untranscribed gene still yields
large z-scores, so the score had been ranking receptors the tumor never transcribes.

**The stress test.** Re-scoring the source criterion recovered seven candidate axes, all
involving pleiotrophin. Kill criteria fixed before that analysis ran withdrew every one,
because 45 to 77% of candidate receptors satisfied the receptor criterion with no ligand
involved. The pipeline invalidated its own positive result through a check written
before the result existed.

**The biological result is a negative, and it stands on its own.** Across 30 solid tumor
types no ligand-receptor axis has both a distinctively cardiac source and a receptor
enriched in tumor. NPR1, the sole floor-passing receptor of the two cardiac-enriched
ligands NPPA and NPPB, is depleted in tumor versus matched adjacent normal in all 15
testable types (median log2FC -1.91, all FDR<0.05).

**Reuse.** The skills apply to any secretome-to-receptor screen: nothing in the scoring
module depends on the source being heart or the receivers being tumors. The floor
threshold is calibrated in RSEM units against GAPDH and must be recalibrated on another
platform — the assertion enforces that a floor was applied, not that this number is
right everywhere.

## The six skills

Each is a `SKILL.md` (instructions) plus a `kernel.py` of plain Python. **No
`kernel.py` imports or references any agent-platform API** — they are pure
stdlib + numpy/pandas/scipy (plus `requests` for `depmap-lineage`), so they run
under any agent harness or as ordinary scripts.

| Skill | Purpose | Guardrail enforced in code |
|---|---|---|
| [`geo-bulk-de`](hf_crosstalk/skills/geo-bulk-de) | GEO series download, metadata parsing, two-group DE | Raises on duplicate sample IDs, groups < 2, all-NA genes |
| [`hpa-secretome`](hf_crosstalk/skills/hpa-secretome) | Protein Atlas annotation, secretion tiering | Raises unless response is a TSV with >10,000 rows (catches error pages) |
| [`cellchat-lr`](hf_crosstalk/skills/cellchat-lr) | Parse CellChatDB, expand receptor complexes | Raises if the interaction table has <1,000 rows |
| [`tcga-pancan`](hf_crosstalk/skills/tcga-pancan) | Memory-bounded PanCanAtlas streaming | Raises if no rows matched or a sample lacks a tumor-type assignment |
| [`crosstalk-score`](hf_crosstalk/skills/crosstalk-score) | Two-factor scoring, induction-based source criteria, conjunction verdicts | **`require_expression_floor()` raises**; null matches on expression decile; `conjunction_verdict()` refuses to call a partial pass a pass; `receptor_pass_rate()` reports the ligand-free base rate |
| [`depmap-lineage`](hf_crosstalk/skills/depmap-lineage) | DepMap resolution, dependency by lineage | Raises if a download is HTML; lineage mapping must be passed, never inferred |

Guardrails are **assertions that raise**, not advice in prose — so they hold
whichever model, or no model, is driving.

## Quick start

```bash
pip install -e .                 # or: pip install git+https://github.com/USER/hf-crosstalk.git
python tests/test_guardrails.py  # verifies the floor validator actually rejects
```

The skills are also usable without installing anything — every `kernel.py` is a
standalone file:

```bash
python -c "
import importlib.util as u
s = u.spec_from_file_location('k', 'hf_crosstalk/skills/crosstalk-score/kernel.py')
m = u.module_from_spec(s); s.loader.exec_module(m)
print([f for f in dir(m) if not f.startswith('_')])"
```

Scoring, standalone:

```python
import sys; import pandas as pd
from hf_crosstalk import crosstalk_score as ck

expr  = pd.read_csv("results/tables/tcga_expr_by_tumor_type.csv", index_col=0)
inter = pd.read_csv("results/tables/cellchatdb_receptor_subunits.csv")
de    = pd.read_csv("results/tables/hf_signature_replicated.csv")

# NOTE: join on gene_symbol, not gene (which holds Ensembl IDs)
avail = ck.ligand_availability(de, inter["ligand_subunit_gene"].unique(),
                               gene_col="gene_symbol")
floor = ck.housekeeping_floor(expr, "GAPDH", 0.005)
res   = ck.crosstalk_score(expr, inter, avail, floor=floor, enforce_floor=False)
print(res["score"].head())
```

## Repository layout

```
hf_crosstalk/      importable package
  skills/          six skills (SKILL.md + kernel.py), no platform dependency
results/tables/    77 CSVs — DE tables, network, ranking, sensitivity sweeps,
                   protein tests, and validation_summary.csv (161 controls)
results/figures/   9 figures (7 main + 2 supplementary)
results/summaries/ JSON provenance for each pipeline stage
PAPER.md           manuscript source (10 sections; architecture first,
                   screen as case study)
SUPPLEMENT.md      S1-S7: additional methods, orthogonal checks, protein
                   detail, exploratory NPR1, sensitivity tables, limma
                   validation, analysis timeline
hf-crosstalk-manuscript.pdf   typeset two-column PDF, 16 pages, 57 references
PRESPECIFICATION.md           criteria for the corrected specificity phase
PRESPECIFICATION_PTN.md       kill criteria for the PTN adversarial phase
PRESPECIFICATION_REPLICATION.md  cross-platform replication thresholds
ANALYSIS_TIMELINE.md          artifact-store timestamps for every phase
PTN_VERDICT.md                the withdrawal, with all three kill criteria
REPRODUCTION.md               endpoints, seeds, and the traps
tests/             guardrail regression tests (run with pytest)
```

## Two traps that fail silently

1. **`hf_signature_replicated.csv` has both `gene` (Ensembl) and `gene_symbol`.**
   Join on **`gene_symbol`**. Joining on `gene` returns zero rows against
   HPA/CellChatDB/TCGA — it looks like a biological result and is not.
2. **TCGA PanCanAtlas values are RSEM normalized counts, not TPM.** GAPDH's
   median here is 69,134. A "2 TPM" threshold is meaningless; use
   `housekeeping_floor()`.

TCGA and GTEx also cannot be subtracted or ratioed — use rank-based comparison,
or the unit-matched adjacent-normal route.

## Honesty notes

`results/tables/validation_summary.csv` records **153 controls, pre-specified
verdicts and audit corrections**: what each expected, what was observed, the verdict, and what it
bears on. Failed controls are reported, not hidden — several are the most
informative results here.

Two audits were run. The first found six numeric errors (rows 26–31); the
second used a deliberately different method — recomputing from source tables
rather than checking prose against summaries — and found five more (rows 51–58),
all from quoting a delegated summary sentence instead of the column it describes.
Both are listed in Appendix A of `PAPER.md` rather than silently fixed. One made
the primary ranking claim *weaker* than first written (4 of 30 tumor types beat
the specificity null, not 6); another corrected the adjacent-normal count from
710 to 680.

The PTN phase is the clearest illustration of why the criteria are pre-specified:
a relaxation recovered seven axes, and a base-rate check written before the tests
ran withdrew them. Any screen that loosens a criterion to recover hits should
report the rate at which the *unmodified* half of the test passes.

This study is hypothesis-generating and **does not** show that heart failure
causes cancer. No dataset here links HF status to cancer incidence.

## Data sources

GSE116250, GSE141910, Human Protein Atlas, CellChatDB, GDC PanCanAtlas
(expression + Thorsson immune landscape), DepMap 24Q4 Public (figshare article
27993248), GTEx v8. Endpoints, access notes, and seeds in `REPRODUCTION.md`.
No raw source downloads are vendored — every table here is derived.

## Citation

If the skills or the floor/null methodology are useful, please cite this
repository and the underlying data sources listed in `REPRODUCTION.md`.

## License

MIT for the code (`skills/`, `tests/`). Derived data tables in `results/` remain
subject to the terms of their upstream sources (TCGA, GTEx, DepMap, HPA,
CellChatDB).
