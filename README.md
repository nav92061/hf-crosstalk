# hf-crosstalk

**Does the failing heart meet a tumour-specific receptor?** A pan-cancer
ligand–receptor compatibility screen across 30 solid tumour types — and six
reusable, platform-independent analysis skills.

**Headline result: no.** Of 42 heart-failure-upregulated secreted ligands, only
NPPA and NPPB are cardiac-enriched and abundant. Their sole floor-passing
receptor, NPR1, is *depleted* in tumour versus adjacent normal in all 15
unit-matched tumour types (median log₂FC −1.91, all FDR<0.05). **Zero of 160
candidate axes survive both specificity checks.**

**Headline method: the expression floor is load-bearing.** A receptor
compatibility score built by z-standardising expression will rank untranscribed
receptors highly, because z-scoring noise still yields large z-values.

| Null configuration | Mean ρ vs observed ranking | Draws with ρ>0.5 |
|---|---|---|
| No expression floor | 0.363 | 21.5% |
| **Floor at 0.005×GAPDH** | **0.107** | **1.7%** |
| Unmatched (naive) random pool | 0.323 | 17.0% |

The "global expression axis" confound is substantially *caused by* scoring
receptors that are not transcribed. Remove them and the ranking becomes
receptor-specific. The null pool must also be expression-matched — a naive pool
gives ρ=0.323 and would wrongly condemn a specific score.

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
| [`tcga-pancan`](hf_crosstalk/skills/tcga-pancan) | Memory-bounded PanCanAtlas streaming | Raises if no rows matched or a sample lacks a tumour-type assignment |
| [`crosstalk-score`](hf_crosstalk/skills/crosstalk-score) | Two-factor compatibility scoring | **`require_expression_floor()` raises**; null matches on expression decile |
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
  skills/          six skills (SKILL.md + kernel.py)
results/tables/    29 CSVs — DE tables, network, ranking, validation
results/figures/   5 publication figures
results/summaries/ JSON provenance for each pipeline stage
PAPER.md           manuscript, with a supported/not-claimed section
REPRODUCTION.md    data endpoints, seeds, and the traps below
tests/             guardrail regression test
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

`results/tables/validation_summary.csv` records **31 controls and audit
corrections**: what each expected, what was observed, the verdict, and what it
bears on. Failed controls are reported, not hidden — several are the most
informative results here.

Six numeric errors found during an internal audit are logged as rows 26–31 and
listed in §6 of `PAPER.md` rather than silently fixed. One of them made the
primary ranking claim *weaker* than first written (4 of 30 tumour types beat the
specificity null, not 6).

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
