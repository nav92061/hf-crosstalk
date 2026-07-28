# Reproduction guide

Everything in `PAPER.md` traces to a saved artifact produced by the six skills
below. This document records what was run, against which data, with which seeds,
and what to expect when re-running.

## 1. Skills

Six skills, each `skills/<name>/SKILL.md` + `kernel.py`, all published to the
catalog and loadable with `skill({skill: "<name>"})`.

| Skill | Purpose | Key guardrail |
|---|---|---|
| `geo-bulk-de` | GEO series download, metadata parsing, two-group DE | Raises on duplicate sample IDs, groups < 2, all-NA genes |
| `hpa-secretome` | HPA protein-class annotation, secretion tiering | Raises unless the response is a TSV with >10,000 rows (catches error pages) |
| `cellchat-lr` | Parse `CellChatDB.human.rda`, expand receptor complexes | Raises if the interaction table has <1,000 rows |
| `tcga-pancan` | Memory-bounded PanCanAtlas streaming, barcode mapping | Raises if no gene rows matched or a sample lacks a tumour-type assignment |
| `crosstalk-score` | Two-factor compatibility scoring | `require_expression_floor()` **raises**; `matched_random_null()` matches on expression decile |
| `depmap-lineage` | DepMap release resolution, dependency by lineage | Raises if a download is HTML; lineage mapping must be passed, never inferred |

**Portability.** No `kernel.py` references any platform API (`host.*`). Each is
pure stdlib plus numpy/pandas/scipy (plus `requests` for `depmap-lineage`), so
they run under a different agent harness, or with no agent:

```bash
pip install numpy pandas scipy statsmodels requests
python -c "
import sys; sys.path.insert(0, 'skills/crosstalk-score')
import kernel as ck, pandas as pd
expr  = pd.read_csv('tcga_expr_by_tumor_type.csv', index_col=0)
inter = pd.read_csv('cellchatdb_receptor_subunits.csv')
de    = pd.read_csv('hf_signature_replicated.csv')
avail = ck.ligand_availability(de, inter['ligand_subunit_gene'].unique(),
                               gene_col='gene_symbol')
print(ck.crosstalk_score(expr, inter, avail, enforce_floor=False)['score'])
"
```

Guardrails are code-level assertions rather than prose instructions
specifically so they survive being driven by a weaker model, or by no model.

## 2. Data sources and access

| Source | Endpoint | Notes |
|---|---|---|
| GSE116250 | `ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/` | **RPKM only**, no counts. 65 genes carried a sentinel `999999999` and were removed. |
| GSE141910 | same host, `GSE141910_RAW.tar` | Per-sample CSVs, log-scale normalized values (MAGNet). |
| HPA secretome | `www.proteinatlas.org/api/search_download.php` | Protein-class strings parsed; no hardcoded gene lists. |
| CellChatDB | `raw.githubusercontent.com/sqjin/CellChat` | `CellChatDB.human.rda`, parsed in Python. |
| TCGA PanCanAtlas | `api.gdc.cancer.gov/data/3586c0da-…` | ~1.9 GB; **RSEM normalized_count, not TPM** (GAPDH median ≈63,656). Streamed row-wise. |
| Thorsson immune landscape | `api.gdc.cancer.gov/data/1a7d7be8-…` | Substituted for TIMER2. |
| DepMap 24Q4 Public | figshare article `27993248` | `depmap.org` serves a bot-verification page; resolve via the figshare API. |
| GTEx v8 | expression connector | Median TPM, 54 tissues (52 after excluding cultured-cell entries). |

**Network grants required** (approved during this work): `api.gdc.cancer.gov`,
`depmap.org`, `api.figshare.com`, `ndownloader.figshare.com`, `compbio.cn`.

On TIMER2, correcting an earlier error in this project's notes: `compbio.cn` (the
host the user supplied) **is** grantable, and was granted. An initial probe had
tested a different TIMER host, `timer.comp-genomics.org`, and recorded TIMER2 as
"not grantable" without ever requesting the real domain. Having requested it, the
site responds — but TIMER2 is an interactive Shiny application with no bulk
download route (`/timer2/data/`, `/timer2/estimation.php` → 404), so pan-cancer
deconvolution estimates cannot be retrieved programmatically. The Thorsson
immune-landscape table is used instead, and the paper now states this reason
rather than claiming the domain was unreachable.

## 3. Pipeline order

1. **HF signature** — `geo-bulk-de` on GSE116250 (failing vs non-failing; DCM-only and ICM-only as secondary), then GSE141910. Replication statistics plus a label-permutation control (seed 42).
2. **Annotation** — `hpa-secretome` genome-wide; `cellchat-lr` parsed and expanded.
3. **Tumour axis** — `tcga-pancan` streams receptor, ligand, housekeeping, pathway-target and compartment-marker rows; 9,538 primary solid tumours in 30 types (LAML, DLBC, THYM excluded); 396 expression-matched background genes for the null.
4. **Scoring** — `crosstalk-score`. Floor = 0.005×GAPDH = 345.7. Primary, floor-filtered and double-centred variants; bootstrap 400; matched-random null 300 iterations (seed 20260727).
5. **Validation** — `depmap-lineage` (81 receptors × 22 lineages, with common-essential positive controls); GTEx tissue specificity; TCGA adjacent-normal unit-matched comparison (710 samples, 15 types).

## 4. Joins that will silently return nothing if done wrong

- **`hf_signature_replicated.csv` has both `gene` (Ensembl) and `gene_symbol`.** Join on **`gene_symbol`**. Joining on `gene` returns zero rows against HPA/CellChatDB/TCGA — it looks like a biological result and is not.
- **TCGA values are RSEM normalized counts.** Never apply a TPM-scaled threshold. Use `housekeeping_floor()`.
- **TCGA and GTEx cannot be subtracted or ratioed.** Rank-based comparison, or the unit-matched adjacent-normal route.
- **The hypergeometric overlap test is near-saturated here** (p=2.1×10⁻³⁵ at 1.07× enrichment). Report directional concordance and effect-size correlation.
- **Two correlation values are both correct** and must not be conflated: ρ=0.541 across all 18,061 shared genes; ρ=0.759 among the 5,717 significant in both.

## 5. Seeds and expected variation

| Step | Seed | Iterations |
|---|---|---|
| DE label permutation | 42 | 1 |
| Matched-random specificity null | 20260727 | 300 |
| Bootstrap CIs | 20260727 | 400 |
| TCGA background gene sample | recorded in `tcga_background_genes.json` | 396 genes |

The null and bootstrap are stochastic; expect ρ_mean and CI bounds to move in
the third decimal. The qualitative results — floor null ≈0.11 versus ≈0.36
without, zero axes surviving both specificity checks, zero DepMap hits after
FDR — are not marginal and should reproduce.

## 6. Environment

`vibemed` (conda): pandas, numpy, scipy, statsmodels, matplotlib, seaborn,
scikit-learn, requests, pyreadr, rdata, openpyxl. 8 GiB RAM, 8 cores, no GPU, no
remote compute — which is why the PanCanAtlas matrix is streamed and DepMap
files are read with `usecols`.

## 7. Outputs

- `TABLE_INDEX.csv` — 22 numbered tables with artifact version IDs
- `FIGURE_INDEX.csv` — 5 figures with captions
- `validation_summary.csv` — 25 controls plus 6 audit corrections: expected, observed, verdict, and what each bears on
- `literature_citations.csv` / `literature_gaps.json` — 107 verified PMIDs, and claims the literature does *not* support
- `PAPER.md` — manuscript

## 8. Note on a discarded prior draft

An earlier `PAPER.md` in this project described a twelve-skill pipeline and cited
31 artifacts. An audit found no registered skills and none of the cited figures
or tables in the artifact store, in this project or any other; its numbers had no
backing data or lineage. At the user's direction it was discarded, and nothing
from it is reused here. Its one methodological point — that a co-expression
ranking can be reproduced by random receptor sets — was treated as a hypothesis
and tested; §2.3 of the paper reports what we actually found, which is that the
confound is caused by scoring untranscribed receptors and is largely removed by
an expression floor.
