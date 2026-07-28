---
name: depmap-lineage
description: Resolve a DepMap public release via the figshare API (depmap.org serves a bot-verification interstitial), download CRISPR gene-effect / expression / model-annotation files with memory-safe column subsetting, compute per-lineage mean gene effect, correlate it against an external per-lineage score with BH-FDR, and select high/low-receptor cell-line panels. Load when asked whether a gene is a lineage-specific dependency, when correlating dependency with an external tumour-type score, or when building a testable cell-line panel.
---

# depmap-lineage

CRISPR dependency and expression from DepMap public releases, keyed to lineage,
with an explicit bridge to externally derived per-tumour-type scores.

## Why figshare, not depmap.org

`depmap.org` (including `depmap.org/portal/api/download/files`) returns a
bot-verification HTML page to non-browser clients, so a naive download succeeds
with HTTP 200 and writes an HTML document where a CSV was expected. DepMap
deposits each quarterly public release on figshare+ as a single article titled
`DepMap <YYQN> Public`, and `ndownloader.figshare.com` serves the files
directly. `resolve_depmap_release` therefore goes through
`https://api.figshare.com/v2/articles/search` (POST, `search_for`), keeps
articles whose title matches `^DepMap (\d\dQ\d) Public$`, picks the highest
quarter tag (or the one you asked for), then reads the file list from
`GET /v2/articles/{id}`. `fetch_depmap_file` asserts the payload is not HTML,
so an interstitial fails loudly instead of contaminating the analysis.

Do not hardcode a release. Resolve it and record what you got.

## Functions

```python
resolve_depmap_release(release=None, cache_path=None, session=None, timeout=60)
    -> {"release","article_id","doi","title","published_date","files","sizes","resolved_via"}
fetch_depmap_file(url, dest, expected_min_bytes=1000, session=None, chunk=4<<20, timeout=120) -> path
parse_gene_columns(header) -> {"TP53 (7157)": "TP53", ...}
load_gene_effect(path, genes=None) -> DataFrame  (ModelID x symbol, Chronos)
load_expression(path, genes=None) -> DataFrame   (ModelID x symbol, log2(TPM+1))
load_model_annotation(path) -> DataFrame indexed by ModelID
lineage_dependency(gene_effect, model_df, genes=None, lineage_col="OncotreeLineage", min_lines=1)
    -> long [gene, lineage, mean_gene_effect, median_gene_effect, n_lines, n_dependent_lines]
correlate_with_external(per_lineage_effect, external_score, lineage_map, min_lineages=5,
                        effect_col="mean_gene_effect")
    -> [gene, spearman_rho, p_value, n_lineages, fdr]
select_panel(expr, gene_effect, model_df, receptor, cofactors=(), lineage=None,
             high_pct=85, low_pct=15, activity_genes=(), exclude_subtypes=(),
             lineage_col="OncotreeLineage", subtype_col="OncotreeSubtype", n_per_arm=4)
    -> panel DataFrame (arm, cell_line, ModelID, subtype, expr, percentile, cofactors,
       baseline activity, receptor gene effect)
bh_fdr(pvals) -> np.ndarray of BH-adjusted p-values
```

## Memory discipline (load-bearing on small machines)

`CRISPRGeneEffect.csv` is ~430 MB and
`OmicsExpressionProteinCodingGenesTPMLogp1.csv` ~510 MB: ~1800 rows by ~19 000
gene columns. Loading either whole costs several GB of RAM. Always pass
`genes=[...]`; the loaders resolve `"SYMBOL (ENTREZ)"` headers to symbols and
read only those columns via `usecols`, keeping peak RSS in the tens of MB. Pass
`genes=None` only if you know you have the memory.

## Preconditions that raise

- `resolve_depmap_release`: the search returns at least one `DepMap <YYQN>
  Public` article; the requested release exists; the article carries
  `CRISPRGeneEffect.csv`, `OmicsExpressionProteinCodingGenesTPMLogp1.csv` and
  `Model.csv`.
- `fetch_depmap_file`: payload does not start with an HTML doctype/`<html>` and
  is at least `expected_min_bytes`.
- `load_gene_effect` / `load_expression`: every requested gene appears in the
  file header (the error names the missing symbols).
- `lineage_dependency`: `lineage_col` exists; ModelIDs overlap; no selected line
  has a null lineage.
- `correlate_with_external`: `lineage_map` is supplied non-empty — the
  DepMap-lineage to external-type mapping is a disclosed judgement call and is
  never inferred.
- `select_panel`: receptor present in the expression matrix; pool large enough
  for both arms; neither arm empty; every selected line has a lineage.

## Worked example

```python
import os, json, pandas as pd
from kernel import (resolve_depmap_release, fetch_depmap_file, load_gene_effect,
                    load_expression, load_model_annotation, lineage_dependency,
                    correlate_with_external, select_panel)

man = resolve_depmap_release(cache_path="data/depmap_manifest.json")
print(man["release"], man["article_id"])          # e.g. 24Q4 27993248

paths = {}
for name in ("CRISPRGeneEffect.csv",
             "OmicsExpressionProteinCodingGenesTPMLogp1.csv", "Model.csv"):
    paths[name] = fetch_depmap_file(man["files"][name], os.path.join("data", name),
                                    expected_min_bytes=500_000)

genes = ["FZD1", "LRP6", "RPL3", "POLR2A", "PLK1"]   # controls included on purpose
ge  = load_gene_effect(paths["CRISPRGeneEffect.csv"], genes=genes)
ex  = load_expression(paths["OmicsExpressionProteinCodingGenesTPMLogp1.csv"], genes=genes)
mod = load_model_annotation(paths["Model.csv"])

per_lin = lineage_dependency(ge, mod, genes=genes, min_lines=5)
# positive controls must be strongly negative, else parsing is wrong:
assert per_lin.query("gene == 'RPL3'")["mean_gene_effect"].max() < -0.5

lmap = pd.DataFrame({"depmap_lineage": ["Ovary/Fallopian Tube", "Prostate"],
                     "tcga_type":      ["OV",                   "PRAD"]})
res = correlate_with_external(per_lin, {"OV": 0.43, "PRAD": 0.27}, lmap, min_lineages=5)

panel = select_panel(ex, ge, mod, receptor="FZD1", cofactors=["LRP6"],
                     lineage="Ovary/Fallopian Tube",
                     activity_genes=["AXIN2", "LGR5", "NKD1"], n_per_arm=4)
```

## Interpretation guardrails

- Chronos gene effect is scaled so median common-essential = -1 and
  non-essential = 0. `< -0.5` is the conventional dependency call; the module
  exposes it as `DEPENDENCY_THRESHOLD`.
- Always score known common essentials (`RPL3`, `POLR2A`, `PLK1`) alongside your
  genes of interest. If they are not strongly negative, the header parse or the
  join is wrong and every other number is untrustworthy.
- A receptor that is expressed but has gene effect near 0 is a candidate
  signalling conduit, not a drug target. Do not describe such a gene as a
  vulnerability.
- `correlate_with_external` correlates lineage-level aggregates, so n is the
  number of mapped types (typically 10-20). Power is low; report FDR and the
  null result plainly rather than reading through to nominal p-values.

## Standalone use (no agent)

`kernel.py` has no platform dependencies — stdlib plus pandas, numpy, scipy and
requests. Copy it next to your script and `from kernel import ...` (or
`sys.path.insert(0, "skills/depmap-lineage")`). Network access to
`api.figshare.com` and `ndownloader.figshare.com` is required for resolution and
download; once the three CSVs are on disk, everything else works offline.

```bash
python - <<'PY'
import sys; sys.path.insert(0, "skills/depmap-lineage")
from kernel import resolve_depmap_release
m = resolve_depmap_release(cache_path="depmap_manifest.json")
print(m["release"], m["doi"], len(m["files"]), "files")
PY
```
