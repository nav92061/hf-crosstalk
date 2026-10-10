---
name: hpa-secretome
description: Retrieve genome-wide Human Protein Atlas protein-class and secretome annotations and classify gene symbols into secreted / plasma-detectable / structural-ECM tiers, then score tissue-of-origin specificity against GTEx medians. Use when you need to decide which genes in a differential-expression result encode proteins that could plausibly leave their tissue of origin and travel through blood, or which tissue actually sources a secreted protein - e.g. filtering a tissue signature down to a candidate endocrine secretome. Encodes the GTEx v2 API silent-empty-response trap, the enrichment-ratio artifact that passes liver-dominant proteins as tissue-enriched, and the Fisher base-rate guard a candidate screen needs before its hits are interpreted.
---

# hpa-secretome

## Purpose

Given any list of human gene symbols (or the whole genome), answer: **is the
protein secreted, is it detectable in plasma, and is it a signalling factor
rather than structural matrix?**

The resource is the Human Protein Atlas download API, which returns a
`Protein class` field containing curated/predicted class tokens, plus an
independent `Secretome location` / `Secretome function` axis. All classification
is done by parsing those strings — no gene lists are hardcoded anywhere in this
skill except the documented structural-ECM symbol families.

## Data source

```
https://www.proteinatlas.org/api/search_download.php?search=<query>&format=tsv&columns=<codes>&compress=no
```

`search=` empty returns the full genome-wide table (20,162 genes as of the
version fetched here). **The API silently drops unrecognised column codes** —
always inspect the header of what came back rather than trusting the request.
Verified working codes are in `HPA_COLUMN_CODES`; the default request is
`g,eg,pc,secl,scml,rnats,blconcms` which yields these columns:

`Gene`, `Ensembl`, `Protein class`, `Secretome location`, `Secretome function`,
`Subcellular main location`, `RNA tissue specificity`,
`Blood concentration - Conc. blood MS [pg/L]`.

Note `secl` is a single code that expands into **two** columns (Secretome
location + Secretome function).

## Classification scheme

Derived from exact token matches inside the comma-separated `Protein class`
string:

| column | HPA token |
|---|---|
| `is_secreted` | `Predicted secreted proteins` |
| `is_plasma` | `Plasma proteins` |
| `is_membrane` | `Predicted membrane proteins` |
| `is_intracellular` | `Predicted intracellular proteins` |

`is_secreted_to_blood` comes from the independent `Secretome location` field
containing "blood".

`secretion_tier`:

- **`core`** — `is_secreted AND is_plasma`. Predicted secreted *and* observed as
  a plasma protein, so there is direct evidence the protein reaches circulation.
  This is the tier to use for a blood-borne / endocrine mechanism claim.
- **`extended`** — `is_secreted` only. Secretion is predicted but circulation is
  not evidenced; may act locally or be below plasma detection.
- **`none`** — everything else.

The classes are not mutually exclusive in HPA (a protein can be annotated both
secreted and intracellular, e.g. alternative isoforms or unconventional
secretion), so the booleans are reported independently rather than collapsed.

## Structural ECM flag

Structural extracellular-matrix proteins are secreted but act as scaffold, not
as diffusible signals. **Excluding them is a judgment call and must be disclosed**,
so `flag_structural_ecm` returns an explicit exclusion log alongside the flags.

The rule (also available at runtime as `STRUCTURAL_ECM_RULE_TEXT`) requires
**both**:

1. the HGNC symbol matches one of the documented family regexes in
   `STRUCTURAL_ECM_FAMILIES` — collagens (`COL\d+[A-Z]\d*`), laminins,
   fibrillins, fibulins/EFEMP, elastin+microfibril (ELN/MFAP/EMILIN/LTBP),
   fibronectin, nidogens, tenascins, structural proteoglycans (AGRN, HSPG2,
   BGN, DCN, LUM, VCAN, ACAN, OGN, PRELP, FMOD), cartilage matrilins/COMP; **and**
2. HPA independently supports extracellular localisation — `Protein class`
   contains `Predicted secreted proteins` **or** `Secretome location` contains
   "extracellular matrix".

Requirement 2 prevents symbol homonyms being thrown away. It does real work:
`MFAP1` and `MFAP3` match the elastin-microfibril pattern but HPA gives them no
secreted/ECM support (they are nuclear/intracellular), so they are **kept** and
recorded in the log with `hpa_supports_extracellular=False`.

Keratin-associated proteins (`KRTAP*`) were evaluated and deliberately left out
of the rule: none of the 95 KRTAP symbols in HPA carries secreted/ECM support,
so they are hair-keratin structural proteins rather than matrix.

## Functions

```python
fetch_hpa_annotations(dest_path, columns="g,eg,pc,secl,scml,rnats,blconcms",
                      query="", force=False, timeout=300) -> pd.DataFrame
classify_secretion(hpa_df) -> pd.DataFrame
flag_structural_ecm(genes, hpa_df=None) -> (pd.Series, pd.DataFrame)
annotate_gene_list(genes, hpa_df) -> pd.DataFrame
build_secretome_table(hpa_df) -> (pd.DataFrame, pd.DataFrame)
```

- `fetch_hpa_annotations` caches to `dest_path` and re-reads it on later calls.
- `classify_secretion` returns one row per HPA row: `gene`, `protein_class` (raw
  string), `secretome_location`, `secretome_function`, the four booleans,
  `is_secreted_to_blood`, `secretion_tier`.
- `flag_structural_ecm` returns `(flags, log)`. `flags` is a boolean Series
  indexed by unique input symbol. `log` has one row per family-regex match with
  `matched_family`, `matched_pattern`, `hpa_supports_extracellular`, `excluded`,
  `reason`. **With `hpa_df=None` nothing is flagged** — the rule is deliberately
  not applied on symbol pattern alone.
- `annotate_gene_list` returns exactly one row per input gene, in input order,
  including `found_in_hpa=False` rows for symbols absent from HPA (aliases,
  non-coding, deprecated symbols) so nothing is silently dropped.
- `build_secretome_table` is the genome-wide convenience wrapper.

## Worked example

```python
import kernel as hpa

df = hpa.fetch_hpa_annotations("data/hpa_full.tsv")     # genome-wide, cached
table, ecm_log = hpa.build_secretome_table(df)

table["secretion_tier"].value_counts()
# none 18260, extended 1256, core 646

# annotate a differential-expression hit list
hits = ["NPPB", "NPPA", "GDF15", "COL1A1", "ACTB", "TTN", "IL6"]
ann = hpa.annotate_gene_list(hits, df)
ann.loc[ann.secretion_tier == "core", "gene"].tolist()
# ['NPPB', 'NPPA', 'GDF15', 'COL1A1', 'IL6']   (ACTB, TTN are tier 'none')

# candidate circulating signalling factors: core tier, minus structural matrix
signal = ann[(ann.secretion_tier == "core") & ~ann.is_structural_ecm]
# COL1A1 drops out here (is_structural_ecm=True)

# always report the exclusion log in the methods
ecm_log[ecm_log.excluded].matched_family.value_counts()
```

## Tissue-specificity companion (GTEx)

Deciding whether a secreted protein is *sourced* from a given tissue needs an
expression axis HPA alone does not provide. These helpers add it, and encode two
failure modes that cost real debugging time.

### GTEx API silent no-op

`GTEX_API_GOTCHA` at runtime. The v2 median-expression endpoint requires **both**
`datasetId=gtex_v8` **and** gencode IDs versioned to that release. Omit the
dataset ID, or pass a bare unversioned `ENSG00000123456`, and the API returns
**HTTP 200 with an empty `data` array** and `totalNumberOfItems=0`. There is no
error — the response is indistinguishable from "this gene is not expressed
anywhere", which is a wrong answer that looks like a real one.

`gtex_resolve_gencode_ids` resolves symbols against `gencodeVersion=v26` /
`genomeBuild=GRCh38/hg38`, and `gtex_fetch_tissue_medians` asserts a non-zero row
count rather than returning an empty frame.

### The enrichment-ratio artifact

`RATIO_METRIC_CAVEAT` at runtime. A ratio of `tissue TPM / median across tissues`
is the conventional enrichment metric, and it has a specific failure mode: it
rewards genes that are near-zero in **most** tissues regardless of where their
actual maximum sits. Two real examples from a cardiac-source screen:

- **APOA1** scored a 17-fold "cardiac enrichment" (heart 12.2 TPM vs cross-tissue
  median 0.71) while liver sits at 5,952 TPM — 490-fold higher than heart.
- **PI16** scored 8-fold with heart atrial appendage ranking only **4th**, behind
  ectocervix (179), vagina (140) and tibial nerve (93 TPM).

Both pass a ratio threshold; neither is plausibly cardiac in source. `tissue_specificity`
therefore returns `best_tissue_rank`, `top_tissue`, `top_tissue_tpm` and
`fold_below_top` alongside the ratio, requires a top-`max_rank` placement for
`passes`, and sets `ratio_only_artifact=True` on exactly this pattern. **Never
report a ratio without the rank companion.**

### Base-rate guard

`base_rate_guard(hit_flags, background_flags)` runs the Fisher test that decides
whether a screen found anything. Apply the *same* criterion to the candidate set
and to a matched non-candidate set from the same cohorts. In the screen this
skill was extended for, 11/295 (3.7%) HF-upregulated secreted proteins passed a
cardiac-source test versus 4/86 (4.7%) of non-upregulated ones — OR 0.79,
p=0.75. Eleven named hits with a working positive control, and no enrichment
whatsoever. Run this guard **before** interpreting a candidate list, not after a
reviewer asks.

```python
ids = hpa.gtex_resolve_gencode_ids(candidate_symbols)
med = hpa.gtex_fetch_tissue_medians(list(ids.values()))
mat = hpa.gtex_tissue_matrix(med)
spec = hpa.tissue_specificity(mat, ["Heart_Left_Ventricle", "Heart_Atrial_Appendage"])

spec[spec.ratio_only_artifact]          # inspect before trusting any ratio
bg = hpa.tissue_specificity(mat_background, [...])
hpa.base_rate_guard(spec.passes, bg.passes)
```

Functions: `gtex_resolve_gencode_ids`, `gtex_fetch_tissue_medians`,
`gtex_tissue_matrix`, `tissue_specificity`, `base_rate_guard`.

## Preconditions / postconditions

Preconditions enforced as raising assertions:

- `columns` must include `g`; otherwise the table has no gene key.
- The fetched TSV must contain a `Gene` column — catches an HTML error page
  being cached as if it were data.
- A genome-wide fetch (`query=""`) must return ≥ `MIN_GENOME_ROWS` (10,000)
  rows — catches a truncated or filtered response.
- No empty gene symbols in the fetched table.
- `classify_secretion` requires `Gene` and `Protein class`.

Postconditions asserted:

- `secretion_tier` only ever takes `core` / `extended` / `none`.
- `core` implies `is_secreted`.
- `flag_structural_ecm` never flags a gene absent from its input.
- `annotate_gene_list` returns exactly `len(genes)` rows.

Caveats the caller must handle:

- HPA's genome-wide table contains a handful of duplicated gene symbols (11 in
  the version fetched here); `annotate_gene_list` takes the first match per
  symbol and `flag_structural_ecm` collapses duplicates to one flag.
- `Predicted secreted proteins` is a prediction, not proteomic proof. Only
  `is_plasma` / `is_secreted_to_blood` carry observational weight for
  circulation.
- Absence from HPA is not evidence of non-secretion; check `found_in_hpa`
  before interpreting a `none` tier.

## Standalone use (no agent)

`kernel.py` depends only on the standard library and pandas, and references no
platform APIs. Copy it anywhere and use it directly:

```bash
python - <<'PY'
import sys; sys.path.insert(0, "skills/hpa-secretome")
import kernel as hpa
df = hpa.fetch_hpa_annotations("hpa_full.tsv")
tab, log = hpa.build_secretome_table(df)
tab.to_csv("hpa_secretome_annotation.csv", index=False)
log.to_csv("structural_ecm_exclusion_log.csv", index=False)
print(tab.secretion_tier.value_counts())
PY
```

Requires network access to `www.proteinatlas.org` on first run only; afterwards
the cached TSV is sufficient and the skill is fully offline.
