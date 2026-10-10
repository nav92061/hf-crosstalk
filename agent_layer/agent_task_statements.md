# Sub-agent task statements (verbatim)

The complete task statement dispatched to each of the 23 sub-agent frames, recovered from the platform's frame records. Each sub-agent received this text plus a context summary, and nothing else from the orchestrator's conversation.


---

## Wave 1 — HF signature

`model=claude-opus-5`, `reasoning_effort=high`, status=failed, duration=28.0 min

### Task statement

```

## Shared engineering rules (apply to all your work)

**Environment**: use `environment="vibemed"` on every python/bash call. It has
pandas, numpy, scipy, statsmodels, seaborn, matplotlib, scikit-learn, requests,
pyreadr, rdata. Install anything else with manage_packages(environment="vibemed").

**Machine limits**: 8 GiB RAM, 8 cores, no GPU, no remote compute configured.
Never load a multi-GB file fully into memory; stream it.

**Network**: these domains are granted and reachable (already probed):
ftp.ncbi.nlm.nih.gov, www.ncbi.nlm.nih.gov, www.proteinatlas.org,
raw.githubusercontent.com, api.gdc.cancer.gov, api.figshare.com,
ndownloader.figshare.com, bioconductor.org.
compbio.cn / timer.comp-genomics.org are BLOCKED and NOT grantable - do not try.
If you hit a new block, call request_network_access; never fabricate data as a
workaround and never silently substitute a different source without saying so.

**Skill authoring rules** (critical - these skills must be portable):
- Write the skill as a directory `skills/<skill-name>/` containing `SKILL.md`
  and `kernel.py`. DO NOT publish it - the orchestrator publishes centrally later.
- `kernel.py` MUST NOT reference `host.*` anything. No host.llm, no
  host.artifacts, no save_artifacts. Pure stdlib + pandas/numpy/scipy only.
  Functions take file paths / DataFrames and return DataFrames or dicts.
  This is so the skill runs outside this platform (e.g. under a different model
  or with no agent at all).
- `kernel.py` must contain only functions, imports and literal constants at top
  level - no top-level classes, no decorators (the loader requires this).
- Guardrails belong in CODE as assertions that RAISE, not as advice in prose.
- `SKILL.md` starts with YAML frontmatter (`name:`, `description:`), then
  purpose, the exact function signatures, a worked example, preconditions/
  postconditions, and a "Standalone use (no agent)" section.
- Docstring every public function with argument types and return shape.

**Scientific integrity - this is the whole point of the project**:
- Every number you report must come from code you actually ran on real data.
- NEVER invent, estimate, or carry over a number from any prior draft.
- There WAS a previous PAPER.md in this project full of numbers with no backing
  data. It has been DISCARDED. Do not look for it, do not reuse its values
  (e.g. "528 genes", "OPRK1 1.00 TPM", specific cell line names). Derive
  everything fresh.
- Report null, weak and contradictory results plainly. A failed control is a
  finding, not something to hide or work around.
- If something cannot be done with available data, say so explicitly in your
  structured output rather than approximating it.

**Artifacts**: save real output files with `save_artifacts(files=[...],
language="python")`. Use the EXACT filenames requested. Report every
artifact filename + version_id in your structured output.

# Your task: the `geo-bulk-de` skill + the human HF signature

## Step 1 - Author the `geo-bulk-de` skill (plan step: "Author geo-bulk-de skill")

Create `skills/geo-bulk-de/SKILL.md` + `skills/geo-bulk-de/kernel.py`: a
reusable skill for differential expression from a GEO series.

`kernel.py` should provide roughly:
- `fetch_geo_supplementary(gse, dest_dir)` -> list of downloaded paths. Build the
  FTP path from the GSE number (e.g. GSE116250 -> series/GSE116nnn/GSE116250/suppl/).
- `fetch_series_matrix(gse, dest_dir)` -> path to the series_matrix.txt.gz.
- `parse_series_matrix_metadata(path)` -> DataFrame of per-sample characteristics
  (GSM id, title, and all `!Sample_characteristics_ch1` fields parsed into columns).
  This is how you get the disease group labels - parse them, do not hardcode.
- `load_expression(path, ...)` -> genes x samples DataFrame.
- `run_de(expr, groups, design="limma_trend"|"voom_like", ...)` -> DataFrame with
  columns gene, log2FC, t, p_value, fdr, mean_expr, n_case, n_ctrl.
  Implement in Python/scipy/statsmodels (no R dependency): for counts use a
  log-CPM + variance-trend moderated t (a limma-trend equivalent); for RPKM/FPKM
  use log2(x+1) then moderated t with empirical-Bayes variance shrinkage toward
  the mean-variance trend. Implement the eBayes shrinkage explicitly and document
  the formula you used in SKILL.md - be honest that this is a limma-trend
  equivalent, not a call into limma itself.
- `replication_stats(de_a, de_b, fdr=0.05)` -> dict with n_shared_genes,
  n_sig_a, n_sig_b, n_replicated, hypergeometric p for the overlap,
  directional_concordance (fraction agreeing in sign among genes significant in
  both), and Spearman rho of log2FC on shared genes.

PRECONDITIONS that raise: sample IDs unique; >=2 samples per group; expression
matrix non-empty with no all-NA genes; group vector aligned to matrix columns.

## Step 2 - Derive the GSE116250 signature (plan step: "Derive GSE116250 HF signature")

GSE116250: LV RNA-seq, non-failing vs dilated (DCM) vs ischemic (ICM)
cardiomyopathy. Supplementary file is
`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116250/suppl/GSE116250_rpkm.txt.gz`
- NOTE this is RPKM, not counts. Use the RPKM-appropriate path and state the
limitation. Get group labels by parsing the series matrix.
Primary contrast: failing (DCM + ICM pooled) vs non-failing. ALSO run DCM-only
and ICM-only vs non-failing as secondary contrasts.
Save `hf_de_GSE116250.csv` (full table, all genes).

## Step 3 - Validate in GSE141910 (plan step: "Validate in GSE141910")

GSE141910 is the larger cohort (MAGNet). Supplementary is `GSE141910_RAW.tar`
(~60 MB) containing one `GSM*_C*.csv.gz` per sample - these are COUNTS, so use
the counts path. Download the tar, extract, and assemble the count matrix
(inspect one CSV first to learn its column layout before writing the loop).
Group labels come from the series matrix characteristics (expect
Non-Failing Donor / Dilated cardiomyopathy / Hypertrophic / Peripartal etc).
Primary contrast: failing (DCM + HCM, i.e. all cardiomyopathy) vs non-failing
donor. Save `hf_de_GSE141910.csv`.
Then compute replication statistics between the two cohorts with
`replication_stats`. Be careful: gene identifiers may differ between cohorts
(symbols vs Ensembl) - harmonize explicitly and report how many genes mapped.

## Step 4 - Save artifacts (plan step: "Save signature artifacts")

Save with these exact filenames:
- `hf_de_GSE116250.csv`, `hf_de_GSE141910.csv` - full DE tables
- `hf_signature_replicated.csv` - genes significant (FDR<0.05) and same-sign in
  BOTH cohorts. Columns: gene, log2FC_discovery, fdr_discovery, log2FC_validation,
  fdr_validation, direction. Include BOTH up and down genes with a `direction`
  column (downstream needs the up set, but the down set is needed for controls).
- `hf_replication_stats.json` - the replication_stats dict plus per-cohort sample
  counts per group, and the gene-ID harmonization counts.
- the two skill files.

Call update_step_status for each of the four plan steps as you go, using the
exact titles given above in parentheses.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "n_up_replicated": {
      "type": "integer",
      "description": "HF-UPregulated genes replicated in both cohorts (FDR<0.05 both, positive sign both)"
    },
    "n_down_replicated": {
      "type": "integer"
    },
    "replication": {
      "type": "object",
      "description": "the replication_stats dict: overlap p, directional concordance, spearman rho, n_shared_genes"
    },
    "cohort_samples": {
      "type": "object",
      "description": "per-cohort per-group sample counts actually used"
    },
    "top_up_genes": {
      "type": "array",
      "items": {
        "type": "string"
      },
      "description": "top 40 HF-upregulated replicated genes by combined evidence"
    },
    "gene_id_harmonization": {
      "type": "string",
      "description": "how IDs were harmonized and how many genes mapped"
    },
    "method_caveats": {
      "type": "array",
      "items": {
        "type": "string"
      },
      "description": "honest limitations: RPKM vs counts, eBayes reimplementation, group pooling"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object"
      }
    },
    "skill_path": {
      "type": "string"
    }
  },
  "required": [
    "n_up_replicated",
    "replication",
    "cohort_samples",
    "top_up_genes",
    "artifacts",
    "method_caveats"
  ]
}
```


---

## Wave 1 — Secretome and ligands

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=28.4 min

### Task statement

```

## Shared engineering rules (apply to all your work)

**Environment**: use `environment="vibemed"` on every python/bash call. It has
pandas, numpy, scipy, statsmodels, seaborn, matplotlib, scikit-learn, requests,
pyreadr, rdata. Install anything else with manage_packages(environment="vibemed").

**Machine limits**: 8 GiB RAM, 8 cores, no GPU, no remote compute configured.
Never load a multi-GB file fully into memory; stream it.

**Network**: these domains are granted and reachable (already probed):
ftp.ncbi.nlm.nih.gov, www.ncbi.nlm.nih.gov, www.proteinatlas.org,
raw.githubusercontent.com, api.gdc.cancer.gov, api.figshare.com,
ndownloader.figshare.com, bioconductor.org.
compbio.cn / timer.comp-genomics.org are BLOCKED and NOT grantable - do not try.
If you hit a new block, call request_network_access; never fabricate data as a
workaround and never silently substitute a different source without saying so.

**Skill authoring rules** (critical - these skills must be portable):
- Write the skill as a directory `skills/<skill-name>/` containing `SKILL.md`
  and `kernel.py`. DO NOT publish it - the orchestrator publishes centrally later.
- `kernel.py` MUST NOT reference `host.*` anything. No host.llm, no
  host.artifacts, no save_artifacts. Pure stdlib + pandas/numpy/scipy only.
  Functions take file paths / DataFrames and return DataFrames or dicts.
  This is so the skill runs outside this platform (e.g. under a different model
  or with no agent at all).
- `kernel.py` must contain only functions, imports and literal constants at top
  level - no top-level classes, no decorators (the loader requires this).
- Guardrails belong in CODE as assertions that RAISE, not as advice in prose.
- `SKILL.md` starts with YAML frontmatter (`name:`, `description:`), then
  purpose, the exact function signatures, a worked example, preconditions/
  postconditions, and a "Standalone use (no agent)" section.
- Docstring every public function with argument types and return shape.

**Scientific integrity - this is the whole point of the project**:
- Every number you report must come from code you actually ran on real data.
- NEVER invent, estimate, or carry over a number from any prior draft.
- There WAS a previous PAPER.md in this project full of numbers with no backing
  data. It has been DISCARDED. Do not look for it, do not reuse its values
  (e.g. "528 genes", "OPRK1 1.00 TPM", specific cell line names). Derive
  everything fresh.
- Report null, weak and contradictory results plainly. A failed control is a
  finding, not something to hide or work around.
- If something cannot be done with available data, say so explicitly in your
  structured output rather than approximating it.

**Artifacts**: save real output files with `save_artifacts(files=[...],
language="python")`. Use the EXACT filenames requested. Report every
artifact filename + version_id in your structured output.

# Your task: the `hpa-secretome` and `cellchat-lr` skills + their tables

## Step 1 - Author `hpa-secretome` (plan step: "Author hpa-secretome skill")

Create `skills/hpa-secretome/SKILL.md` + `kernel.py`.
Human Protein Atlas has a working download API (verified):
`https://www.proteinatlas.org/api/search_download.php?search=<query>&format=tsv&columns=g,gs,pc,scl&compress=no`
Column codes: g=Gene, gs=gene synonym, pc=Protein class, scl=subcellular location.
Investigate what other useful columns exist (e.g. secretome location/annotation,
tissue specificity) and use them if available - inspect the returned header rather
than assuming.

`kernel.py` should provide:
- `fetch_hpa_annotations(dest_path, columns=..., query=...)` -> DataFrame keyed on
  gene symbol with the raw protein-class string. Cache to disk; re-read if present.
- `classify_secretion(hpa_df)` -> DataFrame with boolean/tier columns derived from
  the Protein class field: `is_secreted` (predicted secreted), `is_plasma`
  (plasma protein), `is_membrane`, `is_intracellular`, and a `secretion_tier`
  of "core" (secreted AND plasma-detectable), "extended" (secreted only),
  or "none". Parse the class strings - do not hardcode gene lists.
- `flag_structural_ecm(genes)` -> boolean Series flagging structural
  extracellular-matrix genes (collagens COL*, elastin, fibrillins, laminins,
  etc.) which are secreted but are structural rather than signaling. Use a
  documented, rule-based definition (gene-symbol families + HPA class) and
  return an explicit exclusion log, since excluding these is a judgment call
  the paper must disclose.
- `annotate_gene_list(genes, hpa_df)` -> tidy DataFrame, one row per input gene,
  with all of the above plus whether the gene was found in HPA at all.

PRECONDITION that raises: HPA fetch returned a parseable TSV with a Gene column
and >10000 rows (guards against silently saving an error page).

## Step 2 - Author `cellchat-lr` (plan step: "Author cellchat-lr skill")

Create `skills/cellchat-lr/SKILL.md` + `kernel.py`.
CellChatDB human is an R .rda (verified reachable, 756 KB):
`https://raw.githubusercontent.com/sqjin/CellChat/master/data/CellChatDB.human.rda`
Parse it in Python with `rdata` or `pyreadr` (both installed). The object is a
list with elements `interaction`, `complex`, `cofactor`, `geneInfo`. You will
need to inspect the parsed structure carefully - the `interaction` table has
columns like interaction_name, pathway_name, ligand, receptor, annotation,
evidence. The `complex` table maps a complex name to subunit_1..subunit_n.

`kernel.py` should provide:
- `fetch_cellchatdb(dest_path)` -> path (cache-aware).
- `load_cellchatdb(path)` -> dict of tidy DataFrames: interaction, complex,
  cofactor, geneInfo.
- `expand_receptor_complexes(interaction_df, complex_df)` -> long DataFrame with
  one row per (interaction, receptor_subunit_gene), plus `n_subunits` and
  `is_complex`. A receptor that is not in the complex table is its own single
  subunit. Same treatment for multi-subunit ligands.
- `map_ligands_to_receptors(ligand_genes, db)` -> DataFrame of interactions whose
  ligand is in the input list, with receptor subunits expanded, pathway_name,
  annotation (Secreted Signaling / ECM-Receptor / Cell-Cell Contact) and
  evidence. Include an `annotation` filter argument defaulting to
  "Secreted Signaling" since only secreted signaling is mechanistically
  plausible at a distance - but keep the others retrievable and document why.

PRECONDITION that raises: the parsed interaction table has >1000 rows and the
required columns are present.

## Step 3 - Emit the tables (plan step: "Emit secretome and L-R tables")

Save with these exact filenames:
- `hpa_secretome_annotation.csv` - genome-wide: gene, protein_class raw string,
  is_secreted, is_plasma, secretion_tier, is_structural_ecm.
- `cellchatdb_interactions.csv` - full parsed interaction table.
- `cellchatdb_complexes.csv` - the complex -> subunits mapping, tidied long.
- `cellchatdb_receptor_subunits.csv` - the fully expanded per-subunit table for
  ALL interactions (not filtered to any ligand list), so the downstream scoring
  track can join to it directly.
- `secretome_lr_resource_summary.json` - counts: HPA rows, n secreted, n plasma,
  n core tier, n structural-ECM flagged; CellChatDB n interactions, n pathways,
  n complexes, n multi-subunit receptors, breakdown by annotation category.
- the four skill files.

Do NOT filter to the HF gene list - that join happens downstream and the HF list
is being produced in parallel. Your job is the complete, reusable resource.

Call update_step_status for each of the three plan steps using the exact titles
in parentheses above.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "hpa": {
      "type": "object",
      "description": "row counts: total genes, n_secreted, n_plasma, n_core_tier, n_structural_ecm"
    },
    "cellchatdb": {
      "type": "object",
      "description": "n_interactions, n_pathways, n_complexes, n_multisubunit_receptors, counts by annotation category"
    },
    "hpa_columns_available": {
      "type": "array",
      "items": {
        "type": "string"
      },
      "description": "what HPA columns you actually retrieved"
    },
    "structural_ecm_rule": {
      "type": "string",
      "description": "the exact rule used to flag structural ECM genes"
    },
    "parsing_notes": {
      "type": "string",
      "description": "how the .rda was parsed and any surprises in its structure"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object"
      }
    },
    "skill_paths": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "hpa",
    "cellchatdb",
    "structural_ecm_rule",
    "artifacts",
    "skill_paths"
  ]
}
```


---

## Wave 1 — Tumor axis

`model=claude-opus-5`, `reasoning_effort=high`, status=failed, duration=32.9 min

### Task statement

```

## Shared engineering rules (apply to all your work)

**Environment**: use `environment="vibemed"` on every python/bash call. It has
pandas, numpy, scipy, statsmodels, seaborn, matplotlib, scikit-learn, requests,
pyreadr, rdata. Install anything else with manage_packages(environment="vibemed").

**Machine limits**: 8 GiB RAM, 8 cores, no GPU, no remote compute configured.
Never load a multi-GB file fully into memory; stream it.

**Network**: these domains are granted and reachable (already probed):
ftp.ncbi.nlm.nih.gov, www.ncbi.nlm.nih.gov, www.proteinatlas.org,
raw.githubusercontent.com, api.gdc.cancer.gov, api.figshare.com,
ndownloader.figshare.com, bioconductor.org.
compbio.cn / timer.comp-genomics.org are BLOCKED and NOT grantable - do not try.
If you hit a new block, call request_network_access; never fabricate data as a
workaround and never silently substitute a different source without saying so.

**Skill authoring rules** (critical - these skills must be portable):
- Write the skill as a directory `skills/<skill-name>/` containing `SKILL.md`
  and `kernel.py`. DO NOT publish it - the orchestrator publishes centrally later.
- `kernel.py` MUST NOT reference `host.*` anything. No host.llm, no
  host.artifacts, no save_artifacts. Pure stdlib + pandas/numpy/scipy only.
  Functions take file paths / DataFrames and return DataFrames or dicts.
  This is so the skill runs outside this platform (e.g. under a different model
  or with no agent at all).
- `kernel.py` must contain only functions, imports and literal constants at top
  level - no top-level classes, no decorators (the loader requires this).
- Guardrails belong in CODE as assertions that RAISE, not as advice in prose.
- `SKILL.md` starts with YAML frontmatter (`name:`, `description:`), then
  purpose, the exact function signatures, a worked example, preconditions/
  postconditions, and a "Standalone use (no agent)" section.
- Docstring every public function with argument types and return shape.

**Scientific integrity - this is the whole point of the project**:
- Every number you report must come from code you actually ran on real data.
- NEVER invent, estimate, or carry over a number from any prior draft.
- There WAS a previous PAPER.md in this project full of numbers with no backing
  data. It has been DISCARDED. Do not look for it, do not reuse its values
  (e.g. "528 genes", "OPRK1 1.00 TPM", specific cell line names). Derive
  everything fresh.
- Report null, weak and contradictory results plainly. A failed control is a
  finding, not something to hide or work around.
- If something cannot be done with available data, say so explicitly in your
  structured output rather than approximating it.

**Artifacts**: save real output files with `save_artifacts(files=[...],
language="python")`. Use the EXACT filenames requested. Report every
artifact filename + version_id in your structured output.

# Your task: the `tcga-pancan` skill + tumor expression and TME tables

The TIMER2 web server is blocked and not grantable. Use the GDC-hosted Thorsson
immune-landscape table instead (it is the published source for TCGA
immune/stromal state and covers the same quantities) and say so explicitly.

## Step 1 - Author `tcga-pancan` (plan step: "Author tcga-pancan skill")

Create `skills/tcga-pancan/SKILL.md` + `kernel.py`.

These GDC endpoints are VERIFIED reachable and need no token:
- PanCanAtlas EBPlusPlus RNA-seq (batch-corrected, ~2 GB uncompressed TSV):
  `https://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611`
  Layout: first column `gene_id` formatted `SYMBOL|ENTREZ`, remaining columns are
  TCGA ALIQUOT barcodes (e.g. TCGA-OR-A5J1-01A-11R-A29S-07). Values are
  RSEM normalized counts (NOT TPM - you must state this; do not call them TPM).
- Thorsson immune landscape: `https://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf`
  Has patient_barcode, aliquot_barcode, `cancer type`, immune subtype, leukocyte
  fraction, stromal fraction and related scores.
- TCGA-CDR clinical (xlsx): `https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81`

**Memory discipline is mandatory** - 8 GiB RAM total. `kernel.py` must provide:
- `stream_expression_rows(url_or_path, wanted_genes, dest_csv, chunk_bytes=...)`:
  stream the TSV line by line (gzip/plain, straight from HTTP or a cached local
  file), keep only rows whose symbol is in `wanted_genes`, and write a tidy
  gene x sample CSV. NEVER read the whole matrix into memory. Report how many
  requested genes were found vs missing. Cache the download to disk first if
  that is more robust, but still parse it streaming.
- `parse_tcga_barcode(barcode)` -> dict with patient, sample_type_code,
  is_primary_tumor. Sample type code is barcode positions 13-14 ("01" primary
  solid tumor, "11" solid normal, "06" metastatic, "03" peripheral blood).
- `build_sample_map(thorsson_df)` -> DataFrame mapping aliquot barcode ->
  tumor_type (TCGA study abbreviation) -> and a full human-readable tumor name.
  Provide a GDC `cases` API fallback (`https://api.gdc.cancer.gov/cases`, paginated,
  fields=submitter_id,project.project_id) for barcodes Thorsson does not cover.
- `filter_solid_primary(expr, sample_map, exclude_studies=("LAML","DLBC","THYM"))`
  -> expression restricted to primary solid tumors, heme/thymoma excluded.
  Document why: this project asks about solid tumors receiving cardiac signals.
- `summarize_by_tumor_type(expr, sample_map, stat="median")` -> gene x tumor_type
  matrix, plus a per-tumor-type n table.
- `fetch_thorsson_immune(dest)` / `fetch_tcga_cdr(dest)` loaders.

PRECONDITIONS that raise: requested gene list non-empty; the streamed matrix
yielded >0 matching rows; every retained sample has a tumor-type assignment;
no duplicate aliquot barcodes.

## Step 2 - Extract expression (plan step: "Extract receptor and reference expression")

Build the wanted-gene superset:
(a) ALL receptor subunit genes in CellChatDB - fetch and parse
    `https://raw.githubusercontent.com/sqjin/CellChat/master/data/CellChatDB.human.rda`
    yourself with the `rdata` or `pyreadr` package (do not wait for the parallel
    track; a rough parse to get the union of receptor genes is enough here),
(b) all CellChatDB ligand genes too,
(c) housekeeping reference genes: ACTB, GAPDH, TUBB, RPLP0, B2M, PGK1, PPIA,
    TBP, RPL13A, SDHA, UBC, YWHAZ, HPRT1, EEF1A1,
(d) canonical pathway target genes for downstream activity scoring: AXIN2, LEF1,
    TCF7, CCND1, MYC, NKD1, RNF43, SP5, NOTUM (Wnt); SERPINE1, SMAD7, JUNB,
    SKIL (TGFb/GDF); NPPA, NPPB, NR4A1, FOS, EGR1 (cAMP/PKG-ish);
    plus CD8A, CD68, PTPRC, PECAM1, COL1A1, ACTA2, VIM, EPCAM as
    compartment markers so bulk-tissue composition can be assessed.
Deduplicate. Then stream-extract those rows.

Also extract a set of ~200 randomly chosen genes (fixed seed, recorded) spanning
the full expression range, to serve as the background pool for the downstream
matched-random specificity null. Save the seed and the gene list.

## Step 3 - TME scores (plan step: "Fetch tumor microenvironment scores")

Download the Thorsson table; extract per-sample leukocyte fraction, stromal
fraction, immune subtype, and any TGF-beta / wound-healing / proliferation
signature scores present. Aggregate to per-tumor-type medians too.
Inspect the real column names before selecting - report what you actually found.

## Step 4 - Save (plan step: "Save tumor expression artifacts")

Exact filenames:
- `tcga_expression_subset.csv.gz` - gene x sample, primary solid tumors only,
  gzipped (this is the big one; keep it compressed).
- `tcga_expr_by_tumor_type.csv` - gene x tumor_type medians.
- `tcga_sample_map.csv` - aliquot barcode, patient, tumor_type, tumor_name, n.
- `tcga_tme_scores.csv` - per-sample TME/immune scores.
- `tcga_tme_by_tumor_type.csv` - per-tumor-type medians.
- `tcga_background_genes.json` - the random background gene list + seed.
- `tcga_extraction_summary.json` - n samples, n tumor types with n each, n genes
  requested/found/missing (list the missing ones), the expression unit, and the
  housekeeping gene median values per tumor type (needed to calibrate the floor).
- the two skill files.

If `tcga_expression_subset.csv.gz` exceeds ~500 MB, note it and pass
destination={"tcga_expression_subset.csv.gz": "working_data"} to save_artifacts.

Call update_step_status for each of the four plan steps using the exact titles
in parentheses above.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "n_samples": {
      "type": "integer"
    },
    "n_tumor_types": {
      "type": "integer"
    },
    "tumor_type_counts": {
      "type": "object"
    },
    "expression_unit": {
      "type": "string",
      "description": "the true unit of the values - do not call RSEM values TPM"
    },
    "genes_requested": {
      "type": "integer"
    },
    "genes_found": {
      "type": "integer"
    },
    "genes_missing": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "housekeeping_medians": {
      "type": "object",
      "description": "median value of each housekeeping gene across all samples, for floor calibration"
    },
    "tme_columns_found": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "tme_substitution_note": {
      "type": "string",
      "description": "explicit statement that Thorsson replaced blocked TIMER2 and what that changes"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object"
      }
    },
    "skill_paths": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "n_samples",
    "n_tumor_types",
    "expression_unit",
    "genes_found",
    "housekeeping_medians",
    "artifacts",
    "skill_paths"
  ]
}
```


---

## Wave 1 — Literature base

`model=claude-opus-5`, `reasoning_effort=high`, status=failed, duration=18.6 min

### Task statement

```

## Shared engineering rules (apply to all your work)

**Environment**: use `environment="vibemed"` on every python/bash call. It has
pandas, numpy, scipy, statsmodels, seaborn, matplotlib, scikit-learn, requests,
pyreadr, rdata. Install anything else with manage_packages(environment="vibemed").

**Machine limits**: 8 GiB RAM, 8 cores, no GPU, no remote compute configured.
Never load a multi-GB file fully into memory; stream it.

**Network**: these domains are granted and reachable (already probed):
ftp.ncbi.nlm.nih.gov, www.ncbi.nlm.nih.gov, www.proteinatlas.org,
raw.githubusercontent.com, api.gdc.cancer.gov, api.figshare.com,
ndownloader.figshare.com, bioconductor.org.
compbio.cn / timer.comp-genomics.org are BLOCKED and NOT grantable - do not try.
If you hit a new block, call request_network_access; never fabricate data as a
workaround and never silently substitute a different source without saying so.

**Skill authoring rules** (critical - these skills must be portable):
- Write the skill as a directory `skills/<skill-name>/` containing `SKILL.md`
  and `kernel.py`. DO NOT publish it - the orchestrator publishes centrally later.
- `kernel.py` MUST NOT reference `host.*` anything. No host.llm, no
  host.artifacts, no save_artifacts. Pure stdlib + pandas/numpy/scipy only.
  Functions take file paths / DataFrames and return DataFrames or dicts.
  This is so the skill runs outside this platform (e.g. under a different model
  or with no agent at all).
- `kernel.py` must contain only functions, imports and literal constants at top
  level - no top-level classes, no decorators (the loader requires this).
- Guardrails belong in CODE as assertions that RAISE, not as advice in prose.
- `SKILL.md` starts with YAML frontmatter (`name:`, `description:`), then
  purpose, the exact function signatures, a worked example, preconditions/
  postconditions, and a "Standalone use (no agent)" section.
- Docstring every public function with argument types and return shape.

**Scientific integrity - this is the whole point of the project**:
- Every number you report must come from code you actually ran on real data.
- NEVER invent, estimate, or carry over a number from any prior draft.
- There WAS a previous PAPER.md in this project full of numbers with no backing
  data. It has been DISCARDED. Do not look for it, do not reuse its values
  (e.g. "528 genes", "OPRK1 1.00 TPM", specific cell line names). Derive
  everything fresh.
- Report null, weak and contradictory results plainly. A failed control is a
  finding, not something to hide or work around.
- If something cannot be done with available data, say so explicitly in your
  structured output rather than approximating it.

**Artifacts**: save real output files with `save_artifacts(files=[...],
language="python")`. Use the EXACT filenames requested. Report every
artifact filename + version_id in your structured output.

# Your task: the literature evidence base

Use the PubMed connector from the `repl` tool: `host.mcp("pubmed", <method>, ...)`.
First run `search_skills(query="pubmed literature search")` and load the
`mcp-pubmed` skill doc to get exact method names and signatures. Related useful
connectors: `mcp-biorxiv` for preprints. Do the MCP calls in `repl` cells and
write results to `handoff/*.json`, then analyze in a `python` cell
(environment="vibemed"). MCP calls do not work in the python kernel.

**Every citation must be real.** For each, record PMID, first author, year,
journal, title, and the specific finding. Never write a PMID you did not
retrieve from the connector. If you cannot verify a claim, mark it
`unsupported` rather than attaching an approximate citation. Spot-check a few
PMIDs by re-fetching them and confirming the title matches.

## Step 1 - Reverse cardio-oncology (plan step: "Survey reverse cardio-oncology literature")

Search and synthesize on:
- heart failure promoting tumour growth / tumourigenesis (the foundational
  experimental work - e.g. HF mouse models with implanted tumours, and the
  identification of specific HF-secreted factors driving tumour growth)
- epidemiology: cancer incidence in heart-failure cohorts, and the confounding
  problem (shared risk factors, surveillance bias, competing risk of death)
- the reverse direction and shared-substrate view: inflammation, clonal
  haematopoiesis (CHIP) linking cardiovascular disease and cancer
- "reverse cardio-oncology" as a named field, reviews defining it
For each, note study design and whether the evidence is mechanistic, associative,
or a review. Explicitly capture the strongest arguments AGAINST a causal HF->cancer
link.

## Step 2 - Candidate ligand and receptor biology (plan step: "Survey candidate ligand biology")

The HF-upregulated secreted ligand list is being computed in parallel, so work
from the well-established HF secretome / cardiac stress genes and cover the
canonical candidates: NPPA, NPPB (natriuretic peptides) and NPR1/NPR2/NPR3;
PENK/proenkephalin and opioid receptors OPRK1/OPRD1/OPRM1; POSTN (periostin);
SFRP/WNT family members and FZD/LRP receptors; GDF15; SERPINA3; CTGF/CCN2;
FSTL1; THBS4; COMP; galectin-3 (LGALS3); IGFBP7; CXCL/CCL chemokines.
For each: is it established as elevated in human HF (plasma or myocardium), and
what is known about its receptor in cancer? Record the evidence type
(transcript / protein / functional / clinical) - this distinction matters
because the paper's analysis is transcript-only.
Also search specifically for whether opioid receptors are expressed in solid
tumours, and whether FZD2 is expressed/functional in gynaecologic cancers -
these are likely to be axes the analysis nominates.

## Step 3 - Emit the evidence table (plan step: "Emit citation evidence table")

Exact filenames:
- `literature_citations.csv` - one row per citation: pmid, first_author, year,
  journal, title, topic (reverse_cardio_oncology / epidemiology / ligand_biology /
  receptor_cancer / counter_evidence / methods), claim_supported (one sentence),
  evidence_type (mechanistic_animal / human_observational / human_interventional /
  in_vitro / review / meta_analysis), evidence_strength (strong/moderate/weak),
  direction (supports / contradicts / neutral).
- `ligand_evidence.csv` - one row per ligand: gene, established_in_human_HF
  (yes/no/limited), hf_evidence_type, best_pmids (semicolon separated),
  receptor(s), receptor_cancer_evidence, notes.
- `literature_synthesis.md` - prose suitable for the paper: an introduction
  section on what is and is not established about HF->cancer, a paragraph on
  the confounding problem, a paragraph on why a ligand-receptor compatibility
  design is a reasonable but limited way to generate hypotheses, and a
  discussion section on the leading candidate axes. Cite inline as
  (Author Year, PMID: xxxxxxx). Be explicit about where evidence is thin.
- `literature_gaps.json` - claims the paper might want to make that the
  literature does NOT support, so the writing track avoids them.

Aim for 40-70 verified citations - depth and verification over volume.

Call update_step_status for each of the three plan steps using the exact titles
in parentheses above.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "n_citations": {
      "type": "integer"
    },
    "n_verified": {
      "type": "integer",
      "description": "how many PMIDs you re-fetched and confirmed"
    },
    "key_foundational_pmids": {
      "type": "array",
      "items": {
        "type": "object"
      },
      "description": "the handful of must-cite papers with pmid+finding"
    },
    "counter_evidence_summary": {
      "type": "string",
      "description": "strongest arguments against a causal HF->cancer link"
    },
    "ligands_established_in_hf": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "ligands_weak_evidence": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "opioid_receptor_tumor_evidence": {
      "type": "string",
      "description": "what the literature actually says about opioid receptor expression in solid tumours"
    },
    "unsupported_claims": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object"
      }
    }
  },
  "required": [
    "n_citations",
    "n_verified",
    "key_foundational_pmids",
    "counter_evidence_summary",
    "unsupported_claims",
    "artifacts"
  ]
}
```


---

## Wave 1 — DepMap

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=46.2 min

### Task statement

```

## Shared engineering rules

**Environment**: use `environment="vibemed"` on every python/bash call (has pandas,
numpy, scipy, statsmodels, seaborn, matplotlib, scikit-learn, requests, pyreadr, rdata).

**Machine**: 8 GiB RAM, 8 cores, no GPU, no remote compute. Stream large files;
never load multi-GB into memory.

**Network granted**: ftp.ncbi.nlm.nih.gov, www.ncbi.nlm.nih.gov, www.proteinatlas.org,
raw.githubusercontent.com, api.gdc.cancer.gov, api.figshare.com,
ndownloader.figshare.com, bioconductor.org, depmap.org (NOTE: depmap.org itself
serves a bot-verification interstitial - use figshare for actual DepMap files).
compbio.cn / timer.comp-genomics.org are BLOCKED and NOT grantable.
If blocked on something needed, call request_network_access. NEVER fabricate data.

**Skill authoring**: write `skills/<name>/SKILL.md` + `skills/<name>/kernel.py`.
Do NOT publish (orchestrator publishes centrally). `kernel.py` MUST NOT use
`host.*` - pure stdlib + pandas/numpy/scipy so it runs outside this platform.
Only functions/imports/literal constants at top level (no top-level classes or
decorators). Guardrails = assertions that RAISE, not prose. SKILL.md needs YAML
frontmatter (name, description), signatures, worked example, pre/postconditions,
and a "Standalone use (no agent)" section.

**Scientific integrity**: every number must come from code you ran on real data.
Never invent or carry over a number. A previous PAPER.md in this project was full
of unbacked numbers and has been DISCARDED - do not reuse its values or gene/cell
line names. Report null and contradictory results plainly.

**Artifacts**: save with `save_artifacts(files=[...], language="python")` using
EXACT filenames requested; report filename + version_id in structured output.

## What the project has established so far (all from real data, verified)

Hypothesis: heart failure may promote cancer via secreted factors; susceptibility
depends on ligand availability (source side) AND receiver capacity (tumor side).

Phase 0-1 results you should build on:
- HF signature from GSE116250 (discovery, n=50 failing vs 14 non-failing) and
  GSE141910 (validation, MAGNet). 3084 HF-UPregulated genes replicated in both
  cohorts (FDR<0.05, same sign). Directional concordance 0.881, log2FC Spearman
  0.541. Label-permutation control gave 0 significant genes.
- 295 of those are HPA-annotated secreted; 72 are CellChatDB ligands.
- 188 Secreted-Signaling interactions from 51 HF-up secreted ligands, 38 pathways.
- TCGA: 9538 primary solid tumor samples, 30 tumor types (LAML/DLBC/THYM excluded).
  Expression unit is RSEM normalized_count (EBPlusPlus batch-corrected),
  **NOT TPM** - GAPDH median is ~63,656. Any floor must be relative to
  housekeeping, never an absolute "2 TPM".
- Expression floor = 0.005 x GAPDH median = 345.67 RSEM units. It drops 28 of 188
  interactions, INCLUDING all three PENK->opioid receptor interactions
  (PENK_OPRK1, PENK_OPRD1, PENK_OPRM1). OPRD1 max across tumor types is 4e-5 of
  GAPDH; OPRK1 1.4e-3. The opioid axis is dead on receptor expression grounds.
- Top ligands by availability: NPPA (3.05), NPPB (3.03), PENK (2.85), WNT9A (2.32),
  CCL5 (1.59), CXCL9 (1.34).
- Primary ranking (floor-filtered) top 6: UCS 0.679 (95% CI 0.400-0.914),
  OV 0.432, KIRC 0.377, THCA 0.285, PRAD 0.272, PAAD 0.244.
  Double-centered ranking reorders to: ACC, OV, UCEC, PRAD, PCPG, UCS.
  Primary vs double-centered Spearman = 0.318 (correction-sensitive!).
- KEY METHODOLOGICAL FINDING: the matched-random specificity null PASSES when the
  floor is applied (rho_mean 0.107, 1.7% of random sets exceed rho 0.5) but FAILS
  without it (rho_mean 0.363, 21.5%). An unmatched/naive random pool gives 0.323.
  So the "global expression axis confound" is substantially CAUSED by scoring
  untranscribed receptors, and the null pool must be expression-matched.
- Receiver capacity vs tumor microenvironment: only 1 of 29 Thorsson TME features
  survives FDR (B-cell receptor score rho=0.57, FDR=0.028). The "exploit /
  microenvironment" half of the hypothesis is largely UNSUPPORTED.

Useful artifact version_ids (use host.artifact_path(vid) or {{artifact:vid}}):
- hf_signature_replicated.csv = 76685c52-1489-4bd0-9e51-60c282afbf19
  (IMPORTANT: join on the `gene_symbol` column, NOT `gene` which is Ensembl)
- hpa_secretome_annotation.csv = efb93594-4a96-4a02-8bab-e3030adece36
- cellchatdb_receptor_subunits.csv = 4aabe49f-60f3-4179-a333-ba2cc7c548a0
- tcga_expr_by_tumor_type.csv = 07e2a94c-6093-4cfb-b343-56f042647e20
- crosstalk_network.csv = 8086ab6e-ce8d-466a-804e-31b39df5e0e1
- pancancer_susceptibility_ranking.csv = c8e165e9-4061-43cb-b270-af809942d9d6
- twofactor_decomposition.csv = 23bff2cf-a7f4-42d8-b520-ffcf9b279745
- literature_citations.csv = 940dda9f-ce28-4f2a-a1b6-bcc9d5cc26a8
- ligand_evidence.csv = 8f8d6101-f192-48d8-80ba-10e96f46e295

# Your task: the `depmap-lineage` skill + dependency testing + cell-line panel

## Step 1 - Author `depmap-lineage` (plan step: "Author depmap-lineage skill")

Create `skills/depmap-lineage/SKILL.md` + `skills/depmap-lineage/kernel.py`.

depmap.org serves a bot-verification page, so resolve files via the figshare API
instead. DepMap public releases are figshare depositions. Use
`https://api.figshare.com/v2/` - search for the DepMap Public release (e.g.
query "DepMap Public 24Q4" or browse the institution's articles) to get the
article's file list and download URLs on ndownloader.figshare.com.
Investigate and report which release you actually resolved - do not assume a
version number. If figshare resolution fails entirely, say so explicitly in your
structured output rather than inventing data, and fall back to trying
`https://depmap.org/portal/api/download/files` (may be blocked by the interstitial).

Files you want (names vary slightly by release):
- `CRISPRGeneEffect.csv` (Chronos gene effect, cell lines x genes)
- `OmicsExpressionProteinCodingGenesTPMLogp1.csv` (log2(TPM+1))
- `Model.csv` (cell line annotation incl. OncotreeLineage / primary disease)

`kernel.py` should provide roughly:
- `resolve_depmap_release(release=None)` -> dict of {filename: download_url},
  via the figshare API. Document how it resolves and cache the manifest.
- `fetch_depmap_file(url, dest, expected_min_bytes=...)` -> path, cache-aware,
  with a precondition that the download is not an HTML error page.
- `load_gene_effect(path, genes=None)` -> DataFrame. **Memory discipline**:
  CRISPRGeneEffect.csv is large and wide (genes as columns, ~1000 lines). Read
  with `usecols` on the requested genes rather than loading everything.
- `load_expression(path, genes=None)` -> DataFrame, same discipline.
- `parse_gene_columns(header)` -> map "SYMBOL (ENTREZ)" column names to symbols.
- `lineage_dependency(gene_effect, model_df, genes, lineage_col=...)` ->
  per-lineage mean gene effect per gene, plus n lines per lineage.
- `correlate_with_external(per_lineage_effect, external_score, lineage_map)` ->
  DataFrame gene, spearman rho, p, n_lineages, fdr (Benjamini-Hochberg).
  Requires an explicit lineage->tumor-type mapping argument; do NOT guess
  silently - the mapping must be passed in and recorded.
- `select_panel(expr, gene_effect, model_df, receptor, cofactors=(), lineage=None,
  high_pct=85, low_pct=15, activity_genes=(), exclude_subtypes=())` ->
  DataFrame of candidate high/low lines with expression percentile, cofactor
  level, baseline pathway activity, and gene effect.

PRECONDITIONS that raise: requested genes present in the file header; downloads
are not HTML; every selected line has a lineage annotation.

## Step 2 - Test lineage vulnerabilities (plan step: "Test lineage vulnerabilities")

Take the floor-PASSING receptor subunits from `crosstalk_network.csv`
(version 8086ab6e-ce8d-466a-804e-31b39df5e0e1, use rows where passed_floor is
True). Build an explicit, documented mapping from DepMap OncotreeLineage (or
primary disease) to TCGA tumor-type abbreviations used in
`pancancer_susceptibility_ranking.csv` (c8e165e9-4061-43cb-b270-af809942d9d6).
Record the mapping as a table - it is a judgment call the paper must disclose,
and some TCGA types will have no DepMap counterpart.

Then, for each receptor gene, correlate per-lineage mean CRISPR gene effect
against that lineage's compatibility score. Report FDR-corrected results.
ALSO report, separately and plainly: which receptors are dependencies at all
(mean gene effect < -0.5 in any lineage), and what the gene effect is for the
specific receptors in the top-ranked axes. A receptor that is expressed but not
essential means the axis is a candidate signalling conduit, NOT a drug target -
state that distinction explicitly.

Sanity check your pipeline: include known common-essential genes (e.g. RPL3,
POLR2A, PLK1) as positive controls - they MUST show strongly negative gene
effect. If they do not, your parsing is wrong. Report the control values.

## Step 3 - Cell-line panel (plan step: "Build a testable cell-line panel")

The leading axis after the floor filter should be determined from the data, not
assumed - inspect `crosstalk_network.csv` for the highest-availability ligands
whose receptors pass the floor, and pick the axis to operationalize. (WNT9A ->
FZD/LRP receptors is a likely candidate given WNT9A availability 2.32 and that
the WNT pathway is the largest single contributor to the ranking, but verify
against the data.)

Select a responder arm (high receptor, cofactor subunit present, LOW baseline
pathway activity so induction is measurable) and a lineage-matched control arm
(low receptor). Exclude subtypes that would confound receptor level with lineage.
Report: fold difference in receptor expression between arms, whether receptor
level tracks library size / housekeeping (it must not), and gene effect for the
receptor in those lines (to confirm a knockdown arm is interpretable).

## Step 4 - Save artifacts

Exact filenames:
- `depmap_lineage_dependency.csv` - per receptor x lineage mean gene effect, n lines.
- `depmap_compatibility_correlation.csv` - gene, rho, p, fdr, n_lineages, verdict.
- `depmap_lineage_map.csv` - the DepMap lineage -> TCGA type mapping you used,
  including unmapped types.
- `cell_line_panel.csv` - arm, cell line name, DepMap ID, subtype, receptor
  expression, panel percentile, cofactor expression, baseline activity, gene effect.
- `depmap_summary.json` - release resolved, file URLs, n lines, n genes,
  positive-control gene effects, n receptors tested, n significant after FDR,
  and an explicit statement of the expression-vs-dependency distinction.
- the two skill files.

Call update_step_status for each of the three plan steps using the exact titles
in parentheses.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "depmap_release_resolved": {
      "type": "string"
    },
    "n_cell_lines": {
      "type": "integer"
    },
    "positive_control_gene_effects": {
      "type": "object",
      "description": "common-essential genes; must be strongly negative or parsing is wrong"
    },
    "n_receptors_tested": {
      "type": "integer"
    },
    "n_significant_after_fdr": {
      "type": "integer"
    },
    "any_receptor_essential": {
      "type": "string",
      "description": "which receptors if any have mean gene effect < -0.5, and in which lineage"
    },
    "axis_operationalized": {
      "type": "string",
      "description": "which ligand->receptor axis you built the panel for and why, from the data"
    },
    "panel_fold_difference": {
      "type": "string"
    },
    "expression_vs_dependency_statement": {
      "type": "string"
    },
    "lineage_mapping_caveats": {
      "type": "string",
      "description": "which TCGA types had no DepMap counterpart"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object"
      }
    },
    "skill_paths": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "depmap_release_resolved",
    "positive_control_gene_effects",
    "n_receptors_tested",
    "n_significant_after_fdr",
    "expression_vs_dependency_statement",
    "artifacts",
    "skill_paths"
  ]
}
```


---

## Wave 1 — Normal tissue

`model=claude-opus-5`, `reasoning_effort=high`, status=cancelled, duration=198.9 min

### Task statement

```

## Shared engineering rules

**Environment**: use `environment="vibemed"` on every python/bash call (has pandas,
numpy, scipy, statsmodels, seaborn, matplotlib, scikit-learn, requests, pyreadr, rdata).

**Machine**: 8 GiB RAM, 8 cores, no GPU, no remote compute. Stream large files;
never load multi-GB into memory.

**Network granted**: ftp.ncbi.nlm.nih.gov, www.ncbi.nlm.nih.gov, www.proteinatlas.org,
raw.githubusercontent.com, api.gdc.cancer.gov, api.figshare.com,
ndownloader.figshare.com, bioconductor.org, depmap.org (NOTE: depmap.org itself
serves a bot-verification interstitial - use figshare for actual DepMap files).
compbio.cn / timer.comp-genomics.org are BLOCKED and NOT grantable.
If blocked on something needed, call request_network_access. NEVER fabricate data.

**Skill authoring**: write `skills/<name>/SKILL.md` + `skills/<name>/kernel.py`.
Do NOT publish (orchestrator publishes centrally). `kernel.py` MUST NOT use
`host.*` - pure stdlib + pandas/numpy/scipy so it runs outside this platform.
Only functions/imports/literal constants at top level (no top-level classes or
decorators). Guardrails = assertions that RAISE, not prose. SKILL.md needs YAML
frontmatter (name, description), signatures, worked example, pre/postconditions,
and a "Standalone use (no agent)" section.

**Scientific integrity**: every number must come from code you ran on real data.
Never invent or carry over a number. A previous PAPER.md in this project was full
of unbacked numbers and has been DISCARDED - do not reuse its values or gene/cell
line names. Report null and contradictory results plainly.

**Artifacts**: save with `save_artifacts(files=[...], language="python")` using
EXACT filenames requested; report filename + version_id in structured output.

## What the project has established so far (all from real data, verified)

Hypothesis: heart failure may promote cancer via secreted factors; susceptibility
depends on ligand availability (source side) AND receiver capacity (tumor side).

Phase 0-1 results you should build on:
- HF signature from GSE116250 (discovery, n=50 failing vs 14 non-failing) and
  GSE141910 (validation, MAGNet). 3084 HF-UPregulated genes replicated in both
  cohorts (FDR<0.05, same sign). Directional concordance 0.881, log2FC Spearman
  0.541. Label-permutation control gave 0 significant genes.
- 295 of those are HPA-annotated secreted; 72 are CellChatDB ligands.
- 188 Secreted-Signaling interactions from 51 HF-up secreted ligands, 38 pathways.
- TCGA: 9538 primary solid tumor samples, 30 tumor types (LAML/DLBC/THYM excluded).
  Expression unit is RSEM normalized_count (EBPlusPlus batch-corrected),
  **NOT TPM** - GAPDH median is ~63,656. Any floor must be relative to
  housekeeping, never an absolute "2 TPM".
- Expression floor = 0.005 x GAPDH median = 345.67 RSEM units. It drops 28 of 188
  interactions, INCLUDING all three PENK->opioid receptor interactions
  (PENK_OPRK1, PENK_OPRD1, PENK_OPRM1). OPRD1 max across tumor types is 4e-5 of
  GAPDH; OPRK1 1.4e-3. The opioid axis is dead on receptor expression grounds.
- Top ligands by availability: NPPA (3.05), NPPB (3.03), PENK (2.85), WNT9A (2.32),
  CCL5 (1.59), CXCL9 (1.34).
- Primary ranking (floor-filtered) top 6: UCS 0.679 (95% CI 0.400-0.914),
  OV 0.432, KIRC 0.377, THCA 0.285, PRAD 0.272, PAAD 0.244.
  Double-centered ranking reorders to: ACC, OV, UCEC, PRAD, PCPG, UCS.
  Primary vs double-centered Spearman = 0.318 (correction-sensitive!).
- KEY METHODOLOGICAL FINDING: the matched-random specificity null PASSES when the
  floor is applied (rho_mean 0.107, 1.7% of random sets exceed rho 0.5) but FAILS
  without it (rho_mean 0.363, 21.5%). An unmatched/naive random pool gives 0.323.
  So the "global expression axis confound" is substantially CAUSED by scoring
  untranscribed receptors, and the null pool must be expression-matched.
- Receiver capacity vs tumor microenvironment: only 1 of 29 Thorsson TME features
  survives FDR (B-cell receptor score rho=0.57, FDR=0.028). The "exploit /
  microenvironment" half of the hypothesis is largely UNSUPPORTED.

Useful artifact version_ids (use host.artifact_path(vid) or {{artifact:vid}}):
- hf_signature_replicated.csv = 76685c52-1489-4bd0-9e51-60c282afbf19
  (IMPORTANT: join on the `gene_symbol` column, NOT `gene` which is Ensembl)
- hpa_secretome_annotation.csv = efb93594-4a96-4a02-8bab-e3030adece36
- cellchatdb_receptor_subunits.csv = 4aabe49f-60f3-4179-a333-ba2cc7c548a0
- tcga_expr_by_tumor_type.csv = 07e2a94c-6093-4cfb-b343-56f042647e20
- crosstalk_network.csv = 8086ab6e-ce8d-466a-804e-31b39df5e0e1
- pancancer_susceptibility_ranking.csv = c8e165e9-4061-43cb-b270-af809942d9d6
- twofactor_decomposition.csv = 23bff2cf-a7f4-42d8-b520-ffcf9b279745
- literature_citations.csv = 940dda9f-ce28-4f2a-a1b6-bcc9d5cc26a8
- ligand_evidence.csv = 8f8d6101-f192-48d8-80ba-10e96f46e295

# Your task: normal-tissue and cardiac-specificity comparison

GTEx is available through the expression MCP connector. From the `repl` tool run
`search_skills(query="GTEx tissue expression")` then load the `mcp-expression`
skill doc for exact method names/signatures (methods include
`gtex_median_expression`, `gtex_gene_expression`). MCP calls ONLY work in the
`repl` tool - write results to `handoff/*.json` and analyze in a `python` cell
(environment="vibemed").

Note the unit mismatch you must handle honestly: TCGA here is RSEM
normalized_count (EBPlusPlus), GTEx is TPM (or the connector's stated unit -
check it). You therefore CANNOT subtract them directly. Use rank-based or
within-dataset-standardized comparisons, and state the limitation explicitly.
Do not manufacture a fake common scale.

## Step 1 - Receptors vs matched normal tissue (plan step: "Compare receptors to normal tissue")

Take the floor-PASSING receptor subunits from `crosstalk_network.csv`
(8086ab6e-ce8d-466a-804e-31b39df5e0e1, rows with passed_floor True), focusing on
the receptors in the top-contributing interactions.
Build an explicit, documented TCGA-tumor-type -> GTEx-tissue mapping (e.g.
OV->Ovary, UCS/UCEC->Uterus, KIRC->Kidney Cortex, PRAD->Prostate, THCA->Thyroid,
PAAD->Pancreas, BRCA->Breast, LUAD/LUSC->Lung, ...). Save the mapping; note
tumor types with no clean normal counterpart (e.g. UCS, MESO, TGCT may be
imperfect) and treat those results as weaker.
For each (receptor, tumor type) pair, compare tumor expression rank against the
matched normal tissue rank and report a direction (enriched / depleted /
comparable) with the caveat that this is rank-based, not a fold change.

Also do a within-TCGA check that does not need GTEx at all: for tumor types
where TCGA has adjacent-normal samples ('11' sample type code), you could
compare - but the extracted matrix was filtered to primary tumors only, so note
this as a limitation rather than re-downloading the 1.9 GB matrix. Use judgment
about whether a targeted re-extraction is worth it; if you do it, stream, do not
load the whole matrix.

## Step 2 - Ligand cardiac specificity (plan step: "Verify ligand cardiac specificity")

For the top HF-up secreted ligands (NPPA, NPPB, PENK, WNT9A, CCL5, CXCL9, POSTN,
GDF15, THBS4, COMP, IGF1, INHBA, MDK, PTN and the rest of the 51 in
crosstalk_network.csv), pull GTEx median expression across all tissues and
compute a tissue-specificity metric (e.g. heart LV expression relative to the
median across tissues, and the rank of heart among tissues).
Classify each ligand as cardiac-enriched / cardiac-expressed-but-not-enriched /
not-cardiac. This directly tests whether the failing heart is a plausibly
DISTINCTIVE source. Expect NPPA/NPPB to be strongly cardiac and expect several
others NOT to be - report that plainly, since it bounds the hypothesis.

Cross-reference `ligand_evidence.csv` (8f8d6101-f192-48d8-80ba-10e96f46e295)
which records whether each ligand is established as elevated in human HF at the
protein level.

## Step 3 - Emit tables (plan step: "Emit normal-comparison table")

Exact filenames:
- `receptor_tumor_vs_normal.csv` - receptor, tumor_type, gtex_tissue,
  tumor_rank_metric, normal_rank_metric, direction, caveat.
- `tumor_normal_tissue_map.csv` - the mapping used, with a quality flag per row.
- `ligand_tissue_specificity.csv` - gene, heart_lv_median, median_across_tissues,
  ratio, heart_rank_among_tissues, n_tissues, classification,
  established_in_human_HF (joined from ligand_evidence.csv).
- `specificity_summary.json` - the unit-mismatch limitation stated explicitly,
  n receptors enriched vs comparable vs depleted, n ligands cardiac-enriched,
  and which top-ranked axes survive both checks.

Call update_step_status for each of the three plan steps using the exact titles
in parentheses.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "gtex_unit": {
      "type": "string"
    },
    "unit_mismatch_handling": {
      "type": "string",
      "description": "how you handled TCGA RSEM vs GTEx TPM without faking a common scale"
    },
    "n_receptors_compared": {
      "type": "integer"
    },
    "receptors_enriched_in_tumor": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "receptors_not_enriched": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "ligands_cardiac_enriched": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "ligands_not_cardiac_enriched": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "cardiac_specificity_verdict": {
      "type": "string",
      "description": "does the failing heart look like a distinctive source of these ligands?"
    },
    "axes_surviving_both_checks": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "mapping_caveats": {
      "type": "string"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object"
      }
    }
  },
  "required": [
    "gtex_unit",
    "unit_mismatch_handling",
    "ligands_cardiac_enriched",
    "ligands_not_cardiac_enriched",
    "cardiac_specificity_verdict",
    "artifacts"
  ]
}
```


---

## Wave 1 — Normal tissue

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=16.6 min

### Task statement

```

## Shared engineering rules

**Environment**: use `environment="vibemed"` on every python/bash call (pandas,
numpy, scipy, statsmodels, seaborn, matplotlib, scikit-learn, requests, pyreadr, rdata).

**Machine**: 8 GiB RAM, 8 cores. DISK IS TIGHT (~10 GiB free). Stream large files.

**Network granted**: ftp.ncbi.nlm.nih.gov, www.ncbi.nlm.nih.gov, www.proteinatlas.org,
raw.githubusercontent.com, api.gdc.cancer.gov, api.figshare.com,
ndownloader.figshare.com, bioconductor.org, depmap.org.
compbio.cn / timer.comp-genomics.org are BLOCKED and NOT grantable.
If blocked on something needed, call request_network_access. NEVER fabricate data.

**Scientific integrity**: every number must come from code you ran on real data.
Never invent or carry over a number. A previous PAPER.md in this project was full
of unbacked numbers and has been DISCARDED - do not reuse its values. Report null
and deflationary results plainly; they are the point of your task.

**Artifacts**: save with `save_artifacts(files=[...], language="python")` using
EXACT filenames requested. Report filename + version_id in structured output.

## Additional established results (phase 2, DepMap arm - completed)

- DepMap 24Q4 Public resolved via figshare article 27993248. 1178 cell lines.
  Positive controls confirm correct parsing: RPL3 -2.533, PLK1 -2.693,
  POLR2B -2.418 (POLR2A absent from 24Q4 matrix).
- 81 floor-passing receptors x 22 mapped lineages: **0 significant at FDR<0.05**
  for compatibility-correlated dependency (3 nominal only).
- NPR1 (receptor for NPPA/NPPB, the two highest-availability ligands) pan-lineage
  mean gene effect -0.163. All WNT9A receptor subunits between LRP5 -0.163 and
  FZD6 +0.146. NO WNT/FZD/LRP or NPR1 lineage cell crosses -0.5.
- Conclusion recorded: the nominated axes are candidate signalling conduits, NOT
  drug targets.
- Cell-line panel built for WNT9A->FZD7 in DepMap Uterus (maps to UCS + UCEC),
  4v4 arms, 48.7x linear TPM fold difference, housekeeping-independent
  (rho=0.131, p=0.447).
- Note: ACC (rank 1 under double-centering) has only 1 CRISPR-screened adrenal
  line, so it is untestable in DepMap.

## IMPORTANT operational note

A previous attempt at your task failed with 'database or disk is full' because
another track saved two ~500 MB raw source files as artifacts. Those have been
deleted and space reclaimed, but disk is tight (about 10 GiB free).
**Do NOT save large raw downloads as artifacts.** Save only the small derived
tables listed in your task. Keep any bulk download in the workspace and let it
be discarded. If you must cache something big, keep it under ~50 MB.

# Your task: normal-tissue and cardiac-specificity comparison

GTEx is available through an MCP connector. From the `repl` tool run
`search_skills(query="GTEx tissue expression")` then load the `mcp-expression`
skill doc for exact method names/signatures (methods include
`gtex_median_expression`, `gtex_gene_expression`). MCP calls ONLY work in the
`repl` tool - write results to `handoff/*.json` and analyze in a `python` cell
(environment="vibemed"). If the GTEx connector is unavailable, say so explicitly
in your structured output rather than substituting a different source silently.

Unit mismatch you MUST handle honestly: TCGA here is RSEM normalized_count
(EBPlusPlus), GTEx is TPM (check the connector's stated unit). You CANNOT
subtract or ratio them directly. Use rank-based or within-dataset-standardized
comparisons and state the limitation. Do not manufacture a fake common scale.

## Step 1 - Receptors vs matched normal tissue (plan step: "Compare receptors to normal tissue")

Input: `crosstalk_network.csv` = artifact version 8086ab6e-ce8d-466a-804e-31b39df5e0e1
(use rows where passed_floor is True). Focus on receptors in the top-contributing
interactions - at minimum NPR1 and the WNT9A receptor subunits (FZD1-10, LRP5, LRP6).
Per-tumor-type receptor expression: `tcga_expr_by_tumor_type.csv` =
07e2a94c-6093-4cfb-b343-56f042647e20 (genes x tumor type medians).

Build an explicit, documented TCGA-type -> GTEx-tissue mapping (OV->Ovary,
UCS/UCEC->Uterus, KIRC->Kidney Cortex, PRAD->Prostate, THCA->Thyroid,
PAAD->Pancreas, BRCA->Breast, LUAD/LUSC->Lung, etc). Flag tumor types with no
clean normal counterpart (UCS, MESO, TGCT are imperfect) and treat those as weaker.
For each (receptor, tumor type): compare tumor expression rank vs matched normal
tissue rank, report direction (enriched / comparable / depleted) with the caveat
that this is rank-based, not a fold change.

## Step 2 - Ligand cardiac specificity (plan step: "Verify ligand cardiac specificity")

For the HF-up secreted ligands in crosstalk_network.csv (51 distinct; at minimum
NPPA, NPPB, PENK, WNT9A, CCL5, CXCL9, POSTN, GDF15, THBS4, COMP, IGF1, INHBA,
MDK, PTN), pull GTEx median expression across ALL tissues and compute:
heart-LV expression relative to the median across tissues, and heart's RANK
among tissues. Classify each as cardiac-enriched / cardiac-expressed-but-not-
enriched / not-cardiac. This directly tests whether the failing heart is a
plausibly DISTINCTIVE source. Expect NPPA/NPPB strongly cardiac and several
others NOT - report plainly, it bounds the hypothesis.

Cross-reference `ligand_evidence.csv` = 8f8d6101-f192-48d8-80ba-10e96f46e295
which records whether each ligand is established as elevated in human HF at
protein level.

## Step 3 - Emit tables (plan step: "Emit normal-comparison table")

Exact filenames (all SMALL derived tables - do not save bulk downloads):
- `receptor_tumor_vs_normal.csv` - receptor, tumor_type, gtex_tissue,
  tumor_rank_metric, normal_rank_metric, direction, caveat.
- `tumor_normal_tissue_map.csv` - the mapping used, with a quality flag per row.
- `ligand_tissue_specificity.csv` - gene, heart_lv_median, median_across_tissues,
  ratio, heart_rank_among_tissues, n_tissues, classification,
  established_in_human_HF (joined from ligand_evidence.csv).
- `specificity_summary.json` - the unit-mismatch limitation stated explicitly,
  n receptors enriched vs comparable vs depleted, n ligands cardiac-enriched,
  and which top-ranked axes survive both checks.

Call update_step_status for each of the three plan steps using the exact titles
in parentheses.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "gtex_unit": {
      "type": "string"
    },
    "gtex_source_used": {
      "type": "string",
      "description": "which connector/method actually supplied GTEx data, or an explicit statement it was unavailable"
    },
    "unit_mismatch_handling": {
      "type": "string"
    },
    "n_receptors_compared": {
      "type": "integer"
    },
    "receptors_enriched_in_tumor": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "receptors_not_enriched": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "ligands_cardiac_enriched": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "ligands_not_cardiac_enriched": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "cardiac_specificity_verdict": {
      "type": "string"
    },
    "axes_surviving_both_checks": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "mapping_caveats": {
      "type": "string"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object"
      }
    }
  },
  "required": [
    "gtex_unit",
    "gtex_source_used",
    "unit_mismatch_handling",
    "ligands_cardiac_enriched",
    "ligands_not_cardiac_enriched",
    "cardiac_specificity_verdict",
    "artifacts"
  ]
}
```


---

## Wave 2 — Failing-heart source

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=9.9 min

### Task statement

```

ENVIRONMENT: use environment="vibemed" for every python cell (pandas/numpy/scipy/statsmodels
installed). If a package is missing, install it with manage_packages into vibemed.

THE PRE-SPECIFICATION IS BINDING. Read it first: /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/cd788e4a-50a9-4e62-a4db-e7a74c54100c/v733b9494_PRESPECIFICATION.md
Every threshold, direction, seed and pass/fail rule is fixed there. You may NOT adjust a
threshold, swap an outcome variable, change a covariate set, or substitute a correction method
after seeing a result. A near-miss (e.g. p=0.06 against p<0.05) is a FAIL. Report FAILs with
exactly the same prominence as PASSes. If a criterion cannot be computed at all, say so
explicitly rather than substituting a proxy.

GLOBAL PROHIBITION P1: never place TCGA RSEM, GTEx TPM, GSE116250 RPKM and GSE141910
normalised counts on a common numeric scale. Cross-dataset comparisons are ordinal only.

Print only computed values, not narration. Save every table you produce with save_artifacts
and report the version_ids in your structured output.

TRACK A — Re-test cardiac ligand distinctiveness using a FAILURE-INDUCED definition.

The parent screen asked whether each ligand is enriched in GTEx *healthy* left ventricle.
That is the wrong question for a hypothesis about the FAILING heart: a ligand unremarkable in
healthy myocardium but strongly induced in failure, and unchanged elsewhere, is a distinctive
failing-cardiac source, and the original filter discarded it. Correct that, then re-run axis
survival. The parent result was ZERO surviving axes out of 160 floor-passing.

INPUT ARTIFACTS (use the literal markers):
- discovery DE:        /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/a06e693e-0e8a-47c1-a618-3f29ee78d4db/v6a2d301f_hf_de_GSE116250.csv   (hf_de_GSE116250.csv; RPKM-derived; cols gene_symbol, log2FC, fdr, mean_expr)
- validation DE:       /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/f76fe119-05ac-415b-a1f1-d3b80bf512c8/vcac6c07c_hf_de_GSE141910.csv   (hf_de_GSE141910.csv; counts-derived; same cols)
- GTEx ligand ranks:   /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/d58b9768-3152-4c89-8682-586508e35b55/va82f4063_ligand_tissue_specificity.csv   (ligand_tissue_specificity.csv; 42 rows; cols include best_heart_rank, heart_rank_among_tissues, n_tissues, classification, top_tissue)
- crosstalk network:   /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/890af4bd-6338-410b-8744-e13de45260eb/v8086ab6e_crosstalk_network.csv   (crosstalk_network.csv; 296 rows; cols ligand_subunit_gene, receptor, receptor_subunit_gene, pathway_name, passed_floor, ligand_availability, interaction_name)
- receptor vs normal:  /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/eaf40910-43f9-4d69-87a5-fc359ef800ff/v9c318e4e_receptor_tumor_vs_normal.csv   (receptor_tumor_vs_normal.csv; 2430 rows; cols receptor, tumor_type, delta_percentile_bg_corrected, direction_gtex_rank, passes_tcga_floor_in_type, n_tcga_adjacent_normals, adjnormal_median_rsem, tumor_median_rsem, mapping_quality)
- parent axis survival:/Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/bbbb5618-f77e-4ce2-a4fe-4046c6b48263/v2d0e9667_axis_survival.csv   (axis_survival.csv; 146 rows; the result you are superseding)

Note both mean_expr columns are already on a log2 scale (GSE116250 median 2.85, max 18.58;
GSE141910 median 10.49, max 21.03), so criterion A2's linear ratio needs 2**mean_expr.
GAPDH mean_expr is 10.942161 in GSE116250 and 18.569934 in GSE141910 — verify these from the
files rather than trusting these numbers.

STEPS:
1. Apply criteria A1, A2, A3 exactly as pre-specified to all 51 ligand subunits in the network
   (not just the 42 in the GTEx table — report how many lack a GTEx row and how you treat them;
   a ligand with no GTEx row cannot satisfy A3 and must be recorded as failing A3, not skipped).
   Produce ligand_specificity_failure_based.csv with one row per ligand, a column per criterion
   (boolean plus the underlying value), and a final pass column. Report per-criterion attrition.
2. Re-run axis survival: an axis survives if its ligand passes A1+A2+A3 AND its receptor side
   passes UNCHANGED parent criteria (tumour-enriched vs matched adjacent normal AND clears the
   0.005xGAPDH TCGA floor in that type). Do not modify the receptor side in any way.
   Produce axis_survival_corrected.csv with the same shape as the parent axis_survival.csv plus
   the ligand-criterion columns, and state n_surviving.
3. Score against the pre-specified pass criterion: PASS if n_surviving >= 1, FAIL if 0,
   CRITERION UNINFORMATIVE if > 20. Apply that third verdict honestly if it triggers — a
   relaxation admitting >20 of 160 axes is not a discovery.
4. Report, whatever the outcome: the count under the ORIGINAL criterion for contrast; every
   ligand newly admitted by the relaxed A3 and what its top GTEx tissue is; and which of A1/A2/A3
   is the binding constraint for each of the 51.

Do NOT build final publication figures — the synthesis phase does that. A quick diagnostic plot
for your own checking is fine but do not spend time styling it.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "n_surviving_corrected": {
      "type": "integer"
    },
    "n_surviving_original": {
      "type": "integer"
    },
    "verdict": {
      "type": "string",
      "enum": [
        "PASS",
        "FAIL",
        "CRITERION_UNINFORMATIVE"
      ]
    },
    "n_ligands_pass_all": {
      "type": "integer"
    },
    "attrition": {
      "type": "object",
      "description": "how many of the 51 ligands each criterion removes"
    },
    "newly_admitted_ligands": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "surviving_axes": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "binding_constraint_counts": {
      "type": "object"
    },
    "artifacts": {
      "type": "object",
      "description": "filename -> version_id"
    },
    "caveats": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "n_surviving_corrected",
    "n_surviving_original",
    "verdict",
    "n_ligands_pass_all",
    "attrition",
    "artifacts"
  ]
}
```


---

## Wave 2 — NPR1 exploratory

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=30.4 min

### Task statement

```

ENVIRONMENT: use environment="vibemed" for every python cell (pandas/numpy/scipy/statsmodels
installed). If a package is missing, install it with manage_packages into vibemed.

THE PRE-SPECIFICATION IS BINDING. Read it first: /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/cd788e4a-50a9-4e62-a4db-e7a74c54100c/v733b9494_PRESPECIFICATION.md
Every threshold, direction, seed and pass/fail rule is fixed there. You may NOT adjust a
threshold, swap an outcome variable, change a covariate set, or substitute a correction method
after seeing a result. A near-miss (e.g. p=0.06 against p<0.05) is a FAIL. Report FAILs with
exactly the same prominence as PASSes. If a criterion cannot be computed at all, say so
explicitly rather than substituting a proxy.

GLOBAL PROHIBITION P1: never place TCGA RSEM, GTEx TPM, GSE116250 RPKM and GSE141910
normalised counts on a common numeric scale. Cross-dataset comparisons are ordinal only.

Print only computed values, not narration. Save every table you produce with save_artifacts
and report the version_ids in your structured output.

TRACK B — EXPLORATORY: is NPR1 a lost antiproliferative receptor?

STATUS LABEL, which must appear on every table and in your summary: this hypothesis was
generated BY the refutation, not predicted by the original screen. The original hypothesis was
that a cardiac ligand engages a tumour receptor. This is close to the opposite — that tumours
benefit from SHEDDING that receptor. Never write it up as a confirmed prediction of the
original design.

THE CONJUNCTION RULE IS THE POINT. B1 AND B2 AND B3 must ALL pass AND agree in direction for
the track to pass. One or two passing is a FAIL for the track, reported as a partial result
naming the failing conjunct. Do not report "2 of 3 passed" as a success.

INPUT ARTIFACTS:
- receptor vs normal:  /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/eaf40910-43f9-4d69-87a5-fc359ef800ff/v9c318e4e_receptor_tumor_vs_normal.csv   (receptor_tumor_vs_normal.csv — identifies the 15 tumour types with NPR1 depletion and unit-matched adjacent normals)
- TCGA sample map:     /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/3cc44965-0c60-49d5-aea1-2610a3734b98/v056092b3_tcga_sample_map.csv   (tcga_sample_map.csv; 9538 rows; cols aliquot_barcode, patient, tumor_type, sample_type_code)
- TME / proliferation: /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/7509ec6f-2fd5-4088-a193-1cbd1e2a292d/vfbfd0a30_tcga_tme_scores.csv   (tcga_tme_scores.csv; 9532 rows; col sample_barcode, Module11_Prolif_score)
- per-type expression: /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/51bd31c5-0e36-4709-ac73-f5ffb358e9ce/v07e2a94c_tcga_expr_by_tumor_type.csv   (tcga_expr_by_tumor_type.csv; NPR1 present; per-type medians only — NOT per-sample)

For PER-SAMPLE NPR1 you need tcga_expression_subset.csv.gz from the parent project. Resolve it
with host.artifacts(filename="tcga_expression_subset.csv.gz", exact=True) in a python cell and
read it via host.artifact_path(vid). If NPR1 is absent from that subset, re-stream just NPR1
and GAPDH from the PanCanAtlas matrix using the published `tcga-pancan` skill
(skill(skill="tcga-pancan")) — it streams row-wise and is memory-bounded for an 8 GiB machine.
Do not load the full ~1.9 GB matrix into memory.

B1 — SURVIVAL. TCGA-CDR clinical: GDC file id 1b5f413e-a8d1-4d10-92eb-7c4ae739ed81
(https://api.gdc.cancer.gov/data/<id>) — confirmed reachable, arrives as
TCGA-CDR-SupplementalTableS1.xlsx (needs openpyxl). Endpoint OS / OS.time ONLY; progression-free
interval is NOT a permitted fallback (prohibition P3). Cox PH, NPR1 z-scored WITHIN tumour type,
stratified by tumour type, covariates age at diagnosis and sex. Pre-specified direction HR < 1.
PASS = pooled p < 0.05 two-sided AND HR < 1. Also report per-type Cox with BH-FDR and the count
of types with HR < 1. Use lifelines if available, else statsmodels PHReg — say which.

B2 — PROLIFERATION. Spearman NPR1 vs Module11_Prolif_score WITHIN each tumour type.
Pre-specified direction: negative rho. PASS = median within-type rho < 0 AND >= 8 of 15 types
individually significant at BH-FDR < 0.05 with negative rho. Join tcga_tme_scores.sample_barcode
to the expression barcodes carefully — barcode granularity differs between files (aliquot vs
sample level); state exactly how you matched them and how many samples matched.

B3 — GENOMIC LOSS vs SILENCING. TCGA ABSOLUTE allele-specific copy number:
GDC 0f4f5701-7b61-41ae-bda9-2805d1ca9781 -> TCGA_mastercalls.abs_segtabs.fixed.txt
(cols Sample, Chromosome, Start, End, Num_Probes, Length, Modal_HSCN_1, Modal_HSCN_2,
Modal_Total_CN, Subclonal_HSCN_a1 ...). Confirmed reachable, hg19/GRCh37.
NPR1 GRCh37 locus, already resolved: chr1:153,651,113-153,666,468 (ENSG00000169418).
Per sample take Modal_Total_CN of overlapping segment(s); call NPR1 "lost" if its copy number is
below that sample's modal ploidy (derive modal ploidy per sample from its own segment profile,
or use the ABSOLUTE purity/ploidy table GDC 4f277128-f793-4354-a13d-30cc7fe9f6b5 which has
per-sample purity and ploidy — say which you used).
NULL: 1000 random loci matched on chromosome arm 1q and on segment-length distribution,
seed 20260728. PASS = NPR1 loss frequency above the 95th percentile of that null.
Methylation is a CONDITIONAL SECONDARY read and is NOT part of the conjunction. The merged
PanCanAtlas HumanMethylation27/450 beta matrix is GDC d82e2c44-89eb-43d9-b6d3-712732bf6a53
(reachable, but it is very large and probe-to-gene mapping is not in hand). Only attempt it if
B1+B2+B3 all pass and you can get a probe mapping from an allowlisted source; otherwise record
"methylation not assessed" as a limitation. Absence of methylation data is NOT grounds to relax B3.

DELIVERABLES: npr1_survival.csv, npr1_proliferation.csv, npr1_copynumber.csv, each carrying an
explicit exploratory-status column or header comment. State the conjunction verdict plainly.

Do NOT build final publication figures — the synthesis phase does that.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "conjunction_verdict": {
      "type": "string",
      "enum": [
        "PASS",
        "FAIL"
      ]
    },
    "b1_survival": {
      "type": "object",
      "description": "pass bool, HR, p, n, model used, n_types_HR_below_1"
    },
    "b2_proliferation": {
      "type": "object",
      "description": "pass bool, median rho, n_types_sig_negative, n_types_tested, barcode match method and n matched"
    },
    "b3_copynumber": {
      "type": "object",
      "description": "pass bool, NPR1 loss freq, null 95th pct, ploidy source"
    },
    "failing_conjuncts": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "methylation_assessed": {
      "type": "boolean"
    },
    "artifacts": {
      "type": "object",
      "description": "filename -> version_id"
    },
    "caveats": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "conjunction_verdict",
    "b1_survival",
    "b2_proliferation",
    "b3_copynumber",
    "failing_conjuncts",
    "artifacts"
  ]
}
```


---

## Wave 2 — Receptor discrimination

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=7.5 min

### Task statement

```

ENVIRONMENT: use environment="vibemed" for every python cell. Install into vibemed via
manage_packages if something is missing.

THE PRE-SPECIFICATION IS BINDING. Read it FIRST: /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/c615b111-e38e-4459-a6b5-9e8191e0d33c/v1babf914_PRESPECIFICATION_PTN.md
Thresholds, directions, seeds and kill criteria are fixed there. You may NOT adjust a threshold,
swap a cohort/endpoint/statistic because the pre-specified one gave an unwelcome answer, or promote
a near-miss. Report FAILs with exactly the same prominence as PASSes.

CONTEXT. The parent screen found 0 surviving ligand-receptor axes. A corrected specificity
criterion (scoring cardiac distinctiveness by FAILURE INDUCTION rather than healthy-heart tissue
rank) produced 7 surviving axes, ALL with PTN (pleiotrophin) as the ligand: PTN->SDC1 (enriched in
6 tumour types), PTN->NCL (5), PTN->ITGAV (5), PTN->SDC3 (4), PTN->SDC4 (4), PTN->PTPRZ1 (2),
PTN->SDC2 (1) — 27 (receptor, tumour-type) pairs in total. PTN was admitted partly BY a relaxation
I introduced (heart rank 23 of 52 GTEx tissues, upper-half rather than first). YOUR JOB IS TO TRY
TO KILL THIS RESULT, not to support it. The pre-registered expectation is that PTN fails.

GLOBAL PROHIBITION P1: never place TCGA RSEM, GTEx TPM, RPKM and normalised counts on a common
numeric scale. Cross-dataset comparisons are ordinal only.

Print only computed values, not narration. Save tables with save_artifacts and report version_ids.
Do NOT build publication figures — synthesis does that.

TRACK 1 — Does the receptor side discriminate at all, or would almost any receptor pass?

INPUTS:
- receptor vs normal:  /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/eaf40910-43f9-4d69-87a5-fc359ef800ff/v9c318e4e_receptor_tumor_vs_normal.csv  (receptor_tumor_vs_normal.csv; 2430 rows; cols receptor, tumor_type, direction_within_tcga, direction_gtex_rank, passes_tcga_floor_in_type, delta_percentile_bg_corrected, tumor_median_rsem, adjnormal_median_rsem, n_tcga_adjacent_normals)
- crosstalk network:   /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/890af4bd-6338-410b-8744-e13de45260eb/v8086ab6e_crosstalk_network.csv  (crosstalk_network.csv; col passed_floor identifies the 81 floor-passing receptor subunits)
- corrected axes:      /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/c934590c-e27c-4e64-a2f0-8eef0c5b2e66/vfbcaf78d_axis_survival_corrected.csv  (axis_survival_corrected.csv; the 7 surviving axes)
- per-type expression: /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/51bd31c5-0e36-4709-ac73-f5ffb358e9ce/v07e2a94c_tcga_expr_by_tumor_type.csv  (tcga_expr_by_tumor_type.csv; 1436 genes x 30 tumour types, RSEM medians — use for T1.3 co-expression)

The receptor criterion as applied in the parent work is:
(direction_within_tcga == 'enriched') AND (passes_tcga_floor_in_type == True), evaluated per
(receptor, tumour type).

T1.1 BASE-RATE AUDIT. For each of the 81 floor-passing receptors compute the number of tumour types
it passes in, ligand-independent. Report the full distribution and the fraction passing in >=1 type.
KILL CRITERION: >=60% passing in >=1 type means the 7-axis count is UNINFORMATIVE as evidence for
PTN specifically. I measured 74.1% (60/81) when setting thresholds — CONFIRM OR REFUTE that under
the exact per-axis logic, and quantify how far PTN's 7 receptors sit above the base rate
(percentile of PTN's mean n_types within the all-receptor distribution).

T1.2 MATCHED-RANDOM NULL. Substitute PTN's 7 receptors with expression-decile-matched random
receptor sets of the same size (7), 1000 permutations, SEED 20260729. Statistic: total number of
(receptor, tumour type) enriched-and-floor-passing pairs, observed = 27 (verify this from the data).
Match deciles on the same expression measure the parent analysis used for its null — check the
crosstalk-score skill (skill(skill="crosstalk-score")) for matched_random_null() and reuse its
decile-matching logic so this null is comparable to the parent's.
PASS = observed exceeds the 95th percentile of the null. FAIL otherwise.

T1.3 FAMILY INDEPENDENCE. Spearman-correlate per-tumour-type expression of SDC1, SDC2, SDC3, SDC4,
NCL, ITGAV, PTPRZ1 across the 30 types; hierarchically cluster. Report the mean pairwise rho.
If mean pairwise rho >= 0.7 the 7 axes must be reported as ONE observation, not seven.

DELIVERABLES: receptor_base_rate_audit.csv, ptn_receptor_null.csv, ptn_receptor_correlation.csv.
State the Track 1 verdict: PASS requires T1.2 to pass AND T1.1 not to trigger the uninformative rule.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "track1_verdict": {
      "type": "string",
      "enum": [
        "PASS",
        "FAIL"
      ]
    },
    "t1_1": {
      "type": "object",
      "description": "base rate, n passing >=1 type, PTN mean n_types, PTN percentile in distribution, uninformative_triggered bool"
    },
    "t1_2": {
      "type": "object",
      "description": "observed pairs, null p95, null median, percentile, pass bool, seed"
    },
    "t1_3": {
      "type": "object",
      "description": "mean pairwise rho, cluster structure, single_block bool"
    },
    "artifacts": {
      "type": "object"
    },
    "caveats": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "track1_verdict",
    "t1_1",
    "t1_2",
    "t1_3",
    "artifacts"
  ]
}
```


---

## Wave 2 — Injury generality

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=36.2 min

### Task statement

```

ENVIRONMENT: use environment="vibemed" for every python cell. Install into vibemed via
manage_packages if something is missing.

THE PRE-SPECIFICATION IS BINDING. Read it FIRST: /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/c615b111-e38e-4459-a6b5-9e8191e0d33c/v1babf914_PRESPECIFICATION_PTN.md
Thresholds, directions, seeds and kill criteria are fixed there. You may NOT adjust a threshold,
swap a cohort/endpoint/statistic because the pre-specified one gave an unwelcome answer, or promote
a near-miss. Report FAILs with exactly the same prominence as PASSes.

CONTEXT. The parent screen found 0 surviving ligand-receptor axes. A corrected specificity
criterion (scoring cardiac distinctiveness by FAILURE INDUCTION rather than healthy-heart tissue
rank) produced 7 surviving axes, ALL with PTN (pleiotrophin) as the ligand: PTN->SDC1 (enriched in
6 tumour types), PTN->NCL (5), PTN->ITGAV (5), PTN->SDC3 (4), PTN->SDC4 (4), PTN->PTPRZ1 (2),
PTN->SDC2 (1) — 27 (receptor, tumour-type) pairs in total. PTN was admitted partly BY a relaxation
I introduced (heart rank 23 of 52 GTEx tissues, upper-half rather than first). YOUR JOB IS TO TRY
TO KILL THIS RESULT, not to support it. The pre-registered expectation is that PTN fails.

GLOBAL PROHIBITION P1: never place TCGA RSEM, GTEx TPM, RPKM and normalised counts on a common
numeric scale. Cross-dataset comparisons are ordinal only.

Print only computed values, not narration. Save tables with save_artifacts and report version_ids.
Do NOT build publication figures — synthesis does that.

TRACK 2 — Is PTN a cardiac signal or a generic tissue-injury signal?

Cardiac reference values (from this project, verified): PTN replicated induction min log2FC = 1.5023
across GSE116250 (discovery) and GSE141910 (validation); abundance 2.4-2.6% of GAPDH in failing
myocardium; best heart rank 23 of 52 GTEx tissues, top tissue pituitary.

T2.1 INDUCTION OUTSIDE THE HEART.
Use the published `geo-bulk-de` skill — load it with skill(skill="geo-bulk-de") and use its
functions; it handles GEO supplementary matrices, parses group labels from the series matrix
metadata rather than hardcoding them, and runs an empirical-Bayes moderated t test. It also has a
scale-check helper: some GEO series ship log2 values while labelled RAW, so check before analysing.
Process public bulk RNA-seq cohorts for NON-CARDIAC organ failure or fibrosis. Candidate contexts:
idiopathic pulmonary fibrosis, cirrhotic liver / NASH fibrosis, chronic kidney disease or renal
fibrosis. Choose series yourself; state GSE accessions and group definitions. Keep each cohort's
normalisation separate (P1) — compare effect sizes, never absolute values across cohorts.
AT LEAST TWO non-cardiac cohorts must be successfully processed for a verdict. If you cannot obtain
two, report the track INCONCLUSIVE and say which failed and why — do not pad with a marginal series.
KILL CRITERION: PTN is a generic injury signal if it is significantly induced (FDR<0.05, same
direction) in >=2 non-cardiac cohorts with a MEDIAN effect size >=50% of the cardiac effect
(i.e. median log2FC >= 0.751).
PASS = induction absent or substantially weaker outside the heart.
Report the full PTN row per cohort (log2FC, FDR, n per group), and as a positive control also report
a canonical fibrosis gene (e.g. COL1A1) per cohort to show the cohort's injury signal was detected
at all — a cohort where COL1A1 is flat is a cohort that did not capture fibrosis and must be flagged.

T2.2 CELL OF ORIGIN IN FAILING MYOCARDIUM.
Find a PUBLIC HUMAN heart-failure single-nucleus or single-cell RNA-seq dataset and determine which
cell population expresses PTN and whether that population's PTN rises in failure. Known candidates
include Chaffin et al. 2022 (Nature, DCM/HCM snRNA-seq) and Koenig et al. 2022; the Broad Single
Cell Portal, CELLxGENE and GEO are all plausible sources — CELLxGENE is allowlisted. Marker-level or
pseudobulk analysis is fine on 8 GiB; state which you used and do not attempt to load a full
multi-GB matrix into memory.
KILL CRITERION: if PTN is predominantly FIBROBLAST-derived and tracks a general fibrotic programme,
report PTN as a fibroblast wound-healing signal — this does not make the failing heart a distinctive
systemic source, and the caveat must reach the abstract. Cardiomyocyte-dominant expression would
instead support a cardiac-source reading.
If no suitable HUMAN dataset is obtainable, say so plainly. Do NOT substitute mouse data for a human
claim without labelling it as such.

DELIVERABLES: ptn_cross_tissue_induction.csv, ptn_cell_of_origin.csv.
State the Track 2 verdict.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "track2_verdict": {
      "type": "string",
      "enum": [
        "PASS",
        "FAIL",
        "INCONCLUSIVE"
      ]
    },
    "t2_1": {
      "type": "object",
      "description": "per-cohort GSE, tissue, PTN log2FC, FDR, n; COL1A1 control; median non-cardiac effect; cardiac reference; kill_triggered bool"
    },
    "t2_2": {
      "type": "object",
      "description": "dataset used, method (marker/pseudobulk), dominant cell type, fibroblast_dominant bool, whether PTN rises in failure in that population"
    },
    "cohorts_attempted": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "artifacts": {
      "type": "object"
    },
    "caveats": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "track2_verdict",
    "t2_1",
    "t2_2",
    "cohorts_attempted",
    "artifacts"
  ]
}
```


---

## Wave 2 — Plasma evidence

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=44.1 min

### Task statement

```

ENVIRONMENT: use environment="vibemed" for every python cell. Install into vibemed via
manage_packages if something is missing.

THE PRE-SPECIFICATION IS BINDING. Read it FIRST: /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/c615b111-e38e-4459-a6b5-9e8191e0d33c/v1babf914_PRESPECIFICATION_PTN.md
Thresholds, directions, seeds and kill criteria are fixed there. You may NOT adjust a threshold,
swap a cohort/endpoint/statistic because the pre-specified one gave an unwelcome answer, or promote
a near-miss. Report FAILs with exactly the same prominence as PASSes.

CONTEXT. The parent screen found 0 surviving ligand-receptor axes. A corrected specificity
criterion (scoring cardiac distinctiveness by FAILURE INDUCTION rather than healthy-heart tissue
rank) produced 7 surviving axes, ALL with PTN (pleiotrophin) as the ligand: PTN->SDC1 (enriched in
6 tumour types), PTN->NCL (5), PTN->ITGAV (5), PTN->SDC3 (4), PTN->SDC4 (4), PTN->PTPRZ1 (2),
PTN->SDC2 (1) — 27 (receptor, tumour-type) pairs in total. PTN was admitted partly BY a relaxation
I introduced (heart rank 23 of 52 GTEx tissues, upper-half rather than first). YOUR JOB IS TO TRY
TO KILL THIS RESULT, not to support it. The pre-registered expectation is that PTN fails.

GLOBAL PROHIBITION P1: never place TCGA RSEM, GTEx TPM, RPKM and normalised counts on a common
numeric scale. Cross-dataset comparisons are ordinal only.

Print only computed values, not narration. Save tables with save_artifacts and report version_ids.
Do NOT build publication figures — synthesis does that.

TRACK 3 — Has circulating pleiotrophin ever been measured in heart failure?

This is the question the whole endocrine premise rests on. The parent project already built a
verified literature base of 107 PubMed records using the literature connectors; use the same
approach here.

Discover the available literature tools FIRST with search_skills(prefix="mcp-") and load the
relevant connector doc, then call methods via host.mcp(...) from the repl tool. Note MCP calls only
work in the repl tool, never in python/r — write per-query loops in ONE repl cell and json.dump the
results to ./handoff/, then analyse in a python cell.

SEARCHES (run several, don't rely on one phrasing):
 a. circulating / plasma / serum pleiotrophin (PTN) in heart failure, cardiomyopathy, myocardial
    infarction, cardiac injury
 b. pleiotrophin as a biomarker in any cardiovascular context
 c. PTN as a secreted factor in cancer — tumour growth, angiogenesis, the syndecan and PTPRZ1
    receptors — to establish whether the cancer-side biology is real
 d. PTN in cardiac fibrosis or cardiac remodelling specifically

RULES:
- P4: every record must be RETRIEVED and VERIFIED. No citation written from memory. Record PMID/DOI,
  title, year, journal for each.
- Classify each record: supporting / contradicting / irrelevant to a CARDIAC PLASMA SOURCE.
- KILL CRITERION: if published data show plasma PTN is NOT elevated in heart failure, that directly
  contradicts the endocrine premise and kills the axis. Report it prominently if found.
- An ABSENCE of measurement is NOT evidence against — it is an informative gap, and if plasma PTN has
  never been measured in heart failure you must say so explicitly, since that is then the single
  gap the entire axis rests on.
- Do not overstate: a study measuring PTN in tumour tissue is not a study measuring plasma PTN in
  heart failure.

DELIVERABLES: ptn_literature.csv (one row per verified record with identifier, title, year, journal,
classification, and what it actually measured), and a short ptn_literature_summary.md.
State plainly whether plasma PTN in heart failure has ever been measured.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "plasma_ptn_measured_in_hf": {
      "type": "boolean"
    },
    "contradiction_found": {
      "type": "boolean",
      "description": "published evidence that plasma PTN is NOT elevated in HF"
    },
    "n_records_verified": {
      "type": "integer"
    },
    "n_supporting": {
      "type": "integer"
    },
    "n_contradicting": {
      "type": "integer"
    },
    "cancer_biology_established": {
      "type": "boolean",
      "description": "is PTN's secreted role in cancer well supported"
    },
    "key_records": {
      "type": "array",
      "items": {
        "type": "object"
      },
      "description": "identifier, year, what it measured, classification"
    },
    "artifacts": {
      "type": "object"
    },
    "caveats": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "plasma_ptn_measured_in_hf",
    "contradiction_found",
    "n_records_verified",
    "cancer_biology_established",
    "artifacts"
  ]
}
```


---

## Wave 3 — Proteome

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=1880.3 min

### Task statement

```
Read the pre-specification /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/165bf2e4-7ab3-41ed-a408-421bc1508de6/vf15b5995_PRESPECIFICATION_REPLICATION.md before anything else. Sections 0, 2 and 3 govern you.

STEP 1 - Author and publish a `pdc-proteome` skill.
Write SKILL.md + kernel.py for the PDC GraphQL endpoint using the verified facts above. The
kernel must include, as functions that RAISE on violation:
  - a testability guard implementing pre-spec section 2.1 (gene present in matrix; >=10 non-missing
    aliquots per arm; >=60% non-missing per arm). It must return UNTESTABLE rather than a p-value,
    and must refuse to be talked out of it.
  - a plex-confound checker: cross-tabulate plex x sample_type, chi-square for uniformity, and
    report how many case-pairs share a plex. This is the fatal-if-present risk.
  - an arm-dependent-missingness test (Fisher exact on missing x arm).
  - a genome-wide base-rate function per pre-spec 2.4.
Document in SKILL.md the traps: introspection disabled; acceptDUA accepted on some fields and
rejected on others; quantDataMatrix takes no subfield selection and returns list-of-lists;
header columns are case_id:aliquot_submitter_id; GAPDH is NOT arm-invariant in protein
(0.295 tumour vs -0.427 normal in ccRCC) so the transcript's 0.005xGAPDH floor MUST NOT be
ported; absent gene != depleted gene.
Publish it with host.skills.* (load skill 'customize' for the API).

STEP 2 - Download. For the CPTAC Proteome studies covering the tumour types the screen ranked -
ccRCC, OV, BRCA, PAAD, UCEC, LUAD, HNSC, GBM - pull log2_ratio quantDataMatrix and
biospecimenPerStudy, and studyExperimentalDesign for plex structure. Pick the largest
discovery-cohort Proteome study per tumour type; record which pdc_study_id you chose and why.
Run the plex-confound check on EVERY study - do not assume ccRCC's clean design generalises.
Save: cptac_study_selection.csv, cptac_aliquot_map.csv, cptac_plex_confound.csv.

STEP 3 - NPR1 (pre-spec 2.2-2.5). Paired-within-plex Wilcoxon signed-rank on tumour-normal
differences, BH-FDR across types, median paired difference as effect size. Report per type:
n pairs, effect, FDR, testable/untestable with reason, and whether pairing was within-plex
(batch_free) or blocked (batch_corrected). Run ALL FOUR sensitivity analyses in 2.5.
COMPUTE THE BASE RATE FIRST (2.4) and state NPR1's genome-wide percentile - if the base rate
exceeds 40%, NPR1's result is UNINFORMATIVE no matter what its p-value is.
Then apply 2.3: REPLICATES (<0, FDR<0.05, >=5 of 7 blind types) / REVERSES (>0, FDR<0.05,
>=3 of 7 blind types) / INCONCLUSIVE. The 7 BLIND types are OV, BRCA, PAAD, UCEC, LUAD, HNSC, GBM.
ccRCC is UNBLINDED (I already saw -0.153) - report it separately and give the headline verdict
both with and without it. Apply the 2.6 kill criterion: <4 of 7 testable means ABANDON as
underpowered.
Save: npr1_protein_tumor_normal.csv, npr1_protein_sensitivity.csv, protein_base_rate.csv.

STEP 4 - PTN panel (pre-spec 3). SDC1, SDC2, SDC3, SDC4, NCL, ITGAV, PTPRZ1, NPR2, NPR3. Same
testability and pairing rules. Abundance-decile-matched random null, 300 draws, seed 20260727 -
an unmatched null is not acceptable. Apply the 3 kill criterion (95th percentile of matched null).
PTN itself: report if quantified, but state that CPTAC measures tumour tissue not plasma, so its
absence is uninformative about a circulating cardiac ligand.
Save: ptn_panel_protein.csv, ptn_panel_protein_null.csv.

Return in structured_output the verdict for each track, the numbers behind it, every artifact id,
and any pre-spec criterion you could not apply (with the reason).
```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "npr1_verdict": {
      "type": "string",
      "enum": [
        "REPLICATES",
        "REVERSES",
        "INCONCLUSIVE",
        "UNINFORMATIVE",
        "ABANDONED_UNDERPOWERED"
      ]
    },
    "npr1_verdict_excluding_ccrcc": {
      "type": "string"
    },
    "npr1_per_type": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "tumor_type": {
            "type": "string"
          },
          "pdc_study_id": {
            "type": "string"
          },
          "n_pairs": {
            "type": "integer"
          },
          "median_paired_diff": {
            "type": [
              "number",
              "null"
            ]
          },
          "fdr": {
            "type": [
              "number",
              "null"
            ]
          },
          "testable": {
            "type": "boolean"
          },
          "reason_if_untestable": {
            "type": "string"
          },
          "batch_handling": {
            "type": "string"
          },
          "blind": {
            "type": "boolean"
          }
        }
      }
    },
    "base_rate_fraction_depleted": {
      "type": "number"
    },
    "base_rate_exceeds_40pct": {
      "type": "boolean"
    },
    "npr1_genomewide_percentile": {
      "type": [
        "number",
        "null"
      ]
    },
    "n_blind_types_testable": {
      "type": "integer"
    },
    "plex_confound_any_study": {
      "type": "string"
    },
    "ptn_panel_verdict": {
      "type": "string"
    },
    "ptn_panel_observed_vs_null": {
      "type": "string"
    },
    "sensitivity_agreement": {
      "type": "string"
    },
    "overturns_paper": {
      "type": "boolean"
    },
    "artifact_ids": {
      "type": "object"
    },
    "prespec_criteria_not_applied": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "skill_published": {
      "type": "string"
    }
  },
  "required": [
    "npr1_verdict",
    "npr1_per_type",
    "base_rate_fraction_depleted",
    "base_rate_exceeds_40pct",
    "n_blind_types_testable",
    "overturns_paper",
    "artifact_ids",
    "skill_published"
  ]
}
```


---

## Wave 3 — Platforms

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=1870.7 min

### Task statement

```
Read the pre-specification /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/165bf2e4-7ab3-41ed-a408-421bc1508de6/vf15b5995_PRESPECIFICATION_REPLICATION.md before anything else. Section 4 governs you.

STEP 1 - Extend the `geo-bulk-de` skill for microarray.
Load the existing skill (host.skills.read('geo-bulk-de', path=...)), then add:
  - Affymetrix series-matrix parsing (probe-level, GPL96 and GPL11532).
  - Probe-to-gene collapsing by MAXIMUM MEAN INTENSITY across samples. Add an assertion that
    RAISES on an ambiguous mapping (a probe mapping to multiple genes) rather than silently
    taking the first - silent first-pick is exactly the class of bug this project has been
    catching.
  - Log-scale DETECTION rather than assumption: infer whether values are already log2 and record
    which branch fired. Do not assume.
  - Keep the kernel pure stdlib + numpy/pandas/scipy. Preserve every existing function and its
    behaviour - this skill is published and has downstream users; do not break its API.
Publish the updated skill (load skill 'customize' for the API) and say what you changed.

STEP 2 - Run DE on GSE5406, GSE1869, GSE57345.
Parse failing vs non-failing labels FROM the series matrix characteristics. Do NOT hardcode group
membership. GSE5406 and GSE1869 both distinguish ischaemic from non-ischaemic cardiomyopathy -
pool them as "failing" for the primary contrast, and report the ischaemic-only contrast separately.
Run the label-permutation control on every cohort: pre-spec 4.1 requires 0 significant genes at
FDR<0.05. If a permutation control yields hits, STOP and report that - it means the pipeline is
broken, and no downstream number can be trusted.
Verify the four canonical markers in each cohort and report their direction and FDR.
Save: hf_de_GSE5406.csv, hf_de_GSE1869.csv, hf_de_GSE57345.csv, hf_permutation_controls.csv.

STEP 3 - Test platform-independence (pre-spec 4.1-4.3).
Against the existing replicated signature (artifact hf_signature_replicated.csv - find it with
host.artifacts(filename='hf_signature_replicated.csv', exact=True)), compute for each new cohort:
  - directional concordance among genes significant in both
  - log2FC Spearman rho, reported TWICE and never conflated: across ALL shared genes, and among
    genes significant in both. Conflating these two was a prior audit finding in this project.
  - the specific induction of NPPA, NPPB and PTN.
Apply the 4.1 thresholds (concordance >=75%, rho >=0.40 among sig-in-both, NPPA and NPPB up at
FDR<0.05 in >=2 of 3 cohorts) and the 4.3 kill criterion (below threshold in >=2 of 3 cohorts =
signature is PLATFORM-DEPENDENT, and downstream claims must be re-scoped).
Also answer directly: did the shared-platform weakness of the original pair materially affect the
signature? Quantify rather than assert.
Save: hf_cross_platform_replication.csv, hf_platform_independence_summary.json.

Return in structured_output the verdict, all correlation values kept separate, and every artifact id.
```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "signature_verdict": {
      "type": "string",
      "enum": [
        "PLATFORM_INDEPENDENT",
        "PLATFORM_DEPENDENT",
        "INCONCLUSIVE"
      ]
    },
    "per_cohort": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "gse": {
            "type": "string"
          },
          "platform": {
            "type": "string"
          },
          "technology": {
            "type": "string"
          },
          "n_failing": {
            "type": "integer"
          },
          "n_nonfailing": {
            "type": "integer"
          },
          "n_sig_fdr05": {
            "type": "integer"
          },
          "concordance_sig_in_both": {
            "type": [
              "number",
              "null"
            ]
          },
          "rho_all_shared": {
            "type": [
              "number",
              "null"
            ]
          },
          "rho_sig_in_both": {
            "type": [
              "number",
              "null"
            ]
          },
          "n_shared_genes": {
            "type": [
              "integer",
              "null"
            ]
          },
          "permutation_control_n_sig": {
            "type": "integer"
          },
          "NPPA_log2fc": {
            "type": [
              "number",
              "null"
            ]
          },
          "NPPA_fdr": {
            "type": [
              "number",
              "null"
            ]
          },
          "NPPB_log2fc": {
            "type": [
              "number",
              "null"
            ]
          },
          "NPPB_fdr": {
            "type": [
              "number",
              "null"
            ]
          },
          "PTN_log2fc": {
            "type": [
              "number",
              "null"
            ]
          },
          "PTN_fdr": {
            "type": [
              "number",
              "null"
            ]
          },
          "log_scale_branch": {
            "type": "string"
          },
          "passes_thresholds": {
            "type": "boolean"
          }
        }
      }
    },
    "canonical_markers_correct": {
      "type": "string"
    },
    "shared_platform_effect": {
      "type": "string"
    },
    "skill_updated": {
      "type": "string"
    },
    "artifact_ids": {
      "type": "object"
    },
    "prespec_criteria_not_applied": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "signature_verdict",
    "per_cohort",
    "canonical_markers_correct",
    "shared_platform_effect",
    "artifact_ids",
    "skill_updated"
  ]
}
```


---

## Wave 4 — Prior art

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=70.7 min

### Task statement

```

CONTEXT — project background you must respect.

This project has already done extensive work on reverse cardio-oncology (the hypothesis that
cardiac disease promotes cancer). Two doors are already closed by prior work IN THIS PROJECT:

1. The endocrine ligand-receptor crosstalk hypothesis was tested and WITHDRAWN after adversarial
   testing. Seven candidate axes (all pleiotrophin/PTN) failed all three pre-specified kill
   criteria: 74% of candidate receptors passed the receptor test with no ligand involved; PTN's
   receptor set UNDERperformed expression-matched random draws (p=0.688); PTN is a fibroblast
   wound-healing transcript (71% of counts) with no per-cell cardiac induction; and circulating PTN
   has never been measured in a heart-failure cohort.
2. An attempt to run Mendelian randomization with HEART FAILURE as the exposure was BLOCKED at a
   feasibility gate. The critical finding: the 11 independent genome-wide-significant HF loci are
   canonical SHARED-RISK-FACTOR loci — PITX2 (atrial fibrillation), 9p21 (coronary artery disease),
   LPA, SORT1, FTO (BMI). A positive-control MR (HF -> NPPB protein) also failed (p=0.209).

That second finding is the seed of the current work: if HF's genetic instruments are actually
AF/CAD/BMI loci, then published MR studies using HF as an exposure for cancer may be estimating
the effect of those shared risk factors, not of heart failure.

TOOLS. MCP calls ONLY work in the `repl` tool via host.mcp(server, method, **kwargs) — never in
`python`/`r`. Pass data between kernels via ./handoff/*.json files. Relevant servers:
  - "pubmed": search_articles, get_article_metadata, convert_article_ids, get_full_text_article
  - "human-genetics": gwas_search_traits, gwas_associations_for_trait, gwas_associations_for_variant,
    gwas_get_study
Load the mcp-pubmed / mcp-human-genetics skill docs if you need exact method schemas.
Use fetch_article_fulltext(doi=...) for non-PMC papers; it saves to articles/ for read_file.

RULES.
- Ground every claim in a retrieved record. Never write a PMID, DOI, effect estimate or sample size
  from memory — copy it from the tool result. If you cannot retrieve something, say so explicitly.
- Do NOT fabricate data if something is inaccessible. Record the gap.
- Save all deliverables with save_artifacts and report version_ids.

YOUR TASK — Track: Prior art. Establish the prior-art boundary so later phases cannot duplicate
published work. Two plan steps.

STEP 1 "Establish what is already published".
Retrieve and read IN FULL (not just abstracts) these papers, recording for each: the exact trait
pairs analysed, the methods used, the data sources, the sample sizes, and the headline findings.

Already-published cross-trait genetic analyses (these define what is DONE):
  a) doi 10.1161/ATVBAHA.125.322433 — PMID 40874306 — genome-wide pleiotropy analysis, coronary
     artery disease vs breast/colorectal/lung/prostate cancer. Reported 60 colocalized loci, 43
     novel, 35 lead SNPs with OPPOSING effects. Get the full text; record precisely which trait
     pairs and which methods (pleiotropy, colocalization, multitrait colocalization, PheWAS).
  b) doi 10.3389/fgene.2025.1450259 — PMID 40201568 — breast cancer vs atrial fibrillation shared
     architecture (LDSC, HDL, MiXeR, LAVA, GWAS-PW, CPASSOC, colocalization, TWAS, MR).

Foundational reverse cardio-oncology primary papers, and the expert critique of one of them:
  c) PMID 29459363 — Meijers, Circulation 2018, "Heart Failure Stimulates Tumor Growth by
     Circulating Factors" (heterotopic heart transplantation).
  d) PMID 32475164 — Avraham, Circulation 2020, transverse aortic constriction promotes tumour
     growth and metastasis.
  e) PMID 32661390 — Koelwyn, Nat Med 2020, "Myocardial infarction accelerates breast cancer via
     innate immune reprogramming".
  f) doi 10.1161/CIRCULATIONAHA.120.048655 — Moslehi, Zhang & Moore, Circulation 2020 editorial
     "Crosstalk Between the Heart and Cancer: Beyond Drug Toxicity" — this is an EXPERT CRITICAL
     APPRAISAL of (d) and is load-bearing for a later benchmark. Extract its enumerated criticisms
     VERBATIM into a structured table, one row per distinct criticism, with the exact wording. A
     later phase scores an automated appraiser against this list, so fidelity matters more than
     summary. Save as `editorial_criticisms.csv` with columns:
     criticism_id, target_pmid, criticism_verbatim, category (one of: dose/physiological-plausibility,
     model-system-limitation, missing-causal-experiment, immune-confound, human-data-limitation,
     epidemiological-bias, scope-of-inference), detectable_from_fulltext (yes/no — could an
     appraiser reading ONLY the Avraham full text plausibly flag this?).
  g) PMID 38139195 — the later paper showing periostin-knockout mice STILL showed
     remodelling-dependent tumour promotion, and integrin beta-1 loss in cancer cells abolished it.
     This is the overturned-claim fixture. Record exactly what was claimed vs what was overturned.

STEP 2 "Write the prior-art boundary document".
Write `PRIOR_ART_BOUNDARY.md` stating, in prose with inline DOI links:
  - Which cross-trait genetic analyses are ALREADY DONE and must not be duplicated (be specific
    about trait pairs and methods, so a later phase can check itself against this).
  - What the foundational animal-model evidence does and does not establish in humans.
  - Which questions remain genuinely OPEN. Note specifically that cardiomyopathy vs cancer has not
    been covered by the published pleiotropy work.
  - An explicit "do not attempt" list.
Also save `prior_art_papers.csv` (pmid, doi, year, journal, first_author, design, traits_analysed,
methods, key_finding, fulltext_retrieved).

Deliverables: PRIOR_ART_BOUNDARY.md, prior_art_papers.csv, editorial_criticisms.csv.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "summary": {
      "type": "string"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "filename": {
            "type": "string"
          },
          "version_id": {
            "type": "string"
          },
          "description": {
            "type": "string"
          }
        },
        "required": [
          "filename",
          "version_id"
        ]
      }
    },
    "key_numbers": {
      "type": "object",
      "description": "load-bearing counts/values computed, as name->value"
    },
    "gaps": {
      "type": "array",
      "items": {
        "type": "string"
      },
      "description": "what could not be retrieved or computed, and why"
    }
  },
  "required": [
    "summary",
    "artifacts",
    "key_numbers",
    "gaps"
  ]
}
```


---

## Wave 4 — MR corpus

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=77.6 min

### Task statement

```

CONTEXT — project background you must respect.

This project has already done extensive work on reverse cardio-oncology (the hypothesis that
cardiac disease promotes cancer). Two doors are already closed by prior work IN THIS PROJECT:

1. The endocrine ligand-receptor crosstalk hypothesis was tested and WITHDRAWN after adversarial
   testing. Seven candidate axes (all pleiotrophin/PTN) failed all three pre-specified kill
   criteria: 74% of candidate receptors passed the receptor test with no ligand involved; PTN's
   receptor set UNDERperformed expression-matched random draws (p=0.688); PTN is a fibroblast
   wound-healing transcript (71% of counts) with no per-cell cardiac induction; and circulating PTN
   has never been measured in a heart-failure cohort.
2. An attempt to run Mendelian randomization with HEART FAILURE as the exposure was BLOCKED at a
   feasibility gate. The critical finding: the 11 independent genome-wide-significant HF loci are
   canonical SHARED-RISK-FACTOR loci — PITX2 (atrial fibrillation), 9p21 (coronary artery disease),
   LPA, SORT1, FTO (BMI). A positive-control MR (HF -> NPPB protein) also failed (p=0.209).

That second finding is the seed of the current work: if HF's genetic instruments are actually
AF/CAD/BMI loci, then published MR studies using HF as an exposure for cancer may be estimating
the effect of those shared risk factors, not of heart failure.

TOOLS. MCP calls ONLY work in the `repl` tool via host.mcp(server, method, **kwargs) — never in
`python`/`r`. Pass data between kernels via ./handoff/*.json files. Relevant servers:
  - "pubmed": search_articles, get_article_metadata, convert_article_ids, get_full_text_article
  - "human-genetics": gwas_search_traits, gwas_associations_for_trait, gwas_associations_for_variant,
    gwas_get_study
Load the mcp-pubmed / mcp-human-genetics skill docs if you need exact method schemas.
Use fetch_article_fulltext(doi=...) for non-PMC papers; it saves to articles/ for read_file.

RULES.
- Ground every claim in a retrieved record. Never write a PMID, DOI, effect estimate or sample size
  from memory — copy it from the tool result. If you cannot retrieve something, say so explicitly.
- Do NOT fabricate data if something is inaccessible. Record the gap.
- Save all deliverables with save_artifacts and report version_ids.

YOUR TASK — Track: MR corpus. Assemble and characterise the corpus of Mendelian randomization
papers that test cardiac disease as an EXPOSURE for cancer. Three plan steps.

Background measured already: PubMed counts for `Mendelian randomization AND "heart failure" AND
cancer` by year are 2018:1, 2019:4, 2020:3, 2021:8, 2022:3, 2023:4, 2024:20, 2025:16, 2026:5
(66 total). Verify these yourself and treat the sharp post-2023 rise as a feature to document.

STEP 1 "Assemble the reverse cardio-oncology MR corpus".
Search systematically. Use SEVERAL query formulations (heart failure, cardiomyopathy, atrial
fibrillation, coronary artery disease, myocardial infarction, each vs cancer / neoplasm /
site-specific cancers) so the corpus is not an artifact of one query. Record the query log.
INCLUSION: two-sample or one-sample MR where a CARDIAC/cardiovascular disease or trait is the
exposure and a CANCER is an outcome. Bidirectional studies count — record both directions.
EXCLUSION: cardiotoxicity-of-cancer-therapy papers (that is forward cardio-oncology, not this),
and MR studies where neither exposure nor outcome is cardiac.
For each included paper get the PMID, DOI, year, journal, and attempt full text: convert PMID to
PMCID and use get_full_text_article; where there is no PMC record try fetch_article_fulltext(doi).
Record retrieval status honestly per paper. Save `mr_corpus.csv` and `mr_corpus_query_log.csv`.

STEP 2 "Extract instruments, exposures and claims".
For every paper WITH full text, extract into `mr_extraction.csv` (one row per exposure-outcome
pair; a bidirectional paper gives >=2 rows):
  pmid, doi, direction (cardiac_to_cancer / cancer_to_cardiac), exposure_trait,
  exposure_gwas_source (consortium/biobank and accession if stated), exposure_sample_size,
  outcome_trait, outcome_gwas_source, outcome_sample_size, instrument_pvalue_threshold,
  instrument_r2_clumping, n_instruments, primary_method, effect_estimate, ci_lower, ci_upper,
  p_value, reported_conclusion_verbatim (the causal claim from the abstract/conclusions, quoted),
  instrument_snps_listed (yes/no/partial — did the paper or supplement give the rsIDs?),
  instrument_snps (semicolon-separated rsIDs if recoverable).
Quote effect estimates exactly as printed. Where a field is absent from the paper, write NOT_REPORTED
— that absence is itself the finding for step 3. Never infer a value.

STEP 3 "Flag reporting and design gaps".
Build `mr_reporting_gaps.csv` — one row per paper, columns as yes/no/NOT_REPORTED:
  reported_instrument_strength (F-statistic or R2), reported_variance_explained,
  tested_sample_overlap (exposure and outcome GWAS drawing on the same cohort, e.g. both UK
  Biobank — a major bias source), reported_heterogeneity (Cochran Q), reported_egger_intercept,
  reported_leave_one_out, reported_steiger_or_reverse_causation, reported_power_calculation,
  corrected_for_multiple_testing, cites_STROBE_MR, exposure_is_binary_disease (relevant: MR with a
  binary disease exposure has known interpretational problems).
Then compute and save `mr_reporting_summary.csv`: for each column, the count and percentage of the
corpus satisfying it. The published comparator (do NOT recompute it, just record it for later use)
is doi 10.1093/ije/dyac074: across 87 MR-Base studies, mean 53% of reporting items satisfied, 44%
sufficient for the relevance assumption, 31% for independence.

Also write a short `CORPUS_NOTES.md`: how the corpus was built, what the growth curve looks like,
which journals dominate, how many papers had no accessible full text, and any patterns you noticed
(e.g. repeated boilerplate methods, same GWAS reused, suspiciously uniform structure).

Deliverables: mr_corpus.csv, mr_corpus_query_log.csv, mr_extraction.csv, mr_reporting_gaps.csv,
mr_reporting_summary.csv, CORPUS_NOTES.md.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "summary": {
      "type": "string"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "filename": {
            "type": "string"
          },
          "version_id": {
            "type": "string"
          },
          "description": {
            "type": "string"
          }
        },
        "required": [
          "filename",
          "version_id"
        ]
      }
    },
    "key_numbers": {
      "type": "object",
      "description": "load-bearing counts/values computed, as name->value"
    },
    "gaps": {
      "type": "array",
      "items": {
        "type": "string"
      },
      "description": "what could not be retrieved or computed, and why"
    }
  },
  "required": [
    "summary",
    "artifacts",
    "key_numbers",
    "gaps"
  ]
}
```


---

## Wave 4 — Claim audit

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=81.2 min

### Task statement

```

CONTEXT.

This project audits the reverse cardio-oncology literature (the hypothesis that cardiac disease
promotes cancer). Phase 0 is complete. Its findings that constrain your work:

- A 33-paper corpus of Mendelian randomization (MR) studies using a cardiac disease as EXPOSURE and
  cancer as OUTCOME was assembled; 29 have accessible full text; 192 exposure-outcome rows extracted.
- CRITICAL: instrument rsIDs are unrecoverable for 186 of 192 rows (96.9%) and the instrument count
  is NOT_REPORTED for 130 of 192. Published papers cite supplementary tables absent from PMC text.
  Therefore instruments must be RE-DERIVED from the named GWAS accessions, not read off the papers.
- The dominant heart-failure exposure accession in the corpus is GCST009541, verified as Shah 2020
  (PMID 31919418), 47,309 European cases / 930,014 controls.
- Prior work IN THIS PROJECT found heart failure's genome-wide-significant instruments are canonical
  SHARED-RISK-FACTOR loci: PITX2 (atrial fibrillation), 9p21 (coronary artery disease), LPA, SORT1,
  FTO (BMI). A positive control (HF -> NPPB protein) also failed at p=0.209. That is the hypothesis
  under test: published HF-exposure MR may be estimating shared risk factors, not heart failure.
- Published comparator for MR reporting quality (RECORD it, never recompute): doi 10.1093/ije/dyac074,
  87 MR-Base studies, mean 53% of reporting items satisfied, 44% sufficient for the relevance
  assumption, 31% for independence. This corpus scored 49.7%, 75.9% and 10.3% respectively.

PHASE 0 ARTIFACTS you may read with read_file(version_id=...) or host.artifact_path(vid):
  PRIOR_ART_BOUNDARY.md    541bf3fa-db40-4338-9e6a-24de912fb9f9   (read this FIRST; it has a
                                                                   do-not-attempt list you must obey)
  editorial_criticisms.csv 115007c7-85a5-4a7c-9f1c-cedba48dda5e
  prior_art_papers.csv     b97a6c2d-b19a-45c7-a7e5-c7e174beb61f
  mr_corpus.csv            587aff0f-f85d-45a1-8474-34acebcaec01
  mr_extraction.csv        644688bd-5ed5-4783-866b-f8e481c5e14c
  mr_reporting_gaps.csv    bc8b4de2-dc68-46e0-b562-abc2ab18cb27
  mr_reporting_summary.csv d931e9eb-6e8a-42c9-8f68-944ffdc35f71
  CORPUS_NOTES.md          e79c496a-e4d9-420f-81ea-40af5ae202d5

TOOLS. MCP calls work ONLY in the `repl` tool via host.mcp(server, method, **kwargs) — never in
`python`/`r`. Pass data between kernels via ./handoff/*.json. Servers: "pubmed" (search_articles,
get_article_metadata, convert_article_ids, get_full_text_article), "human-genetics"
(gwas_search_traits, gwas_associations_for_trait, gwas_associations_for_variant, gwas_get_study).
Full texts are NOT in your workspace — re-retrieve from PMC using the PMCIDs in mr_corpus.csv.

SKILL AUTHORING. You are publishing a real reusable skill. Load skill('customize') for the
host.skills.* API. Requirements, non-negotiable:
  - kernel.py at the skill root: pure numpy/pandas/scipy/statsmodels — NO host.* calls, NO network,
    NO MCP inside kernel.py, so it runs standalone. Functions + imports + literal constants only;
    no top-level classes or decorators (wrap those in a factory function).
  - Guardrails are ASSERTIONS that RAISE, not prose a reader can skip. This is the house convention
    from this project's crosstalk-score skill: require_* validators raise AssertionError by default
    and return an audit table so anything dropped is disclosed rather than silently removed.
  - SKILL.md must state what the skill does NOT establish, and must credit prior art rather than
    implying novelty where none exists.
  - Ship regression tests (test_*.py) that reproduce real measured numbers, and confirm they pass.
  - Publish with host.skills.publish(name) and report the skill name.

RULES.
- Never write a PMID, DOI, effect estimate, sample size or accession from memory — copy it from a
  tool result. If you cannot retrieve something, record the gap explicitly.
- Never fabricate data to fill a gap.
- Do not duplicate anything on the PRIOR_ART_BOUNDARY.md do-not-attempt list.
- Save deliverables with save_artifacts and report version_ids.

YOUR TASK — Track: Claim audit. Build, benchmark and publish the `claim-audit` skill. Five steps.

STEP 1 "Build claim-audit kernel".
Implement kernel.py providing full-text appraisal. Required capabilities:
  - segment_sections(text) -> dict: split an article into abstract, introduction, methods,
    statistical_analysis, results, discussion, limitations, conclusions. Real PMC text is messy;
    handle numbered headings, ALL-CAPS headings, and missing sections. Return what was found and
    record what was absent — never silently return an empty section as if it were empty in the paper.
  - score_strobe_mr(sections) -> DataFrame: score a paper against the STROBE-MR checklist
    (doi 10.1001/jama.2021.18236), an EXISTING published standard — you are encoding it, not
    inventing a rubric. Score the three core instrumental-variable assumptions explicitly and
    separately: RELEVANCE (instrument associated with exposure: F-statistic, R2, variance explained),
    INDEPENDENCE/exchangeability (no instrument-outcome confounding: pleiotropy tests, Egger
    intercept, Cochran Q, sample overlap), EXCLUSION RESTRICTION (no effect except through the
    exposure). For every item return the item, a yes/no/not_reported verdict, AND the verbatim
    evidence span from the text that justifies it — a score with no quotable evidence is not usable.
  - detect_causal_overclaim(sections) -> DataFrame: find causal language ("causal effect",
    "causes", "leads to", "protective effect", "confers risk") in abstract/conclusions and pair it
    with whether the corresponding validity checks were actually reported. The output claim is a
    MISMATCH between claim strength and evidence, not a judgement of truth.
  - design_strength(design_string) -> ordered score: place a design in the evidence hierarchy
    (RCT > MR with validated instruments > prospective cohort with lag > cohort > cross-sectional >
    animal model > in vitro).
Sanity-check on real retrieved text as you go, not on synthetic strings.

STEP 2 "Encode computability triage as validators".
Implement triage classifying each claim into exactly one of:
  - COMPUTABLE_FROM_PUBLIC_DATA
  - REQUIRES_UNAVAILABLE_MEASUREMENT (the measurement is possible in principle but has never been
    made, so no public dataset can settle it)
  - UNTESTABLE_IN_PRINCIPLE (no observational or public-data design could settle it)
Make require_computable() RAISE by default when a claim in the UNTESTABLE or UNAVAILABLE class is
about to be treated as testable, returning an audit table of what was flagged. Seed the triage with
measurement gaps this project already established: circulating pleiotrophin has never been measured
in a heart-failure cohort against controls; human plasma periostin concentrations in heart failure
patients are unmeasured; no human study demonstrates engagement of a heart-failure-derived ligand
with its receptor on a tumour; tissue mRNA cannot be converted to a plasma concentration without
translation efficiency, secretion rate, clearance and volume of distribution.

STEP 3 "Benchmark against expert appraisal".
This is the step that decides whether the skill is usable. Score the appraiser against the 14
verbatim criticisms in editorial_criticisms.csv (version 115007c7-85a5-4a7c-9f1c-cedba48dda5e),
which are an INDEPENDENT expert critique by Moslehi, Zhang and Moore of Avraham 2020 (PMID 32475164).
  IMPORTANT CONSTRAINT: the Avraham full text is PAYWALLED and could not be retrieved (all five
  fetch routes failed). You therefore CANNOT score recall against the full text. Do this instead,
  and be explicit about it:
   (a) Score the appraiser on the Avraham material that IS available — the PubMed abstract — and
       report recall against the subset of criticisms that are in principle detectable from an
       abstract alone. State plainly that this is a harder and narrower test than the full-text test.
   (b) Additionally run the appraiser on a paper whose full text IS available and which the same
       editorial context applies to: Koelwyn 2020 (PMID 32661390, doi 10.1038/s41591-020-0964-7),
       retrievable from PMC. Report what it flags.
   (c) Report BOTH recall and FALSE POSITIVES. An appraiser that flags everything has recall 1.0 and
       is useless. Report precision alongside recall.
  Do NOT edit the wording in editorial_criticisms.csv — it is the fixture.
  Save `appraiser_benchmark.csv` (criticism_id, category, detectable_from_abstract, flagged_by_appraiser,
  matching_appraiser_output, verdict) and `appraiser_benchmark_summary.csv` with recall, precision and
  the denominators used.

STEP 4 "Test the overturned-claim case".
Two fixtures.
  (a) OVERTURNED: Avraham 2020 claimed periostin is the necessary circulating mediator of
      remodelling-driven tumour promotion. Langier Goncalves 2023 (PMID 38139195,
      doi 10.3390/ijms242417367) overturned it — in POSTN-knockout mice TAC STILL promoted tumour
      growth, and causality moved to cancer-cell Integrin beta-1. Confirm the appraiser flags the
      original claim as establishing SUFFICIENCY but not NECESSITY (an in-vitro depletion plus
      exogenous-addition experiment cannot establish necessity in vivo). This is a specific,
      checkable logical gap — verify the appraiser catches it.
  (b) RETRACTED: these PMIDs are retracted and must never be cited as support — 23903826 and
      26549838 (delta-opioid-receptor hepatocellular papers), 28601389 (periostin/tumour
      microenvironment, retraction notice PMID 36990904). Confirm the appraiser refuses to treat a
      retracted record as supporting evidence. Verify retraction status via PubMed rather than
      trusting this list.
  Save `overturned_claim_fixtures.csv`.

STEP 5 "Add regression tests and publish".
Freeze the benchmark and fixture results as regression tests (test_claim_audit.py) so a later change
cannot silently lower recall. Then apply the finished appraiser to the 29 full-text corpus papers and
save `corpus_claim_audit.csv` (per paper: STROBE-MR items satisfied, the three IV assumptions,
causal-overclaim flags, computability class of the headline claim). Compare the corpus against the
published base rates (53% / 44% / 31%, doi 10.1093/ije/dyac074) so a poor score is read against the
field norm, not presented as unique to cardio-oncology.
Publish the skill as `claim-audit`. SKILL.md must document: the failure modes encoded, the measured
recall AND precision with their denominators, the paywall limitation on the benchmark, and an explicit
statement that the appraiser flags WEAK SUPPORT and REPORTING GAPS — it does not adjudicate whether a
claim is true.

Deliverables: published `claim-audit` skill (SKILL.md + kernel.py + tests), appraiser_benchmark.csv,
appraiser_benchmark_summary.csv, overturned_claim_fixtures.csv, corpus_claim_audit.csv.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "summary": {
      "type": "string"
    },
    "skill_published": {
      "type": "string"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "filename": {
            "type": "string"
          },
          "version_id": {
            "type": "string"
          },
          "description": {
            "type": "string"
          }
        },
        "required": [
          "filename",
          "version_id"
        ]
      }
    },
    "key_numbers": {
      "type": "object"
    },
    "tests_passing": {
      "type": "string"
    },
    "gaps": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "summary",
    "skill_published",
    "artifacts",
    "key_numbers",
    "tests_passing",
    "gaps"
  ]
}
```


---

## Wave 4 — Instrument audit

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=58.6 min

### Task statement

```

CONTEXT.

This project audits the reverse cardio-oncology literature (the hypothesis that cardiac disease
promotes cancer). Phase 0 is complete. Its findings that constrain your work:

- A 33-paper corpus of Mendelian randomization (MR) studies using a cardiac disease as EXPOSURE and
  cancer as OUTCOME was assembled; 29 have accessible full text; 192 exposure-outcome rows extracted.
- CRITICAL: instrument rsIDs are unrecoverable for 186 of 192 rows (96.9%) and the instrument count
  is NOT_REPORTED for 130 of 192. Published papers cite supplementary tables absent from PMC text.
  Therefore instruments must be RE-DERIVED from the named GWAS accessions, not read off the papers.
- The dominant heart-failure exposure accession in the corpus is GCST009541, verified as Shah 2020
  (PMID 31919418), 47,309 European cases / 930,014 controls.
- Prior work IN THIS PROJECT found heart failure's genome-wide-significant instruments are canonical
  SHARED-RISK-FACTOR loci: PITX2 (atrial fibrillation), 9p21 (coronary artery disease), LPA, SORT1,
  FTO (BMI). A positive control (HF -> NPPB protein) also failed at p=0.209. That is the hypothesis
  under test: published HF-exposure MR may be estimating shared risk factors, not heart failure.
- Published comparator for MR reporting quality (RECORD it, never recompute): doi 10.1093/ije/dyac074,
  87 MR-Base studies, mean 53% of reporting items satisfied, 44% sufficient for the relevance
  assumption, 31% for independence. This corpus scored 49.7%, 75.9% and 10.3% respectively.

PHASE 0 ARTIFACTS you may read with read_file(version_id=...) or host.artifact_path(vid):
  PRIOR_ART_BOUNDARY.md    541bf3fa-db40-4338-9e6a-24de912fb9f9   (read this FIRST; it has a
                                                                   do-not-attempt list you must obey)
  editorial_criticisms.csv 115007c7-85a5-4a7c-9f1c-cedba48dda5e
  prior_art_papers.csv     b97a6c2d-b19a-45c7-a7e5-c7e174beb61f
  mr_corpus.csv            587aff0f-f85d-45a1-8474-34acebcaec01
  mr_extraction.csv        644688bd-5ed5-4783-866b-f8e481c5e14c
  mr_reporting_gaps.csv    bc8b4de2-dc68-46e0-b562-abc2ab18cb27
  mr_reporting_summary.csv d931e9eb-6e8a-42c9-8f68-944ffdc35f71
  CORPUS_NOTES.md          e79c496a-e4d9-420f-81ea-40af5ae202d5

TOOLS. MCP calls work ONLY in the `repl` tool via host.mcp(server, method, **kwargs) — never in
`python`/`r`. Pass data between kernels via ./handoff/*.json. Servers: "pubmed" (search_articles,
get_article_metadata, convert_article_ids, get_full_text_article), "human-genetics"
(gwas_search_traits, gwas_associations_for_trait, gwas_associations_for_variant, gwas_get_study).
Full texts are NOT in your workspace — re-retrieve from PMC using the PMCIDs in mr_corpus.csv.

SKILL AUTHORING. You are publishing a real reusable skill. Load skill('customize') for the
host.skills.* API. Requirements, non-negotiable:
  - kernel.py at the skill root: pure numpy/pandas/scipy/statsmodels — NO host.* calls, NO network,
    NO MCP inside kernel.py, so it runs standalone. Functions + imports + literal constants only;
    no top-level classes or decorators (wrap those in a factory function).
  - Guardrails are ASSERTIONS that RAISE, not prose a reader can skip. This is the house convention
    from this project's crosstalk-score skill: require_* validators raise AssertionError by default
    and return an audit table so anything dropped is disclosed rather than silently removed.
  - SKILL.md must state what the skill does NOT establish, and must credit prior art rather than
    implying novelty where none exists.
  - Ship regression tests (test_*.py) that reproduce real measured numbers, and confirm they pass.
  - Publish with host.skills.publish(name) and report the skill name.

RULES.
- Never write a PMID, DOI, effect estimate, sample size or accession from memory — copy it from a
  tool result. If you cannot retrieve something, record the gap explicitly.
- Never fabricate data to fill a gap.
- Do not duplicate anything on the PRIOR_ART_BOUNDARY.md do-not-attempt list.
- Save deliverables with save_artifacts and report version_ids.

YOUR TASK — Track: Instrument audit. Build and publish the `mr-instrument-audit` skill. Three steps.

The purpose: given a cardiac exposure trait, RE-DERIVE the instrument set a published MR study would
have used, then test whether those instruments are specific to that trait or are really
shared-risk-factor loci. This is the machinery that decides whether the corpus's causal claims mean
what they say.

STEP 1 "Build instrument re-derivation kernel".
Implement in kernel.py:
  - clump_by_distance(assoc_df, window_kb=1000, p_threshold=5e-8) -> DataFrame: collapse
    genome-wide-significant associations into independent loci by genomic distance, keeping the
    lead variant per locus. Distance-based clumping is a documented APPROXIMATION to LD-based
    clumping — say so in the docstring and in SKILL.md; do not imply it is equivalent.
  - derive_instruments(assoc_df, ...) -> DataFrame: reconstruct the instrument set for a trait,
    handling the GWAS Catalog's multi-trait annotation so one locus is not counted twice, and
    handling rows with missing rsIDs or missing positions explicitly rather than dropping silently.
  - Data retrieval itself goes in the SKILL.md workflow using host.mcp("human-genetics",
    "gwas_associations_for_trait", efo_id=..., max_records=...), NOT inside kernel.py (kernel.py
    must stay dependency-free and offline-runnable).
Verify against a trait with known architecture before trusting it.

STEP 2 "Implement pleiotropy annotation and specificity scoring".
  - annotate_variant_traits(variant_trait_table) -> DataFrame: given each instrument variant's full
    trait-association profile (retrieved via host.mcp("human-genetics",
    "gwas_associations_for_variant", rs_id=...) in the workflow, passed in as data), classify the
    variant as CARDIAC_SPECIFIC or SHARED_RISK_FACTOR according to whether it is also associated
    with body mass index, smoking, blood lipids, blood pressure, type 2 diabetes, atrial
    fibrillation or coronary artery disease. Make the shared-risk-factor trait list an explicit,
    inspectable constant — never buried in a regex — and support EFO-id matching as well as label
    substring matching, because catalog labels vary ("blood VLDL cholesterol amount",
    "low density lipoprotein cholesterol measurement", "coronary artery disorder").
  - specificity_score(...) -> dict: the fraction of an instrument set that is shared-risk-factor,
    with the per-variant table so every classification is auditable.
  IMPORTANT — credit prior art. Variant-level phenome lookup ALREADY EXISTS as PhenoScanner V2
  (doi 10.1093/bioinformatics/btz469; over 150 million variants, 65 billion associations, EFO-mapped,
  with 1000 Genomes LD proxy search). Your contribution is programmatic corpus-scale application over
  the GWAS Catalog, NOT a new lookup method. SKILL.md must say this plainly.

STEP 3 "Add matched null, base-rate guard, tests and publish".
  - matched_null(...) -> dict: sample random variant sets MATCHED on the number of trait
    associations per variant (a highly pleiotropic variant is more likely to hit any trait list by
    chance, so an unmatched null is misleading — this is the matched-null lesson from this project's
    crosstalk-score skill, where an unmatched pool gave rho=0.323 and wrongly implied confounding).
    Return the null distribution and an empirical p-value for the observed shared-risk-factor fraction.
  - base_rate_guard(...) -> dict: the rate at which ANY variant set of comparable size and
    pleiotropy would be classified shared-risk-factor. If most variants qualify regardless, a high
    fraction says nothing. Make this RAISE or return an explicit uninformative verdict when the base
    rate is so high the test cannot discriminate.
  - Regression tests must reproduce the real finding that heart-failure instruments are dominated by
    shared-risk-factor loci. Derive HF instruments from the GWAS Catalog (heart failure is
    MONDO_0005252; note the catalog also carries congestive MONDO_0005009, diastolic MONDO_0006727,
    systolic MONDO_0006993) and confirm the known confounded loci appear: PITX2 near chr4:110.7Mb,
    9p21 near chr9:22.1Mb, LPA chr6:160.6Mb, SORT1 chr1:109.3Mb, FTO chr16:53.8Mb. Report the actual
    measured shared-risk-factor fraction with its matched-null p-value and base rate — whatever the
    numbers turn out to be. Do not tune anything to reach a desired answer; if the guard says the
    test is uninformative, that is the result and it must be reported as such.
Also save `hf_instrument_audit.csv` (per instrument variant: rsID, locus, mapped genes, n_traits,
trait labels, classification) and `instrument_audit_validation.csv` (the null and base-rate numbers).
Publish as `mr-instrument-audit`. SKILL.md must document the distance-clumping approximation, the
PhenoScanner prior art, what a shared-risk-factor fraction does NOT prove (it flags that an estimate
may not be exposure-specific; it does not establish the true causal effect), and that GWAS Catalog
curated top-hits are not full summary statistics.

Deliverables: published `mr-instrument-audit` skill (SKILL.md + kernel.py + tests),
hf_instrument_audit.csv, instrument_audit_validation.csv.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "summary": {
      "type": "string"
    },
    "skill_published": {
      "type": "string"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "filename": {
            "type": "string"
          },
          "version_id": {
            "type": "string"
          },
          "description": {
            "type": "string"
          }
        },
        "required": [
          "filename",
          "version_id"
        ]
      }
    },
    "key_numbers": {
      "type": "object"
    },
    "tests_passing": {
      "type": "string"
    },
    "gaps": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "summary",
    "skill_published",
    "artifacts",
    "key_numbers",
    "tests_passing",
    "gaps"
  ]
}
```


---

## Wave 4 — Instrument verdicts

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=132.0 min

### Task statement

```

CONTEXT.

This project audits the reverse cardio-oncology literature — the hypothesis that cardiac disease
promotes cancer. Phases 0 and 1 are complete. You run Phase 2: the audit itself.

WHAT PHASE 0 ESTABLISHED (all numbers copied from retrieved records; treat as given):
- 33-paper corpus of Mendelian randomization (MR) studies with a cardiac disease as EXPOSURE and
  cancer as OUTCOME. 29 have accessible full text. 192 exposure-outcome rows.
- Instrument rsIDs are UNRECOVERABLE for 186 of 192 rows (96.9%); n_instruments NOT_REPORTED for
  130 of 192. Instruments therefore MUST be re-derived from named GWAS accessions.
- Corpus exposure mix (cardiac_to_cancer rows): coronary/CAD/CHD/IHD/angina 55, atrial fibrillation
  43, myocardial infarction 20, heart failure 10, composite CVD 2. NOTE: heart failure is the LEAST
  studied cardiac exposure in a literature nominally about heart failure and cancer.
- Three data sources dominate (of 384 source fields): FinnGen 103, IEU OpenGWAS 97, UK Biobank 83.
  33 papers reusing three datasets are not 33 independent tests.
- 50 extraction rows across 8 papers draw BOTH exposure and outcome from UK Biobank/IEU; 7 rows
  across 2 papers draw both from FinnGen. Only 4 of 29 papers raise sample overlap at all.
- PUBLISHED PRIOR ART that constrains you (do not duplicate):
  * Shui 2025, Human Genomics, PMID 41366697, doi 10.1186/s40246-025-00888-6 — HF vs 11 cancers.
    LDSC rg: HF-skin 0.381 (P-FDR 1.12e-26), HF-endometrial 0.375 (P-FDR 5.11e-6), HF-lung 0.216
    (P-FDR 2.571e-4); NULL for breast (0.011), prostate (0.035), colorectal (-0.120). Yet MR was
    NULL in all three carried-forward pairs: HF->lung OR 1.05 (0.83-1.34), HF->skin OR 1.07
    (0.97-1.19), HF->endometrial OR 1.05 (0.76-1.45).
  * Mu 2026, Cardio-Oncology, PMID 41527136, doi 10.1186/s40959-025-00438-y — lung cancer x HF:
    ALL 48 pleiotropic SNPs at one locus, 15q25.1 (CHRNA3/CHRNA5/HYKK/PSMA4, the smoking locus);
    HF->lung IVW beta=-0.049 P=0.688. DO NOT re-derive the lung pair; it is published.
  * Yarmolinsky 2025 ATVB doi 10.1161/ATVBAHA.125.322433 — CAD x breast/colorectal/lung/prostate
    pleiotropy: 60 colocalized loci, 43 novel, 35 opposing-direction. Do not redo this.
  * Yang 2025 doi 10.3389/fgene.2025.1450259 — AF x breast cancer, full method stack, MR null both
    directions. Do not redo this.

WHAT PHASE 1 BUILT — two published skills you MUST use rather than reimplement:
  * `mr-instrument-audit` — load with skill('mr-instrument-audit'). Provides clump_by_distance,
    derive_instruments, annotate_variant_traits, specificity_score, matched_null, base_rate_guard,
    audit_instrument_set, flatten_locus_traits. Measured results already in hand:
      - HF accession GCST009541 (Shah, PMID 31919418, 47,309 cases / 930,014 controls) re-derives to
        10 independent loci, recovering PITX2 rs17042102, 9p21 rs1556516, LPA rs55730499,
        SORT1/CELSR2 rs660240, FTO rs56094641 (all within 48 kb of expected).
      - Lead-SNP level, that 10-locus set: 60.0% shared-risk-factor, matched-null mean 0.4348,
        empirical p=0.177, base rate 0.332 -> verdict NOT_DISTINGUISHABLE_FROM_MATCHED_NULL. n=10 is
        UNDERPOWERED. Do not overstate it.
      - Pooled cross-study top-40 HF locus set: 79.5% (31/39), p=2.0e-4 -> SHARED_RISK_FACTOR_DOMINATED.
      - Gene-region-level annotation: base rate 96.4% trips base_rate_guard -> UNINFORMATIVE_BASE_RATE.
        Lead-SNP annotation is the informative resolution; gene-region is not.
  * `claim-audit` — load with skill('claim-audit'). Corpus audit already produced: 55.6% of
    STROBE-MR items satisfied (published base rate 53.0%), 69.0% sufficient for relevance (vs 44.0%),
    3.4% sufficient for independence (vs 31.0%), and 24 of 29 papers with a causal-overclaim mismatch.

PHASE 0/1 ARTIFACTS — read with read_file(version_id=...) or host.artifact_path(vid):
  PRIOR_ART_BOUNDARY.md         541bf3fa-db40-4338-9e6a-24de912fb9f9  (READ FIRST — do-not-attempt list)
  mr_corpus.csv                 587aff0f-f85d-45a1-8474-34acebcaec01
  mr_extraction.csv             644688bd-5ed5-4783-866b-f8e481c5e14c  (192 rows; exposures/outcomes/estimates)
  mr_reporting_gaps.csv         bc8b4de2-dc68-46e0-b562-abc2ab18cb27
  mr_reporting_summary.csv      d931e9eb-6e8a-42c9-8f68-944ffdc35f71
  corpus_claim_audit.csv        1e154133-67ce-455f-984d-291cd8483608
  corpus_vs_published_base_rate.csv 847cd346-8754-471a-8a88-17ece73fa081
  hf_instrument_audit.csv       02fffede-69ce-47a5-b2fc-7bce9c12dc94
  instrument_audit_validation.csv 422d9b78-709c-4287-8545-e3e250e3dcfa
  editorial_criticisms.csv      115007c7-85a5-4a7c-9f1c-cedba48dda5e

TOOLS. MCP calls work ONLY in the `repl` tool via host.mcp(server, method, **kwargs) — never in
`python`/`r`. Pass data between kernels via ./handoff/*.json. Servers: "human-genetics"
(gwas_search_traits, gwas_search_studies, gwas_get_study, gwas_associations_for_trait,
gwas_associations_for_variant), "pubmed". The GWAS Catalog API is rate-limited and caps at ~1000
records per variant query — record when truncated=True rather than ignoring it.

RULES.
- Never write a PMID, DOI, accession, effect estimate or sample size from memory. Copy from a tool
  result. If you cannot retrieve it, record the gap.
- NEVER tune a threshold or trait list to reach a desired answer. If a test comes out underpowered
  or uninformative, that IS the result and you report it as such.
- Report every triggered kill criterion without post-hoc adjustment, including marginal triggers.
- Do not duplicate anything on the do-not-attempt list.

YOUR TASK — Phase 2, five steps in order. STEP 1 MUST BE COMPLETE AND SAVED BEFORE ANY VERDICT IS
COMPUTED.

STEP 1 "Fix kill criteria before the audit runs".
Write and SAVE `PRESPECIFICATION_AUDIT.md` before computing anything. It must fix, numerically and
in advance:
  - The primary statistic (shared-risk-factor fraction of a re-derived instrument set at lead-SNP
    resolution) and the exact shared-risk-factor trait family list used.
  - The matched null (pleiotropy-count-binned permutation), permutation count, and seed.
  - The base-rate guard threshold above which a result is declared UNINFORMATIVE.
  - The minimum instrument count below which a set is declared UNDERPOWERED and no significance
    claim is made (justify it from the n=10 p=0.177 result already measured).
  - Thresholds separating verdicts: INSTRUMENTS_ADEQUATE / PARTIALLY_CONFOUNDED /
    SHARED_RISK_FACTOR_DOMINATED / UNDERPOWERED / UNINFORMATIVE_BASE_RATE.
  - What would FALSIFY the project's hypothesis: state explicitly what result would show published
    cardiac-exposure MR instruments ARE exposure-specific. If the audit returns that, you report it.
  - The positive and negative control pass conditions from step 4, fixed in advance.

STEP 2 "Audit instrument validity across the corpus".
The corpus is dominated by CAD (55 rows) and AF (43 rows), NOT heart failure (10 rows) — and the
HF single-study set is already known underpowered. So audit the exposures the corpus actually uses,
which also have more loci and therefore more power:
  - Resolve each distinct cardiac exposure in mr_extraction.csv to its GWAS Catalog trait id
    (gwas_search_traits). Cover at minimum: coronary artery disease, atrial fibrillation,
    myocardial infarction, heart failure, and any cardiomyopathy exposure present.
  - For each, derive instruments and run the full audit (specificity_score + matched_null +
    base_rate_guard) at LEAD-SNP resolution using the `mr-instrument-audit` skill.
  - Where a paper named a specific accession (e.g. GCST009541 for HF), audit that accession's
    instrument set specifically, and separately audit the pooled trait-level set. Report both and
    label which is which — they answer different questions.
  - Save `corpus_instrument_verdicts.csv`: one row per (exposure, source) with exposure_trait,
    efo_id, accession_or_pooled, n_loci, n_classified, srf_fraction, matched_null_mean, empirical_p,
    base_rate, verdict, srf_family_counts, n_corpus_rows_relying_on_this_exposure.
  - Cross-link to the corpus: for each paper/row, attach the verdict for its exposure, and record
    from mr_reporting_gaps.csv whether that paper's OWN reported sensitivity analyses (Egger, Q,
    leave-one-out, Steiger) could in principle have detected instrument non-specificity. Note the
    key methodological point: MR-Egger and Cochran Q detect DIRECTIONAL pleiotropy and heterogeneity,
    but instruments that are uniformly shared-risk-factor loci can produce a HOMOGENEOUS, non-
    directionally-pleiotropic, entirely confounded estimate that passes both tests. State whether
    that is the case here. Save `corpus_row_verdicts.csv`.

STEP 3 "Re-examine the strongest published positive claim".
The corpus's strongest positive HF->cancer claim is doi 10.3389/fimmu.2024.1454021 (PMID 39346905),
which reported HF associated with ~30% increased colorectal cancer risk in discovery (OR 1.32,
95% CI 1.03-1.69, P=0.025), non-significant in validation (OR 1.19, 0.97-1.46, P=0.090), and a
pooled OR 1.24 (95% CI 1.06-1.25, P=0.007). VERIFY those numbers against the retrieved full text
(PMC11427256) before using them — and note the printed pooled CI 1.06-1.25 does not contain its own
point estimate 1.24 correctly at the upper bound, so check whether the paper misprinted it.
Then: audit the instruments for its HF exposure, determine what fraction are shared-risk-factor,
and assess whether the claim survives. Prior project work found an analogous ACE->colorectal signal
was carried entirely by a single variant (rs4363) and collapsed on leave-one-out, so single-variant
dependence is specifically in scope. Report whether the paper's own sensitivity analyses could have
caught what you find. Save `strongest_claim_reexamination.csv` and a short prose verdict.

STEP 4 "Run positive and negative controls".
An auditor that cannot pass a valid instrument set cannot be trusted to fail an invalid one.
  - POSITIVE CONTROL (must pass as exposure-specific-ish): instruments for a trait whose genetics
    are well characterised and specific — e.g. LDL cholesterol instruments audited against a trait
    list that EXCLUDES blood lipids, or a similarly constructed self-consistent control. Design it
    so the correct answer is known in advance and state that design in the pre-specification.
  - NEGATIVE CONTROL (must fail): a set known to be dominated by a shared risk factor — the
    15q25.1 smoking locus for lung cancer is the published exemplar (Mu 2026).
  - Also run a NULL-TRAIT control: an exposure with no plausible relation to the shared-risk-factor
    list, to check the classifier does not label everything shared.
  Save `audit_controls.csv` with the pre-specified pass condition and whether each control met it.
  If a control FAILS, say so plainly and treat the corresponding audit results as unreliable.

STEP 5 "Scope the residual open question".
Do NOT run an underpowered cardiomyopathy-vs-cancer overlap test. Instead quantify the ceiling:
  - For each of the nine cardiomyopathy EFO traits in the catalog (cardiomyopathy MONDO_0004994,
    dilated MONDO_0005021, hypertrophic MONDO_0005045, idiopathic dilated EFO_0009094, ischemic
    EFO_0001425, nonischemic EFO_0009881, peripartum MONDO_0018920, Tako-tsubo MONDO_0019018,
    Chagas MONDO_0005491): count catalog associations, count independent loci after clumping, and
    pull the largest study's case/control counts via gwas_search_studies / gwas_get_study.
  - Do the same for the cancer traits, and quantify the ancestry and sample-size asymmetry between
    the cardiac and cancer sides.
  - State how many independent instruments a cardiomyopathy exposure would need for the audit to be
    informative given the n=10 -> p=0.177 result, and whether any cardiomyopathy trait reaches it.
  - ALSO address the genuinely open question Phase 0 identified: Shui 2025 reports rg=0.381
    (HF-skin) and rg=0.375 (HF-endometrial) yet null MR. Explain the mechanism of that discordance
    in terms of what you measured: genetic correlation is driven by genome-wide shared polygenic
    signal including shared risk factors, while MR uses only genome-wide-significant instruments.
    Say clearly whether your instrument audit supports reframing those published nulls as
    UNINFORMATIVE (instruments not exposure-specific, so the null does not test the HF hypothesis)
    rather than as evidence of no causal effect. Be careful and honest: if your data do not support
    that reframing, say so.
  Save `residual_question_scoping.csv` and `RESIDUAL_SCOPING.md`.

Deliverables: PRESPECIFICATION_AUDIT.md, corpus_instrument_verdicts.csv, corpus_row_verdicts.csv,
strongest_claim_reexamination.csv, audit_controls.csv, residual_question_scoping.csv,
RESIDUAL_SCOPING.md.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "summary": {
      "type": "string"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "filename": {
            "type": "string"
          },
          "version_id": {
            "type": "string"
          },
          "description": {
            "type": "string"
          }
        },
        "required": [
          "filename",
          "version_id"
        ]
      }
    },
    "key_numbers": {
      "type": "object"
    },
    "verdicts": {
      "type": "object",
      "description": "per-exposure verdict summary"
    },
    "controls_passed": {
      "type": "string"
    },
    "kill_criteria_triggered": {
      "type": "array",
      "items": {
        "type": "string"
      }
    },
    "gaps": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "summary",
    "artifacts",
    "key_numbers",
    "verdicts",
    "controls_passed",
    "kill_criteria_triggered",
    "gaps"
  ]
}
```


---

## Wave 4 — Synthesis

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=51.4 min

### Task statement

```

CONTEXT. You write the FINAL SYNTHESIS (Phase 3) of a project auditing the reverse cardio-oncology
literature — the hypothesis that cardiac disease promotes cancer. Phases 0-2 are complete. Your job
is adjudication, number-auditing and the final report. You are NOT running new analyses.

THE PROJECT IN ONE LINE: published Mendelian randomization (MR) studies use cardiac diseases as
exposures for cancer; this project tested whether their genetic instruments are specific to those
cardiac diseases or are really shared-risk-factor loci (BMI, smoking, lipids, blood pressure, type 2
diabetes, atrial fibrillation, coronary artery disease).

=== ALL PRIOR ARTIFACTS. Read with read_file(version_id=...) or host.artifact_path(vid). ===
Phase 0 (evidence base):
  PRIOR_ART_BOUNDARY.md              541bf3fa-db40-4338-9e6a-24de912fb9f9
  prior_art_papers.csv               b97a6c2d-b19a-45c7-a7e5-c7e174beb61f
  editorial_criticisms.csv           115007c7-85a5-4a7c-9f1c-cedba48dda5e
  mr_corpus.csv                      587aff0f-f85d-45a1-8474-34acebcaec01
  mr_corpus_query_log.csv            b584fdc8-8e0a-4da2-8066-55791e742e6d
  mr_extraction.csv                  644688bd-5ed5-4783-866b-f8e481c5e14c
  mr_reporting_gaps.csv              bc8b4de2-dc68-46e0-b562-abc2ab18cb27
  mr_reporting_summary.csv           d931e9eb-6e8a-42c9-8f68-944ffdc35f71
  CORPUS_NOTES.md                    e79c496a-e4d9-420f-81ea-40af5ae202d5
Phase 1 (skills + benchmarks):
  appraiser_benchmark.csv            e3280ac3-fa3d-4884-b9c2-53be0c4f1388
  appraiser_benchmark_summary.csv    e17453f6-39c0-4ac6-90e2-34a8389443c0
  overturned_claim_fixtures.csv      452abf42-a43f-4241-a00e-9670099dae32
  corpus_claim_audit.csv             1e154133-67ce-455f-984d-291cd8483608
  corpus_vs_published_base_rate.csv  847cd346-8754-471a-8a88-17ece73fa081
  appraiser_phase0_concordance.csv   a3210548-7283-4c2c-ae7d-a885664d96c4
  hf_instrument_audit.csv            02fffede-69ce-47a5-b2fc-7bce9c12dc94
  instrument_audit_validation.csv    422d9b78-709c-4287-8545-e3e250e3dcfa
Phase 2 (the audit):
  PRESPECIFICATION_AUDIT.md          51d2b2da-2e2e-406a-9b97-462f28e9e974
  corpus_instrument_verdicts.csv     90a16f36-d295-4113-9e46-a82123152741
  corpus_row_verdicts.csv            92445070-55e9-409a-ac28-7d275adf5a0e
  strongest_claim_reexamination.csv  47d15364-3373-4562-9a5a-cc578fefd369
  STRONGEST_CLAIM_VERDICT.md         1be4589f-e808-4204-a1b4-de34cdd8eae1
  strongest_claim_leave_one_out.csv  9b6f7c61-0733-471a-bb95-546f7dfab7b4
  audit_controls.csv                 d4dcdb04-18c5-483b-8ef0-cebe354bd6a4
  residual_question_scoping.csv      82427ff3-ab5e-40c2-9e8f-a54ff8ebd0a7
  RESIDUAL_SCOPING.md                37b260fe-aa6d-4a2b-812f-158df9a2e817
  PROTOCOL_DEVIATIONS.md             b2dd60c1-34df-4d95-9030-16a01b9dfba8
  instrument_audit_summary.png       cd6d999f-160a-4e81-95d0-390e79d68c31

Two published skills exist and are what this project shipped: `claim-audit` and
`mr-instrument-audit`. Load them if you need their documented semantics
(skill('claim-audit'), skill('mr-instrument-audit')).

=== HEADLINE RESULTS TO ADJUDICATE (verify every one against the tables; do not trust this summary) ===
Instrument verdicts at lead-SNP resolution, pooled top-40 sets unless stated:
  heart failure 0.875 (p=2e-4) SHARED_RISK_FACTOR_DOMINATED; coronary atherosclerosis 0.750 (p=0.0082)
  DOMINATED; coronary artery disease 0.675 (p=0.04099) DOMINATED (marginal band); myocardial infarction
  0.600 (p=0.0048) DOMINATED (exactly on the 0.60 boundary); angina 0.575 (p=0.0064)
  PARTIALLY_CONFOUNDED; MI accession GCST011365 0.5256 (p=2e-4, n=78) PARTIALLY_CONFOUNDED.
  UNDERPOWERED (n<22 floor, no significance claim): HF GCST009541 n=10 0.600 p=0.158; AF GCST006414
  n=18; AF GCST006061 n=17.
  ATRIAL FIBRILLATION pooled: 0.250, p=0.733, n=40 -> INSTRUMENTS_ADEQUATE. This was the
  PRE-SPECIFIED FALSIFICATION RESULT for the corpus's most-used exposure (43 rows). It was then
  judged UNRELIABLE_annotation_depth_confound.
  88 of 150 cardiac-to-cancer rows rest on a flagged exposure; 58 rows across 11 papers report clean
  MR-Egger AND Cochran Q while relying on flagged instruments; 80 of 88 flagged rows never address
  sample overlap.
Reporting quality (claim-audit over 29 papers) vs published base rates from doi 10.1093/ije/dyac074
(87 MR-Base studies; RECORDED, never recomputed): 55.6% of STROBE-MR items satisfied vs 53.0%;
69.0% sufficient for relevance vs 44.0%; 3.4% sufficient for independence vs 31.0%. 24 of 29 papers
show a causal-overclaim mismatch.
Appraiser benchmark: Avraham 2020 full text is PAYWALLED and recall against it is UNMEASURABLE.
On the abstract: 6/6 abstract-detectable (1.000), 6/13 all-mapped (0.462), precision 0.857.
Koelwyn full-text arm: 5/13 all-mapped (0.385), 5/9 material-present (0.556), precision 0.833.
Independent concordance with Phase 0 hand scoring: 233/261 = 89.3%.

=== THREE CONSTRAINTS YOU MUST HONOUR. Do not soften these. ===
1. TWO OF THREE CONTROLS FAILED (audit_controls.csv). Positive control P-B failed at the boundary
   (0.600 not < 0.60). The NEGATIVE control failed on fraction (0.433 vs required 0.80) though its
   modal-family clause passed — meaning the auditor UNDER-calls shared status where per-variant
   annotation is sparse. Only the null-trait control (vitiligo, 0.1765, p=0.755, n=68) passed
   outright, and it is the control that protects the positive findings from the "labels everything
   shared" failure mode. A report that presents the four DOMINATED verdicts without these failures
   prominently attached is not acceptable.
2. ANNOTATION-DEPTH CONFOUND. Classification tracks how deeply a variant is annotated, not only its
   true pleiotropy: r=0.240, p=3.2e-09 across 593 variants; shared-risk-factor rate 0.406 in the
   lowest annotation quartile rising to 0.923 in the highest; AF instruments carry a median of 2.0
   catalog trait annotations versus 13.5 for CAD and 16.0 for heart failure. This cuts BOTH ways —
   it is why the AF falsification was withdrawn, AND it is a bias toward the positive findings for
   well-studied exposures. PROTOCOL_DEVIATIONS.md records that the top-40 pooled sets are biased
   toward most-significant hence most-annotated loci, a bias running TOWARD higher fractions. State
   this as a limitation on the DOMINATED verdicts, not only on AF.
3. THE SHUI REFRAMING IS NOT SUPPORTED. Shui 2025 (PMID 41366697) used exposure accession
   GCST90399713 (retrieved: PMID 36777996, 52,496 cases / 967,945 controls), which had ZERO
   association rows in the retrieved set — so the instruments behind ITS published nulls were never
   audited. The general mechanism (genetic correlation uses genome-wide polygenic signal, MR uses
   only genome-wide-significant instruments) may be explained, but you MUST NOT claim this project
   showed Shui's specific nulls are uninformative. Report it as a gap.

=== YOUR TASK — three steps ===

STEP 1 "Adjudicate against pre-specified criteria".
Read PRESPECIFICATION_AUDIT.md and compare EVERY result to the thresholds fixed in advance. Produce
`ADJUDICATION.csv`: criterion_id, criterion_as_prespecified, observed_value, verdict
(MET / NOT_MET / MARGINAL / UNINFORMATIVE / NOT_APPLICABLE), and consequence for the project's claims.
Include every triggered criterion, every marginal trigger (CAD p=0.04099 in the 0.04-0.06 band; MI
fraction exactly 0.6000), both control failures, and the returned-then-withdrawn AF falsification.
Report marginal triggers as triggered, with NO post-hoc adjustment.

STEP 2 "Audit every number against source tables".
This is a verification step and it has caught real errors in this project before — a prior audit found
three numerical errors, all from quoting a summary instead of the column it described. So: recompute
each numeric claim you intend to publish DIRECTLY from the saved result tables, never from any prose
summary (including the summary in this task message). Check in particular:
  - every instrument verdict fraction, p-value, n_loci and corpus-row count against
    corpus_instrument_verdicts.csv and corpus_row_verdicts.csv;
  - the reporting percentages against corpus_claim_audit.csv / corpus_vs_published_base_rate.csv, and
    confirm the 53/44/31 comparator is recorded from doi 10.1093/ije/dyac074 and NOT recomputed;
  - the benchmark recall/precision figures and their DENOMINATORS against appraiser_benchmark*.csv;
  - the 88/150, 58-rows-across-11-papers and 80/88 cross-tabulations by recomputing them from
    corpus_row_verdicts.csv;
  - the strongest-claim numbers against strongest_claim_reexamination.csv and
    strongest_claim_leave_one_out.csv, including the misprint finding (printed pooled CI upper 1.25
    vs recomputed 1.453, point estimate 1.241, lower 1.060, P=0.0072).
Save `NUMBER_AUDIT.csv`: claim_location, claim_as_stated, recomputed_value, source_table, verdict
(VERIFIED / ERROR / INCONSISTENT / STALE / AMBIGUOUS), and note. If you find an error, CORRECT it in
the report and record it — do not quietly drop it.

STEP 3 "Write verdict report and figures".
Write `VERDICT_REPORT.md` in the register of a methods-and-results section: connected prose
paragraphs, each opening on a synthetic claim and then spending evidence, with inline DOI links as
[Author Year](https://doi.org/...). No bulleted bibliography. It must state plainly:
  (a) WHICH published causal claims survive instrument scrutiny and which do not, naming papers and
      DOIs, and being explicit that a flagged instrument set means an estimate may not be
      exposure-specific — it does NOT establish that the true causal effect is zero, nor that any
      paper's arithmetic is wrong.
  (b) WHAT IS COMPUTATIONALLY UNPROVABLE from public data and the specific measurement each would
      require. Draw on the computability triage and this project's earlier catalogued gaps
      (circulating pleiotrophin never measured in a heart-failure cohort; human plasma periostin in
      heart failure unmeasured; no human study demonstrating a heart-failure-derived ligand engaging
      a tumour receptor; tissue mRNA not convertible to plasma concentration without translation
      efficiency, secretion rate, clearance and volume of distribution).
  (c) The reporting-quality finding, with the independence assumption (3.4% vs 31.0%) as the centre
      of gravity, and the mechanism: MR-Egger detects DIRECTIONAL pleiotropy and Cochran Q detects
      HETEROGENEITY, so instruments that are uniformly shared-risk-factor loci produce a homogeneous,
      non-directionally-pleiotropic, entirely confounded estimate that passes both tests. Support
      this with the measured diffuseness (each dominated set spreading across 6-7 families with modal
      share 0.19-0.30) and the 58-rows/11-papers count.
  (d) A LIMITATIONS section carrying the three constraints above at full strength, plus: distance
      clumping is not LD clumping; GWAS Catalog curated top-hits are not full summary statistics;
      API truncation (CAD 600 of 4,081; AF 600 of 3,124) makes all locus counts lower bounds;
      instrument rsIDs were unrecoverable for 186 of 192 corpus rows so every set is a re-derivation
      and never a paper's actual list; 4 of 33 papers had no accessible full text.
  (e) What a stronger version of this audit would require (full summary statistics, an LD reference
      panel, PhenoScanner-style LD-proxy annotation, institutional access to the paywalled papers).
  (f) Credit prior art properly: PhenoScanner V2 (doi 10.1093/bioinformatics/btz469) for
      variant-phenome lookup; STROBE-MR (doi 10.1001/jama.2021.18236) as the checklist encoded;
      Mu 2026 (doi 10.1186/s40959-025-00438-y) as the published precedent showing HF-lung pleiotropy
      reduces to the 15q25.1 smoking locus. This project's contribution is corpus-scale application
      and the reframing of instrument non-specificity as a reporting-and-validity problem.

FIGURES — build with matplotlib. Load skill('figure-style') FIRST and call apply_figure_style()
before rendering; these are deliverables. Produce at least:
  fig1_corpus_growth.png — the corpus growth curve by year (from mr_corpus.csv; 28 of 33 papers are
    2024 or later), ideally alongside the seed-query year counts.
  fig2_instrument_specificity.png — per-exposure shared-risk-factor fraction with the matched-null
    mean and the empirical p, ordered, with verdict encoded and the underpowered/unreliable sets
    visually distinguished (do NOT present AF as a clean adequate result).
  fig3_reporting_vs_baserate.png — corpus reporting percentages against the published 53/44/31 base
    rates, with the independence collapse as the visual focus.
Every figure must use fig, ax = plt.subplots() and fig.savefig(...) — never plt.savefig. Caption each
in the report. Do not invent data for a figure; if a number is not in the tables, leave it out.

Deliverables: ADJUDICATION.csv, NUMBER_AUDIT.csv, VERDICT_REPORT.md, fig1_corpus_growth.png,
fig2_instrument_specificity.png, fig3_reporting_vs_baserate.png.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "summary": {
      "type": "string"
    },
    "artifacts": {
      "type": "array",
      "items": {
        "type": "object",
        "properties": {
          "filename": {
            "type": "string"
          },
          "version_id": {
            "type": "string"
          },
          "description": {
            "type": "string"
          }
        },
        "required": [
          "filename",
          "version_id"
        ]
      }
    },
    "number_audit": {
      "type": "object",
      "description": "counts of VERIFIED/ERROR/etc and any errors found"
    },
    "adjudication": {
      "type": "object"
    },
    "key_numbers": {
      "type": "object"
    },
    "gaps": {
      "type": "array",
      "items": {
        "type": "string"
      }
    }
  },
  "required": [
    "summary",
    "artifacts",
    "number_audit",
    "adjudication",
    "key_numbers",
    "gaps"
  ]
}
```


---

## Wave 5 — limma

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=21.1 min

### Task statement

```
Validate a custom differential-expression implementation against canonical limma. This is the
decisive item in a reviewer critique: the manuscript implemented an "empirical-Bayes moderated t-test
equivalent to limma-trend" in NumPy/SciPy instead of calling limma, and a reviewer asked whether the
5,036-gene signature partly reflects implementation choices.

FIX THE CRITERION BEFORE YOU RUN ANYTHING. Write it into limma_validation.csv as its own row:
  VALIDATED requires BOTH (a) Pearson r >= 0.99 between custom and limma log2 fold changes on shared
  genes, AND (b) >= 95% Jaccard-style overlap of the FDR<0.05 gene sets.
  Anything less is NOT validated and must be reported as a caveat on the signature. Do not soften this
  after seeing the numbers, and do not pick the better of several comparison variants post hoc.

DATA. Re-fetch both cohorts from GEO (allowlisted; use the geo-bulk-de skill, load it first):
  - GSE116250: left-ventricular RNA-seq, submitter-supplied RPKM. 37 dilated + 13 ischaemic
    cardiomyopathy vs 14 non-failing donors; the primary contrast pools the 50 failing against the 14
    non-failing. The original analysis used a limma-trend-equivalent moderated t on log2(RPKM+1).
    65 genes carrying a sentinel value of 999999999 were detected and removed before analysis - you
    must reproduce that exclusion or your gene set will not match.
  - GSE141910: MAGNet left ventricle, counts. 166 dilated + 28 hypertrophic cardiomyopathy vs 166
    non-failing donors, with 6 peripartum cardiomyopathy samples EXCLUDED, giving 194 vs 166.
  Group labels must be parsed from the series-matrix characteristics fields, not hardcoded.

CANONICAL COMPARISON. In the R environment named `limma-check` (limma 3.66.0, edgeR 4.8.2 already
installed - use environment="limma-check" on r tool calls):
  - GSE116250: lmFit on log2(RPKH+1) with a two-group design, then eBayes(fit, trend=TRUE). This is
    limma-trend.
  - GSE141910: counts, so run BOTH (i) edgeR TMM normalisation + voom + lmFit + eBayes, and (ii)
    log2-CPM + eBayes(trend=TRUE). Report both. Say which one the original implementation is closer to
    and why; do not just report the flattering one.
  Extract logFC, t, P.Value and adj.P.Val (BH) per gene via topTable(..., number=Inf, sort.by="none").

WHAT TO COMPARE, per cohort:
  - Pearson AND Spearman on log2FC over shared genes; same for the moderated t statistic.
  - FDR<0.05 gene-set overlap: n in both, n custom-only, n limma-only, Jaccard.
  - The empirical-Bayes hyperparameters side by side: limma's prior degrees of freedom (fit$df.prior)
    and prior variance (fit$s2.prior) against whatever the custom implementation estimated. If the
    custom code's values are not recoverable from its outputs, say so rather than guessing.
  - Whether the 5,036-gene replicated signature is recovered when the SAME replication rule (FDR<0.05
    and concordant direction in both cohorts) is applied to the limma results. Report the recovered
    count and the overlap with the published 5,036.
  - The canonical markers: NPPA, NPPB, MYH6, ATP2A2 - logFC and FDR under both implementations.

CUSTOM-SIDE VALUES to compare against are in these artifacts (read them, do not recompute):
  hf_de_GSE116250.csv  /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/a06e693e-0e8a-47c1-a618-3f29ee78d4db/v6a2d301f_hf_de_GSE116250.csv
  hf_de_GSE141910.csv  /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/f76fe119-05ac-415b-a1f1-d3b80bf512c8/vcac6c07c_hf_de_GSE141910.csv
  hf_signature_replicated.csv  /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/a379a891-5d69-4bdd-a668-9b71abbe27c1/v76685c52_hf_signature_replicated.csv
Inspect their columns first; do not assume names.

BINDING CONSTRAINTS:
- If GEO is unreachable or a cohort cannot be rebuilt, report that as BLOCKED and say which comparison
  could not be made. NEVER simulate, synthesise or approximate an expression matrix, and never
  fabricate a limma result. A blocked comparison reported honestly is the correct outcome; an invented
  one is not.
- Report disagreement at full strength if you find it. A finding that the custom implementation
  diverges is a legitimate and important result - the manuscript would be corrected, not defended.
- Any figure you build must be a final deliverable: load the figure-style skill and call
  apply_figure_style() first, save with fig.savefig (never plt.savefig), and verify the rendered
  output before saving.
- Every number you report must come from a computation in this session.

OUTPUTS (save as artifacts):
  1. limma_validation.csv - one row per comparison metric per cohort per limma variant, plus the
     pre-fixed criterion row and a final VERDICT row reading VALIDATED or NOT VALIDATED.
  2. limma_de_GSE116250.csv and limma_de_GSE141910.csv - the canonical limma per-gene results.
  3. fig9_limma_concordance.png - log2FC scatter (custom vs limma) per cohort with the correlation
     annotated, and a panel showing the FDR<0.05 set overlap. Only build it if the comparison ran.
  4. A short markdown note limma_validation_note.md stating the verdict in prose, the exact
     hyperparameters from both sides, and any caveat a reviewer should be told.

Report in structured output: verdict, per-cohort logFC Pearson and Spearman, gene-set Jaccard,
recovered signature size vs 5036, df.prior and s2.prior from limma, which GSE141910 pipeline you
treated as primary and why, canonical marker agreement, and anything BLOCKED.
```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "verdict": {
      "type": "string"
    },
    "gse116250": {
      "type": "object"
    },
    "gse141910": {
      "type": "object"
    },
    "signature_recovery": {
      "type": "object"
    },
    "hyperparameters": {
      "type": "object"
    },
    "gse141910_primary_pipeline": {
      "type": "string"
    },
    "canonical_markers": {
      "type": "string"
    },
    "criterion_as_fixed": {
      "type": "string"
    },
    "blocked": {
      "type": "string"
    },
    "artifact_ids": {
      "type": "object"
    }
  },
  "required": [
    "verdict",
    "gse116250",
    "gse141910",
    "signature_recovery",
    "criterion_as_fixed",
    "artifact_ids"
  ]
}
```


---

## Wave 5 — Sensitivity

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=49.7 min

### Task statement

```
Answer six methodological challenges raised by a reviewer of a ligand-receptor compatibility
screen. Each is a "why is this choice right" question and each needs a computed answer, not an argument.
Load the crosstalk-score skill first - it implements the scoring, the expression floor and the
matched-random null, and its guardrails raise on violation.

INPUT ARTIFACTS (inspect columns before use; do not assume names):
  crosstalk_network.csv                /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/890af4bd-6338-410b-8744-e13de45260eb/v8086ab6e_crosstalk_network.csv
  pancancer_susceptibility_ranking.csv /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/c32f62f4-275f-4aac-8b97-a84fe1dc6031/vc8e165e9_pancancer_susceptibility_ranking.csv
  tcga_expr_by_tumor_type.csv          /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/51bd31c5-0e36-4709-ac73-f5ffb358e9ce/v07e2a94c_tcga_expr_by_tumor_type.csv
  cellchatdb_receptor_subunits.csv     /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/9e2362f9-8f0a-4331-8d10-ff0e6074cd6c/v4aabe49f_cellchatdb_receptor_subunits.csv
  axis_survival.csv                    /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/bbbb5618-f77e-4ce2-a4fe-4046c6b48263/v2d0e9667_axis_survival.csv
  floor_audit.csv                      /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/2536d056-b29a-4b53-9fca-5dadb9ff351e/v89e852c0_floor_audit.csv
  twofactor_decomposition.csv          /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/35a3ae86-3e5b-4d75-9b98-ef09470e3331/v23bff2cf_twofactor_decomposition.csv
  hf_de_GSE116250.csv                  /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/a06e693e-0e8a-47c1-a618-3f29ee78d4db/v6a2d301f_hf_de_GSE116250.csv
  hf_de_GSE141910.csv                  /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/f76fe119-05ac-415b-a1f1-d3b80bf512c8/vcac6c07c_hf_de_GSE141910.csv
Published values you are testing against: floor = 0.005 x GAPDH = 345.7 in RSEM normalized-count units
(GAPDH median across the 30 per-type medians is 69,134); matched-random null mean rho 0.363 with no
floor, 0.107 with the floor, 0.323 for an unmatched naive pool; 4 of 30 tumour types beat the
specificity null; seed 20260727; 300 permutations; 188 interactions from 51 ligands; 103 receptor
subunits; 146 ligand-receptor-subunit axes; 7 surviving PTN axes; 74% receiver base rate.

TASK 1 - Is 0.005 x GAPDH tuned or on a plateau? (-> floor_sensitivity.csv)
Recompute the ranking and matched-random null across floor fractions 0.001, 0.0025, 0.005, 0.01, 0.025,
0.05 of GAPDH, plus two alternatives: an absolute-percentile floor (receptor above the Nth percentile of
expressed genes) and no floor. Per setting report null mean rho, fraction of draws above 0.5,
interactions excluded, and how many of 30 types beat the null. State plainly whether 0.005 sits on a
plateau or whether headline conclusions move with the choice. If they move, that is the finding.

TASK 2 - Is "heart ranks first among 52 GTEx tissues" the right specificity definition?
(-> specificity_definition_sensitivity.csv)
Re-score ligand cardiac specificity under: rank 1 (published), top 3, top 5, upper half, and fold
enrichment (heart median vs all-tissue median >= 2x and >= 5x). Report how many of the 42 HF-upregulated
secreted ligands pass each, and whether NPPA and NPPB remain the only ligands both cardiac-enriched and
abundant under every definition. Say which definitions are more permissive and whether any admits a
ligand the paper excluded. GTEx v8 median TPM by tissue is fetchable; if it is not reachable, use the
ligand_tissue_specificity table if you can locate it via host.artifacts(search=...) and say so.

TASK 3 - Does 300 permutations resolve the reported p-values? (-> permutation_resolution.csv)
Reported empirical p-values go as low as 0.003 = 1/300, at the resolution limit. Re-run the principal
specificity null at 300, 1000, 5000 and 20000 draws and report per tumour type the p-value at each depth,
the Monte Carlo standard error at each depth, and whether the set of types passing p<0.05 changes.
Conclude whether 300 sufficed; if the pass set changes, give the corrected set.

TASK 4 - Does the matched null pool need to be wider? (-> null_pool_width.csv)
Vary the expression-decile-matched pool width and report null mean rho and pass count for several
widths, and whether the conclusion is stable.

TASK 5 - Does expanding receptor complexes into subunits double-count? (-> complex_independence.csv)
The screen expanded heteromeric complexes into 103 subunits. Recompute axis-survival tallies at COMPLEX
level (complex counts once, using the published min-across-subunits rule, and separately the geometric
mean) rather than per subunit. Report how many complexes contribute more than one subunit, the effective
number of independent tests, and whether any headline count changes - specifically the 7 surviving PTN
axes and the 74% receiver base rate. If subunit expansion inflated a count, say which.

TASK 6 - Does pooling HF etiologies mask an etiology-specific signal? (-> etiology_stratified.csv)
GSE116250 pools 37 dilated + 13 ischaemic vs 14 non-failing; GSE141910 pools 166 dilated + 28
hypertrophic vs 166 non-failing. Re-derive differential expression within each etiology separately
against the same controls (re-fetch from GEO with the geo-bulk-de skill; group labels parsed from
series-matrix characteristics, and GSE116250 has 65 genes carrying a sentinel 999999999 that must be
removed). Then ask of the SECRETED ligands: is any strongly induced in one etiology and not the other,
and would it have passed the cardiac-specificity and abundance criteria NPPA and NPPB passed? Report
per-etiology log2FC and FDR for the 42 HF-upregulated secreted ligands. Note the power asymmetry (13
ischaemic, 28 hypertrophic are small) and do not read a null in a small stratum as absence.

OUTPUTS: the six CSVs, plus sensitivity_summary.json with a one-line verdict per task, plus
fig10_sensitivity.png IF AND ONLY IF at least one sweep shows a headline conclusion is choice-dependent.
If everything is stable, say so and skip the figure rather than plotting null findings.

Report in structured output: per-task verdict (STABLE / CHOICE-DEPENDENT / BLOCKED), whether 0.005 is on
a plateau, whether 300 draws suffice and the corrected pass set if not, whether subunit expansion
inflated any count and which, whether any etiology-specific ligand would have changed the conclusion,
and all artifact ids.

BINDING CONSTRAINTS (project convention, non-negotiable):
- Fix any threshold or criterion BEFORE you compute the thing it judges, and write it into the output
  table as its own row. Do not adjust it after seeing results.
- NEVER simulate, synthesise or approximate data. If a resource is unreachable or a quantity is not
  recoverable, report BLOCKED or UNTESTABLE and say exactly which comparison could not be made.
  A blocked result reported honestly is correct; an invented one is not.
- Report findings that undermine the manuscript at full strength. The manuscript gets corrected, not
  defended. Several results in it are already withdrawals and negative findings.
- Every number must come from a computation or a fetch in this session. Nothing from memory.
- Any figure must be a final deliverable: load the figure-style skill, call apply_figure_style() first,
  use fig.savefig (never plt.savefig), and inspect the render before saving.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "task_verdicts": {
      "type": "object"
    },
    "floor_plateau": {
      "type": "string"
    },
    "permutation_sufficient": {
      "type": "string"
    },
    "corrected_pass_set": {
      "type": "string"
    },
    "subunit_inflation": {
      "type": "string"
    },
    "etiology_finding": {
      "type": "string"
    },
    "specificity_definition_effect": {
      "type": "string"
    },
    "null_pool_effect": {
      "type": "string"
    },
    "blocked": {
      "type": "string"
    },
    "artifact_ids": {
      "type": "object"
    }
  },
  "required": [
    "task_verdicts",
    "floor_plateau",
    "permutation_sufficient",
    "subunit_inflation",
    "etiology_finding",
    "artifact_ids"
  ]
}
```


---

## Wave 5 — Citations

`model=claude-opus-5`, `reasoning_effort=high`, status=completed, duration=20.5 min

### Task statement

```
Build a complete, verified reference list for a manuscript. A reviewer noted it carries only 17
references while spanning reverse cardio-oncology, HF molecular biology, ligand-receptor databases,
transcriptomic statistics, TCGA, GTEx, DepMap, CPTAC, single-cell atlases, tumour microenvironment
scoring and plasma proteomics - and that Methods cite dataset NAMES without their foundational papers.

Fetch REAL records. Check search_skills for an mcp-* PubMed connector and use host.mcp in a repl cell if
one exists; otherwise NCBI E-utilities esearch/esummary over HTTPS (allowlisted). For a contact email
call host.get_user_email() and catch host.ContactEmailUnavailable - if it raises, omit the email
parameter entirely. Never fabricate an address.

RESOURCES NEEDING THEIR FOUNDATIONAL PAPER:
  - limma (Ritchie et al. Nucleic Acids Res 2015) and the voom/limma-trend methodology (Law et al. 2014)
  - Benjamini-Hochberg FDR (1995, JRSS-B). Likely not in PubMed; if absent, cite it correctly from the
    journal and mark the record non-PubMed rather than inventing a PMID.
  - TCGA PanCanAtlas flagship, and the EBPlusPlus batch-corrected expression resource if separately cited
  - GTEx v8 (2020 Science flagship)
  - Human Protein Atlas (Uhlen et al. tissue-based map; plus the secretome paper if distinct)
  - CellChat / CellChatDB (Jin et al. Nat Commun 2021)
  - DepMap resource citation and the Chronos gene-effect algorithm (Dempster et al.)
  - CPTAC and the Proteomic Data Commons
  - Thorsson et al. immune landscape of cancer (Immunity 2018)
  - The standard human LV single-nucleus atlases (Litvinukova et al. Nature 2020; and determine whether
    Koenig et al. or Chaffin et al. is the second standard one). If you cannot establish which two the
    analysis used, cite the standard ones and state that caveat explicitly.
  - Cox proportional hazards (Cox 1972) for the survival models
  - GEO itself (Barrett et al.), and the cohort papers for GSE116250 and GSE141910 (MAGNet)

VERIFICATION IS THE POINT. For every record confirm the returned title and year actually match the work
intended, not merely that a search returned something. Report a per-record verified flag and its basis
(title similarity). If a search returns an unrelated paper, mark UNRESOLVED rather than accepting it. No
bibliographic field may be written from memory - authors, journal, volume, pages, DOI all come from the
fetched record.

EXISTING BIBLIOGRAPHY (13 of the current 17; reuse these records rather than re-fetching blindly):
  literature_citations.csv /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/c2280f91-5b72-41b2-a325-04eec06c950b/v940dda9f_literature_citations.csv
The 4 fetched later and absent from that file are PMIDs 39523983, 21685905, 35861824, 28885394.
Manuscript source, to see where each citation must land:
  PAPER.md /Users/nicolasa/.claude-science/orgs/ee427824-fe4d-4ed6-bcf9-29e31fe9ea33/artifacts/proj_272902b9683e/c60aea5b-67e0-447e-83f0-d9e025273d90/v9574fef3_PAPER.md

OUTPUTS:
  1. references_expanded.csv - every reference with pmid (or blank plus a source note for non-PubMed),
     authors, title, journal, year, volume, issue, pages, doi, verified flag, verification basis, and
     which Methods subsection or claim it supports.
  2. citation_placements.csv - for each new reference, the exact unique anchor text in PAPER.md where the
     citation should be inserted, with enough surrounding text to be unambiguous.
  3. references_note.md - anything UNRESOLVED, any record you could not verify, and any citation you
     believe the manuscript needs but could not locate.

Do NOT edit PAPER.md - a later step does placement and renumbering.

Report in structured output: total references, how many newly fetched, how many verified and by what
basis, any UNRESOLVED, the non-PubMed records and how they are cited, and all artifact ids.

BINDING CONSTRAINTS (project convention, non-negotiable):
- Fix any threshold or criterion BEFORE you compute the thing it judges, and write it into the output
  table as its own row. Do not adjust it after seeing results.
- NEVER simulate, synthesise or approximate data. If a resource is unreachable or a quantity is not
  recoverable, report BLOCKED or UNTESTABLE and say exactly which comparison could not be made.
  A blocked result reported honestly is correct; an invented one is not.
- Report findings that undermine the manuscript at full strength. The manuscript gets corrected, not
  defended. Several results in it are already withdrawals and negative findings.
- Every number must come from a computation or a fetch in this session. Nothing from memory.
- Any figure must be a final deliverable: load the figure-style skill, call apply_figure_style() first,
  use fig.savefig (never plt.savefig), and inspect the render before saving.

```

### Required structured output schema

```json
{
  "type": "object",
  "properties": {
    "total_references": {
      "type": "integer"
    },
    "newly_fetched": {
      "type": "integer"
    },
    "verified_count": {
      "type": "integer"
    },
    "verification_basis": {
      "type": "string"
    },
    "unresolved": {
      "type": "string"
    },
    "non_pubmed_records": {
      "type": "string"
    },
    "artifact_ids": {
      "type": "object"
    }
  },
  "required": [
    "total_references",
    "newly_fetched",
    "verified_count",
    "artifact_ids"
  ]
}
```
