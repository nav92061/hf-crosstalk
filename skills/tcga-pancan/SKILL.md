---
name: tcga-pancan
description: Memory-bounded extraction of TCGA PanCanAtlas pan-cancer expression and tumor microenvironment state from the GDC. Streams the ~1.9 GB EBPlusPlus RNA-seq matrix row-wise to pull selected genes without loading it, maps TCGA aliquot barcodes to tumor type, filters to primary solid tumors, and assembles per-sample leukocyte/stromal/purity and PanImmune signature scores. Use when you need per-tumor-type expression of a gene panel (receptors, ligands, housekeeping references) in true published units, or TCGA immune/stromal composition, on a machine with only a few GiB of RAM.
---

# tcga-pancan

Pull pan-cancer expression and tumor microenvironment state out of the TCGA
PanCanAtlas / PanImmune data products hosted on the GDC, on a laptop-sized
machine.

## Purpose

The PanCanAtlas expression matrix (`EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv`)
is ~1.9 GB of plain TSV with 20,531 genes x 11,069 aliquots. Loading it with
`pd.read_csv` needs well over 8 GiB. This skill streams it line by line, keeps
only the rows you asked for, and never holds more than the retained submatrix.

Two properties matter for downstream scoring and are handled explicitly:

1. **Units are real, not z-scores.** Values are RSEM `normalized_count`
   (upper-quartile scaled, EBPlusPlus batch-corrected). They are **not TPM and
   not FPKM** — `kernel.EXPRESSION_UNIT` carries the exact string; do not
   relabel them. Because the scale is absolute, a gene that is not transcribed
   is distinguishable from one that is merely low relative to other tumors.
2. **Housekeeping references travel with the panel.** `housekeeping_medians()`
   reports medians for 14 housekeeping genes in the same unit, so an absolute
   expression floor can be calibrated (e.g. "below 1% of the GAPDH median")
   rather than inferred from a within-gene z-score, which cannot express
   "absent".

## GDC endpoints used (all open access, no token)

| Quantity | File | UUID |
|---|---|---|
| Expression | `EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv` | `3586c0da-64d0-4b74-a449-5ff4d9136611` |
| Aliquot → cancer type + QC | `merged_sample_quality_annotations.tsv` | `1a7d7be8-675d-4e60-a105-19d4121bdebf` |
| Leukocyte fraction (Thorsson) | `TCGA_all_leuk_estimate.masked.20170107.tsv` | `6f75c9d7-5134-4ed1-b8f3-72856c98a4e8` |
| 160 immune signature scores | `Scores_160_Signatures.tsv.gz` | `80a82092-161d-4615-9d96-e858f113618d` |
| ABSOLUTE purity/ploidy | `TCGA_mastercalls.abs_tables_JSedit.fixed.txt` | `4f277128-f793-4354-a13d-30cc7fe9f6b5` |
| Clinical endpoints | TCGA-CDR `.xlsx` | `1b5f413e-a8d1-4d10-92eb-7c4ae739ed81` |

Base URL is `https://api.gdc.cancer.gov/data/<uuid>`. These are the files listed
at `https://gdc.cancer.gov/about-data/publications/panimmune`, i.e. the
published Thorsson et al. immune-landscape data products. **The TIMER2 web
server is not used** — where TIMER2 would supply immune/stromal deconvolution,
this skill uses the Thorsson leukocyte fraction plus ABSOLUTE purity (and their
difference, the stromal fraction) and the PanImmune signature scores. These are
the source data behind the published immune landscape, but they are not
identical to TIMER2 output: they are methylation-based (leukocyte fraction) and
DNA-copy-number-based (purity), not TIMER2's expression deconvolution, and they
do not include TIMER2's per-cell-type abundance estimates.

## Functions

```python
# Downloads (all cache to disk; re-download only with force=True)
download_to_disk(url, dest, force=False, timeout=120, chunk_bytes=8<<20) -> str
fetch_thorsson_immune(dest, force=False) -> DataFrame        # aliquot -> cancer type + QC
fetch_leukocyte_fraction(dest, force=False) -> DataFrame      # tumor_type, aliquot_barcode, leukocyte_fraction
fetch_absolute_purity(dest, force=False) -> DataFrame         # ABSOLUTE purity / ploidy
fetch_signature_scores(dest, wanted_signatures=None, force=False) -> DataFrame  # samples x signatures
fetch_tcga_cdr(dest, force=False, sheet=0) -> DataFrame       # clinical endpoints

# Streaming expression extraction (memory-bounded)
stream_expression_rows(url_or_path, wanted_genes, dest_csv, chunk_bytes=1<<20,
                       gene_id_sep="|", progress_every=5000, verbose=True) -> dict
load_expression_csv(path, index_col="gene") -> DataFrame
stream_gene_medians(url_or_path, keep_samples=None, gene_id_sep="|",
                    min_symbol_len=2, progress_every=0, verbose=False) -> Series
sample_background_genes(url_or_path, n=200, seed=0, exclude=(), min_symbol_len=2) -> dict
sample_background_genes_stratified(gene_medians, n=200, seed=0, exclude=(), n_strata=10) -> dict

# Symbol vocabulary (the matrix is frozen at ~2016 HGNC)
resolve_symbols_to_matrix(url_or_path, symbols, gene_info_path, gene_id_sep="|") -> dict

# Barcodes and sample mapping
parse_tcga_barcode(barcode) -> dict
build_sample_map(thorsson_df, aliquot_col=None, study_col=None) -> DataFrame
lookup_tumor_types_via_gdc(barcodes, page_size=400, timeout=90, verbose=True) -> dict

# Filtering / summarisation
filter_solid_primary(expr, sample_map, exclude_studies=("LAML","DLBC","THYM"),
                     require_all_mapped=True) -> (DataFrame, DataFrame)
summarize_by_tumor_type(expr, sample_map, stat="median") -> (DataFrame, DataFrame)
housekeeping_medians(expr, genes=HOUSEKEEPING_GENES, sample_map=None) -> dict

# TME assembly
build_tme_table(leukocyte_df, purity_df=None, signature_df=None, sample_map=None) -> DataFrame
summarize_tme_by_tumor_type(tme_df, tumor_type_col="tumor_type", stat="median") -> DataFrame

# Gene sets
parse_cellchatdb_genes(rda_path) -> dict   # needs `rdata`; expands complexes/cofactors
```

Constants: `EXPRESSION_UNIT`, `GDC_*_URL`, `TCGA_SAMPLE_TYPES`,
`TCGA_STUDY_NAMES`, `DEFAULT_EXCLUDED_STUDIES`, `HOUSEKEEPING_GENES`.

## Worked example

```python
import kernel as K

# 1. gene panel: CellChatDB receptors + ligands, plus housekeeping references
cc = K.parse_cellchatdb_genes("CellChatDB.human.rda")
wanted = sorted(set(cc["union"]) | set(K.HOUSEKEEPING_GENES))

# 2. stream the 1.9 GB matrix; only `wanted` rows are held in memory
info = K.stream_expression_rows(K.GDC_EXPRESSION_URL, wanted, "expr_subset.csv.gz")
print(info["genes_found"][:5], len(info["genes_missing"]), info["expression_unit"])
expr = K.load_expression_csv("expr_subset.csv.gz")

# 3. map aliquots to tumor type and keep primary solid tumors only
ann   = K.fetch_thorsson_immune("sample_annotations.tsv")
smap  = K.build_sample_map(ann)
expr, smap = K.filter_solid_primary(expr, smap, require_all_mapped=False)

# 4. per-tumor-type medians + the absolute floor calibration
by_type, n_table = K.summarize_by_tumor_type(expr, smap, stat="median")
hk = K.housekeeping_medians(expr, sample_map=smap)
floor = 0.01 * hk["overall"]["GAPDH"]          # example absolute floor
silent = by_type.index[(by_type < floor).all(axis=1)]

# 5. tumor microenvironment state
lf   = K.fetch_leukocyte_fraction("leuk.tsv")
pur  = K.fetch_absolute_purity("absolute.txt")
sig  = K.fetch_signature_scores("sig160.tsv.gz",
                                wanted_signatures=["TGFB_score_21050467",
                                                   "CHANG_CORE_SERUM_RESPONSE_UP",
                                                   "Module11_Prolif_score"])
tme  = K.build_tme_table(lf, purity_df=pur, signature_df=sig, sample_map=smap)
tme_by_type = K.summarize_tme_by_tumor_type(tme)
```

## Preconditions (enforced as raising assertions)

- `stream_expression_rows`: `wanted_genes` must be non-empty; the header must
  have sample columns; aliquot barcodes in the header must be unique; **at least
  one requested gene must match** — a silent zero-row result is impossible.
- `parse_tcga_barcode`: barcode must have ≥3 hyphen-delimited fields.
- `build_sample_map`: an aliquot column and a cancer-type column must be
  locatable; the resulting index must contain no duplicate aliquot barcodes.
- `filter_solid_primary`: every expression column must have a tumor-type
  assignment (or `require_all_mapped=False` to drop the unmapped ones); at
  least one primary solid tumor must survive.
- `summarize_by_tumor_type`: `stat` must be one of median/mean/q25/q75/frac_pos;
  every expression column must be present in the sample map.
- `parse_cellchatdb_genes`: must yield >0 receptor genes.
- `sample_background_genes`: the eligible symbol universe must be ≥ `n`.

## Postconditions

- The CSV written by `stream_expression_rows` has one row per matched gene
  (symbol in column `gene`), one column per aliquot barcode, values in
  `EXPRESSION_UNIT`. Symbol duplicates in the source keep the first occurrence
  and are listed in `duplicate_symbols`.
- `filter_solid_primary` returns `expr` and `sample_map` aligned on identical,
  identically ordered aliquot barcodes.
- `summarize_by_tumor_type` returns tumor-type columns sorted alphabetically,
  and an `n_table` whose `n_samples` sum equals the number of expression columns.
- `build_tme_table` is indexed by 15-character TCGA **sample** barcode, not
  aliquot: ABSOLUTE calls come from DNA aliquots while expression comes from RNA
  aliquots of the same sample, so joins must happen at sample depth.
  `stromal_fraction = 1 - tumor_purity - leukocyte_fraction`, clipped at 0 and
  NaN where purity is missing.

## The symbol-drift trap (read before trusting a "gene absent" result)

The matrix carries ~2016 HGNC symbols. A modern gene list will report genes as
absent that are in fact present under an old name — in a 1,053-gene CellChatDB
panel, 23 of the 36 apparent misses were renames, including **CXCL8 stored as
IL8, ACKR1–4 as DARC/CCBP2/CXCR7/CCRL1, NECTIN1–4 as PVRL1–4, VSIR as
C10orf54, VEGFD as FIGF**. Falsely calling a chemokine receptor "not
transcribed" is exactly the error an absolute-floor analysis must not make, so
run the resolver before interpreting misses:

```python
res = K.resolve_symbols_to_matrix(K.GDC_EXPRESSION_URL, wanted,
                                  "Homo_sapiens.gene_info.gz")
# res["resolved"]: current symbol -> matrix symbol, matched on exact Entrez GeneID
wanted2 = sorted(set(wanted) | set(res["resolved"].values()))
info = K.stream_expression_rows(K.GDC_EXPRESSION_URL, wanted2, "expr.csv.gz")
expr = K.load_expression_csv("expr.csv.gz").rename(
    index={v: k for k, v in res["resolved"].items()})   # back to current symbols
```

`gene_info` comes from
`https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz`.
Matching is on GeneID, not string similarity, so it cannot produce a wrong pairing.

## Background pools for a specificity null

`sample_background_genes` draws uniformly over the transcriptome. That pool is
dominated by low- and mid-expressed genes: in the retained solid-tumor cohort a
200-gene uniform draw topped out at a median of ~5,700 while ACTB sits near
94,000, so it offers **no candidates at housekeeping level** and an
expression-matched null for a highly expressed receptor has nothing to match
against. Use `stream_gene_medians` + `sample_background_genes_stratified`
(equal-count expression deciles) when the null must be expression-matched.

## Caveats that change interpretation

- **RSEM normalized_count, not TPM.** Values are comparable across samples
  within this matrix but are not transcript-per-million; do not feed them to a
  tool expecting TPM without saying so.
- **Bulk tissue.** A tumor-type median mixes malignant, immune, stromal and
  endothelial cells. Carry compartment markers (PTPRC, PECAM1, COL1A1, ACTA2,
  EPCAM) in the panel so a "high receptor" call can be checked against the
  possibility that it is simply a leukocyte- or stroma-rich cohort.
- **LAML/DLBC/THYM excluded by default.** LAML and DLBC tumor aliquots are
  blood/marrow, so their expression reflects circulating leukocytes rather than
  a perfused solid tumor; THYM is dominated by non-neoplastic thymocytes.
- **`merged_sample_quality_annotations.tsv` carries QC flags** (`Do_not_use`,
  `AWG_excluded_because_of_pathology`) that this skill does not apply
  automatically — filter on them yourself if a stricter cohort is wanted.
- **Per-sample immune subtype labels (C1–C6) are not in the GDC PanImmune file
  set.** The signature scores that define them are (`Scores_160_Signatures`),
  but the mclust subtype calls are published only as a journal supplementary
  table. If you need C1–C6 labels, supply that table yourself.

## Standalone use (no agent)

`kernel.py` depends only on the standard library plus pandas/numpy (and `rdata`
for `parse_cellchatdb_genes`, `openpyxl` for `fetch_tcga_cdr`). It makes no
reference to any agent host, platform SDK or artifact store. Copy the file
anywhere and run it:

```bash
pip install pandas numpy rdata openpyxl
python - <<'PY'
import kernel as K
info = K.stream_expression_rows(K.GDC_EXPRESSION_URL, ["GAPDH","ACTB","OXTR"],
                                "small.csv.gz")
print(info["genes_found"], info["genes_missing"], info["n_samples"])
PY
```

Every function takes file paths or DataFrames and returns DataFrames or plain
dicts, so the whole pipeline can be scripted from a shell, a Makefile, or
another language via CSV/JSON hand-off. Peak memory for a 1,000-gene panel over
11,000 aliquots is a few hundred MB.
