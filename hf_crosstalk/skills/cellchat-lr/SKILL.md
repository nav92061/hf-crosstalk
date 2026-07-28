---
name: cellchat-lr
description: Fetch and parse CellChatDB (human) ligand-receptor interactions from the CellChat GitHub .rda, expand multi-subunit receptor and ligand complexes to per-gene subunits, and map a ligand gene list to the receptors that could receive it. Use when you have candidate secreted ligands and need the receptor genes, pathway names and interaction annotation (Secreted Signaling / ECM-Receptor / Cell-Cell Contact) to test receptor expression in a target cell population.
---

# cellchat-lr

## Purpose

Turn a list of candidate ligand genes into a table of **which receptor genes
could receive them**, with pathway names, evidence strings, and the interaction's
mechanistic class. Multi-subunit receptors are expanded to one row per subunit
gene so the table joins directly against an expression matrix keyed on gene
symbol.

## Data source

```
https://raw.githubusercontent.com/sqjin/CellChat/master/data/CellChatDB.human.rda
```

(~756 KB, bzip2-compressed R serialization; mouse equivalent at
`CellChatDB.mouse.rda`, exposed as `CELLCHATDB_MOUSE_URL`.)

Parsed in pure Python with the `rdata` package — no R installation is required.
The `.rda` holds a single top-level named object which is an R **list** of four
data.frames:

| element | shape (human, version fetched here) | notes |
|---|---|---|
| `interaction` | 1939 × 11 | one row per L-R interaction |
| `complex` | 157 × 4 (+1 after load) | wide: `subunit_1..subunit_4`; `load_cellchatdb` promotes the R row names to a `complex_name` column, so the returned frame is 157 × 5 |
| `cofactor` | 31 × 16 (+1 after load) | wide: `cofactor1..cofactor16`; returned as 31 × 17 with the row-name column |
| `geneInfo` | 41787 × 6 (+1 after load) | HGNC symbol / Entrez / Ensembl mapping; returned as 41787 × 7 |

**Structural notes that cost time if you rediscover them:**

- `complex` and `cofactor` carry their key **only as R row names**, not as a
  column. `load_cellchatdb` promotes those to `complex_name` / `cofactor_name`
  (and `geneInfo` row names, which are HGNC IDs, to `hgnc_id`).
- Unused subunit/cofactor slots are **empty strings, not NA**, so a naive
  `dropna()` keeps them. `tidy_complexes` filters on empty strings.
- Complex *names* are free text, not gene symbols — `TGFbR1_R2`, `Activin AB`,
  `ACVR1_TGFbR`. Never try to split a receptor name on `_` to recover subunits;
  join to the complex table.
- `interaction.receptor` / `.ligand` hold either a plain gene symbol or a complex
  name. A name absent from the complex table is its own single subunit.
- Interaction row names duplicate `interaction_name`, which already exists as a
  column.

## Annotation categories and why the default filter is `Secreted Signaling`

| annotation | interactions | mechanism |
|---|---|---|
| `Secreted Signaling` | 1199 | diffusible ligand → receptor |
| `ECM-Receptor` | 421 | matrix-deposited ligand → receptor |
| `Cell-Cell Contact` | 319 | requires physical juxtaposition |

`map_ligands_to_receptors` defaults to `annotation="Secreted Signaling"` because
only that class supports a ligand acting **at a distance** from its tissue of
origin (a circulating-factor mechanism). ECM-Receptor requires the ligand to be
deposited in the local matrix around the receiving cell, and Cell-Cell Contact
requires the two cells to touch — neither is available to a factor released into
blood by a distant organ.

Both remain fully retrievable: pass `annotation=None` for everything, or a
string / list of category names. Retrieve them when the question is local
signalling (e.g. within the failing myocardium, or within the tumour
microenvironment), and disclose which class the claim rests on.

## Functions

```python
fetch_cellchatdb(dest_path, url=CELLCHATDB_HUMAN_URL, force=False, timeout=300) -> str
load_cellchatdb(path) -> dict[str, pd.DataFrame]
tidy_complexes(complex_df) -> pd.DataFrame
expand_receptor_complexes(interaction_df, complex_df, side="receptor") -> pd.DataFrame
expand_all_subunits(interaction_df, complex_df) -> pd.DataFrame
map_ligands_to_receptors(ligand_genes, db, annotation="Secreted Signaling",
                         expand_ligand_complexes=True) -> pd.DataFrame
summarize_db(db) -> dict
```

- `fetch_cellchatdb` is cache-aware and returns the path.
- `load_cellchatdb` returns `{'interaction', 'complex', 'cofactor', 'geneInfo'}`,
  all values coerced to stripped strings.
- `tidy_complexes` melts wide → long: `complex_name`, `subunit_index` (1-based),
  `subunit_gene`, `n_subunits`.
- `expand_receptor_complexes` accepts the complex table in **either** wide or
  tidy form and adds `<side>_subunit_index`, `<side>_subunit_gene`,
  `n_subunits`, `is_complex`. `side` is `"receptor"` or `"ligand"`.
- `expand_all_subunits` expands both sides, renaming counts to
  `n_ligand_subunits` / `n_receptor_subunits` and
  `ligand_is_complex` / `receptor_is_complex`.
- `map_ligands_to_receptors` matches a ligand if the interaction's ligand name is
  in the input list **or** (default) any subunit of a ligand complex is. This
  matters: `INHBA` is the ligand name of 2 interactions directly, but reaches 7
  once `Activin AB` and `Inhibin A` complex membership is expanded.
- `summarize_db` returns descriptive counts including `by_annotation`.

## Worked example

```python
import kernel as cc

path = cc.fetch_cellchatdb("data/CellChatDB.human.rda")
db = cc.load_cellchatdb(path)
cc.summarize_db(db)["by_annotation"]
# {'Secreted Signaling': 1199, 'ECM-Receptor': 421, 'Cell-Cell Contact': 319}

# per-subunit table for ALL interactions - join target for expression matrices
subunits = cc.expand_receptor_complexes(db["interaction"], db["complex"])
len(subunits)          # 2884 rows from 1939 interactions

# candidate circulating ligands -> receivable receptors
hits = cc.map_ligands_to_receptors(["NPPB", "IL6", "TGFB1", "GDF15", "INHBA"], db)
hits[["interaction_name", "pathway_name", "matched_ligand_gene",
      "receptor", "receptor_subunit_gene", "n_receptor_subunits"]].head()
# TGFB1_TGFBR1_TGFBR2 | TGFb | TGFB1 | TGFbR1_R2 | TGFBR1 | 2
# TGFB1_TGFBR1_TGFBR2 | TGFb | TGFB1 | TGFbR1_R2 | TGFBR2 | 2

# a complex receptor is only "expressed" if ALL subunits are - test per interaction
expr = set(tumour_expressed_genes)
ok = (hits.groupby("interaction_name")["receptor_subunit_gene"]
          .apply(lambda s: set(s) <= expr))
usable = ok[ok].index

# local (non-endocrine) signalling instead:
cc.map_ligands_to_receptors(["COL1A1"], db, annotation="ECM-Receptor")
cc.map_ligands_to_receptors(["COL1A1"], db, annotation=None)
```

## Preconditions / postconditions

Preconditions enforced as raising assertions:

- Downloaded file is > 100 KB and begins with compressed-R / RDX magic bytes —
  catches an HTML error page saved as `.rda`.
- The parsed object is an R list containing all four of `interaction`,
  `complex`, `cofactor`, `geneInfo`; each is a data.frame.
- The interaction table has ≥ `MIN_INTERACTION_ROWS` (1000) rows and carries all
  of `interaction_name`, `pathway_name`, `ligand`, `receptor`, `annotation`,
  `evidence`.
- `complex_df` has a `complex_name` column and at least one `subunit*` column.
- `side` is `"receptor"` or `"ligand"`.
- An `annotation` value not present in the database raises rather than silently
  returning an empty frame.

Postconditions asserted:

- Subunit expansion never drops an interaction (`interaction_name.nunique()`
  is preserved) and never emits an empty subunit gene.
- Melted complexes contain no empty `subunit_gene`.

Interpretation caveats:

- CellChatDB is curated from literature/KEGG; absence of a pair is not evidence
  against it. `evidence` holds the supporting citation (e.g. `KEGG: hsa04350`).
- A multi-subunit receptor requires **all** subunits present to signal. After
  expansion, aggregate per `interaction_name` before calling a receptor
  expressed — do not treat single-subunit rows as independent hits.
- Interactions sharing a ligand and pathway are partly redundant (many TGFb rows
  differ only in receptor pairing); count pathways, not rows, when describing
  breadth.

## Standalone use (no agent)

`kernel.py` depends only on the standard library, pandas, and `rdata`
(`pip install rdata`), and references no platform APIs.

```bash
python - <<'PY'
import sys; sys.path.insert(0, "skills/cellchat-lr")
import kernel as cc
db = cc.load_cellchatdb(cc.fetch_cellchatdb("CellChatDB.human.rda"))
db["interaction"].to_csv("cellchatdb_interactions.csv", index=False)
cc.tidy_complexes(db["complex"]).to_csv("cellchatdb_complexes.csv", index=False)
cc.expand_receptor_complexes(db["interaction"], db["complex"]).to_csv(
    "cellchatdb_receptor_subunits.csv", index=False)
print(cc.summarize_db(db))
PY
```

Requires network access to `raw.githubusercontent.com` on first run only.
`pyreadr` is an alternative reader but does **not** handle nested R lists
cleanly; `rdata` is the working path and is what `load_cellchatdb` uses.
