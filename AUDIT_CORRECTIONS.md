# Numeric audit corrections

Removed from the manuscript body during final drafting: these are internal quality-control
records, not results. They are retained here in full because each correction changed a
number that had already been written down, and several made the paper's claims weaker.
Cross-referenced to `validation_summary.csv`.

Every numeric claim in this manuscript was machine-checked against the saved
artifacts. Two audits were run, and the second used a deliberately different method
from the first: rather than checking prose against the summary artifacts, it
recomputed from the underlying source tables wherever an artifact stored only a
summary. That method difference is what found the second batch: all five of the
numerical errors below arose from quoting a summary sentence instead of the column
the sentence describes.

**First audit: six errors, logged as rows 26–31 of `validation_summary.csv`.**

1. The floor drops four PENK→opioid-receptor axes (OPRD1, OPRK1, OPRL1, OPRM1), not three. The whole OPIOID pathway is non-viable.
2. Four of 30 tumour types beat the specificity null at p<0.05, not six. Ovarian carcinoma ranks second by score but does not pass (p=0.060); rank and receptor-specificity are not the same thing.
3. GAPDH's median across per-tumour-type medians is 69,134, not ≈64,000. The floor itself (345.7) was always computed rather than typed, so no result changes.
4. `compbio.cn` is grantable; an earlier note claimed TIMER2 was "not grantable" having probed a different host and never requested the real one. TIMER2 is nonetheless unusable here for a different reason (Shiny app, no bulk endpoint).
5. A documented scale-check example in `geo-bulk-de` reported output values that a real run does not produce (min/max 1.00/22.39 versus the actual 1.82/21.03). The qualitative conclusion it supports, that this "RAW" series ships log₂-scale values rather than counts, is unaffected and was independently confirmed.
6. The claim that the skill catalog held "only Anthropic built-ins" was overstated; one unrelated personal skill pre-existed. None of the six analysis skills did.

**Second audit: 61 claims re-checked, 49 verified exactly, five numerical errors
plus seven attribution or currency findings, logged as rows 51–58:**

7. Adjacent-normal samples: 680, not 710. The larger figure counts samples mapped to a tumour type *before* the ≥10-normals-per-type filter; 30 of them fall in types subsequently dropped.
8. NPR1 fails the expression floor in 24 of 30 types, not 26. It passes in KIRC, KIRP, MESO, OV, PCPG and SARC. The conclusion is unaffected: NPR1 still fails in the large majority.
9. DepMap positive controls are now reported with their aggregation named. The paper's −2.53/−2.69 are pan-cell-line means; the per-lineage means that enter the correlation are −2.55/−2.67.
10. CellChatDB ligand overlap is 70, not 72: 70 across all annotation classes, 50 within `Secreted Signaling`. The downstream chain is unaffected: the Secreted-Signaling restriction gives 50 ligands, expanding to the 51 network subunits.
11. WNT contributes 40 of the 160 floor-passing interactions, not 80. The 80 is a row count at interaction × receptor-subunit granularity. WNT remains the largest pathway (BMP 21, CCL 15).
12. Two further findings were real numbers attached to the wrong population, now restated: the survival denominator (146 ligand–receptor-subunit axes, not 160 complex-level interactions) and the FZD7 housekeeping-independence statistic (the 36-line DepMap Uterus pool, not the 16-line panel).

A reproducibility defect was also fixed rather than logged as a numeric error: three
filenames existed as two distinct artifacts each, and `TABLE_INDEX.csv` pointed at the
earlier version, which lacked the columns Section 3.7's claims depend on. The index now names
the versions the results actually use.

None of these alter the study's conclusions. Two make claims weaker than first written
(first-audit item 2 on the ranking, and second-audit item 7 on sample count), and one
strengthens the paper's methodological claim by adding the base-rate guard that Section 3.9
required.
