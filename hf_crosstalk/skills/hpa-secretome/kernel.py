"""Human Protein Atlas secretome annotation.

Retrieves genome-wide protein-class / secretome annotations from the Human
Protein Atlas download API and classifies genes by secretion route.

Pure stdlib + pandas. No platform dependencies: every function takes file paths
or DataFrames and returns DataFrames / dicts / Series.
"""

import os
import re
import urllib.request

import pandas as pd

HPA_API_URL = "https://www.proteinatlas.org/api/search_download.php"

# Column codes accepted by the HPA download API. Verified against the live API
# (invalid codes are silently dropped from the response, so always inspect the
# returned header rather than trusting the request).
HPA_COLUMN_CODES = {
    "g": "Gene",
    "eg": "Ensembl",
    "gs": "Gene synonym",
    "gd": "Gene description",
    "up": "Uniprot",
    "chr": "Chromosome",
    "pc": "Protein class",
    "secl": "Secretome location + Secretome function",
    "scl": "Subcellular location",
    "scml": "Subcellular main location",
    "scal": "Subcellular additional location",
    "up_mf": "Molecular function",
    "upbp": "Biological process",
    "rnats": "RNA tissue specificity",
    "rnatd": "RNA tissue distribution",
    "blconcms": "Blood concentration MS [pg/L]",
    "blconcia": "Blood concentration IM [pg/L]",
    "di": "Disease involvement",
    "pe": "Evidence",
    "relih": "Reliability (IH)",
}

DEFAULT_HPA_COLUMNS = "g,eg,pc,secl,scml,rnats,blconcms"

# HPA "Protein class" tokens used for secretion classification.
CLASS_SECRETED = "Predicted secreted proteins"
CLASS_PLASMA = "Plasma proteins"
CLASS_MEMBRANE = "Predicted membrane proteins"
CLASS_INTRACELLULAR = "Predicted intracellular proteins"

# Minimum genome-wide row count. HPA ships ~20k protein-coding genes; anything
# far below this means the API returned an error page or a filtered subset.
MIN_GENOME_ROWS = 10000

# ---------------------------------------------------------------------------
# Structural extracellular-matrix definition.
#
# These are gene-symbol families whose products are secreted but act as
# structural scaffold rather than as diffusible signaling ligands. Flagging
# them lets a downstream analysis exclude them from "endocrine signal" claims;
# the exclusion is a judgment call and must be disclosed, hence the log.
#
# Each entry: family label -> compiled regex over HGNC gene symbols.
# ---------------------------------------------------------------------------
STRUCTURAL_ECM_FAMILIES = {
    "collagen": r"^COL\d+[A-Z]\d*$",
    "laminin": r"^LAM[ABCG]\d+$",
    "fibrillin": r"^FBN\d+$",
    "fibulin": r"^(FBLN\d+|EFEMP\d+)$",
    "elastin": r"^(ELN|MFAP\d+|EMILIN\d+)$",
    "fibronectin": r"^FN1$",
    "nidogen": r"^NID\d+$",
    "tenascin": r"^TN[CRNXW]?$|^TNXB$",
    "elastin_microfibril": r"^(LTBP\d+|FBLN5)$",
    # NOTE: keratin-associated proteins (KRTAP*) and desmosomal plakins were
    # evaluated and deliberately EXCLUDED from this rule: none of the 95 KRTAP
    # symbols in HPA carries secreted/ECM support, so they are hair-keratin
    # structural proteins rather than extracellular matrix.
    "proteoglycan_structural": r"^(AGRN|HSPG2|BGN|DCN|LUM|VCAN|ACAN|OGN|PRELP|FMOD)$",
    "matrilin_cartilage": r"^(MATN\d+|COMP|CHAD|EPYC)$",
    "collagen_like_other": r"^(EMID\d+|MULTIMERIN\d*|MMRN\d+|VWA\d+)$",
}


STRUCTURAL_ECM_RULE_TEXT = (
    "A gene is flagged is_structural_ecm if (a) its HGNC symbol matches one of "
    "the documented structural-ECM family regexes in STRUCTURAL_ECM_FAMILIES "
    "(collagens, laminins, fibrillins, fibulins, elastin/microfibril, "
    "fibronectin, nidogens, tenascins, structural proteoglycans, cartilage "
    "matrilins) AND (b) HPA independently supports an extracellular localisation "
    "for it, i.e. Protein class contains 'Predicted secreted proteins' OR "
    "Secretome location contains 'extracellular matrix'. Requirement (b) "
    "prevents intracellular symbol homonyms from being excluded. Genes matching "
    "(a) but failing (b) are recorded in the exclusion log with "
    "hpa_supports_extracellular=False and are NOT flagged."
)


def h_hpa_url(query, columns):
    from urllib.parse import urlencode

    params = {
        "search": query,
        "format": "tsv",
        "columns": columns,
        "compress": "no",
    }
    return HPA_API_URL + "?" + urlencode(params)


def fetch_hpa_annotations(
    dest_path,
    columns=None,
    query="",
    force=False,
    timeout=300,
):
    """Download (or re-read) the HPA annotation table.

    Args:
        dest_path (str): local TSV cache path. If it exists and force is False,
            it is re-read instead of re-downloaded.
        columns (str): comma-separated HPA column codes. See HPA_COLUMN_CODES.
            Must include 'g' (Gene). Default retrieves gene, Ensembl id,
            protein class, secretome location+function, subcellular main
            location, RNA tissue specificity and MS blood concentration.
        query (str): HPA search string. The empty string returns the full
            genome-wide table (~20k protein-coding genes).
        force (bool): re-download even if dest_path exists.
        timeout (int): socket timeout in seconds.

    Returns:
        pandas.DataFrame: one row per gene, indexed 0..n-1, with a 'Gene'
            column plus whichever columns HPA actually returned. All values are
            strings; missing values are the empty string.

    Raises:
        AssertionError: if the response has no 'Gene' column, or (for a
            genome-wide query) fewer than MIN_GENOME_ROWS rows. This guards
            against silently caching an HTML error page.
    """
    if columns is None:
        columns = DEFAULT_HPA_COLUMNS
    assert "g" in [c.strip() for c in columns.split(",")], (
        "columns must include 'g' (Gene); got %r" % columns
    )

    if force or not os.path.exists(dest_path) or os.path.getsize(dest_path) == 0:
        parent = os.path.dirname(os.path.abspath(dest_path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        url = h_hpa_url(query, columns)
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            payload = resp.read()
        assert payload[:1000].lstrip()[:6].lower() != b"<html>", (
            "HPA returned an HTML page, not a TSV, for url=%s" % url
        )
        with open(dest_path, "wb") as fh:
            fh.write(payload)

    df = pd.read_csv(dest_path, sep="\t", dtype=str)
    df = df.fillna("")

    assert "Gene" in df.columns, (
        "HPA response has no 'Gene' column (columns=%s); cache at %s is probably "
        "an error page - delete it and retry" % (list(df.columns), dest_path)
    )
    if query == "":
        assert len(df) >= MIN_GENOME_ROWS, (
            "genome-wide HPA fetch returned only %d rows (<%d); refusing to use "
            "a truncated table" % (len(df), MIN_GENOME_ROWS)
        )
    assert df["Gene"].str.len().gt(0).all(), "HPA table contains empty gene symbols"
    return df.reset_index(drop=True)


def h_has_class(series, token):
    """Exact token match inside HPA's comma-separated Protein class string."""
    return series.fillna("").apply(
        lambda s: token in [t.strip() for t in s.split(",")]
    )


def classify_secretion(hpa_df):
    """Derive secretion booleans and a secretion tier from HPA protein class.

    Args:
        hpa_df (pandas.DataFrame): output of fetch_hpa_annotations. Must have
            'Gene' and 'Protein class'. 'Secretome location' and
            'Secretome function' are used when present.

    Returns:
        pandas.DataFrame: one row per input row with columns
            gene (str), protein_class (str, raw HPA string),
            secretome_location (str, '' if column absent),
            secretome_function (str, '' if column absent),
            is_secreted (bool)      - HPA 'Predicted secreted proteins',
            is_plasma (bool)        - HPA 'Plasma proteins',
            is_membrane (bool)      - HPA 'Predicted membrane proteins',
            is_intracellular (bool) - HPA 'Predicted intracellular proteins',
            is_secreted_to_blood (bool) - Secretome location says blood,
            secretion_tier (str)    - 'core' if is_secreted and is_plasma,
                                      'extended' if is_secreted only,
                                      'none' otherwise.

    Raises:
        AssertionError: if required columns are missing.
    """
    for col in ("Gene", "Protein class"):
        assert col in hpa_df.columns, "hpa_df is missing required column %r" % col

    pc = hpa_df["Protein class"].fillna("")
    out = pd.DataFrame(
        {
            "gene": hpa_df["Gene"].astype(str),
            "protein_class": pc,
            "secretome_location": hpa_df.get(
                "Secretome location", pd.Series("", index=hpa_df.index)
            ).fillna(""),
            "secretome_function": hpa_df.get(
                "Secretome function", pd.Series("", index=hpa_df.index)
            ).fillna(""),
        }
    )
    out["is_secreted"] = h_has_class(pc, CLASS_SECRETED).values
    out["is_plasma"] = h_has_class(pc, CLASS_PLASMA).values
    out["is_membrane"] = h_has_class(pc, CLASS_MEMBRANE).values
    out["is_intracellular"] = h_has_class(pc, CLASS_INTRACELLULAR).values
    out["is_secreted_to_blood"] = (
        out["secretome_location"].str.contains("blood", case=False, na=False).values
    )
    out["secretion_tier"] = "none"
    out.loc[out["is_secreted"], "secretion_tier"] = "extended"
    out.loc[out["is_secreted"] & out["is_plasma"], "secretion_tier"] = "core"

    assert set(out["secretion_tier"]) <= {"core", "extended", "none"}
    assert not (
        (out["secretion_tier"] == "core") & ~out["is_secreted"]
    ).any(), "core tier must imply is_secreted"
    return out.reset_index(drop=True)


def flag_structural_ecm(genes, hpa_df=None):
    """Flag structural (scaffold) extracellular-matrix genes.

    Structural ECM proteins are secreted but act as scaffold rather than as
    diffusible signals, so an endocrine-signalling analysis usually excludes
    them. The rule is symbol-family based AND requires independent HPA support
    for extracellular localisation; see STRUCTURAL_ECM_RULE_TEXT.

    Args:
        genes (iterable of str): gene symbols.
        hpa_df (pandas.DataFrame, optional): output of fetch_hpa_annotations,
            used for requirement (b). If None, requirement (b) cannot be
            evaluated and NO gene is flagged (the function then returns an
            all-False Series and an exclusion log whose
            hpa_supports_extracellular column is None) - this is deliberate:
            the rule is not applied on symbol pattern alone.

    Returns:
        (pandas.Series, pandas.DataFrame):
            flags: boolean Series indexed by the input gene symbols.
            log: exclusion log, one row per gene matching a family regex, with
                columns gene, matched_family, matched_pattern,
                hpa_supports_extracellular (bool or None), hpa_protein_class,
                hpa_secretome_location, excluded (bool), reason (str).

    Raises:
        AssertionError: if a flagged gene is absent from the input list.
    """
    gene_list = [str(g) for g in genes]
    # Unique index (first-appearance order): HPA's genome-wide table contains a
    # small number of duplicated gene symbols, and a duplicated index makes the
    # Series unusable for .map()/.get(). Duplicates collapse to a single flag.
    uniq = list(dict.fromkeys(gene_list))
    flags = pd.Series(False, index=uniq, dtype=bool)

    support = {}
    pcmap = {}
    slmap = {}
    if hpa_df is not None:
        assert "Gene" in hpa_df.columns, "hpa_df must have a 'Gene' column"
        pc = hpa_df["Protein class"].fillna("") if "Protein class" in hpa_df else ""
        sl = (
            hpa_df["Secretome location"].fillna("")
            if "Secretome location" in hpa_df
            else pd.Series("", index=hpa_df.index)
        )
        secreted = h_has_class(
            hpa_df.get("Protein class", pd.Series("", index=hpa_df.index)),
            CLASS_SECRETED,
        )
        ecm_loc = sl.str.contains("extracellular matrix", case=False, na=False)
        sup = (secreted | ecm_loc).values
        for g, s, p, l in zip(hpa_df["Gene"].astype(str), sup, pc, sl):
            support[g] = bool(s)
            pcmap[g] = p
            slmap[g] = l

    rows = []
    for g in uniq:
        fam = None
        pat = None
        for name, rx in {name: re.compile(pat) for name, pat in STRUCTURAL_ECM_FAMILIES.items()}.items():
            if rx.match(g):
                fam, pat = name, rx.pattern
                break
        if fam is None:
            continue
        if hpa_df is None:
            sup = None
            excluded = False
            reason = "family match but no hpa_df supplied; rule not applied"
        else:
            sup = support.get(g)
            if sup is True:
                excluded = True
                reason = "structural ECM family %s + HPA extracellular support" % fam
            elif sup is False:
                excluded = False
                reason = (
                    "family %s matched but HPA gives no secreted/ECM support; kept"
                    % fam
                )
            else:
                sup = None
                excluded = False
                reason = "family %s matched but gene absent from HPA; kept" % fam
        rows.append(
            {
                "gene": g,
                "matched_family": fam,
                "matched_pattern": pat,
                "hpa_supports_extracellular": sup,
                "hpa_protein_class": pcmap.get(g, ""),
                "hpa_secretome_location": slmap.get(g, ""),
                "excluded": excluded,
                "reason": reason,
            }
        )
        if excluded:
            flags[g] = True

    log = pd.DataFrame(
        rows,
        columns=[
            "gene",
            "matched_family",
            "matched_pattern",
            "hpa_supports_extracellular",
            "hpa_protein_class",
            "hpa_secretome_location",
            "excluded",
            "reason",
        ],
    )
    assert set(flags.index[flags]) <= set(gene_list), "flagged a gene not in input"
    return flags, log


def annotate_gene_list(genes, hpa_df):
    """Annotate an arbitrary gene list against HPA secretome classification.

    Args:
        genes (iterable of str): gene symbols to annotate (duplicates are kept
            in input order; the first HPA match is used per symbol).
        hpa_df (pandas.DataFrame): output of fetch_hpa_annotations.

    Returns:
        pandas.DataFrame: one row per input gene, in input order, with columns
            gene, found_in_hpa (bool), protein_class, secretome_location,
            secretome_function, is_secreted, is_plasma, is_membrane,
            is_intracellular, is_secreted_to_blood, secretion_tier,
            is_structural_ecm. Genes absent from HPA get found_in_hpa=False,
            all booleans False and secretion_tier='none'.

    Raises:
        AssertionError: if the output row count differs from the input length.
    """
    gene_list = [str(g) for g in genes]
    cls = classify_secretion(hpa_df)
    cls = cls.drop_duplicates(subset="gene", keep="first").set_index("gene")

    flags, _log = flag_structural_ecm(gene_list, hpa_df=hpa_df)

    recs = []
    for g in gene_list:
        if g in cls.index:
            r = cls.loc[g].to_dict()
            r["found_in_hpa"] = True
        else:
            r = {
                "protein_class": "",
                "secretome_location": "",
                "secretome_function": "",
                "is_secreted": False,
                "is_plasma": False,
                "is_membrane": False,
                "is_intracellular": False,
                "is_secreted_to_blood": False,
                "secretion_tier": "none",
                "found_in_hpa": False,
            }
        r["gene"] = g
        r["is_structural_ecm"] = bool(flags.get(g, False))
        recs.append(r)

    out = pd.DataFrame(recs)[
        [
            "gene",
            "found_in_hpa",
            "protein_class",
            "secretome_location",
            "secretome_function",
            "is_secreted",
            "is_plasma",
            "is_membrane",
            "is_intracellular",
            "is_secreted_to_blood",
            "secretion_tier",
            "is_structural_ecm",
        ]
    ]
    assert len(out) == len(gene_list), "annotate_gene_list changed row count"
    return out.reset_index(drop=True)


def build_secretome_table(hpa_df):
    """Convenience: genome-wide secretome table with the structural-ECM flag.

    Args:
        hpa_df (pandas.DataFrame): output of fetch_hpa_annotations (genome-wide).

    Returns:
        (pandas.DataFrame, pandas.DataFrame): (table, ecm_log) where table has
            one row per HPA gene with gene, protein_class, secretome_location,
            secretome_function, is_secreted, is_plasma, is_membrane,
            is_intracellular, is_secreted_to_blood, secretion_tier,
            is_structural_ecm; and ecm_log is the flag_structural_ecm log.
    """
    cls = classify_secretion(hpa_df)
    flags, log = flag_structural_ecm(cls["gene"].tolist(), hpa_df=hpa_df)
    cls["is_structural_ecm"] = cls["gene"].map(flags).fillna(False).astype(bool)
    return cls, log

# ---------------------------------------------------------------------------
# GTEx tissue-specificity companion (added after the HF-secretome screen)
# ---------------------------------------------------------------------------

GTEX_MEDIAN_URL = "https://gtexportal.org/api/v2/expression/medianGeneExpression"
GTEX_GENE_URL = "https://gtexportal.org/api/v2/reference/gene"
GTEX_DATASET_ID = "gtex_v8"
GTEX_GENCODE_VERSION = "v26"
GTEX_GENOME_BUILD = "GRCh38/hg38"
GTEX_CELL_LINES = ("Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes")
GTEX_API_GOTCHA = (
    "GTEx API v2 requires BOTH datasetId=gtex_v8 AND gencode IDs versioned to "
    "that release (resolve via /reference/gene with gencodeVersion=v26 and "
    "genomeBuild=GRCh38/hg38). Omitting datasetId, or passing a bare unversioned "
    "ENSG, returns HTTP 200 with an EMPTY data array and totalNumberOfItems=0 -- "
    "a silent no-op that looks like 'gene not expressed' rather than an error. "
    "Always assert a non-zero row count per gene before interpreting."
)
RATIO_METRIC_CAVEAT = (
    "An enrichment ratio of tissue-of-interest TPM divided by the cross-tissue "
    "MEDIAN rewards genes that are near-zero in most tissues, regardless of where "
    "their actual maximum sits. A liver-dominant protein can score a high cardiac "
    "ratio purely because its median across tissues is ~0. Never report a ratio "
    "without the companion absolute-rank and fold-vs-top-tissue columns."
)


def gtex_resolve_gencode_ids(symbols, batch=50, timeout=120, retries=4):
    """Map HGNC symbols -> release-matched GTEx gencode IDs (see GTEX_API_GOTCHA)."""
    import json
    import time
    import urllib.parse
    out = {}
    syms = [str(s).upper() for s in symbols]
    for i in range(0, len(syms), batch):
        chunk = syms[i:i + batch]
        q = [("geneId", s) for s in chunk] + [
            ("gencodeVersion", GTEX_GENCODE_VERSION),
            ("genomeBuild", GTEX_GENOME_BUILD),
            ("itemsPerPage", "500"),
        ]
        url = GTEX_GENE_URL + "?" + urllib.parse.urlencode(q)
        for attempt in range(retries):
            try:
                req = urllib.request.Request(url, headers={"Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    payload = json.loads(resp.read().decode())
                for rec in payload.get("data", []):
                    sym = str(rec.get("geneSymbol", "")).upper()
                    if sym and sym not in out:
                        out[sym] = rec.get("gencodeId")
                break
            except Exception:
                if attempt == retries - 1:
                    raise
                time.sleep(2 * (attempt + 1))
    return out


def gtex_fetch_tissue_medians(gencode_ids, per_request=25, timeout=180, retries=4):
    """Fetch median TPM per tissue. Asserts the request was not a silent no-op."""
    import json
    import time
    import urllib.parse
    ids = [g for g in gencode_ids if g]
    assert ids, "no gencode IDs supplied -- resolve symbols first"
    rows = []
    for i in range(0, len(ids), per_request):
        chunk = ids[i:i + per_request]
        q = [("gencodeId", g) for g in chunk] + [
            ("datasetId", GTEX_DATASET_ID),
            ("itemsPerPage", "5000"),
        ]
        url = GTEX_MEDIAN_URL + "?" + urllib.parse.urlencode(q)
        for attempt in range(retries):
            try:
                req = urllib.request.Request(url, headers={"Accept": "application/json"})
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    payload = json.loads(resp.read().decode())
                rows.extend(payload.get("data", []))
                break
            except Exception:
                if attempt == retries - 1:
                    raise
                time.sleep(3 * (attempt + 1))
    assert rows, (
        "GTEx returned zero rows for every gene. " + GTEX_API_GOTCHA
    )
    frame = pd.DataFrame(rows)
    n_genes = frame["geneSymbol"].nunique() if "geneSymbol" in frame.columns else 0
    assert n_genes > 0, "GTEx response carried no geneSymbol column"
    return frame


def gtex_tissue_matrix(medians_df, drop_cell_lines=True):
    """Pivot long-form GTEx medians into a gene x tissue TPM matrix."""
    need = {"geneSymbol", "tissueSiteDetailId", "median"}
    missing = need - set(medians_df.columns)
    assert not missing, "medians_df missing columns: %s" % sorted(missing)
    mat = medians_df.pivot_table(
        index="geneSymbol", columns="tissueSiteDetailId", values="median", aggfunc="median"
    )
    if drop_cell_lines:
        mat = mat.drop(columns=[c for c in GTEX_CELL_LINES if c in mat.columns])
    return mat


def tissue_specificity(matrix, tissues, min_ratio=3.0, max_rank=5, min_tpm=5.0,
                       max_fold_below_top=3.0):
    """Score enrichment for `tissues` with the rank guard RATIO_METRIC_CAVEAT requires.

    Returns one row per gene with the ratio-vs-median metric AND the companion
    columns that catch its failure mode: absolute rank of the best target tissue,
    the top tissue overall, and fold difference between top tissue and target.
    `passes` requires enrichment, a top-`max_rank` placement, AND abundance.
    `ratio_only_artifact` is True whenever a gene clears the ratio threshold but
    the target tissue is not the global maximum, or sits >=`max_fold_below_top`
    below it -- the APOA1/PI16 pattern in RATIO_METRIC_CAVEAT. Inspect those rows
    by hand; a True here means the ratio is arithmetically real but the
    tissue-of-origin reading it invites is not.
    """
    import numpy as np
    tissues = [t for t in tissues if t in matrix.columns]
    assert tissues, "none of the requested tissues are columns in the matrix"
    out = []
    for gene in matrix.index:
        row = matrix.loc[gene].dropna()
        present = [t for t in tissues if t in row.index]
        if not present:
            continue
        ranked = row.sort_values(ascending=False)
        best_t = max(present, key=lambda t: float(row[t]))
        best = float(row[best_t])
        med = float(row.median())
        top_t = ranked.index[0]
        top_v = float(ranked.iloc[0])
        ratio = best / med if med > 0 else np.inf
        rank = int(ranked.index.get_loc(best_t)) + 1
        enriched = (ratio >= min_ratio) and (rank <= max_rank)
        out.append(dict(
            gene=gene, best_tissue=best_t, best_tpm=best,
            median_across_tissues=med, ratio_vs_median=ratio,
            best_tissue_rank=rank, top_tissue=top_t, top_tissue_tpm=top_v,
            fold_below_top=(top_v / best) if best > 0 else np.inf,
            n_tissues=len(row), enriched=enriched, abundant=best >= min_tpm,
            passes=bool(enriched and best >= min_tpm),
            ratio_only_artifact=bool(
                ratio >= min_ratio
                and (rank > max_rank
                     or rank > 1
                     or (top_v / best if best > 0 else np.inf) >= max_fold_below_top)
            ),
        ))
    frame = pd.DataFrame(out)
    if len(frame):
        assert frame.best_tissue_rank.min() >= 1, "rank must be 1-based"
    return frame


def base_rate_guard(hit_flags, background_flags):
    """Fisher test: is the pass rate in the hit set above the background base rate?

    Both args are boolean Series/arrays of the SAME criterion applied to the
    candidate set and to a matched non-candidate set. A screen whose criterion
    passes at background rate has found nothing, however plausible the hits look.
    """
    from scipy.stats import fisher_exact
    hit = [bool(x) for x in hit_flags]
    bg = [bool(x) for x in background_flags]
    a = sum(hit)
    b = len(hit) - a
    c = sum(bg)
    d = len(bg) - c
    odds, pval = fisher_exact([[a, b], [c, d]])
    return dict(
        hit_pass=a, hit_n=len(hit), hit_rate=(a / len(hit) if hit else float("nan")),
        bg_pass=c, bg_n=len(bg), bg_rate=(c / len(bg) if bg else float("nan")),
        odds_ratio=float(odds), p_value=float(pval),
        enriched_above_base_rate=bool(pval < 0.05 and odds > 1),
    )

