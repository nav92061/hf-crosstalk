"""depmap-lineage: resolve, fetch and analyse DepMap public release data.

Pure stdlib + pandas/numpy/scipy/requests. No platform (`host.*`) calls, so this
module runs unchanged in a plain Python session.

Why figshare and not depmap.org: depmap.org serves a bot-verification
interstitial to non-browser clients, so its download API returns HTML rather
than CSV. DepMap public quarterly releases are also deposited on figshare+ as
one article per release ("DepMap 24Q4 Public"), and ndownloader.figshare.com
serves the files directly. All resolution here goes through
https://api.figshare.com/v2/.
"""

import json
import os
import re
import time

import numpy as np
import pandas as pd
import requests
from scipy import stats

FIGSHARE_API = "https://api.figshare.com/v2"
RELEASE_TITLE_RE = r"^DepMap\s+(\d\dQ\d)\s+Public$"
CORE_FILES = (
    "CRISPRGeneEffect.csv",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "Model.csv",
)
GENE_COL_RE = r"^(?P<symbol>.+?)\s*\((?P<entrez>\d+)\)$"
HTML_SNIFF = (b"<!DOCTYPE", b"<!doctype", b"<html", b"<HTML")
DEFAULT_LINEAGE_COL = "OncotreeLineage"
DEPENDENCY_THRESHOLD = -0.5


def release_sort_key(tag):
    """'24Q4' -> (24, 4) so releases sort chronologically."""
    m = re.match(r"^(\d\d)Q(\d)$", tag)
    assert m is not None, "release tag %r is not of the form YYQN" % (tag,)
    return (int(m.group(1)), int(m.group(2)))


def resolve_depmap_release(release=None, cache_path=None, session=None,
                           timeout=60):
    """Resolve a DepMap public release to {filename: download_url} via figshare.

    Parameters
    ----------
    release : str or None
        Quarter tag such as "24Q4". None picks the most recent release
        discoverable through the figshare search index.
    cache_path : str or None
        JSON file to read/write the resolved manifest. When it exists and
        matches the requested release, no network call is made.
    session : requests.Session or None
    timeout : float

    Returns
    -------
    dict
        {"release": "24Q4", "article_id": int, "doi": str,
         "files": {filename: download_url}, "sizes": {filename: int}}

    Postcondition: every name in CORE_FILES is present in ``files``.
    """
    if cache_path and os.path.exists(cache_path):
        with open(cache_path) as fh:
            cached = json.load(fh)
        if release is None or cached.get("release") == release:
            return cached

    sess = session or requests.Session()
    query = "DepMap %s Public" % release if release else "DepMap Public"
    resp = sess.post(
        FIGSHARE_API + "/articles/search",
        json={"search_for": query, "page_size": 100,
              "order": "published_date", "order_direction": "desc"},
        headers={"Accept": "application/json"}, timeout=timeout,
    )
    resp.raise_for_status()
    candidates = {}
    for art in resp.json():
        m = re.match(RELEASE_TITLE_RE, art.get("title", "").strip())
        if m:
            candidates[m.group(1)] = art
    assert candidates, (
        "figshare search for %r returned no article titled 'DepMap <YYQN> "
        "Public'; cannot resolve a release" % query)
    if release is not None:
        assert release in candidates, (
            "release %r not found on figshare; available: %s"
            % (release, sorted(candidates)))
        tag = release
    else:
        tag = max(candidates, key=release_sort_key)

    art_id = candidates[tag]["id"]
    detail = sess.get("%s/articles/%d" % (FIGSHARE_API, art_id),
                      timeout=timeout)
    detail.raise_for_status()
    detail = detail.json()
    files = {f["name"]: f["download_url"] for f in detail.get("files", [])}
    sizes = {f["name"]: f["size"] for f in detail.get("files", [])}
    missing = [n for n in CORE_FILES if n not in files]
    assert not missing, (
        "DepMap %s figshare article %d is missing expected files: %s"
        % (tag, art_id, missing))

    manifest = {"release": tag, "article_id": art_id,
                "doi": detail.get("doi"), "title": detail.get("title"),
                "published_date": detail.get("published_date"),
                "files": files, "sizes": sizes,
                "resolved_via": "figshare articles/search + articles/{id}"}
    if cache_path:
        with open(cache_path, "w") as fh:
            json.dump(manifest, fh, indent=1)
    return manifest


def fetch_depmap_file(url, dest, expected_min_bytes=1000, session=None,
                      chunk=4194304, timeout=120):
    """Download ``url`` to ``dest`` unless a large-enough copy already exists.

    Precondition (raises): the downloaded payload does not begin with an HTML
    document sniff, and is at least ``expected_min_bytes`` long. This catches
    bot-verification interstitials and error pages served with HTTP 200.
    """
    if os.path.exists(dest) and os.path.getsize(dest) >= expected_min_bytes:
        return dest
    sess = session or requests.Session()
    tmp = dest + ".part"
    os.makedirs(os.path.dirname(os.path.abspath(dest)) or ".", exist_ok=True)
    with sess.get(url, stream=True, timeout=timeout) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for block in r.iter_content(chunk_size=chunk):
                if block:
                    fh.write(block)
    with open(tmp, "rb") as fh:
        head = fh.read(512)
    assert not any(head.lstrip().startswith(s) for s in HTML_SNIFF), (
        "download from %s is an HTML page, not data (first bytes: %r). "
        "depmap.org serves a bot-verification interstitial; use figshare."
        % (url, head[:120]))
    size = os.path.getsize(tmp)
    assert size >= expected_min_bytes, (
        "download from %s is %d bytes, below expected_min_bytes=%d"
        % (url, size, expected_min_bytes))
    os.replace(tmp, dest)
    return dest


def parse_gene_columns(header):
    """Map DepMap 'SYMBOL (ENTREZ)' column names to bare symbols.

    Returns {column_name: symbol} for every column that matches the pattern.
    Columns that do not match (e.g. the row-id column) are omitted.
    """
    out = {}
    for col in header:
        m = re.match(GENE_COL_RE, str(col))
        if m:
            out[col] = m.group("symbol").strip()
    return out


def read_header(path):
    return pd.read_csv(path, nrows=0).columns.tolist()


def load_wide(path, genes, index_name):
    """Shared loader for the wide (models x genes) DepMap matrices.

    Memory discipline: these files are 400-500 MB with ~19k gene columns. When
    ``genes`` is given only those columns are parsed via ``usecols``, so peak
    RSS stays in the tens of MB instead of several GB.
    """
    header = _read_header(path)
    id_col = header[0]
    sym_by_col = parse_gene_columns(header)
    if genes is None:
        df = pd.read_csv(path, index_col=0, low_memory=False)
        df.columns = [sym_by_col.get(c, c) for c in df.columns]
    else:
        genes = list(dict.fromkeys(genes))
        col_by_sym = {}
        for col, sym in sym_by_col.items():
            col_by_sym.setdefault(sym, col)
        missing = [g for g in genes if g not in col_by_sym]
        assert not missing, (
            "requested genes absent from %s header: %s"
            % (os.path.basename(path), missing))
        cols = [col_by_sym[g] for g in genes]
        df = pd.read_csv(path, index_col=0, usecols=[id_col] + cols,
                         low_memory=False)
        df = df[cols]
        df.columns = [sym_by_col[c] for c in cols]
    df.index.name = index_name
    return df


def load_gene_effect(path, genes=None):
    """Load CRISPRGeneEffect.csv (Chronos). Rows = ModelID, cols = symbols."""
    return _load_wide(path, genes, "ModelID")


def load_expression(path, genes=None):
    """Load OmicsExpressionProteinCodingGenesTPMLogp1.csv, log2(TPM+1)."""
    return _load_wide(path, genes, "ModelID")


def load_model_annotation(path):
    """Load Model.csv indexed by ModelID."""
    df = pd.read_csv(path, low_memory=False)
    assert "ModelID" in df.columns, "Model.csv has no ModelID column"
    return df.set_index("ModelID")


def lineage_dependency(gene_effect, model_df, genes=None,
                       lineage_col=None, min_lines=1):
    """Per-lineage mean CRISPR gene effect.

    Parameters
    ----------
    gene_effect : DataFrame  (ModelID x gene symbol)
    model_df : DataFrame indexed by ModelID, carrying ``lineage_col``
    genes : iterable or None -- subset of columns; None uses all
    min_lines : int -- drop lineage/gene cells scored on fewer lines than this

    Returns
    -------
    long DataFrame [gene, lineage, mean_gene_effect, median_gene_effect,
                    n_lines, n_dependent_lines]
    where n_dependent_lines counts lines with gene effect < -0.5.

    Precondition (raises): every cell line kept has a non-null lineage.
    """
    if lineage_col is None:
        lineage_col = DEFAULT_LINEAGE_COL
    assert lineage_col in model_df.columns, (
        "model annotation has no column %r; available: %s"
        % (lineage_col, list(model_df.columns)[:25]))
    cols = list(gene_effect.columns) if genes is None else [
        g for g in dict.fromkeys(genes) if g in gene_effect.columns]
    assert cols, "none of the requested genes are columns of gene_effect"
    lin = model_df[lineage_col]
    shared = gene_effect.index.intersection(lin.dropna().index)
    assert len(shared) > 0, "no ModelID overlap between gene_effect and model_df"
    ge = gene_effect.loc[shared, cols]
    lab = lin.loc[shared]
    assert lab.notna().all(), "internal: unannotated lineage survived the join"

    rows = []
    for lineage, idx in lab.groupby(lab).groups.items():
        block = ge.loc[idx]
        for g in cols:
            v = block[g].dropna()
            if len(v) < min_lines:
                continue
            rows.append({
                "gene": g, "lineage": lineage,
                "mean_gene_effect": float(v.mean()),
                "median_gene_effect": float(v.median()),
                "n_lines": int(len(v)),
                "n_dependent_lines": int((v < DEPENDENCY_THRESHOLD).sum()),
            })
    out = pd.DataFrame(rows)
    return out.sort_values(["gene", "mean_gene_effect"]).reset_index(drop=True)


def bh_fdr(pvals):
    """Benjamini-Hochberg adjusted p-values. NaN-safe, order preserving."""
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = int(ok.sum())
    if n == 0:
        return out
    idx = np.where(ok)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order] * n / np.arange(1, n + 1)
    out[order] = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    return out


def correlate_with_external(per_lineage_effect, external_score, lineage_map,
                            min_lineages=5, effect_col="mean_gene_effect"):
    """Spearman-correlate per-lineage gene effect against an external score.

    Parameters
    ----------
    per_lineage_effect : long DataFrame from :func:`lineage_dependency`
    external_score : Series or dict keyed by EXTERNAL type label (e.g. TCGA
        tumour-type abbreviation) -> score
    lineage_map : dict or DataFrame
        Explicit DepMap-lineage -> external-type mapping. REQUIRED: the mapping
        is a judgement call that must be recorded and disclosed, so this
        function never infers it. A DataFrame must have columns
        ``depmap_lineage`` and ``tcga_type``.
    min_lineages : int -- genes with fewer paired lineages are skipped

    Returns
    -------
    DataFrame [gene, spearman_rho, p_value, n_lineages, fdr]
    """
    if isinstance(lineage_map, pd.DataFrame):
        assert {"depmap_lineage", "tcga_type"} <= set(lineage_map.columns), (
            "lineage_map DataFrame needs columns depmap_lineage, tcga_type")
        lm = lineage_map.dropna(subset=["tcga_type"])
        lineage_map = dict(zip(lm["depmap_lineage"], lm["tcga_type"]))
    assert isinstance(lineage_map, dict) and lineage_map, (
        "lineage_map must be a non-empty dict or DataFrame -- pass it "
        "explicitly; it is never inferred")
    ext = pd.Series(external_score, dtype=float).dropna()

    df = per_lineage_effect.copy()
    df["tcga_type"] = df["lineage"].map(lineage_map)
    df = df.dropna(subset=["tcga_type"])
    df["external_score"] = df["tcga_type"].map(ext)
    df = df.dropna(subset=["external_score", effect_col])

    rows = []
    for gene, block in df.groupby("gene"):
        block = block.groupby("tcga_type").agg(
            {effect_col: "mean", "external_score": "first"})
        n = len(block)
        if n < min_lineages:
            rows.append({"gene": gene, "spearman_rho": np.nan,
                         "p_value": np.nan, "n_lineages": n})
            continue
        rho, p = stats.spearmanr(block[effect_col], block["external_score"])
        rows.append({"gene": gene, "spearman_rho": float(rho),
                     "p_value": float(p), "n_lineages": n})
    out = pd.DataFrame(rows)
    out["fdr"] = bh_fdr(out["p_value"].values)
    return out.sort_values("p_value", na_position="last").reset_index(drop=True)


def select_panel(expr, gene_effect, model_df, receptor, cofactors=(),
                 lineage=None, high_pct=85, low_pct=15, activity_genes=(),
                 exclude_subtypes=(), lineage_col=None,
                 subtype_col="OncotreeSubtype", n_per_arm=4):
    """Pick high- and low-receptor cell lines for a perturbation experiment.

    The candidate pool is restricted to ``lineage`` (when given), minus any line
    whose ``subtype_col`` value is in ``exclude_subtypes``. Receptor percentiles
    are computed WITHIN that pool, so they describe the panel rather than the
    whole of DepMap. ``activity_genes`` gives a baseline pathway-activity
    readout (mean log2(TPM+1) of downstream targets): the responder arm should
    be low here, or an induction assay has no dynamic range.

    Returns
    -------
    DataFrame [arm, cell_line, ModelID, lineage, subtype, receptor,
               receptor_expr, receptor_panel_pct, <cofactor cols>,
               baseline_activity, receptor_gene_effect]

    Precondition (raises): every selected line has a non-null lineage.
    """
    if lineage_col is None:
        lineage_col = DEFAULT_LINEAGE_COL
    assert receptor in expr.columns, (
        "receptor %r not in expression matrix columns" % receptor)
    ann = model_df
    pool = expr.index.intersection(ann.index)
    if lineage is not None:
        want = {lineage} if isinstance(lineage, str) else set(lineage)
        pool = pool[ann.loc[pool, lineage_col].isin(want)]
    if len(exclude_subtypes) and subtype_col in ann.columns:
        pool = pool[~ann.loc[pool, subtype_col].isin(set(exclude_subtypes))]
    pool = pool[expr.loc[pool, receptor].notna()]
    assert len(pool) >= 2 * n_per_arm, (
        "candidate pool has %d lines, too few for %d per arm"
        % (len(pool), n_per_arm))

    rec = expr.loc[pool, receptor]
    pct = rec.rank(pct=True) * 100
    hi = rec[pct >= high_pct].sort_values(ascending=False).index[:n_per_arm]
    lo = rec[pct <= low_pct].sort_values().index[:n_per_arm]
    assert len(hi) and len(lo), (
        "high_pct=%s / low_pct=%s produced an empty arm on %d lines"
        % (high_pct, low_pct, len(pool)))

    act_use = [g for g in activity_genes if g in expr.columns]
    rows = []
    for arm, ids in (("responder_high_receptor", hi),
                     ("control_low_receptor", lo)):
        for mid in ids:
            r = {
                "arm": arm,
                "cell_line": ann.loc[mid].get("StrippedCellLineName",
                                              ann.loc[mid].get("CellLineName")),
                "ModelID": mid,
                "lineage": ann.loc[mid, lineage_col],
                "subtype": ann.loc[mid].get(subtype_col),
                "receptor": receptor,
                "receptor_expr_log2tpm1": float(expr.loc[mid, receptor]),
                "receptor_panel_pct": float(pct.loc[mid]),
            }
            for cf in cofactors:
                r["cofactor_%s_log2tpm1" % cf] = (
                    float(expr.loc[mid, cf]) if cf in expr.columns else np.nan)
            r["baseline_activity_log2tpm1"] = (
                float(expr.loc[mid, act_use].mean()) if act_use else np.nan)
            r["n_activity_genes"] = len(act_use)
            r["receptor_gene_effect"] = (
                float(gene_effect.loc[mid, receptor])
                if (receptor in gene_effect.columns
                    and mid in gene_effect.index)
                else np.nan)
            rows.append(r)
    out = pd.DataFrame(rows)
    assert out["lineage"].notna().all(), (
        "selected lines lack a lineage annotation: %s"
        % out.loc[out["lineage"].isna(), "ModelID"].tolist())
    return out
