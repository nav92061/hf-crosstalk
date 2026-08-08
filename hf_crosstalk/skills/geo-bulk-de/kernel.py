"""geo-bulk-de: fetch a GEO series, parse its sample metadata, and run
two-group differential expression entirely in Python (no R / no limma call).

Pure stdlib + pandas/numpy/scipy/statsmodels. No platform-specific imports.
Every public function takes file paths / DataFrames and returns DataFrames or
dicts, so the module is usable from a plain Python session with no agent.
"""

import gzip
import io
import os
import re
import tarfile
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import digamma, polygamma
from statsmodels.nonparametric.smoothers_lowess import lowess

GEO_FTP_ROOT = "https://ftp.ncbi.nlm.nih.gov/geo"
GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series"
DEFAULT_TIMEOUT = 300
CPM_PRIOR_COUNT = 0.5
LOG2_PRIOR_COUNT = 1.0


# ---------------------------------------------------------------------------
# GEO download helpers
# ---------------------------------------------------------------------------

def gse_ftp_dir(gse):
    """Return the GEO FTP/HTTPS directory URL for a series.

    Args:
        gse (str): series accession, e.g. "GSE116250" (case-insensitive).

    Returns:
        str: e.g. "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE116nnn/GSE116250"
    """
    gse = str(gse).strip().upper()
    m = re.fullmatch(r"GSE(\d+)", gse)
    if m is None:
        raise ValueError("not a GSE accession: %r" % (gse,))
    digits = m.group(1)
    stub = "GSE%snnn" % (digits[:-3] if len(digits) > 3 else "",)
    return "%s/%s/%s" % (GEO_BASE, stub, gse)


def h_http_get(url, timeout=None):
    if timeout is None:
        timeout = None
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        return fh.read()


def h_download(url, dest_path, timeout=None, overwrite=False):
    if timeout is None:
        timeout = None
    if os.path.exists(dest_path) and not overwrite and os.path.getsize(dest_path) > 0:
        return dest_path
    tmp = dest_path + ".part"
    with urllib.request.urlopen(url, timeout=timeout) as src, open(tmp, "wb") as dst:
        while True:
            chunk = src.read(1 << 20)
            if not chunk:
                break
            dst.write(chunk)
    os.replace(tmp, dest_path)
    return dest_path


def list_geo_supplementary(gse, timeout=None):
    """List supplementary filenames available for a GEO series.

    Args:
        gse (str): series accession.
        timeout (int): socket timeout in seconds.

    Returns:
        list[str]: filenames (not URLs), directory order.
    """
    if timeout is None:
        timeout = None
    url = gse_ftp_dir(gse) + "/suppl/"
    html = h_http_get(url, timeout=timeout).decode("utf-8", "replace")
    names = re.findall(r'href="([^"?/][^"]*)"', html)
    seen, out = set(), []
    for n in names:
        if n in seen or n.startswith(("http", "..")):
            continue
        seen.add(n)
        out.append(n)
    return out


def fetch_geo_supplementary(gse, dest_dir, pattern=None, timeout=None,
                            overwrite=False):
    """Download supplementary files of a GEO series.

    Args:
        gse (str): series accession, e.g. "GSE116250".
        dest_dir (str): directory to write into (created if absent).
        pattern (str | None): regex; if given, only filenames matching it
            (re.search, case-insensitive) are downloaded. None = all files.
        timeout (int): socket timeout in seconds.
        overwrite (bool): re-download even if a non-empty file already exists.

    Returns:
        list[str]: local paths of downloaded files, in download order.

    Raises:
        FileNotFoundError: if no supplementary file matches `pattern`.
    """
    if timeout is None:
        timeout = DEFAULT_TIMEOUT
    if timeout is None:
        timeout = DEFAULT_TIMEOUT
    os.makedirs(dest_dir, exist_ok=True)
    base = gse_ftp_dir(gse) + "/suppl/"
    names = list_geo_supplementary(gse, timeout=timeout)
    if pattern is not None:
        rx = re.compile(pattern, re.I)
        names = [n for n in names if rx.search(n)]
    if not names:
        raise FileNotFoundError(
            "no supplementary files for %s matching %r" % (gse, pattern))
    return [h_download(base + n, os.path.join(dest_dir, n), timeout=timeout,
                      overwrite=overwrite) for n in names]


def fetch_series_matrix(gse, dest_dir, timeout=None, overwrite=False,
                        platform=None):
    """Download the series matrix file (sample metadata) for a GEO series.

    Args:
        gse (str): series accession.
        dest_dir (str): directory to write into (created if absent).
        timeout (int): socket timeout in seconds.
        overwrite (bool): re-download even if present.
        platform (str | None): for a multi-platform series, the GPL accession
            whose matrix to fetch, giving `<GSE>-<GPL>_series_matrix.txt.gz`.
            None (default) fetches the single-platform
            `<GSE>_series_matrix.txt.gz` - which does not exist for a
            multi-platform submission; call `list_series_matrix_files` when
            unsure.

    Returns:
        str: local path to the downloaded series matrix.
    """
    if timeout is None:
        timeout = DEFAULT_TIMEOUT
    os.makedirs(dest_dir, exist_ok=True)
    gse = str(gse).strip().upper()
    if platform is None:
        name = "%s_series_matrix.txt.gz" % gse
    else:
        name = "%s-%s_series_matrix.txt.gz" % (gse, str(platform).strip().upper())
    url = "%s/matrix/%s" % (gse_ftp_dir(gse), name)
    return h_download(url, os.path.join(dest_dir, name), timeout=timeout,
                     overwrite=overwrite)


def extract_tar(tar_path, dest_dir, pattern=None):
    """Extract (optionally filtered) members of a tar archive.

    Args:
        tar_path (str): path to a .tar / .tar.gz archive.
        dest_dir (str): destination directory (created if absent).
        pattern (str | None): regex filter on member names (re.search, case-insens.).

    Returns:
        list[str]: local paths of extracted regular files, sorted.
    """
    os.makedirs(dest_dir, exist_ok=True)
    rx = re.compile(pattern, re.I) if pattern else None
    out = []
    with tarfile.open(tar_path) as tf:
        for m in tf.getmembers():
            if not m.isfile():
                continue
            name = os.path.basename(m.name)
            if not name or name.startswith("."):
                continue
            if rx is not None and not rx.search(name):
                continue
            target = os.path.join(dest_dir, name)
            if not os.path.abspath(target).startswith(os.path.abspath(dest_dir)):
                raise ValueError("unsafe tar member path: %r" % (m.name,))
            if not os.path.exists(target) or os.path.getsize(target) == 0:
                src = tf.extractfile(m)
                with open(target, "wb") as dst:
                    dst.write(src.read())
            out.append(target)
    return sorted(out)


# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

def h_open_text(path):
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def h_split_meta_line(line):
    parts = line.rstrip("\n").rstrip("\r").split("\t")
    key = parts[0].lstrip("!")
    vals = [p.strip().strip('"') for p in parts[1:]]
    return key, vals


def parse_series_matrix_metadata(path):
    """Parse per-sample metadata out of a GEO series_matrix.txt.gz.

    Every `!Sample_characteristics_ch1` row is split on the first ":" into a
    key/value pair; the key becomes a column (deduplicated with `_2`, `_3`
    suffixes when the same key appears on several rows). `!Sample_title`,
    `!Sample_source_name_ch1`, `!Sample_description` and
    `!Sample_geo_accession` are also returned. `description` matters for older
    series (e.g. GSE1869) that carry the disease group there and have no
    `!Sample_characteristics_ch1` line at all - check it before concluding a
    series has no usable labels.

    Args:
        path (str): path to the (optionally gzipped) series matrix file.

    Returns:
        pandas.DataFrame: one row per sample, indexed by GSM accession.
            Always contains a 'gsm' column; other columns depend on the series.

    Raises:
        ValueError: if no `!Sample_geo_accession` line is present, or if the
            characteristics rows are ragged relative to the sample count.
    """
    gsms = None
    simple = {}
    char_rows = []
    with h_open_text(path) as fh:
        for line in fh:
            if line.startswith("!Sample_"):
                key, vals = h_split_meta_line(line)
                if key == "Sample_geo_accession":
                    gsms = vals
                elif key == "Sample_characteristics_ch1":
                    char_rows.append(vals)
                elif key in ("Sample_title", "Sample_source_name_ch1",
                             "Sample_organism_ch1", "Sample_library_strategy",
                             "Sample_description"):
                    simple[key.replace("Sample_", "").lower()] = vals
            elif line.startswith("!series_matrix_table_begin"):
                break
    if not gsms:
        raise ValueError("no !Sample_geo_accession line found in %s" % path)
    n = len(gsms)
    df = pd.DataFrame({"gsm": gsms})
    for k, v in simple.items():
        if len(v) != n:
            raise ValueError("metadata row %r has %d values for %d samples"
                             % (k, len(v), n))
        df[k] = v
    used = {}
    for vals in char_rows:
        if len(vals) != n:
            raise ValueError("ragged !Sample_characteristics_ch1 row: %d values "
                             "for %d samples" % (len(vals), n))
        keys = [v.split(":", 1)[0].strip().lower() if ":" in v else "characteristic"
                for v in vals]
        key = pd.Series(keys).mode().iat[0] if keys else "characteristic"
        key = re.sub(r"[^0-9a-z]+", "_", key).strip("_") or "characteristic"
        used[key] = used.get(key, 0) + 1
        col = key if used[key] == 1 else "%s_%d" % (key, used[key])
        df[col] = [v.split(":", 1)[1].strip() if ":" in v else v.strip() for v in vals]
    df = df.set_index("gsm", drop=False)
    if df.index.duplicated().any():
        raise ValueError("duplicate GSM accessions in series matrix")
    return df


# ---------------------------------------------------------------------------
# Expression loading
# ---------------------------------------------------------------------------

def load_expression(path, sep="\t", gene_col=0, drop_cols=None, aggregate="sum"):
    """Load a genes x samples expression matrix from a delimited text file.

    Args:
        path (str): path to a (optionally .gz) text matrix.
        sep (str): field separator.
        gene_col (int | str): column holding gene identifiers.
        drop_cols (iterable[str]): extra annotation columns to drop.
        aggregate ('sum' | 'max' | 'first' | None): how to collapse duplicate
            gene identifiers. None leaves duplicates in place.

    Returns:
        pandas.DataFrame: numeric, genes (index) x samples (columns).

    Raises:
        ValueError: if the matrix is empty, non-numeric after coercion, or a
            gene row is entirely NA.
    """
    if drop_cols is None:
        drop_cols = ()
    df = pd.read_csv(path, sep=sep, low_memory=False)
    gcol = df.columns[gene_col] if isinstance(gene_col, int) else gene_col
    df = df.set_index(gcol)
    drop = [c for c in drop_cols if c in df.columns]
    if drop:
        df = df.drop(columns=drop)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.loc[:, ~df.isna().all(axis=0)]
    df.index = df.index.astype(str)
    if aggregate is not None and df.index.duplicated().any():
        df = getattr(df.groupby(level=0), aggregate)() if aggregate != "first" \
            else df.groupby(level=0).first()
    h_check_matrix(df)
    return df


def load_per_sample_tables(paths, value_col=None, gene_col=0, sep=",",
                           sample_name_regex=r"(GSM\d+)"):
    """Assemble a genes x samples matrix from one file per sample.

    Args:
        paths (iterable[str]): per-sample table paths (optionally .gz).
        value_col (str | int | None): column holding the values. None = the
            first numeric column after the gene column.
        gene_col (int | str): gene identifier column.
        sep (str): field separator.
        sample_name_regex (str): regex applied to the basename; group 1 is the
            sample name. If it does not match, the basename (minus extensions)
            is used.

    Returns:
        pandas.DataFrame: genes x samples, inner-joined on the gene index.

    Raises:
        ValueError: on duplicate sample names, or if the per-file gene indexes
            share no genes.
    """
    rx = re.compile(sample_name_regex)
    cols = {}
    for p in paths:
        base = os.path.basename(p)
        m = rx.search(base)
        name = m.group(1) if m else re.sub(r"\.(csv|tsv|txt)(\.gz)?$", "", base)
        if name in cols:
            raise ValueError("duplicate sample name %r from %s" % (name, p))
        d = pd.read_csv(p, sep=sep, low_memory=False)
        gcol = d.columns[gene_col] if isinstance(gene_col, int) else gene_col
        if value_col is None:
            cand = [c for c in d.columns if c != gcol
                    and pd.api.types.is_numeric_dtype(d[c])]
            if not cand:
                raise ValueError("no numeric value column in %s" % p)
            vcol = cand[0]
        else:
            vcol = d.columns[value_col] if isinstance(value_col, int) else value_col
        s = pd.to_numeric(d.set_index(d[gcol].astype(str))[vcol], errors="coerce")
        s = s.groupby(level=0).sum()
        cols[name] = s
    mat = pd.DataFrame(cols)
    mat = mat.dropna(how="any")
    if mat.empty:
        raise ValueError("per-sample tables share no common genes")
    h_check_matrix(mat)
    return mat


def strip_ensembl_version(index):
    """Strip trailing `.N` version suffixes from Ensembl-style identifiers.

    Args:
        index (iterable[str]): identifiers.

    Returns:
        list[str]: identifiers with `ENSG00000123456.7` -> `ENSG00000123456`.
    """
    return [re.sub(r"^(ENS[A-Z]*[GTP]\d+)\.\d+$", r"\1", str(i)) for i in index]


# ---------------------------------------------------------------------------
# Preconditions
# ---------------------------------------------------------------------------

def h_check_matrix(expr):
    if expr is None or expr.shape[0] == 0 or expr.shape[1] == 0:
        raise ValueError("expression matrix is empty (shape=%r)"
                         % (None if expr is None else expr.shape,))
    if pd.Index(expr.columns).duplicated().any():
        raise ValueError("sample identifiers (columns) are not unique")
    allna = expr.isna().all(axis=1)
    if bool(allna.any()):
        raise ValueError("%d gene rows are entirely NA (e.g. %s)"
                         % (int(allna.sum()), list(expr.index[allna][:5])))
    nonnum = [c for c in expr.columns if not pd.api.types.is_numeric_dtype(expr[c])]
    if nonnum:
        raise ValueError("non-numeric sample columns: %s" % nonnum[:5])


def check_de_inputs(expr, groups, case, control, min_per_group=2):
    """Validate a DE call's inputs; raise on any violated precondition.

    Args:
        expr (pandas.DataFrame): genes x samples.
        groups (pandas.Series | dict | list): group label per sample. A Series
            is aligned by index to `expr.columns`; a list must be positional.
        case (str): group label treated as the numerator of log2FC.
        control (str): group label treated as the denominator.
        min_per_group (int): minimum samples required in each group.

    Returns:
        pandas.Series: group labels reindexed to `expr.columns`.

    Raises:
        ValueError: empty/NA matrix, duplicate sample ids, groups not aligned to
            the matrix columns, missing case/control label, or fewer than
            `min_per_group` samples in either group.
    """
    h_check_matrix(expr)
    if isinstance(groups, pd.Series):
        g = groups.copy()
    elif isinstance(groups, dict):
        g = pd.Series(groups)
    else:
        g = pd.Series(list(groups), index=list(expr.columns))
    g.index = g.index.astype(str)
    cols = pd.Index([str(c) for c in expr.columns])
    missing = [c for c in cols if c not in g.index]
    if missing:
        raise ValueError("group vector is not aligned to the expression matrix: "
                         "%d columns have no label (e.g. %s)"
                         % (len(missing), missing[:5]))
    g = g.reindex(cols)
    for lab in (case, control):
        if lab not in set(g.dropna()):
            raise ValueError("group label %r absent from the group vector "
                             "(present: %s)" % (lab, sorted(set(g.dropna()))))
    n_case = int((g == case).sum())
    n_ctrl = int((g == control).sum())
    if n_case < min_per_group or n_ctrl < min_per_group:
        raise ValueError("need >=%d samples per group; got case=%d control=%d"
                         % (min_per_group, n_case, n_ctrl))
    return g


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def filter_by_expression(expr, min_value=1.0, min_fraction=0.2, groups=None):
    """Drop low-expression genes.

    A gene is kept when `min_value` is exceeded in at least `min_fraction` of
    samples (or, if `groups` is given, in at least that fraction of the smaller
    group's samples).

    Args:
        expr (pandas.DataFrame): genes x samples, on the scale you want to
            threshold (CPM for counts, RPKM/FPKM for normalised data).
        min_value (float): threshold value.
        min_fraction (float): fraction of samples that must exceed it.
        groups (pandas.Series | None): optional group labels; the required
            sample count is based on the smallest group when supplied.

    Returns:
        pandas.DataFrame: filtered matrix (same column order).
    """
    n = expr.shape[1]
    if groups is not None:
        sizes = pd.Series(groups).value_counts()
        n = int(sizes.min()) if len(sizes) else n
    need = max(2, int(np.ceil(min_fraction * n)))
    keep = (expr > min_value).sum(axis=1) >= need
    return expr.loc[keep]


def cpm(counts, log=False, prior_count=None, lib_sizes=None):
    """Counts per million, optionally log2.

    log-CPM uses the edgeR-style offset: log2((count + prior)/(lib + 2*prior) * 1e6).

    Args:
        counts (pandas.DataFrame): genes x samples raw counts.
        log (bool): return log2-CPM instead of CPM.
        prior_count (float): pseudo-count added before logging.
        lib_sizes (pandas.Series | None): per-sample library sizes; default is
            the column sums of `counts`.

    Returns:
        pandas.DataFrame: same shape as `counts`.
    """
    if prior_count is None:
        prior_count = CPM_PRIOR_COUNT
    lib = counts.sum(axis=0) if lib_sizes is None else pd.Series(lib_sizes)
    lib = lib.reindex(counts.columns).astype(float)
    if (lib <= 0).any():
        raise ValueError("library size <= 0 for samples: %s"
                         % list(lib.index[lib <= 0])[:5])
    if not log:
        return counts.divide(lib, axis=1) * 1e6
    return np.log2((counts + prior_count).divide(lib + 2 * prior_count, axis=1) * 1e6)


# ---------------------------------------------------------------------------
# limma-equivalent machinery (implemented here, not a call into limma)
# ---------------------------------------------------------------------------

def trigamma_inverse(x):
    """Invert the trigamma function elementwise (Newton iteration).

    Args:
        x (array-like): positive values.

    Returns:
        numpy.ndarray: y such that trigamma(y) == x.
    """
    x = np.asarray(x, dtype=float)
    y = np.full_like(x, np.nan)
    ok = np.isfinite(x) & (x > 0)
    xs = x[ok]
    out = np.where(xs > 1e7, 1.0 / np.sqrt(xs), np.exp(xs) + 0.5)
    out = np.where(xs < 1e-6, 1.0 / xs, out)
    for _ in range(60):
        tri = polygamma(1, out)
        dif = tri * (1 - tri / xs) / polygamma(2, out)
        out = out + dif
        if np.max(np.abs(dif / np.maximum(out, 1e-12))) < 1e-8:
            break
    y[ok] = out
    return y


def fit_f_dist(s2, df_resid):
    """Smyth (2004) moment estimator of the scaled-F prior on residual variances.

    Matches the first two moments of log(s2) to those of a scaled
    chi-square/F, giving the prior degrees of freedom d0 and prior variance s0^2.

        e      = log(s2) - digamma(df/2) + log(df/2)
        evar   = var(e) - mean(trigamma(df/2))
        d0     = 2 * trigamma^-1(evar)                      (evar > 0)
        s0^2   = exp(mean(e) + digamma(d0/2) - log(d0/2))
        evar<=0 -> d0 = inf, s0^2 = exp(mean(e))            (full shrinkage)

    Args:
        s2 (array-like): per-gene residual variances (> 0).
        df_resid (float | array-like): residual degrees of freedom per gene.

    Returns:
        dict: {'d0': float, 's0_2': float} (d0 may be numpy.inf).
    """
    s2 = np.asarray(s2, dtype=float)
    dfr = np.broadcast_to(np.asarray(df_resid, dtype=float), s2.shape)
    ok = np.isfinite(s2) & (s2 > 0) & (dfr > 0)
    if ok.sum() < 2:
        raise ValueError("fewer than 2 usable residual variances for fit_f_dist")
    z = np.log(s2[ok])
    d = dfr[ok]
    e = z - digamma(d / 2.0) + np.log(d / 2.0)
    ebar = float(np.mean(e))
    evar = float(np.mean((e - ebar) ** 2) * (ok.sum() / (ok.sum() - 1.0))
                 - np.mean(polygamma(1, d / 2.0)))
    if evar > 0:
        d0 = float(2.0 * trigamma_inverse(evar))
        s0_2 = float(np.exp(ebar + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    else:
        d0 = np.inf
        s0_2 = float(np.exp(ebar))
    return {"d0": d0, "s0_2": s0_2}


def h_variance_trend(mean_expr, s2, frac=0.4):
    """lowess trend of residual variance on mean expression (returns fitted s2)."""
    y = np.log(np.maximum(s2, 1e-12))
    ok = np.isfinite(y) & np.isfinite(mean_expr)
    fit = lowess(y[ok], mean_expr[ok], frac=frac, return_sorted=True)
    xs, ys = fit[:, 0], fit[:, 1]
    trend = np.interp(mean_expr, xs, ys, left=ys[0], right=ys[-1])
    return np.exp(trend)


def h_lstsq_fit(Y, X, weights=None):
    """Least-squares (optionally weighted) fit of every gene row of Y on X."""
    G, n = Y.shape
    p = X.shape[1]
    if weights is None:
        W = np.ones((G, n))
    else:
        W = np.asarray(weights, dtype=float)
    XtWX = np.einsum("ni,gn,nj->gij", X, W, X)
    XtWy = np.einsum("ni,gn,gn->gi", X, W, Y)
    beta = np.linalg.solve(XtWX, XtWy[:, :, None])[:, :, 0]
    fitted = beta @ X.T
    resid = Y - fitted
    sse = np.einsum("gn,gn->g", W, resid ** 2)
    df_resid = float(n - p)
    s2 = sse / df_resid
    cov = np.linalg.inv(XtWX)
    return beta, s2, df_resid, cov, fitted


def h_voom_weights(counts, X, lib_sizes, frac=0.5):
    """voom-style precision weights from a counts matrix and a design."""
    lib = np.asarray(lib_sizes, dtype=float)
    y = np.log2((counts + 0.5) / (lib + 1.0) * 1e6)
    beta, s2, dfr, _, fitted = h_lstsq_fit(y, X)
    sd = np.sqrt(np.maximum(s2, 1e-12))
    mean_log_count = np.log2(counts + 0.5).mean(axis=1)
    ok = np.isfinite(sd) & (sd > 0)
    fit = lowess(np.sqrt(sd[ok]), mean_log_count[ok], frac=frac, return_sorted=True)
    xs, ys = fit[:, 0], fit[:, 1]
    # predicted per-observation log2 count from the fitted log-CPM values
    fitted_log_count = fitted + np.log2(lib + 1.0)[None, :] - np.log2(1e6)
    sqrt_sd = np.interp(fitted_log_count, xs, ys, left=ys[0], right=ys[-1])
    w = 1.0 / np.maximum(sqrt_sd, 1e-8) ** 4
    return w


def moderated_t(beta, s2, df_resid, stdev_unscaled, mean_expr, trend=True,
                trend_frac=0.4):
    """Empirical-Bayes moderated t statistics (limma eBayes equivalent).

    Posterior variance:  s~^2 = (d0*s0^2 + df*s^2) / (d0 + df)
    Moderated t:         t = beta / (s~ * stdev_unscaled), df = df + d0
    With `trend=True`, s0^2 is gene-specific: s0^2_g = scale * trend(mean_expr_g),
    where trend() is a lowess fit of log(s^2) on mean expression and `scale`
    comes from `fit_f_dist` applied to s^2 / trend (limma-trend).

    Args:
        beta (array): per-gene contrast estimate (log2 fold change).
        s2 (array): per-gene residual variance.
        df_resid (float): residual degrees of freedom.
        stdev_unscaled (float | array): sqrt of the contrast's variance factor
            from the design, i.e. sqrt(c' (X'X)^-1 c).
        mean_expr (array): per-gene average expression (trend covariate).
        trend (bool): fit a mean-variance trend for the prior (limma-trend) or
            use a single global prior variance.
        trend_frac (float): lowess span for the trend.

    Returns:
        dict: {'t', 'p_value', 's2_post', 'df_total', 'd0', 's0_2'} where 's0_2'
            is an array when `trend=True`, else a float.
    """
    beta = np.asarray(beta, dtype=float)
    s2 = np.asarray(s2, dtype=float)
    mean_expr = np.asarray(mean_expr, dtype=float)
    if trend:
        tr = h_variance_trend(mean_expr, s2, frac=trend_frac)
        fit = fit_f_dist(s2 / tr, df_resid)
        d0, s0_2 = fit["d0"], fit["s0_2"] * tr
    else:
        fit = fit_f_dist(s2, df_resid)
        d0, s0_2 = fit["d0"], np.full_like(s2, fit["s0_2"])
    if np.isinf(d0):
        s2_post = np.asarray(s0_2, dtype=float)
        df_total = np.full_like(s2, np.inf)
    else:
        s2_post = (d0 * s0_2 + df_resid * s2) / (d0 + df_resid)
        df_total = np.full_like(s2, df_resid + d0)
    se = np.sqrt(s2_post) * np.asarray(stdev_unscaled, dtype=float)
    t = beta / np.where(se > 0, se, np.nan)
    with np.errstate(invalid="ignore"):
        p = np.where(np.isinf(df_total),
                     2.0 * stats.norm.sf(np.abs(t)),
                     2.0 * stats.t.sf(np.abs(t), df_total))
    return {"t": t, "p_value": p, "s2_post": s2_post, "df_total": df_total,
            "d0": d0, "s0_2": s0_2}


def bh_fdr(p):
    """Benjamini-Hochberg FDR.

    Args:
        p (array-like): p-values (NaNs propagate as NaN).

    Returns:
        numpy.ndarray: adjusted p-values, same order as input.
    """
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    n = pv.size
    if n == 0:
        return out
    order = np.argsort(pv)
    ranked = pv[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty(n)
    adj[order] = np.minimum(ranked, 1.0)
    out[ok] = adj
    return out


# ---------------------------------------------------------------------------
# Main DE entry point
# ---------------------------------------------------------------------------

def run_de(expr, groups, case, control, design="limma_trend", covariates=None,
           min_per_group=2, filter_min_value=None, filter_min_fraction=0.2,
           trend_frac=0.4, already_logged=False):
    """Two-group differential expression with empirical-Bayes moderated t.

    Args:
        expr (pandas.DataFrame): genes x samples. Raw counts for
            design='voom_like' or 'limma_trend_counts'; RPKM/FPKM/TPM (or
            already-logged values, see `already_logged`) for 'limma_trend'.
        groups (pandas.Series | dict | list): group label per sample; a Series
            is aligned by index to `expr.columns`. Samples labelled with
            anything other than `case`/`control` are dropped.
        case (str): group label used as the log2FC numerator (e.g. "failing").
        control (str): group label used as the denominator.
        design (str): one of
            'limma_trend'         - log2(x+1) of normalised values (RPKM/FPKM/TPM),
                                    moderated t with a mean-variance trend prior;
            'limma_trend_counts'  - log2-CPM from raw counts, same trend prior;
            'voom_like'           - log2-CPM with voom precision weights and a
                                    non-trended prior.
        covariates (pandas.DataFrame | None): numeric per-sample covariates
            (index = sample ids) added to the design matrix.
        min_per_group (int): precondition on group sizes.
        filter_min_value (float | None): low-expression filter threshold applied
            on the CPM (counts designs) or input (limma_trend) scale. None = no
            filtering. Recommended: 1.0.
        filter_min_fraction (float): fraction of the smaller group that must
            exceed `filter_min_value`.
        trend_frac (float): lowess span for the variance trend.
        already_logged (bool): for design='limma_trend', treat `expr` as already
            log-scale (skip the log2(x+1) transform).

    Returns:
        pandas.DataFrame: one row per tested gene, sorted by p_value, with
            columns gene, log2FC, t, p_value, fdr, mean_expr, n_case, n_ctrl,
            and attrs {'design', 'd0', 'n_genes_tested', 'df_resid'}.

    Raises:
        ValueError: on any violated precondition (see `check_de_inputs`) or an
            unknown `design`.
    """
    if design not in ("limma_trend", "limma_trend_counts", "voom_like"):
        raise ValueError("unknown design %r" % (design,))
    expr = expr.copy()
    expr.columns = [str(c) for c in expr.columns]
    g = check_de_inputs(expr, groups, case, control, min_per_group=min_per_group)
    keep_samples = [c for c in expr.columns if g[c] in (case, control)]
    expr = expr[keep_samples]
    g = g[keep_samples]
    n_case = int((g == case).sum())
    n_ctrl = int((g == control).sum())

    counts_like = design in ("limma_trend_counts", "voom_like")
    if counts_like:
        lib = expr.sum(axis=0).astype(float)
        if filter_min_value is not None:
            keep = filter_by_expression(cpm(expr, log=False), filter_min_value,
                                        filter_min_fraction, groups=g).index
            expr = expr.loc[keep]
        y = cpm(expr, log=True)
    else:
        if filter_min_value is not None:
            expr = filter_by_expression(expr, filter_min_value,
                                        filter_min_fraction, groups=g)
        y = expr if already_logged else np.log2(expr + LOG2_PRIOR_COUNT)
        lib = None
    if y.shape[0] == 0:
        raise ValueError("no genes survived the expression filter "
                         "(filter_min_value=%r)" % (filter_min_value,))

    # design matrix: intercept (control) + case indicator + covariates
    X = np.column_stack([np.ones(len(g)), (g == case).astype(float).values])
    cov_names = []
    if covariates is not None:
        cv = covariates.reindex(y.columns)
        if cv.isna().any().any():
            raise ValueError("covariates contain NA after aligning to samples")
        X = np.column_stack([X] + [cv[c].astype(float).values for c in cv.columns])
        cov_names = list(cv.columns)
    if np.linalg.matrix_rank(X) < X.shape[1]:
        raise ValueError("design matrix is rank-deficient (columns: intercept, "
                         "group, %s)" % cov_names)

    Y = y.values.astype(float)
    if design == "voom_like":
        W = h_voom_weights(expr.values.astype(float), X, lib.values)
        beta, s2, dfr, cov, _ = h_lstsq_fit(Y, X, weights=W)
        trend = False
    else:
        beta, s2, dfr, cov, _ = h_lstsq_fit(Y, X)
        trend = True
    stdev_unscaled = np.sqrt(np.maximum(cov[:, 1, 1], 0.0))
    mean_expr = Y.mean(axis=1)
    eb = moderated_t(beta[:, 1], s2, dfr, stdev_unscaled, mean_expr,
                     trend=trend, trend_frac=trend_frac)
    out = pd.DataFrame({
        "gene": list(y.index),
        "log2FC": beta[:, 1],
        "t": eb["t"],
        "p_value": eb["p_value"],
        "mean_expr": mean_expr,
        "n_case": n_case,
        "n_ctrl": n_ctrl,
    })
    out["fdr"] = bh_fdr(out["p_value"].values)
    out = out[["gene", "log2FC", "t", "p_value", "fdr", "mean_expr",
               "n_case", "n_ctrl"]]
    out = out.sort_values("p_value", kind="mergesort").reset_index(drop=True)
    out.attrs["design"] = design
    out.attrs["d0"] = eb["d0"]
    out.attrs["df_resid"] = dfr
    out.attrs["n_genes_tested"] = int(out.shape[0])
    return out


# ---------------------------------------------------------------------------
# Cross-cohort replication
# ---------------------------------------------------------------------------

def replication_stats(de_a, de_b, fdr=0.05, gene_col="gene"):
    """Compare two DE tables on their shared genes.

    Args:
        de_a (pandas.DataFrame): discovery table (needs `gene_col`, 'log2FC', 'fdr').
        de_b (pandas.DataFrame): validation table, same columns.
        fdr (float): significance threshold applied to both tables.
        gene_col (str): identifier column name.

    Returns:
        dict with keys:
            n_genes_a, n_genes_b, n_shared_genes,
            n_sig_a, n_sig_b            - significant among shared genes,
            n_replicated                - significant in both (any sign),
            n_replicated_same_sign,
            hypergeom_p                 - P(overlap >= observed | independence),
            expected_overlap,
            directional_concordance     - fraction with matching log2FC sign
                                          among genes significant in both,
            spearman_rho, spearman_p    - log2FC correlation on all shared genes,
            pearson_r_sig               - log2FC correlation among genes
                                          significant in both,
            fdr_threshold.

    Raises:
        ValueError: if a required column is missing or the tables share no genes.
    """
    for name, d in (("de_a", de_a), ("de_b", de_b)):
        for c in (gene_col, "log2FC", "fdr"):
            if c not in d.columns:
                raise ValueError("%s is missing required column %r" % (name, c))
        if d[gene_col].duplicated().any():
            raise ValueError("%s has duplicate %s values" % (name, gene_col))
    a = de_a.set_index(gene_col)
    b = de_b.set_index(gene_col)
    shared = a.index.intersection(b.index)
    if len(shared) == 0:
        raise ValueError("the two DE tables share no genes - harmonize "
                         "identifiers before calling replication_stats")
    A, B = a.loc[shared], b.loc[shared]
    sig_a = (A["fdr"] < fdr).values
    sig_b = (B["fdr"] < fdr).values
    both = sig_a & sig_b
    same_sign = np.sign(A["log2FC"].values) == np.sign(B["log2FC"].values)
    M, n, N, k = len(shared), int(sig_a.sum()), int(sig_b.sum()), int(both.sum())
    hyp = float(stats.hypergeom.sf(k - 1, M, n, N)) if k > 0 else 1.0
    rho, rho_p = stats.spearmanr(A["log2FC"].values, B["log2FC"].values)
    if both.sum() >= 3:
        pr = float(stats.pearsonr(A["log2FC"].values[both],
                                  B["log2FC"].values[both])[0])
    else:
        pr = float("nan")
    return {
        "n_genes_a": int(de_a.shape[0]),
        "n_genes_b": int(de_b.shape[0]),
        "n_shared_genes": int(M),
        "n_sig_a": n,
        "n_sig_b": N,
        "n_replicated": k,
        "n_replicated_same_sign": int((both & same_sign).sum()),
        "expected_overlap": float(n * N / M) if M else float("nan"),
        "hypergeom_p": hyp,
        "directional_concordance": (float(np.mean(same_sign[both]))
                                    if k > 0 else float("nan")),
        "spearman_rho": float(rho),
        "spearman_p": float(rho_p),
        "pearson_r_sig": pr,
        "fdr_threshold": float(fdr),
    }


def replicated_signature(de_a, de_b, fdr=0.05, gene_col="gene",
                         suffixes=None):
    """Genes significant and same-sign in both cohorts.

    Args:
        de_a (pandas.DataFrame): discovery DE table.
        de_b (pandas.DataFrame): validation DE table.
        fdr (float): significance threshold for both.
        gene_col (str): identifier column.
        suffixes (tuple[str, str]): column-name suffixes for the two cohorts.

    Returns:
        pandas.DataFrame: columns gene, log2FC_<a>, fdr_<a>, log2FC_<b>,
            fdr_<b>, direction ('up' | 'down'), sorted by combined evidence
            (max of the two FDRs, ascending).
    """
    if suffixes is None:
        suffixes = ("discovery", "validation")
    sa, sb = suffixes
    a = de_a.set_index(gene_col)
    b = de_b.set_index(gene_col)
    shared = a.index.intersection(b.index)
    A, B = a.loc[shared], b.loc[shared]
    sel = ((A["fdr"] < fdr) & (B["fdr"] < fdr)
           & (np.sign(A["log2FC"]) == np.sign(B["log2FC"]))
           & (A["log2FC"] != 0)).values
    out = pd.DataFrame({
        "gene": shared[sel],
        "log2FC_%s" % sa: A["log2FC"].values[sel],
        "fdr_%s" % sa: A["fdr"].values[sel],
        "log2FC_%s" % sb: B["log2FC"].values[sel],
        "fdr_%s" % sb: B["fdr"].values[sel],
    })
    out["direction"] = np.where(out["log2FC_%s" % sa] > 0, "up", "down")
    out["_worst_fdr"] = np.maximum(out["fdr_%s" % sa], out["fdr_%s" % sb])
    out = out.sort_values("_worst_fdr").drop(columns="_worst_fdr")
    return out.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Microarray: series-matrix expression tables, platform annotation, probe
# collapsing, and log-scale detection
# ---------------------------------------------------------------------------

MICROARRAY_LOG_MAX = 30.0      # log2 intensities essentially never exceed this
MICROARRAY_LINEAR_MIN = 100.0  # linear intensities essentially always exceed this


def list_series_matrix_files(gse, timeout=None):
    """List the series-matrix filenames GEO publishes for a series.

    A multi-platform (SuperSeries or multi-array) submission has no single
    `<GSE>_series_matrix.txt.gz`; it has one file per platform, named
    `<GSE>-<GPL>_series_matrix.txt.gz`. Call this before `fetch_series_matrix`
    when you do not already know which is the case.

    Args:
        gse (str): series accession.
        timeout (int | None): socket timeout in seconds.

    Returns:
        list[str]: filenames, sorted.

    Raises:
        ValueError: if the series has no matrix directory or it lists no files.
    """
    if timeout is None:
        timeout = DEFAULT_TIMEOUT
    url = "%s/matrix/" % gse_ftp_dir(gse)
    try:
        html = h_http_get(url, timeout=timeout).decode("utf-8", "replace")
    except Exception as exc:                              # noqa: BLE001
        raise ValueError("cannot list matrix directory for %s: %s" % (gse, exc))
    names = sorted(set(re.findall(r'href="([^"?/][^"]*_series_matrix\.txt\.gz)"',
                                  html)))
    if not names:
        raise ValueError("no series_matrix files listed at %s" % url)
    return names


def series_matrix_platform(path):
    """Return the platform accession(s) recorded in a series matrix.

    Args:
        path (str): path to a (optionally gzipped) series matrix file.

    Returns:
        str: the single GPL accession, or `"GPL1|GPL2"` if the samples in the
            file do not all share one platform.

    Raises:
        ValueError: if no `!Sample_platform_id` line is present.
    """
    with h_open_text(path) as fh:
        for line in fh:
            if line.startswith("!Sample_platform_id"):
                _, vals = h_split_meta_line(line)
                uniq = sorted(set(v for v in vals if v))
                return uniq[0] if len(uniq) == 1 else "|".join(uniq)
            if line.startswith("!series_matrix_table_begin"):
                break
    raise ValueError("no !Sample_platform_id line found in %s" % path)


def load_series_matrix_expression(path, aggregate=None):
    """Load the value table embedded in a GEO series matrix file.

    The table lies between `!series_matrix_table_begin` and
    `!series_matrix_table_end`; its first column is the platform's probe /
    feature id (`ID_REF`) and its remaining columns are GSM accessions. For an
    Affymetrix or Gene-ST series this is probe-level data, not gene-level -
    collapse it with `collapse_probes_to_genes` before any cross-platform
    comparison.

    Args:
        path (str): path to a (optionally gzipped) series matrix file.
        aggregate ('sum' | 'max' | 'mean' | 'first' | None): how to collapse
            duplicate feature ids. None (default) raises on duplicates rather
            than guessing - probe ids are meant to be unique.

    Returns:
        pandas.DataFrame: features (index) x samples (GSM columns), float64.

    Raises:
        ValueError: if no table is present, the table has no value rows, feature
            ids are duplicated while `aggregate is None`, or the frame fails
            `h_check_matrix`.
    """
    with h_open_text(path) as fh:
        lines = []
        state = 0
        for line in fh:
            if state == 0:
                if line.startswith("!series_matrix_table_begin"):
                    state = 1
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            lines.append(line)
    if state == 0:
        raise ValueError("no !series_matrix_table_begin marker in %s" % path)
    if len(lines) < 2:
        raise ValueError("series matrix %s carries no value rows (normal for a "
                         "SuperSeries stub - fetch the per-platform matrix "
                         "files instead, see list_series_matrix_files)" % path)
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", low_memory=False)
    df = df.set_index(df.columns[0])
    df.index = [str(i).strip().strip('"') for i in df.index]
    df.index.name = "ID_REF"
    df.columns = [str(c).strip().strip('"') for c in df.columns]
    df = df.apply(pd.to_numeric, errors="coerce").astype("float64")
    df = df.loc[:, ~df.isna().all(axis=0)]
    idx = pd.Index(df.index)
    if idx.duplicated().any():
        dups = sorted(set(idx[idx.duplicated()]))
        if aggregate is None:
            raise ValueError("duplicate feature ids in %s (%d, e.g. %s); pass "
                             "aggregate= to collapse them deliberately"
                             % (path, len(dups), dups[:5]))
        df = (df.groupby(level=0).first() if aggregate == "first"
              else getattr(df.groupby(level=0), aggregate)())
    h_check_matrix(df)
    return df.astype("float64")


def fetch_platform_annotation(gpl, dest_dir, timeout=None, overwrite=False):
    """Download a GEO platform annotation (`GPLxxx.annot.gz`).

    Args:
        gpl (str): platform accession, e.g. "GPL96".
        dest_dir (str): directory to write into (created if absent).
        timeout (int | None): socket timeout in seconds.
        overwrite (bool): re-download even if present.

    Returns:
        str: local path to `<GPL>.annot.gz`.

    Raises:
        ValueError: if `gpl` is not a GPL accession. Platforms with no curated
            annotation file (several sequencing platforms, e.g. GPL9052) raise
            urllib's HTTPError 404 - there is no probe set to annotate.
    """
    if timeout is None:
        timeout = DEFAULT_TIMEOUT
    gpl = str(gpl).strip().upper()
    m = re.fullmatch(r"GPL(\d+)", gpl)
    if m is None:
        raise ValueError("not a GPL accession: %r" % (gpl,))
    digits = m.group(1)
    stub = "GPL%snnn" % (digits[:-3] if len(digits) > 3 else "")
    os.makedirs(dest_dir, exist_ok=True)
    name = "%s.annot.gz" % gpl
    url = "%s/platforms/%s/%s/annot/%s" % (GEO_FTP_ROOT, stub, gpl, name)
    return h_download(url, os.path.join(dest_dir, name), timeout=timeout,
                      overwrite=overwrite)


def parse_platform_annotation(path, id_col="ID", symbol_col="Gene symbol"):
    """Parse the probe -> gene-symbol mapping out of a GEO `.annot.gz` file.

    The mapping is returned verbatim, including GEO's `///`-joined multi-gene
    entries (e.g. `"MIR4640///DDR1"`). Resolving those is
    `collapse_probes_to_genes`' job, under an explicit policy - this function
    never silently picks the first symbol.

    Args:
        path (str): path to a (optionally gzipped) `GPLxxx.annot.gz`.
        id_col (str): probe id column in the platform table.
        symbol_col (str): gene symbol column.

    Returns:
        pandas.Series: index = probe id, values = raw symbol string. Probes with
            an empty symbol are dropped.

    Raises:
        ValueError: if the platform-table markers or requested columns are
            absent, or probe ids are duplicated.
    """
    with h_open_text(path) as fh:
        lines = []
        state = 0
        for line in fh:
            if state == 0:
                if line.startswith("!platform_table_begin"):
                    state = 1
                continue
            if line.startswith("!platform_table_end"):
                break
            lines.append(line)
    if state == 0:
        raise ValueError("no !platform_table_begin marker in %s" % path)
    df = pd.read_csv(io.StringIO("".join(lines)), sep="\t", low_memory=False,
                     dtype=str)
    for c in (id_col, symbol_col):
        if c not in df.columns:
            raise ValueError("annotation %s has no column %r (columns: %s)"
                             % (path, c, list(df.columns)[:12]))
    s = pd.Series(df[symbol_col].values,
                  index=[str(i).strip() for i in df[id_col].values])
    if s.index.duplicated().any():
        raise ValueError("duplicate probe ids in annotation %s" % path)
    s = s.fillna("").map(lambda v: str(v).strip())
    return s[s != ""]


def detect_log_scale(expr, log_max=None, linear_min=None):
    """Decide whether an expression matrix is already on a log scale.

    Detection, not assumption: microarray series matrices ship RMA/GC-RMA log2
    values, MAS5 linear intensities, or occasionally something else, and the
    series metadata frequently does not say which. Getting it wrong applies an
    unwanted monotone transform to every log2FC.

    Rules, in order (the one that fires is reported):

    * any value < 0             -> logged  (`logged:negative_values`)
    * max <= `log_max` (30)      -> logged  (`logged:max_le_30`)
    * max >= `linear_min` (100)  -> linear  (`linear:max_ge_100`)
    * otherwise                  -> raise   (nothing decides it)

    Args:
        expr (pandas.DataFrame): features x samples.
        log_max (float | None): upper bound on a plausible log2 maximum.
        linear_min (float | None): lower bound on a plausible linear maximum.

    Returns:
        dict: `{'already_logged': bool, 'branch': str, 'max': float,
            'min': float, 'median': float, 'n_negative': int}`.

    Raises:
        ValueError: if the matrix has no finite values, or its maximum falls in
            the undecidable band between `log_max` and `linear_min`. In that
            case establish the scale from the series' own documentation and pass
            `already_logged=` explicitly - do not let this function guess.
    """
    if log_max is None:
        log_max = MICROARRAY_LOG_MAX
    if linear_min is None:
        linear_min = MICROARRAY_LINEAR_MIN
    v = np.asarray(expr.values, dtype="float64")
    finite = v[np.isfinite(v)]
    if finite.size == 0:
        raise ValueError("detect_log_scale: matrix has no finite values")
    vmax = float(finite.max())
    n_neg = int((finite < 0).sum())
    info = {"max": vmax, "min": float(finite.min()),
            "median": float(np.median(finite)), "n_negative": n_neg}
    if n_neg > 0:
        info.update(already_logged=True, branch="logged:negative_values")
    elif vmax <= log_max:
        info.update(already_logged=True, branch="logged:max_le_%g" % log_max)
    elif vmax >= linear_min:
        info.update(already_logged=False, branch="linear:max_ge_%g" % linear_min)
    else:
        raise ValueError(
            "detect_log_scale cannot decide the scale: max=%.4g lies between "
            "log_max=%.4g and linear_min=%.4g. Establish the scale from the "
            "series documentation and pass already_logged= explicitly."
            % (vmax, log_max, linear_min))
    return info


def collapse_probes_to_genes(expr, probe_to_gene, method="max_mean",
                             ambiguous="raise", multi_gene_sep="///",
                             gene_name="gene"):
    """Collapse a probe-level matrix to one row per gene.

    Default `method="max_mean"`: among the probes mapping to a gene, keep the
    single probe with the highest mean value across samples. On a log2 intensity
    scale that is the highest-mean-intensity probe - the conventional choice (it
    favours the probe least dominated by background), and it keeps one real
    measurement per gene rather than averaging probes of different affinity.

    Ambiguous probes - a probe whose annotation names several genes, e.g. GEO's
    `"MIR4640///DDR1"` - are **not** resolved by taking the first symbol.
    `ambiguous="raise"` (default) raises; `ambiguous="drop"` excludes them and
    reports the count. Silent first-picking is not offered.

    Selection is deterministic: ties on the ranking statistic are broken by
    probe id, so the result does not depend on input row order.

    Args:
        expr (pandas.DataFrame): probes x samples, numeric.
        probe_to_gene (pandas.Series | dict | pandas.DataFrame): probe id ->
            gene symbol. A DataFrame must have exactly one column.
        method ('max_mean' | 'mean' | 'median' | 'max_var'): collapsing rule.
            'max_mean' / 'max_var' select one probe; 'mean' / 'median' average
            across probes.
        ambiguous ('raise' | 'drop'): what to do with multi-gene probes.
        multi_gene_sep (str): separator marking a multi-gene annotation.
        gene_name (str): name given to the resulting index.

    Returns:
        (pandas.DataFrame, dict): genes x samples (float64, unique index,
            input column order preserved) and a report with keys
            `method, n_probes_in, n_probes_annotated, n_probes_unannotated,
            n_probes_ambiguous, ambiguous_policy, ambiguous_examples,
            n_genes_out, n_genes_multiprobe, max_probes_per_gene`.

    Raises:
        ValueError: on an unknown `method` / `ambiguous`, duplicate probe ids in
            the mapping, a mapping with no overlap with `expr.index`, an
            ambiguous mapping under `ambiguous="raise"`, or an empty result.
    """
    if method not in ("max_mean", "mean", "median", "max_var"):
        raise ValueError("unknown collapsing method %r" % (method,))
    if ambiguous not in ("raise", "drop"):
        raise ValueError("ambiguous must be 'raise' or 'drop', got %r"
                         % (ambiguous,))
    if isinstance(probe_to_gene, pd.DataFrame):
        if probe_to_gene.shape[1] != 1:
            raise ValueError("probe_to_gene DataFrame must have one column, "
                             "got %s" % list(probe_to_gene.columns))
        probe_to_gene = probe_to_gene.iloc[:, 0]
    m = pd.Series(probe_to_gene, dtype=object)
    m.index = [str(i).strip() for i in m.index]
    if m.index.duplicated().any():
        raise ValueError("probe_to_gene has duplicate probe ids")
    X = expr.copy()
    X.index = [str(i).strip() for i in X.index]
    X = X.apply(pd.to_numeric, errors="coerce").astype("float64")

    n_in = int(X.shape[0])
    mapped = m.reindex(X.index)
    mapped = mapped.map(lambda v: "" if v is None
                        or (isinstance(v, float) and np.isnan(v))
                        else str(v).strip())
    annotated = mapped != ""
    if int(annotated.sum()) == 0:
        raise ValueError("probe_to_gene shares no probe ids with the expression "
                         "matrix - check that the annotation platform matches "
                         "the series platform")
    amb_mask = annotated & mapped.map(lambda v: multi_gene_sep in v)
    n_amb = int(amb_mask.sum())
    amb_examples = [(p, mapped[p]) for p in list(mapped.index[amb_mask])[:5]]
    if n_amb and ambiguous == "raise":
        raise ValueError(
            "%d of %d annotated probes map to more than one gene (separator "
            "%r), e.g. %s. Refusing to pick the first symbol silently. Pass "
            "ambiguous='drop' to exclude them, or supply a one-gene-per-probe "
            "mapping." % (n_amb, int(annotated.sum()), multi_gene_sep,
                          amb_examples))
    keep = (annotated & ~amb_mask).values
    Xk = X.loc[keep]
    genes = mapped[keep].values

    if method in ("mean", "median"):
        out = getattr(Xk.groupby(genes), method)()
    else:
        stat = (Xk.mean(axis=1) if method == "max_mean"
                else Xk.var(axis=1, ddof=1))
        order = pd.DataFrame({"gene": genes, "stat": stat.values,
                              "probe": list(Xk.index)}, index=list(Xk.index))
        order = order.sort_values(["gene", "stat", "probe"],
                                  ascending=[True, False, True],
                                  kind="mergesort")
        chosen = order.groupby("gene", sort=True).head(1)
        out = Xk.loc[chosen["probe"].values]
        out.index = chosen["gene"].values
    out = out[[c for c in expr.columns]]
    out.index = pd.Index([str(g) for g in out.index], name=gene_name)
    out = out.astype("float64")
    if out.shape[0] == 0:
        raise ValueError("probe collapsing produced an empty matrix")
    if out.index.duplicated().any():
        raise ValueError("internal error: collapsed matrix has duplicate genes")
    per_gene = pd.Series(genes).value_counts()
    report = {
        "method": method,
        "n_probes_in": n_in,
        "n_probes_annotated": int(annotated.sum()),
        "n_probes_unannotated": int((~annotated).sum()),
        "n_probes_ambiguous": n_amb,
        "ambiguous_policy": ambiguous,
        "ambiguous_examples": amb_examples,
        "n_genes_out": int(out.shape[0]),
        "n_genes_multiprobe": int((per_gene > 1).sum()),
        "max_probes_per_gene": int(per_gene.max()) if len(per_gene) else 0,
    }
    return out, report


def run_de_microarray(expr, groups, case, control, already_logged=None,
                      log_max=None, linear_min=None, filter_min_value=None,
                      filter_min_fraction=0.2, covariates=None,
                      min_per_group=2, trend_frac=0.4):
    """Differential expression on a microarray matrix, with scale detection.

    A thin wrapper over `run_de(design="limma_trend")`. Its only job is to
    decide - or record - whether the input is already log2, and to report which
    branch fired so the choice is auditable rather than implicit.

    Args:
        expr (pandas.DataFrame): genes x samples (collapse probes first).
        groups (pandas.Series | dict | list): group label per sample.
        case (str): log2FC numerator label.
        control (str): denominator label.
        already_logged (bool | None): None = call `detect_log_scale`; True/False
            overrides detection and the branch is recorded as
            `caller:already_logged=...`.
        log_max, linear_min (float | None): passed to `detect_log_scale`.
        filter_min_value (float | None): expression filter on the input scale.
            None (default) = no filtering. An RMA matrix is already
            background-corrected, and thresholding a log2 matrix at 1.0 would
            drop essentially nothing while a threshold near the array's noise
            floor discards real low-expressed genes - so leave this None unless
            you have a specific reason.
        filter_min_fraction (float): fraction of the smaller group that must
            exceed `filter_min_value`.
        covariates (pandas.DataFrame | None): per-sample numeric covariates.
        min_per_group (int): precondition on group sizes.
        trend_frac (float): lowess span for the variance trend.

    Returns:
        pandas.DataFrame: as `run_de`, with `.attrs['log_scale_branch']` and
            `.attrs['log_scale_info']` added.

    Raises:
        ValueError: as `detect_log_scale` and `run_de`.
    """
    if already_logged is None:
        info = detect_log_scale(expr, log_max=log_max, linear_min=linear_min)
        logged = bool(info["already_logged"])
        branch = info["branch"]
    else:
        logged = bool(already_logged)
        try:
            info = detect_log_scale(expr, log_max=log_max,
                                    linear_min=linear_min)
        except ValueError as exc:
            info = {"detection_error": str(exc)}
        branch = "caller:already_logged=%s" % logged
    out = run_de(expr, groups, case, control, design="limma_trend",
                 covariates=covariates, min_per_group=min_per_group,
                 filter_min_value=filter_min_value,
                 filter_min_fraction=filter_min_fraction,
                 trend_frac=trend_frac, already_logged=logged)
    out.attrs["log_scale_branch"] = branch
    out.attrs["log_scale_info"] = info
    return out


def label_permutation_control(expr, groups, case, control, seed, n_perm=1,
                              fdr=0.05, de_fn=None, **de_kwargs):
    """Re-run a DE contrast on shuffled group labels and count the hits.

    Under a correct pipeline a permuted contrast yields no genes at
    FDR<`fdr`. A non-zero count means the variance model, the FDR, or the label
    handling is wrong, and no unpermuted number from the same pipeline can be
    trusted.

    The shuffle is a permutation of the label vector, so both group sizes are
    preserved and the only thing destroyed is the association between label and
    expression.

    Args:
        expr (pandas.DataFrame): genes x samples.
        groups (pandas.Series | dict | list): true labels.
        case (str): case label.
        control (str): control label.
        seed (int): RNG seed - required, not defaulted, so a control is always
            reproducible.
        n_perm (int): number of permutations.
        fdr (float): significance threshold.
        de_fn (callable | None): DE function; default `run_de_microarray`.
        **de_kwargs: passed through to `de_fn`.

    Returns:
        list[dict]: one row per permutation with `permutation, seed, n_sig,
            n_genes_tested, min_fdr, min_p_value, frac_p_below_05, n_case,
            n_ctrl`.

    Raises:
        ValueError: if `n_perm < 1`, or from the underlying DE call.
    """
    if n_perm < 1:
        raise ValueError("n_perm must be >= 1, got %r" % (n_perm,))
    if de_fn is None:
        de_fn = run_de_microarray
    g = groups.copy() if isinstance(groups, pd.Series) else pd.Series(groups)
    g.index = [str(i) for i in g.index]
    g = g.reindex([str(c) for c in expr.columns])
    keep = g[g.isin([case, control])].index.tolist()
    sub = expr[keep]
    labels = np.asarray(g[keep].values).copy()
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(int(n_perm)):
        perm = pd.Series(rng.permutation(labels), index=keep)
        de = de_fn(sub, perm, case, control, **de_kwargs)
        p = np.asarray(de["p_value"].values, dtype="float64")
        rows.append({
            "permutation": i + 1,
            "seed": int(seed),
            "n_sig": int((de["fdr"].values < fdr).sum()),
            "n_genes_tested": int(de.shape[0]),
            "min_fdr": float(np.nanmin(de["fdr"].values)),
            "min_p_value": float(np.nanmin(p)),
            "frac_p_below_05": float(np.mean(p < 0.05)),
            "n_case": int((perm == case).sum()),
            "n_ctrl": int((perm == control).sum()),
        })
    return rows
