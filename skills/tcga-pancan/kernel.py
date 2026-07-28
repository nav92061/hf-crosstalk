"""Memory-bounded access to the TCGA PanCanAtlas / PanImmune data products on GDC.

Pure stdlib + pandas/numpy. No agent host calls, no platform dependencies.
Designed for machines with a few GiB of RAM: the ~1.9 GB EBPlusPlus expression
matrix is never loaded whole, only streamed line by line.
"""

import csv
import gzip
import io
import json
import os
import re
import sys
import time
import urllib.request

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# GDC PanCanAtlas / PanImmune file endpoints (open access, no token required).
# UUIDs taken from https://gdc.cancer.gov/about-data/publications/panimmune
# ---------------------------------------------------------------------------

GDC_DATA_BASE = "https://api.gdc.cancer.gov/data/"

# EBPlusPlusAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.tsv  (~1.9 GB, plain TSV)
# Values are RSEM normalized counts (upper-quartile-scaled, batch-corrected).
# They are NOT TPM and NOT FPKM. Do not relabel them.
GDC_EXPRESSION_URL = "https://api.gdc.cancer.gov/data/3586c0da-64d0-4b74-a449-5ff4d9136611"
EXPRESSION_UNIT = (
    "RSEM normalized_count (EBPlusPlus batch-corrected, upper-quartile scaled); "
    "NOT TPM, NOT FPKM"
)

# merged_sample_quality_annotations.tsv - aliquot barcode -> cancer type + QC flags.
GDC_SAMPLE_ANNOTATION_URL = "https://api.gdc.cancer.gov/data/1a7d7be8-675d-4e60-a105-19d4121bdebf"

# TCGA_all_leuk_estimate.masked.20170107.tsv - Thorsson et al. leukocyte fraction
# (methylation-based). Headerless: cancer_type, aliquot_barcode, leukocyte_fraction.
GDC_LEUKOCYTE_FRACTION_URL = "https://api.gdc.cancer.gov/data/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8"

# Scores_160_Signatures.tsv.gz - PanImmune signature scores, samples in COLUMNS.
GDC_SIGNATURE_SCORES_URL = "https://api.gdc.cancer.gov/data/80a82092-161d-4615-9d96-e858f113618d"

# TCGA_mastercalls.abs_tables_JSedit.fixed.txt - ABSOLUTE purity / ploidy.
GDC_ABSOLUTE_PURITY_URL = "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5"

# TCGA-CDR clinical outcome endpoints (xlsx).
GDC_TCGA_CDR_URL = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81"

GDC_CASES_API = "https://api.gdc.cancer.gov/cases"

# TCGA sample type codes: barcode characters 14-15 (1-based), i.e. positions
# 13-14 zero-based, the two digits of the 4th hyphen-delimited field.
TCGA_SAMPLE_TYPES = {
    "01": "Primary Solid Tumor",
    "02": "Recurrent Solid Tumor",
    "03": "Primary Blood Derived Cancer - Peripheral Blood",
    "05": "Additional New Primary",
    "06": "Metastatic",
    "07": "Additional Metastatic",
    "09": "Primary Blood Derived Cancer - Bone Marrow",
    "10": "Blood Derived Normal",
    "11": "Solid Tissue Normal",
    "12": "Buccal Cell Normal",
    "14": "Bone Marrow Normal",
}

PRIMARY_SOLID_TUMOR_CODE = "01"

# TCGA study abbreviation -> human-readable disease name.
TCGA_STUDY_NAMES = {
    "ACC": "Adrenocortical carcinoma",
    "BLCA": "Bladder urothelial carcinoma",
    "BRCA": "Breast invasive carcinoma",
    "CESC": "Cervical squamous cell carcinoma and endocervical adenocarcinoma",
    "CHOL": "Cholangiocarcinoma",
    "COAD": "Colon adenocarcinoma",
    "DLBC": "Lymphoid neoplasm diffuse large B-cell lymphoma",
    "ESCA": "Esophageal carcinoma",
    "GBM": "Glioblastoma multiforme",
    "HNSC": "Head and neck squamous cell carcinoma",
    "KICH": "Kidney chromophobe",
    "KIRC": "Kidney renal clear cell carcinoma",
    "KIRP": "Kidney renal papillary cell carcinoma",
    "LAML": "Acute myeloid leukemia",
    "LGG": "Brain lower grade glioma",
    "LIHC": "Liver hepatocellular carcinoma",
    "LUAD": "Lung adenocarcinoma",
    "LUSC": "Lung squamous cell carcinoma",
    "MESO": "Mesothelioma",
    "OV": "Ovarian serous cystadenocarcinoma",
    "PAAD": "Pancreatic adenocarcinoma",
    "PCPG": "Pheochromocytoma and paraganglioma",
    "PRAD": "Prostate adenocarcinoma",
    "READ": "Rectum adenocarcinoma",
    "SARC": "Sarcoma",
    "SKCM": "Skin cutaneous melanoma",
    "STAD": "Stomach adenocarcinoma",
    "TGCT": "Testicular germ cell tumors",
    "THCA": "Thyroid carcinoma",
    "THYM": "Thymoma",
    "UCEC": "Uterine corpus endometrial carcinoma",
    "UCS": "Uterine carcinosarcoma",
    "UVM": "Uveal melanoma",
}

# Studies excluded by default from "solid tumor" analyses.
# LAML (acute myeloid leukemia) and DLBC (diffuse large B-cell lymphoma) are
# haematological malignancies: their "tumor" aliquots are blood/marrow derived,
# so bulk expression reflects circulating leukocytes rather than a solid tumor
# receiving signals from perfused interstitium. THYM (thymoma) is a solid mass
# but is overwhelmingly composed of non-neoplastic thymocytes, so receptor
# expression cannot be attributed to a tumor compartment. A question about
# solid tumors exposed to circulating cardiac-derived factors is not answerable
# in these three cohorts, so they are dropped rather than silently mis-scored.
DEFAULT_EXCLUDED_STUDIES = ("LAML", "DLBC", "THYM")

HOUSEKEEPING_GENES = (
    "ACTB", "GAPDH", "TUBB", "RPLP0", "B2M", "PGK1", "PPIA",
    "TBP", "RPL13A", "SDHA", "UBC", "YWHAZ", "HPRT1", "EEF1A1",
)


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

def download_to_disk(url, dest, force=False, timeout=120, chunk_bytes=None):
    """Stream a URL to a local file without holding it in memory.

    Args:
        url (str): source URL.
        dest (str): destination path on disk.
        force (bool): re-download even if dest already exists and is non-empty.
        timeout (int): socket timeout in seconds.
        chunk_bytes (int): read block size.

    Returns:
        str: dest.
    """
    if chunk_bytes is None:
        chunk_bytes = 8 << 20
    if (not force) and os.path.exists(dest) and os.path.getsize(dest) > 0:
        return dest
    tmp = dest + ".part"
    d = os.path.dirname(os.path.abspath(dest))
    if d:
        os.makedirs(d, exist_ok=True)
    with urllib.request.urlopen(url, timeout=timeout) as resp, open(tmp, "wb") as fh:
        while True:
            block = resp.read(chunk_bytes)
            if not block:
                break
            fh.write(block)
    os.replace(tmp, dest)
    assert os.path.getsize(dest) > 0, "downloaded file is empty: %s" % url
    return dest


def h_open_text_stream(url_or_path, timeout=120):
    """Open a local path or URL as a decoded text stream, transparently gunzipping.

    Returns:
        tuple: (text_iterator_context, closer_callable)
    """
    if re.match(r"^https?://", str(url_or_path)):
        raw = urllib.request.urlopen(url_or_path, timeout=timeout)
        head = raw.read(2)
        gz = head[:2] == b"\x1f\x8b"
        body = io.BufferedReader(PrefixedReader(head, raw), buffer_size=1 << 20)
        binary = gzip.GzipFile(fileobj=body) if gz else body
        closer = raw.close
    else:
        with open(url_or_path, "rb") as probe:
            gz = probe.read(2)[:2] == b"\x1f\x8b"
        binary = (gzip.open(url_or_path, "rb") if gz
                  else open(url_or_path, "rb", buffering=1 << 20))
        closer = binary.close
    return io.TextIOWrapper(binary, encoding="utf-8", errors="replace"), closer


def h_prefixed_read_factory(prefix, stream):
    """Internal: build a file-like object that replays `prefix` before `stream`."""
    state = {"buf": bytes(prefix)}

    def read(n=-1):
        if state["buf"]:
            if n is None or n < 0:
                out = state["buf"] + stream.read()
                state["buf"] = b""
                return out
            take = state["buf"][:n]
            state["buf"] = state["buf"][n:]
            if len(take) < n:
                take += stream.read(n - len(take))
            return take
        return stream.read() if (n is None or n < 0) else stream.read(n)

    return read


def PrefixedReader(prefix, stream):  # noqa: N802 - factory, not a class
    """Internal: RawIOBase-compatible shim replaying sniffed bytes.

    Implemented as a factory returning an io.RawIOBase instance so that
    kernel.py keeps no top-level class definitions.
    """
    read = h_prefixed_read_factory(prefix, stream)
    obj = io.RawIOBase()
    obj.readable = lambda: True
    obj.read = read
    obj.readinto = None  # force BufferedReader to use .read()

    def readinto(b):
        data = read(len(b))
        n = len(data)
        b[:n] = data
        return n

    obj.readinto = readinto
    obj.close = lambda: None
    return obj


def h_strip(tok):
    """Internal: strip surrounding double quotes and whitespace from a TSV token."""
    tok = tok.strip()
    if len(tok) >= 2 and tok[0] == '"' and tok[-1] == '"':
        tok = tok[1:-1]
    return tok


# ---------------------------------------------------------------------------
# Streaming expression extraction
# ---------------------------------------------------------------------------

def stream_expression_rows(url_or_path, wanted_genes, dest_csv,
                          chunk_bytes=None, gene_id_sep="|",
                          progress_every=5000, verbose=True):
    """Stream the PanCanAtlas gene x aliquot TSV and keep only requested genes.

    The matrix is read one line at a time; at no point is the full file held in
    memory. Only the retained rows are buffered, so peak RSS scales with
    len(wanted_genes) x n_samples, not with the file size.

    Row format expected: first column `gene_id` as "SYMBOL|ENTREZ", remaining
    columns TCGA aliquot barcodes. Values are RSEM normalized counts
    (see EXPRESSION_UNIT) - they are not TPM.

    Args:
        url_or_path (str): HTTP(S) URL or local path; gzip auto-detected.
        wanted_genes (iterable of str): HGNC symbols to retain (case-sensitive
            match on the symbol part of gene_id, with an upper-case fallback).
        dest_csv (str): output CSV path. Written with the gene symbol as the
            first column ("gene") followed by one column per aliquot barcode.
            If it ends in .gz the output is gzipped.
        chunk_bytes (int): text-stream read buffer hint.
        gene_id_sep (str): separator inside gene_id.
        progress_every (int): print a progress marker every N input lines
            (0 disables).
        verbose (bool): print progress markers to stderr.

    Returns:
        dict: {
          'dest_csv': str,
          'n_input_rows': int,
          'n_rows_kept': int,
          'genes_requested': int,
          'genes_found': sorted list[str],
          'genes_missing': sorted list[str],
          'n_samples': int,
          'samples': list[str],
          'expression_unit': str,
          'duplicate_symbols': sorted list[str],
        }

    Raises:
        AssertionError: if wanted_genes is empty, if the header cannot be
            parsed, if zero requested genes matched, or if duplicate aliquot
            barcodes appear in the header.
    """
    if chunk_bytes is None:
        chunk_bytes = 1 << 20
    wanted = {str(g).strip() for g in wanted_genes if str(g).strip()}
    assert wanted, "wanted_genes is empty - refusing to stream a 1.9 GB matrix for nothing"
    wanted_upper = {g.upper(): g for g in wanted}

    stream, closer = h_open_text_stream(url_or_path)
    n_input = 0
    kept_rows = []
    kept_symbols = []
    dup_symbols = []
    seen = set()
    try:
        header_line = stream.readline()
        assert header_line, "expression stream is empty"
        header = [h_strip(t) for t in header_line.rstrip("\r\n").split("\t")]
        assert len(header) > 1, "header has no sample columns: %r" % header[:3]
        samples = header[1:]
        assert len(samples) == len(set(samples)), (
            "duplicate aliquot barcodes in expression header (%d cols, %d unique)"
            % (len(samples), len(set(samples)))
        )

        for line in stream:
            n_input += 1
            if progress_every and verbose and n_input % progress_every == 0:
                sys.stderr.write("  streamed %d rows, kept %d\n" % (n_input, len(kept_rows)))
                sys.stderr.flush()
            tab = line.find("\t")
            if tab < 0:
                continue
            gene_id = h_strip(line[:tab])
            symbol = gene_id.split(gene_id_sep)[0]
            match = None
            if symbol in wanted:
                match = symbol
            elif symbol.upper() in wanted_upper:
                match = wanted_upper[symbol.upper()]
            if match is None:
                continue
            if match in seen:
                dup_symbols.append(match)
                continue
            seen.add(match)
            vals = line[tab + 1:].rstrip("\r\n").split("\t")
            if len(vals) != len(samples):
                vals = (vals + [""] * len(samples))[:len(samples)]
            kept_symbols.append(match)
            kept_rows.append(vals)
    finally:
        try:
            stream.close()
        finally:
            closer()

    assert kept_rows, (
        "streamed %d rows from %s but matched 0 of %d requested genes - "
        "check symbol casing or gene_id format" % (n_input, url_or_path, len(wanted))
    )

    opener = gzip.open if str(dest_csv).endswith(".gz") else open
    with opener(dest_csv, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene"] + samples)
        for sym, vals in zip(kept_symbols, kept_rows):
            w.writerow([sym] + vals)

    found = sorted(seen)
    missing = sorted(wanted - seen)
    if verbose:
        sys.stderr.write("stream done: %d input rows, %d genes kept, %d missing\n"
                         % (n_input, len(found), len(missing)))
    return {
        "dest_csv": str(dest_csv),
        "n_input_rows": n_input,
        "n_rows_kept": len(kept_rows),
        "genes_requested": len(wanted),
        "genes_found": found,
        "genes_missing": missing,
        "n_samples": len(samples),
        "samples": samples,
        "expression_unit": EXPRESSION_UNIT,
        "duplicate_symbols": sorted(set(dup_symbols)),
    }


def load_expression_csv(path, index_col="gene"):
    """Load a CSV written by stream_expression_rows into a float DataFrame.

    Args:
        path (str): CSV or CSV.GZ path.
        index_col (str): gene column name.

    Returns:
        pandas.DataFrame: genes (rows) x aliquot barcodes (columns), float64.
    """
    df = pd.read_csv(path, index_col=index_col)
    return df.apply(pd.to_numeric, errors="coerce")


# ---------------------------------------------------------------------------
# Barcode parsing and sample mapping
# ---------------------------------------------------------------------------

def parse_tcga_barcode(barcode):
    """Decompose a TCGA barcode into patient / sample-type components.

    Sample type is read from barcode characters 14-15 (zero-based positions
    13-14), i.e. the two leading digits of the 4th hyphen-delimited field:
    "01" primary solid tumor, "11" solid tissue normal, "06" metastatic,
    "03" primary blood derived cancer (peripheral blood).

    Args:
        barcode (str): any TCGA barcode of patient depth or deeper, e.g.
            "TCGA-OR-A5J1-01A-11R-A29S-07".

    Returns:
        dict: {
          'barcode': str, 'patient': str (TCGA-XX-YYYY), 'tss': str,
          'participant': str, 'sample_type_code': str|None, 'vial': str|None,
          'sample_type': str|None (human-readable), 'is_primary_tumor': bool,
          'is_tumor': bool, 'is_normal': bool,
        }

    Raises:
        AssertionError: if the barcode has fewer than 3 hyphen-delimited fields.
    """
    bc = str(barcode).strip()
    parts = bc.split("-")
    assert len(parts) >= 3, "not a TCGA barcode: %r" % bc
    patient = "-".join(parts[:3])
    code = None
    vial = None
    if len(parts) >= 4 and len(parts[3]) >= 2 and parts[3][:2].isdigit():
        code = parts[3][:2]
        vial = parts[3][2:] or None
    st = TCGA_SAMPLE_TYPES.get(code)
    is_normal = code in ("10", "11", "12", "14") if code else False
    return {
        "barcode": bc,
        "patient": patient,
        "tss": parts[1],
        "participant": parts[2],
        "sample_type_code": code,
        "vial": vial,
        "sample_type": st,
        "is_primary_tumor": code == PRIMARY_SOLID_TUMOR_CODE,
        "is_tumor": (code is not None) and (not is_normal),
        "is_normal": is_normal,
    }


def build_sample_map(thorsson_df, aliquot_col=None, study_col=None):
    """Map aliquot barcodes to TCGA study abbreviation and full disease name.

    Args:
        thorsson_df (pandas.DataFrame): the PanCanAtlas
            merged_sample_quality_annotations table (or any frame carrying an
            aliquot barcode column and a cancer-type column). Column names are
            auto-detected when not given.
        aliquot_col (str|None): aliquot barcode column name.
        study_col (str|None): study abbreviation column name (e.g. "cancer type").

    Returns:
        pandas.DataFrame indexed by aliquot_barcode with columns
        ['patient', 'tumor_type', 'tumor_name', 'sample_type_code',
         'sample_type', 'is_primary_tumor'].

    Raises:
        AssertionError: if the required columns cannot be located, or if the
            resulting index contains duplicate aliquot barcodes.
    """
    df = thorsson_df
    cols = {str(c).strip().lower(): c for c in df.columns}
    if aliquot_col is None:
        for cand in ("aliquot_barcode", "aliquotbarcode", "aliquot", "sample_barcode"):
            if cand in cols:
                aliquot_col = cols[cand]
                break
    if study_col is None:
        for cand in ("cancer type", "cancer_type", "study", "disease", "project_id",
                     "cancer.type"):
            if cand in cols:
                study_col = cols[cand]
                break
    assert aliquot_col is not None, "no aliquot barcode column in %s" % list(df.columns)[:12]
    assert study_col is not None, "no cancer-type column in %s" % list(df.columns)[:12]

    sub = df[[aliquot_col, study_col]].copy()
    sub.columns = ["aliquot_barcode", "tumor_type"]
    sub["aliquot_barcode"] = sub["aliquot_barcode"].astype(str).str.strip()
    sub["tumor_type"] = (sub["tumor_type"].astype(str).str.strip().str.upper()
                         .str.replace("^TCGA-", "", regex=True))
    sub = sub[sub["aliquot_barcode"].str.startswith("TCGA-")]
    sub = sub[~sub["tumor_type"].isin(["", "NAN", "NONE"])]
    sub = sub.drop_duplicates(subset=["aliquot_barcode"], keep="first")

    parsed = pd.DataFrame([parse_tcga_barcode(b) for b in sub["aliquot_barcode"]])
    out = pd.DataFrame({
        "aliquot_barcode": sub["aliquot_barcode"].values,
        "patient": parsed["patient"].values,
        "tumor_type": sub["tumor_type"].values,
        "sample_type_code": parsed["sample_type_code"].values,
        "sample_type": parsed["sample_type"].values,
        "is_primary_tumor": parsed["is_primary_tumor"].values,
    })
    out["tumor_name"] = out["tumor_type"].map(TCGA_STUDY_NAMES).fillna(out["tumor_type"])
    out = out.set_index("aliquot_barcode")
    assert out.index.is_unique, "duplicate aliquot barcodes in sample map"
    return out[["patient", "tumor_type", "tumor_name", "sample_type_code",
                "sample_type", "is_primary_tumor"]]


def lookup_tumor_types_via_gdc(barcodes, page_size=400, timeout=90, verbose=True):
    """Fallback: resolve barcodes to TCGA projects via the GDC cases API.

    Queries https://api.gdc.cancer.gov/cases with
    fields=submitter_id,project.project_id, paginated. Barcodes are truncated to
    patient depth (TCGA-XX-YYYY) because `cases.submitter_id` is a patient id.

    Args:
        barcodes (iterable of str): aliquot or patient barcodes.
        page_size (int): records per request.
        timeout (int): socket timeout in seconds.
        verbose (bool): print progress markers to stderr.

    Returns:
        dict: patient barcode -> TCGA study abbreviation (e.g. "BRCA").
    """
    patients = sorted({parse_tcga_barcode(b)["patient"] for b in barcodes})
    resolved = {}
    for start in range(0, len(patients), page_size):
        batch = patients[start:start + page_size]
        payload = {
            "filters": {"op": "in",
                        "content": {"field": "submitter_id", "value": batch}},
            "fields": "submitter_id,project.project_id",
            "format": "JSON",
            "size": str(len(batch) * 2),
        }
        req = urllib.request.Request(
            GDC_CASES_API,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
        )
        for attempt in range(3):
            try:
                with urllib.request.urlopen(req, timeout=timeout) as resp:
                    hits = json.loads(resp.read().decode("utf-8"))["data"]["hits"]
                break
            except Exception:
                if attempt == 2:
                    raise
                time.sleep(2 ** attempt)
        for h in hits:
            proj = (h.get("project") or {}).get("project_id", "")
            if proj.startswith("TCGA-"):
                resolved[h["submitter_id"]] = proj.split("-", 1)[1]
        if verbose:
            sys.stderr.write("  GDC cases: %d/%d patients resolved\n"
                             % (len(resolved), len(patients)))
    return resolved


# ---------------------------------------------------------------------------
# Filtering and summarisation
# ---------------------------------------------------------------------------

def filter_solid_primary(expr, sample_map, exclude_studies=None,
                        require_all_mapped=True):
    """Restrict an expression matrix to primary solid tumors of solid-organ studies.

    Keeps only aliquots whose sample type code is "01" (Primary Solid Tumor) and
    drops the studies in `exclude_studies`. LAML and DLBC are haematological
    (their tumor aliquots are blood/marrow, so bulk expression measures
    circulating leukocytes, not a perfused solid tumor); THYM is dominated by
    non-neoplastic thymocytes. This project asks which SOLID tumors could
    respond to circulating cardiac-derived factors, so those cohorts are
    dropped rather than silently mis-scored.

    Args:
        expr (pandas.DataFrame): genes x aliquot barcodes.
        sample_map (pandas.DataFrame): output of build_sample_map (indexed by
            aliquot barcode, with 'tumor_type').
        exclude_studies (tuple of str): study abbreviations to drop.
        require_all_mapped (bool): if True, raise when an expression column has
            no tumor-type assignment; if False, drop such columns.

    Returns:
        tuple: (expr_filtered, sample_map_filtered) aligned on the same
        aliquot barcodes, in the same column order.

    Raises:
        AssertionError: if columns are unmapped while require_all_mapped, or if
            no samples survive filtering.
    """
    if exclude_studies is None:
        exclude_studies = DEFAULT_EXCLUDED_STUDIES
    cols = pd.Index([str(c).strip() for c in expr.columns])
    expr = expr.copy()
    expr.columns = cols
    assert cols.is_unique, "duplicate aliquot barcodes in expression columns"

    mapped = cols.intersection(sample_map.index)
    unmapped = cols.difference(sample_map.index)
    if len(unmapped):
        assert not require_all_mapped, (
            "%d expression aliquots have no tumor-type assignment (e.g. %s) - "
            "extend the sample map or pass require_all_mapped=False"
            % (len(unmapped), list(unmapped[:5]))
        )
    sm = sample_map.loc[mapped]
    keep = sm.index[sm["is_primary_tumor"].astype(bool)
                    & ~sm["tumor_type"].isin(list(exclude_studies))]
    assert len(keep) > 0, "no primary solid tumor samples survived filtering"
    keep = [c for c in cols if c in set(keep)]
    return expr.loc[:, keep], sample_map.loc[keep]


def summarize_by_tumor_type(expr, sample_map, stat="median"):
    """Collapse a gene x sample matrix to gene x tumor_type summaries.

    Args:
        expr (pandas.DataFrame): genes x aliquot barcodes (numeric).
        sample_map (pandas.DataFrame): indexed by aliquot barcode with
            'tumor_type', 'tumor_name'.
        stat (str): one of 'median', 'mean', 'q25', 'q75', 'frac_pos'
            ('frac_pos' = fraction of samples with value > 0).

    Returns:
        tuple: (by_type, n_table)
          by_type (pandas.DataFrame): genes x tumor_type.
          n_table (pandas.DataFrame): one row per tumor_type with
            ['tumor_type', 'tumor_name', 'n_samples'].

    Raises:
        AssertionError: on unknown stat, or if expression columns are not all
            present in sample_map.
    """
    allowed = ("median", "mean", "q25", "q75", "frac_pos")
    assert stat in allowed, "stat must be one of %s, got %r" % (allowed, stat)
    cols = [str(c).strip() for c in expr.columns]
    missing = [c for c in cols if c not in sample_map.index]
    assert not missing, "%d expression columns absent from sample_map (e.g. %s)" % (
        len(missing), missing[:5])

    groups = sample_map.loc[cols, "tumor_type"].astype(str)
    num = expr.apply(pd.to_numeric, errors="coerce")
    num.columns = cols
    gb = num.T.groupby(groups.values)
    if stat == "median":
        by_type = gb.median().T
    elif stat == "mean":
        by_type = gb.mean().T
    elif stat == "q25":
        by_type = gb.quantile(0.25).T
    elif stat == "q75":
        by_type = gb.quantile(0.75).T
    else:
        by_type = gb.apply(lambda g: (g > 0).mean()).T
    by_type = by_type.sort_index(axis=1)

    n_table = (sample_map.loc[cols].groupby("tumor_type")
               .agg(tumor_name=("tumor_name", "first"), n_samples=("patient", "size"))
               .reset_index().sort_values("tumor_type").reset_index(drop=True))
    return by_type, n_table


# ---------------------------------------------------------------------------
# PanImmune / clinical loaders
# ---------------------------------------------------------------------------

def fetch_thorsson_immune(dest, force=False):
    """Download + load the PanCanAtlas sample annotation table (aliquot -> cancer type).

    This is merged_sample_quality_annotations.tsv, the PanImmune companion table
    published with Thorsson et al. It provides the aliquot barcode -> cancer type
    mapping and QC exclusion flags.

    Args:
        dest (str): local path to cache the TSV.
        force (bool): re-download if already present.

    Returns:
        pandas.DataFrame
    """
    download_to_disk(GDC_SAMPLE_ANNOTATION_URL, dest, force=force)
    return pd.read_csv(dest, sep="\t", low_memory=False)


def fetch_leukocyte_fraction(dest, force=False):
    """Download + load Thorsson et al. methylation-based leukocyte fraction.

    The file (TCGA_all_leuk_estimate.masked.20170107.tsv) is HEADERLESS with
    three columns: cancer type, aliquot barcode, leukocyte fraction.

    Args:
        dest (str): local cache path.
        force (bool): re-download if present.

    Returns:
        pandas.DataFrame with columns
        ['tumor_type', 'aliquot_barcode', 'leukocyte_fraction'].
    """
    download_to_disk(GDC_LEUKOCYTE_FRACTION_URL, dest, force=force)
    df = pd.read_csv(dest, sep="\t", header=None,
                     names=["tumor_type", "aliquot_barcode", "leukocyte_fraction"])
    df["leukocyte_fraction"] = pd.to_numeric(df["leukocyte_fraction"], errors="coerce")
    return df


def fetch_absolute_purity(dest, force=False):
    """Download + load ABSOLUTE tumor purity / ploidy master calls.

    Args:
        dest (str): local cache path.
        force (bool): re-download if present.

    Returns:
        pandas.DataFrame as published (columns include 'array', 'sample',
        'purity', 'ploidy').
    """
    download_to_disk(GDC_ABSOLUTE_PURITY_URL, dest, force=force)
    return pd.read_csv(dest, sep="\t", low_memory=False)


def fetch_signature_scores(dest, wanted_signatures=None, force=False):
    """Download + load PanImmune Scores_160_Signatures (samples in COLUMNS).

    The published layout is signatures x samples: columns 'Source', 'SetName',
    then one column per TCGA aliquot barcode. This function transposes to
    samples x signatures and optionally keeps only selected signatures, so the
    ~9100-column table need not be carried whole.

    Args:
        dest (str): local cache path for the .tsv.gz.
        wanted_signatures (iterable of str|None): SetName values to keep;
            None keeps all 160.
        force (bool): re-download if present.

    Returns:
        pandas.DataFrame indexed by aliquot barcode, one column per signature
        (column name = SetName).
    """
    download_to_disk(GDC_SIGNATURE_SCORES_URL, dest, force=force)
    df = pd.read_csv(dest, sep="\t", low_memory=False)
    assert "SetName" in df.columns, "unexpected signature table columns: %s" % list(
        df.columns)[:5]
    if wanted_signatures is not None:
        want = {str(s) for s in wanted_signatures}
        df = df[df["SetName"].astype(str).isin(want)]
        assert len(df), "none of the requested signatures were present"
    df = df.set_index(df["SetName"].astype(str))
    value_cols = [c for c in df.columns if str(c).startswith("TCGA-")]
    out = df[value_cols].T.apply(pd.to_numeric, errors="coerce")
    out.index.name = "aliquot_barcode"
    return out.loc[:, ~out.columns.duplicated()]


def fetch_tcga_cdr(dest, force=False, sheet=0):
    """Download + load the TCGA-CDR clinical outcome endpoints workbook.

    Args:
        dest (str): local cache path for the .xlsx.
        force (bool): re-download if present.
        sheet (int|str): sheet index or name to read.

    Returns:
        pandas.DataFrame (one row per patient; includes 'bcr_patient_barcode',
        'type', 'OS', 'OS.time', 'PFI', 'PFI.time', stage and demographics).
    """
    download_to_disk(GDC_TCGA_CDR_URL, dest, force=force)
    return pd.read_excel(dest, sheet_name=sheet)


def build_tme_table(leukocyte_df, purity_df=None, signature_df=None,
                   sample_map=None):
    """Assemble a per-sample tumor microenvironment table.

    Stromal fraction follows the Thorsson et al. definition
    `stromal_fraction = 1 - tumor_purity - leukocyte_fraction`, computed only
    where ABSOLUTE purity is available; it is clipped at 0 and NaN elsewhere.
    Purity is matched at SAMPLE depth (first 15 barcode characters) because
    ABSOLUTE calls come from DNA aliquots while expression comes from RNA
    aliquots of the same sample.

    Args:
        leukocyte_df (pandas.DataFrame): output of fetch_leukocyte_fraction.
        purity_df (pandas.DataFrame|None): output of fetch_absolute_purity.
        signature_df (pandas.DataFrame|None): output of fetch_signature_scores
            (samples x signatures, aliquot-indexed).
        sample_map (pandas.DataFrame|None): output of build_sample_map, used to
            attach tumor_type/tumor_name and restrict to mapped samples.

    Returns:
        pandas.DataFrame indexed by 'sample_barcode' (15-char TCGA sample id)
        with leukocyte_fraction, tumor_purity, ploidy, stromal_fraction, any
        signature columns, plus patient/tumor_type/tumor_name when available.
    """
    lf = leukocyte_df.copy()
    lf["sample_barcode"] = lf["aliquot_barcode"].astype(str).str[:15]
    lf_types = None
    if "tumor_type" in lf.columns:
        lf_types = (lf.drop_duplicates(subset=["sample_barcode"])
                      .set_index("sample_barcode")["tumor_type"]
                      .astype(str).str.strip().str.upper())
    lf = (lf.dropna(subset=["leukocyte_fraction"])
            .groupby("sample_barcode", as_index=True)["leukocyte_fraction"].median()
            .to_frame())

    out = lf
    if purity_df is not None:
        pu = purity_df.copy()
        key = "array" if "array" in pu.columns else pu.columns[0]
        pu["sample_barcode"] = pu[key].astype(str).str.strip().str[:15]
        keep = {}
        if "purity" in pu.columns:
            keep["tumor_purity"] = "purity"
        if "ploidy" in pu.columns:
            keep["ploidy"] = "ploidy"
        agg = pu.groupby("sample_barcode").agg(
            **{new: (old, "median") for new, old in keep.items()})
        out = out.join(agg, how="outer")

    if signature_df is not None:
        sg = signature_df.copy()
        sg.index = pd.Index([str(i)[:15] for i in sg.index], name="sample_barcode")
        sg = sg.groupby(level=0).median()
        out = out.join(sg, how="outer")

    if "tumor_purity" in out.columns:
        out["stromal_fraction"] = (
            1.0 - out["tumor_purity"] - out["leukocyte_fraction"]).clip(lower=0.0)

    if sample_map is not None:
        sm = sample_map.copy()
        sm["sample_barcode"] = [str(i)[:15] for i in sm.index]
        sm = (sm.drop_duplicates(subset=["sample_barcode"])
                .set_index("sample_barcode")[["patient", "tumor_type", "tumor_name",
                                              "sample_type_code", "is_primary_tumor"]])
        out = out.join(sm, how="left")
    # Fall back to the leukocyte table's own cancer-type column for samples the
    # annotation table does not cover, so tumor_type is missing only where no
    # published source states it.
    if lf_types is not None:
        if "tumor_type" not in out.columns:
            out["tumor_type"] = np.nan
        fill = out.index.to_series().map(lf_types)
        out["tumor_type"] = out["tumor_type"].fillna(fill)
        if "tumor_name" not in out.columns:
            out["tumor_name"] = np.nan
        out["tumor_name"] = out["tumor_name"].fillna(
            out["tumor_type"].map(TCGA_STUDY_NAMES))
    out.index.name = "sample_barcode"
    return out.sort_index()


def summarize_tme_by_tumor_type(tme_df, tumor_type_col="tumor_type", stat="median"):
    """Aggregate a per-sample TME table to per-tumor-type summaries.

    Args:
        tme_df (pandas.DataFrame): output of build_tme_table (must carry a
            tumor-type column).
        tumor_type_col (str): grouping column.
        stat (str): 'median' or 'mean'.

    Returns:
        pandas.DataFrame indexed by tumor_type: numeric column summaries plus
        'n_samples'.

    Raises:
        AssertionError: if the tumor-type column is absent or stat unknown.
    """
    assert tumor_type_col in tme_df.columns, (
        "no %r column - pass sample_map to build_tme_table" % tumor_type_col)
    assert stat in ("median", "mean"), "stat must be 'median' or 'mean'"
    num = tme_df.select_dtypes(include=[np.number])
    gb = num.groupby(tme_df[tumor_type_col].astype(str))
    agg = gb.median() if stat == "median" else gb.mean()
    agg["n_samples"] = tme_df.groupby(tme_df[tumor_type_col].astype(str)).size()
    agg.index.name = tumor_type_col
    return agg.sort_index()


# ---------------------------------------------------------------------------
# Gene set helpers
# ---------------------------------------------------------------------------

def parse_cellchatdb_genes(rda_path):
    """Extract ligand / receptor / cofactor gene symbols from CellChatDB.human.rda.

    Multi-subunit entries in the interaction table are names of rows in the
    `complex` table; cofactor entries index the `cofactor` table. Both are
    expanded to their constituent gene symbols.

    Requires the `rdata` package (pip install rdata).

    Args:
        rda_path (str): path to CellChatDB.human.rda.

    Returns:
        dict: {'receptors': sorted list[str], 'ligands': sorted list[str],
               'cofactors': sorted list[str], 'union': sorted list[str],
               'n_interactions': int, 'n_pathways': int}

    Raises:
        AssertionError: if no receptor genes could be parsed.
    """
    import rdata as _rdata

    conv = _rdata.conversion.convert(_rdata.parser.parse_file(rda_path))
    db = conv[[k for k in conv if "CellChatDB" in str(k)][0]]
    inter, cplx, cof = db["interaction"], db["complex"], db["cofactor"]

    def _tokens(val):
        s = str(val).strip()
        if s in ("", "nan", "None"):
            return []
        return [t.strip() for t in re.split(r"[,;]", s) if t.strip()]

    def _table_map(tbl):
        return {str(k): [str(v).strip() for v in row if str(v).strip() not in ("", "nan")]
                for k, row in zip(tbl.index, tbl.values)}

    cmap, comap = _table_map(cplx), _table_map(cof)

    def _expand(val, mapping):
        out = []
        for tok in _tokens(val):
            out.extend(mapping.get(tok, [tok]))
        return out

    receptors, ligands, cofactors = set(), set(), set()
    for _, row in inter.iterrows():
        ligands.update(_expand(row["ligand"], cmap))
        receptors.update(_expand(row["receptor"], cmap))
        for c in ("co_A_receptor", "co_I_receptor", "agonist", "antagonist"):
            if c in inter.columns:
                cofactors.update(_expand(row[c], comap))

    ok = lambda s: sorted({g for g in s if re.fullmatch(r"[A-Za-z0-9._-]+", g)})
    receptors, ligands, cofactors = ok(receptors), ok(ligands), ok(cofactors)
    assert receptors, "parsed 0 receptor genes from %s" % rda_path
    return {
        "receptors": receptors,
        "ligands": ligands,
        "cofactors": cofactors,
        "union": sorted(set(receptors) | set(ligands) | set(cofactors)),
        "n_interactions": int(len(inter)),
        "n_pathways": int(inter["pathway_name"].nunique())
                      if "pathway_name" in inter.columns else 0,
    }


def sample_background_genes(url_or_path, n=200, seed=0, exclude=None,
                           min_symbol_len=2):
    """Draw a reproducible random background gene set spanning the matrix.

    Streams only the gene_id column of the expression matrix (values are
    discarded), so this is memory-bounded regardless of matrix size. Intended
    as the null pool for matched-random specificity tests.

    Args:
        url_or_path (str): expression matrix URL or path.
        n (int): number of genes to draw.
        seed (int): numpy RNG seed - record it alongside the gene list.
        exclude (iterable of str): symbols to exclude (e.g. the target set).
        min_symbol_len (int): drop shorter symbols (filters '?' placeholders).

    Returns:
        dict: {'genes': sorted list[str], 'seed': int, 'n_requested': int,
               'n_universe': int}

    Raises:
        AssertionError: if fewer than n eligible symbols exist.
    """
    if exclude is None:
        exclude = None
    ex = {str(g).strip() for g in exclude}
    stream, closer = h_open_text_stream(url_or_path)
    universe = []
    try:
        stream.readline()
        for line in stream:
            tab = line.find("\t")
            if tab < 0:
                continue
            sym = h_strip(line[:tab]).split("|")[0]
            if len(sym) >= min_symbol_len and sym not in ex and re.fullmatch(
                    r"[A-Za-z0-9._-]+", sym):
                universe.append(sym)
    finally:
        try:
            stream.close()
        finally:
            closer()
    universe = sorted(set(universe))
    assert len(universe) >= n, "only %d eligible symbols, need %d" % (len(universe), n)
    rng = np.random.default_rng(seed)
    picks = rng.choice(np.array(universe, dtype=object), size=n, replace=False)
    return {"genes": sorted(str(g) for g in picks), "seed": int(seed),
            "n_requested": int(n), "n_universe": int(len(universe))}


def resolve_symbols_to_matrix(url_or_path, symbols, gene_info_path,
                             gene_id_sep="|"):
    """Map current HGNC symbols onto the matrix's frozen 2016 symbol vocabulary.

    The PanCanAtlas matrix was annotated around 2016, so genes renamed since
    then (CXCL8 was IL8, ACKR1 was DARC, NECTIN2 was PVRL2, ...) will not match
    by symbol and would be silently reported as absent. This resolves them
    exactly, by Entrez GeneID: the matrix's own `gene_id` field is
    "SYMBOL|ENTREZ", and NCBI `gene_info` gives the GeneID for each current
    symbol and each of its synonyms.

    Only the gene_id column of the matrix is streamed, so memory use is
    independent of matrix size.

    Args:
        url_or_path (str): expression matrix URL or path.
        symbols (iterable of str): current symbols to resolve.
        gene_info_path (str): path to NCBI Homo_sapiens.gene_info(.gz), from
            ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/.
        gene_id_sep (str): separator inside gene_id.

    Returns:
        dict: {
          'resolved': {current_symbol: matrix_symbol} for symbols reachable via
              Entrez ID or synonym but under a different name,
          'unresolvable': sorted list[str] genuinely absent from the matrix,
          'by_entrez': {current_symbol: matrix_symbol} subset matched on GeneID,
          'by_synonym': {current_symbol: matrix_symbol} subset matched on synonym only,
        }
    """
    want = {str(s).strip() for s in symbols if str(s).strip()}
    assert want, "symbols is empty"

    # symbol/synonym -> GeneID, from NCBI
    sym2id, syn2id = {}, {}
    opener = gzip.open if str(gene_info_path).endswith(".gz") else open
    with opener(gene_info_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            gid, sym, syns = f[1], f[2], f[4]
            sym2id[sym] = gid
            for s in syns.split("|"):
                if s and s != "-":
                    syn2id.setdefault(s, gid)

    # matrix vocabulary: entrez -> symbol, and symbol set
    stream, closer = h_open_text_stream(url_or_path)
    mat_by_entrez, mat_syms = {}, set()
    try:
        stream.readline()
        for line in stream:
            tab = line.find("\t")
            if tab < 0:
                continue
            gene_id = h_strip(line[:tab])
            parts = gene_id.split(gene_id_sep)
            mat_syms.add(parts[0])
            if len(parts) > 1 and parts[1].isdigit():
                mat_by_entrez.setdefault(parts[1], parts[0])
    finally:
        try:
            stream.close()
        finally:
            closer()

    resolved, by_entrez, by_syn, unresolvable = {}, {}, {}, []
    for s in sorted(want):
        if s in mat_syms:
            continue  # already matches directly
        gid = sym2id.get(s) or syn2id.get(s)
        if gid and gid in mat_by_entrez:
            resolved[s] = mat_by_entrez[gid]
            by_entrez[s] = mat_by_entrez[gid]
            continue
        hit = None
        for syn, sid in syn2id.items():
            if sid == gid and syn in mat_syms:
                hit = syn
                break
        if hit is None and gid:
            for msym, mid in ((m, sym2id.get(m)) for m in mat_syms):
                if mid == gid:
                    hit = msym
                    break
        if hit:
            resolved[s] = hit
            by_syn[s] = hit
        else:
            unresolvable.append(s)
    return {"resolved": resolved, "unresolvable": sorted(unresolvable),
            "by_entrez": by_entrez, "by_synonym": by_syn}


def stream_gene_medians(url_or_path, keep_samples=None, gene_id_sep="|",
                       min_symbol_len=2, progress_every=0, verbose=False):
    """Compute per-gene median expression across the whole matrix, one row at a time.

    Streams the matrix and reduces each row to a single median immediately, so
    peak memory is O(n_genes) floats regardless of matrix size. The resulting
    transcriptome-wide median distribution is what lets an ABSOLUTE expression
    floor be placed on a percentile basis, and what lets a background pool be
    stratified across the real dynamic range.

    Args:
        url_or_path (str): expression matrix URL or path.
        keep_samples (iterable of str|None): restrict to these aliquot barcodes
            (e.g. primary solid tumors only). None uses all columns.
        gene_id_sep (str): separator inside gene_id.
        min_symbol_len (int): drop shorter symbols (filters '?' placeholders).
        progress_every (int): progress marker interval (0 disables).
        verbose (bool): print progress to stderr.

    Returns:
        pandas.Series: gene symbol -> median value, sorted descending. Duplicate
        symbols keep the first occurrence.
    """
    stream, closer = h_open_text_stream(url_or_path)
    out = {}
    n = 0
    try:
        header_line = stream.readline()
        assert header_line, "expression stream is empty"
        header = [h_strip(t) for t in header_line.rstrip("\r\n").split("\t")]
        samples = header[1:]
        if keep_samples is None:
            idx = None
        else:
            want = {str(s).strip() for s in keep_samples}
            idx = np.array([i for i, s in enumerate(samples) if s in want], dtype=int)
            assert idx.size > 0, "keep_samples matched no columns in the header"
        for line in stream:
            n += 1
            if progress_every and verbose and n % progress_every == 0:
                sys.stderr.write("  medians: %d rows\n" % n)
                sys.stderr.flush()
            tab = line.find("\t")
            if tab < 0:
                continue
            sym = h_strip(line[:tab]).split(gene_id_sep)[0]
            if len(sym) < min_symbol_len or sym in out:
                continue
            vals = np.array(line[tab + 1:].rstrip("\r\n").split("\t"), dtype=object)
            if idx is not None:
                if vals.size != len(samples):
                    continue
                vals = vals[idx]
            arr = pd.to_numeric(pd.Series(vals), errors="coerce").to_numpy(dtype=float)
            out[sym] = float(np.nanmedian(arr)) if np.isfinite(arr).any() else np.nan
    finally:
        try:
            stream.close()
        finally:
            closer()
    assert out, "no gene rows parsed from %s" % url_or_path
    return pd.Series(out, name="median").sort_values(ascending=False)


def sample_background_genes_stratified(gene_medians, n=200, seed=0, exclude=(),
                                      n_strata=10):
    """Draw a seeded background gene pool stratified across expression strata.

    A uniform random draw over the transcriptome is dominated by low- and
    mid-expressed genes and yields almost no candidates at housekeeping level,
    so a matched-random specificity null has nothing to match a highly expressed
    receptor against. Stratifying by expression decile guarantees candidates at
    every level of the real dynamic range.

    Args:
        gene_medians (pandas.Series): gene -> median expression, e.g. from
            stream_gene_medians.
        n (int): total genes to draw (split as evenly as possible across strata).
        seed (int): numpy RNG seed - record it with the gene list.
        exclude (iterable of str): symbols to exclude (e.g. the target panel).
        n_strata (int): number of equal-count expression strata.

    Returns:
        dict: {'genes': sorted list[str], 'seed': int, 'n_requested': int,
               'n_strata': int, 'n_universe': int,
               'strata': [{'stratum': int, 'median_low': float,
                           'median_high': float, 'n_drawn': int,
                           'n_available': int}, ...]}

    Raises:
        AssertionError: if the eligible universe is smaller than n.
    """
    if exclude is None:
        exclude = ()
    ex = {str(g).strip() for g in exclude}
    s = gene_medians.dropna()
    s = s[~s.index.isin(ex)]
    assert len(s) >= n, "only %d eligible genes, need %d" % (len(s), n)

    rng = np.random.default_rng(seed)
    order = s.sort_values().index.to_numpy()
    parts = np.array_split(order, n_strata)
    base, extra = divmod(n, n_strata)
    picks, strat_info = [], []
    for i, part in enumerate(parts):
        take = min(base + (1 if i < extra else 0), len(part))
        chosen = rng.choice(part, size=take, replace=False) if take else np.array([])
        picks.extend(str(g) for g in chosen)
        strat_info.append({
            "stratum": i,
            "median_low": float(s[part].min()),
            "median_high": float(s[part].max()),
            "n_drawn": int(take),
            "n_available": int(len(part)),
        })
    return {"genes": sorted(set(picks)), "seed": int(seed), "n_requested": int(n),
            "n_strata": int(n_strata), "n_universe": int(len(s)),
            "strata": strat_info}


def housekeeping_medians(expr, genes=None, sample_map=None):
    """Median housekeeping-gene expression, overall and optionally per tumor type.

    These values calibrate an ABSOLUTE expression floor: a receptor whose value
    sits orders of magnitude below the housekeeping medians in the same unit is
    effectively untranscribed, which a within-gene z-score cannot express.

    Args:
        expr (pandas.DataFrame): genes x samples, in RSEM normalized counts.
        genes (iterable of str): housekeeping symbols.
        sample_map (pandas.DataFrame|None): if given, also compute per-tumor-type
            medians.

    Returns:
        dict: {'overall': {gene: float}, 'by_tumor_type': {tumor_type:
        {gene: float}} or {}, 'genes_absent': list[str]}
    """
    if genes is None:
        genes = HOUSEKEEPING_GENES
    present = [g for g in genes if g in expr.index]
    absent = [g for g in genes if g not in expr.index]
    sub = expr.loc[present].apply(pd.to_numeric, errors="coerce")
    overall = {g: float(v) for g, v in sub.median(axis=1).items()}
    by_type = {}
    if sample_map is not None and len(present):
        bt, _ = summarize_by_tumor_type(sub, sample_map, stat="median")
        by_type = {str(t): {g: float(bt.loc[g, t]) for g in bt.index}
                   for t in bt.columns}
    return {"overall": overall, "by_tumor_type": by_type, "genes_absent": absent}
