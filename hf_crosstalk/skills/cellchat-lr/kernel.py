"""CellChatDB human ligand-receptor database: fetch, parse, expand, map.

Downloads CellChatDB.human.rda from the CellChat GitHub repository, parses the
R serialization in pure Python, and exposes the interaction / complex / cofactor
/ geneInfo tables plus complex-aware receptor and ligand subunit expansion.

Pure stdlib + pandas + rdata/pyreadr. No platform dependencies: every function
takes file paths or DataFrames and returns DataFrames / dicts.
"""

import os
import urllib.request

import pandas as pd

CELLCHATDB_HUMAN_URL = (
    "https://raw.githubusercontent.com/sqjin/CellChat/master/data/CellChatDB.human.rda"
)
CELLCHATDB_MOUSE_URL = (
    "https://raw.githubusercontent.com/sqjin/CellChat/master/data/CellChatDB.mouse.rda"
)

# Elements of the CellChatDB list object.
DB_ELEMENTS = ("interaction", "complex", "cofactor", "geneInfo")

# Columns the interaction table must carry for downstream use.
REQUIRED_INTERACTION_COLUMNS = (
    "interaction_name",
    "pathway_name",
    "ligand",
    "receptor",
    "annotation",
    "evidence",
)

# CellChatDB annotation categories.
ANNOTATION_SECRETED = "Secreted Signaling"
ANNOTATION_ECM = "ECM-Receptor"
ANNOTATION_CONTACT = "Cell-Cell Contact"

# Minimum interaction row count; the human DB ships ~1900 interactions.
MIN_INTERACTION_ROWS = 1000


def fetch_cellchatdb(dest_path, url=None, force=False, timeout=300):
    """Download (or reuse) the CellChatDB .rda file.

    Args:
        dest_path (str): local path for the .rda cache. Reused if it already
            exists and is non-empty, unless force is True.
        url (str): source URL. Defaults to the human database.
        force (bool): re-download even if dest_path exists.
        timeout (int): socket timeout in seconds.

    Returns:
        str: dest_path.

    Raises:
        AssertionError: if the downloaded file is implausibly small (<100 KB) or
            is not a compressed R serialization (bzip2/gzip/xz magic bytes),
            which would indicate an HTML error page was saved instead.
    """
    if url is None:
        url = CELLCHATDB_HUMAN_URL
    if force or not os.path.exists(dest_path) or os.path.getsize(dest_path) == 0:
        parent = os.path.dirname(os.path.abspath(dest_path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            payload = resp.read()
        with open(dest_path, "wb") as fh:
            fh.write(payload)

    size = os.path.getsize(dest_path)
    assert size > 100_000, (
        "CellChatDB download at %s is only %d bytes; expected >100 KB - the "
        "fetch probably returned an error page" % (dest_path, size)
    )
    with open(dest_path, "rb") as fh:
        magic = fh.read(6)
    assert (
        magic[:3] == b"BZh" or magic[:2] == b"\x1f\x8b" or magic[:6] == b"\xfd7zXZ\x00"
        or magic[:5] == b"RDX2\n" or magic[:5] == b"RDX3\n"
    ), (
        "file at %s does not look like an R .rda (magic=%r); refusing to parse"
        % (dest_path, magic)
    )
    return dest_path


def h_clean_str_df(df):
    """Coerce every cell to a stripped python str, empty string for missing."""
    out = df.copy()
    out.columns = [str(c) for c in out.columns]
    for c in out.columns:
        out[c] = out[c].apply(
            lambda v: "" if v is None or (isinstance(v, float) and v != v) else str(v).strip()
        )
    return out


def load_cellchatdb(path):
    """Parse a CellChatDB .rda into tidy DataFrames.

    Parsed with the `rdata` package (pure-Python R serialization reader). The
    top-level R object is a named list; the four elements are returned as
    DataFrames with their R row names promoted to an explicit column:
    'interaction_name' already exists on the interaction table, while the
    complex and cofactor tables carry their key only as row names, so those are
    promoted to 'complex_name' / 'cofactor_name'.

    Args:
        path (str): path to the .rda file (see fetch_cellchatdb).

    Returns:
        dict: {'interaction': DataFrame, 'complex': DataFrame,
               'cofactor': DataFrame, 'geneInfo': DataFrame}. All values are
        strings ('' for missing). The interaction table has one row per
        ligand-receptor interaction with columns including interaction_name,
        pathway_name, ligand, receptor, agonist, antagonist, co_A_receptor,
        co_I_receptor, evidence, annotation, interaction_name_2. The complex
        table has complex_name + subunit_1..subunit_n. The cofactor table has
        cofactor_name + cofactor1..cofactorN. geneInfo has Symbol, Name,
        EntrezGene.ID, Ensembl.Gene.ID, MGI.ID, Gene.group.name.

    Raises:
        AssertionError: if the object is not a list containing the four expected
            elements, if the interaction table has fewer than
            MIN_INTERACTION_ROWS rows, or if any of
            REQUIRED_INTERACTION_COLUMNS is missing.
    """
    import rdata

    parsed = rdata.parser.parse_file(path)
    converted = rdata.conversion.convert(parsed)
    assert isinstance(converted, dict) and len(converted) >= 1, (
        "parsed .rda is not a named-object mapping: %r" % type(converted)
    )
    # The .rda contains exactly one top-level object (CellChatDB.human/mouse).
    obj_name = [k for k in converted.keys()][0]
    db = converted[obj_name]
    assert isinstance(db, dict), (
        "%s is not an R list (got %r); expected a list with elements %s"
        % (obj_name, type(db), list(DB_ELEMENTS))
    )
    keys = {str(k): k for k in db.keys()}
    missing = [e for e in DB_ELEMENTS if e not in keys]
    assert not missing, "CellChatDB object is missing element(s) %s; has %s" % (
        missing,
        list(keys),
    )

    out = {}
    for elem in DB_ELEMENTS:
        raw = db[keys[elem]]
        assert isinstance(raw, pd.DataFrame), (
            "element %r is %r, expected a data.frame" % (elem, type(raw))
        )
        df = h_clean_str_df(raw)
        idx = pd.Index([str(i) for i in raw.index])
        if elem == "complex":
            df.insert(0, "complex_name", idx)
        elif elem == "cofactor":
            df.insert(0, "cofactor_name", idx)
        elif elem == "geneInfo":
            df.insert(0, "hgnc_id", idx)
        elif elem == "interaction" and "interaction_name" not in df.columns:
            df.insert(0, "interaction_name", idx)
        out[elem] = df.reset_index(drop=True)

    inter = out["interaction"]
    assert len(inter) >= MIN_INTERACTION_ROWS, (
        "interaction table has only %d rows (<%d); the parse or the download is "
        "incomplete" % (len(inter), MIN_INTERACTION_ROWS)
    )
    miss = [c for c in REQUIRED_INTERACTION_COLUMNS if c not in inter.columns]
    assert not miss, "interaction table is missing required column(s) %s; has %s" % (
        miss,
        list(inter.columns),
    )
    return out


def tidy_complexes(complex_df):
    """Melt the wide complex table into one row per (complex, subunit).

    Args:
        complex_df (pandas.DataFrame): the 'complex' element from
            load_cellchatdb (complex_name + subunit_1..subunit_n).

    Returns:
        pandas.DataFrame: columns complex_name, subunit_index (int, 1-based),
            subunit_gene (str), n_subunits (int, total non-empty subunits of
            that complex). Empty subunit slots are dropped.

    Raises:
        AssertionError: if no subunit columns are found or any retained
            subunit_gene is empty.
    """
    assert "complex_name" in complex_df.columns, "complex_df needs 'complex_name'"
    sub_cols = [c for c in complex_df.columns if c.startswith("subunit")]
    assert sub_cols, "no subunit_* columns in complex_df: %s" % list(complex_df.columns)

    long = complex_df.melt(
        id_vars=["complex_name"],
        value_vars=sub_cols,
        var_name="subunit_col",
        value_name="subunit_gene",
    )
    long["subunit_gene"] = long["subunit_gene"].fillna("").str.strip()
    long = long[long["subunit_gene"] != ""].copy()
    long["subunit_index"] = (
        long["subunit_col"].str.extract(r"(\d+)$")[0].astype(int)
    )
    long = long.sort_values(["complex_name", "subunit_index"])
    long["n_subunits"] = long.groupby("complex_name")["subunit_gene"].transform("size")
    out = long[
        ["complex_name", "subunit_index", "subunit_gene", "n_subunits"]
    ].reset_index(drop=True)
    assert (out["subunit_gene"].str.len() > 0).all(), "empty subunit survived melt"
    return out


def h_expand_one(names, cmap):
    """Expand a Series of ligand/receptor names into (index, subunit) tuples."""
    recs = []
    for pos, name in enumerate(names):
        nm = "" if name is None else str(name).strip()
        subs = cmap.get(nm)
        if subs:
            for i, g in enumerate(subs, start=1):
                recs.append((pos, nm, i, g, len(subs), True))
        else:
            recs.append((pos, nm, 1, nm, 1, False))
    return recs


def expand_receptor_complexes(interaction_df, complex_df, side="receptor"):
    """Expand multi-subunit receptors (or ligands) to one row per subunit gene.

    A receptor/ligand name that appears in the complex table is expanded to its
    constituent subunit genes; a name that does not appear is treated as its own
    single subunit (is_complex=False, n_subunits=1).

    Args:
        interaction_df (pandas.DataFrame): the 'interaction' element from
            load_cellchatdb.
        complex_df (pandas.DataFrame): the 'complex' element (wide form) or the
            output of tidy_complexes (long form); both are accepted.
        side (str): 'receptor' or 'ligand' - which column to expand.

    Returns:
        pandas.DataFrame: all original interaction columns plus
            <side>_subunit_index (int, 1-based), <side>_subunit_gene (str),
            n_subunits (int), is_complex (bool). One row per
            (interaction, subunit).

    Raises:
        AssertionError: if side is not 'receptor'/'ligand', if the column is
            absent, or if the expansion loses any interaction.
    """
    assert side in ("receptor", "ligand"), "side must be 'receptor' or 'ligand'"
    assert side in interaction_df.columns, (
        "interaction_df has no %r column" % side
    )

    if "subunit_gene" in complex_df.columns:
        tidy = complex_df
    else:
        tidy = tidy_complexes(complex_df)
    cmap = (
        tidy.sort_values("subunit_index")
        .groupby("complex_name")["subunit_gene"]
        .apply(list)
        .to_dict()
    )

    recs = h_expand_one(interaction_df[side].tolist(), cmap)
    positions = [r[0] for r in recs]
    out = interaction_df.iloc[positions].reset_index(drop=True)
    out["%s_subunit_index" % side] = [r[2] for r in recs]
    out["%s_subunit_gene" % side] = [r[3] for r in recs]
    out["n_subunits"] = [r[4] for r in recs]
    out["is_complex"] = [r[5] for r in recs]

    assert out["interaction_name"].nunique() == interaction_df[
        "interaction_name"
    ].nunique(), "expansion dropped interactions"
    assert (out["%s_subunit_gene" % side].str.len() > 0).all(), (
        "expansion produced an empty subunit gene"
    )
    return out


def expand_all_subunits(interaction_df, complex_df):
    """Expand BOTH ligand and receptor sides to one row per subunit pair.

    Args:
        interaction_df (pandas.DataFrame): 'interaction' element.
        complex_df (pandas.DataFrame): 'complex' element (wide or tidy).

    Returns:
        pandas.DataFrame: one row per (interaction, ligand_subunit,
            receptor_subunit) with columns interaction_name, pathway_name,
            ligand, receptor, annotation, evidence, plus
            ligand_subunit_index, ligand_subunit_gene, n_ligand_subunits,
            ligand_is_complex, receptor_subunit_index, receptor_subunit_gene,
            n_receptor_subunits, receptor_is_complex, and any other original
            interaction columns.
    """
    lig = expand_receptor_complexes(interaction_df, complex_df, side="ligand")
    lig = lig.rename(
        columns={"n_subunits": "n_ligand_subunits", "is_complex": "ligand_is_complex"}
    )
    both = expand_receptor_complexes(lig, complex_df, side="receptor")
    both = both.rename(
        columns={
            "n_subunits": "n_receptor_subunits",
            "is_complex": "receptor_is_complex",
        }
    )
    return both.reset_index(drop=True)


def map_ligands_to_receptors(
    ligand_genes,
    db,
    annotation="Secreted Signaling",
    expand_ligand_complexes=True,
):
    """Find interactions whose ligand matches an input gene list.

    A ligand matches if the interaction's ligand name equals a gene symbol in
    ligand_genes, or (when expand_ligand_complexes is True) if any subunit of a
    multi-subunit ligand complex is in ligand_genes.

    The annotation filter defaults to 'Secreted Signaling' because only secreted
    signalling is mechanistically plausible for a factor acting at a distance
    from its tissue of origin (e.g. heart -> distant tumour). 'ECM-Receptor'
    interactions require the ligand to be deposited in the local matrix and
    'Cell-Cell Contact' requires physical juxtaposition, so neither supports a
    circulating-factor mechanism. Both remain retrievable: pass
    annotation=None for everything, or a string / list of category names.

    Args:
        ligand_genes (iterable of str): candidate secreted ligand gene symbols.
        db (dict): output of load_cellchatdb (needs 'interaction' and 'complex').
        annotation (str | list[str] | None): CellChatDB annotation category or
            categories to keep. None keeps all.
        expand_ligand_complexes (bool): also match ligand-complex subunits.

    Returns:
        pandas.DataFrame: one row per (matched interaction, receptor subunit),
            columns: interaction_name, pathway_name, ligand, receptor,
            annotation, evidence, matched_ligand_gene, ligand_is_complex,
            n_ligand_subunits, receptor_subunit_gene, receptor_subunit_index,
            n_receptor_subunits, receptor_is_complex, plus agonist/antagonist/
            co_A_receptor/co_I_receptor when present. Empty DataFrame with the
            same columns if nothing matches.

    Raises:
        AssertionError: if db lacks the required elements, or if annotation
            names a category absent from the database.
    """
    for k in ("interaction", "complex"):
        assert k in db, "db is missing element %r" % k
    inter = db["interaction"]
    cplx = db["complex"]

    if annotation is not None:
        wanted = [annotation] if isinstance(annotation, str) else list(annotation)
        present = set(inter["annotation"].unique())
        bad = [a for a in wanted if a not in present]
        assert not bad, "annotation %s not in database; available: %s" % (
            bad,
            sorted(present),
        )
        inter = inter[inter["annotation"].isin(wanted)].copy()

    genes = {str(g).strip() for g in ligand_genes if str(g).strip()}

    lig = expand_receptor_complexes(inter, cplx, side="ligand") if len(inter) else inter
    if len(inter) == 0:
        return pd.DataFrame(
            columns=list(inter.columns)
            + [
                "matched_ligand_gene",
                "ligand_is_complex",
                "n_ligand_subunits",
                "receptor_subunit_gene",
                "receptor_subunit_index",
                "n_receptor_subunits",
                "receptor_is_complex",
            ]
        )

    if expand_ligand_complexes:
        hit = lig["ligand_subunit_gene"].isin(genes)
    else:
        hit = lig["ligand"].isin(genes)
    lig = lig[hit].copy()
    lig = lig.rename(
        columns={
            "ligand_subunit_gene": "matched_ligand_gene",
            "n_subunits": "n_ligand_subunits",
            "is_complex": "ligand_is_complex",
        }
    ).drop(columns=["ligand_subunit_index"], errors="ignore")

    if len(lig) == 0:
        return lig.assign(
            receptor_subunit_gene=pd.Series(dtype=str),
            receptor_subunit_index=pd.Series(dtype=int),
            n_receptor_subunits=pd.Series(dtype=int),
            receptor_is_complex=pd.Series(dtype=bool),
        ).reset_index(drop=True)

    out = expand_receptor_complexes(lig, cplx, side="receptor")
    out = out.rename(
        columns={
            "n_subunits": "n_receptor_subunits",
            "is_complex": "receptor_is_complex",
        }
    )
    front = [
        c
        for c in (
            "interaction_name",
            "pathway_name",
            "ligand",
            "matched_ligand_gene",
            "ligand_is_complex",
            "n_ligand_subunits",
            "receptor",
            "receptor_subunit_gene",
            "receptor_subunit_index",
            "n_receptor_subunits",
            "receptor_is_complex",
            "annotation",
            "evidence",
        )
        if c in out.columns
    ]
    rest = [c for c in out.columns if c not in front]
    return out[front + rest].reset_index(drop=True)


def summarize_db(db):
    """Compute descriptive counts for a parsed CellChatDB.

    Args:
        db (dict): output of load_cellchatdb.

    Returns:
        dict: n_interactions, n_pathways, n_ligands, n_receptors, n_complexes,
            n_multisubunit_receptor_interactions, n_multisubunit_receptor_names,
            n_cofactors, n_geneinfo_symbols, by_annotation (dict),
            receptor_subunit_rows.
    """
    inter = db["interaction"]
    cplx = db["complex"]
    tidy = tidy_complexes(cplx)
    exp = expand_receptor_complexes(inter, tidy, side="receptor")
    return {
        "n_interactions": int(len(inter)),
        "n_pathways": int(inter["pathway_name"].nunique()),
        "n_ligands": int(inter["ligand"].nunique()),
        "n_receptors": int(inter["receptor"].nunique()),
        "n_complexes": int(tidy["complex_name"].nunique()),
        "n_multisubunit_receptor_interactions": int(
            exp.loc[exp["n_subunits"] > 1, "interaction_name"].nunique()
        ),
        "n_multisubunit_receptor_names": int(
            exp.loc[exp["n_subunits"] > 1, "receptor"].nunique()
        ),
        "n_cofactors": int(len(db["cofactor"])),
        "n_geneinfo_symbols": int(db["geneInfo"]["Symbol"].nunique()),
        "by_annotation": {
            str(k): int(v) for k, v in inter["annotation"].value_counts().items()
        },
        "receptor_subunit_rows": int(len(exp)),
    }
