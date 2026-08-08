"""Two-factor ligand-receptor crosstalk scoring with mandatory validators.

Scores how compatible each receiver context (e.g. tumor type) is with a set of
ligands released by a source tissue, decomposed into:

    availability  - how strongly the source releases the ligand
    capacity      - how well the receiver expresses the receptor complex

The module's reason for existing is the *validators*. A relative score
(z-scored across receivers) will happily rank a receptor highly when that
receptor is not transcribed anywhere, because z-scoring a noise-level gene
still yields large z-values. Two guards are therefore enforced in code, not
prose:

    require_expression_floor()  - raises if a scored receptor never clears an
                                  absolute floor in any receiver.
    matched_random_null()       - permutation null over expression-matched
                                  receptor sets; a score that a random matched
                                  set reproduces is not receptor-specific.

Portable: pure stdlib + numpy/pandas/scipy. No platform APIs, no network.
"""

import json
import numpy as np
import pandas as pd
from scipy import stats

# Receptor complex aggregation rules. "min" is the default because a
# heteromeric receptor cannot signal if any obligate subunit is absent.
COMPLEX_RULES = ("min", "geomean", "mean")

DEFAULT_FLOOR_FRACTION = 0.005
DEFAULT_FLOOR_REFERENCE = "GAPDH"


def housekeeping_floor(expr, reference=None,
                       fraction=None):
    """Absolute expression floor expressed relative to a housekeeping gene.

    Using a housekeeping reference makes the floor unit-agnostic: it works for
    RSEM normalized counts, TPM or FPKM without hardcoding a value that is
    only meaningful in one unit.

    Args:
        expr: DataFrame, genes (index) x receivers (columns).
        reference: housekeeping gene symbol present in expr.index.
        fraction: floor as a fraction of the reference's median level.

    Returns:
        float: the absolute floor in the same unit as expr.

    Raises:
        KeyError: reference gene absent from expr.
    """
    if reference is None:
        reference = DEFAULT_FLOOR_REFERENCE
    if fraction is None:
        fraction = DEFAULT_FLOOR_FRACTION
    if reference not in expr.index:
        raise KeyError(
            "floor reference gene %r absent from expression matrix; pass a "
            "reference that was extracted (e.g. one of ACTB/GAPDH/RPLP0)"
            % reference)
    ref_median = float(np.nanmedian(expr.loc[reference].astype(float)))
    if not np.isfinite(ref_median) or ref_median <= 0:
        raise ValueError("reference gene %r has non-positive median %r"
                         % (reference, ref_median))
    return ref_median * float(fraction)


def receptor_capacity(expr, interactions, rule="min",
                      receptor_col="receptor_subunit_gene",
                      key_col="interaction_name"):
    """Aggregate subunit expression into a per-interaction receptor capacity.

    Args:
        expr: DataFrame genes x receivers (absolute units, non-negative).
        interactions: DataFrame with one row per (interaction, subunit).
        rule: one of COMPLEX_RULES. "min" = obligate-subunit semantics.
        receptor_col: column holding the subunit gene symbol.
        key_col: column identifying the interaction.

    Returns:
        DataFrame interactions (index) x receivers, the aggregated receptor
        level in the SAME absolute unit as expr (no standardization yet).

    Raises:
        ValueError: unknown rule, or no subunit found in expr.
    """
    if rule not in COMPLEX_RULES:
        raise ValueError("rule must be one of %r, got %r" % (COMPLEX_RULES, rule))
    rows = {}
    dropped = {}
    for key, grp in interactions.groupby(key_col):
        genes = [g for g in grp[receptor_col].astype(str).unique() if g in expr.index]
        missing = [g for g in grp[receptor_col].astype(str).unique()
                   if g not in expr.index]
        if not genes:
            dropped[key] = missing
            continue
        block = expr.loc[genes].astype(float)
        if rule == "min":
            agg = block.min(axis=0)
        elif rule == "mean":
            agg = block.mean(axis=0)
        else:
            agg = np.exp(np.log(block.clip(lower=1e-9)).mean(axis=0))
        rows[key] = agg
        if missing:
            dropped.setdefault(key, missing)
    if not rows:
        raise ValueError("no interaction had any receptor subunit present in "
                         "the expression matrix")
    out = pd.DataFrame(rows).T
    out.attrs["dropped_subunits"] = dropped
    out.attrs["rule"] = rule
    return out


def require_expression_floor(capacity, floor, raise_on_violation=True):
    """MANDATORY VALIDATOR. Reject receptors never expressed above `floor`.

    A receptor complex whose capacity is below `floor` in *every* receiver is
    not transcribed anywhere in the cohort and cannot mediate signalling,
    however high it scores after standardization. This is the guard whose
    absence makes relative co-expression rankings unsafe.

    Args:
        capacity: DataFrame interactions x receivers (absolute units).
        floor: absolute threshold, e.g. from housekeeping_floor().
        raise_on_violation: if True raise when violations exist; if False
            return them for the caller to drop.

    Returns:
        (passing_capacity, audit) where audit is a DataFrame with columns
        interaction, max_capacity, floor, passes.

    Raises:
        AssertionError: violations present and raise_on_violation=True.
    """
    mx = capacity.max(axis=1).astype(float)
    audit = pd.DataFrame({
        "interaction": mx.index,
        "max_capacity": mx.values,
        "floor": float(floor),
    })
    audit["passes"] = audit["max_capacity"] >= float(floor)
    audit = audit.sort_values("max_capacity").reset_index(drop=True)
    bad = audit.loc[~audit["passes"], "interaction"].tolist()
    if bad and raise_on_violation:
        raise AssertionError(
            "expression floor violated by %d interaction(s): %s. These "
            "receptors are below %.4g in every receiver and must not be "
            "scored. Pass raise_on_violation=False to drop them instead."
            % (len(bad), bad[:12], float(floor)))
    return capacity.loc[audit.loc[audit["passes"], "interaction"]], audit


def ligand_availability(de_table, ligand_genes, gene_col="gene_symbol",
                        effect_col="log2FC_discovery",
                        effect_col2="log2FC_validation",
                        tier_map=None, tier_weights=None):
    """Per-ligand availability from source-tissue differential expression.

    Availability = mean replicated effect size, optionally scaled by a
    secretion-tier weight (a plasma-detectable secreted protein is a more
    plausible long-range signal than a merely predicted-secreted one).

    Args:
        de_table: DataFrame of source DE results.
        ligand_genes: iterable of ligand gene symbols to score.
        gene_col, effect_col, effect_col2: column names. If effect_col2 is
            present both cohorts are averaged.
        tier_map: optional dict gene -> tier string.
        tier_weights: optional dict tier -> multiplier (default core=1.0,
            extended=0.5).

    Returns:
        Series indexed by ligand gene, availability >= 0 (non-positive effects
        clipped to 0 since only upregulation is a plausible HF-driven signal).
    """
    if tier_weights is None:
        tier_weights = {"core": 1.0, "extended": 0.5, "none": 0.25}
    d = de_table.drop_duplicates(subset=[gene_col]).set_index(gene_col)
    keep = [g for g in ligand_genes if g in d.index]
    eff = d.loc[keep, effect_col].astype(float)
    if effect_col2 and effect_col2 in d.columns:
        eff = (eff + d.loc[keep, effect_col2].astype(float)) / 2.0
    eff = eff.clip(lower=0.0)
    if tier_map:
        w = pd.Series({g: tier_weights.get(str(tier_map.get(g, "none")), 0.25)
                       for g in keep})
        eff = eff * w
    eff.name = "availability"
    return eff


def double_center(mat):
    """Remove additive row and column effects (confound adjustment).

    A per-receiver global expression axis (RNA content, normalization,
    cellularity) inflates every gene in some receivers. Double-centering
    removes additive row and column means so what remains is receptor-by-
    receiver *specificity*.

    Args:
        mat: DataFrame, numeric.

    Returns:
        DataFrame, same shape, grand-mean-restored double-centered values.
    """
    m = mat.astype(float)
    return m.sub(m.mean(axis=1), axis=0).sub(m.mean(axis=0), axis=1) + m.values.mean()


def crosstalk_score(expr, interactions, availability, rule="min",
                    floor=None, floor_reference=None,
                    floor_fraction=None,
                    adjust="none", enforce_floor=True,
                    ligand_col="ligand_subunit_gene",
                    receptor_col="receptor_subunit_gene",
                    key_col="interaction_name"):
    """Score each receiver's compatibility with a source ligand set.

    Pipeline: aggregate receptor subunits -> enforce the absolute expression
    floor -> optionally double-center -> z-standardize each receptor across
    receivers -> average across interactions weighted by ligand availability.

    Args:
        expr: DataFrame genes x receivers, absolute units.
        interactions: long DataFrame, one row per (interaction, subunit), with
            ligand_col, receptor_col, key_col and optionally pathway_name.
        availability: Series ligand gene -> availability weight.
        rule: receptor complex aggregation rule.
        floor: absolute floor; if None derived from floor_reference/fraction.
        adjust: "none" or "double_center".
        enforce_floor: if True a below-floor receptor raises; if False it is
            dropped and recorded.

    Returns:
        dict with keys:
          score           Series receiver -> compatibility
          capacity        DataFrame interactions x receivers (absolute)
          zscore          DataFrame interactions x receivers (standardized)
          contributions   DataFrame per-interaction weighted contribution
          floor_audit     DataFrame from require_expression_floor
          floor           float used
          n_interactions  int scored
          dropped         list of interactions dropped by the floor
    """
    if floor_reference is None:
        floor_reference = DEFAULT_FLOOR_REFERENCE
    if floor_fraction is None:
        floor_fraction = DEFAULT_FLOOR_FRACTION
    if floor is None:
        floor = housekeeping_floor(expr, floor_reference, floor_fraction)
    cap = receptor_capacity(expr, interactions, rule=rule,
                            receptor_col=receptor_col, key_col=key_col)
    cap_ok, audit = require_expression_floor(cap, floor,
                                             raise_on_violation=enforce_floor)
    dropped = audit.loc[~audit["passes"], "interaction"].tolist()

    mat = double_center(cap_ok) if adjust == "double_center" else cap_ok
    # z-standardize each receptor across receivers
    sd = mat.std(axis=1).replace(0.0, np.nan)
    z = mat.sub(mat.mean(axis=1), axis=0).div(sd, axis=0).dropna(how="all")

    lig_of = (interactions.drop_duplicates(subset=[key_col])
              .set_index(key_col)[ligand_col].astype(str))
    w = pd.Series({k: float(availability.get(lig_of.get(k, ""), 0.0))
                   for k in z.index})
    w = w.reindex(z.index).fillna(0.0)
    if w.sum() <= 0:
        raise ValueError("all ligand availability weights are zero; check that "
                         "availability index matches the interaction ligands")
    contrib = z.mul(w, axis=0)
    score = contrib.sum(axis=0) / w.sum()
    return {
        "score": score.sort_values(ascending=False),
        "capacity": cap_ok,
        "zscore": z,
        "weights": w,
        "contributions": contrib,
        "floor_audit": audit,
        "floor": float(floor),
        "n_interactions": int(z.shape[0]),
        "dropped": dropped,
        "rule": rule,
        "adjust": adjust,
    }


def matched_random_null(expr, interactions, availability, background_genes,
                        n_iter=1000, seed=0, n_bins=10, tolerance=0.5,
                        **score_kwargs):
    """MANDATORY VALIDATOR. Is the ranking receptor-specific, or generic?

    Replaces each real receptor subunit with a random gene matched on mean
    expression decile, re-scores, and correlates the null ranking with the
    observed one. If expression-matched random receptors reproduce the
    ranking (high rho), the score reflects a global per-receiver expression
    axis rather than the nominated receptors.

    Args:
        expr, interactions, availability: as crosstalk_score.
        background_genes: candidate pool for substitution.
        n_iter: permutations.
        seed: RNG seed (recorded for reproducibility).
        n_bins: expression deciles used for matching.
        tolerance: unused placeholder for API stability.

    Returns:
        dict with observed_score, null_rho (list), rho_mean, rho_p95,
        frac_rho_above_0_5, and per-receiver empirical p-values for being
        ranked as high as observed.
    """
    rng = np.random.default_rng(seed)
    obs = crosstalk_score(expr, interactions, availability, **score_kwargs)
    obs_score = obs["score"]

    pool = [g for g in background_genes if g in expr.index]
    if len(pool) < 20:
        raise ValueError("background pool too small after intersecting with "
                         "the expression matrix: %d genes" % len(pool))
    mean_expr = expr.mean(axis=1).astype(float)
    # decile bins over the union of pool and real receptors
    bins = pd.qcut(mean_expr.rank(method="first"), n_bins, labels=False)
    pool_by_bin = {}
    for g in pool:
        pool_by_bin.setdefault(int(bins[g]), []).append(g)

    rhos, higher = [], pd.Series(0, index=obs_score.index, dtype=float)
    kw = dict(score_kwargs)
    kw["enforce_floor"] = False  # random draws may land below floor
    for _ in range(int(n_iter)):
        sub = {}
        for g in interactions["receptor_subunit_gene"].astype(str).unique():
            if g not in expr.index:
                continue
            b = int(bins[g])
            cands = pool_by_bin.get(b) or pool
            sub[g] = cands[int(rng.integers(len(cands)))]
        perm = interactions.copy()
        perm["receptor_subunit_gene"] = (perm["receptor_subunit_gene"]
                                        .astype(str).map(lambda x: sub.get(x, x)))
        try:
            null = crosstalk_score(expr, perm, availability, **kw)
        except (ValueError, AssertionError):
            continue
        common = obs_score.index.intersection(null["score"].index)
        if len(common) > 3:
            rhos.append(float(stats.spearmanr(obs_score[common],
                                              null["score"][common]).statistic))
        higher = higher.add(
            (null["score"].reindex(obs_score.index) >= obs_score).astype(float),
            fill_value=0.0)
    n_ok = max(len(rhos), 1)
    return {
        "observed_score": obs_score,
        "null_rho": rhos,
        "rho_mean": float(np.mean(rhos)) if rhos else float("nan"),
        "rho_p95": float(np.percentile(rhos, 95)) if rhos else float("nan"),
        "frac_rho_above_0_5": float(np.mean([r > 0.5 for r in rhos])) if rhos else float("nan"),
        "per_receiver_p": ((higher + 1) / (n_ok + 1)).sort_values(),
        "n_iter_ok": n_ok,
        "seed": int(seed),
    }


def bootstrap_ci(expr, interactions, availability, n_boot=1000, seed=0,
                 **score_kwargs):
    """Bootstrap CIs and rank stability by resampling interactions.

    Returns:
        DataFrame indexed by receiver with columns score, ci_low, ci_high,
        rank_mean, rank_sd.
    """
    rng = np.random.default_rng(seed)
    obs = crosstalk_score(expr, interactions, availability, **score_kwargs)
    keys = list(obs["zscore"].index)
    scores, ranks = [], []
    kw = dict(score_kwargs); kw["enforce_floor"] = False
    for _ in range(int(n_boot)):
        pick = [keys[i] for i in rng.integers(0, len(keys), len(keys))]
        sub = interactions[interactions["interaction_name"].isin(set(pick))]
        try:
            s = crosstalk_score(expr, sub, availability, **kw)["score"]
        except (ValueError, AssertionError):
            continue
        scores.append(s)
        ranks.append(s.rank(ascending=False))
    S = pd.DataFrame(scores)
    R = pd.DataFrame(ranks)
    out = pd.DataFrame({
        "score": obs["score"],
        "ci_low": S.quantile(0.025),
        "ci_high": S.quantile(0.975),
        "rank_mean": R.mean(),
        "rank_sd": R.std(),
    }).sort_values("score", ascending=False)
    return out


def leave_one_out(expr, interactions, availability, by="ligand_subunit_gene",
                  **score_kwargs):
    """Stability control: drop each ligand/pathway in turn, correlate ranking.

    Returns:
        DataFrame with columns dropped, spearman_rho vs the full ranking.
    """
    kw = dict(score_kwargs); kw["enforce_floor"] = False
    full = crosstalk_score(expr, interactions, availability, **score_kwargs)["score"]
    rows = []
    for val in interactions[by].astype(str).unique():
        sub = interactions[interactions[by].astype(str) != val]
        if sub.empty:
            continue
        try:
            s = crosstalk_score(expr, sub, availability, **kw)["score"]
        except (ValueError, AssertionError):
            continue
        common = full.index.intersection(s.index)
        rows.append({"dropped": val,
                     "spearman_rho": float(stats.spearmanr(full[common],
                                                           s[common]).statistic)})
    return pd.DataFrame(rows).sort_values("spearman_rho")


def decompose(availability, capacity, interactions,
              ligand_col="ligand_subunit_gene", key_col="interaction_name"):
    """Split compatibility into availability vs capacity per receiver.

    Tests the two-factor hypothesis explicitly: a receiver can rank highly
    because the source releases strong ligands (availability, which is
    constant across receivers) or because it expresses the receptors
    (capacity, which varies). Only capacity can explain receiver differences,
    so this quantifies how much of the ranking is receiver-intrinsic.

    Returns:
        DataFrame indexed by receiver with mean_capacity_z,
        availability_weighted_capacity, and the availability total.
    """
    lig_of = (interactions.drop_duplicates(subset=[key_col])
              .set_index(key_col)[ligand_col].astype(str))
    sd = capacity.std(axis=1).replace(0.0, np.nan)
    z = capacity.sub(capacity.mean(axis=1), axis=0).div(sd, axis=0)
    w = pd.Series({k: float(availability.get(lig_of.get(k, ""), 0.0))
                   for k in capacity.index}).reindex(capacity.index).fillna(0.0)
    return pd.DataFrame({
        "mean_capacity_z": z.mean(axis=0),
        "availability_weighted_capacity": z.mul(w, axis=0).sum(axis=0) / max(w.sum(), 1e-9),
        "median_capacity_absolute": capacity.median(axis=0),
    }).sort_values("availability_weighted_capacity", ascending=False)


def write_manifest(path, inputs, params, outputs, notes=""):
    """Write a JSON manifest recording inputs, parameters and outputs."""
    man = {"inputs": inputs, "parameters": params, "outputs": outputs,
           "notes": notes}
    with open(path, "w") as fh:
        json.dump(man, fh, indent=1, default=str)
    return path


# ---------------------------------------------------------------------------
# Source-specificity by induction, and multi-part conjunction verdicts.
#
# Added after a screen returned zero surviving axes because source-tissue
# specificity was scored in HEALTHY tissue. A ligand that is unremarkable at
# baseline but strongly induced in disease is a plausible disease-source
# ligand, and a healthy-baseline filter discards it. Correcting that recovered
# axes that the original criterion had hidden.
#
# The correction also showed why a relaxation needs a guard: the recovered hits
# ran on receptors that 74% of all floor-passing receptors could have matched.
# induction_specificity() therefore returns the criteria SEPARATELY so the
# binding constraint is visible, and conjunction_verdict() refuses to call a
# partial pass a pass.
# ---------------------------------------------------------------------------


def induction_specificity(disease_de, baseline_rank, housekeeping_gene="GAPDH",
                          min_log2fc=1.0, max_fdr=0.05, min_housekeeping_ratio=0.005,
                          max_baseline_rank_frac=0.5, log_scale=True):
    """Score ligands on disease induction rather than healthy-tissue specificity.

    Three criteria, each computed WITHIN a single dataset so no cross-platform
    absolute comparison is implied (see SKILL.md on unit mismatches):

      A1  replicated induction: log2FC >= min_log2fc and FDR < max_fdr in EVERY
          cohort supplied. Requiring every cohort is the point; a single-cohort
          hit is not replication.
      A2  abundance guard: the ligand reaches min_housekeeping_ratio of the
          housekeeping gene INSIDE each disease cohort. This is the same
          housekeeping-relative device as housekeeping_floor(), applied to the
          source tissue instead of the receiver.
      A3  baseline non-exclusion: the source tissue ranks within the top
          max_baseline_rank_frac of tissues at baseline. This is deliberately
          weaker than "ranks first" -- it is the relaxation. Keep it explicit so
          a reader can see the hit depended on it.

    Args:
        disease_de: {cohort_name: DataFrame} with columns gene_symbol, log2FC,
            fdr, mean_expr. mean_expr is assumed log2 when log_scale is True.
        baseline_rank: DataFrame with columns gene, best_source_rank, n_tissues.
            A ligand absent here FAILS A3 -- it is never silently skipped.
        housekeeping_gene: reference for the A2 ratio.
        min_log2fc: A1 threshold, applied to the MINIMUM across cohorts.
        max_fdr: A1 significance threshold in every cohort.
        min_housekeeping_ratio: A2 threshold as a fraction of housekeeping.
        max_baseline_rank_frac: A3 threshold as a fraction of n_tissues.
        log_scale: whether mean_expr needs 2** before ratios are taken.

    Returns:
        DataFrame, one row per ligand, with per-criterion booleans, the value
        behind each, pass_all, and binding_constraint naming which criteria
        failed. Inspect binding_constraint before believing pass_all: if one
        criterion removes nearly everything, the screen is testing that
        criterion and not the biology.

    Raises:
        ValueError: if a cohort lacks the housekeeping gene or required columns.
    """
    if not disease_de:
        raise ValueError("disease_de is empty; supply at least one cohort")

    req = {"gene_symbol", "log2FC", "fdr", "mean_expr"}
    hk = {}
    for name, df in disease_de.items():
        missing = req - set(df.columns)
        if missing:
            raise ValueError(f"cohort {name!r} missing columns: {sorted(missing)}")
        row = df.loc[df["gene_symbol"] == housekeeping_gene, "mean_expr"]
        if row.empty:
            raise ValueError(
                f"housekeeping gene {housekeeping_gene!r} absent from cohort {name!r}; "
                "the A2 abundance guard cannot be computed without it")
        hk[name] = (2.0 ** float(row.mean())) if log_scale else float(row.mean())

    genes = sorted(set().union(*[set(df["gene_symbol"].dropna()) for df in disease_de.values()]))
    brank = baseline_rank.set_index("gene") if "gene" in baseline_rank.columns else baseline_rank

    rows = []
    for g in genes:
        fcs, fdrs, ratios, present = [], [], [], True
        for name, df in disease_de.items():
            sub = df.loc[df["gene_symbol"] == g]
            if sub.empty:
                present = False
                break
            fcs.append(float(sub["log2FC"].mean()))
            fdrs.append(float(sub["fdr"].min()))
            me = float(sub["mean_expr"].mean())
            ratios.append(((2.0 ** me) if log_scale else me) / hk[name])
        if not present:
            rows.append({"ligand": g, "A1_pass": False, "A2_pass": False, "A3_pass": False,
                         "min_log2FC": np.nan, "min_ratio_housekeeping": np.nan,
                         "baseline_rank": np.nan, "pass_all": False,
                         "binding_constraint": "absent_from_a_cohort"})
            continue

        a1 = bool(min(fcs) >= min_log2fc and max(fdrs) < max_fdr and min(fcs) > 0)
        a2 = bool(min(ratios) >= min_housekeeping_ratio)

        if g in brank.index:
            r = brank.loc[g]
            rank = float(r["best_source_rank"])
            ntis = float(r.get("n_tissues", np.nan))
            a3 = bool(np.isfinite(rank) and np.isfinite(ntis)
                      and rank <= max_baseline_rank_frac * ntis)
        else:
            rank, a3 = np.nan, False   # absent from baseline table => fails A3

        failed = [n for n, ok in (("A1", a1), ("A2", a2), ("A3", a3)) if not ok]
        rows.append({
            "ligand": g, "A1_pass": a1, "A2_pass": a2, "A3_pass": a3,
            "min_log2FC": min(fcs), "min_ratio_housekeeping": min(ratios),
            "baseline_rank": rank, "pass_all": not failed,
            "binding_constraint": "+".join(failed) if failed else "none",
        })
    return pd.DataFrame(rows).set_index("ligand").sort_values(
        ["pass_all", "min_log2FC"], ascending=[False, False])


def conjunction_verdict(results, require_all=True):
    """Score several sub-analyses as a conjunction, not a best-of.

    Where one hypothesis implies several independent predictions, testing each
    and reporting whichever clears significance inflates the false-positive
    rate by exactly the number of predictions. Requiring simultaneous agreement
    is the multiplicity control, and it has to be fixed BEFORE the tests run --
    deciding afterwards which conjuncts counted is the failure this guards.

    Args:
        results: {name: {"pass": bool, "direction_ok": bool, "detail": str}}.
            direction_ok defaults to the pass value when absent, but supply it
            explicitly: a conjunct can be significant in the WRONG direction,
            which is evidence against, not for.
        require_all: keep True unless a rule fixed in advance said otherwise.

    Returns:
        {"verdict": "PASS"|"FAIL", "n_passed", "n_total", "failing_conjuncts",
         "wrong_direction", "detail"}. A FAIL carries the failing names so the
        write-up cannot quietly become "n of m passed".

    Raises:
        ValueError: if results is empty or a conjunct lacks a "pass" key.
    """
    if not results:
        raise ValueError("results is empty; a conjunction needs at least one conjunct")

    passed, failing, wrong_dir = [], [], []
    for name, r in results.items():
        if "pass" not in r:
            raise ValueError(f"conjunct {name!r} has no 'pass' key")
        ok = bool(r["pass"])
        dir_ok = bool(r.get("direction_ok", ok))
        if not dir_ok:
            wrong_dir.append(name)
        if ok and dir_ok:
            passed.append(name)
        else:
            failing.append(name)

    verdict = "PASS" if (not failing if require_all else bool(passed)) else "FAIL"
    return {
        "verdict": verdict,
        "n_passed": len(passed),
        "n_total": len(results),
        "passing_conjuncts": passed,
        "failing_conjuncts": failing,
        "wrong_direction": wrong_dir,
        "detail": {k: v.get("detail", "") for k, v in results.items()},
    }


def receptor_pass_rate(receptor_table, receptor_col="receptor",
                       pass_col="recpass", min_types=1):
    """Base rate at which receptors satisfy the receiver criterion, ligand-free.

    Ask this BEFORE reading an axis count as evidence for a particular ligand.
    If most receptors in the universe pass with some tumour type regardless of
    ligand, an axis count is measuring the permissiveness of the receiver test.
    In the screen this was written for the rate was 74% -- which made a
    seven-axis result uninformative on its own.

    Returns:
        {"base_rate", "n_passing", "n_receptors", "distribution",
         "mean_types_per_receptor"} where base_rate is the fraction passing in
        at least min_types contexts.
    """
    if pass_col not in receptor_table.columns:
        raise ValueError(f"receptor_table has no {pass_col!r} column")
    per = receptor_table.groupby(receptor_col)[pass_col].sum()
    return {
        "base_rate": float((per >= min_types).mean()),
        "n_passing": int((per >= min_types).sum()),
        "n_receptors": int(per.size),
        "distribution": {int(k): int(v) for k, v in per.value_counts().sort_index().items()},
        "mean_types_per_receptor": float(per.mean()),
    }

