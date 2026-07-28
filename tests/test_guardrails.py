"""Regression tests for the guardrails that make this analysis trustworthy.

Run:  python tests/test_guardrails.py
Exits non-zero on any failure. No network access required.
"""

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "skills", "crosstalk-score"))
import kernel as ck  # noqa: E402

FAILURES = []


def check(name, condition, detail=""):
    status = "PASS" if condition else "FAIL"
    print("[%s] %s%s" % (status, name, ("  -- " + detail) if detail else ""))
    if not condition:
        FAILURES.append(name)


def toy():
    """A 4-gene toy cohort: NPR1 and FZD7 expressed, OPRK1 at noise level."""
    expr = pd.DataFrame(
        {"A": [1000.0, 50000.0, 2.0, 900.0],
         "B": [1200.0, 52000.0, 1.5, 1100.0],
         "C": [800.0, 48000.0, 3.0, 700.0]},
        index=["NPR1", "GAPDH", "OPRK1", "FZD7"])
    inter = pd.DataFrame({
        "interaction_name": ["NPPA_NPR1", "PENK_OPRK1", "WNT9A_FZD7"],
        "ligand_subunit_gene": ["NPPA", "PENK", "WNT9A"],
        "receptor_subunit_gene": ["NPR1", "OPRK1", "FZD7"]})
    avail = pd.Series({"NPPA": 3.0, "PENK": 2.8, "WNT9A": 2.3})
    return expr, inter, avail


def test_floor_is_unit_agnostic():
    expr, _, _ = toy()
    floor = ck.housekeeping_floor(expr, "GAPDH", 0.005)
    check("floor derives from the housekeeping reference, not a hardcoded unit",
          abs(floor - 250.0) < 1e-9, "floor=%.2f" % floor)
    try:
        ck.housekeeping_floor(expr, "NOT_A_GENE", 0.005)
        check("missing floor reference raises", False)
    except KeyError:
        check("missing floor reference raises KeyError", True)


def test_floor_rejects_untranscribed_receptor():
    """The central guardrail: a receptor at noise level must not be scored."""
    expr, inter, avail = toy()
    floor = ck.housekeeping_floor(expr, "GAPDH", 0.005)
    try:
        ck.crosstalk_score(expr, inter, avail, floor=floor, enforce_floor=True)
        check("below-floor receptor raises by default", False,
              "OPRK1 was scored despite being at noise level")
    except AssertionError as exc:
        check("below-floor receptor raises by default", "OPRK1" in str(exc),
              "raised naming the offender")

    res = ck.crosstalk_score(expr, inter, avail, floor=floor,
                             enforce_floor=False)
    check("enforce_floor=False drops rather than scores it",
          res["dropped"] == ["PENK_OPRK1"] and res["n_interactions"] == 2,
          "dropped=%s scored=%d" % (res["dropped"], res["n_interactions"]))
    check("floor audit accounts for every candidate interaction",
          len(res["floor_audit"]) == 3)


def test_complex_rule_is_obligate_subunit():
    """min-of-subunits: a complex is limited by its least-expressed subunit."""
    expr = pd.DataFrame({"A": [1000.0, 10.0], "B": [1000.0, 10.0]},
                        index=["FZD1", "LRP5"])
    inter = pd.DataFrame({
        "interaction_name": ["WNT9A_FZD1_LRP5", "WNT9A_FZD1_LRP5"],
        "ligand_subunit_gene": ["WNT9A", "WNT9A"],
        "receptor_subunit_gene": ["FZD1", "LRP5"]})
    cap = ck.receptor_capacity(expr, inter, rule="min")
    check("min rule takes the limiting subunit",
          abs(float(cap.loc["WNT9A_FZD1_LRP5", "A"]) - 10.0) < 1e-9)
    cap_gm = ck.receptor_capacity(expr, inter, rule="geomean")
    check("geomean rule differs from min (sensitivity alternative)",
          float(cap_gm.loc["WNT9A_FZD1_LRP5", "A"]) > 10.0)
    try:
        ck.receptor_capacity(expr, inter, rule="not_a_rule")
        check("unknown complex rule raises", False)
    except ValueError:
        check("unknown complex rule raises ValueError", True)


def test_mismatched_availability_raises():
    """A silent all-zero weighting would produce a meaningless score."""
    expr, inter, _ = toy()
    floor = ck.housekeeping_floor(expr, "GAPDH", 0.005)
    bad = pd.Series({"NOT_A_LIGAND": 1.0})
    try:
        ck.crosstalk_score(expr, inter, bad, floor=floor, enforce_floor=False)
        check("all-zero availability weights raise", False)
    except ValueError:
        check("all-zero availability weights raise ValueError", True)


def test_double_centering_removes_additive_effects():
    m = pd.DataFrame(np.arange(12, dtype=float).reshape(3, 4))
    dc = ck.double_center(m)
    check("double-centering zeroes row and column deviations",
          np.allclose(dc.sub(dc.mean(axis=1), axis=0)
                        .sub(dc.mean(axis=0), axis=1).values,
                      dc.values - dc.values.mean(), atol=1e-9))


def test_published_results_are_self_consistent():
    """Guard the two documented traps against the shipped tables."""
    tdir = os.path.join(ROOT, "results", "tables")
    sig = os.path.join(tdir, "hf_signature_replicated.csv")
    if os.path.exists(sig):
        cols = pd.read_csv(sig, nrows=1).columns
        check("signature table still carries both id columns (join trap)",
              "gene" in cols and "gene_symbol" in cols)
    rank = os.path.join(tdir, "pancancer_susceptibility_ranking.csv")
    if os.path.exists(rank):
        rk = pd.read_csv(rank)
        n_sig = int((rk["null_p_primary"] < 0.05).sum())
        check("4 of 30 tumour types beat the specificity null",
              n_sig == 4, "observed %d" % n_sig)
        check("ranking is sorted by score",
              rk["score"].is_monotonic_decreasing)
    net = os.path.join(tdir, "crosstalk_network.csv")
    if os.path.exists(net):
        nt = pd.read_csv(net)
        op = nt[nt["pathway_name"].astype(str).str.upper() == "OPIOID"]
        check("all four PENK-opioid axes are below the floor",
              len(set(op["interaction_name"])) == 4
              and not op["passed_floor"].any(),
              "%d axes, any passing: %s"
              % (len(set(op["interaction_name"])), bool(op["passed_floor"].any())))


if __name__ == "__main__":
    for fn in [test_floor_is_unit_agnostic,
               test_floor_rejects_untranscribed_receptor,
               test_complex_rule_is_obligate_subunit,
               test_mismatched_availability_raises,
               test_double_centering_removes_additive_effects,
               test_published_results_are_self_consistent]:
        print("\n== %s ==" % fn.__name__)
        fn()
    print("\n%s" % ("-" * 56))
    if FAILURES:
        print("FAILED: %s" % ", ".join(FAILURES))
        sys.exit(1)
    print("All guardrail tests passed.")
