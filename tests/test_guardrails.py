"""Regression tests for the guardrails that make this analysis trustworthy.

Run:  python tests/test_guardrails.py
Exits non-zero on any failure. No network access required.
"""

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from hf_crosstalk import crosstalk_score as ck  # noqa: E402

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


def test_conjunction_refuses_partial_pass():
    """A hypothesis with three predictions needs all three, in direction."""
    one_of_three = {
        "survival": {"pass": False, "direction_ok": False},
        "proliferation": {"pass": True, "direction_ok": True},
        "copy_number": {"pass": False, "direction_ok": False},
    }
    v = ck.conjunction_verdict(one_of_three)
    check("1 of 3 conjuncts is a FAIL, not a discovery",
          v["verdict"] == "FAIL" and v["n_passed"] == 1,
          "verdict=%s n_passed=%s" % (v["verdict"], v["n_passed"]))
    check("failing conjuncts are named",
          set(v["failing_conjuncts"]) == {"survival", "copy_number"},
          "failing=%s" % v["failing_conjuncts"])
    all_three = {k: {"pass": True, "direction_ok": True}
                 for k in ("a", "b", "c")}
    check("all three passing is a PASS",
          ck.conjunction_verdict(all_three)["verdict"] == "PASS")
    wrong_dir = {"a": {"pass": True, "direction_ok": False},
                 "b": {"pass": True, "direction_ok": True}}
    check("significant in the wrong direction does not count as passed",
          ck.conjunction_verdict(wrong_dir)["verdict"] == "FAIL")
    try:
        ck.conjunction_verdict({})
        check("empty conjunction raises", False)
    except ValueError:
        check("empty conjunction raises ValueError", True)


def test_receptor_base_rate_guard():
    """Before reading an axis count as evidence, ask the ligand-free rate."""
    tbl = pd.DataFrame({
        "receptor": ["R1", "R1", "R2", "R2", "R3", "R3"],
        "recpass": [True, True, True, False, False, False]})
    out = ck.receptor_pass_rate(tbl, min_types=1)
    check("base rate counts receptors passing with no ligand involved",
          abs(out["base_rate"] - 2 / 3) < 1e-9 and out["n_receptors"] == 3,
          "base_rate=%.3f n=%d" % (out["base_rate"], out["n_receptors"]))
    check("mean types per receptor is reported",
          abs(out["mean_types_per_receptor"] - 1.0) < 1e-9)
    try:
        ck.receptor_pass_rate(tbl.drop(columns=["recpass"]))
        check("missing pass column raises", False)
    except ValueError:
        check("missing pass column raises ValueError", True)


def test_induction_criteria_are_within_dataset():
    """Source specificity by induction: three criteria, reported separately."""
    de = pd.DataFrame({
        "gene_symbol": ["PTN", "QUIET", "GAPDH"],
        "log2FC": [1.6, 0.2, 0.0],
        "fdr": [1e-4, 0.5, 0.9],
        "mean_expr": [np.log2(1600.0), np.log2(3.0), np.log2(69000.0)]})
    base = pd.DataFrame({"gene": ["PTN", "QUIET"],
                         "best_source_rank": [23, 2],
                         "n_tissues": [52, 52]})
    out = ck.induction_specificity({"c1": de, "c2": de}, base)
    check("criteria are returned separately so the binding one is visible",
          all(c in out.columns for c in ("A1_pass", "A2_pass", "A3_pass",
                                        "pass_all")))
    check("an induced, abundant, upper-half ligand passes all three",
          bool(out.loc["PTN", "pass_all"]))
    check("an unchanged low-abundance ligand fails induction and abundance",
          not bool(out.loc["QUIET", "A1_pass"])
          and not bool(out.loc["QUIET", "A2_pass"]))
    try:
        ck.induction_specificity({"c1": de.drop(columns=["fdr"])}, base)
        check("missing required column raises", False)
    except ValueError:
        check("missing required column raises ValueError", True)
    try:
        ck.induction_specificity({"c1": de[de.gene_symbol != "GAPDH"]}, base)
        check("missing housekeeping reference raises", False)
    except ValueError:
        check("missing housekeeping reference raises ValueError", True)


def test_ptn_withdrawal_is_recorded():
    """The adversarial phase's verdicts must survive in the shipped tables."""
    tdir = os.path.join(ROOT, "results", "tables")
    bra = os.path.join(tdir, "receptor_base_rate_audit.csv")
    if os.path.exists(bra):
        b = pd.read_csv(bra)
        rate = float(b["passes_ge1_type"].mean())
        check("receptor base rate is 74.1%, above the 60% uninformative bar",
              abs(rate - 0.7407) < 0.001, "rate=%.4f" % rate)
    nl = os.path.join(tdir, "ptn_receptor_null.csv")
    if os.path.exists(nl):
        draws = pd.read_csv(nl)["total_pairs_primary"].values
        check("PTN's 27 observed pairs do NOT exceed the null 95th percentile",
              27 <= np.percentile(draws, 95),
              "obs=27 p95=%.1f" % np.percentile(draws, 95))
    vs = os.path.join(tdir, "validation_summary.csv")
    if os.path.exists(vs):
        v = pd.read_csv(vs)
        row = v[v["n"] == 50]
        check("the PTN retention rule is recorded as FAIL",
              len(row) == 1 and row.iloc[0]["verdict"] == "FAIL",
              "verdict=%s" % (row.iloc[0]["verdict"] if len(row) else "missing"))


if __name__ == "__main__":
    for fn in [test_floor_is_unit_agnostic,
               test_floor_rejects_untranscribed_receptor,
               test_complex_rule_is_obligate_subunit,
               test_mismatched_availability_raises,
               test_double_centering_removes_additive_effects,
               test_published_results_are_self_consistent,
               test_conjunction_refuses_partial_pass,
               test_receptor_base_rate_guard,
               test_induction_criteria_are_within_dataset,
               test_ptn_withdrawal_is_recorded]:
        print("\n== %s ==" % fn.__name__)
        fn()
    print("\n%s" % ("-" * 56))
    if FAILURES:
        print("FAILED: %s" % ", ".join(FAILURES))
        sys.exit(1)
    print("All guardrail tests passed.")
