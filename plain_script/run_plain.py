#!/usr/bin/env python3
"""Plain-script reproduction of the pan-cancer compatibility analysis.

No agent, no language-model call and no platform SDK is involved. The released
crosstalk_score module is imported as an ordinary Python file and run on the
released input tables. Each guardrail is then exercised on deliberately invalid
input to confirm that it halts the script.

Usage: python run_plain.py <inputs_dir> <module_path> <out_dir>
"""
import sys, os, json, time, platform, hashlib, importlib.util, traceback
import numpy as np, pandas as pd

inp, mod_path, out = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out, exist_ok=True)
spec = importlib.util.spec_from_file_location("crosstalk_score", mod_path)
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)

def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
expr = pd.read_csv(f"{inp}/tcga_expr_by_tumor_type.csv", index_col=0)
net = pd.read_csv(f"{inp}/crosstalk_network.csv")
pub = pd.read_csv(f"{inp}/pancancer_susceptibility_ranking.csv")
bg = json.load(open(f"{inp}/tcga_background_genes.json"))
avail = net.drop_duplicates("ligand_subunit_gene").set_index("ligand_subunit_gene")["ligand_availability"].astype(float)
inter = net[["interaction_name", "ligand_subunit_gene", "receptor_subunit_gene"]].copy()

log = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__,
       "module_sha256_16": sha(mod_path), "inputs_sha256_16": {f: sha(f"{inp}/{f}") for f in sorted(os.listdir(inp))}}
t0 = time.time()

# 1. Primary score: floor enforced, default settings.
res = cs.crosstalk_score(expr, inter, avail, rule="min", enforce_floor=False)
score = res["score"]
floor_audit = res["floor_audit"]
scored = floor_audit[floor_audit.passes].interaction.tolist()
res = cs.crosstalk_score(expr, inter[inter.interaction_name.isin(scored)], avail, rule="min", enforce_floor=True)
score = res["score"]
log["floor"] = res["floor"]; log["n_interactions_scored"] = res["n_interactions"]
log["n_interactions_below_floor"] = int((~floor_audit.passes).sum())

# 2. Matched-random null on the floor-passing network.
B = int(os.environ.get("N_NULL", "300"))
# Configuration A, as in the published floor-sensitivity analysis: full network passed,
# floor applied inside scoring without raising. Configuration B: only floor-passing
# interactions passed, so null draws replace only receptors that are actually scored.
nul = cs.matched_random_null(expr, inter, avail, bg["genes"], n_iter=B, seed=20260727, enforce_floor=False)
nulB = cs.matched_random_null(expr, inter[inter.interaction_name.isin(scored)], avail,
                              bg["genes"], n_iter=B, seed=20260727)
log["null_n_iter"] = B
for k, v in (("A", nul), ("B", nulB)):
    log[f"null_{k}_rho_mean"] = v["rho_mean"]; log[f"null_{k}_frac_rho_gt_0_5"] = v["frac_rho_above_0_5"]

# 3. Same analysis with the floor removed (the failure mode the floor prevents).
nofloor = cs.crosstalk_score(expr, inter, avail, rule="min", enforce_floor=False, floor=0.0)["score"]

tab = pd.DataFrame({"tumor_type": score.index, "score_plain": score.values})
tab["rank_plain"] = tab.score_plain.rank(ascending=False, method="first").astype(int)
tab["null_p_plain_A"] = tab.tumor_type.map(nul["per_receiver_p"])
tab["null_p_plain_B"] = tab.tumor_type.map(nulB["per_receiver_p"])
tab["score_no_floor_plain"] = tab.tumor_type.map(nofloor)
tab = tab.merge(pub[["tumor_type", "rank", "score", "null_p_primary", "score_no_floor"]], on="tumor_type", how="outer")
tab.to_csv(f"{out}/plain_script_ranking.csv", index=False)
log["runtime_s_primary"] = round(time.time() - t0, 1)

# 4. Guardrail challenge tests: each must raise on invalid input.
def expect_raise(name, fn, exc):
    try:
        fn(); return {"check": name, "halted": False, "exception": None, "message": None}
    except exc as e:
        return {"check": name, "halted": True, "exception": type(e).__name__, "message": str(e).splitlines()[0][:160]}
    except Exception as e:
        return {"check": name, "halted": True, "exception": type(e).__name__ + " (unexpected type)", "message": str(e).splitlines()[0][:160]}

below = floor_audit[~floor_audit.passes].interaction.tolist()
tests = [
  expect_raise("Expression floor: scoring includes below-floor interactions",
               lambda: cs.crosstalk_score(expr, inter, avail, enforce_floor=True), AssertionError),
  expect_raise("Floor reference gene missing from matrix",
               lambda: cs.crosstalk_score(expr.drop(index="GAPDH"), inter[inter.interaction_name.isin(scored)], avail), KeyError),
  expect_raise("Ligand availability does not match network ligands",
               lambda: cs.crosstalk_score(expr, inter[inter.interaction_name.isin(scored)], pd.Series({"NOT_A_GENE": 1.0})), ValueError),
  expect_raise("Null background pool too small",
               lambda: cs.matched_random_null(expr, inter[inter.interaction_name.isin(scored)], avail, bg["genes"][:5], n_iter=2), ValueError),
  expect_raise("Unknown receptor-complex rule",
               lambda: cs.crosstalk_score(expr, inter[inter.interaction_name.isin(scored)], avail, rule="max"), ValueError),
  expect_raise("Empty conjunction",
               lambda: cs.conjunction_verdict([]), ValueError),
]
tests.append({"check": "Control: valid input runs without halting", "halted": False, "exception": None,
              "message": f"{log['n_interactions_scored']} interactions scored"})
pd.DataFrame(tests).to_csv(f"{out}/plain_script_guardrail_tests.csv", index=False)
log["n_below_floor_examples"] = below[:5]
json.dump(log, open(f"{out}/plain_script_run_log.json", "w"), indent=1, default=str)
print("done", round(time.time() - t0, 1), "s")
