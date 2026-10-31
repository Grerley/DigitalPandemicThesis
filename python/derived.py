"""Derived results added to outputs/results.json after run_all.py (thesis v2).
Recovery-study relative bias from outputs/recovery.csv, overshoot share in per cent, and a targeting
experiment with run-level spread (paired differences)."""
import csv, json, numpy as np, math
import run_all as ra
from dpm.model import Params
from dpm import network as nw

R = json.load(open("outputs/results.json"))
def put(k, v, note=""): R[k] = {"value": float(v) if isinstance(v, (np.floating, float)) else v, "note": note}

for row in csv.DictReader(open("outputs/recovery.csv")):
    tr = row["transition"].replace("→", "_")
    put(f"rec_{row['noise']}_{row['n']}_{tr}_relbias", 100 * float(row["bias"]) / float(row["true"]), "per cent")
    put(f"rec_{row['noise']}_{row['n']}_{tr}_relrmse_pct", 100 * float(row["rmse"]) / float(row["true"]), "per cent")
put("P_share_overshoot_pct", 100 * R["P_share_overshoot"]["value"])
put("int_n_draws_eff", R["int_n_draws"]["value"])

# targeting with paired runs: same graph and seed for each strategy
ranges = json.load(open("data/parameter_ranges.json"))
pc0, pc1, da, x18 = ra.central(ranges)
from dpm.model import simulate
_, y = simulate(pc1, x18, 4, t_eval=[4]); x22 = y[:, -1]
beta = nw.beta_from_b(pc1, 10)
reps, n = 60, 2000
diffs_r, diffs_h, base = [], [], []
for s in range(reps):
    G = nw.make_graph("BA", n, 10, seed=ra.SEED + 500 + s)
    deg = np.array([d for _, d in G.degree()]); rng = np.random.default_rng(ra.SEED + 600 + s)
    m = int(0.10 * n)
    w0 = np.ones(n); wr = np.ones(n); wr[rng.choice(n, m, replace=False)] = 0.2
    wh = np.ones(n); wh[np.argsort(-deg)[:m]] = 0.2
    seed = ra.SEED + 700 + s
    a = ra.simulate_weighted(G, pc1, beta, 240, x22, w0, seed)
    b = ra.simulate_weighted(G, pc1, beta, 240, x22, wr, seed)
    c = ra.simulate_weighted(G, pc1, beta, 240, x22, wh, seed)
    base.append(a); diffs_r.append(a - b); diffs_h.append(a - c)
dr, dh = np.array(diffs_r) * 100, np.array(diffs_h) * 100
put("tgt_reps", reps); put("tgt_base_mean", 100 * np.mean(base))
put("tgt_random_reduction_mean", dr.mean()); put("tgt_random_reduction_se", dr.std(ddof=1) / math.sqrt(reps))
put("tgt_hub_reduction_mean", dh.mean()); put("tgt_hub_reduction_se", dh.std(ddof=1) / math.sqrt(reps))
dd = dh - dr
put("tgt_hub_minus_random_mean", dd.mean()); put("tgt_hub_minus_random_se", dd.std(ddof=1) / math.sqrt(reps))
json.dump(R, open("outputs/results.json", "w"), indent=1, ensure_ascii=False)
print({k: round(R[k]["value"], 3) for k in R if k.startswith("tgt_")})
