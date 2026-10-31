"""Derived results added to outputs/results.json after run_all.py (thesis).
Recovery-study relative bias and error from outputs/recovery.csv, and the overshoot share in per cent."""
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

# do the calibration failures depend on the solver bounds? re-solve rejected draws with bounds ten times wider
from dpm import calibrate as cal
from dpm.model import equilibrium
rows = json.load(open("outputs/calibration_draws.json"))
bad = [r for r in rows if r.get("b0") is not None and not r["ok"]]
resP = resT = 0
for r in bad:
    d = {k: r[k] for k in cal.PARAM_ORDER}
    p0 = Params(b=r["b0"], **d); x18 = equilibrium(p0)
    if r["b1"] is None and cal.solve_b1(p0, x18, r["t22"], bmax=2000.0) is not None:
        resP += 1
    if r["d_alpha"] is None and cal.solve_dalpha(p0, x18, r["t22"], dmax=200.0) is not None:
        resT += 1
put("cal_n_failP", sum(1 for r in bad if r["b1"] is None)); put("cal_n_failT", sum(1 for r in bad if r["d_alpha"] is None))
put("cal_n_rescued_P_x10", resP); put("cal_n_rescued_T_x10", resT)

json.dump(R, open("outputs/results.json", "w"), indent=1, ensure_ascii=False)
print("derived results added")
