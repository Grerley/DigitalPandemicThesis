"""Numbers specific to Paper 2 (calibration under data scarcity).

The likelihood ratios that a published 2026 HBSC figure would give for the persistent (P) against the temporary (T)
explanation, computed as in ext_calibration.sec5 but under three priors on the solved quantities: the implicit prior
of the exact-solution calibration, a flat prior (with b1 <= 20 per year in Scenario P) and a log-uniform prior.
The draws and the re-weighting are those of ext_bayes.py; each scenario is re-weighted over every draw that solves it.

Writes outputs/results_p2.json (keys p2_...).
"""
from __future__ import annotations

import json, pathlib
import numpy as np
from scipy.stats import gaussian_kde

import run_all as ra
import ext_bayes as eb

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"
R: dict = {}
YEARS = [9, 10, 11, 12, 13, 14]


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["p2_" + k] = {"value": v, "note": note}


def scenario_rows(sc, allrows, cache):
    field = "b1" if sc == "P" else "d_alpha"
    rows = [r for r in allrows if r.get(field) is not None]
    I26, J, sol, b0v = [], [], [], []
    for r in rows:
        k = (r["eta"], r["sigma"], r["gamma1"], r["gamma2"], r["delta"], r["nu"], r["alpha"])
        if k in cache:
            ip, it, dd = cache[k]
            d1, d2 = dd[0], (dd[1] if sc == "P" else dd[2])
            path = ip if sc == "P" else it
        else:
            d1, d2 = eb.derivs_one(r, sc)
            path = ra.project(r, sc)[2]
        I26.append(100 * float(np.asarray(path)[2026 - 2018])); J.append(d1 * d2); sol.append(r[field]); b0v.append(r["b0"])
    return np.array(I26), np.array(J), np.array(sol), np.array(b0v)


def main():
    ok0 = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]]
    PP = np.loadtxt(OUT / "projection_P.csv", delimiter=","); PT = np.loadtxt(OUT / "projection_T.csv", delimiter=",")
    extra = json.load(open(OUT / "bayx_extra_draws.json"))
    okx = [r for r in extra if r["ok"]]
    ok = ok0 + okx
    PPk = np.vstack([PP, np.array([r["IP"] for r in okx])]); PTk = np.vstack([PT, np.array([r["IT"] for r in okx])])
    Dk = np.vstack([np.array([eb.derivs(r) for r in ok0]), np.array([r["D"] for r in okx])])
    key = lambda r: (r["eta"], r["sigma"], r["gamma1"], r["gamma2"], r["delta"], r["nu"], r["alpha"])
    cache = {key(r): (PPk[i], PTk[i], Dk[i]) for i, r in enumerate(ok)}
    allrows = json.load(open(OUT / "calibration_draws.json")) + extra

    samples = {}
    # implicit prior: draws kept in both scenarios, unweighted (the main calibration plus the extra draws)
    samples["implicit"] = {"P": (100 * PPk[:, 2026 - 2018], None), "T": (100 * PTk[:, 2026 - 2018], None)}
    xP, JP, solP, b0P = scenario_rows("P", allrows, cache)
    xT, JT, solT, b0T = scenario_rows("T", allrows, cache)
    samples["flat"] = {"P": (xP, np.where(solP <= 20, 1 / JP, 0.0)), "T": (xT, 1 / JT)}
    samples["log"] = {"P": (xP, 1 / (JP * b0P * solP)), "T": (xT, 1 / (JT * b0T * solT))}
    put("n_implicit", len(ok)); put("n_P_solvable", len(xP)); put("n_T_solvable", len(xT))

    for nm, d in samples.items():
        kd = {}
        for sc in "PT":
            x, w = d[sc]
            if w is not None:
                keep = w > 0; x, w = x[keep], w[keep] / w[keep].sum()
                put(f"{nm}_{sc}_ess", float(1 / np.sum(w ** 2)))
            kd[sc] = gaussian_kde(x, weights=w)
            put(f"{nm}_{sc}2026_median", float(eb.wq(x, w if w is not None else np.ones_like(x), 0.5)))
        for y in YEARS:
            pP = kd["P"].integrate_box_1d(y - 0.5, y + 0.5); pT = kd["T"].integrate_box_1d(y - 0.5, y + 0.5)
            put(f"{nm}_lr_{y}", float(pP / pT) if pT > 1e-12 else float("inf"))
            put(f"{nm}_post_P_{y}_pct", float(100 * pP / (pP + pT)))
            put(f"{nm}_pP_{y}", float(pP)); put(f"{nm}_pT_{y}", float(pT))
    json.dump(R, open(OUT / "results_p2.json", "w"), indent=1)
    for k, v in R.items():
        print(k, round(v["value"], 3) if isinstance(v["value"], float) else v["value"])


if __name__ == "__main__":
    main()
