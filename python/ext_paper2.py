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
    # weighted medians of the solved quantities under each prior
    for nm, (wP_, wT_) in {"flat": (samples["flat"]["P"][1], samples["flat"]["T"][1]),
                           "log": (samples["log"]["P"][1], samples["log"]["T"][1])}.items():
        put(f"{nm}_b1_median", eb.wq(solP, wP_, 0.5)); put(f"{nm}_dalpha_median", eb.wq(solT, wT_, 0.5))
    put("implicit_b1_median", float(np.median([r["b1"] for r in ok])))
    put("implicit_dalpha_median", float(np.median([r["d_alpha"] for r in ok])))
    # flat priors with bounds on the solved quantity of each scenario: sensitivity of the next-wave reading
    def lr_pair(wP_, wT_):
        kP = gaussian_kde(xP[wP_ > 0], weights=wP_[wP_ > 0] / wP_.sum())
        kT = gaussian_kde(xT[wT_ > 0], weights=wT_[wT_ > 0] / wT_.sum())
        out = {}
        for y in YEARS:
            pP = kP.integrate_box_1d(y - 0.5, y + 0.5); pT = kT.integrate_box_1d(y - 0.5, y + 0.5)
            out[y] = (pP, pT, pP / pT if pT > 1e-12 else float("inf"))
        return out
    for cb in (20, 10, 5):
        wP_ = np.where(solP <= cb, 1 / JP, 0.0)
        put(f"cap_b{cb}_P2026_median", eb.wq(xP, wP_, 0.5)); put(f"cap_b{cb}_b1_median", eb.wq(solP, wP_, 0.5))
        put(f"cap_b{cb}_P_ess", float(wP_.sum() ** 2 / np.sum(wP_ ** 2)))
    for cd in (20, 10, 5, 2):
        wT_ = np.where(solT <= cd, 1 / JT, 0.0)
        put(f"cap_d{cd}_T2026_median", eb.wq(xT, wT_, 0.5)); put(f"cap_d{cd}_dalpha_median", eb.wq(solT, wT_, 0.5))
        put(f"cap_d{cd}_T_ess", float(wT_.sum() ** 2 / np.sum(wT_ ** 2)))
    for cb, cd in [(20, 20), (20, 5), (10, 10), (5, 5), (5, 2)]:
        o = lr_pair(np.where(solP <= cb, 1 / JP, 0.0), np.where(solT <= cd, 1 / JT, 0.0))
        for y in YEARS:
            put(f"cap_b{cb}_d{cd}_lr_{y}", float(o[y][2])); put(f"cap_b{cb}_d{cd}_pT_{y}", float(o[y][1]))
            put(f"cap_b{cb}_d{cd}_pP_{y}", float(o[y][0]))
    # bootstrap of the log-uniform ratios at 10% and 11%
    g = np.random.default_rng(2026)
    wPl, wTl = samples["log"]["P"][1], samples["log"]["T"][1]
    bs = {10: [], 11: []}
    for _ in range(200):
        iP = g.integers(0, len(xP), len(xP)); iT = g.integers(0, len(xT), len(xT))
        kP = gaussian_kde(xP[iP], weights=wPl[iP] / wPl[iP].sum()); kT = gaussian_kde(xT[iT], weights=wTl[iT] / wTl[iT].sum())
        for y in bs:
            bs[y].append(kP.integrate_box_1d(y - 0.5, y + 0.5) / kT.integrate_box_1d(y - 0.5, y + 0.5))
    for y in bs:
        put(f"log_lr_{y}_boot_p05", float(np.percentile(bs[y], 5))); put(f"log_lr_{y}_boot_p95", float(np.percentile(bs[y], 95)))
    json.dump(R, open(OUT / "results_p2.json", "w"), indent=1)
    monotonicity(allrows)
    json.dump(R, open(OUT / "results_p2.json", "w"), indent=1)
    for k, v in R.items():
        print(k, round(v["value"], 3) if isinstance(v["value"], float) else v["value"])


def monotonicity(allrows, n=100, m=60):
    """Is the 2022 prevalence increasing in the solved quantity? Log-spaced grids over the whole search range,
    for draws that solve each scenario (including draws that solve only one)."""
    import dpm.calibrate as cal
    from dpm.model import Params
    g = np.random.default_rng(7)
    for sc, field in (("P", "b1"), ("T", "d_alpha")):
        rows = [r for r in allrows if r.get(field) is not None]
        pick = [rows[i] for i in g.choice(len(rows), n, replace=False)]
        good = 0; worst = np.inf
        for r in pick:
            d = {k: r[k] for k in cal.PARAM_ORDER}; p0 = Params(b=r["b0"], **d); x18 = np.array(r["x18"])
            if sc == "P":
                grid = np.geomspace(r["b0"], 200.0, m)
                v = [cal.I_after(p0.with_(b=b), x18, 4.0) for b in grid]
            else:
                grid = np.concatenate([[0.0], np.geomspace(1e-4, 20.0, m - 1)])
                v = [cal.I_after(p0, x18, 4.0, alpha_t=cal.shock(p0.alpha, z)) for z in grid]
            dv = np.diff(v); good += bool(np.all(dv > 0)); worst = min(worst, float(dv.min()))
        put(f"mono_{sc}_n", n); put(f"mono_{sc}_ok", good); put(f"mono_{sc}_grid", m); put(f"mono_{sc}_min_step", worst)


if __name__ == "__main__":
    main()
