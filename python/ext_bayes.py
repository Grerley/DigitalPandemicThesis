"""A Bayesian reading of the calibration (Chapter 7, section "The calibration as a Bayesian procedure").

Solving (b0, b1) exactly for targets drawn uniformly within the rounding intervals samples the solved
quantities with density proportional to the Jacobian determinant of the map from (b0, b1) to the two
targets. The procedure is therefore a Bayesian calibration with a rounding likelihood and an implicit
prior on (b0, b1) proportional to that determinant. Re-weighting each kept draw by the inverse
determinant gives the posterior under a flat prior on (b0, b1); a further factor 1/(b0 b1) gives a
log-uniform prior. The same holds for (b0, d_alpha) in Scenario T.

Writes outputs/results_bayx.json (keys bayx_...) and figures/fig_bayx_weights.png.
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt
from dpm.model import Params, equilibrium, R0
from dpm import calibrate as cal

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R: dict = {}
N_EXTRA = 8000
SEED_EXTRA = 20270101


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["bayx_" + k] = {"value": v, "note": note}


def wq(x, w, q):
    """Weighted quantile (q in [0, 1]) by linear interpolation of the weighted empirical CDF."""
    o = np.argsort(x); x, w = np.asarray(x)[o], np.asarray(w)[o]
    c = (np.cumsum(w) - 0.5 * w) / w.sum()
    return float(np.interp(q, c, x))


def ess(w):
    return float(w.sum() ** 2 / (w ** 2).sum())


def derivs(r, h=1e-5):
    d = {k: r[k] for k in cal.PARAM_ORDER}
    p0 = Params(b=r["b0"], **d)
    # d I_eq / d b0 (central difference on the log scale, converted)
    f = lambda b: equilibrium(p0.with_(b=b))[2]
    dI18_db0 = (f(r["b0"] * (1 + h)) - f(r["b0"] * (1 - h))) / (2 * h * r["b0"])
    x18 = np.array(r["x18"])
    g = lambda b: cal.I_after(p0.with_(b=b), x18, 4.0)
    dI22_db1 = (g(r["b1"] * (1 + h)) - g(r["b1"] * (1 - h))) / (2 * h * r["b1"])
    s = lambda da: cal.I_after(p0, x18, 4.0, alpha_t=cal.shock(p0.alpha, da))
    hd = h * max(r["d_alpha"], 1e-3)
    dI22_dda = (s(r["d_alpha"] + hd) - s(r["d_alpha"] - hd)) / (2 * hd)
    return dI18_db0, dI22_db1, dI22_dda


def derivs_one(r, sc, h=1e-5):
    """Derivatives for one scenario only (for draws that solve one scenario but not the other)."""
    d = {k: r[k] for k in cal.PARAM_ORDER}
    p0 = Params(b=r["b0"], **d)
    f = lambda b: equilibrium(p0.with_(b=b))[2]
    d1 = (f(r["b0"] * (1 + h)) - f(r["b0"] * (1 - h))) / (2 * h * r["b0"])
    x18 = np.array(r["x18"])
    if sc == "P":
        g = lambda b: cal.I_after(p0.with_(b=b), x18, 4.0)
        d2 = (g(r["b1"] * (1 + h)) - g(r["b1"] * (1 - h))) / (2 * h * r["b1"])
    else:
        s_ = lambda z: cal.I_after(p0, x18, 4.0, alpha_t=cal.shock(p0.alpha, z))
        hd = h * max(r["d_alpha"], 1e-3)
        d2 = (s_(r["d_alpha"] + hd) - s_(r["d_alpha"] - hd)) / (2 * hd)
    return d1, d2


def main():
    t0 = time.time()
    ok0 = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]]
    PP = np.loadtxt(OUT / "projection_P.csv", delimiter=",")
    PT = np.loadtxt(OUT / "projection_T.csv", delimiter=",")
    assert PP.shape[0] == len(ok0)
    D0 = np.array([derivs(r) for r in ok0])
    # unweighted summaries of the main calibration, for reference
    put("main_n", len(ok0))
    # additional draws from the same ranges (different seed) to give the re-weighting an adequate effective size
    extra_path = OUT / "bayx_extra_draws.json"
    if extra_path.exists():
        extra = json.load(open(extra_path))
    else:
        ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
        extra = []
        for blk in range(N_EXTRA // 500):
            rows = cal.calibrate(ranges, n=500, seed=SEED_EXTRA + blk)
            for r in rows:
                if r["ok"]:
                    yP, yT = ra.project(r, "P"), ra.project(r, "T")
                    r["IP"] = [float(v) for v in yP[2]]; r["IT"] = [float(v) for v in yT[2]]
                    r["D"] = list(derivs(r))
                extra.append(r)
            print("block", blk, time.time() - t0, flush=True)
        json.dump(extra, open(extra_path, "w"))
    okx = [r for r in extra if r["ok"]]
    put("extra_n_draws", len(extra)); put("extra_n_ok", len(okx))
    ok = ok0 + okx
    put("n", len(ok))
    # ---- as calibrated: draws kept in both scenarios, unweighted
    PPk = np.vstack([PP, np.array([r["IP"] for r in okx])]); PTk = np.vstack([PT, np.array([r["IT"] for r in okx])])
    Dk = np.vstack([D0, np.array([r["D"] for r in okx])])
    JPk = Dk[:, 0] * Dk[:, 1]
    put("JP_ratio_p95_p05", float(np.percentile(JPk, 95) / np.percentile(JPk, 5)))
    put("JT_ratio_p95_p05", float(np.percentile(Dk[:, 0] * Dk[:, 2], 95) / np.percentile(Dk[:, 0] * Dk[:, 2], 5)))
    put("dI18_db0_median", float(np.median(Dk[:, 0]))); put("dI22_db1_median", float(np.median(Dk[:, 1])))
    for sc, I in [("P", PPk), ("T", PTk)]:
        for yr in [2026, 2042]:
            jj = yr - 2018
            put(f"{sc}{yr}_unw_median", 100 * float(np.median(I[:, jj])))
            put(f"{sc}{yr}_unw_p05", 100 * float(np.percentile(I[:, jj], 5)))
            put(f"{sc}{yr}_unw_p95", 100 * float(np.percentile(I[:, jj], 95)))
        put(f"{sc}_gamma2_unw_median", float(np.median([r["gamma2"] for r in ok])))
    # ---- re-weighting: each scenario over every draw that solves that scenario (not only draws kept in both)
    allrows = json.load(open(OUT / "calibration_draws.json")) + extra
    key = lambda r: (r["eta"], r["sigma"], r["gamma1"], r["gamma2"], r["delta"], r["nu"], r["alpha"])
    cache = {key(r): (r, PPk[i], PTk[i], Dk[i]) for i, r in enumerate(ok)}
    for sc in ["P", "T"]:
        field = "b1" if sc == "P" else "d_alpha"
        rows = [r for r in allrows if r.get(field) is not None]
        I, J, sol, b0v = [], [], [], []
        n_new = 0
        for r in rows:
            if key(r) in cache:
                _, ip, it, dd = cache[key(r)]
                d1, d2 = dd[0], (dd[1] if sc == "P" else dd[2])
                path = ip if sc == "P" else it
            else:
                n_new += 1
                dd = derivs_one(r, sc)
                d1, d2 = dd
                path = ra.project(r, sc)[2]
            I.append(np.asarray(path)); J.append(d1 * d2); sol.append(r[field]); b0v.append(r["b0"])
        I = np.array(I); J = np.array(J); sol = np.array(sol); b0v = np.array(b0v)
        assert np.all(J > 0)
        gam = np.array([r["gamma2"] for r in rows])
        put(f"{sc}_n_solvable", len(rows)); put(f"{sc}_n_solvable_only", n_new)
        W = {"flat": 1 / J, "log": 1 / (J * b0v * sol)}
        if sc == "P":
            for cap in [20, 50]:
                W[f"flat{cap}"] = np.where(sol <= cap, 1 / J, 0.0)
            put("P_b1_flat_median", wq(sol, W["flat"], 0.5)); put("P_b1_flat_p95", wq(sol, W["flat"], 0.95))
            put("P_b1_unw_median", float(np.median([r["b1"] for r in ok])))
            put("corr_logwPflat_loggamma2", float(np.corrcoef(np.log(W["flat"]), np.log(gam))[0, 1]))
            put("P_flat_weight_share_Ponly_pct", float(100 * W["flat"][[key(r) not in cache for r in rows]].sum() / W["flat"].sum()))
            WP_fig = (I[:, -1], W)
        else:
            WT_fig = (I[:, -1], W)
        for nm, w in W.items():
            put(f"ess_{sc}_{nm}", ess(w)); put(f"maxw_{sc}_{nm}_pct", float(100 * w.max() / w.sum()))
            put(f"{sc}_gamma2_{nm}_median", wq(gam, w, 0.5))
            for yr in [2026, 2042]:
                jj = yr - 2018
                put(f"{sc}{yr}_{nm}_median", 100 * wq(I[:, jj], w, 0.5))
                put(f"{sc}{yr}_{nm}_p05", 100 * wq(I[:, jj], w, 0.05))
                put(f"{sc}{yr}_{nm}_p95", 100 * wq(I[:, jj], w, 0.95))
        print(sc, "done", n_new, time.time() - t0, flush=True)
    # ---- is the 2022 prevalence monotone in the solved quantity? (injectivity of the map)
    g = np.random.default_rng(SEED_EXTRA)
    sample = [ok[i] for i in g.choice(len(ok), 100, replace=False)]
    mono_P = mono_T = 0
    for r in sample:
        d = {k: r[k] for k in cal.PARAM_ORDER}; p0 = Params(b=r["b0"], **d); x18 = np.array(r["x18"])
        vb = [cal.I_after(p0.with_(b=b), x18, 4.0) for b in np.linspace(r["b0"], 200, 40)]
        vd = [cal.I_after(p0, x18, 4.0, alpha_t=cal.shock(p0.alpha, z)) for z in np.linspace(0, 20, 40)]
        mono_P += bool(np.all(np.diff(vb) > 0)); mono_T += bool(np.all(np.diff(vd) > 0))
    put("mono_n", len(sample)); put("mono_P_n", mono_P); put("mono_T_n", mono_T)

    # figure
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    for a_, (xw, W), Ik, title in [(ax[0], WP_fig, PPk, "Scenario P"), (ax[1], WT_fig, PTk, "Scenario T")]:
        xk = 100 * Ik[:, -1]; x = 100 * xw
        bins = np.linspace(np.percentile(xk, 0.5), np.percentile(xk, 99.5), 40)
        a_.hist(xk, bins=bins, density=True, color=GREY, alpha=0.5, label="as calibrated")
        a_.hist(x, bins=bins, weights=W["flat"], density=True, histtype="step", lw=1.8, color=INK, label="flat prior on solved quantity")
        a_.hist(x, bins=bins, weights=W["log"], density=True, histtype="step", lw=1.8, ls="--", color=C1, label="log-uniform prior")
        a_.set_title(title, fontsize=10); a_.set_xlabel("Problematic use in 2042 (%)"); a_.grid(alpha=0.3)
    ax[0].set_ylabel("Density"); ax[0].legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "fig_bayx_weights.png", dpi=200); plt.close(fig)

    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_bayx.json", "w"), indent=1, ensure_ascii=False)
    print("done", len(R), (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
