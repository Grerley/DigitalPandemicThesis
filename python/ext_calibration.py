"""Extended calibration analysis for Chapter 7.

Writes outputs/results_calx.json (keys prefixed calx_) and figures/fig_calx_*.png.

  1. Kept against sampled parameter distributions (figure).
  2. Two further scenarios that also fit the two HBSC figures: a gradual rise in exposure between 2018 and
     2022 (Scenario G) and a shock confined to 2021 (Scenario T1).
  3. Sensitivity of the projections to how the parameter ranges are sampled.
  4. Plausibility weighting against the entry-rate bound (3.2).
  5. What the 2025/2026 HBSC wave could tell: likelihood ratios for Scenario P against Scenario T.

    python ext_calibration.py
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np
from scipy.optimize import brentq
from scipy.stats import gaussian_kde

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt
from dpm.model import Params, equilibrium, simulate
from dpm import calibrate as cal

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R = {}
SEED = ra.SEED + 2100
YEARS = np.arange(2018, 2043)


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["calx_" + k] = {"value": v, "note": note}


def q3(v):
    v = np.asarray(v)
    return float(np.median(v)), float(np.percentile(v, 5)), float(np.percentile(v, 95))


# ---------------------------------------------------------------- 1. kept vs sampled
def sec1(rows, ok):
    labels = {"eta": "η", "alpha": "α", "sigma": "σ", "gamma1": "γ₁", "gamma2": "γ₂", "delta": "δ"}
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.0))
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    for ax, k in zip(axs.flat, labels):
        lo, hi = ranges[k]["low"], ranges[k]["high"]
        logsc = ranges[k].get("log", False)
        bins = np.geomspace(lo, hi, 16) if logsc else np.linspace(lo, hi, 16)
        ax.hist([r[k] for r in rows], bins=bins, color=GREY, alpha=0.6, label="Sampled")
        ax.hist([r[k] for r in ok], bins=bins, color=C1, alpha=0.8, label="Kept")
        if logsc:
            from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter
            ax.set_xscale("log")
            ticks = [0.002, 0.01, 0.05] if k == "alpha" else [0.15, 0.3, 0.6, 1.2]
            ax.xaxis.set_major_locator(FixedLocator(ticks)); ax.xaxis.set_major_formatter(FixedFormatter([str(t) for t in ticks]))
            ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_title(labels[k], fontsize=9); ax.tick_params(labelsize=7)
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(FIG / "fig_calx_kept.png"); plt.close(fig)
    # share of kept by gamma2 tercile of its range
    g2 = np.array([r["gamma2"] for r in rows]); kept = np.array([r["ok"] for r in rows])
    edges = np.geomspace(0.15, 1.25, 4)
    for i in range(3):
        m = (g2 >= edges[i]) & (g2 < edges[i + 1] if i < 2 else g2 <= edges[i + 1])
        put(f"keep_rate_gamma2_t{i + 1}_pct", 100 * float(kept[m].mean()))
        put(f"gamma2_tercile_edge_{i + 1}", float(edges[i + 1]))


# ---------------------------------------------------------------- 2. further scenarios
def gradual_b(b0, b1):
    return lambda t: b0 + (b1 - b0) * min(max(t, 0.0), 4.0) / 4.0


def sim_gradual(p0, x18, bG, years, t_eval):
    f = gradual_b(p0.b, bG)
    from scipy.integrate import solve_ivp
    from dpm.model import rhs as mrhs
    sol = solve_ivp(lambda t, x: mrhs(t, x, p0.with_(b=f(t))), (0, years), x18, t_eval=t_eval,
                    rtol=1e-9, atol=1e-11, max_step=0.1)
    return sol.y


def solve_gradual(p0, x18, target):
    g = lambda b: sim_gradual(p0, x18, b, 4.0, [4.0])[2, -1] - target
    if g(p0.b) > 0 or g(400.0) < 0:
        return None
    return brentq(g, p0.b, 400.0, xtol=1e-9)


def shock1(alpha0, d):
    return lambda t: alpha0 + (d if 3.0 <= t < 4.0 else 0.0)


def solve_shock1(p0, x18, target):
    g = lambda d: simulate(p0, x18, 4.0, t_eval=[4.0], alpha_t=shock1(p0.alpha, d), rtol=1e-9, atol=1e-11)[1][2, -1] - target
    if g(0.0) > 0 or g(60.0) < 0:
        return None
    return brentq(g, 0.0, 60.0, xtol=1e-10)


def sec2(ok):
    t = YEARS - 2018.0
    res = {"G": [], "T1": []}; fails = {"G": 0, "T1": 0}
    for r in ok:
        d = {k: r[k] for k in cal.PARAM_ORDER}
        p0 = Params(b=r["b0"], **d); x18 = np.array(r["x18"])
        bG = solve_gradual(p0, x18, r["t22"])
        if bG is None:
            fails["G"] += 1
        else:
            y = sim_gradual(p0, x18, bG, 24.0, t)
            res["G"].append((bG, y[2]))
        d1 = solve_shock1(p0, x18, r["t22"])
        if d1 is None:
            fails["T1"] += 1
        else:
            _, y = simulate(p0, x18, 24.0, t_eval=t, alpha_t=shock1(p0.alpha, d1))
            res["T1"].append((d1, y[2]))
    for sc in res:
        put(f"{sc}_n_ok", len(res[sc])); put(f"{sc}_n_fail", fails[sc])
        par = [x[0] for x in res[sc]]
        m, lo, hi = q3(par)
        put(f"{sc}_par_median", m); put(f"{sc}_par_p05", lo); put(f"{sc}_par_p95", hi)
        I = np.array([x[1] for x in res[sc]])
        for yr in [2026, 2032, 2042]:
            m, lo, hi = q3(100 * I[:, yr - 2018])
            put(f"{sc}_{yr}_median", m); put(f"{sc}_{yr}_p05", lo); put(f"{sc}_{yr}_p95", hi)
    # figure: four scenarios
    P = np.loadtxt(OUT / "projection_P.csv", delimiter=","); T = np.loadtxt(OUT / "projection_T.csv", delimiter=",")
    G = np.array([x[1] for x in res["G"]]); T1 = np.array([x[1] for x in res["T1"]])
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for arr, col, lab, ls in [(P, C1, "P: persistent shift in 2018", "-"), (G, C3, "G: gradual rise 2018 to 2022", ":"),
                              (T, C2, "T: shock 2020 to 2021", "-."), (T1, INK, "T1: shock in 2021 only", "--")]:
        ax.plot(YEARS, np.median(arr, 0) * 100, color=col, lw=2, ls=ls, label=lab)
        ax.fill_between(YEARS, np.percentile(arr, 5, 0) * 100, np.percentile(arr, 95, 0) * 100, color=col, alpha=0.08, lw=0)
    ax.errorbar([2018, 2022], [7, 11], yerr=0.5, fmt="o", color=GREY, capsize=3)
    ax.set_xlabel("Year"); ax.set_ylabel("Problematic social media use (%)"); ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=7, loc="upper left"); fig.savefig(FIG / "fig_calx_scenarios.png"); plt.close(fig)
    return res


# ---------------------------------------------------------------- 3. sampling sensitivity
def sec3(ranges):
    variants = {
        "uniform_g2_alpha": {k: dict(v, log=False) if k in ("gamma2", "alpha") else v for k, v in ranges.items() if not k.startswith("_")},
        "narrow": {k: (dict(v, low=v["central"] - 0.5 * (v["central"] - v["low"]), high=v["central"] + 0.5 * (v["high"] - v["central"]))
                       if k != "nu" else v) for k, v in ranges.items() if not k.startswith("_")},
        "wide_g2": {k: (dict(v, low=0.10, high=2.0) if k == "gamma2" else v) for k, v in ranges.items() if not k.startswith("_")},
    }
    t = YEARS - 2018.0
    for i, (name, rg) in enumerate(variants.items()):
        rows = cal.calibrate(rg, n=300, seed=SEED + 10 * (i + 1))
        ok = [r for r in rows if r["ok"]]
        put(f"sens_{name}_n_ok", len(ok)); put(f"sens_{name}_n", len(rows))
        for sc in ["P", "T"]:
            I = []
            for r in ok:
                d = {k: r[k] for k in cal.PARAM_ORDER}
                p0 = Params(b=r["b0"], **d); x18 = np.array(r["x18"])
                if sc == "P":
                    _, y = simulate(p0.with_(b=r["b1"]), x18, 24, t_eval=[24])
                else:
                    _, y = simulate(p0, x18, 24, t_eval=[24], alpha_t=cal.shock(p0.alpha, r["d_alpha"]))
                I.append(100 * y[2, -1])
            m, lo, hi = q3(I)
            put(f"sens_{name}_{sc}_2042_median", m); put(f"sens_{name}_{sc}_2042_p05", lo); put(f"sens_{name}_{sc}_2042_p95", hi)


# ---------------------------------------------------------------- 4. plausibility weighting
def sec4(ok):
    """Weight each draw by w = 1 if the entry-rate h <= 0.4, else exp(-((h - 0.4)/s)^2 / 2), with s = 0.2 and 0.4.
    h is the hazard of leaving S in 2022 (Scenario P) or at the peak of the shock (Scenario T)."""
    P = np.loadtxt(OUT / "projection_P.csv", delimiter=","); T = np.loadtxt(OUT / "projection_T.csv", delimiter=",")
    hP, hT = [], []
    for r in ok:
        d = {k: r[k] for k in cal.PARAM_ORDER}
        y = simulate(Params(b=r["b1"], **d), np.array(r["x18"]), 4, t_eval=[4])[1][:, -1]
        hP.append(r["b1"] * (y[2] + d["eta"] * y[1]) + d["alpha"])
        p0 = Params(b=r["b0"], **d)
        y = simulate(p0, np.array(r["x18"]), 3, t_eval=[3], alpha_t=cal.shock(p0.alpha, r["d_alpha"]))[1][:, -1]
        hT.append(r["b0"] * (y[2] + d["eta"] * y[1]) + d["alpha"] + r["d_alpha"])
    hP, hT = np.array(hP), np.array(hT)
    for s in (0.2, 0.4):
        wP = np.where(hP <= 0.4, 1.0, np.exp(-0.5 * ((hP - 0.4) / s) ** 2))
        wT = np.where(hT <= 0.4, 1.0, np.exp(-0.5 * ((hT - 0.4) / s) ** 2))
        tag = f"s{int(s * 10)}"
        put(f"w_{tag}_meanP", float(wP.mean())); put(f"w_{tag}_meanT", float(wT.mean()))
        put(f"w_{tag}_ratio_P_over_T", float(wP.mean() / wT.mean()))
        put(f"w_{tag}_essP", float(wP.sum() ** 2 / (wP ** 2).sum())); put(f"w_{tag}_essT", float(wT.sum() ** 2 / (wT ** 2).sum()))
        for sc, w, arr in [("P", wP, P), ("T", wT, T)]:
            v = arr[:, 2042 - 2018] * 100
            o = np.argsort(v); cw = np.cumsum(w[o]) / w.sum()
            put(f"w_{tag}_{sc}_2042_median", float(v[o][np.searchsorted(cw, 0.5)]))
            put(f"w_{tag}_{sc}_2042_p05", float(v[o][np.searchsorted(cw, 0.05)]))
            put(f"w_{tag}_{sc}_2042_p95", float(v[o][min(np.searchsorted(cw, 0.95), len(v) - 1)]))
            g2 = np.array([r["gamma2"] for r in ok])[o]
            put(f"w_{tag}_{sc}_gamma2_median", float(g2[np.searchsorted(cw, 0.5)]))


# ---------------------------------------------------------------- 5. the next HBSC wave
def sec5():
    P = np.loadtxt(OUT / "projection_P.csv", delimiter=",")[:, 2026 - 2018] * 100
    T = np.loadtxt(OUT / "projection_T.csv", delimiter=",")[:, 2026 - 2018] * 100
    kP, kT = gaussian_kde(P), gaussian_kde(T)
    for y in [9, 10, 11, 12, 13, 14]:
        # a published whole-per-cent figure y covers [y - 0.5, y + 0.5)
        pP = kP.integrate_box_1d(y - 0.5, y + 0.5); pT = kT.integrate_box_1d(y - 0.5, y + 0.5)
        put(f"lr_{y}", float(pP / pT) if pT > 0 else None)
        put(f"pP_{y}", float(pP)); put(f"pT_{y}", float(pT))
        put(f"post_P_{y}_pct", float(100 * pP / (pP + pT)) if pP + pT > 0 else None)
    put("P2026_share_ge12_pct", float(100 * np.mean(P >= 11.5))); put("T2026_share_ge12_pct", float(100 * np.mean(T >= 11.5)))
    put("P2026_share_le10_pct", float(100 * np.mean(P < 10.5))); put("T2026_share_le10_pct", float(100 * np.mean(T < 10.5)))
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    xs = np.linspace(5, 20, 300)
    ax.plot(xs, kP(xs), color=C1, lw=2, label="Scenario P"); ax.plot(xs, kT(xs), color=C2, lw=2, label="Scenario T")
    ax.set_xlabel("Projected problematic use in 2026 (%)"); ax.set_ylabel("Density across draws")
    ax.legend(frameon=False, fontsize=8); fig.savefig(FIG / "fig_calx_2026.png"); plt.close(fig)


def main():
    t0 = time.time()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    rows = json.load(open(OUT / "calibration_draws.json")); ok = [r for r in rows if r["ok"]]
    sec1(rows, ok); print("1", time.time() - t0)
    sec2(ok); print("2", time.time() - t0)
    sec4(ok); print("4", time.time() - t0)
    sec5(); print("5", time.time() - t0)
    sec3(ranges); print("3", time.time() - t0)
    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_calx.json", "w"), indent=1, ensure_ascii=False)
    print(f"done, {len(R)} results, {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
