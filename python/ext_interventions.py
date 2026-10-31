"""Extended intervention analysis for Chapter 9.

Writes outputs/results_intx.json (keys prefixed intx_) and figures/fig_intx_*.png.

  1. Response curves: reduction in 2042 against assumed efficacy, central parameters, Scenario P.
  2. Timing: start year 2024, 2028 or 2032; reduction in 2042 and person-years of problematic use averted.
  3. Ramp-up and waning of effects.
  4. What drives the uncertainty in each effect: partial rank correlations with efficacies and model rates.
  5. Interventions on networks: platform regulation and treatment access on Barabasi-Albert graphs.

    python ext_interventions.py
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np
from scipy.integrate import solve_ivp

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt, INTERVENTIONS, apply
from dpm.model import Params, equilibrium, simulate, rhs
from dpm import calibrate as cal
from dpm import network as nw

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R = {}
SEED = ra.SEED + 2500


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["intx_" + k] = {"value": v, "note": note}


def q3(v):
    v = np.asarray(v)
    return float(np.median(v)), float(np.percentile(v, 5)), float(np.percentile(v, 95))


def state_at(p1, x18, year):
    return simulate(p1, x18, year - 2018.0, t_eval=[year - 2018.0])[1][:, -1]


def run_from(p, x0, start, end=2042, mult_t=None):
    """Integrate from `start` to `end`. mult_t(t) -> dict of parameter multipliers at time t after start."""
    t_eval = np.arange(0, end - start + 1, 1.0)
    if mult_t is None:
        return simulate(p, x0, end - start, t_eval=t_eval)[1]
    def f(t, x):
        return rhs(t, x, apply(p, mult_t(t)))
    sol = solve_ivp(f, (0, end - start), x0, t_eval=t_eval, rtol=1e-9, atol=1e-11, max_step=0.1)
    return sol.y


# ---------------------------------------------------------------- 1. response curves
def sec1(pc0, pc1, x18):
    x24 = state_at(pc1, x18, 2024)
    base = run_from(pc1, x24, 2024)[2, -1]
    curves = {
        "platform": (np.linspace(0, 0.9, 19), lambda e: {"b": 1 - e}, "Cut in exposure rate b"),
        "school": (np.linspace(0, 1.5, 16), lambda e: {"gamma1": 1 + e}, "Increase in remission γ₁"),
        "minimum_age": (np.linspace(0, 1.0, 11), lambda e: {"b": 1 - 0.4 * e, "alpha": 1 - 0.4 * e}, "Compliance with minimum age"),
        "treatment_g2": (np.linspace(0, 2.0, 21), lambda e: {"gamma2": 1 + e}, "Increase in recovery γ₂"),
        "treatment_d": (np.linspace(0, 0.9, 10), lambda e: {"delta": 1 - e}, "Cut in relapse δ"),
    }
    out = {}
    for k, (grid, mk, lab) in curves.items():
        red = [100 * (1 - run_from(apply(pc1, mk(e)), x24, 2024)[2, -1] / base) for e in grid]
        out[k] = (grid, red, lab)
        put(f"curve_{k}_grid", [float(g) for g in grid]); put(f"curve_{k}_red", [float(r) for r in red])
    # summary points used in the text
    put("curve_platform_red_at50", float(np.interp(0.5, out["platform"][0], out["platform"][1])))
    put("curve_platform_red_at90", float(out["platform"][1][-1]))
    put("curve_g2_red_at100", float(np.interp(1.0, out["treatment_g2"][0], out["treatment_g2"][1])))
    put("curve_g2_red_at200", float(out["treatment_g2"][1][-1]))
    put("curve_d_red_at50", float(np.interp(0.5, out["treatment_d"][0], out["treatment_d"][1])))
    put("curve_school_red_at100", float(np.interp(1.0, out["school"][0], out["school"][1])))
    floor = equilibrium(pc1.with_(b=0.0))[2]
    put("curve_platform_max_possible_pct", 100 * (1 - floor / base))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    ax = axs[0]
    for k, col, lst in [("platform", C1, "-"), ("minimum_age", C3, ":"), ("treatment_d", C2, "-.")]:
        g, r, lab = out[k]; ax.plot(g * 100, r, color=col, lw=2, ls=lst, label=lab)
    ax.axhline(100 * (1 - floor / base), color=GREY, lw=0.8, ls="--")
    ax.text(2, 100 * (1 - floor / base) + 1.5, "limit for exposure-only policies", fontsize=7, color=GREY)
    ax.set_xlabel("Efficacy or compliance (%)"); ax.set_ylabel("Reduction in 2042 (%)"); ax.legend(frameon=False, fontsize=7, loc="center left")
    ax = axs[1]
    for k, col in [("treatment_g2", C2), ("school", INK)]:
        g, r, lab = out[k]; ax.plot(g * 100, r, color=col, lw=2, label=lab)
    ax.set_xlabel("Increase in rate (%)"); ax.legend(frameon=False, fontsize=7)
    for ax in axs:
        ax.set_ylim(0, 105)
    fig.tight_layout(); fig.savefig(FIG / "fig_intx_curves.png"); plt.close(fig)


# ---------------------------------------------------------------- 2 and 3. timing, ramp-up, waning
CENTRAL_EFF = {"platform": {"b": 0.75}, "school": {"gamma1": 1.30}, "minimum_age": {"b": 0.8, "alpha": 0.8},
               "treatment": {"gamma2": 1.60, "delta": 0.85}}


def person_years(I_series):
    return float(np.trapezoid(I_series, dx=1.0))


def sec2(draws, start_years=(2024, 2028, 2032)):
    res = {}
    for r in draws:
        d = {k: r[k] for k in cal.PARAM_ORDER}
        p1 = Params(b=r["b1"], **d); x18 = np.array(r["x18"])
        base24 = run_from(p1, state_at(p1, x18, 2024), 2024)[2]
        for name, mult in CENTRAL_EFF.items():
            for sy in start_years:
                xs = state_at(p1, x18, sy)
                y = run_from(apply(p1, mult), xs, sy)[2]
                b_s = base24[sy - 2024:]
                red42 = 1 - y[-1] / b_s[-1]
                averted = person_years(b_s - y)            # fraction of population x years, sy to 2042
                res.setdefault((name, sy), []).append((red42, averted))
            # ramp-up over 3 years and waning with a 5-year half-life (start 2024)
            x24 = state_at(p1, x18, 2024)
            ramp = lambda t, m=mult: {k: 1 + (v - 1) * min(t / 3.0, 1.0) for k, v in m.items()}
            wane = lambda t, m=mult: {k: 1 + (v - 1) * 0.5 ** (t / 5.0) for k, v in m.items()}
            for tag, mt in [("ramp", ramp), ("wane", wane)]:
                y = run_from(p1, x24, 2024, mult_t=mt)[2]
                res.setdefault((name, tag), []).append((1 - y[-1] / base24[-1], person_years(base24 - y)))
            y = run_from(apply(p1, mult), x24, 2024)[2]
            res.setdefault((name, "full"), []).append((1 - y[-1] / base24[-1], person_years(base24 - y)))
    for (name, key), v in res.items():
        v = np.array(v)
        m, lo, hi = q3(100 * v[:, 0])
        put(f"tim_{name}_{key}_red42_median", m); put(f"tim_{name}_{key}_red42_p05", lo); put(f"tim_{name}_{key}_red42_p95", hi)
        # person-years per 1,000 young people
        m, lo, hi = q3(1000 * v[:, 1])
        put(f"tim_{name}_{key}_py_median", m); put(f"tim_{name}_{key}_py_p05", lo); put(f"tim_{name}_{key}_py_p95", hi)
    put("tim_n_draws", len(draws))
    put("central_eff_platform_cut_pct", 25.0); put("central_eff_school_up_pct", 30.0); put("central_eff_minage_compliance_pct", 50.0)
    put("central_eff_g2_up_pct", 60.0); put("central_eff_delta_cut_pct", 15.0)


# ---------------------------------------------------------------- 4. uncertainty decomposition
def prcc(X, y):
    from scipy.stats import rankdata
    Xr = np.column_stack([rankdata(c) for c in X.T]); yr = rankdata(y)
    out = []
    for j in range(Xr.shape[1]):
        others = np.delete(Xr, j, axis=1); A = np.column_stack([np.ones(len(yr)), others])
        rx = Xr[:, j] - A @ np.linalg.lstsq(A, Xr[:, j], rcond=None)[0]
        ry = yr - A @ np.linalg.lstsq(A, yr, rcond=None)[0]
        out.append(np.corrcoef(rx, ry)[0, 1])
    return np.array(out)


def sec4(draws):
    g = np.random.default_rng(SEED + 4)
    rows = {k: [] for k in INTERVENTIONS}
    for r in draws:
        d = {k: r[k] for k in cal.PARAM_ORDER}
        p1 = Params(b=r["b1"], **d); x18 = np.array(r["x18"])
        x24 = state_at(p1, x18, 2024)
        base = run_from(p1, x24, 2024)[2, -1]
        for k in INTERVENTIONS:
            m = INTERVENTIONS[k][1](g)
            y = run_from(apply(p1, m), x24, 2024)[2, -1]
            keys = ["b"] if k == "minimum_age" else sorted(m)      # b and alpha share one compliance draw
            effs = [m[x] for x in keys]
            rows[k].append(([d[x] for x in cal.PARAM_ORDER] + [r["b1"]] + effs, 1 - y / base, keys))
    names = cal.PARAM_ORDER + ["b1"]
    for k, v in rows.items():
        X = np.array([x[0] for x in v]); y = np.array([x[1] for x in v])
        pc = prcc(X, y)
        effn = ["eff_" + e for e in v[0][2]]
        for nm, val in zip(names + effn, pc):
            put(f"prcc_{k}_{nm}", float(val))
        # share of variance explained by efficacies alone (R^2 of rank regression on efficacies)
        from scipy.stats import rankdata
        Er = np.column_stack([rankdata(c) for c in X[:, len(names):].T]); yr = rankdata(y)
        A = np.column_stack([np.ones(len(yr)), Er]); res = yr - A @ np.linalg.lstsq(A, yr, rcond=None)[0]
        put(f"r2eff_{k}_pct", float(100 * (1 - res.var() / yr.var())))
        Pr = np.column_stack([rankdata(c) for c in X[:, :len(names)].T])
        A = np.column_stack([np.ones(len(yr)), Pr]); res = yr - A @ np.linalg.lstsq(A, yr, rcond=None)[0]
        put(f"r2par_{k}_pct", float(100 * (1 - res.var() / yr.var())))
    put("prcc_n", len(draws))


# ---------------------------------------------------------------- 5. networks
def sec5(pc1, x18, reps=30, n=2000, k=10):
    _, y = simulate(pc1, x18, 6, t_eval=[6]); x24 = y[:, -1]
    beta = nw.beta_from_b(pc1, k)
    months = 18 * 12
    out = {}
    for name, mult in [("none", {}), ("platform", CENTRAL_EFF["platform"]), ("treatment", CENTRAL_EFF["treatment"])]:
        v = []
        for s in range(reps):
            G = nw.make_graph("BA", n, k, seed=SEED + 100 + s)
            p = apply(pc1, mult) if mult else pc1
            bta = nw.beta_from_b(p, k)
            traj = nw.simulate(G, p, bta, months, x24, seed=SEED + 300 + s, record_every=months)
            v.append(traj[-1, 2])
        out[name] = np.array(v)
    for name in ["platform", "treatment"]:
        red = 100 * (1 - out[name] / out["none"])
        put(f"net_{name}_red_mean", float(np.mean(1 - out[name].mean() / out["none"].mean()) * 100))
        put(f"net_{name}_red_se", float(red.std(ddof=1) / math.sqrt(reps)))
    # mean-field counterparts with the same central efficacies
    for name in ["platform", "treatment"]:
        yb = run_from(pc1, x24, 2024)[2, -1]; yi = run_from(apply(pc1, CENTRAL_EFF[name]), x24, 2024)[2, -1]
        put(f"mf_{name}_red", 100 * (1 - yi / yb))
    put("net_reps", reps); put("net_nodes", n)


def main():
    t0 = time.time()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)
    ok = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]]
    sec1(pc0, pc1, x18); print("1", time.time() - t0)
    sec2(ok[:150]); print("2", time.time() - t0)
    sec4(ok[:400]); print("4", time.time() - t0)
    sec5(pc1, x18); print("5", time.time() - t0)
    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_intx.json", "w"), indent=1, ensure_ascii=False)
    print(f"done, {len(R)} results, {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
