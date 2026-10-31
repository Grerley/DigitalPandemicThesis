"""Age- and sex-structured model: calibration, projections and the minimum-age scenario.

Writes outputs/results_agesex.json (keys prefixed asx_) and figures/fig_asx_*.png.

Calibration targets (Boniel-Nissim et al., 2024): HBSC international average problematic social media use
7% in 2018 (pp. vi, 6); in 2022 girls 13%, boys 9%, 11-year-olds 9% (p. 7). Each is rounded, so each target is
drawn within +/- 0.5 percentage points. The HBSC figures for 13-year-olds (12%) and the 2022 total (11%) are
not fitted and serve as checks. Four quantities are solved for each draw:
  Scenario P: b0 (exposure before 2018), b1 (from 2018), phi (girls' vulnerability), psi (10-to-12 vulnerability).
  Scenario T: b0, d_alpha (2020-2021 shock), phi, psi.

    python ext_agesex.py
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np
from scipy.optimize import least_squares

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt
from dpm.model import Params, simulate as sim_h
from dpm import agesex as AS
from dpm import calibrate as cal

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R = {}
SEED = ra.SEED + 1300
START = 2024


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["asx_" + k] = {"value": v, "note": note}


def pct(x):
    return 100.0 * x


def shock(alpha0, da):
    return lambda t: alpha0 + (da if 2.0 <= t < 4.0 else 0.0)


def fit(p_base: Params, st0: AS.Struct, targets, scenario="P", z0=None):
    """Solve the four unknowns by least squares on log scale. Returns dict or None."""
    t18, tg, tb, ta = targets

    def unpack(z):
        b0, k2, phi, psi = np.exp(z)
        return b0, k2, phi, psi

    def resid(z):
        b0, k2, phi, psi = unpack(z)
        st = st0.with_(phi=phi, psi=psi)
        p0 = p_base.with_(b=b0)
        x18 = AS.equilibrium(p0, st)
        o18 = AS.hbsc_obs(x18)
        if scenario == "P":
            _, y = AS.simulate(p0.with_(b=k2), st, x18, 4, t_eval=[4])
        else:
            _, y = AS.simulate(p0, st, x18, 4, t_eval=[2, 4], alpha_t=shock(p0.alpha, k2))
        o22 = AS.hbsc_obs(y[:, -1])
        return np.array([o18["total"] - t18, o22["girls"] - tg, o22["boys"] - tb, o22["age11"] - ta]) * 100

    if z0 is None:
        z0 = np.log([1.0, 2.0 if scenario == "P" else 0.4, 1.4, 1.2])
    res = least_squares(resid, z0, method="lm", xtol=1e-12, ftol=1e-12, max_nfev=400)
    if not res.success or np.max(np.abs(res.fun)) > 1e-4:
        return None
    b0, k2, phi, psi = unpack(res.x)
    if scenario == "P" and k2 < b0:
        pass                                   # allowed: reported, not imposed
    return dict(b0=b0, k2=k2, phi=phi, psi=psi, z=res.x)


def run_scenario(p_base, st, sol, scenario, years=24, access_t=None, t_eval=None):
    p0 = p_base.with_(b=sol["b0"])
    x18 = AS.equilibrium(p0, st)
    if t_eval is None:
        t_eval = np.arange(0, years + 1, 1.0)
    if scenario == "P":
        return AS.simulate(p0.with_(b=sol["k2"]), st, x18, years, t_eval=t_eval, access_t=access_t)
    return AS.simulate(p0, st, x18, years, t_eval=t_eval, alpha_t=shock(p0.alpha, sol["k2"]), access_t=access_t)


def obs_series(y):
    out = {k: [] for k in ["total", "girls", "boys", "age11", "age13_15", "whole"]}
    for j in range(y.shape[1]):
        o = AS.hbsc_obs(y[:, j])
        for k in ["total", "girls", "boys", "age11", "age13_15"]:
            out[k].append(o[k])
        out["whole"].append(AS.whole_prev(y[:, j]))
    return {k: np.array(v) for k, v in out.items()}


def band_prev(x, a):
    yv = AS.group_prev(x)
    return 0.5 * (yv[AS.gid(0, a)] + yv[AS.gid(1, a)])


# ---------------------------------------------------------------- central point
def central_analysis(pc_base, st):
    targets = (0.07, 0.13, 0.09, 0.09)
    out = {}
    for sc in ["P", "T"]:
        sol = fit(pc_base, st, targets, sc)
        out[sc] = sol
        stc = st.with_(phi=sol["phi"], psi=sol["psi"])
        put(f"c_{sc}_b0", sol["b0"]); put(f"c_{sc}_k2", sol["k2"]); put(f"c_{sc}_phi", sol["phi"]); put(f"c_{sc}_psi", sol["psi"])
        put(f"c_{sc}_R0_b0", AS.R0_struct(pc_base.with_(b=sol["b0"], alpha=0.0), stc))
        if sc == "P":
            put("c_P_R0_b1", AS.R0_struct(pc_base.with_(b=sol["k2"], alpha=0.0), stc))
        t, y = run_scenario(pc_base, stc, sol, sc)
        s = obs_series(y)
        yrs = 2018 + t
        for yr in [2022, 2026, 2032, 2042]:
            j = int(yr - 2018)
            for k in s:
                put(f"c_{sc}_{k}_{yr}", pct(s[k][j]))
            for a, nm in enumerate(["b1019", "b1315", "b1619", "b2024"]):
                put(f"c_{sc}_{nm}_{yr}", pct(band_prev(y[:, j], a)))
            yv = AS.group_prev(y[:, j])
            put(f"c_{sc}_girls_whole_{yr}", pct(np.sum(yv[:4] * AS.N_G[:4]) / AS.N_G[:4].sum()))
            put(f"c_{sc}_boys_whole_{yr}", pct(np.sum(yv[4:] * AS.N_G[4:]) / AS.N_G[4:].sum()))
        out[sc + "_series"] = (yrs, s, y)
    # homogeneous model at the same central rates, for comparison (whole population, 2042)
    pc0 = pc_base.with_(b=cal.solve_b0(pc_base, 0.07)); x18 = ra.equilibrium(pc0)
    b1 = cal.solve_b1(pc0, x18, 0.11)
    _, yh = sim_h(pc0.with_(b=b1), x18, 24, t_eval=[24])
    put("c_homog_P_2042", pct(yh[2, -1]))
    # figure: HBSC-age series by sex and the 2022 targets
    yrs, s, _ = out["P_series"]
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.2), sharey=True)
    ax = axs[0]
    ax.plot(yrs, s["girls"] * 100, color=C2, lw=2, ls="-", label="Girls (HBSC ages)")
    ax.plot(yrs, s["boys"] * 100, color=C1, lw=2, ls="-.", label="Boys (HBSC ages)")
    ax.plot(yrs, s["whole"] * 100, color=INK, lw=1, ls="--", label="Whole population 10 to 24")
    ax.errorbar([2022, 2022], [13, 9], yerr=0.5, fmt="o", color=GREY, capsize=3)
    ax.errorbar([2018], [7], yerr=0.5, fmt="s", color=GREY, capsize=3)
    ax.set_title("Scenario P, by sex", fontsize=9); ax.set_ylabel("Problematic social media use (%)"); ax.legend(frameon=False, fontsize=7)
    ax = axs[1]
    _, _, y = out["P_series"]
    for a, col, lab in [(0, C3, "10 to 12"), (1, C1, "13 to 15"), (2, C2, "16 to 19"), (3, INK, "20 to 24")]:
        ax.plot(yrs, [band_prev(y[:, j], a) * 100 for j in range(y.shape[1])], color=col, lw=2, label=lab)
    ax.errorbar([2022], [9], yerr=0.5, fmt="o", color=GREY, capsize=3)
    ax.set_title("Scenario P, by age band (both sexes)", fontsize=9); ax.legend(frameon=False, fontsize=7, title="Age band", title_fontsize=7)
    for ax in axs:
        ax.set_xlabel("Year"); ax.set_ylim(0, None); ax.set_xticks([2018, 2022, 2026, 2030, 2034, 2038, 2042])
        ax.tick_params(axis="x", labelsize=7)
    fig.tight_layout(); fig.savefig(FIG / "fig_asx_central.png"); plt.close(fig)
    return out


# ---------------------------------------------------------------- uncertainty across draws
def draws_analysis(st, n):
    g = np.random.default_rng(SEED)
    ok = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]][:n]
    rows = {"P": [], "T": []}
    fails = {"P": 0, "T": 0}
    for r in ok:
        d = {k: r[k] for k in cal.PARAM_ORDER}
        base = Params(b=0.0, **d)
        stx = st.with_(eps_a=g.uniform(0.5, 0.9), eps_s=g.uniform(0.3, 0.7))
        targets = (g.uniform(0.065, 0.075), g.uniform(0.125, 0.135), g.uniform(0.085, 0.095), g.uniform(0.085, 0.095))
        for sc in ["P", "T"]:
            sol = fit(base, stx, targets, sc)
            if sol is None:
                fails[sc] += 1; continue
            stc = stx.with_(phi=sol["phi"], psi=sol["psi"])
            t, y = run_scenario(base, stc, sol, sc)
            s = obs_series(y)
            rows[sc].append(dict(sol=sol, base=base, st=stx, phi=sol["phi"], psi=sol["psi"], b0=sol["b0"], k2=sol["k2"],
                                 total42=s["total"][-1], girls42=s["girls"][-1], boys42=s["boys"][-1],
                                 whole42=s["whole"][-1], age13_22=s["age13_15"][4], eps_a=stx.eps_a, eps_s=stx.eps_s))
    put("d_n", len(ok))
    for sc in ["P", "T"]:
        put(f"d_{sc}_n_ok", len(rows[sc])); put(f"d_{sc}_n_fail", fails[sc])
        for k in ["phi", "psi", "total42", "girls42", "boys42", "whole42", "age13_22"]:
            v = np.array([x[k] for x in rows[sc]])
            scale = 100.0 if k not in ("phi", "psi") else 1.0
            put(f"d_{sc}_{k}_median", scale * float(np.median(v)))
            put(f"d_{sc}_{k}_p05", scale * float(np.percentile(v, 5)))
            put(f"d_{sc}_{k}_p95", scale * float(np.percentile(v, 95)))
        gap = np.array([x["girls42"] - x["boys42"] for x in rows[sc]]) * 100
        put(f"d_{sc}_gap42_median", float(np.median(gap)))
        put(f"d_{sc}_gap42_p05", float(np.percentile(gap, 5))); put(f"d_{sc}_gap42_p95", float(np.percentile(gap, 95)))
        # how strongly do the assumed mixing parameters move phi and psi? rank correlations
        from scipy.stats import spearmanr
        for k in ["phi", "psi"]:
            for e in ["eps_a", "eps_s"]:
                put(f"d_{sc}_rho_{k}_{e}", float(spearmanr([x[k] for x in rows[sc]], [x[e] for x in rows[sc]])[0]))
    return rows


# ---------------------------------------------------------------- minimum age of 16, structured
def min_age_access(compliance, start=START - 2018):
    """From the start year a share `compliance` of young people in the 10-12 and 13-15 bands is off social media."""
    def f(t):
        m = np.ones(AS.G)
        if t >= start:
            for s in range(2):
                for a in (0, 1):
                    m[AS.gid(s, a)] = 1 - compliance
        return m
    return f


def run_min_age(p_base, st, sol, compliance):
    """Scenario P with the policy from START; integration split at the policy start."""
    stc = st.with_(phi=sol["phi"], psi=sol["psi"])
    p0 = p_base.with_(b=sol["b0"]); x18 = AS.equilibrium(p0, stc)
    p1 = p0.with_(b=sol["k2"])
    _, y = AS.simulate(p1, stc, x18, START - 2018, t_eval=[START - 2018])
    x24 = y[:, -1]
    acc = min_age_access(compliance, start=0.0)
    _, y2 = AS.simulate(p1, stc, x24, 2042 - START, t_eval=np.arange(0, 2042 - START + 1, 1.0), access_t=acc)
    _, y0 = AS.simulate(p1, stc, x24, 2042 - START, t_eval=np.arange(0, 2042 - START + 1, 1.0))
    return y0, y2


def min_age_analysis(pc_base, st, central, rows_n):
    sol = central["P"]
    out = {}
    for c in [0.3, 0.5, 0.7]:
        y0, y2 = run_min_age(pc_base, st, sol, c)
        x0, x2 = y0[:, -1], y2[:, -1]
        tag = f"c{int(c * 100)}"
        put(f"ma_{tag}_whole_red_pct", 100 * (1 - AS.whole_prev(x2) / AS.whole_prev(x0)))
        put(f"ma_{tag}_hbsc_red_pct", 100 * (1 - AS.hbsc_obs(x2)["total"] / AS.hbsc_obs(x0)["total"]))
        for a, nm in enumerate(["b1012", "b1315", "b1619", "b2024"]):
            put(f"ma_{tag}_{nm}_red_pct", 100 * (1 - band_prev(x2, a) / band_prev(x0, a)))
        out[c] = (y0, y2)
    # figure: reduction by age band over time at 50% compliance
    y0, y2 = out[0.5]
    yrs = START + np.arange(y0.shape[1])
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    for a, col, lab in [(0, C3, "10 to 12"), (1, C1, "13 to 15"), (2, C2, "16 to 19"), (3, INK, "20 to 24")]:
        red = [100 * (1 - band_prev(y2[:, j], a) / band_prev(y0[:, j], a)) for j in range(y0.shape[1])]
        ax.plot(yrs, red, color=col, lw=2, label=lab)
    whole = [100 * (1 - AS.whole_prev(y2[:, j]) / AS.whole_prev(y0[:, j])) for j in range(y0.shape[1])]
    ax.plot(yrs, whole, color=GREY, lw=1.5, ls="--", label="Whole population")
    ax.set_xlabel("Year"); ax.set_ylabel("Reduction in problematic use (%)"); ax.legend(frameon=False, fontsize=7, title="Age band", title_fontsize=7)
    ax.set_xticks(np.arange(START, 2043, 3))
    fig.savefig(FIG / "fig_asx_minage.png"); plt.close(fig)
    # across draws (reusing the Scenario P fits), compliance drawn as in Chapter 9 (30% to 70%)
    g = np.random.default_rng(SEED + 5)
    red_w, red_h = [], []
    for x in rows_n["P"]:
        c = g.uniform(0.3, 0.7)
        y0, y2 = run_min_age(x["base"], x["st"], x["sol"], c)
        red_w.append(100 * (1 - AS.whole_prev(y2[:, -1]) / AS.whole_prev(y0[:, -1])))
        red_h.append(100 * (1 - AS.hbsc_obs(y2[:, -1])["total"] / AS.hbsc_obs(y0[:, -1])["total"]))
    for nm, v in [("whole", red_w), ("hbsc", red_h)]:
        v = np.array(v)
        put(f"ma_d_{nm}_median", float(np.median(v))); put(f"ma_d_{nm}_p05", float(np.percentile(v, 5)))
        put(f"ma_d_{nm}_p95", float(np.percentile(v, 95)))
    put("ma_d_n", len(red_w))


def checks(pc_base, st, central):
    """Numerical checks reported in Chapter 5: population conservation, uniqueness of the equilibrium from
    random starts, R0 threshold behaviour."""
    sol = central["P"]
    stc = st.with_(phi=sol["phi"], psi=sol["psi"])
    p1 = pc_base.with_(b=sol["k2"])
    g = np.random.default_rng(SEED + 9)
    xs = []
    for _ in range(40):
        X = np.zeros((AS.G, AS.K))
        for gg in range(AS.G):
            X[gg] = g.dirichlet(np.ones(4)) * AS.N_G[gg]
        _, y = AS.simulate(p1, stc, X.ravel(), 400, t_eval=[400])
        xs.append(y[:, -1])
    xs = np.array(xs)
    put("chk_unique_maxdev", float(np.max(np.abs(xs - xs.mean(0)))))
    put("chk_mass_maxdev", float(np.max(np.abs(xs.sum(1) - 1))))
    # alpha = 0 threshold: R0_struct slightly below and above one
    q = pc_base.with_(alpha=0.0)
    from scipy.optimize import brentq
    bstar = brentq(lambda b: AS.R0_struct(q.with_(b=b), stc) - 1, 1e-4, 50)
    put("chk_bstar", bstar)
    for tag, f in [("below", 0.9), ("above", 1.1)]:
        _, y = AS.simulate(q.with_(b=bstar * f), stc, AS.initial_state(0.05), 600, t_eval=[600])
        put(f"chk_{tag}_whole_pct", pct(AS.whole_prev(y[:, -1])))


def main():
    t0 = time.time()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)
    base = pc1.with_(b=0.0)
    st = AS.Struct()
    put("eps_a_central", st.eps_a); put("eps_s_central", st.eps_s); put("start_year", START)
    central = central_analysis(base, st); print("central", time.time() - t0)
    checks(base, st, central); print("checks", time.time() - t0)
    rows = draws_analysis(st, 100); print("draws", time.time() - t0)
    min_age_analysis(base, st, central, rows); print("minage", time.time() - t0)
    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_agesex.json", "w"), indent=1, ensure_ascii=False)
    print(f"done, {len(R)} results, {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
