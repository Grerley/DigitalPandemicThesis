"""Extended mathematical analysis for Chapters 4 and 5 (thesis expansion).

Writes outputs/results_mathx.json (keys prefixed mathx_) and figures/fig_mathx_*.png.

Sections
  1. Closed-form sensitivity indices of R0 for every parameter, checked against numerical derivatives.
  2. Closed-form sensitivity indices of the endemic prevalence I* (alpha = 0) and implicit ones (alpha > 0).
  3. Time scales: eigenvalues at the endemic equilibrium; approach times across calibrated draws.
  4. Stochastic invasion: multitype branching-process extinction probability, checked by Gillespie simulation.
  5. Stochastic fluctuations with alpha > 0: linear noise approximation, checked by simulation.
  6. Profile likelihoods for an idealised twenty-year annual series.
  7. Two-parameter maps: equilibrium prevalence over (b, gamma2) and the floor over (alpha, gamma2).
  8. Structural alternatives: simpler models calibrated to the same two figures, and their projections.

    python ext_math.py
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np
from scipy.optimize import brentq, least_squares
from scipy.linalg import solve_continuous_lyapunov, eigvals
from scipy.integrate import solve_ivp

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt
from dpm.model import Params, R0, equilibrium, simulate, jacobian, elasticity, relapse_multiplier, rhs as rhs_m
from dpm import calibrate as cal

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R = {}
SEED = ra.SEED + 900


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["mathx_" + k] = {"value": v, "note": note}


# ---------------------------------------------------------------- 1. R0 sensitivity indices in closed form
def blocks(p):
    a = p.sigma + p.gamma1 + p.nu; c = p.gamma2 + p.nu; d = p.delta + p.nu
    D = a * c * d - p.sigma * p.gamma2 * p.delta
    return a, c, d, D


def R0_indices(p):
    """Closed-form normalised sensitivity indices of R0 = b d (eta c + sigma) / D."""
    a, c, d, D = blocks(p)
    e, s, g1, g2, de, nu = p.eta, p.sigma, p.gamma1, p.gamma2, p.delta, p.nu
    return {
        "b": 1.0,
        "eta": e * c / (e * c + s),
        "sigma": s * (1 / (e * c + s) - (c * d - g2 * de) / D),
        "gamma1": -g1 * c * d / D,
        "gamma2": g2 * (e / (e * c + s) - (a * d - s * de) / D),
        "delta": s * g2 * de * nu / (d * D),
        "nu": nu * (1 / d + e / (e * c + s) - (c * d + a * d + a * c) / D),
    }


def sec1(pc1):
    q = pc1.with_(alpha=0.0)
    cf = R0_indices(q)
    worst = 0.0
    for k, v in cf.items():
        num = elasticity(lambda p: R0(p.with_(alpha=0.0)), q, k)
        worst = max(worst, abs(num - v))
        put(f"R0idx_{k}", v)
    put("R0idx_max_abs_error", worst, "closed form against central differences")
    # random-parameter check of the closed forms and of the sign of the gamma2 index
    g = np.random.default_rng(SEED)
    worst = 0.0; neg = 0; n = 2000
    for _ in range(n):
        p = Params(b=g.uniform(0.2, 6), eta=g.uniform(0, 1), alpha=0.0, sigma=g.uniform(0.05, 1),
                   gamma1=g.uniform(0.05, 1), gamma2=g.uniform(0.1, 1.5), delta=g.uniform(0.01, 0.5),
                   nu=g.uniform(0.03, 0.1))
        cf = R0_indices(p)
        for k in ["sigma", "gamma2", "nu"]:
            worst = max(worst, abs(elasticity(R0, p, k) - cf[k]))
        neg += cf["gamma2"] < 0
    put("R0idx_random_n", n); put("R0idx_random_max_abs_error", worst); put("R0idx_gamma2_negative_share", neg / n)


# ---------------------------------------------------------------- 2. indices of I*
def Istar_indices_alpha0(p):
    """I* = (sigma/c)(1 - 1/R0)/m, so Y^I = Y^(sigma/c) - Y^m + Y^R0 / (R0 - 1)."""
    a, c, d, D = blocks(p)
    s, g2, de, nu = p.sigma, p.gamma2, p.delta, p.nu
    m = 1 + s / c + s * g2 / (c * d)
    r0 = R0(p)
    ri = R0_indices(p)
    # derivatives of m
    dm = {"b": 0.0, "eta": 0.0, "gamma1": 0.0,
          "sigma": 1 / c + g2 / (c * d),
          "gamma2": -s / c ** 2 + s * (c * d - g2 * d) / (c * d) ** 2,
          "delta": -s * g2 / (c * d ** 2),
          "nu": -s / c ** 2 - s * g2 * (c + d) / (c * d) ** 2}
    dlog_sc = {"b": 0.0, "eta": 0.0, "gamma1": 0.0, "delta": 0.0,
               "sigma": 1 / s, "gamma2": -1 / c, "nu": -1 / c}
    th = {"b": p.b, "eta": p.eta, "sigma": s, "gamma1": p.gamma1, "gamma2": g2, "delta": de, "nu": nu}
    return {k: th[k] * dlog_sc[k] - th[k] * dm[k] / m + ri[k] / (r0 - 1) for k in th}


def sec2(pc1):
    q = pc1.with_(alpha=0.0)
    cf = Istar_indices_alpha0(q)
    worst = 0.0
    for k, v in cf.items():
        num = elasticity(lambda p: equilibrium(p)[2], q, k)
        worst = max(worst, abs(num - v))
        put(f"Iidx0_{k}", v)
    put("Iidx0_max_abs_error", worst)
    # alpha > 0: compare the indices with alpha = 0 at the same point
    for k in ["b", "gamma2", "sigma", "delta", "alpha"]:
        put(f"Iidx_{k}", elasticity(lambda p: equilibrium(p)[2], pc1, k))
    # amplification: Y^I_b = Y^R0_b / (R0 - 1) when alpha = 0
    put("Iidx0_b_formula", 1 / (R0(q) - 1))


# ---------------------------------------------------------------- 3. time scales
def reduced_jac(p, x):
    J = jacobian(p, x)
    # eliminate S = 1 - A - I - R
    T = np.array([[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    return (J @ T)[1:, :]


def sec3(pc0, pc1, ok):
    for tag, p in [("b0", pc0), ("b1", pc1)]:
        x = equilibrium(p)
        ev = np.linalg.eigvals(reduced_jac(p, x))
        ev = ev[np.argsort(-ev.real)]
        slow = -ev[0].real
        put(f"eig_{tag}_re", [float(e.real) for e in ev]); put(f"eig_{tag}_im", [float(e.imag) for e in ev])
        put(f"eig_{tag}_slowest_rate", slow); put(f"eig_{tag}_halflife_years", math.log(2) / slow)
        put(f"eig_{tag}_complex", bool(np.any(np.abs(ev.imag) > 1e-9)))
    # critical slowing down near the threshold (alpha = 0): slowest rate against R0
    rows = []
    q = pc1.with_(alpha=0.0)
    for r0t in [1.05, 1.2, 1.5, 2.0, 3.0, 4.0]:
        bt = brentq(lambda b: R0(q.with_(b=b)) - r0t, 1e-6, 100)
        x = equilibrium(q.with_(b=bt))
        ev = np.linalg.eigvals(reduced_jac(q.with_(b=bt), x))
        rows.append((r0t, -max(ev.real)))
    put("csd_R0", [r[0] for r in rows]); put("csd_rate", [r[1] for r in rows])
    # approach to the post-2018 equilibrium across calibrated draws (Scenario P)
    P = np.loadtxt(OUT / "projection_P.csv", delimiter=",")
    t90, peakyr = [], []
    for i, r in enumerate(ok):
        d = {k: r[k] for k in cal.PARAM_ORDER}
        Ist = equilibrium(Params(b=r["b1"], **d))[2]
        I18 = r["x18"][2]
        _, y = simulate(Params(b=r["b1"], **d), np.array(r["x18"]), 60, t_eval=np.linspace(0, 60, 601))
        frac = (y[2] - I18) / (Ist - I18)
        # first time the trajectory stays within 10% of the remaining gap
        within = np.abs(1 - frac) <= 0.10
        idx = next((j for j in range(len(within)) if within[j:].all()), len(within) - 1)
        t90.append(idx / 10.0)
        peakyr.append(2018 + np.argmax(y[2]) / 10.0)
    t90 = np.array(t90)
    put("t90_median", float(np.median(t90))); put("t90_p05", float(np.percentile(t90, 5))); put("t90_p95", float(np.percentile(t90, 95)))
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot([r[0] for r in rows], [r[1] for r in rows], "o-", color=C1)
    ax.set_xlabel("$\\mathcal{R}_0$ (exposure rate varied, $\\alpha = 0$, other rates central)")
    ax.set_ylabel("Slowest return rate (per year)"); ax.set_ylim(0, None)
    fig.savefig(FIG / "fig_mathx_slowing.png"); plt.close(fig)


# ---------------------------------------------------------------- 4. branching process and Gillespie
def extinction_probs(p):
    """Minimal fixed point of the offspring generating functions of the three-type branching process
    (types A, I, R) linearised at the problem-free state, alpha = 0."""
    a, c, d, _ = blocks(p)
    bA, bI = p.b * p.eta, p.b

    def G(q):
        qA, qI, qR = q
        return np.array([(bA * qA * qA + p.sigma * qI + (p.gamma1 + p.nu)) / (bA + a),
                         (bI * qI * qA + p.gamma2 * qR + p.nu) / (bI + c),
                         (p.delta * qA + p.nu) / d])
    q = np.zeros(3)
    for _ in range(200000):
        qn = G(q)
        if np.max(np.abs(qn - q)) < 1e-14:
            break
        q = qn
    return q


def gillespie_invasion(p, N, start, reps, threshold, rng):
    """Exact stochastic simulation of the individual-level model (alpha = 0) from `start` = (A, I, R) counts.
    Returns the share of runs in which the number ever having entered A reaches `threshold` (a major outbreak)."""
    big = 0
    for _ in range(reps):
        nA, nI, nR = start; nS = N - nA - nI - nR; ever = nA + nI
        while nA + nI + nR > 0 and ever < threshold:
            lam = p.b * (nI + p.eta * nA) / N
            rates = np.array([lam * nS, p.sigma * nA, (p.gamma1 + p.nu) * nA, p.gamma2 * nI, p.nu * nI,
                              p.delta * nR, p.nu * nR])
            tot = rates.sum()
            k = np.searchsorted(np.cumsum(rates), rng.random() * tot)
            if k == 0: nS -= 1; nA += 1; ever += 1
            elif k == 1: nA -= 1; nI += 1
            elif k == 2: nA -= 1; nS += 1
            elif k == 3: nI -= 1; nR += 1
            elif k == 4: nI -= 1; nS += 1
            elif k == 5: nR -= 1; nA += 1
            else: nR -= 1; nS += 1
        big += ever >= threshold
    return big / reps


def sec4(pc0, pc1):
    rng = np.random.default_rng(SEED + 1)
    for tag, p in [("b0", pc0), ("b1", pc1)]:
        q = p.with_(alpha=0.0)
        ext = extinction_probs(q)
        put(f"ext_{tag}_qA", float(ext[0])); put(f"ext_{tag}_qI", float(ext[1])); put(f"ext_{tag}_qR", float(ext[2]))
        put(f"ext_{tag}_major_from_A", 1 - float(ext[0])); put(f"ext_{tag}_major_from_I", 1 - float(ext[1]))
        reps = 2000
        sim = gillespie_invasion(q, 20000, (1, 0, 0), reps, 200, rng)
        put(f"ext_{tag}_sim_major_from_A", sim); put(f"ext_{tag}_sim_se", math.sqrt(sim * (1 - sim) / reps))
        put(f"ext_{tag}_R0", R0(q))
    put("ext_sim_N", 20000); put("ext_sim_threshold", 200); put("ext_sim_reps", 2000)
    # extinction probability against R0
    q = pc1.with_(alpha=0.0)
    r0s = np.linspace(1.01, 5, 40); pm = []
    for r0t in r0s:
        bt = brentq(lambda b: R0(q.with_(b=b)) - r0t, 1e-6, 100)
        pm.append(1 - extinction_probs(q.with_(b=bt))[0])
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot(r0s, pm, color=C1, lw=2, label="Starting from one At-Risk young person")
    ax.plot(r0s, 1 - 1 / r0s, color=GREY, lw=1, ls="--", label="$1 - 1/\\mathcal{R}_0$ (single-type reference)")
    ax.set_xlabel("$\\mathcal{R}_0$"); ax.set_ylabel("Probability of a major outbreak"); ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8); fig.savefig(FIG / "fig_mathx_invasion.png"); plt.close(fig)


# ---------------------------------------------------------------- 5. linear noise approximation
def stoich_and_rates(p, x):
    S, A, I, Rr = x
    lam = p.b * (I + p.eta * A)
    # transitions: S->A, A->I, A->S (remission), I->R, R->A, turnover A->S, I->S, R->S
    v = np.array([[-1, 1, 0, 0], [0, -1, 1, 0], [1, -1, 0, 0], [0, 0, -1, 1], [0, 1, 0, -1],
                  [1, -1, 0, 0], [1, 0, -1, 0], [1, 0, 0, -1]], float)
    w = np.array([(lam + p.alpha) * S, p.sigma * A, p.gamma1 * A, p.gamma2 * I, p.delta * Rr,
                  p.nu * A, p.nu * I, p.nu * Rr])
    return v, w


def lna_sd_I(p):
    """Stationary covariance of sqrt(N)(x - x*) from J Sigma + Sigma J' + B = 0 on (A, I, R)."""
    x = equilibrium(p)
    Jr = reduced_jac(p, x)
    v, w = stoich_and_rates(p, x)
    vr = v[:, 1:]
    B = (vr.T * w) @ vr
    Sig = solve_continuous_lyapunov(Jr, -B)
    return math.sqrt(Sig[1, 1]), x


def tau_leap(p, N, years, dt, x0, rng):
    n = np.round(np.asarray(x0) * N).astype(np.int64); n[0] = N - n[1:].sum()
    steps = int(years / dt); out = np.empty(steps)
    for t in range(steps):
        v, w = stoich_and_rates(p, n / N)
        k = rng.poisson(np.maximum(w * N * dt, 0))
        src = np.argmax(v < 0, axis=1)
        for j in range(len(k)):                        # never move more than are present
            k[j] = min(k[j], n[src[j]])
            n = n + (v[j] * k[j]).astype(np.int64)
        out[t] = n[2] / N
    return out


def sec5(pc1):
    sd, x = lna_sd_I(pc1)
    put("lna_sd_scaled", sd)
    rng = np.random.default_rng(SEED + 2)
    for N in [1000, 10000, 100000]:
        put(f"lna_N{N}_sd_pp", 100 * sd / math.sqrt(N))
        put(f"lna_N{N}_cv_pct", 100 * sd / math.sqrt(N) / x[2])
    for N in [1000, 10000]:
        traj = tau_leap(pc1, N, 400, 1 / 52, x, rng)[52 * 50:]     # drop 50 years of burn-in
        put(f"lna_N{N}_sim_sd_pp", 100 * float(traj.std())); put(f"lna_N{N}_sim_mean_pct", 100 * float(traj.mean()))
    put("lna_I_star_pct", 100 * x[2])


# ---------------------------------------------------------------- 6. profile likelihood, idealised annual series
NAMES = ["b", "eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"]


def sec6(pc1, x18):
    t = np.arange(0, 21, 1.0)
    def traj(p):
        sol = solve_ivp(lambda tt, x: rhs_m(tt, x, p), (0, 20), x18, t_eval=t, method="LSODA", rtol=1e-8, atol=1e-10)
        return sol.y[2]
    data = traj(pc1)                     # expected data; profiles of the expected likelihood
    sig = 0.005
    th0 = np.array([getattr(pc1, k) for k in NAMES])

    def sse(logth):
        p = Params(**dict(zip(NAMES, np.exp(logth))))
        return (traj(p) - data) / sig

    grid = np.linspace(-1.2, 1.2, 9)          # log multiplier of the profiled parameter
    prof = {}
    mid = len(grid) // 2
    hits = [0]
    for j, nm in enumerate(NAMES):
        free = [i for i in range(len(NAMES)) if i != j]
        # free parameters may move at most a factor of 20 from their central values, which keeps the
        # optimiser away from stiff regions without binding at any reported solution
        lo = np.log(th0)[free] - 3.0; hi = np.log(th0)[free] + 3.0
        if 1 in free:
            hi[free.index(1)] = min(hi[free.index(1)], 0.0)   # eta <= 1
        vals = [None] * len(grid)
        for order in (range(mid, len(grid)), range(mid, -1, -1)):   # walk outwards, warm-starting
            z = np.minimum(np.log(th0)[free], hi - 1e-9)
            for gi in order:
                fixed = math.log(th0[j]) + grid[gi]
                if nm == "eta" and fixed > 0:
                    vals[gi] = float("nan"); continue
                def resid(zz):
                    full = np.empty(len(NAMES)); full[free] = zz; full[j] = fixed
                    return sse(full)
                res = least_squares(resid, z, method="trf", bounds=(lo, hi), max_nfev=40)
                hits[0] += int(np.any(np.isclose(res.x, lo, atol=1e-6) | np.isclose(res.x, hi, atol=1e-6) & (hi != 0.0)))
                z = res.x
                v = float(2 * res.cost)
                vals[gi] = v if vals[gi] is None else min(vals[gi], v)
            print("profile", nm, flush=True)
        prof[nm] = vals
        put(f"prof_{nm}_chi2", [None if (v is None or v != v) else float(v) for v in vals])
        v = np.array(vals, float)
        inside = grid[np.nan_to_num(v, nan=np.inf) <= 3.84]
        put(f"prof_{nm}_ci_lo_mult", float(math.exp(inside.min())) if len(inside) else None)
        put(f"prof_{nm}_ci_hi_mult", float(math.exp(inside.max())) if len(inside) else None)
        put(f"prof_{nm}_flat", bool(np.nanmax(v) < 3.84))
        put(f"prof_{nm}_max_chi2", float(np.nanmax(v)))
    put("prof_bound_hits", hits[0])
    put("prof_grid_mult_lo", float(math.exp(grid[0]))); put("prof_grid_mult_hi", float(math.exp(grid[-1])))
    put("prof_n_flat", int(sum(R[f"mathx_prof_{nm}_flat"]["value"] for nm in NAMES)))
    labels = {"b": "b", "eta": "η", "alpha": "α", "sigma": "σ", "gamma1": "γ₁", "gamma2": "γ₂", "delta": "δ", "nu": "ν"}
    fig, axs = plt.subplots(2, 4, figsize=(7.2, 3.8), sharey=True)
    for ax, nm in zip(axs.flat, NAMES):
        ax.plot(np.exp(grid), np.minimum(prof[nm], 40), color=C1, marker="o", ms=2)
        ax.axhline(3.84, color=C2, lw=0.8, ls="--"); ax.set_xscale("log"); ax.set_title(labels[nm], fontsize=9)
        from matplotlib.ticker import NullFormatter, FixedLocator, FixedFormatter
        ax.xaxis.set_major_locator(FixedLocator([0.3, 1, 3])); ax.xaxis.set_major_formatter(FixedFormatter(["0.3", "1", "3"]))
        ax.xaxis.set_minor_formatter(NullFormatter()); ax.tick_params(axis="x", labelsize=7)
    for ax in axs[:, 0]:
        ax.set_ylabel("Δχ²", fontsize=8)
    for ax in axs[1, :]:
        ax.set_xlabel("multiple of central value", fontsize=8)
    fig.tight_layout(); fig.savefig(FIG / "fig_mathx_profiles.png"); plt.close(fig)


# ---------------------------------------------------------------- 7. two-parameter maps
def sec7(pc1):
    bs = np.linspace(0.05, 6, 80); g2s = np.linspace(0.1, 1.5, 70)
    Z = np.array([[equilibrium(pc1.with_(b=b, gamma2=g))[2] * 100 for b in bs] for g in g2s])
    R0line = []
    for g in g2s:
        q = pc1.with_(gamma2=g, alpha=0.0)
        R0line.append(brentq(lambda b: R0(q.with_(b=b)) - 1, 1e-6, 1000))
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    cs = ax.contourf(bs, g2s, Z, levels=np.arange(0, 32, 2), cmap="Blues")
    ax.contour(bs, g2s, Z, levels=[7, 11], colors=[INK], linewidths=0.8, linestyles=["--", "-"])
    ax.plot(R0line, g2s, color=C2, lw=1.5, label="$\\mathcal{R}_0 = 1$ (for $\\alpha = 0$)")
    ax.plot([pc1.b], [pc1.gamma2], "o", color=C2)
    ax.set_xlim(bs[0], bs[-1]); ax.set_xlabel("Exposure rate $b$ (per year)"); ax.set_ylabel("Recovery rate $\\gamma_2$ (per year)")
    fig.colorbar(cs, ax=ax, label="Equilibrium problematic use (%)"); ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.savefig(FIG / "fig_mathx_map.png"); plt.close(fig)
    # floor over alpha and gamma2
    als = np.geomspace(0.002, 0.05, 60)
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    for g, col in [(0.25, C2), (0.55, C1), (1.0, C3)]:
        ax.plot(als, [equilibrium(pc1.with_(b=0.0, alpha=a, gamma2=g))[2] * 100 for a in als], color=col, lw=2,
                label=f"$\\gamma_2 = {g}$")
    ax.set_xscale("log"); ax.set_xlabel("Exogenous acquisition $\\alpha$ (per year)"); ax.set_ylabel("Floor $I^*_{\\min}$ (%)")
    ax.legend(frameon=False, fontsize=8); fig.savefig(FIG / "fig_mathx_floor.png"); plt.close(fig)
    put("floor_alpha_002_pct", equilibrium(pc1.with_(b=0.0, alpha=0.002))[2] * 100)
    put("floor_alpha_05_pct", equilibrium(pc1.with_(b=0.0, alpha=0.05))[2] * 100)


# ---------------------------------------------------------------- 8. structural alternatives
def alt_rhs(kind, p, b):
    """Right-hand sides of simpler structures, written on the same four-state vector for convenience."""
    def f(t, x):
        S, A, I, Rr = x
        if kind == "SIS":          # S <-> I only: problematic use acquired directly and lost directly
            lam = b * I; inc = (lam + p.alpha) * S
            return [p.nu - inc + p.gamma2 * I - p.nu * S, 0.0, inc - (p.gamma2 + p.nu) * I, 0.0]
        if kind == "SIRS":         # no At-Risk state
            lam = b * I; inc = (lam + p.alpha) * S
            return [p.nu - inc - p.nu * S, 0.0, inc + p.delta * Rr - (p.gamma2 + p.nu) * I,
                    p.gamma2 * I - (p.delta + p.nu) * Rr]
        lam = b * (I + (0.0 if kind == "eta0" else p.eta) * A)
        dl = 0.0 if kind == "norelapse" else p.delta
        inc = (lam + p.alpha) * S
        return [p.nu - inc + p.gamma1 * A - p.nu * S, inc - (p.sigma + p.gamma1 + p.nu) * A + dl * Rr,
                p.sigma * A - (p.gamma2 + p.nu) * I, p.gamma2 * I - (dl + p.nu) * Rr]
    return f


def alt_I(kind, p, b, x0, years):
    sol = solve_ivp(alt_rhs(kind, p, b), (0, years), x0, t_eval=[years], rtol=1e-9, atol=1e-11)
    return sol.y[:, -1]


def alt_eq(kind, p, b):
    x = np.array([0.9, 0.05 if kind not in ("SIS", "SIRS") else 0.0, 0.05, 0.0])
    return alt_I(kind, p, b, x, 400)


def sec8(pc1):
    p = pc1
    out = {}
    for kind in ["full", "eta0", "norelapse", "SIRS", "SIS"]:
        if kind == "full":
            f_eq = lambda b: equilibrium(p.with_(b=b))[2] - 0.07
            b0 = brentq(f_eq, 0.0, 200); x18 = equilibrium(p.with_(b=b0))
            b1 = brentq(lambda b: simulate(p.with_(b=b), x18, 4, t_eval=[4])[1][2, -1] - 0.11, b0, 200)
            I42 = simulate(p.with_(b=b1), x18, 24, t_eval=[24])[1][2, -1]
            Ieq = equilibrium(p.with_(b=b1))[2]
        else:
            b0 = brentq(lambda b: alt_eq(kind, p, b)[2] - 0.07, 1e-6, 200)
            x18 = alt_eq(kind, p, b0)
            put(f"alt_{kind}_b0", b0)
            put(f"alt_{kind}_R2018_pct", 100 * x18[3]); put(f"alt_{kind}_max2022_pct", 100 * alt_I(kind, p, 500.0, x18, 4)[2])
            if alt_I(kind, p, 500.0, x18, 4)[2] < 0.11:
                put(f"alt_{kind}_reachable", False)
                continue
            put(f"alt_{kind}_reachable", True)
            b1 = brentq(lambda b: alt_I(kind, p, b, x18, 4)[2] - 0.11, b0, 500)
            I42 = alt_I(kind, p, b1, x18, 24)[2]
            Ieq = alt_eq(kind, p, b1)[2]
        out[kind] = (b0, b1, I42, Ieq)
        put(f"alt_{kind}_b0", b0); put(f"alt_{kind}_b1", b1)
        put(f"alt_{kind}_I2042_pct", 100 * I42); put(f"alt_{kind}_Ieq_pct", 100 * Ieq)
    return out


# ---------------------------------------------------------------- 9. stability, robustness of indices, age structure
def routh_hurwitz(J3):
    """Coefficients of the characteristic polynomial l^3 + a1 l^2 + a2 l + a3 of a 3x3 matrix and the
    Routh-Hurwitz conditions a1 > 0, a3 > 0, a1 a2 > a3."""
    a1 = -np.trace(J3)
    a2 = 0.5 * (np.trace(J3) ** 2 - np.trace(J3 @ J3))
    a3 = -np.linalg.det(J3)
    return a1, a2, a3, (a1 > 0) and (a3 > 0) and (a1 * a2 > a3)


def sec9(pc1, ok):
    g = np.random.default_rng(SEED + 9)
    n = 20000; good = 0; good0 = 0; n0 = 0
    a3_err, a1_slack, margin = [], [], []
    for i in range(n):
        alpha = 0.0 if i % 2 == 0 else g.uniform(1e-4, 0.2)
        p = Params(b=g.uniform(0.2, 8), eta=g.uniform(0, 1), alpha=alpha, sigma=g.uniform(0.05, 1.5),
                   gamma1=g.uniform(0.05, 1.5), gamma2=g.uniform(0.1, 2), delta=g.uniform(0, 0.6), nu=g.uniform(0.03, 0.1))
        if alpha == 0.0 and R0(p) <= 1:
            continue
        x = equilibrium(p)
        a1_, a2_, a3_, ok_ = routh_hurwitz(reduced_jac(p, x))
        S_, A_, I_, R_ = x
        Lam = p.b * (I_ + p.eta * A_) + p.alpha
        a, c, d = p.sigma + p.gamma1 + p.nu, p.gamma2 + p.nu, p.delta + p.nu
        D = a * c * d - p.sigma * p.gamma2 * p.delta
        m = 1 + p.sigma / c + p.sigma * p.gamma2 / (c * d)
        a3_formula = D * p.alpha / Lam + Lam * c * d * m
        a1_bound = c + d + Lam
        a3_err.append(abs(a3_ - a3_formula) / a3_)
        a1_slack.append(a1_ - a1_bound)
        margin.append((a1_ * a2_ - a3_) / (a1_ * a2_))
        if alpha == 0.0:
            n0 += 1; good0 += ok_
        else:
            good += ok_
    put("rh_a3_formula_max_rel_err", float(max(a3_err)))
    put("rh_a1_bound_min_slack", float(min(a1_slack)))
    put("rh_hurwitz_margin_min", float(min(margin))); put("rh_hurwitz_margin_min_c", float(1 - min(margin))); put("rh_hurwitz_margin_median", float(np.median(margin)))
    put("rh_n_alpha_pos", n // 2); put("rh_ok_alpha_pos", int(good))
    put("rh_n_alpha0_endemic", n0); put("rh_ok_alpha0_endemic", int(good0))
    a1, a2, a3, _ = routh_hurwitz(reduced_jac(pc1, equilibrium(pc1)))
    put("rh_central_a1", a1); put("rh_central_a2", a2); put("rh_central_a3", a3)
    # sensitivity indices of I* across the calibrated draws (post-2018, Scenario P equilibrium)
    names = ["b", "eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"]
    rows = []
    for r in ok:
        p = Params(b=r["b1"], **{k: r[k] for k in cal.PARAM_ORDER})
        rows.append([elasticity(lambda q: equilibrium(q)[2], p, k) for k in names])
    rows = np.array(rows)
    ranks = np.argsort(np.argsort(-np.abs(rows), axis=1), axis=1)
    for j, k in enumerate(names):
        m, lo, hi = np.median(rows[:, j]), np.percentile(rows[:, j], 5), np.percentile(rows[:, j], 95)
        put(f"idxd_{k}_median", float(m)); put(f"idxd_{k}_p05", float(lo)); put(f"idxd_{k}_p95", float(hi))
        put(f"idxd_{k}_top_pct", float(100 * np.mean(ranks[:, j] == 0)))
    put("idxd_n", len(ok))
    put("idxd_b_below_gamma2_pct", float(100 * np.mean(np.abs(rows[:, 0]) < np.abs(rows[:, 5]))))
    # exponential turnover against ageing through bands: the structured model with no heterogeneity
    from dpm import agesex as AS
    st0 = AS.Struct(phi=1.0, psi=1.0, eps_a=0.0, eps_s=0.0)
    q0 = pc1.with_(alpha=0.0)
    put("ageing_R0_exponential", R0(q0)); put("ageing_R0_bands", AS.R0_struct(q0, st0))
    x = AS.equilibrium(pc1, st0)
    put("ageing_I_exponential_pct", 100 * equilibrium(pc1)[2]); put("ageing_I_bands_pct", 100 * AS.whole_prev(x))
    put("ageing_I_bands_hbsc_pct", 100 * AS.hbsc_obs(x)["total"])
    gp = AS.group_prev(x)
    for k in range(4):
        put(f"ageing_I_band{k + 1}_pct", 100 * float(gp[k]))


def main():
    t0 = time.time()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)
    ok = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]]
    save = lambda: json.dump(R, open(OUT / "results_mathx.json", "w"), indent=1, ensure_ascii=False)
    sec1(pc1); print("1", time.time() - t0)
    sec2(pc1); print("2", time.time() - t0)
    sec3(pc0, pc1, ok); print("3", time.time() - t0)
    sec4(pc0, pc1); print("4", time.time() - t0)
    sec5(pc1); print("5", time.time() - t0, flush=True); save()
    sec7(pc1); print("7", time.time() - t0, flush=True); save()
    sec8(pc1); print("8", time.time() - t0, flush=True); save()
    sec9(pc1, ok); print("9", time.time() - t0, flush=True); save()
    sec6(pc1, x18); print("6", time.time() - t0, flush=True)
    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_mathx.json", "w"), indent=1, ensure_ascii=False)
    print(f"done, {len(R)} results, {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--only9":
        R.update(json.load(open(OUT / "results_mathx.json")))
        pc0, pc1, da, x18 = ra.central(json.load(open(HERE / "data" / "parameter_ranges.json")))
        ok = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]]
        sec9(pc1, ok); json.dump(R, open(OUT / "results_mathx.json", "w"), indent=1, ensure_ascii=False); print("done 9")
    elif len(sys.argv) > 1 and sys.argv[1] == "--only8":
        R.update(json.load(open(OUT / "results_mathx.json")))
        pc0, pc1, da, x18 = ra.central(json.load(open(HERE / "data" / "parameter_ranges.json")))
        sec8(pc1); json.dump(R, open(OUT / "results_mathx.json", "w"), indent=1, ensure_ascii=False)
    elif len(sys.argv) > 1 and sys.argv[1] == "--from8":
        # resume: keep sections 1 to 7 already saved, run 8 and 6
        R.update(json.load(open(OUT / "results_mathx.json")))
        t0 = time.time()
        ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
        pc0, pc1, da, x18 = ra.central(ranges)
        sec8(pc1); json.dump(R, open(OUT / "results_mathx.json", "w"), indent=1, ensure_ascii=False); print("8", flush=True)
        sec6(pc1, x18); put("runtime_minutes_resume", (time.time() - t0) / 60)
        json.dump(R, open(OUT / "results_mathx.json", "w"), indent=1, ensure_ascii=False); print("done", flush=True)
    else:
        main()
