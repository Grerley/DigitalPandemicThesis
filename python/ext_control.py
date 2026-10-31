"""Optimal allocation of intervention effort over time (thesis Chapter 9, section "Optimal allocation of effort over time").

Writes outputs/results_ctl.json (keys prefixed ctl_) and figures/fig_ctl_*.png.

Problem. On [0, T], T = 20 years from the intervention start year of Chapter 9, minimise

    J(u) = int_0^T [ I(t) + (w1/2) u1(t)^2 + (w2/2) u2(t)^2 ] dt

subject to the S-A-I-R model (4.2) with b -> b (1 - U1MAX u1) and gamma2 -> gamma2 (1 + U2MAX u2),
u1, u2 in [0, 1]: each control is the share of its assumed efficacy ceiling in use (U1MAX: largest cut
in exposure, U2MAX: largest proportional rise in recovery, the upper ends of Table 9.1).
Costs enter only through the relative weights w1, w2 (cost of one year at the full ceiling, in units of
prevalence-years); no cost figures are used.

Method. Pontryagin's maximum principle; forward-backward sweep (Lenhart and Workman 2007, ch. 4) with
fourth-order Runge-Kutta on a fixed grid, convex-combination control updates and the relative
convergence test of Lenhart and Workman. Checks: step-halving, independent direct optimisation over
piecewise-constant controls, different initial guesses, and the uncontrolled limit.

Sections
  1. Central parameters, Scenario P: reference solution, numerical checks, comparison with no control
     and with the best constant controls of the same total cost.
  2. Grid of weights: ratio w1/w2 over 1/256 to 16 at three overall levels.
  3. Robustness: calibrated draws (Scenario P) and Scenario T at the central point.
  0. Weight-free marginal values of the two controls along the uncontrolled path (gradient of the
     prevalence-years at zero control), checked against finite differences.

    python ext_control.py
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np
from scipy.optimize import minimize, minimize_scalar

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt
from dpm.model import Params, simulate
from dpm import calibrate as cal

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R = {}
SEED = ra.SEED + 3100

T_HOR = 20.0          # years of control
N_STEPS = 2000        # RK4 steps on [0, T] (h = 0.01 year); halved in the step-size check
U1MAX = 0.40          # upper end of the assumed platform-design range (Table 9.1)
U2MAX = 1.00          # upper end of the assumed rise in recovery (Table 9.1)
W_REF = 0.10          # reference weight level: one year at a full ceiling counts as W/2 prevalence-years
RATIOS = [2.0 ** k for k in range(-8, 5)]     # w1/w2 from 1/256 to 16
LEVELS = {"low": W_REF / 4, "ref": W_REF, "high": W_REF * 4}
TOL = 1e-4            # relative convergence tolerance of the sweep


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["ctl_" + k] = {"value": v, "note": note}


# ---------------------------------------------------------------- state and adjoint systems (scalar code for speed)
def f_state(x, u1, u2, p):
    S, A, I, Rr = x
    bt = p.b * (1.0 - U1MAX * u1); g2 = p.gamma2 * (1.0 + U2MAX * u2)
    inc = (bt * (I + p.eta * A) + p.alpha) * S
    return (p.nu - inc + p.gamma1 * A - p.nu * S,
            inc - (p.sigma + p.gamma1 + p.nu) * A + p.delta * Rr,
            p.sigma * A - (g2 + p.nu) * I,
            g2 * I - (p.delta + p.nu) * Rr)


def f_adj(lam, x, u1, u2, p):
    """d lambda / dt = - dH/dx  (equations (9.5))."""
    lS, lA, lI, lR = lam
    S, A, I, Rr = x
    bt = p.b * (1.0 - U1MAX * u1); g2 = p.gamma2 * (1.0 + U2MAX * u2)
    Lam = bt * (I + p.eta * A) + p.alpha
    dSA = lS - lA
    return (Lam * dSA + p.nu * lS,
            bt * p.eta * S * dSA - p.gamma1 * dSA + p.sigma * (lA - lI) + p.nu * lA,
            -1.0 + bt * S * dSA + g2 * (lI - lR) + p.nu * lI,
            p.delta * (lR - lA) + p.nu * lR)


def _add(a, b, s):
    return (a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2], a[3] + s * b[3])


def forward(x0, u1, u2, p, n):
    """RK4 for the state on a uniform grid; controls at midpoints are averages of neighbouring nodes."""
    h = T_HOR / n
    X = np.empty((n + 1, 4)); x = tuple(float(v) for v in x0); X[0] = x
    for i in range(n):
        a1, a2 = u1[i], u2[i]; m1, m2 = 0.5 * (u1[i] + u1[i + 1]), 0.5 * (u2[i] + u2[i + 1]); e1, e2 = u1[i + 1], u2[i + 1]
        k1 = f_state(x, a1, a2, p)
        k2 = f_state(_add(x, k1, h / 2), m1, m2, p)
        k3 = f_state(_add(x, k2, h / 2), m1, m2, p)
        k4 = f_state(_add(x, k3, h), e1, e2, p)
        x = tuple(x[j] + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(4))
        X[i + 1] = x
    return X


def backward(X, u1, u2, p, n):
    """RK4 backwards in time for the adjoint from lambda(T) = 0."""
    h = T_HOR / n
    L = np.zeros((n + 1, 4)); lam = (0.0, 0.0, 0.0, 0.0)
    for i in range(n, 0, -1):
        xe, xm, xs = tuple(X[i]), tuple(0.5 * (X[i] + X[i - 1])), tuple(X[i - 1])
        e1, e2 = u1[i], u2[i]; m1, m2 = 0.5 * (u1[i] + u1[i - 1]), 0.5 * (u2[i] + u2[i - 1]); s1, s2 = u1[i - 1], u2[i - 1]
        k1 = f_adj(lam, xe, e1, e2, p)
        k2 = f_adj(_add(lam, k1, -h / 2), xm, m1, m2, p)
        k3 = f_adj(_add(lam, k2, -h / 2), xm, m1, m2, p)
        k4 = f_adj(_add(lam, k3, -h), xs, s1, s2, p)
        lam = tuple(lam[j] - h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]) for j in range(4))
        L[i - 1] = lam
    return L


def sweep_horizon(T, *args, **kw):
    """Run the sweep on [0, T] instead of [0, T_HOR] (used only for the horizon check)."""
    global T_HOR
    keep = T_HOR
    try:
        T_HOR = T
        return sweep(*args, **kw)
    finally:
        T_HOR = keep


def characterise(X, L, p, w1, w2):
    """Projection formulas (9.7): minimisers of the Hamiltonian over the control box."""
    S, A, I = X[:, 0], X[:, 1], X[:, 2]
    phi1 = U1MAX * p.b * (I + p.eta * A) * S * (L[:, 1] - L[:, 0])     # marginal value of exposure effort
    phi2 = U2MAX * p.gamma2 * I * (L[:, 2] - L[:, 3])                   # marginal value of recovery effort
    return np.clip(phi1 / w1, 0, 1.0), np.clip(phi2 / w2, 0, 1.0), phi1, phi2


def trap(y, n):
    h = T_HOR / n
    return float(h * (y.sum() - 0.5 * (y[0] + y[-1])))


def sweep(x0, p, w1, w2, n=N_STEPS, u_init=(0.0, 0.0), tol=TOL, max_iter=2000):
    """Forward-backward sweep. Returns dict with grid, state, adjoint, controls, iterations and J components."""
    u1 = np.full(n + 1, float(u_init[0])); u2 = np.full(n + 1, float(u_init[1]))
    X = forward(x0, u1, u2, p, n); L = backward(X, u1, u2, p, n)
    for it in range(1, max_iter + 1):
        o1, o2, oX, oL = u1.copy(), u2.copy(), X.copy(), L.copy()
        X = forward(x0, u1, u2, p, n)
        L = backward(X, u1, u2, p, n)
        n1, n2, _, _ = characterise(X, L, p, w1, w2)
        u1 = 0.5 * (n1 + o1); u2 = 0.5 * (n2 + o2)
        rel = [np.abs(a - b).sum() / max(np.abs(a).sum(), 1e-12) for a, b in
               [(u1, o1), (u2, o2), (X, oX), (L, oL)] if np.abs(a).sum() > 1e-12]
        if it > 2 and max(rel + [0.0]) < tol:
            break
    else:
        raise RuntimeError("sweep did not converge")
    X = forward(x0, u1, u2, p, n)
    t = np.linspace(0, T_HOR, n + 1)
    cost = 0.5 * w1 * u1 ** 2 + 0.5 * w2 * u2 ** 2
    return dict(t=t, X=X, L=backward(X, u1, u2, p, n), u1=u1, u2=u2, iters=it, n=n,
                PY=trap(X[:, 2], n), C=trap(cost, n), C1=trap(0.5 * w1 * u1 ** 2, n), C2=trap(0.5 * w2 * u2 ** 2, n))


def prev_years(x0, p, c1, c2, n=N_STEPS):
    """Integral of I over [0, T] under constant controls (c1, c2)."""
    X = forward(x0, np.full(n + 1, c1), np.full(n + 1, c2), p, n)
    return trap(X[:, 2], n), X


def best_constant(x0, p, w1, w2, C, n=N_STEPS, only=None):
    """Constant controls with total cost T (w1 c1^2 + w2 c2^2)/2 = C that minimise the integral of I.

    Parametrise c1 = r1 cos(theta), c2 = r2 sin(theta), r_i = sqrt(2C / (T w_i)), restricted to the box.
    only = 'u1' or 'u2' gives the best single-lever constant policy of the same cost (if feasible)."""
    r1 = math.sqrt(2 * C / (T_HOR * w1)); r2 = math.sqrt(2 * C / (T_HOR * w2))
    if only == "u1":
        return (None if r1 > 1 else (prev_years(x0, p, r1, 0.0, n)[0], r1, 0.0))
    if only == "u2":
        return (None if r2 > 1 else (prev_years(x0, p, 0.0, r2, n)[0], 0.0, r2))
    lo = 0.0 if r1 <= 1 else math.acos(1 / r1)
    hi = math.pi / 2 if r2 <= 1 else math.asin(1 / r2)
    if lo > hi:
        return None

    def f(th):
        return prev_years(x0, p, r1 * math.cos(th), r2 * math.sin(th), n)[0]
    grid = np.linspace(lo, hi, 13)
    vals = [f(th) for th in grid]
    k = int(np.argmin(vals))
    a, b = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
    res = minimize_scalar(f, bounds=(a, b), method="bounded", options={"xatol": 1e-6})
    th = res.x if res.fun <= vals[k] else grid[k]
    return f(th), r1 * math.cos(th), r2 * math.sin(th)


def compare(x0, p, w1, w2, sol):
    """Prevalence-years under no control, the optimal control and the best constant control of equal cost."""
    PY0, X0 = prev_years(x0, p, 0.0, 0.0, sol["n"])
    bc = best_constant(x0, p, w1, w2, sol["C"], sol["n"])
    PYc, c1, c2 = bc
    av_opt = PY0 - sol["PY"]; av_c = PY0 - PYc
    return dict(PY0=PY0, X0=X0, PYc=PYc, c1=c1, c2=c2, av_opt=av_opt, av_c=av_c,
                gain_pct=100 * (av_opt - av_c) / av_opt if av_opt > 0 else float("nan"),
                const_share_pct=100 * av_c / av_opt if av_opt > 0 else float("nan"))


def summarise(sol, w1, w2):
    u1, u2, n = sol["u1"], sol["u2"], sol["n"]
    return dict(u1_mean=trap(u1, n) / T_HOR, u2_mean=trap(u2, n) / T_HOR,
                u1_share_cost_pct=100 * sol["C1"] / sol["C"] if sol["C"] > 0 else float("nan"),
                u1_sat_pct=100 * float(np.mean(u1 >= 1 - 1e-3)), u2_sat_pct=100 * float(np.mean(u2 >= 1 - 1e-3)))


def state_at_start(p_post, x18, start, alpha_t=None):
    return simulate(p_post, x18, start - 2018.0, t_eval=[start - 2018.0], alpha_t=alpha_t)[1][:, -1]


# ---------------------------------------------------------------- 0. weight-free marginal values
def marginal(x0, p, n=N_STEPS):
    """Switching functions along the uncontrolled path. Their integrals are the reductions in prevalence-years
    per unit of a constant control held over [0, T], at zero control (first-order effect)."""
    z = np.zeros(n + 1)
    X = forward(x0, z, z, p, n); L = backward(X, z, z, p, n)
    _, _, p1, p2 = characterise(X, L, p, 1.0, 1.0)
    return p1, p2, trap(p1, n), trap(p2, n), X


def sec0(pc1, x0):
    p1, p2, P1, P2, X = marginal(x0, pc1)
    PY0 = trap(X[:, 2], N_STEPS); eps = 1e-4
    fd1 = (PY0 - prev_years(x0, pc1, eps, 0.0)[0]) / eps
    fd2 = (PY0 - prev_years(x0, pc1, 0.0, eps)[0]) / eps
    put("mv_Phi1", P1); put("mv_Phi2", P2); put("mv_ratio", P2 / P1); put("mv_ratio_inv", P1 / P2)
    put("mv_Phi1_per1000", 1000 * P1); put("mv_Phi2_per1000", 1000 * P2)
    put("mv_fd_rel_err", max(abs(fd1 - P1) / P1, abs(fd2 - P2) / P2), "adjoint gradient against finite differences")
    # per unit proportional change in the rate (1% cut in b against 1% rise in gamma2)
    put("mv_ratio_per_unit_rate", (P2 / U2MAX) / (P1 / U1MAX))
    Rm = json.load(open(OUT / "results.json"))
    put("mv_elast_ratio", abs(Rm["elast_I_gamma2"]["value"]) / Rm["elast_I_b"]["value"],
        "|index of I* to gamma2| / index of I* to b, Chapter 5")
    put("mv_Phi1_pct_PY0", 100 * P1 / PY0); put("mv_Phi2_pct_PY0", 100 * P2 / PY0)
    # shape over time
    n = N_STEPS
    for yr in (0, 10, 19):
        j = int(round(yr / T_HOR * n))
        put(f"mv_phi1_at{yr}", p1[j]); put(f"mv_phi2_at{yr}", p2[j])
    return p1, p2


# ---------------------------------------------------------------- 1. reference solution and checks
def sec1(pc1, x0, start):
    w1 = w2 = W_REF
    t0 = time.time()
    sol = sweep(x0, pc1, w1, w2)
    put("ref_runtime_s", time.time() - t0)
    put("ref_iters", sol["iters"])
    cmp_ = compare(x0, pc1, w1, w2, sol)
    s = summarise(sol, w1, w2)
    n = sol["n"]
    put("T_years", T_HOR); put("start_year", start); put("end_year", start + T_HOR)
    put("u1max_pct", 100 * U1MAX); put("u2max_pct", 100 * U2MAX); put("w_ref", W_REF,
        "weight level: running one control at value 1 for a year counts as w/2 prevalence-years")
    put("n_steps", N_STEPS); put("h_years", T_HOR / N_STEPS); put("tol", TOL)
    put("ref_I0_pct", 100 * x0[2])
    put("ref_PY0_per1000", 1000 * cmp_["PY0"]); put("ref_PYopt_per1000", 1000 * sol["PY"])
    put("ref_PYconst_per1000", 1000 * cmp_["PYc"])
    put("ref_averted_opt_per1000", 1000 * cmp_["av_opt"]); put("ref_averted_const_per1000", 1000 * cmp_["av_c"])
    put("ref_averted_opt_pct", 100 * cmp_["av_opt"] / cmp_["PY0"])
    put("ref_averted_const_pct", 100 * cmp_["av_c"] / cmp_["PY0"])
    put("ref_gain_pct", cmp_["gain_pct"], "extra prevalence-years averted by the optimal over the best equal-cost constant, % of optimal")
    put("ref_const_share_pct", cmp_["const_share_pct"])
    put("ref_const_c1_pct", 100 * cmp_["c1"]); put("ref_const_c2_pct", 100 * cmp_["c2"])
    for k, v in s.items():
        put("ref_" + k, v)
    put("ref_u1_mean_pct", 100 * s["u1_mean"]); put("ref_u2_mean_pct", 100 * s["u2_mean"])
    put("w_ref_half_per1000", 1000 * W_REF / 2, "one year at a full ceiling, in person-years per 1,000")
    put("w_low_half_per1000", 1000 * LEVELS["low"] / 2); put("w_high_half_per1000", 1000 * LEVELS["high"] / 2)
    # profile features
    u1, u2, t = sol["u1"], sol["u2"], sol["t"]
    put("ref_u1_at0_pct", 100 * u1[0]); put("ref_u2_at0_pct", 100 * u2[0])
    put("ref_u1_max_pct", 100 * u1.max()); put("ref_u2_max_pct", 100 * u2.max())
    j5 = int(round(5 / T_HOR * n)); j10 = int(round(10 / T_HOR * n)); j15 = int(round(15 / T_HOR * n))
    put("ref_u1_at5_pct", 100 * u1[j5]); put("ref_u2_at5_pct", 100 * u2[j5])
    put("ref_u1_at10_pct", 100 * u1[j10]); put("ref_u2_at10_pct", 100 * u2[j10])
    put("ref_u1_at15_pct", 100 * u1[j15]); put("ref_u2_at15_pct", 100 * u2[j15])
    # time at which u2 leaves its ceiling (if it starts there)
    put("ref_u2_settle_years", float(t[int(np.argmax(u2 <= 1.02 * u2[j10]))]),
        "first time u2 is within 2% of its mid-horizon value")
    put("ref_u2_ratio_start_to_mid", u2[0] / u2[j10]); put("ref_u1_ratio_end_to_mid", u1[j15] / u1[j10])
    # length of the terminal taper: time before T at which each control falls below half its mid-horizon value
    for nm, u in (("u1", u1), ("u2", u2)):
        k = n - int(np.argmax((u[::-1] >= 0.5 * u[j10])))
        put(f"ref_{nm}_taper_years", float(T_HOR - t[k]))
    j18 = int(round((2042 - start) / T_HOR * n))   # the year 2042 used throughout Chapter 9
    put("ref_I_end_none_pct", 100 * cmp_["X0"][-1, 2]); put("ref_I_end_opt_pct", 100 * sol["X"][-1, 2])
    put("ref_I_2042_none_pct", 100 * cmp_["X0"][j18, 2]); put("ref_I_2042_opt_pct", 100 * sol["X"][j18, 2])
    _, Xc = prev_years(x0, pc1, cmp_["c1"], cmp_["c2"])
    put("ref_I_2042_const_pct", 100 * Xc[j18, 2]); put("ref_I_end_const_pct", 100 * Xc[-1, 2])
    put("ref_I_min_opt_pct", 100 * sol["X"][:, 2].min())
    # gain measured only to 2042 (year 18), away from the terminal taper
    py18 = lambda X: (T_HOR / n) * (X[:j18 + 1, 2].sum() - 0.5 * (X[0, 2] + X[j18, 2]))
    a_opt = py18(cmp_["X0"]) - py18(sol["X"]); a_c = py18(cmp_["X0"]) - py18(Xc)
    put("ref_gain_to2042_pct", 100 * (a_opt - a_c) / a_opt)
    # longer horizon: the same problem on [0, 40]; controls over the first 20 years
    s40 = sweep_horizon(40.0, x0, pc1, w1, w2, n=2 * N_STEPS)
    k20 = N_STEPS; k10 = N_STEPS // 2
    put("T40_u1_at10_pct", 100 * s40["u1"][k10]); put("T40_u2_at10_pct", 100 * s40["u2"][k10])
    put("T40_u1_at20_pct", 100 * s40["u1"][k20]); put("T40_u2_at20_pct", 100 * s40["u2"][k20])
    put("T40_max_du_first10", float(max(np.abs(s40["u1"][:k10 + 1] - sol["u1"][:j10 + 1]).max(),
                                        np.abs(s40["u2"][:k10 + 1] - sol["u2"][:j10 + 1]).max())))
    # single-lever constant policies of the same cost
    for lever in ("u1", "u2"):
        r = best_constant(x0, pc1, w1, w2, sol["C"], only=lever)
        if r is None:
            put(f"ref_only_{lever}_feasible", 0)
        else:
            put(f"ref_only_{lever}_feasible", 1)
            put(f"ref_only_{lever}_level_pct", 100 * (r[1] if lever == "u1" else r[2]))
            put(f"ref_only_{lever}_averted_pct_of_opt", 100 * (cmp_["PY0"] - r[0]) / cmp_["av_opt"])
            put(f"ref_only_{lever}_PY_per1000", 1000 * r[0])
            _, Xl = prev_years(x0, pc1, r[1], r[2])
            put(f"ref_only_{lever}_I_2042_pct", 100 * Xl[int(round((2042 - start) / T_HOR * n)), 2])
    # adjoint signs: the marginal values that drive the controls are positive throughout
    L, X = sol["L"], sol["X"]
    put("ref_lam_ratio_I_over_A_at0", float(L[0, 2] / L[0, 1]))
    put("ref_min_lA_minus_lS", float((L[:-1, 1] - L[:-1, 0]).min()))
    put("ref_min_lI_minus_lR", float((L[:-1, 2] - L[:-1, 3]).min()))
    put("ref_lI_at0", float(L[0, 2])); put("ref_lA_at0", float(L[0, 1]))
    # ---------------- numerical checks
    # (a) step halving
    sol2 = sweep(x0, pc1, w1, w2, n=2 * N_STEPS)
    u1h = sol2["u1"][::2]; u2h = sol2["u2"][::2]
    put("chk_halving_max_du1", float(np.abs(u1h - sol["u1"]).max()))
    put("chk_halving_max_du2", float(np.abs(u2h - sol["u2"]).max()))
    put("chk_halving_rel_dPY", abs(sol2["PY"] - sol["PY"]) / sol["PY"])
    put("chk_halving_rel_dJ", abs((sol2["PY"] + sol2["C"]) - (sol["PY"] + sol["C"])) / (sol["PY"] + sol["C"]))
    # (b) other initial guesses
    worst = 0.0
    for ui in [(1.0, 1.0), (1.0, 0.0), (0.0, 1.0), (0.5, 0.5)]:
        s3 = sweep(x0, pc1, w1, w2, u_init=ui)
        worst = max(worst, float(np.abs(s3["u1"] - sol["u1"]).max()), float(np.abs(s3["u2"] - sol["u2"]).max()))
    put("chk_init_max_du", worst, "largest difference in controls across four other initial guesses")
    # (c) direct optimisation over piecewise-constant controls on two-year pieces (independent of the adjoint)
    Jopt = sol["PY"] + sol["C"]
    m = 10; n_d = 1000
    edges = np.linspace(0, T_HOR, m + 1); tg = np.linspace(0, T_HOR, n_d + 1)
    idx = np.minimum(np.searchsorted(edges, tg, side="right") - 1, m - 1)

    def Jpc(z):
        a1 = z[:m][idx]; a2 = z[m:][idx]
        X = forward(x0, a1, a2, pc1, n_d)
        return trap(X[:, 2], n_d) + T_HOR / m * (0.5 * w1 * (z[:m] ** 2).sum() + 0.5 * w2 * (z[m:] ** 2).sum())
    # start from the piece averages of the PMP solution, and separately from zero
    z0 = np.concatenate([[sol["u1"][(sol["t"] >= edges[k]) & (sol["t"] < edges[k + 1])].mean() for k in range(m)],
                         [sol["u2"][(sol["t"] >= edges[k]) & (sol["t"] < edges[k + 1])].mean() for k in range(m)]])
    bnds = [(0, 1)] * (2 * m)
    best = None
    for zs in (z0, np.zeros(2 * m)):
        r = minimize(Jpc, zs, method="L-BFGS-B", bounds=bnds, options={"maxiter": 500, "eps": 1e-7, "ftol": 1e-13})
        if best is None or r.fun < best.fun:
            best = r
    put("chk_direct_J", float(best.fun)); put("chk_pmp_J", float(Jopt))
    put("chk_direct_minus_pmp_rel", (best.fun - Jopt) / Jopt,
        "direct piecewise-constant optimum minus PMP optimum, relative; non-negative expected")
    put("chk_direct_max_dpiece", float(np.abs(best.x - z0).max()),
        "largest difference between direct pieces and two-year averages of the PMP controls")
    # (d) uncontrolled limit: huge weights reproduce the uncontrolled trajectory
    s4 = sweep(x0, pc1, 1e6, 1e6)
    PY0 = cmp_["PY0"]
    put("chk_hugeweights_rel_dPY", abs(s4["PY"] - PY0) / PY0)
    # (e) uncontrolled grid trajectory against the adaptive solver of the main pipeline
    tt = np.linspace(0, T_HOR, 21)
    _, Y = simulate(pc1, x0, T_HOR, t_eval=tt)
    put("chk_rk4_vs_adaptive_maxabs", float(np.abs(cmp_["X0"][:: N_STEPS // 20, :].T - Y).max()))
    return sol, cmp_


# ---------------------------------------------------------------- 2. weight grid
def sec2(pc1, x0):
    rows = {}
    for lev, W in LEVELS.items():
        for rho in RATIOS:
            w2 = W / math.sqrt(rho); w1 = W * math.sqrt(rho)   # geometric mean of the weights fixed at W
            sol = sweep(x0, pc1, w1, w2)
            alt = sweep(x0, pc1, w1, w2, u_init=(1.0, 1.0))     # second start: both controls at their ceilings
            MULTI.append(max(float(np.abs(alt["u1"] - sol["u1"]).max()), float(np.abs(alt["u2"] - sol["u2"]).max())))
            if alt["PY"] + alt["C"] < sol["PY"] + sol["C"]:
                sol = alt
            c = compare(x0, pc1, w1, w2, sol)
            s = summarise(sol, w1, w2)
            rows[(lev, rho)] = dict(sol=sol, cmp=c, s=s, w1=w1, w2=w2)
            tag = f"{lev}_r{_rtag(rho)}"
            put(f"grid_{tag}_u1_mean_pct", 100 * s["u1_mean"]); put(f"grid_{tag}_u2_mean_pct", 100 * s["u2_mean"])
            put(f"grid_{tag}_u1_cost_share_pct", s["u1_share_cost_pct"])
            put(f"grid_{tag}_averted_pct", 100 * c["av_opt"] / c["PY0"])
            put(f"grid_{tag}_gain_pct", c["gain_pct"])
            put(f"grid_{tag}_u2_sat_pct", s["u2_sat_pct"]); put(f"grid_{tag}_u1_sat_pct", s["u1_sat_pct"])
            print(f"  {lev} rho={rho:g}: u1 {s['u1_mean']:.3f} u2 {s['u2_mean']:.3f} share1 {s['u1_share_cost_pct']:.1f} "
                  f"avert {100 * c['av_opt'] / c['PY0']:.1f}% gain {c['gain_pct']:.2f}% it {sol['iters']}", flush=True)
    for lev in LEVELS:
        g = [rows[(lev, r)]["cmp"]["gain_pct"] for r in RATIOS]
        put(f"grid_{lev}_gain_max_pct", float(max(g))); put(f"grid_{lev}_gain_min_pct", float(min(g)))
        # ratio at which both controls are used at the same share of their ceilings (log-linear interpolation)
        dv = np.array([rows[(lev, r)]["s"]["u1_mean"] - rows[(lev, r)]["s"]["u2_mean"] for r in RATIOS])
        lr = np.log(RATIOS)
        cross = [k for k in range(len(RATIOS) - 1) if dv[k] * dv[k + 1] <= 0]
        k = cross[0]
        x = lr[k] - dv[k] * (lr[k + 1] - lr[k]) / (dv[k + 1] - dv[k])
        put(f"grid_{lev}_equal_intensity_ratio", float(math.exp(x)))
        put(f"grid_{lev}_equal_intensity_inv", float(math.exp(-x)))
        sh = np.array([rows[(lev, r)]["s"]["u1_share_cost_pct"] for r in RATIOS])
        put(f"grid_{lev}_u1_cost_share_max_pct", float(sh.max()))
    allg = [rows[k]["cmp"]["gain_pct"] for k in rows]
    put("grid_gain_max_pct", float(max(allg))); put("grid_gain_min_pct", float(min(allg)))
    put("grid_multistart_max_du", float(max(MULTI)), "largest control difference between two starts, all grid points")
    put("grid_n", len(rows)); put("grid_ratio_min_inv", 1 / min(RATIOS)); put("grid_ratio_max", max(RATIOS))
    return rows


def _rtag(rho):
    k = int(round(math.log2(rho)))
    return f"m{-k}" if k < 0 else f"{k}"       # ratio 2^k: r_m8 is 1/256, r_0 is 1, r_4 is 16


# ---------------------------------------------------------------- 3. robustness
def sec3(ok, pcT, x0T, n_draws=20, n_mv=200):
    g = np.random.default_rng(SEED)
    # weight-free marginal values at a random sample of calibrated draws
    pick_mv = sorted(g.choice(len(ok), size=n_mv, replace=False).tolist())
    ratios = []
    for j in pick_mv:
        r = ok[j]
        p1 = Params(b=r["b1"], **{k: r[k] for k in cal.PARAM_ORDER})
        x0 = state_at_start(p1, np.array(r["x18"]), START)
        _, _, P1, P2, _ = marginal(x0, p1, n=1000)
        ratios.append(P2 / P1)
    ratios = np.array(ratios); MV_RATIOS[:] = ratios.tolist()
    put("rob_mv_n", n_mv)
    put("rob_mv_ratio_median", float(np.median(ratios)))
    put("rob_mv_ratio_p05", float(np.percentile(ratios, 5))); put("rob_mv_ratio_p95", float(np.percentile(ratios, 95)))
    put("rob_mv_ratio_min", float(ratios.min())); put("rob_mv_ratio_above1_pct", 100 * float(np.mean(ratios > 1)))
    pick = sorted(g.choice(len(ok), size=n_draws, replace=False).tolist())
    res = []
    for j in pick:
        r = ok[j]
        p1 = Params(b=r["b1"], **{k: r[k] for k in cal.PARAM_ORDER})
        x0 = state_at_start(p1, np.array(r["x18"]), START)
        sol = sweep(x0, p1, W_REF, W_REF)
        c = compare(x0, p1, W_REF, W_REF, sol)
        s = summarise(sol, W_REF, W_REF)
        res.append(dict(gain=c["gain_pct"], share1=s["u1_share_cost_pct"], u1=s["u1_mean"], u2=s["u2_mean"],
                        av=100 * c["av_opt"] / c["PY0"]))
    for k in ["gain", "share1", "u1", "u2", "av"]:
        v = np.array([d[k] for d in res])
        sc = 100 if k in ("u1", "u2") else 1
        put(f"rob_{k}_median", sc * float(np.median(v))); put(f"rob_{k}_min", sc * float(v.min()))
        put(f"rob_{k}_max", sc * float(v.max()))
    put("rob_n_draws", n_draws)
    put("rob_share1_below50_n", int(sum(d["share1"] < 50 for d in res)))
    # Scenario T, central parameters
    _, _, P1, P2, _ = marginal(x0T, pcT)
    put("T_mv_ratio", P2 / P1)
    sol = sweep(x0T, pcT, W_REF, W_REF)
    c = compare(x0T, pcT, W_REF, W_REF, sol)
    s = summarise(sol, W_REF, W_REF)
    put("T_I0_pct", 100 * x0T[2])
    put("T_gain_pct", c["gain_pct"]); put("T_u1_cost_share_pct", s["u1_share_cost_pct"])
    put("T_u1_mean_pct", 100 * s["u1_mean"]); put("T_u2_mean_pct", 100 * s["u2_mean"])
    put("T_averted_pct", 100 * c["av_opt"] / c["PY0"]); put("T_averted_per1000", 1000 * c["av_opt"])
    return res, sol


# ---------------------------------------------------------------- figures
def fig_profiles(sol, cmp_, start):
    t = sol["t"] + start
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.5))
    a = ax[0]
    a.plot(t, 100 * sol["u1"], color=C2, lw=1.8, label=f"$u_1$: exposure (ceiling: $b$ cut {100 * U1MAX:.0f}%)")
    a.plot(t, 100 * sol["u2"], color=C1, lw=1.8, label=f"$u_2$: recovery (ceiling: $\\gamma_2$ up {100 * U2MAX:.0f}%)")
    a.axhline(100 * cmp_["c1"], color=C2, lw=1, ls="--")
    a.axhline(100 * cmp_["c2"], color=C1, lw=1, ls="--")
    a.set_ylabel("control (% of its efficacy ceiling)"); a.set_xlabel("year")
    a.set_ylim(0, 105); a.set_title("Controls (dashed: best constant, equal cost)", fontsize=9, loc="left")
    a.legend(frameon=False, fontsize=8, loc="upper right")
    b = ax[1]
    b.plot(t, 100 * cmp_["X0"][:, 2], color=GREY, lw=1.8, label="no control")
    b.plot(t, 100 * cmp_["Xc"][:, 2], color=INK, lw=1.4, ls="--", label="best constant, equal cost")
    b.plot(t, 100 * sol["X"][:, 2], color=C3, lw=1.8, label="optimal control")
    b.set_ylabel("problematic use $I(t)$ (%)"); b.set_xlabel("year"); b.set_ylim(0, None)
    b.set_title("Prevalence of problematic use", fontsize=9, loc="left")
    b.legend(frameon=False, fontsize=8, loc="lower left")
    fig.tight_layout(); fig.savefig(FIG / "fig_ctl_profiles.png"); plt.close(fig)


def fig_grid(rows):
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.5))
    styles = {"low": (":", "cheaper effort (weights / 4)"), "ref": ("-", "reference weights"),
              "high": ("--", "dearer effort (weights x 4)")}
    for lev, (ls, lab) in styles.items():
        sh = [rows[(lev, r)]["s"]["u1_share_cost_pct"] for r in RATIOS]
        ax[0].plot(RATIOS, sh, color=INK, ls=ls, marker="o", ms=3, label=lab)
        m1 = [100 * rows[(lev, r)]["s"]["u1_mean"] for r in RATIOS]
        m2 = [100 * rows[(lev, r)]["s"]["u2_mean"] for r in RATIOS]
        ax[1].plot(RATIOS, m1, color=C2, ls=ls, marker="o", ms=3)
        ax[1].plot(RATIOS, m2, color=C1, ls=ls, marker="o", ms=3)
    for a in ax:
        a.set_xscale("log", base=2); a.set_xlabel("weight ratio $B_1/B_2$")
        a.axvline(1, color=GREY, lw=0.8)
    ax[0].axhline(50, color=GREY, lw=0.8)
    ax[0].set_ylabel("share of control cost on $u_1$ (%)"); ax[0].set_ylim(0, 100)
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].set_title("Split of effort between the two controls", fontsize=9, loc="left")
    ax[1].set_ylabel("time-averaged control (% of ceiling)"); ax[1].set_ylim(0, 105)
    ax[1].set_title("Exposure ($u_1$, orange) and recovery ($u_2$, blue)", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(FIG / "fig_ctl_weights.png"); plt.close(fig)


def fig_marginal(p1, p2, start):
    t = np.linspace(0, T_HOR, N_STEPS + 1) + start
    fig, ax = plt.subplots(1, 2, figsize=(9.0, 3.4))
    a = ax[0]
    a.plot(t, 1000 * p2, color=C1, lw=1.8, label="recovery effort, $\\chi_2$")
    a.plot(t, 1000 * p1, color=C2, lw=1.8, label="exposure effort, $\\chi_1$")
    a.set_xlabel("year"); a.set_ylabel("person-years averted per 1,000\nper unit of control per year")
    a.set_ylim(0, None); a.legend(frameon=False, fontsize=8)
    a.set_title("Marginal value along the uncontrolled path", fontsize=9, loc="left")
    b = ax[1]
    v = np.array(MV_RATIOS)
    b.hist(v, bins=np.logspace(np.log10(max(v.min() * 0.9, 0.1)), np.log10(v.max() * 1.1), 25), color=GREY, edgecolor="white")
    b.set_xscale("log"); b.axvline(1, color=INK, lw=1, ls="--")
    b.set_xlabel("ratio of marginal values, recovery to exposure ($\\bar{\\chi}_2/\\bar{\\chi}_1$)"); b.set_ylabel("calibrated draws")
    b.set_title("Across calibrated draws, Scenario P", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(FIG / "fig_ctl_marginal.png"); plt.close(fig)


START = None
MV_RATIOS = []
MULTI = []


def main():
    global START
    t0 = time.time()
    START = int(json.load(open(OUT / "results.json"))["int_start_year"]["value"])
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)
    x0 = state_at_start(pc1, x18, START)
    p1, p2 = sec0(pc1, x0)
    sol, cmp_ = sec1(pc1, x0, START)
    _, cmp_["Xc"] = prev_years(x0, pc1, cmp_["c1"], cmp_["c2"])
    fig_profiles(sol, cmp_, START)
    print("1", time.time() - t0, flush=True)
    rows = sec2(pc1, x0)
    fig_grid(rows)
    print("2", time.time() - t0, flush=True)
    ok = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]]
    x0T = state_at_start(pc0, x18, START, alpha_t=cal.shock(pc0.alpha, da))
    sec3(ok, pc0, x0T)
    fig_marginal(p1, p2, START)
    print("3", time.time() - t0, flush=True)
    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_ctl.json", "w"), indent=1, ensure_ascii=False)
    print(f"done, {len(R)} results, {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
