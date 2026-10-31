"""Non-exponential stage durations (Chapter 5, section "Stage durations that are not exponential").

The At-Risk, Addicted and Recovered spells are given Erlang (gamma with integer shape) durations with the
same means as in the main model, by the linear chain trick. Ageing out remains exponential at rate nu.
The script

1. checks the Laplace-transform formulas for R0 and the equilibrium against the stage ODE,
2. re-calibrates the two exposure quantities (Scenario P) and the shock (Scenario T) for each shape,
3. compares projections, time scales and the floor across shapes.

Writes outputs/results_durx.json (keys durx_...) and figures/fig_durx_*.png.
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt
from dpm.model import Params, R0 as R0_exp, equilibrium as eq_exp
from dpm import calibrate as cal

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R: dict = {}


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["durx_" + k] = {"value": v, "note": note}


# ---------------------------------------------------------------- Laplace-transform quantities
def stage_quantities(p: Params, kA: int, kI: int, kR: int):
    """Probabilities of completing each spell before ageing out (Laplace transform of the Erlang
    duration at nu) and expected time spent in each spell, for spells with means 1/(sigma+gamma1),
    1/gamma2 and 1/delta."""
    rA, rI = p.sigma + p.gamma1, p.gamma2
    pA = (kA * rA / (kA * rA + p.nu)) ** kA
    pI = (kI * rI / (kI * rI + p.nu)) ** kI
    if p.delta > 0:
        pR = (kR * p.delta / (kR * p.delta + p.nu)) ** kR
    else:
        pR = 0.0
    tA, tI, tR = (1 - pA) / p.nu, (1 - pI) / p.nu, (1 - pR) / p.nu
    piA = pA * p.sigma / rA          # probability an At-Risk spell ends in problematic use
    q = piA * pI * pR                # probability a career loops back to At-Risk
    K = p.eta * tA + piA * tI        # exposure per unit b produced by one entry into A, single spell
    M = tA + piA * tI + piA * pI * tR  # expected time in A, I, R per entry into A, single spell
    return dict(pA=pA, pI=pI, pR=pR, tA=tA, tI=tI, tR=tR, piA=piA, q=q, K=K, M=M)


def R0_dur(p, kA, kI, kR):
    Q = stage_quantities(p, kA, kI, kR)
    return p.b * Q["K"] / (1 - Q["q"])


def eq_flux(p, kA, kI, kR):
    """Equilibrium entry flux J into the At-Risk state: root of b K M J^2 - (b K - alpha M - (1-q)) J - alpha = 0."""
    Q = stage_quantities(p, kA, kI, kR)
    kap, M, om = p.b * Q["K"], Q["M"], 1 - Q["q"]
    if p.alpha == 0:
        r0 = kap / om
        return 0.0 if r0 <= 1 else (1 - 1 / r0) / M, Q
    if kap == 0:
        return p.alpha / (p.alpha * M + om), Q
    B = kap - p.alpha * M - om
    return (B + math.sqrt(B * B + 4 * kap * M * p.alpha)) / (2 * kap * M), Q


def I_eq_dur(p, kA, kI, kR):
    J, Q = eq_flux(p, kA, kI, kR)
    return J * Q["piA"] * Q["tI"]


# ---------------------------------------------------------------- stage ODE (linear chain trick)
class Chain:
    def __init__(self, kA, kI, kR):
        self.k = (kA, kI, kR)
        self.n = 1 + kA + kI + kR
        self.iA = slice(1, 1 + kA); self.iI = slice(1 + kA, 1 + kA + kI); self.iR = slice(1 + kA + kI, self.n)

    def rhs(self, p, b_t=None, alpha_t=None):
        kA, kI, kR = self.k
        rA, rI, rR = kA * (p.sigma + p.gamma1), kI * p.gamma2, kR * p.delta
        fs = p.sigma / (p.sigma + p.gamma1)
        iA, iI, iR = self.iA, self.iI, self.iR

        def f(t, x):
            b = p.b if b_t is None else b_t(t)
            al = p.alpha if alpha_t is None else alpha_t(t)
            S = x[0]; A = x[iA]; I = x[iI]; Rr = x[iR]
            lam = b * (I.sum() + p.eta * A.sum())
            inc = (lam + al) * S
            dx = np.empty_like(x)
            outA = rA * A[-1]
            dx[0] = p.nu - inc - p.nu * S + (1 - fs) * outA
            dA = -(rA + p.nu) * A
            dA[0] += inc + rR * Rr[-1]
            dA[1:] += rA * A[:-1]
            dI = -(rI + p.nu) * I
            dI[0] += fs * outA
            dI[1:] += rI * I[:-1]
            dR = -(rR + p.nu) * Rr
            dR[0] += rI * I[-1]
            dR[1:] += rR * Rr[:-1]
            dx[iA] = dA; dx[iI] = dI; dx[iR] = dR
            return dx
        return f

    def eq_state(self, p):
        """Equilibrium of the stage ODE built from the flux solution."""
        kA, kI, kR = self.k
        J, Q = eq_flux(p, kA, kI, kR)
        rA, rI, rR = kA * (p.sigma + p.gamma1), kI * p.gamma2, kR * p.delta
        x = np.zeros(self.n)
        def chain(Jin, r, k):
            out = np.empty(k); v = Jin / (r + p.nu)
            for j in range(k):
                out[j] = v; v = v * r / (r + p.nu)
            return out
        x[self.iA] = chain(J, rA, kA)
        JI = J * Q["piA"]; x[self.iI] = chain(JI, rI, kI)
        JR = JI * Q["pI"]; x[self.iR] = chain(JR, rR, kR)
        x[0] = 1 - x[1:].sum()
        return x

    def I(self, x):
        return x[self.iI].sum(axis=0)

    def sim(self, p, x0, T, t_eval, b_t=None, alpha_t=None):
        sol = solve_ivp(self.rhs(p, b_t, alpha_t), (0, T), x0, t_eval=t_eval, method="LSODA",
                        rtol=1e-9, atol=1e-11, max_step=0.25)
        assert sol.success
        return sol.t, sol.y


def jac_num(f, x, h=1e-7):
    n = len(x); J = np.empty((n, n)); f0 = f(0, x)
    for j in range(n):
        e = np.zeros(n); e[j] = h
        J[:, j] = (f(0, x + e) - f(0, x - e)) / (2 * h)
    return J


# ---------------------------------------------------------------- calibration for one shape
def calibrate_shape(base: Params, k, t18=0.07, t22=0.11):
    ch = Chain(*k)
    b0 = brentq(lambda b: I_eq_dur(base.with_(b=b), *k) - t18, 0.0, 200, xtol=1e-11)
    p0 = base.with_(b=b0)
    x18 = ch.eq_state(p0)
    def I4(p, alpha_t=None):
        _, y = ch.sim(p, x18, 4.0, [4.0], alpha_t=alpha_t)
        return float(ch.I(y)[-1])
    b1 = brentq(lambda b: I4(p0.with_(b=b)) - t22, b0, 200, xtol=1e-10)
    da = brentq(lambda d: I4(p0, cal.shock(p0.alpha, d)) - t22, 0.0, 20, xtol=1e-11)
    return ch, p0, x18, b1, da


def analyse_shape(base, k, tag, store_path=False):
    ch, p0, x18, b1, da = calibrate_shape(base, k)
    p1 = p0.with_(b=b1)
    tf = np.linspace(0, 24, 24 * 12 + 1); yrs = 2018 + tf
    _, yP = ch.sim(p1, x18, 24.0, tf)
    _, yT = ch.sim(p0, x18, 24.0, tf, alpha_t=cal.shock(p0.alpha, da))
    IPf, ITf = ch.I(yP), ch.I(yT)
    IP, IT = IPf[::12], ITf[::12]
    Ieq1 = I_eq_dur(p1, *k)
    # time from 2018 until within 10% of the remaining distance to the new equilibrium (Scenario P)
    tt = np.linspace(0, 60, 6001)
    _, yy = ch.sim(p1, x18, 60.0, tt)
    dev = np.abs(ch.I(yy) - Ieq1) / abs(Ieq1 - 0.07)
    viol = np.where(dev > 0.1)[0]
    t90 = float(tt[viol[-1] + 1]) if len(viol) else 0.0
    peak = float(ch.I(yy).max())
    J = jac_num(ch.rhs(p1), ch.eq_state(p1))
    # eliminate S = 1 - (sum of the other states): this removes the mode of the total population (rate nu)
    T = np.vstack([-np.ones((1, ch.n - 1)), np.eye(ch.n - 1)])
    ev = np.linalg.eigvals((J @ T)[1:, :])
    slow = float(-np.max(ev.real))
    floor = I_eq_dur(p1.with_(b=0.0), *k)
    r0 = (R0_dur(p0.with_(alpha=0.0), *k), R0_dur(p1.with_(alpha=0.0), *k))
    out = dict(b0=p0.b, b1=b1, da=da, R0_b0=r0[0], R0_b1=r0[1], P2042=100 * IP[-1], T2042=100 * IT[-1],
               P2026=100 * IP[8], T2026=100 * IT[8], gap2026=100 * (IP[8] - IT[8]),
               T_peak=100 * ITf.max(), T_peak_year=float(yrs[np.argmax(ITf)]),
               Ieq_b1=100 * Ieq1, t90=t90, peak=100 * peak, overshoot=100 * (peak - Ieq1), slow=slow,
               halflife=math.log(2) / slow, floor=100 * floor, q=stage_quantities(p1, *k)["q"])
    for key, v in out.items():
        put(f"{tag}_{key}", float(v))
    if store_path:
        return out, yrs, IPf, ITf
    return out


def main():
    t0 = time.time()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da_c, x18_c = ra.central(ranges)
    base = pc0.with_(b=0.0)

    # 1. checks: k = 1 reproduces the main model; flux equilibrium is a zero of the stage ODE; long-run convergence
    g = np.random.default_rng(31)
    err_R0 = err_I = err_fix = err_conv = 0.0
    for i in range(200):
        p = Params(b=g.uniform(0.2, 6), eta=g.uniform(0, 1), alpha=(0.0 if i % 2 == 0 else g.uniform(1e-3, 0.1)),
                   sigma=g.uniform(0.05, 1.5), gamma1=g.uniform(0.05, 1.5), gamma2=g.uniform(0.1, 2),
                   delta=g.uniform(0.01, 0.6), nu=g.uniform(0.03, 0.1))
        err_R0 = max(err_R0, abs(R0_dur(p, 1, 1, 1) - R0_exp(p)) / R0_exp(p))
        err_I = max(err_I, abs(I_eq_dur(p, 1, 1, 1) - eq_exp(p)[2]))
        k = tuple(int(v) for v in g.integers(1, 9, size=3))
        ch = Chain(*k)
        x = ch.eq_state(p)
        err_fix = max(err_fix, float(np.max(np.abs(ch.rhs(p)(0, x)))))
        if i < 40 and (p.alpha > 0 or R0_dur(p, *k) > 1.2):
            x0 = np.zeros(ch.n); x0[0] = 0.9; x0[ch.iA.start] = 0.05; x0[ch.iI.start] = 0.05
            _, y = ch.sim(p, x0, 800.0, [800.0])
            err_conv = max(err_conv, abs(float(ch.I(y)[-1]) - I_eq_dur(p, *k)))
    put("chk_R0_k1_max_rel_err", err_R0); put("chk_I_k1_max_abs_err", err_I)
    put("chk_fixed_point_max_resid", err_fix); put("chk_convergence_max_abs_err", err_conv)
    print("checks", time.time() - t0, flush=True)

    # 2. completion probabilities at the central point for shapes 1 and infinity
    for nm, mean in [("A", 1 / (pc1.sigma + pc1.gamma1)), ("I", 1 / pc1.gamma2), ("R", 1 / pc1.delta)]:
        put(f"mean_{nm}_years", mean)
        put(f"pcomplete_{nm}_k1", 1 / (1 + pc1.nu * mean))
        put(f"pcomplete_{nm}_kinf", math.exp(-pc1.nu * mean))

    # 3. shapes applied to all three spells
    shapes = [1, 2, 3, 5, 10, 20]
    rows = {}
    paths = {}
    for s in shapes:
        out, yrs, IP, IT = analyse_shape(base, (s, s, s), f"all{s}", store_path=True)
        rows[s] = out; paths[s] = (IP, IT)
        print("shape", s, time.time() - t0, flush=True)
    # deterministic limit for R0 and equilibria (no simulation)
    big = 2000
    b0 = brentq(lambda b: I_eq_dur(base.with_(b=b), big, big, big) - 0.07, 0.0, 200, xtol=1e-11)
    put("allinf_b0", b0); put("allinf_R0_b0", R0_dur(base.with_(b=b0, alpha=0.0), big, big, big))
    # one spell at a time
    for nm, k in [("A", (10, 1, 1)), ("I", (1, 10, 1)), ("R", (1, 1, 10))]:
        analyse_shape(base, k, f"only{nm}10")
        print("only", nm, time.time() - t0, flush=True)
    # summary ranges across shapes
    for key in ["P2042", "T2042", "t90", "halflife", "floor", "R0_b1", "overshoot", "gap2026", "da", "T_peak"]:
        v = [rows[s][key] for s in shapes]
        put(f"range_{key}_min", min(v)); put(f"range_{key}_max", max(v))
    put("shapes", shapes)

    # figure: trajectories under P and T for shapes 1, 3, 20
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
    cols = {1: INK, 3: C1, 20: C2}
    for s in [1, 3, 20]:
        IP, IT = paths[s]
        ax[0].plot(yrs, 100 * IP, color=cols[s], lw=1.8, label=f"shape {s}")
        ax[1].plot(yrs, 100 * IT, color=cols[s], lw=1.8, label=f"shape {s}")
    for a_, title in zip(ax, ["Scenario P (persistent shift)", "Scenario T (shock in 2020 and 2021)"]):
        a_.scatter([2018, 2022], [7, 11], color=GREY, zorder=5, s=25)
        a_.set_title(title, fontsize=10); a_.set_xlabel("Year"); a_.grid(alpha=0.3)
        a_.set_xticks([2018, 2022, 2026, 2030, 2034, 2038, 2042])
    ax[0].set_ylabel("Problematic use (%)"); ax[0].legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(FIG / "fig_durx_paths.png", dpi=200); plt.close(fig)

    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_durx.json", "w"), indent=1, ensure_ascii=False)
    print("done", len(R), (time.time() - t0) / 60)


if __name__ == "__main__":
    main()
