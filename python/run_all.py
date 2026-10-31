"""Thesis v2 pipeline: regenerates every number, table and figure used in the thesis.

    python run_all.py            # full run (about 15 minutes on a laptop)
    python run_all.py --quick    # smoke test

Outputs: outputs/results.json (every reported number, keyed), outputs/*.csv, figures/*.png.
The thesis text is built from results.json, so no number is typed by hand.
"""
from __future__ import annotations

import argparse, json, math, pathlib, platform, time
import numpy as np
import scipy, networkx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from dpm.model import Params, R0, R0_ngm, equilibrium, simulate, elasticity, relapse_multiplier, stable
from dpm import calibrate as cal
from dpm import network as nw
from dpm import recovery as rec

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
OUT.mkdir(exist_ok=True); FIG.mkdir(exist_ok=True)
SEED = 20261031
R = {}  # results registry: key -> value


def put(key, value, note=""):
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    R[key] = {"value": value, "note": note}


def pct(x):
    return 100.0 * x


plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 150, "savefig.bbox": "tight", "font.family": "DejaVu Sans"})
INK = "#1f2933"; C1 = "#2f6db5"; C2 = "#c2571a"; C3 = "#4b8b3b"; GREY = "#9aa5b1"


# ---------------------------------------------------------------- 1. calibration
def run_calibration(ranges, n):
    rows = cal.calibrate(ranges, n=n, seed=SEED)
    cal.save(rows, OUT / "calibration_draws.json")
    ok = [r for r in rows if r["ok"]]
    put("cal_n_draws", len(rows)); put("cal_n_ok", len(ok))
    put("cal_n_rejected_2018", sum(1 for r in rows if r.get("reason") == "2018 level not reachable"))
    for k in ["b0", "b1", "d_alpha", "R0_b0", "R0_b1"]:
        v = np.array([r[k] for r in ok])
        put(f"cal_{k}_median", float(np.median(v)))
        put(f"cal_{k}_p05", float(np.percentile(v, 5)))
        put(f"cal_{k}_p95", float(np.percentile(v, 95)))
    put("cal_share_R0_b0_above1", float(np.mean([r["R0_b0"] > 1 for r in ok])))
    put("cal_share_R0_b0_above1_pct", 100 * float(np.mean([r["R0_b0"] > 1 for r in ok])))
    return ok


def central(ranges):
    d = {k: ranges[k]["central"] for k in cal.PARAM_ORDER}
    base = Params(b=0.0, **d)
    b0 = cal.solve_b0(base, 0.07)
    p0 = base.with_(b=b0)
    x18 = equilibrium(p0)
    b1 = cal.solve_b1(p0, x18, 0.11)
    da = cal.solve_dalpha(p0, x18, 0.11)
    return p0, p0.with_(b=b1), da, x18


# ---------------------------------------------------------------- 2. projections
YEARS = np.arange(2018, 2043)


def project(r, scenario):
    d = {k: r[k] for k in cal.PARAM_ORDER}
    p0 = Params(b=r["b0"], **d)
    x18 = np.array(r["x18"])
    t = YEARS - 2018.0
    if scenario == "P":
        _, y = simulate(p0.with_(b=r["b1"]), x18, t[-1], t_eval=t)
    else:
        _, y = simulate(p0, x18, t[-1], t_eval=t, alpha_t=cal.shock(p0.alpha, r["d_alpha"]))
    return y


def run_projections(ok):
    res = {}
    for sc in ["P", "T"]:
        I = np.array([project(r, sc)[2] for r in ok])
        res[sc] = I
        np.savetxt(OUT / f"projection_{sc}.csv", I, delimiter=",", fmt="%.6f")
        for yr in [2022, 2026, 2027, 2032, 2042]:
            j = yr - 2018
            put(f"proj_{sc}_{yr}_median", pct(np.median(I[:, j])))
            put(f"proj_{sc}_{yr}_p05", pct(np.percentile(I[:, j], 5)))
            put(f"proj_{sc}_{yr}_p95", pct(np.percentile(I[:, j], 95)))
    # implied annual incidence of reaching I among those not in I, 2022 (external check, not fitted)
    inc = []
    for r in ok:
        d = {k: r[k] for k in cal.PARAM_ORDER}
        y = project(r, "P")[:, 2022 - 2018]
        inc.append(d["sigma"] * y[1] / (1 - y[2]))
    put("implied_incidence_median", float(np.median(inc))); put("implied_incidence_p05", float(np.percentile(inc, 5)))
    put("implied_incidence_p95", float(np.percentile(inc, 95)))
    put("proj_gap_2042", R["proj_P_2042_median"]["value"] - R["proj_T_2042_median"]["value"])
    # equilibrium under persistent shift
    Istar = np.array([equilibrium(Params(b=r["b1"], **{k: r[k] for k in cal.PARAM_ORDER}))[2] for r in ok])
    put("eqP_median", pct(np.median(Istar))); put("eqP_p05", pct(np.percentile(Istar, 5))); put("eqP_p95", pct(np.percentile(Istar, 95)))
    # does the persistent scenario overshoot (peak before equilibrium)?
    over = np.mean([I.max() > Istar[i] + 1e-4 for i, I in enumerate(res["P"])])
    put("P_share_overshoot", float(over))

    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    for sc, col, lab in [("P", C1, "Persistent shift in exposure from 2018"), ("T", C2, "Temporary shock, 2020 to 2021")]:
        I = res[sc] * 100
        ax.fill_between(YEARS, np.percentile(I, 5, 0), np.percentile(I, 95, 0), color=col, alpha=0.18, lw=0)
        ax.plot(YEARS, np.median(I, 0), color=col, lw=2, label=lab)
    ax.errorbar([2018, 2022], [7, 11], yerr=0.5, fmt="o", color=INK, capsize=3, label="HBSC published value (rounded)")
    ax.set_xlabel("Year"); ax.set_ylabel("Problematic social media use (%)")
    ax.set_ylim(0, None); ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.savefig(FIG / "fig_calibration_scenarios.png"); plt.close(fig)
    return res


# ---------------------------------------------------------------- 3. mathematics at the central point
def run_mathematics(pc0, pc1):
    put("central_b0", pc0.b); put("central_b1", pc1.b)
    for tag, p in [("b0", pc0), ("b1", pc1)]:
        q = p.with_(alpha=0.0)
        put(f"R0_{tag}", R0(q)); put(f"R0_{tag}_ngm_check", R0_ngm(q))
        put(f"relapse_mult_{tag}", relapse_multiplier(q))
        x = equilibrium(p)
        put(f"eq_{tag}_S", pct(x[0])); put(f"eq_{tag}_A", pct(x[1])); put(f"eq_{tag}_I", pct(x[2])); put(f"eq_{tag}_R", pct(x[3]))
        put(f"eq_{tag}_stable", stable(p, x))
        xa0 = equilibrium(q)
        put(f"eq_{tag}_I_alpha0", pct(xa0[2]))
    # elasticities at the persistent-shift central point
    names = ["b", "eta", "sigma", "gamma1", "gamma2", "delta", "nu", "alpha"]
    el = {}
    for nm in names:
        eR = elasticity(lambda p: R0(p.with_(alpha=0.0)), pc1, nm) if nm != "alpha" else float("nan")
        eI = elasticity(lambda p: equilibrium(p)[2], pc1, nm)
        el[nm] = (eR, eI)
        if nm != "alpha":
            put(f"elast_R0_{nm}", eR)
        put(f"elast_I_{nm}", eI)
    put("elast_I_b_x10", 10 * el["b"][1])
    rankI = sorted(names, key=lambda n: -abs(el[n][1]))
    put("elast_I_rank", rankI)
    labels = {"b": "b exposure", "eta": "η at-risk influence", "sigma": "σ progression", "gamma1": "γ₁ remission",
              "gamma2": "γ₂ recovery", "delta": "δ relapse", "nu": "ν turnover", "alpha": "α exogenous"}
    order = sorted(names, key=lambda n: -abs(el[n][1]))
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    yv = np.arange(len(order))
    ax.barh(yv + 0.2, [el[n][0] if not math.isnan(el[n][0]) else 0 for n in order], 0.4, color=C1, label="R₀ (α = 0)")
    ax.barh(yv - 0.2, [el[n][1] for n in order], 0.4, color=C2, label="Equilibrium prevalence I*")
    ax.set_yticks(yv, [labels[n] for n in order]); ax.invert_yaxis(); ax.axvline(0, color=INK, lw=0.8)
    ax.set_xlabel("Normalised sensitivity (per cent change per 1 per cent change in parameter)")
    ax.legend(frameon=False, fontsize=8); fig.savefig(FIG / "fig_elasticities.png"); plt.close(fig)

    # bifurcation diagram: I* against R0 for alpha = 0 and alpha > 0
    bs = np.linspace(0.01, pc1.b * 2.5, 300)
    r0s = [R0(pc1.with_(b=b, alpha=0.0)) for b in bs]
    I0 = [equilibrium(pc1.with_(b=b, alpha=0.0))[2] * 100 for b in bs]
    Ia = [equilibrium(pc1.with_(b=b))[2] * 100 for b in bs]
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.plot(r0s, I0, color=C1, lw=2, label="α = 0: threshold at R₀ = 1")
    ax.plot(r0s, Ia, color=C2, lw=2, ls="--", label=f"α = {pc1.alpha:.3f} per year: no threshold")
    ax.axvline(1, color=GREY, lw=0.8); ax.set_xlabel("R₀ (exposure rate b varied, other parameters central)")
    ax.set_ylabel("Equilibrium problematic use I* (%)"); ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "fig_bifurcation.png"); plt.close(fig)
    # minimum reachable prevalence when social exposure is removed entirely
    put("I_star_b_zero", pct(equilibrium(pc1.with_(b=0.0))[2]))



# ---------------------------------------------------------------- 3b. local identifiability (sensitivity rank)
def run_identifiability(pc1, x18):
    """Singular values of the scaled sensitivity matrix d log I(t) / d log theta, t = 0..20 years,
    under two observation designs: (i) I observed yearly for 20 years from a known 2018 state,
    (ii) only the two HBSC waves (2018 equilibrium and 2022)."""
    names = ["b", "eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"]
    t = np.arange(0, 21, 1.0)
    def traj(p):
        _, y = simulate(p, x18, 20, t_eval=t); return y[2]
    base = traj(pc1)
    cols = []
    for nm in names:
        v = getattr(pc1, nm); h = 1e-4
        up = traj(pc1.with_(**{nm: v * (1 + h)})); dn = traj(pc1.with_(**{nm: v * (1 - h)}))
        cols.append((up - dn) / (2 * h) / base)
    Sm = np.array(cols).T
    sv = np.linalg.svd(Sm, compute_uv=False)
    put("ident_sv_yearly", [float(x) for x in sv])
    put("ident_rank_yearly_1e-3", int((sv / sv[0] > 1e-3).sum()))
    put("ident_n_params", len(names))

def run_convergence_evidence(n):
    g = np.random.default_rng(SEED + 5)
    worst = 0.0; count = 0
    for _ in range(n):
        p = Params(b=g.uniform(0.2, 8), eta=g.uniform(0, 1), alpha=0.0 if g.random() < 0.5 else g.uniform(1e-4, 0.2),
                   sigma=g.uniform(0.05, 1.5), gamma1=g.uniform(0.05, 1.5), gamma2=g.uniform(0.1, 2), delta=g.uniform(0, 0.6),
                   nu=g.uniform(0.03, 0.1))
        if p.alpha == 0 and abs(R0(p) - 1) < 0.2:
            continue
        x0 = g.dirichlet(np.ones(4))
        _, y = simulate(p, x0, 600, t_eval=[600], rtol=1e-8, atol=1e-10)
        worst = max(worst, float(np.abs(y[:, -1] - equilibrium(p)).max())); count += 1
    put("conv_n_tested", count); put("conv_max_deviation", worst)


# ---------------------------------------------------------------- 4. sensitivity (PRCC)
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


def run_prcc(ok, proj):
    X = np.array([[r[k] for k in cal.PARAM_ORDER] for r in ok])
    for sc in ["P", "T"]:
        y = proj[sc][:, 2042 - 2018]
        pc = prcc(X, y)
        for k, v in zip(cal.PARAM_ORDER, pc):
            put(f"prcc_{sc}_{k}", float(v))


# ---------------------------------------------------------------- 5. interventions
INTERVENTIONS = {
    # name: (description, sampler of parameter multipliers given rng)
    "platform": ("Platform design regulation: exposure rate b reduced", lambda g: {"b": 1 - g.uniform(0.10, 0.40)}),
    "school": ("School-based prevention: remission from At-Risk γ₁ increased", lambda g: {"gamma1": 1 + g.uniform(0.10, 0.50)}),
    "minimum_age": ("Minimum age 16: exposure b and α reduced for the covered share",
                    lambda g: (lambda c: {"b": 1 - 0.4 * c, "alpha": 1 - 0.4 * c})(g.uniform(0.3, 0.7))),
    "treatment": ("Treatment access: recovery γ₂ increased, relapse δ reduced",
                  lambda g: {"gamma2": 1 + g.uniform(0.20, 1.00), "delta": 1 - g.uniform(0.0, 0.30)}),
}
START = 2024


def apply(p: Params, mult: dict) -> Params:
    return p.with_(**{k: getattr(p, k) * v for k, v in mult.items()})


def run_with(r, mults):
    d = {k: r[k] for k in cal.PARAM_ORDER}
    p1 = Params(b=r["b1"], **d)
    x18 = np.array(r["x18"])
    _, y = simulate(p1, x18, START - 2018.0, t_eval=[START - 2018.0])
    x24 = y[:, -1]
    q = p1
    for m in mults:
        q = apply(q, m)
    t = np.arange(0, 2042 - START + 1, 1.0)
    _, y2 = simulate(q, x24, t[-1], t_eval=t)
    return y2[2]  # I from 2024 to 2042


def run_interventions(ok, n):
    g = np.random.default_rng(SEED + 7)
    sub = ok[:n]
    names = list(INTERVENTIONS)
    combos = [(a,) for a in names] + [("platform", "minimum_age"), ("school", "treatment"), ("platform", "school"),
                                       tuple(names)]
    out = {c: [] for c in combos}
    base = []
    for r in sub:
        draws = {k: INTERVENTIONS[k][1](g) for k in names}
        b = run_with(r, [])
        base.append(b)
        for c in combos:
            y = run_with(r, [draws[k] for k in c])
            out[c].append(1 - y[-1] / b[-1])         # relative reduction in I in 2042
    table = []
    single = {}
    for c in combos:
        v = np.array(out[c])
        key = "+".join(c)
        put(f"int_{key}_median", pct(np.median(v))); put(f"int_{key}_p05", pct(np.percentile(v, 5)))
        put(f"int_{key}_p95", pct(np.percentile(v, 95)))
        if len(c) == 1:
            single[c[0]] = v
        table.append((key, np.median(v), np.percentile(v, 5), np.percentile(v, 95)))
    # Bliss independence: expected = 1 - prod(1 - e_i); index = observed / expected
    for c in combos:
        if len(c) > 1:
            exp = 1 - np.prod([1 - single[k] for k in c], axis=0)
            obs = np.array(out[c])
            idx = obs / exp
            key = "+".join(c)
            put(f"int_{key}_bliss_expected_median", pct(np.median(exp)))
            put(f"int_{key}_bliss_index_median", float(np.median(idx)))
            put(f"int_{key}_bliss_index_p05", float(np.percentile(idx, 5)))
            put(f"int_{key}_bliss_index_p95", float(np.percentile(idx, 95)))
    put("int_n_draws", len(sub)); put("int_start_year", START)
    with open(OUT / "interventions.csv", "w") as f:
        f.write("combination,median_reduction_2042,p05,p95\n")
        for k, m, lo, hi in table:
            f.write(f"{k},{m:.4f},{lo:.4f},{hi:.4f}\n")
    labels = {"platform": "Platform design", "school": "School prevention", "minimum_age": "Minimum age 16",
              "treatment": "Treatment access"}
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ks = [t[0] for t in table]
    med = [t[1] * 100 for t in table]; lo = [t[2] * 100 for t in table]; hi = [t[3] * 100 for t in table]
    yv = np.arange(len(ks))
    ax.errorbar(med, yv, xerr=[np.array(med) - lo, np.array(hi) - med], fmt="o", color=C1, capsize=3)
    ax.set_yticks(yv, [" + ".join(labels[x] for x in k.split("+")) if k.count("+") < 3 else "All four" for k in ks], fontsize=8)
    ax.invert_yaxis(); ax.set_xlabel("Reduction in problematic use in 2042 relative to no intervention (%)")
    ax.set_xlim(0, None); fig.savefig(FIG / "fig_interventions.png"); plt.close(fig)


# ---------------------------------------------------------------- 6. networks
def run_networks(pc1, x22, reps, n_nodes):
    k = 10
    beta = nw.beta_from_b(pc1, k)
    put("net_beta_per_month", beta); put("net_mean_degree", k); put("net_nodes", n_nodes); put("net_reps", reps)
    put("net_T", nw.transmissibility(pc1, beta))
    months = 20 * 12
    traj = {}
    for kind in ["ER", "WS", "BA"]:
        runs = []
        k1s, ratios, r0n = [], [], []
        for s in range(reps):
            G = nw.make_graph(kind, n_nodes, k, seed=SEED + s)
            k1, k2 = nw.degree_moments(G); k1s.append(k1); ratios.append(k2 / k1)
            r0n.append(nw.R0_network(pc1.with_(alpha=0.0), beta, G))
            runs.append(nw.simulate(G, pc1, beta, months, x22, seed=SEED + 1000 + s, record_every=12)[:, 2])
        runs = np.array(runs)
        traj[kind] = runs
        put(f"net_{kind}_k2_over_k", float(np.mean(ratios)))
        put(f"net_{kind}_R0net", float(np.mean(r0n)))
        put(f"net_{kind}_I2042_mean", pct(runs[:, -1].mean())); put(f"net_{kind}_I2042_p05", pct(np.percentile(runs[:, -1], 5)))
        put(f"net_{kind}_I2042_p95", pct(np.percentile(runs[:, -1], 95)))
        put(f"net_{kind}_peak_mean", pct(runs.max(1).mean()))
    put("net_meanfield_R0", R0(pc1.with_(alpha=0.0)))
    _, y = simulate(pc1, x22, 20, t_eval=np.arange(0, 21))
    put("net_meanfield_I2042", pct(y[2, -1]))
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    yrs = np.arange(2022, 2043)
    for kind, col in [("ER", C1), ("WS", C3), ("BA", C2)]:
        r = traj[kind] * 100
        ax.fill_between(yrs, np.percentile(r, 5, 0), np.percentile(r, 95, 0), color=col, alpha=0.15, lw=0)
        ax.plot(yrs, r.mean(0), color=col, lw=2, label={"ER": "Erdős–Rényi", "WS": "Watts–Strogatz", "BA": "Barabási–Albert"}[kind])
    ax.plot(yrs, y[2] * 100, color=INK, lw=1.2, ls="--", label="Mean-field ODE")
    ax.set_xlabel("Year"); ax.set_ylabel("Problematic use (%)"); ax.set_ylim(0, None); ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIG / "fig_network_topologies.png"); plt.close(fig)

    # targeting on BA networks: protect 10% of nodes (exposure x0.2) chosen by degree or at random
    cover, damp = 0.10, 0.2
    res = {"none": [], "random": [], "hub": []}
    for s in range(reps):
        G = nw.make_graph("BA", n_nodes, k, seed=SEED + s)
        deg = np.array([d for _, d in G.degree()])
        rng = np.random.default_rng(SEED + 2000 + s)
        m = int(cover * n_nodes)
        for strat in res:
            w = np.ones(n_nodes)
            if strat == "random":
                w[rng.choice(n_nodes, m, replace=False)] = damp
            elif strat == "hub":
                w[np.argsort(-deg)[:m]] = damp
            res[strat].append(simulate_weighted(G, pc1, beta, months, x22, w, SEED + 3000 + s))
    for strat, v in res.items():
        v = np.array(v)
        put(f"net_target_{strat}_I2042_mean", pct(v.mean()))
    put("net_target_cover", cover); put("net_target_damp", damp)
    put("net_target_cover_pct", 100 * cover); put("net_target_damp_pct", 100 * damp)


def simulate_weighted(G, p, beta, months, x0, w, seed):
    """Network run in which node i's susceptibility to social exposure is scaled by w[i] (targeting experiment)."""
    import scipy.sparse as sp
    rng = np.random.default_rng(seed)
    n = G.number_of_nodes()
    Adj = networkx.to_scipy_sparse_array(G, format="csr", dtype=float)
    # scaling the edge weights into node i is equivalent to scaling log(1-beta) for node i
    st = rng.choice(4, size=n, p=np.asarray(x0) / np.sum(x0))
    DT = 1 / 12
    am = 1 - math.exp(-p.alpha * DT)
    pA = 1 - math.exp(-(p.sigma + p.gamma1) * DT); shareI = p.sigma / (p.sigma + p.gamma1)
    pI = 1 - math.exp(-p.gamma2 * DT); pR = 1 - math.exp(-p.delta * DT); pexit = 1 - math.exp(-p.nu * DT)
    l = math.log(1 - beta)
    for t in range(months):
        expo = Adj @ ((st == 2) + p.eta * (st == 1))
        pinf = 1 - (1 - am) * np.exp(l * expo * w)
        u = rng.random(n); new = st.copy()
        S, A, I, Rr = st == 0, st == 1, st == 2, st == 3
        new[S & (u < pinf)] = 1
        la = A & (u < pA); new[la & (u < pA * shareI)] = 2; new[la & (u >= pA * shareI)] = 0
        new[I & (u < pI)] = 3; new[Rr & (u < pR)] = 1
        new[rng.random(n) < pexit] = 0
        st = new
    return float(np.mean(st == 2))


# ---------------------------------------------------------------- 7. HMM recovery study
def run_recovery(sizes, reps):
    g = np.random.default_rng(SEED + 11)
    truth = dict(sigma=0.40, gamma1=0.50, gamma2=0.60, delta=0.15, incidence=0.08)
    P = rec.annual_matrix(**truth)
    pi0 = np.array([0.75, 0.12, 0.08, 0.05])
    for k, v in truth.items():
        put(f"rec_truth_{k}", v)
    put("rec_reps", reps)
    put("rec_waves", 6)
    entries = [(0, 1), (1, 0), (1, 2), (2, 3), (3, 1)]
    names = ["S→A", "A→S", "A→I", "I→R", "R→A"]
    rows = []
    for noise, sd in [("low", np.array([0.8, 1.0, 1.0, 0.8])), ("high", np.array([1.5, 1.8, 1.8, 1.5]))]:
        mu = np.array([1.0, 3.5, 7.0, 1.0])
        for n in sizes:
            est = []
            for r in range(reps):
                _, y = rec.generate(P, pi0, mu, sd, n, 6, g)
                f = rec.fit_multistart(y, P, pi0, mu, sd, g, starts=2)
                est.append([f["P"][i, j] for i, j in entries])
            est = np.array(est); true = np.array([P[i, j] for i, j in entries])
            bias = est.mean(0) - true; rmse = np.sqrt(((est - true) ** 2).mean(0))
            for nm, t, b, e in zip(names, true, bias, rmse):
                rows.append((noise, n, nm, t, b, e))
                put(f"rec_{noise}_{n}_{nm.replace(chr(8594), chr(95))}_relrmse", float(e / t))
    with open(OUT / "recovery.csv", "w") as f:
        f.write("noise,n,transition,true,bias,rmse\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for noise, ls in [("low", "-"), ("high", "--")]:
        for nm, col in zip(names, [C1, C2, C3, INK, GREY]):
            xs = [r[1] for r in rows if r[0] == noise and r[2] == nm]
            ys = [r[5] / r[3] * 100 for r in rows if r[0] == noise and r[2] == nm]
            ax.plot(xs, ys, ls=ls, color=col, marker="o", ms=3, label=f"{nm}" if noise == "low" else None)
    ax.set_xscale("log"); ax.set_xlabel("Panel size (young people, six annual waves)")
    ax.set_ylabel("Relative RMSE (%)"); ax.legend(frameon=False, fontsize=8, title="solid: low noise; dashed: high noise", title_fontsize=7)
    fig.savefig(FIG / "fig_recovery.png"); plt.close(fig)


# ---------------------------------------------------------------- 8. South African access scenarios
def run_sa(pc0, pc1, access):
    """Scenario only: exposure (b and alpha) scaled by the connected share. Not an estimate for South Africa."""
    out = {}
    for name, f in access.items():
        q = pc1.with_(b=pc1.b * f, alpha=pc1.alpha * f)
        x = equilibrium(q)
        put(f"sa_{name}_access", f); put(f"sa_{name}_Istar", pct(x[2]))
        put(f"sa_{name}_R0", R0(q.with_(alpha=0.0)))
        out[name] = x[2]
    fs = np.linspace(0.05, 1.0, 60)
    Is = [equilibrium(pc1.with_(b=pc1.b * f, alpha=pc1.alpha * f))[2] * 100 for f in fs]
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot(fs * 100, Is, color=C1, lw=2)
    for name, f in access.items():
        ax.plot(f * 100, out[name] * 100, "o", color=C2)
        ax.annotate(name.replace("_", " "), (f * 100, out[name] * 100), textcoords="offset points", xytext=(4, -10), fontsize=8)
    ax.set_xlabel("Connected share of young people (%)"); ax.set_ylabel("Equilibrium problematic use (%)")
    fig.savefig(FIG / "fig_sa_access.png"); plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); a = ap.parse_args()
    t0 = time.time()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    access = json.load(open(HERE / "data" / "sa_access.json"))
    n_cal = 60 if a.quick else 1000
    ok = run_calibration(ranges, n_cal)
    proj = run_projections(ok)
    pc0, pc1, da_c, x18 = central(ranges)
    put("central_d_alpha", da_c)
    run_mathematics(pc0, pc1)
    run_identifiability(pc1, x18)
    run_convergence_evidence(200 if a.quick else 5000)
    run_prcc(ok, proj)
    run_interventions(ok, 30 if a.quick else 400)
    _, y = simulate(pc1, x18, 4, t_eval=[4]); x22 = y[:, -1]
    put("x2022_central", [float(v) for v in x22])
    run_networks(pc1, x22, 4 if a.quick else 40, 1000 if a.quick else 2000)
    run_recovery([500, 2000] if a.quick else [500, 2000, 8000], 2 if a.quick else 20)
    run_sa(pc0, pc1, {k: v["value"] for k, v in access.items() if not k.startswith("_")})
    put("env_python", platform.python_version()); put("env_numpy", np.__version__); put("env_scipy", scipy.__version__)
    put("env_networkx", networkx.__version__); put("env_matplotlib", matplotlib.__version__); put("seed", SEED)
    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results.json", "w"), indent=1, ensure_ascii=False)
    print(f"done in {(time.time() - t0) / 60:.1f} min, {len(R)} results")


if __name__ == "__main__":
    main()
