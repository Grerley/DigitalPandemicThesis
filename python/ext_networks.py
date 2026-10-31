"""Chapter 8 extension: further network analyses (thesis Chapter 8, sections added in the 2026 expansion).

    python ext_networks.py            # full run (about 10 to 15 minutes on a laptop)
    python ext_networks.py --quick    # smoke test

Writes outputs/results_netx.json (keys prefixed netx_) and figures/fig_netx_*.png. Uses the calibrated
post-2018 central parameter set of run_all.central and the house plotting style of run_all.

Analyses
  1. Threshold sweep: long-run outcome against the heuristic R_net (8.4) and mean-field R0 on ER, WS and BA
     graphs with alpha = 0; finite-size check on BA graphs with 1,000, 2,000 and 8,000 nodes.
  2. Homophily without contagion: ties assortative on a latent trait that raises exogenous acquisition and
     progression, b = 0; clustering and a naive two-wave peer-effect regression, against a contagion
     world on the same graphs.
  3. Complex contagion: social exposure acts only with at least two affected neighbours; invasion from a
     local seed on ER and WS graphs, and long-run level from low and high starting states.
  4. Targeting: paired runs on BA graphs, coverage 5% to 30%, four distances from the threshold, random,
     highest-degree and acquaintance selection.
  5. Mean degree: 2042 level against mean degree, and a decomposition of the exposure deficit of
     Susceptible nodes relative to the mean-field model.
All outputs are model outputs.
"""
from __future__ import annotations

import argparse, json, math, time
import numpy as np
import networkx as nx
from scipy.optimize import brentq

import run_all as ra                                   # house style (rcParams), SEED, central(), _b_for_Rnet
from run_all import INK, C1, C2, C3, GREY, SEED, HERE, OUT, FIG, plt
from dpm import network as nw
from dpm import network_ext as nx2
from dpm.model import R0, equilibrium, simulate

R = {}
NAMES = {"ER": "Erdős–Rényi", "WS": "Watts–Strogatz", "BA": "Barabási–Albert"}
COLS = {"ER": C1, "WS": C3, "BA": C2}
K = 10


def put(key, value, note=""):
    assert key.startswith("netx_"), key
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    R[key] = {"value": value, "note": note}


def tag(v):
    return f"{v:g}".replace(".", "p")


def S(offset):
    """Seed block for this script, disjoint from the offsets used in run_all.py (below 10,000)."""
    return SEED + 100_000 + offset


# ---------------------------------------------------------------- 1. threshold sweep
THR_GRID = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5]


def _sweep(pc1, kind, n, grid, reps, years, seed0):
    """For each target R_net: runs with alpha = 0 from 5% Addicted. Returns per-target summaries."""
    x0 = np.array([0.95, 0.0, 0.05, 0.0])
    out = {}
    graphs = [nw.make_graph(kind, n, K, seed=seed0 + s) for s in range(reps)]
    adjs = [nx.to_scipy_sparse_array(G, format="csr", dtype=float) for G in graphs]
    ratio = np.mean([nw.degree_moments(G)[1] / nw.degree_moments(G)[0] for G in graphs])
    for tgt in grid:
        tot_end, tot_mid, I_end, alive, r0mf, bs = [], [], [], [], [], []
        for s, (G, Adj) in enumerate(zip(graphs, adjs)):
            b = ra._b_for_Rnet(pc1, G, K, tgt)
            q = pc1.with_(b=b, alpha=0.0)
            tr = nx2.simulate_ext(G, q, nw.beta_from_b(q, K), years * 12, x0, seed=seed0 + 50_000 + s,
                                  record_every=12, Adj=Adj)["traj"]
            tot = tr[:, 1] + tr[:, 2] + tr[:, 3]
            tot_end.append(tot[-11:].mean()); tot_mid.append(tot[years // 2 - 5: years // 2 + 6].mean())
            I_end.append(tr[-11:, 2].mean()); alive.append(tot[-1] > 0); r0mf.append(R0(q)); bs.append(b)
        out[tgt] = dict(tot_end=float(np.mean(tot_end)), tot_mid=float(np.mean(tot_mid)), I_end=float(np.mean(I_end)),
                        persist=float(np.mean(alive)), R0mf=float(np.mean(r0mf)), b=float(np.mean(bs)))
    return out, float(ratio)


def eff_threshold(res):
    """Effective threshold: linear extrapolation to zero of the long-run share At-Risk, Addicted or Recovered,
    using the three lowest grid values at which that share is at least 3% and not declining (the mean over
    the last 11 years is at least 95% of the mean over the 11 years around the midpoint)."""
    g = sorted(res)
    ok = [t for t in g if res[t]["tot_end"] >= 0.03 and res[t]["tot_end"] >= 0.95 * res[t]["tot_mid"]][:3]
    if len(ok) < 2:
        return float("nan"), float("nan"), float("nan")
    x = np.array(ok); y = np.array([res[t]["tot_end"] for t in ok])
    slope, icpt = np.polyfit(x, y, 1)
    rstar = -icpt / slope
    r0 = float(np.interp(rstar, g, [res[t]["R0mf"] for t in g]))
    b = float(np.interp(rstar, g, [res[t]["b"] for t in g]))
    return float(rstar), r0, b


def run_threshold(pc1, quick):
    reps = 3 if quick else 16
    years = 60 if quick else 100
    grid = [0.6, 0.9, 1.2, 1.5, 2.0] if quick else THR_GRID
    put("netx_thr_reps", reps); put("netx_thr_years", years); put("netx_thr_grid_min", min(grid)); put("netx_thr_grid_max", max(grid))
    put("netx_thr_grid_n", len(grid))
    allres = {}
    for j, kind in enumerate(["ER", "WS", "BA"]):
        res, ratio = _sweep(pc1, kind, 2000, grid, reps, years, S(1000 * j))
        allres[kind] = res
        for t, r in res.items():
            tg = f"netx_thr_{kind}_{tag(t)}"
            put(f"{tg}_tot_pct", 100 * r["tot_end"]); put(f"{tg}_I_pct", 100 * r["I_end"])
            put(f"{tg}_persist_pct", 100 * r["persist"]); put(f"{tg}_R0mf", r["R0mf"])
        rstar, r0, b = eff_threshold(res)
        put(f"netx_thr_{kind}_Rnet_star", rstar, "effective threshold in R_net units (extrapolation)")
        put(f"netx_thr_{kind}_R0mf_star", r0, "mean-field R0 at the effective threshold")
        put(f"netx_thr_{kind}_bstar_over_b1", b / pc1.b)
        put(f"netx_thr_{kind}_bstar", b)
        put(f"netx_thr_{kind}_k2_over_k", ratio)
    # mean-field reference: ODE threshold b and equilibrium curve
    b_mf1 = pc1.b / R0(pc1.with_(alpha=0.0))
    put("netx_thr_meanfield_b1_over_b1", b_mf1 / pc1.b)

    # figure: long-run problematic use against R_net and against mean-field R0
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3), sharey=True)
    for kind in ["ER", "WS", "BA"]:
        res = allres[kind]; g = sorted(res)
        axs[0].plot(g, [100 * res[t]["I_end"] for t in g], "o-", ms=3, color=COLS[kind], lw=1.5, label=NAMES[kind])
        axs[1].plot([res[t]["R0mf"] for t in g], [100 * res[t]["I_end"] for t in g], "o-", ms=3, color=COLS[kind], lw=1.5)
    r0s = np.linspace(0.2, 4.2, 200)
    Istar = [100 * equilibrium(pc1.with_(b=pc1.b * r / R0(pc1.with_(alpha=0.0)), alpha=0.0))[2] for r in r0s]
    axs[1].plot(r0s, Istar, color=INK, lw=1.2, ls="--", label="Mean-field ODE")
    for ax in axs:
        ax.axvline(1.0, color=GREY, lw=0.8, ls=":")
        ax.set_ylim(0, None)
    axs[0].set_xlabel(r"Heuristic $\mathcal{R}_{\mathrm{net}}$ (8.4)"); axs[1].set_xlabel(r"Mean-field $\mathcal{R}_0$ at the same exposure rate")
    axs[0].set_ylabel("Problematic use, years 90 to 100 (%)")
    axs[0].legend(frameon=False, fontsize=8); axs[1].legend(frameon=False, fontsize=8, loc="lower right")
    axs[0].set_title("A. Against the network heuristic", fontsize=9, loc="left")
    axs[1].set_title("B. Against the mean-field threshold", fontsize=9, loc="left")
    fig.savefig(FIG / "fig_netx_threshold.png"); plt.close(fig)

    # finite-size check on BA graphs
    fs_grid = [0.9, 1.2, 1.5] if quick else [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0]
    fs_reps = 3 if quick else 12
    sizes = [1000, 2000, 8000]
    put("netx_fs_reps", fs_reps)
    for i, n in enumerate(sizes):
        if n == 2000:
            res = {t: allres["BA"][t] for t in fs_grid if t in allres["BA"]}
            ratio = R["netx_thr_BA_k2_over_k"]["value"]
        else:
            res, ratio = _sweep(pc1, "BA", n, fs_grid, fs_reps, years, S(5000 + 1000 * i))
        rstar, r0, b = eff_threshold(res)
        put(f"netx_fs_n{n}_k2_over_k", ratio)
        put(f"netx_fs_n{n}_Rnet_star", rstar); put(f"netx_fs_n{n}_R0mf_star", r0); put(f"netx_fs_n{n}_bstar_over_b1", b / pc1.b)
        # exposure rate at which the heuristic equals one, and its mean-field R0
        Gs = [nw.make_graph("BA", n, K, seed=S(9000 + 100 * i + s)) for s in range(5)]
        b1s = [ra._b_for_Rnet(pc1, G, K, 1.0) for G in Gs]
        put(f"netx_fs_n{n}_b_Rnet1_over_b1", float(np.mean(b1s)) / pc1.b)
        put(f"netx_fs_n{n}_R0mf_at_Rnet1", float(np.mean([R0(pc1.with_(b=b, alpha=0.0)) for b in b1s])))
        for t in [1.0, 1.5]:
            if t in res:
                put(f"netx_fs_n{n}_{tag(t)}_tot_pct", 100 * res[t]["tot_end"])
                put(f"netx_fs_n{n}_{tag(t)}_persist_pct", 100 * res[t]["persist"])
    return allres


# ---------------------------------------------------------------- 2. homophily without contagion
def edge_odds_ratio(G, I):
    e = np.array(G.edges()); a, b = I[e[:, 0]], I[e[:, 1]]
    both = np.sum(a & b); neither = np.sum(~a & ~b); dis = np.sum(a ^ b)
    return 4.0 * both * neither / dis ** 2


def survey_rows(Adj, snaps, t1, t2):
    """Two waves 12 months apart. Sample: the same people at both waves, without problematic use at wave 1.
    Outcome: problematic use at wave 2. Exposure: number of friends with problematic use at wave 1."""
    s1, id1 = snaps[t1]; s2, id2 = snaps[t2]
    I1 = s1 == 2
    nI = np.asarray(Adj @ I1.astype(float))
    sel = (id1 == id2) & ~I1
    return sel, nI, (s2 == 2)


def _fit(y, nI, z=None):
    X = [np.ones(len(y)), nI]
    if z is not None:
        X += [z, z ** 2, z ** 3]
    b, se = nx2.logit_fit(np.column_stack(X), y)
    return b[1], se[1]


def run_homophily(pc1, quick):
    n = 2000; h = 0.6; groups = 20
    reps = 4 if quick else 20
    s_grid = [0.0, 0.5, 1.0] if quick else [0.0, 0.25, 0.5, 0.75, 1.0]
    s_head = 0.5
    t1, t2 = 360, 372
    put("netx_hom_n", n); put("netx_hom_h_pct", 100 * h); put("netx_hom_groups", groups); put("netx_hom_reps", reps)
    put("netx_hom_s_head", s_head); put("netx_hom_burnin_years", t1 // 12)
    c_ = pc1.gamma2 + pc1.nu; d_ = pc1.delta + pc1.nu
    m_ = 1 + pc1.sigma / c_ + pc1.sigma * pc1.gamma2 / (c_ * d_)
    put("netx_hom_I_ceiling_pct", 100 * pc1.sigma / (c_ * m_), "limit of (5.10) as alpha grows, b = 0, central sigma")
    beta = nw.beta_from_b(pc1, K)
    graphs, zs, adjs, rz = [], [], [], []
    for r in range(reps):
        z = np.random.default_rng(S(20_000 + r)).standard_normal(n)
        G = nx2.homophilous_graph(z, K, h, groups, seed=S(21_000 + r))
        graphs.append(G); zs.append(z); adjs.append(nx.to_scipy_sparse_array(G, format="csr", dtype=float))
        e = np.array(G.edges()); rz.append(np.corrcoef(z[e[:, 0]], z[e[:, 1]])[0, 1])
    put("netx_hom_trait_corr_friends", float(np.mean(rz)))
    put("netx_hom_mean_degree", float(np.mean([2 * G.number_of_edges() / n for G in graphs])))

    def collect(world_runs):
        """world_runs: list of (G, Adj, z, snaps). Returns summary dict."""
        Y, NI, Z, ors, prev, est, sig, grad = [], [], [], [], [], [], [], {0: [0, 0], 1: [0, 0], 2: [0, 0]}
        cs_grad = {0: [0, 0], 1: [0, 0], 2: [0, 0]}
        for G, Adj, z, sn in world_runs:
            I1 = sn[t1][0] == 2
            prev.append(I1.mean()); ors.append(edge_odds_ratio(G, I1))
            nI_all = np.asarray(Adj @ I1.astype(float))
            for c in (0, 1, 2):
                m = (nI_all >= 2) if c == 2 else (nI_all == c)
                cs_grad[c][0] += int(I1[m].sum()); cs_grad[c][1] += int(m.sum())
            sel, nI, I2 = survey_rows(Adj, sn, t1, t2)
            y = I2[sel].astype(float); x = nI[sel]
            b, se = _fit(y, x); est.append(b); sig.append(abs(b / se) > 1.96)
            Y.append(y); NI.append(x); Z.append(z[sel])
            for c in (0, 1, 2):
                m = (x >= 2) if c == 2 else (x == c)
                grad[c][0] += int(y[m].sum()); grad[c][1] += int(m.sum())
        y = np.concatenate(Y); x = np.concatenate(NI); zz = np.concatenate(Z)
        b, se = _fit(y, x); ba, sea = _fit(y, x, zz)
        return dict(prev=float(np.mean(prev)), edge_or=float(np.mean(ors)), b=b, se=se, b_adj=ba, se_adj=sea,
                    sig=float(np.mean(sig)), est_sd=float(np.std(est, ddof=1)), inc=float(y.mean()), n_obs=len(y),
                    grad={c: g[0] / g[1] for c, g in grad.items()}, cs={c: g[0] / g[1] for c, g in cs_grad.items()},
                    cs_share={c: g[1] for c, g in cs_grad.items()})

    def store(prefix, d):
        put(f"{prefix}_prev_pct", 100 * d["prev"]); put(f"{prefix}_edge_or", d["edge_or"])
        put(f"{prefix}_or", math.exp(d["b"])); put(f"{prefix}_or_lo", math.exp(d["b"] - 1.96 * d["se"]))
        put(f"{prefix}_or_hi", math.exp(d["b"] + 1.96 * d["se"]))
        put(f"{prefix}_or_adj", math.exp(d["b_adj"])); put(f"{prefix}_or_adj_lo", math.exp(d["b_adj"] - 1.96 * d["se_adj"]))
        put(f"{prefix}_or_adj_hi", math.exp(d["b_adj"] + 1.96 * d["se_adj"]))
        put(f"{prefix}_sig_pct", 100 * d["sig"]); put(f"{prefix}_inc_pct", 100 * d["inc"]); put(f"{prefix}_n_obs", d["n_obs"])
        for c in (0, 1, 2):
            put(f"{prefix}_inc_f{c}_pct", 100 * d["grad"][c]); put(f"{prefix}_cs_f{c}_pct", 100 * d["cs"][c])

    # contagion world: calibrated central parameters on the same graphs; the trait has no effect
    runsC = []
    for r, (G, Adj, z) in enumerate(zip(graphs, adjs, zs)):
        sn = nx2.simulate_ext(G, pc1, beta, t2, equilibrium(pc1), seed=S(22_000 + r), snapshots=(t1, t2), Adj=Adj)["snaps"]
        runsC.append((G, Adj, z, sn))
    dC = collect(runsC); store("netx_hom_C", dC)
    target = dC["prev"]
    put("netx_hom_target_prev_pct", 100 * target)

    # homophily worlds: b = 0; alpha_i = a0 exp(s z_i), sigma_i = sigma exp(s z_i); a0 matches the prevalence
    res = {}
    for s_ in s_grid:
        runsH, a0s, qratio, amed = [], [], [], []
        for r, (G, Adj, z) in enumerate(zip(graphs, adjs, zs)):
            sv = pc1.sigma * np.exp(s_ * z)
            a0 = brentq(lambda a: nx2.stationary_b0(pc1, a * np.exp(s_ * z), sv)[:, 2].mean() - target, 1e-6, 1e4)
            av = a0 * np.exp(s_ * z)
            P = nx2.stationary_discrete_b0(pc1, av, sv)
            rng = np.random.default_rng(S(23_000 + r))
            init = (rng.random(n)[:, None] > P.cumsum(1)).sum(1)
            sn = nx2.simulate_ext(G, pc1.with_(b=0.0), 0.0, t2, seed=S(24_000 + r), init_states=init, alpha_vec=av,
                                  sigma_vec=sv, snapshots=(t1, t2), Adj=Adj)["snaps"]
            runsH.append((G, Adj, z, sn)); a0s.append(a0); amed.append(np.median(av))
            q = np.quantile(z, [0.2, 0.8])
            qratio.append(P[z >= q[1], 2].mean() / P[z <= q[0], 2].mean())
        d = collect(runsH); res[s_] = d
        pre = f"netx_hom_H_s{tag(s_)}"
        store(pre, d)
        put(f"{pre}_alpha_median", float(np.mean(amed))); put(f"{pre}_quintile_ratio", float(np.mean(qratio)))
    put("netx_hom_C_b", pc1.b); put("netx_hom_C_alpha", pc1.alpha)

    # figure: naive and trait-adjusted odds ratio per affected friend against trait strength
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    xs = np.array(s_grid)
    def ci(d, adj=False):
        b, se = (d["b_adj"], d["se_adj"]) if adj else (d["b"], d["se"])
        return math.exp(b), math.exp(b - 1.96 * se), math.exp(b + 1.96 * se)
    for adj, col, lab, off in [(False, C2, "Homophily, no contagion: naive estimate", -0.01),
                               (True, GREY, "Homophily, no contagion: adjusted for the latent trait", 0.01)]:
        v = np.array([ci(res[s_], adj) for s_ in s_grid])
        ax.errorbar(xs + off, v[:, 0], yerr=[v[:, 0] - v[:, 1], v[:, 2] - v[:, 0]], fmt="o-", color=col, ms=4,
                    lw=1.4, capsize=2, label=lab)
    c = ci(dC)
    ax.axhspan(c[1], c[2], color=C1, alpha=0.15, lw=0)
    ax.axhline(c[0], color=C1, lw=1.4, label="Contagion world (calibrated $b$): naive estimate")
    ax.axhline(1.0, color=INK, lw=0.8, ls=":")
    ax.set_xlabel("Strength of the latent trait, $s$")
    ax.set_ylabel("Odds ratio per affected friend")
    ax.legend(frameon=False, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1)
    fig.savefig(FIG / "fig_netx_homophily.png"); plt.close(fig)
    return res, dC


# ---------------------------------------------------------------- 3. complex contagion
def run_complex(pc1, x22, quick):
    n = 2000
    reps = 4 if quick else 20
    mults = [1.0, 2.0, 4.0]
    put("netx_cc_reps", reps); put("netx_cc_min_affected", 2)
    graphs = {kind: [nw.make_graph(kind, n, K, seed=S(30_000 + 100 * j + s)) for s in range(reps)]
              for j, kind in enumerate(["ER", "WS"])}
    adjs = {kind: [nx.to_scipy_sparse_array(G, format="csr", dtype=float) for G in gs] for kind, gs in graphs.items()}
    clus = {kind: float(np.mean([nx.average_clustering(G) for G in gs[:5]])) for kind, gs in graphs.items()}
    for kind in clus:
        put(f"netx_cc_{kind}_clustering", clus[kind])
    # (a) invasion from a local seed (a node and its neighbours Addicted), alpha = 0, 40 years
    years = 40
    put("netx_cc_inv_years", years)
    for mult in mults:
        for kind in ["ER", "WS"]:
            for rule, ma in [("simple", None), ("complex", 2)]:
                est, I10, I40 = [], [], []
                for s, (G, Adj) in enumerate(zip(graphs[kind], adjs[kind])):
                    q = pc1.with_(b=pc1.b * mult, alpha=0.0)
                    init = np.zeros(n, int); init[0] = 2; init[list(G.neighbors(0))] = 2
                    tr = nx2.simulate_ext(G, q, nw.beta_from_b(q, K), years * 12, seed=S(31_000 + s), init_states=init,
                                          min_affected=ma, Adj=Adj)["traj"]
                    est.append(tr[-1, 1] + tr[-1, 2] + tr[-1, 3] >= 0.05); I10.append(tr[10, 2]); I40.append(tr[-1, 2])
                tg = f"netx_cc_inv_{kind}_{rule}_m{tag(mult)}"
                put(f"{tg}_est_pct", 100 * float(np.mean(est)))
                put(f"{tg}_I10_pct", 100 * float(np.mean(I10))); put(f"{tg}_I40_pct", 100 * float(np.mean(I40)))
    # (b) alpha at its calibrated value: long-run level from a problem-free start and from the 2022 state
    years2 = 30 if quick else 60
    hreps = 3 if quick else 8
    put("netx_cc_hys_years", years2); put("netx_cc_hys_reps", hreps)
    hm = [0.3, 0.6, 1.0] if quick else [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.25]
    put("netx_cc_hys_mult_min", min(hm)); put("netx_cc_hys_mult_max", max(hm))
    hys = {}
    for kind in ["ER", "WS"]:
        for rule, ma, starts in [("simple", None, ["low"]), ("complex", 2, ["low", "high"])]:
            for start in starts:
                vals = []
                for mult in hm:
                    v = []
                    for s in range(hreps):
                        G, Adj = graphs[kind][s], adjs[kind][s]
                        q = pc1.with_(b=pc1.b * mult)
                        kw = dict(init_states=np.zeros(n, int)) if start == "low" else dict(x0=x22)
                        tr = nx2.simulate_ext(G, q, nw.beta_from_b(q, K), years2 * 12, seed=S(32_000 + s), min_affected=ma,
                                              Adj=Adj, **kw)["traj"]
                        v.append(tr[-6:, 2].mean())
                    vals.append(float(np.mean(v)))
                    put(f"netx_cc_hys_{kind}_{rule}_{start}_m{tag(mult)}_I_pct", 100 * vals[-1])
                hys[(kind, rule, start)] = vals
    # width of the bistable window: multipliers at which high and low starts differ by more than 2 points
    for kind in ["ER", "WS"]:
        lo, hi = np.array(hys[(kind, "complex", "low")]), np.array(hys[(kind, "complex", "high")])
        gap = 100 * (hi - lo)
        bis = [m for m, g_ in zip(hm, gap) if g_ > 2.0]
        put(f"netx_cc_hys_{kind}_bistable_min", min(bis) if bis else float("nan"))
        put(f"netx_cc_hys_{kind}_bistable_max", max(bis) if bis else float("nan"))
        put(f"netx_cc_hys_{kind}_maxgap_pts", float(gap.max()))
    # is the start dependence lasting? rerun the most start-dependent multiplier for 150 years
    long_years = 60 if quick else 150
    put("netx_cc_long_years", long_years)
    for kind in ["ER", "WS"]:
        lo, hi = np.array(hys[(kind, "complex", "low")]), np.array(hys[(kind, "complex", "high")])
        m_ = hm[int(np.argmax(hi - lo))]
        put(f"netx_cc_long_{kind}_mult", m_)
        for start in ["low", "high"]:
            trs = []
            for s in range(hreps):
                G, Adj = graphs[kind][s], adjs[kind][s]
                q = pc1.with_(b=pc1.b * m_)
                kw = dict(init_states=np.zeros(n, int)) if start == "low" else dict(x0=x22)
                trs.append(nx2.simulate_ext(G, q, nw.beta_from_b(q, K), long_years * 12, seed=S(33_000 + s), min_affected=2,
                                            Adj=Adj, **kw)["traj"][:, 2])
            trs = 100 * np.array(trs).mean(0)
            put(f"netx_cc_long_{kind}_{start}_y60", float(trs[55:61].mean()))
            put(f"netx_cc_long_{kind}_{start}_yend", float(trs[-6:].mean()))
            put(f"netx_cc_long_{kind}_{start}_y10", float(trs[10]))
    ode = [100 * equilibrium(pc1.with_(b=pc1.b * m))[2] for m in hm]
    for m, v in zip(hm, ode):
        put(f"netx_cc_hys_ode_m{tag(m)}_I_pct", v)
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3), sharey=True)
    for ax, kind in zip(axs, ["ER", "WS"]):
        ax.plot(hm, ode, color=INK, ls="--", lw=1.1, label="Mean-field ODE (simple)")
        ax.plot(hm, 100 * np.array(hys[(kind, "simple", "low")]), "o-", color=GREY, ms=3, lw=1.3, label="Simple contagion")
        ax.plot(hm, 100 * np.array(hys[(kind, "complex", "low")]), "o-", color=COLS[kind], ms=3, lw=1.5,
                label="Complex, problem-free start")
        ax.plot(hm, 100 * np.array(hys[(kind, "complex", "high")]), "s--", color=COLS[kind], ms=3, lw=1.5, mfc="white",
                label="Complex, start at 2022 state")
        ax.axvline(1.0, color=GREY, lw=0.8, ls=":")
        ax.set_title(f"{'A' if kind == 'ER' else 'B'}. {NAMES[kind]}", fontsize=9, loc="left")
        ax.set_xlabel("Exposure rate relative to calibrated $b_1$")
        ax.set_ylim(0, None)
    axs[0].set_ylabel(f"Problematic use after {years2} years (%)")
    axs[1].legend(frameon=False, fontsize=7, loc="lower right")
    fig.savefig(FIG / "fig_netx_complex.png"); plt.close(fig)
    return hys


# ---------------------------------------------------------------- 4. targeting
def run_targeting(pc1, x22, quick):
    n = 2000
    reps = 6 if quick else 40
    covers = [0.05, 0.10, 0.20, 0.30]
    damp = 0.2
    put("netx_tgt_reps", reps); put("netx_tgt_damp_pct", 100 * damp)
    G0 = nw.make_graph("BA", n, K, seed=S(40_000))
    dists = [("r1p5", 1.5), ("r2p5", 2.5), ("r4", 4.0), ("cal", None)]
    settings = []
    for lab, tgt in dists:
        if tgt is None:
            p = pc1; x0 = x22
        else:
            b = ra._b_for_Rnet(pc1, G0, K, tgt); p = pc1.with_(b=b); x0 = equilibrium(p)
        put(f"netx_tgt_{lab}_b_ratio", p.b / pc1.b)
        put(f"netx_tgt_{lab}_Rnet", nw.R0_network(p.with_(alpha=0.0), nw.beta_from_b(p, K), G0))
        settings.append((lab, p, x0))
    strategies = ["random", "hub", "acq"]
    res = {}
    acq_deg = []
    for lab, p, x0 in settings:
        beta = nw.beta_from_b(p, K)
        base = []; red = {(st, c): [] for st in strategies for c in covers}
        for s in range(reps):
            G = nw.make_graph("BA", n, K, seed=S(41_000 + s)); Adj = nx.to_scipy_sparse_array(G, format="csr", dtype=float)
            deg = np.array([d for _, d in G.degree()])
            seed = S(42_000 + s)
            a = nx2.simulate_ext(G, p, beta, 240, x0, seed=seed, Adj=Adj, record_every=240)["traj"][-1, 2]
            base.append(a)
            for c in covers:
                m = int(round(c * n))
                rng = np.random.default_rng(S(43_000 + 100 * s + int(100 * c)))
                sel = {"random": rng.choice(n, m, replace=False), "hub": np.argsort(-deg, kind="stable")[:m],
                       "acq": nx2.acquaintance_nodes(G, m, rng)}
                if lab == "cal" and c == 0.10:
                    acq_deg.append(deg[sel["acq"]].mean())
                for st in strategies:
                    w = np.ones(n); w[sel[st]] = damp
                    v = nx2.simulate_ext(G, p, beta, 240, x0, seed=seed, Adj=Adj, w=w, record_every=240)["traj"][-1, 2]
                    red[(st, c)].append(a - v)
        base = 100 * np.array(base)
        put(f"netx_tgt_{lab}_base", float(base.mean()))
        for (st, c), v in red.items():
            v = 100 * np.array(v)
            tg = f"netx_tgt_{lab}_{st}_c{int(round(100 * c))}"
            put(f"{tg}_red", float(v.mean())); put(f"{tg}_se", float(v.std(ddof=1) / math.sqrt(reps)))
            put(f"{tg}_rel_pct", float(100 * v.mean() / base.mean()))
            res[(lab, st, c)] = (float(100 * v.mean() / base.mean()), float(100 * v.std(ddof=1) / math.sqrt(reps) / base.mean()))
        for c in covers:
            dh = 100 * (np.array(red[("hub", c)]) - np.array(red[("random", c)]))
            da = 100 * (np.array(red[("acq", c)]) - np.array(red[("random", c)]))
            tg = f"netx_tgt_{lab}_c{int(round(100 * c))}"
            put(f"{tg}_hub_minus_random", float(dh.mean())); put(f"{tg}_hub_minus_random_se", float(dh.std(ddof=1) / math.sqrt(reps)))
            put(f"{tg}_acq_minus_random", float(da.mean())); put(f"{tg}_acq_minus_random_se", float(da.std(ddof=1) / math.sqrt(reps)))
    put("netx_tgt_acq_mean_degree", float(np.mean(acq_deg)))
    G = nw.make_graph("BA", n, K, seed=S(41_000))
    deg = np.array([d for _, d in G.degree()])
    put("netx_tgt_hub_mean_degree_c10", float(np.sort(deg)[::-1][:200].mean()))
    put("netx_tgt_hub_min_degree_c10", int(np.sort(deg)[::-1][199]))

    fig, axs = plt.subplots(1, 4, figsize=(8.6, 3.0), sharey=True)
    titles = {"r1p5": r"A. $\mathcal{R}_{\mathrm{net}} = 1.5$", "r2p5": r"B. $\mathcal{R}_{\mathrm{net}} = 2.5$",
              "r4": r"C. $\mathcal{R}_{\mathrm{net}} = 4$", "cal": "D. Calibrated $b_1$"}
    for ax, (lab, _, _) in zip(axs, settings):
        for st, col, nm in [("random", GREY, "Random"), ("acq", C3, "Acquaintance"), ("hub", C2, "Highest degree")]:
            v = np.array([res[(lab, st, c)] for c in covers])
            ax.errorbar(100 * np.array(covers), v[:, 0], yerr=1.96 * v[:, 1], fmt="o-", color=col, ms=3.5, lw=1.4,
                        capsize=2, label=nm)
        ax.set_title(titles[lab], fontsize=9, loc="left"); ax.set_xlabel("Coverage (%)")
        ax.set_xticks([5, 10, 20, 30]); ax.set_ylim(0, None)
    axs[0].set_ylabel("Reduction in problematic use\nafter 20 years (%)")
    axs[3].legend(frameon=False, fontsize=7.5, loc="upper left")
    fig.savefig(FIG / "fig_netx_targeting.png"); plt.close(fig)
    return res


# ---------------------------------------------------------------- 5. mean degree and the exposure deficit
def deficit(Adj, deg, st, eta):
    """Ratio of the social exposure of Susceptible nodes on the network to its mean-field counterpart, split
    into a degree factor and a neighbourhood factor:
      rho = mean_S(m_I + eta m_A) / (kbar (I + eta A))
          = [mean_S k / kbar] x [mean_S(m_I + eta m_A) / mean_S k / (I + eta A)]."""
    S_ = st == 0
    f = (st == 2) + eta * (st == 1)
    m = np.asarray(Adj @ f.astype(float))
    kbar = deg.mean(); pop = f.mean()
    degf = deg[S_].mean() / kbar
    locf = m[S_].mean() / deg[S_].mean() / pop
    return degf * locf, degf, locf


def run_degree(pc1, x22, quick):
    n = 2000
    reps = 3 if quick else 24
    ks = [4, 10, 20] if quick else [4, 6, 8, 10, 14, 20, 30, 40]
    months = 240
    put("netx_deg_reps", reps); put("netx_deg_kmin", min(ks)); put("netx_deg_kmax", max(ks))
    _, y = simulate(pc1, x22, 20, t_eval=np.arange(0, 21))
    put("netx_deg_meanfield_I2042", 100 * y[2, -1])
    out = {}
    for j, kind in enumerate(["ER", "WS", "BA"]):
        vals = []
        for kk in ks:
            beta = nw.beta_from_b(pc1, kk)
            v, rho, dg, lc = [], [], [], []
            for s in range(reps):
                G = nw.make_graph(kind, n, kk, seed=S(50_000 + 1000 * j + 10 * kk + s))
                Adj = nx.to_scipy_sparse_array(G, format="csr", dtype=float)
                deg = np.array([d for _, d in G.degree()], float)
                snaps = tuple(range(months - 48, months + 1, 12))
                r = nx2.simulate_ext(G, pc1, beta, months, x22, seed=S(51_000 + 1000 * j + 10 * kk + s), Adj=Adj,
                                     record_every=months, snapshots=snaps)
                v.append(r["traj"][-1, 2])
                dd = np.array([deficit(Adj, deg, r["snaps"][t][0], pc1.eta) for t in snaps])
                rho.append(dd[:, 0].mean()); dg.append(dd[:, 1].mean()); lc.append(dd[:, 2].mean())
            vals.append(100 * float(np.mean(v)))
            put(f"netx_deg_{kind}_k{kk}_I2042", vals[-1])
            put(f"netx_deg_{kind}_k{kk}_I2042_se", 100 * float(np.std(v, ddof=1) / math.sqrt(reps)))
            put(f"netx_deg_{kind}_k{kk}_rho", float(np.mean(rho)))
            put(f"netx_deg_{kind}_k{kk}_degfactor", float(np.mean(dg)))
            put(f"netx_deg_{kind}_k{kk}_locfactor", float(np.mean(lc)))
            # first-order check: mean-field ODE with the exposure rate scaled by the measured deficit
            _, yr = simulate(pc1.with_(b=pc1.b * float(np.mean(rho))), x22, 20, t_eval=np.arange(0, 21))
            put(f"netx_deg_{kind}_k{kk}_ode_rho_I2042", 100 * yr[2, -1])
            out[(kind, kk)] = (vals[-1], float(np.mean(rho)), 100 * float(np.std(v, ddof=1) / math.sqrt(reps)))
    fig, axs = plt.subplots(1, 2, figsize=(7.6, 3.2), gridspec_kw={"wspace": 0.4})
    mf = 100 * y[2, -1]
    for j, kind in enumerate(["ER", "WS", "BA"]):
        off = np.exp((j - 1) * 0.03)
        axs[0].errorbar(np.array(ks) * off, [out[(kind, kk)][0] - mf for kk in ks], yerr=[1.96 * out[(kind, kk)][2] for kk in ks],
                        fmt="o-", color=COLS[kind], ms=3.5, lw=1.5, capsize=2, label=NAMES[kind])
        axs[1].plot(ks, [out[(kind, kk)][1] for kk in ks], "o-", color=COLS[kind], ms=3.5, lw=1.5, label=NAMES[kind])
    axs[0].axhline(0.0, color=INK, ls="--", lw=1.1, label="Mean-field ODE")
    axs[1].axhline(1.0, color=INK, ls="--", lw=1.1)
    for ax in axs:
        ax.set_xscale("log"); ax.set_xticks(ks); ax.set_xticklabels([str(k) for k in ks]); ax.minorticks_off()
        ax.set_xlabel("Mean degree")
    axs[0].set_ylabel("Difference from mean-field level\nin 2042 (percentage points)")
    axs[1].set_ylabel("Exposure of Susceptible nodes\nrelative to mean field")
    axs[0].legend(frameon=False, fontsize=7.5, loc="lower right")
    axs[0].set_title("A. Level in 2042 against the ODE", fontsize=9, loc="left"); axs[1].set_title("B. Exposure deficit", fontsize=9, loc="left")
    fig.savefig(FIG / "fig_netx_degree.png"); plt.close(fig)
    return out


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true"); a = ap.parse_args()
    t0 = time.time()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)
    _, y = simulate(pc1, x18, 4, t_eval=[4]); x22 = y[:, -1]
    put("netx_seed", SEED + 100_000)
    for name, fn in [("threshold", lambda: run_threshold(pc1, a.quick)),
                     ("homophily", lambda: run_homophily(pc1, a.quick)),
                     ("complex", lambda: run_complex(pc1, x22, a.quick)),
                     ("targeting", lambda: run_targeting(pc1, x22, a.quick)),
                     ("degree", lambda: run_degree(pc1, x22, a.quick))]:
        t1 = time.time(); fn()
        print(f"{name}: {(time.time() - t1) / 60:.1f} min", flush=True)
    put("netx_runtime_minutes", (time.time() - t0) / 60)
    name = "quick_netx.json" if a.quick else "results_netx.json"
    json.dump(R, open(OUT / name, "w"), indent=1, ensure_ascii=False)
    print(f"done in {(time.time() - t0) / 60:.1f} min, {len(R)} results -> {name}")


if __name__ == "__main__":
    main()
