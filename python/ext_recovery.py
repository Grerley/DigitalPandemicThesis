"""Extended parameter-recovery study for thesis Chapter 6 (workstream recx).

    python ext_recovery.py            # full run (about 25 minutes on two cores)
    python ext_recovery.py --quick    # smoke test

Writes outputs/results_recx.json (keys prefixed recx_), outputs/recx_estimates.csv (every replicate)
and figures/fig_recx_*.png. Synthetic data with known truth only: nothing here is an estimate for any
population, and nothing enters the calibrated parameter set of Chapter 7.

Analyses
  1. Design grid: waves (3, 6, 10) x panel size (500, 2,000, 8,000) x noise (low, high), at the
     design-study rates and at the calibrated central rates (Gaussian score, as in run_all.run_recovery).
  2. Incidence bias: the original EM stopping rule against strict convergence, each with and without a
     baseline item on past problematic use (perfect, and sensitivity 0.8 / specificity 0.9).
  3. Measurement: a 0 to 9 symptom count with binomial or beta-binomial emissions; misclassification
     implied by the cut-off of 6 or more; recovery from the full count, three bands, or the binary
     classification only; misspecified emission families.
  4. Design table: smallest panel reaching 20% and 10% relative RMSE, by transition, waves and noise.
"""
from __future__ import annotations

import os
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse, json, math, multiprocessing as mp, time

import numpy as np
from scipy.stats import norm

import run_all as ra                    # house colours, rcParams, paths, seed
from dpm import recovery as rec
from dpm import recovery_ext as rx

plt = ra.plt
OUT, FIG = ra.OUT, ra.FIG
SEED = ra.SEED + 700
RX = {}

ENTRIES = [(0, 1), (1, 0), (1, 2), (2, 3), (3, 1)]
TR = ["S_A", "A_S", "A_I", "I_R", "R_A"]
TRL = ["S→A", "A→S", "A→I", "I→R", "R→A"]
COLS = [ra.C1, ra.C2, ra.C3, ra.INK, ra.GREY]
MU = np.array([1.0, 3.5, 7.0, 1.0])
SD = {"low": np.array([0.8, 1.0, 1.0, 0.8]), "high": np.array([1.5, 1.8, 1.8, 1.5])}
STRICT = dict(tol=1e-9, iters=1500, accel=True, rule="abs")      # converged fits (incidence-bias test)
GRIDFIT = dict(tol=1e-9, iters=1500, accel=True, rule="abs")     # design grid and measurement arms
ORIG = dict(tol=1e-7, iters=500, accel=False, rule="rel")        # rule of dpm.recovery.baum_welch
Q_BIN = np.array([0.06, 0.36, 0.75, 0.06])                       # P(item endorsed | S, A, I, R)
RHO = 0.15                                                       # beta-binomial item correlation
H_ITEM = {"perfect": (1.0, 1.0), "imperfect": (0.8, 0.9)}        # (sensitivity, specificity)
H_AI = 0.4                                                       # P(z = 1) for A and I at entry
# (arm, fitting rule, item, item error rates fixed at their stated values)
ITEM_ARMS = [("orig_none", ORIG, None, True), ("orig_perfect", ORIG, "perfect", True),
             ("strict_none", STRICT, None, True), ("strict_perfect", STRICT, "perfect", True),
             ("strict_imperfect", STRICT, "imperfect", True), ("strict_perfectfree", STRICT, "perfect", False)]


def put(key, value, note=""):
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    RX[f"recx_{key}"] = {"value": value, "note": note}


def rate_sets():
    """Design-study rates (as run_all.run_recovery) and calibrated central rates (as stored in results.json
    by run_all.run_recovery_central)."""
    R = json.load(open(OUT / "results.json"))
    design = dict(sigma=0.40, gamma1=0.50, gamma2=0.60, delta=0.15, incidence=0.08)
    central = {k: R[f"recc_truth_{k}"]["value"] for k in design}
    pi_c = np.array([R[f"recc_pi0_{s}"]["value"] for s in "SAIR"]) / 100
    return {"design": (rec.annual_matrix(**design), np.array([0.75, 0.12, 0.08, 0.05])),
            "central": (rec.annual_matrix(**central), pi_c / pi_c.sum())}


RATES = rate_sets()


def est(f):
    return [f["P"][i, j] for i, j in ENTRIES]


# ---------------------------------------------------------------- tasks (one replicate each)
def task_grid(a):
    rates, noise, T, n, seed = a
    P, pi0 = RATES[rates]
    g = np.random.default_rng(seed)
    s = rx.simulate_states(P, pi0, n, T, g)
    y = rx.emit_gauss(s, MU, SD[noise], g)
    f = rx.fit(y, "gauss", g, starts=2, **GRIDFIT)
    return {"grid": est(f) + [f["pi"][3], f["iters"]]}


def task_item(a):
    noise, n, seed = a
    P, pi0 = RATES["design"]
    g = np.random.default_rng(seed)
    s = rx.simulate_states(P, pi0, n, 6, g)
    y = rx.emit_gauss(s, MU, SD[noise], g)
    z = {}
    for k, (se, sp) in H_ITEM.items():
        h = np.array([1 - sp, H_AI, H_AI, se])
        z[k] = rx.baseline_item(s[:, 0], h, g)
    out = {}
    sd0 = int(g.integers(0, 2**31))             # the transition and emission starts are shared by all arms
    for arm, kw, item, fixed in ITEM_ARMS:
        if arm == "strict_perfectfree" and noise == "high":
            continue                            # costly near-flat fits; run at low noise only
        hf = None
        if item is not None and fixed:
            se, sp = H_ITEM[item]; hf = {0: 1 - sp, 3: se}
        f = rx.fit(y, "gauss", np.random.default_rng(sd0), starts=2, z=None if item is None else z[item],
                   h_fix=hf, **kw)
        out[arm] = est(f) + [f["pi"][3], f["iters"], f["loglik"]]
    return out


def hazards(Pm):
    """Invert dpm.recovery.annual_matrix: yearly hazards from a fitted annual transition matrix."""
    lam = -np.log(1 - Pm[0, 1]); g2 = -np.log(1 - Pm[2, 3]); de = -np.log(1 - Pm[3, 1])
    leave = Pm[1, 2] + Pm[1, 0]; hA = -np.log(1 - leave)
    return [lam, hA * Pm[1, 0] / leave, hA * Pm[1, 2] / leave, g2, de]     # order of TR


HZ = ["incidence", "gamma1", "sigma", "gamma2", "delta"]                       # order of TR


def truth_hazards(rates):
    R = json.load(open(OUT / "results.json"))
    if rates == "design":
        return dict(sigma=0.40, gamma1=0.50, gamma2=0.60, delta=0.15, incidence=0.08)
    return {k: R[f"recc_truth_{k}"]["value"] for k in HZ}


def exact_matrix(t):
    """Exact one-year transition matrix of the continuous-time chain (no turnover within the panel)."""
    from scipy.linalg import expm
    Q = np.zeros((4, 4))
    Q[0, 1] = t["incidence"]; Q[1, 0] = t["gamma1"]; Q[1, 2] = t["sigma"]; Q[2, 3] = t["gamma2"]; Q[3, 1] = t["delta"]
    np.fill_diagonal(Q, -Q.sum(1))
    return expm(Q)


def task_exact(a):
    rates, gen, n, seed = a
    _, pi0 = RATES[rates]
    t = truth_hazards(rates)
    P = exact_matrix(t) if gen == "exact" else rec.annual_matrix(**t)
    g = np.random.default_rng(seed)
    s = rx.simulate_states(P, pi0, n, 6, g)
    y = rx.emit_gauss(s, MU, SD["low"], g)
    f = rx.fit(y, "gauss", g, starts=2, **GRIDFIT)
    out = {"exact": hazards(f["P"]) + [f["pi"][3], f["iters"]]}
    if gen == "exact":                          # the remedy: fit the continuous-time model itself
        c = rx.fit_ct(y, f)
        out["ct"] = list(c["hz"]) + [c["pi"][3], float(c["converged"])]
    return out


def run_exact(q, t0):
    """Annual waves of a continuous-time process: fit the one-transition-per-year model to panels generated
    from the exact one-year matrix exp(Q), and from the one-transition matrix as a control."""
    reps = 2 if q else 10
    n = 2000
    jobs = []
    seeds = np.random.SeedSequence(SEED + 2).spawn(4 * reps)
    i = 0
    for rates in RATES:
        for gen in ["onestep", "exact"]:
            for r in range(reps):
                jobs.append(("exact", i, (rates, gen, n, int(seeds[i].generate_state(1)[0])))); i += 1
    with mp.get_context("fork").Pool(2) as pool:
        out = list(pool.imap(run_task, jobs))
    put("exact_reps", reps); put("exact_n", n)
    for rates in RATES:
        t = truth_hazards(rates)
        E1 = exact_matrix(t)
        forb = np.where(rec.ALLOWED, 0.0, E1).sum(1)
        for k, nm in enumerate("SAIR"):
            put(f"exact_{rates}_forbidden_{nm}_pct", 100 * float(forb[k]),
                "one-year probability of a transition the fitted model sets to zero")
        put(f"exact_{rates}_AR_pct", 100 * float(E1[1, 3])); put(f"exact_{rates}_SI_pct", 100 * float(E1[0, 2]))
        true = np.array([t[h] for h in HZ])
        Ec = np.array([o[2]["ct"] for o, j in zip(out, jobs) if j[2][0] == rates and j[2][1] == "exact"])
        put(f"exact_{rates}_ct_converged_pct", 100 * float(Ec[:, 6].mean()))
        for gen in ["onestep", "exact", "ct"]:
            if gen == "ct":
                E = Ec
            else:
                E = np.array([o[2]["exact"] for o, j in zip(out, jobs) if j[2][0] == rates and j[2][1] == gen])
            b, sb, m, sm = summarise(E[:, :5], true)
            for kk, tr in enumerate(TR):
                pre = f"exact_{rates}_{gen}_{tr}"
                put(f"{pre}_hz_relbias_pct", float(b[kk])); put(f"{pre}_hz_relbias_mcse_pct", float(sb[kk]))
                put(f"{pre}_hz_relrmse_pct", float(m[kk]))
    return sum(o[3] for o in out)


def bands(y):
    return np.where(y >= 6, 2, np.where(y >= 2, 1, 0))


def task_meas(a):
    n, seed = a
    P, pi0 = RATES["design"]
    g = np.random.default_rng(seed)
    s = rx.simulate_states(P, pi0, n, 6, g)
    yb = rx.emit_binom(s, Q_BIN, g)
    ybb = rx.emit_betabinom(s, Q_BIN, RHO, g)
    out = {}
    arms = [("bin_binom", yb, "binom", 10), ("bin_cat", yb, "cat", 10), ("bin_band3", bands(yb), "cat", 3),
            ("bin_bin6", (yb >= 6).astype(int), "cat", 2),
            ("bb_binom", ybb, "binom", 10), ("bb_cat", ybb, "cat", 10), ("bb_gauss", ybb, "gauss", 10)]
    for name, y, fam, M in arms:
        U, _, w = rx.compress(y)
        f = rx.fit(U.astype(float) if fam == "gauss" else U, fam, g, starts=2, wts=w, M=M, **GRIDFIT)
        out[name] = est(f) + [f["pi"][3], f["iters"]]
    return out


def run_task(job):
    kind, idx, args = job
    t = time.time()
    r = {"grid": task_grid, "item": task_item, "meas": task_meas, "exact": task_exact}[kind](args)
    return kind, idx, r, time.time() - t


# ---------------------------------------------------------------- summaries
def summarise(E, true):
    """Relative bias and RMSE (per cent of truth) with Monte Carlo standard errors."""
    e = np.asarray(E) - true
    R = len(e)
    bias = e.mean(0); mse = (e ** 2).mean(0); rmse = np.sqrt(mse)
    se_bias = e.std(0, ddof=1) / np.sqrt(R)
    se_rmse = (e ** 2).std(0, ddof=1) / np.sqrt(R) / (2 * rmse)
    return 100 * bias / true, 100 * se_bias / true, 100 * rmse / true, 100 * se_rmse / true


def needed_n(sizes, rr, target):
    """Smallest panel with relative RMSE <= target, by log-log interpolation between grid sizes."""
    if rr[0] <= target:
        return sizes[0], "le"
    for k in range(1, len(sizes)):
        if rr[k] <= target:
            x0, x1 = math.log(sizes[k - 1]), math.log(sizes[k])
            y0, y1 = math.log(rr[k - 1]), math.log(rr[k])
            x = x0 + (math.log(target) - y0) * (x1 - x0) / (y1 - y0)
            return math.exp(x), "in"
    return None, "gt"


def ceil2(x):
    """Round up to two significant figures."""
    m = 10 ** (math.floor(math.log10(x)) - 1)
    return int(math.ceil(x / m) * m)


def design_table(grid, WAVES, SIZES):
    """Smallest panel reaching 20% and 10% relative RMSE, by log-log interpolation between grid sizes.
    A dagger marks cells where the relative bias at the largest size exceeds half the target."""
    for target in [20, 10]:
        for rates in RATES:
            for noise in SD:
                for T in WAVES:
                    for k, tr in enumerate(TR):
                        rr = [grid[(rates, noise, T, n)][2][k] for n in SIZES]
                        bias_last = abs(grid[(rates, noise, T, SIZES[-1])][0][k])
                        nn, how = needed_n(SIZES, rr, target)
                        if how == "le":
                            txt = f"≤ {SIZES[0]:,}"; val = SIZES[0]
                        elif how == "in":
                            val = ceil2(nn); txt = f"{val:,}"
                        else:
                            val = None
                            txt = f"> {SIZES[-1]:,}" + ("†" if bias_last > target / 2 else "")
                        put(f"need{target}_{rates}_{noise}_T{T}_{tr}", txt,
                            f"smallest panel for relative RMSE <= {target}%; numeric {val}")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--quick", action="store_true")
    ap.add_argument("--exact-only", action="store_true", help="run only the exact-process check and merge it")
    ap.add_argument("--figures-only", action="store_true", help="redraw figures from outputs/recx_estimates.csv")
    a = ap.parse_args()
    if a.figures_only:                          # redraw figures and rebuild the design table from saved estimates
        RX.update(json.load(open(OUT / "results_recx.json")))
        grid, _, _ = from_csv()
        design_table(grid, sorted({k[2] for k in grid}), sorted({k[3] for k in grid}))
        figures_only()
        json.dump(RX, open(OUT / "results_recx.json", "w"), indent=1, ensure_ascii=False); return
    t0 = time.time()
    q = a.quick
    if a.exact_only:
        RX.update(json.load(open(OUT / "results_recx.json")))
        cpu = run_exact(q, t0)
        put("cpu_minutes_exact", cpu / 60); put("runtime_minutes_exact", (time.time() - t0) / 60)
        json.dump(RX, open(OUT / "results_recx.json", "w"), indent=1, ensure_ascii=False)
        print(f"exact check done in {(time.time() - t0) / 60:.1f} min"); return
    WAVES = [3, 6, 10]; SIZES = [500, 2000, 8000]
    REPS = {500: 3, 2000: 2, 8000: 1} if q else {500: 32, 2000: 16, 8000: 8}
    ITEM_REPS = {"low": 2, "high": 1} if q else {"low": 8, "high": 6}     # panels of 8,000, six waves
    SIZES_MEAS = {2000: 2, 8000: 1} if q else {2000: 12, 8000: 5}
    N_ITEM = 2000 if q else 8000
    if q:
        SIZES = [500, 2000]; SIZES_MEAS = {2000: 1}

    jobs = []
    for rates in RATES:
        for noise in SD:
            for T in WAVES:
                for n in SIZES:
                    for r in range(REPS[n]):
                        jobs.append(("grid", (rates, noise, T, n), [rates, noise, T, n]))
    for noise, reps in ITEM_REPS.items():
        for r in range(reps):
            jobs.append(("item", (noise, N_ITEM), [noise, N_ITEM]))
    for n, reps in SIZES_MEAS.items():
        for r in range(reps):
            jobs.append(("meas", (n,), [n]))
    seeds = np.random.SeedSequence(SEED).spawn(len(jobs))
    full = []
    for i, (kind, cell, args) in enumerate(jobs):
        full.append((kind, i, tuple(args) + (int(seeds[i].generate_state(1)[0]),)))
    # largest jobs first for load balance; results are keyed by job index, so order does not matter
    cost = {"grid": 1.0, "item": 4.0, "meas": 3.0}
    def jc(j):
        args = j[2]
        if j[0] == "grid":
            return args[3] * args[2]
        return cost[j[0]] * args[-2] * 6
    order = sorted(full, key=jc, reverse=True)
    res = {}; tsum = {"grid": 0.0, "item": 0.0, "meas": 0.0}
    with mp.get_context("fork").Pool(2) as pool:
        for k, (kind, idx, r, dt) in enumerate(pool.imap_unordered(run_task, order)):
            res[idx] = r; tsum[kind] += dt
            if k % 25 == 0:
                print(f"{k}/{len(order)} jobs, {(time.time() - t0) / 60:.1f} min", flush=True)

    # collect per cell
    cells = {}
    for i, (kind, cell, _) in enumerate(jobs):
        cells.setdefault((kind, cell), []).append(res[i])
    rows = []

    # ---------------- 1. grid
    put("grid_waves", WAVES); put("grid_sizes", SIZES)
    for n in SIZES:
        put(f"grid_reps_{n}", REPS[n], "replicates per cell")
    put("noise_low_sd_A", float(SD["low"][1])); put("noise_high_sd_A", float(SD["high"][1]))
    put("noise_low_sd_S", float(SD["low"][0])); put("noise_high_sd_S", float(SD["high"][0]))
    put("strict_tol", STRICT["tol"]); put("strict_iters", STRICT["iters"])
    put("gridfit_tol", GRIDFIT["tol"]); put("gridfit_iters", GRIDFIT["iters"])
    put("orig_tol", ORIG["tol"]); put("orig_iters", ORIG["iters"])
    grid = {}
    for rates in RATES:
        P, _ = RATES[rates]
        true = np.array([P[i, j] for i, j in ENTRIES])
        for noise in SD:
            for T in WAVES:
                for n in SIZES:
                    E = np.array([r["grid"] for r in cells[("grid", (rates, noise, T, n))]])
                    b, sb, m, sm = summarise(E[:, :5], true)
                    grid[(rates, noise, T, n)] = (b, sb, m, sm)
                    for k, tr in enumerate(TR):
                        pre = f"g_{rates}_{noise}_T{T}_n{n}_{tr}"
                        put(f"{pre}_relbias_pct", float(b[k])); put(f"{pre}_relbias_mcse_pct", float(sb[k]))
                        put(f"{pre}_relrmse_pct", float(m[k])); put(f"{pre}_relrmse_mcse_pct", float(sm[k]))
                    put(f"g_{rates}_{noise}_T{T}_n{n}_piR_mean_pct", 100 * float(E[:, 5].mean()))
                    put(f"g_{rates}_{noise}_T{T}_n{n}_capped_pct", 100 * float((E[:, 6] >= GRIDFIT["iters"]).mean()),
                        "share of fits whose best start reached the iteration cap")
                    for e in E:
                        rows.append(["grid", rates, noise, T, n] + list(e))
        for k, tr in enumerate(TR):
            put(f"true_{rates}_{tr}", float(true[k]), "annual transition probability")
    # largest Monte Carlo standard error relative to its RMSE, over the grid
    rel_mcse = [grid[c][3][k] / grid[c][2][k] for c in grid for k in range(5)]
    put("grid_mcse_rel_median_pct", 100 * float(np.median(rel_mcse)))
    put("grid_mcse_rel_max_pct", 100 * float(np.max(rel_mcse)))

    # ---------------- 4. design table (from the grid)
    design_table(grid, WAVES, SIZES)

    # ---------------- 2. baseline item
    P, pi0 = RATES["design"]
    true = np.array([P[i, j] for i, j in ENTRIES])
    put("item_h_AI", H_AI)
    for k, (se, sp) in H_ITEM.items():
        put(f"item_{k}_sens", se); put(f"item_{k}_spec", sp)
    arms = [a[0] for a in ITEM_ARMS]
    for noise in ["low", "high"]:
        for n in [N_ITEM]:
            key = ("item", (noise, n))
            if key not in cells:
                continue
            reps = cells[key]
            put(f"item_{noise}_{n}_reps", len(reps))
            for arm in arms:
                if arm not in reps[0]:
                    continue
                E = np.array([r[arm] for r in reps])
                b, sb, m, sm = summarise(E[:, :5], true)
                for kk, tr in enumerate(TR):
                    pre = f"item_{noise}_{n}_{arm}_{tr}"
                    put(f"{pre}_relbias_pct", float(b[kk])); put(f"{pre}_relbias_mcse_pct", float(sb[kk]))
                    put(f"{pre}_relrmse_pct", float(m[kk])); put(f"{pre}_relrmse_mcse_pct", float(sm[kk]))
                put(f"item_{noise}_{n}_{arm}_piR_mean_pct", 100 * float(E[:, 5].mean()))
                put(f"item_{noise}_{n}_{arm}_iters_median", float(np.median(E[:, 6])))
                for e in E:
                    rows.append(["item", arm, noise, 6, n] + list(e[:7]))
            # paired log-likelihood gain of strict over original rule (no item)
            gain = np.array([r["strict_none"][7] - r["orig_none"][7] for r in reps])
            put(f"item_{noise}_{n}_llgain_median", float(np.median(gain)))
            put(f"item_{noise}_{n}_llgain_min", float(np.min(gain)))
            put(f"item_{noise}_{n}_llgain_max", float(np.max(gain)))
    put("item_true_piR_pct", 100 * float(pi0[3]))

    # ---------------- 3. measurement
    # misclassification implied by the cut-off of 6 or more symptoms
    pb = rx.pmf_binom(Q_BIN); pbb = rx.pmf_betabinom(Q_BIN, RHO)
    put("q_S", float(Q_BIN[0])); put("q_A", float(Q_BIN[1])); put("q_I", float(Q_BIN[2])); put("rho", RHO)
    for nm, pm in [("bin", pb), ("bb", pbb)]:
        put(f"mis_{nm}_SR_pos_pct", 100 * float(pm[0, 6:].sum()))
        put(f"mis_{nm}_A_pos_pct", 100 * float(pm[1, 6:].sum()))
        put(f"mis_{nm}_I_neg_pct", 100 * float(pm[2, :6].sum()))
        put(f"mis_{nm}_A_band_pct", 100 * float(pm[1, 2:6].sum()))
        put(f"mis_{nm}_SR_band0_pct", 100 * float(pm[0, :2].sum()))
        put(f"mis_{nm}_mean_A", float((pm[1] * np.arange(10)).sum()))
    for nm in SD:
        put(f"mis_gauss_{nm}_SR_pos_pct", 100 * float(norm.sf(5.5, MU[0], SD[nm][0])))
        put(f"mis_gauss_{nm}_A_pos_pct", 100 * float(norm.sf(5.5, MU[1], SD[nm][1])))
        put(f"mis_gauss_{nm}_I_neg_pct", 100 * float(norm.cdf(5.5, MU[2], SD[nm][2])))
    # what misclassification does to a cross-section and to a naive two-wave remission estimate
    gm = np.random.default_rng(SEED + 1)
    for rates in RATES:
        P, pi0 = RATES[rates]
        n = 400_000
        s = rx.simulate_states(P, pi0, n, 2, gm)
        put(f"mis_{rates}_trueI_pct", 100 * float(pi0[2]))
        put(f"mis_{rates}_true_IR", float(P[2, 3]))
        for nm in ["bin", "bb", "gauss_low", "gauss_high"]:
            if nm == "bin":
                y = rx.emit_binom(s, Q_BIN, gm); c = y >= 6
            elif nm == "bb":
                y = rx.emit_betabinom(s, Q_BIN, RHO, gm); c = y >= 6
            else:
                y = rx.emit_gauss(s, MU, SD[nm.split("_")[1]], gm); c = y >= 5.5
            put(f"mis_{rates}_{nm}_obs_prev_pct", 100 * float(c[:, 0].mean()))
            put(f"mis_{rates}_{nm}_ppv_pct", 100 * float((s[c[:, 0], 0] == 2).mean()))
            exit_naive = float((~c[c[:, 0], 1]).mean())
            put(f"mis_{rates}_{nm}_naive_exit", exit_naive)
            put(f"mis_{rates}_{nm}_naive_exit_ratio", exit_naive / float(P[2, 3]))
    # recovery under each measurement arm
    marms = ["bin_binom", "bin_cat", "bin_band3", "bin_bin6", "bb_binom", "bb_cat", "bb_gauss"]
    P, pi0 = RATES["design"]
    true = np.array([P[i, j] for i, j in ENTRIES])
    for n in SIZES_MEAS:
        reps = cells[("meas", (n,))]
        put(f"meas_{n}_reps", len(reps))
        for arm in marms:
            E = np.array([r[arm] for r in reps])
            b, sb, m, sm = summarise(E[:, :5], true)
            for kk, tr in enumerate(TR):
                pre = f"meas_{n}_{arm}_{tr}"
                put(f"{pre}_relbias_pct", float(b[kk])); put(f"{pre}_relbias_mcse_pct", float(sb[kk]))
                put(f"{pre}_relrmse_pct", float(m[kk])); put(f"{pre}_relrmse_mcse_pct", float(sm[kk]))
            for e in E:
                rows.append(["meas", arm, "", 6, n] + list(e))

    # ---------------- figures
    cpu_exact = run_exact(q, t0)
    put("cpu_minutes_exact", cpu_exact / 60)
    with open(OUT / "recx_estimates.csv", "w") as f:
        f.write("analysis,setting,noise,waves,n,S_A,A_S,A_I,I_R,R_A,piR,iters,extra\n")
        for r in rows:
            f.write(",".join(str(x) for x in r) + "\n")
    figures_only()                              # figures are drawn from the saved replicate estimates
    put("runtime_minutes", (time.time() - t0) / 60)
    put("cpu_minutes_grid", tsum["grid"] / 60); put("cpu_minutes_item", tsum["item"] / 60)
    put("cpu_minutes_meas", tsum["meas"] / 60)
    put("n_fits_total", len(jobs))
    put("seed", SEED)
    json.dump(RX, open(OUT / "results_recx.json", "w"), indent=1, ensure_ascii=False)
    print(f"done in {(time.time() - t0) / 60:.1f} min, {len(RX)} results")


# ---------------------------------------------------------------- figures
ITEM_LABELS = {"orig_none": "Original rule, no item", "orig_perfect": "Original rule, perfect item",
               "strict_none": "Converged, no item", "strict_perfect": "Converged, perfect item",
               "strict_imperfect": "Converged, imperfect item",
               "strict_perfectfree": "Converged, perfect item,\nerror rates estimated"}
MEAS_LABELS = {"bin_binom": "Count, binomial fit", "bin_cat": "Count, categorical fit",
               "bin_band3": "Three bands (0 to 1, 2 to 5, 6 to 9)", "bin_bin6": "Classification only (6 or more)",
               "bb_binom": "Count, binomial fit (misspecified)", "bb_cat": "Count, categorical fit",
               "bb_gauss": "Count, Gaussian fit (misspecified)"}


def from_csv():
    """Rebuild the per-replicate estimates from outputs/recx_estimates.csv (for --figures-only)."""
    grid_E, item_E, meas_E = {}, {}, {}
    for line in open(OUT / "recx_estimates.csv").read().splitlines()[1:]:
        f = line.split(",")
        kind, setting, noise, T, n = f[0], f[1], f[2], int(f[3]), int(f[4])
        v = [float(x) for x in f[5:10]]
        if kind == "grid":
            grid_E.setdefault((setting, noise, T, n), []).append(v)
        elif kind == "item":
            item_E.setdefault((noise, n), {}).setdefault(setting, []).append(v)
        elif kind == "meas":
            meas_E.setdefault(n, {}).setdefault(setting, []).append(v)
    grid = {}
    for (rates, noise, T, n), E in grid_E.items():
        P, _ = RATES[rates]
        grid[(rates, noise, T, n)] = summarise(np.array(E), np.array([P[i, j] for i, j in ENTRIES]))
    return grid, item_E, meas_E


def fig_grid(grid, WAVES, SIZES):
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 5.0), sharex=True, sharey=True)
    for i, rates in enumerate(["design", "central"]):
        for j, T in enumerate(WAVES):
            ax = axs[i, j]
            for noise, ls in [("low", "-"), ("high", "--")]:
                for k, (nm, col) in enumerate(zip(TRL, COLS)):
                    ys = [grid[(rates, noise, T, n)][2][k] for n in SIZES]
                    ax.plot(SIZES, ys, ls=ls, color=col, marker="o", ms=2.5, lw=1.2,
                            label=nm if (noise == "low" and i == 0 and j == 0) else None)
            ax.axhline(10, color=ra.GREY, lw=0.6, ls=":")
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xticks(SIZES); ax.set_xticklabels([f"{n:,}" for n in SIZES], fontsize=8)
            ax.minorticks_off()
            ax.set_yticks([1, 3, 10, 30, 100, 300]); ax.set_yticklabels(["1", "3", "10", "30", "100", "300"], fontsize=8)
            ax.set_title(f"{'Design-study' if rates == 'design' else 'Calibrated central'} rates, {T} waves", fontsize=8)
            if j == 0:
                ax.set_ylabel("Relative RMSE (%)")
            if i == 1:
                ax.set_xlabel("Panel size")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=8, loc="lower center", ncol=5, bbox_to_anchor=(0.5, -0.01))
    fig.text(0.5, -0.045, "Solid lines: low noise. Dashed lines: high noise. Dotted line: 10%.", ha="center", fontsize=7.5)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(FIG / "fig_recx_grid.png"); plt.close(fig)


def fig_item(item_E, arms):
    """item_E: {(noise, n): {arm: [[five estimates], ...]}}"""
    P, _ = RATES["design"]
    true = np.array([P[i, j] for i, j in ENTRIES])
    panels = sorted(item_E, key=lambda k: (k[0] != "low", -k[1]))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    for ax, (tr_i, title) in zip(axs, [(0, "Incidence S→A"), (4, "Relapse R→A")]):
        for pi_, key in enumerate(panels):
            noise, n = key
            for a_i, arm in enumerate(arms):
                if arm not in item_E[key]:
                    continue
                b, sb, _, _ = summarise(np.array(item_E[key][arm]), true)
                yv = a_i + (pi_ - (len(panels) - 1) / 2) * 0.22
                col = [ra.C1, ra.C2, ra.C3][pi_]
                ax.errorbar(b[tr_i], yv, xerr=1.96 * sb[tr_i], fmt="o", ms=3.5, color=col, lw=1, capsize=0,
                            label=f"{'Low' if noise == 'low' else 'High'} noise, {n:,}" if (a_i == 0 and tr_i == 0) else None)
        ax.axvline(0, color=ra.INK, lw=0.6)
        ax.set_yticks(range(len(arms))); ax.set_yticklabels([ITEM_LABELS[a] for a in arms], fontsize=7.5)
        ax.set_title(title, fontsize=9); ax.set_xlabel("Relative bias (%)", fontsize=8)
    axs[0].invert_yaxis()
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, frameon=False, fontsize=7.5, loc="lower center", ncol=2, bbox_to_anchor=(0.6, -0.02))
    fig.tight_layout(rect=(0, 0.06, 1, 1)); fig.savefig(FIG / "fig_recx_baseline_item.png"); plt.close(fig)


def fig_emissions(pb, pbb):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), sharey=True)
    k = np.arange(10)
    for ax, pm, title in [(axs[0], pb, "Binomial (items independent within state)"),
                          (axs[1], pbb, f"Beta-binomial (item correlation {RHO})")]:
        for st, nm, col, off in [(0, "S and R", ra.GREY, -0.27), (1, "A", ra.C3, 0.0), (2, "I", ra.C2, 0.27)]:
            ax.bar(k + off, pm[st], width=0.26, color=col, label=nm)
        ax.axvline(5.5, color=ra.INK, lw=0.8, ls="--")
        ax.text(5.65, 0.45, "cut-off:\n6 or more", fontsize=7, ha="left", va="center")
        ax.set_xticks(k); ax.set_xlabel("Symptom count (0 to 9)"); ax.set_title(title, fontsize=8.5)
    axs[0].set_ylabel("Probability")
    axs[1].legend(frameon=False, fontsize=7.5, loc="upper right")
    fig.tight_layout(); fig.savefig(FIG / "fig_recx_emissions.png"); plt.close(fig)


def fig_meas(meas_E, marms):
    """meas_E: {n: {arm: [[five estimates], ...]}}; plots the largest panel."""
    P, _ = RATES["design"]
    true = np.array([P[i, j] for i, j in ENTRIES])
    n = max(meas_E)
    groups = [("Binomial data", [a for a in marms if a.startswith("bin_")]),
              ("Beta-binomial data", [a for a in marms if a.startswith("bb_")])]
    pos, labels, ypos = 0, [], {}
    heads = []
    for head, arms in groups:
        heads.append(pos); labels.append(head); pos += 1
        for a in arms:
            ypos[a] = pos; labels.append(MEAS_LABELS[a]); pos += 1
    fig, ax = plt.subplots(figsize=(7.2, 3.9))
    for arm in marms:
        _, _, m, sm = summarise(np.array(meas_E[n][arm]), true)
        for k, (nm, col) in enumerate(zip(TRL, COLS)):
            ax.errorbar(m[k], ypos[arm] + (k - 2) * 0.13, xerr=1.96 * sm[k], fmt="o", ms=3.2, color=col, lw=0.8,
                        capsize=0, label=nm if arm == marms[0] else None)
    ax.set_xscale("log"); ax.axvline(10, color=ra.GREY, lw=0.6, ls=":")
    ax.set_yticks(range(pos)); ax.set_yticklabels(labels, fontsize=7.5)
    for hpos in heads:
        ax.get_yticklabels()[hpos].set_fontweight("bold")
        ax.axhline(hpos, color=ra.GREY, lw=0.4, alpha=0.5)
    ax.set_ylim(pos - 0.5, -0.6)
    ax.set_xlabel(f"Relative RMSE (%), panel of {n:,}, six waves (log scale)")
    ax.legend(frameon=False, fontsize=7.5, loc="upper center", ncol=5, bbox_to_anchor=(0.4, 1.1))
    fig.tight_layout(); fig.savefig(FIG / "fig_recx_measurement.png"); plt.close(fig)


def figures_only():
    grid, item_E, meas_E = from_csv()
    waves = sorted({k[2] for k in grid}); sizes = sorted({k[3] for k in grid})
    fig_grid(grid, waves, sizes)
    fig_item(item_E, [a[0] for a in ITEM_ARMS])
    fig_emissions(rx.pmf_binom(Q_BIN), rx.pmf_betabinom(Q_BIN, RHO))
    fig_meas(meas_E, list(MEAS_LABELS))


if __name__ == "__main__":
    main()
