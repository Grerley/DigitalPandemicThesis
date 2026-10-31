"""Misclassification by the scale (Chapter 7, section "Misclassification by the scale").

If the Social Media Disorder Scale classifies a young person with problematic use as such with
probability Se (sensitivity) and one without it as such with probability 1 - Sp, the published prevalence
y relates to the model's prevalence I by y = Se I + (1 - Sp)(1 - I). The script inverts this for the two
HBSC figures over a grid of (Se, Sp), re-calibrates both scenarios at the central rates and on a subset of
the calibrated draws, and reports projections on both the true and the observed scale.

Writes outputs/results_misx.json (keys misx_...) and figures/fig_misx_grid.png.
"""
from __future__ import annotations

import json, pathlib, time
import numpy as np

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt
from dpm.model import Params, equilibrium, simulate, R0
from dpm import calibrate as cal

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R: dict = {}
SE = [1.0, 0.9, 0.8]
SP = [1.0, 0.99, 0.98, 0.97, 0.96]
N_SUB = 200
SEED = 20261105


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["misx_" + k] = {"value": v, "note": note}


def true_prev(y, se, sp):
    return (y - (1 - sp)) / (se - (1 - sp))


def observed(I, se, sp):
    return se * I + (1 - sp) * (1 - I)


def tag(se, sp):
    return f"se{int(round(100 * se))}_sp{int(round(100 * sp))}"


def calib(base: Params, t18, t22):
    b0 = cal.solve_b0(base, t18)
    if b0 is None:
        return None
    p0 = base.with_(b=b0); x18 = equilibrium(p0)
    b1 = cal.solve_b1(p0, x18, t22); da = cal.solve_dalpha(p0, x18, t22)
    if b1 is None or da is None:
        return None
    t = np.arange(0, 25.0)
    _, yP = simulate(p0.with_(b=b1), x18, 24.0, t_eval=t)
    _, yT = simulate(p0, x18, 24.0, t_eval=t, alpha_t=cal.shock(p0.alpha, da))
    return dict(b0=b0, b1=b1, da=da, R0_b1=R0(p0.with_(b=b1, alpha=0.0)), P2042=yP[2, -1], T2042=yT[2, -1],
                P2026=yP[2, 8], T2026=yT[2, 8], eqP=equilibrium(p0.with_(b=b1))[2])


def figure():
    Rr = R if R else json.load(open(OUT / "results_misx.json"))
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    cols = {1.0: INK, 0.9: C1, 0.8: C2}
    fp = [100 * (1 - s) for s in SP]
    for se in SE:
        ax[0].plot(fp, [100 * true_prev(0.07, se, s) for s in SP], "o-", color=cols[se], label=f"2018, Se = {se:.1f}")
        ax[0].plot(fp, [100 * true_prev(0.11, se, s) for s in SP], "s--", color=cols[se], label=f"2022, Se = {se:.1f}")
        ax[1].plot(fp, [Rr[f"misx_{tag(se, s)}_central_R0_b1"]["value"] for s in SP], "o-", color=cols[se], label=f"Se = {se:.1f}")
    ax[0].set_xlabel("False-positive rate, 1 - Sp (%)"); ax[0].set_ylabel("Implied model prevalence (%)")
    ax[0].legend(frameon=False, fontsize=7, ncol=3, loc="lower left", bbox_to_anchor=(0, 1.0)); ax[0].grid(alpha=0.3)
    ax[1].set_xlabel("False-positive rate, 1 - Sp (%)"); ax[1].set_ylabel(r"$\mathcal{R}_0$ after 2018 (with $\alpha = 0$)")
    ax[1].legend(frameon=False, fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(FIG / "fig_misx_grid.png", dpi=200); plt.close(fig)


def main():
    t0 = time.time()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da_c, x18_c = ra.central(ranges)
    base = pc0.with_(b=0.0)
    ok = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]]
    rng = np.random.default_rng(SEED)
    sub = [ok[i] for i in rng.choice(len(ok), size=N_SUB, replace=False)]
    put("se_grid", SE); put("sp_grid", SP); put("n_sub", N_SUB)
    grid_P = np.full((len(SE), len(SP)), np.nan); grid_ratio = np.full((len(SE), len(SP)), np.nan)
    for i, se in enumerate(SE):
        for j, sp in enumerate(SP):
            tg = tag(se, sp)
            t18, t22 = true_prev(0.07, se, sp), true_prev(0.11, se, sp)
            put(f"{tg}_true2018_pct", 100 * t18); put(f"{tg}_true2022_pct", 100 * t22)
            put(f"{tg}_ratio", t22 / t18)
            grid_ratio[i, j] = t22 / t18
            c = calib(base, t18, t22)
            put(f"{tg}_central_ok", c is not None)
            if c is not None:
                for k in ["b0", "b1", "da", "R0_b1"]:
                    put(f"{tg}_central_{k}", c[k])
                for k in ["P2042", "T2042", "P2026", "T2026", "eqP"]:
                    put(f"{tg}_central_{k}_true_pct", 100 * c[k])
                    put(f"{tg}_central_{k}_obs_pct", 100 * observed(c[k], se, sp))
                put(f"{tg}_central_rel_rise_P_true", c["P2042"] / t18)
                grid_P[i, j] = 100 * observed(c["P2042"], se, sp)
            # across a subset of calibrated draws (their own rates, re-solved exposure)
            if (se, sp) in [(1.0, 1.0), (0.9, 0.98), (0.8, 0.96), (0.9, 0.96)]:
                rows = []
                for r in sub:
                    b = Params(b=0.0, **{k: r[k] for k in cal.PARAM_ORDER})
                    c2 = calib(b, t18, t22)
                    if c2 is not None:
                        rows.append(c2)
                put(f"{tg}_draws_n_ok", len(rows))
                for k in ["P2042", "T2042"]:
                    v = np.array([observed(x[k], se, sp) for x in rows]) * 100
                    put(f"{tg}_draws_{k}_obs_median", float(np.median(v)))
                    put(f"{tg}_draws_{k}_obs_p05", float(np.percentile(v, 5)))
                    put(f"{tg}_draws_{k}_obs_p95", float(np.percentile(v, 95)))
                v = np.array([x["R0_b1"] for x in rows])
                put(f"{tg}_draws_R0_b1_median", float(np.median(v)))
                print("draws", tg, time.time() - t0, flush=True)
    figure()
    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_misx.json", "w"), indent=1, ensure_ascii=False)
    print("done", len(R), (time.time() - t0) / 60)


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--figure":
        figure()
    else:
        main()
