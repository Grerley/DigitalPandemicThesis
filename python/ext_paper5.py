"""Figure for Paper 5 (age and sex): girls' and boys' prevalence at HBSC ages, 2018 to 2042, under the
constant-vulnerability fit (Scenario P, central structure) and under the variant in which girls' vulnerability
changes in 2018, with the 2022 published figures and the 2018 values derived from the pooled odds ratios.

Writes figures(_bw)/fig_p5_sex and outputs/results_p5.json (keys p5_...).
"""
from __future__ import annotations

import json, pathlib
import numpy as np

import run_all as ra
from run_all import INK, C1, GREY, plt
from dpm import agesex as AS
import ext_agesex_shift as sh

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R: dict = {}


def put(k, v):
    R["p5_" + k] = {"value": float(v), "note": ""}


def series(y):
    o = [AS.hbsc_obs(y[:, j]) for j in range(y.shape[1])]
    return np.array([100 * x["girls"] for x in o]), np.array([100 * x["boys"] for x in o])


def ratios_2018_across_draws(n=100):
    """2018 girls-to-boys and age ratios of the constant-vulnerability fit, across the draws of ext_agesex."""
    import ext_agesex as E
    from dpm import calibrate as cal
    from dpm.model import Params
    st = AS.Struct()
    g = np.random.default_rng(E.SEED)
    ok = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]][:n]
    res = {"P": [], "T": []}
    for r in ok:
        d = {k: r[k] for k in cal.PARAM_ORDER}; base = Params(b=0.0, **d)
        stx = st.with_(eps_a=g.uniform(0.5, 0.9), eps_s=g.uniform(0.3, 0.7))
        targets = (g.uniform(0.065, 0.075), g.uniform(0.125, 0.135), g.uniform(0.085, 0.095), g.uniform(0.085, 0.095))
        for sc in "PT":
            sol = E.fit(base, stx, targets, sc)
            if sol is None:
                continue
            o = AS.hbsc_obs(AS.equilibrium(base.with_(b=sol["b0"]), stx.with_(phi=sol["phi"], psi=sol["psi"])))
            res[sc].append((o["girls"] / o["boys"], o["age11"] / o["age13_15"]))
    for sc, rows in res.items():
        a = np.array(rows)
        put(f"{sc}_n", len(a))
        for j, nm in enumerate(("sex", "age")):
            put(f"{sc}_ratio18_{nm}_median", np.median(a[:, j]))
            put(f"{sc}_ratio18_{nm}_p05", np.percentile(a[:, j], 5)); put(f"{sc}_ratio18_{nm}_p95", np.percentile(a[:, j], 95))
            put(f"{sc}_ratio18_{nm}_min", a[:, j].min()); put(f"{sc}_ratio18_{nm}_max", a[:, j].max())


def main():
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)
    base = pc1.with_(b=0.0); st = AS.Struct()
    asx = json.load(open(OUT / "results_agesex.json")); asxs = json.load(open(OUT / "results_asxs.json"))
    v = lambda d, k: d[k]["value"]
    # constant vulnerability, Scenario P, central structure
    stc = st.with_(phi=v(asx, "asx_c_P_phi"), psi=v(asx, "asx_c_P_psi"))
    p0 = base.with_(b=v(asx, "asx_c_P_b0"))
    xc18 = AS.equilibrium(p0, stc)
    tt = np.arange(0, 24.001, 0.25)
    _, yc = AS.simulate(p0.with_(b=v(asx, "asx_c_P_k2")), stc, xc18, 24, t_eval=tt)
    gc, bc = series(yc)
    # vulnerability of girls changes in 2018
    sol = {k: v(asxs, f"asxs_c_{k}") for k in ("b0", "b1", "phi0", "phi1", "psi")}
    st_a = st.with_(phi=sol["phi0"], psi=sol["psi"]); st_b = st.with_(phi=sol["phi1"], psi=sol["psi"])
    q0 = base.with_(b=sol["b0"]); xs18 = AS.equilibrium(q0, st_a)
    _, ys = AS.simulate(q0.with_(b=sol["b1"]), st_b, xs18, 24, t_eval=tt)
    gs, bs = series(ys)
    put("const_2022_girls", gc[16]); put("const_2022_boys", bc[16]); put("shift_2022_girls", gs[16]); put("shift_2022_boys", bs[16])
    d18g, d18b = v(asxs, "asxs_derived18_girls"), v(asxs, "asxs_derived18_boys")
    yrs = 2018 + tt
    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.2), sharey=True)
    for ax, (g_, b_), title in [(axs[0], (gc, bc), "A. Constant vulnerability of girls"),
                                (axs[1], (gs, bs), "B. Vulnerability of girls changes in 2018")]:
        ax.plot(yrs, g_, color=INK, ls="-", lw=1.4, label="Girls, model")
        ax.plot(yrs, b_, color=C1, ls="--", lw=1.4, label="Boys, model")
        ax.errorbar([2022], [13], yerr=[[0.5], [0.5]], fmt="o", color=INK, ms=4, capsize=2, label="Girls, HBSC 2022")
        ax.errorbar([2022], [9], yerr=[[0.5], [0.5]], fmt="s", color=C1, ms=4, capsize=2, label="Boys, HBSC 2022")
        ax.plot([2018], [d18g], marker="D", ls="none", color=INK, ms=5, mfc="white", label="Girls, 2018 (derived)")
        ax.plot([2018], [d18b], marker="^", ls="none", color=C1, ms=5, mfc="white", label="Boys, 2018 (derived)")
        ax.set_title(title, fontsize=9, loc="left"); ax.set_xlabel("Year"); ax.set_ylim(0, 18)
    axs[0].set_ylabel("Problematic use at HBSC ages (%)")
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, frameon=False, fontsize=7.5)
    fig.tight_layout(rect=(0, 0.14, 1, 1)); fig.savefig(FIG / "fig_p5_sex.png", dpi=200); plt.close(fig)
    ratios_2018_across_draws()
    json.dump(R, open(OUT / "results_p5.json", "w"), indent=1)
    print({k: round(x["value"], 3) for k, x in R.items()})


if __name__ == "__main__":
    main()
