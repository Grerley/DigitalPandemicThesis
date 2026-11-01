"""Numbers and one figure specific to Paper 3 (networks), from the saved threshold sweep of ext_networks.py.

1. Sensitivity of the extrapolated effective threshold to the cut-off on the long-run share (1%, 3%, 5%), and an
   alternative summary: the lowest grid value at which at least half the runs persisted for 100 years.
2. The heuristic without the back-edge correction, T <k^2>/<k>, evaluated at the effective thresholds.
3. Figure: long-run share At-Risk, Addicted or Recovered against R_net and against the population R0, with the
   population-model share 1 - 1/R0 (the quantity reported in the threshold table).

Writes outputs/results_p3.json (keys p3_...) and figures(_bw)/fig_p3_threshold.
"""
from __future__ import annotations

import json, os, pathlib
import numpy as np

import run_all as ra
from run_all import INK, C1, C2, GREY, plt

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R: dict = {}
GRID = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5]
NAMES = {"ER": "Erdős–Rényi", "WS": "Watts–Strogatz", "BA": "Barabási–Albert"}


def put(k, v):
    R["p3_" + k] = {"value": float(v) if isinstance(v, (np.floating, float)) else v, "note": ""}


def tag(t):
    s = f"{t:g}".replace(".", "p")
    return s


def main():
    res = json.load(open(OUT / "results_netx.json"))
    val = lambda k: res[k]["value"]
    curves = {}
    for kind in NAMES:
        tot = np.array([val(f"netx_thr_{kind}_{tag(t)}_tot_pct") for t in GRID]) / 100
        per = np.array([val(f"netx_thr_{kind}_{tag(t)}_persist_pct") for t in GRID])
        r0 = np.array([val(f"netx_thr_{kind}_{tag(t)}_R0mf") for t in GRID])
        curves[kind] = (tot, per, r0)
        for cut in (0.01, 0.03, 0.05):
            ok = [i for i, t in enumerate(GRID) if tot[i] >= cut][:3]
            x = np.array(GRID)[ok]; y = tot[ok]
            slope, icpt = np.polyfit(x, y, 1)
            rs = -icpt / slope
            put(f"{kind}_Rnet_star_cut{int(cut * 100)}", rs)
            put(f"{kind}_R0mf_star_cut{int(cut * 100)}", float(np.interp(rs, GRID, r0)))
        half = next(t for i, t in enumerate(GRID) if per[i] >= 50)
        put(f"{kind}_half_persist_Rnet", half)
        k = val(f"netx_thr_{kind}_k2_over_k")
        put(f"{kind}_noback_factor", k / (k - 1))
        put(f"{kind}_noback_at_star", val(f"netx_thr_{kind}_Rnet_star") * k / (k - 1))
    put("check_cut3_matches_ER", abs(R["p3_ER_Rnet_star_cut3"]["value"] - val("netx_thr_ER_Rnet_star")))
    put("check_cut3_matches_WS", abs(R["p3_WS_Rnet_star_cut3"]["value"] - val("netx_thr_WS_Rnet_star")))
    put("check_cut3_matches_BA", abs(R["p3_BA_Rnet_star_cut3"]["value"] - val("netx_thr_BA_Rnet_star")))

    fig, axs = plt.subplots(1, 2, figsize=(7.4, 3.3), sharey=True)
    styles = {"ER": ("o", "-"), "WS": ("s", ":"), "BA": ("^", "-.")}
    for kind, (tot, per, r0) in curves.items():
        mk, ls = styles[kind]
        axs[0].plot(GRID, 100 * tot, marker=mk, ls=ls, ms=3, color=INK, lw=1.2, label=NAMES[kind])
        axs[1].plot(r0, 100 * tot, marker=mk, ls=ls, ms=3, color=INK, lw=1.2)
    r0s = np.linspace(1.0, 4.2, 200)
    axs[1].plot(np.concatenate([[0.2, 1.0], r0s]), 100 * np.concatenate([[0, 0], 1 - 1 / r0s]), color=INK, lw=1.0,
                ls="--", label="Population model")
    for ax in axs:
        ax.axvline(1.0, color=GREY, lw=0.8, ls=":")
        ax.set_ylim(0, None)
    axs[0].set_xlabel(r"Heuristic $\mathcal{R}_{\mathrm{net}}$")
    axs[1].set_xlabel(r"Population-model $\mathcal{R}_0$ at the same exposure rate")
    axs[0].set_ylabel("At-Risk, Addicted or Recovered,\nyears 90 to 100 (%)")
    axs[0].legend(frameon=False, fontsize=8); axs[1].legend(frameon=False, fontsize=8, loc="lower right")
    axs[0].set_title("A. Against the network heuristic", fontsize=9, loc="left")
    axs[1].set_title("B. Against the population threshold", fontsize=9, loc="left")
    fig.tight_layout(); fig.savefig(FIG / "fig_p3_threshold.png", dpi=200); plt.close(fig)
    json.dump(R, open(OUT / "results_p3.json", "w"), indent=1)
    for k, v in R.items():
        print(k, round(v["value"], 3))


if __name__ == "__main__":
    main()
