"""Interpretation of the model's parameters at the calibrated central point (thesis Chapter 4).

Writes outputs/results_modx.json (keys prefixed modx_).

For one young person, the model's rates define a continuous-time Markov chain on {S, A, I, R} with an
absorbing exit at age 25 (rate nu). At a constant social exposure lambda* (the post-2018 equilibrium value)
this chain gives durations, exit probabilities, the expected number of problematic-use spells and the
probability of ever having problematic use between 10 and 25.

    python ext_model.py
"""
from __future__ import annotations

import json, pathlib
import numpy as np
from scipy.linalg import expm

import run_all as ra
from dpm.model import Params, equilibrium

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"
R = {}


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["modx_" + k] = {"value": v, "note": note}


def individual_generator(p: Params, hazard_SA: float, absorbing_I: bool = False):
    """Generator on (S, A, I, R, Exit) for one person. With absorbing_I, entry to I is absorbing
    (used for the probability of ever having problematic use)."""
    Q = np.zeros((5, 5))
    Q[0, 1] = hazard_SA; Q[0, 4] = p.nu
    Q[1, 2] = p.sigma; Q[1, 0] = p.gamma1; Q[1, 4] = p.nu
    if not absorbing_I:
        Q[2, 3] = p.gamma2; Q[2, 4] = p.nu
        Q[3, 1] = p.delta; Q[3, 4] = p.nu
    for i in range(5):
        Q[i, i] = -Q[i].sum()
    return Q


def main():
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)
    p = pc1
    a = p.sigma + p.gamma1 + p.nu; c = p.gamma2 + p.nu; d = p.delta + p.nu
    put("dur_A_years", 1 / a); put("dur_I_years", 1 / c); put("dur_R_years", 1 / d)
    put("p_A_to_I", p.sigma / a); put("p_A_to_S", p.gamma1 / a); put("p_A_exit", p.nu / a)
    put("p_I_to_R", p.gamma2 / c); put("p_R_to_A", p.delta / d)
    q = (p.sigma / a) * (p.gamma2 / c) * (p.delta / d)
    put("q_loop", q); put("expected_A_spells_per_entry", 1 / (1 - q))
    put("expected_I_spells_per_entry", (p.sigma / a) / (1 - q))
    x = equilibrium(p)
    lam = p.b * (x[2] + p.eta * x[1])
    hz = lam + p.alpha
    put("lambda_star", lam); put("hazard_SA_star", hz)
    put("share_SA_social_pct", 100 * lam / hz)
    for tag, h in [("b1", hz), ("b0", pc0.b * (equilibrium(pc0)[2] + pc0.eta * equilibrium(pc0)[1]) + pc0.alpha)]:
        # probability of ever having problematic use between 10 and 25, starting in S at 10, at constant exposure;
        # a young person observed from age 10 to 25: no random exit within the window (nu = 0)
        p_life = p.with_(nu=0.0)
        Q = individual_generator(p_life, h, absorbing_I=True)
        P15 = expm(Q * 15.0)
        put(f"ever_I_by25_{tag}_pct", 100 * P15[0, 2])
        Qall = individual_generator(p_life, h)
        P = expm(Qall * 15.0)
        put(f"in_I_at25_{tag}_pct", 100 * P[0, 2])
        # probability of ever having problematic use by each age (for the text)
        for age in (12, 13, 16, 20):
            put(f"ever_I_by{age}_{tag}_pct", 100 * expm(Q * (age - 10.0))[0, 2])
        # expected years spent in I between 10 and 25
        grid = np.linspace(0, 15, 301)
        occ = [expm(Qall * t)[0, 2] for t in grid]
        put(f"years_in_I_{tag}", float(np.trapezoid(occ, grid)))
        put(f"hazard_SA_{tag}", h)
    # mean time from entering A to first entering I, among those who do (competing risks, first spell)
    put("mean_A_to_I_given_progress_years", 1 / a)
    json.dump(R, open(OUT / "results_modx.json", "w"), indent=1, ensure_ascii=False)
    print({k: round(v["value"], 4) for k, v in R.items()})


if __name__ == "__main__":
    main()
