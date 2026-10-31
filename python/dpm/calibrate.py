"""Calibration of the S-A-I-R model to the published HBSC problematic social media use series
(thesis Chapter 7).

Targets (Boniel-Nissim et al., 2024, pp. iv and 6): problematic social media use was 7% in 2018 and
11% in 2022 across HBSC countries and regions. Both figures are published rounded to whole per cent,
so each target is drawn uniformly within +/- 0.5 percentage points.

Two aggregate figures can identify at most two parameters. All other parameters are drawn from
literature-based ranges (Latin hypercube) and, for each draw, two quantities are solved exactly:

  Scenario P (persistent shift): the population is at equilibrium in 2018 with exposure rate b0;
      from 2018 the exposure rate is b1. Solve b0 so that I*(b0) = target_2018 and b1 so that
      I(2022) = target_2022.
  Scenario T (transient shock): exposure stays at b0; exogenous acquisition rises by d_alpha during
      2020 and 2021 only (the COVID-19 period) and then returns. Solve d_alpha so that I(2022) = target_2022.

A draw is kept only if both scenarios can reproduce both targets. The solver bounds (b up to 200 per
year, d_alpha up to 20 per year) are far above any plausible value; a draw fails when, for example,
recovery is so fast that even unbounded exposure cannot raise prevalence from 7% to 11% in four years.
Kept draws reproduce the data exactly. They differ in what they imply after 2022, which is the
identifiability point the thesis makes.
"""
from __future__ import annotations

import json
import numpy as np
from scipy.optimize import brentq
from scipy.stats import qmc

from .model import Params, equilibrium, simulate, R0

PARAM_ORDER = ["eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"]


def lhs_draws(ranges: dict, n: int, seed: int):
    sampler = qmc.LatinHypercube(d=len(PARAM_ORDER), seed=seed)
    u = sampler.random(n)
    out = []
    for row in u:
        d = {}
        for j, name in enumerate(PARAM_ORDER):
            lo, hi = ranges[name]["low"], ranges[name]["high"]
            if ranges[name].get("log", False):
                d[name] = float(np.exp(np.log(lo) + row[j] * (np.log(hi) - np.log(lo))))
            else:
                d[name] = float(lo + row[j] * (hi - lo))
        out.append(d)
    return out


def I_eq(p: Params) -> float:
    return float(equilibrium(p)[2])


def solve_b0(base: Params, target: float, bmax: float = 200.0):
    f = lambda b: I_eq(base.with_(b=b)) - target
    if f(0.0) > 0 or f(bmax) < 0:
        return None
    return brentq(f, 0.0, bmax, xtol=1e-10)


def I_after(p: Params, x0, years: float, alpha_t=None) -> float:
    _, y = simulate(p, x0, years, t_eval=[years], alpha_t=alpha_t, rtol=1e-9, atol=1e-11)
    return float(y[2, -1])


def solve_b1(p0: Params, x0, target: float, years: float = 4.0, bmax: float = 200.0):
    f = lambda b: I_after(p0.with_(b=b), x0, years) - target
    if f(p0.b) > 0 or f(bmax) < 0:
        return None
    return brentq(f, p0.b, bmax, xtol=1e-9)


def shock(alpha0: float, d_alpha: float, start: float = 2.0, end: float = 4.0):
    """Exogenous acquisition raised by d_alpha between 2020 (t=2) and 2022 (t=4)."""
    return lambda t: alpha0 + (d_alpha if start <= t < end else 0.0)


def solve_dalpha(p0: Params, x0, target: float, dmax: float = 20.0):
    f = lambda d: I_after(p0, x0, 4.0, alpha_t=shock(p0.alpha, d)) - target
    if f(0.0) > 0 or f(dmax) < 0:
        return None
    return brentq(f, 0.0, dmax, xtol=1e-10)


def calibrate(ranges: dict, n: int = 1000, seed: int = 20261031):
    rng = np.random.default_rng(seed)
    draws = lhs_draws(ranges, n, seed)
    rows = []
    for d in draws:
        t18 = rng.uniform(0.065, 0.075)
        t22 = rng.uniform(0.105, 0.115)
        base = Params(b=0.0, **d)
        b0 = solve_b0(base, t18)
        if b0 is None:
            rows.append(dict(d, t18=t18, t22=t22, ok=False, reason="2018 level not reachable")); continue
        p0 = base.with_(b=b0)
        x18 = equilibrium(p0)
        b1 = solve_b1(p0, x18, t22)
        da = solve_dalpha(p0, x18, t22)
        reason = None
        if b1 is None and da is None:
            reason = "2022 level not reachable in either scenario"
        elif b1 is None:
            reason = "2022 level not reachable in Scenario P"
        elif da is None:
            reason = "2022 level not reachable in Scenario T"
        rows.append(dict(d, t18=t18, t22=t22, ok=reason is None, reason=reason, b0=b0, b1=b1, d_alpha=da,
                         R0_b0=R0(p0.with_(alpha=0.0)),
                         R0_b1=R0(p0.with_(b=b1, alpha=0.0)) if b1 else None,
                         x18=x18.tolist()))
    return rows


def save(rows, path):
    with open(path, "w") as f:
        json.dump(rows, f, indent=1)
