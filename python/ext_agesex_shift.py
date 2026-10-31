"""Age- and sex-structured model: out-of-sample checks against 2018 and a variant in which girls'
vulnerability changes in 2018 (thesis Chapter 7).

Writes outputs/results_asxs.json (keys prefixed asxs_).

The 2018 HBSC report gives no prevalence by sex or age. Boer et al. (2022a, Table 2) report, for the pooled
2017/2018 sample (n = 222,532 in 44 countries), adjusted odds ratios of problematic use of 1.208 for girls
against boys, and 1.484 and 1.612 for 13- and 15-year-olds against 11-year-olds, with mean prevalence 7.64%.
Converting these to prevalences, with equal numbers in each sex-by-age cell and no interaction, gives an
approximate 2018 girls-to-boys prevalence ratio of 1.19 and an age-11-to-ages-13-and-15 ratio of 0.66. These are
derived quantities, used here only (1) to check the constant-vulnerability model out of sample and (2) as one
extra target for a variant in which girls' vulnerability can change in 2018.

    python ext_agesex_shift.py
"""
from __future__ import annotations

import json, math, pathlib, time
import numpy as np
from scipy.optimize import least_squares
from scipy.special import expit, logit
from scipy.optimize import brentq

import run_all as ra
from run_all import INK, C1, C2, C3, GREY, plt
from dpm.model import Params
from dpm import agesex as AS
from dpm import calibrate as cal

HERE = pathlib.Path(__file__).parent
OUT = HERE / "outputs"; FIG = HERE / "figures"
R = {}
SEED = ra.SEED + 1700


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["asxs_" + k] = {"value": v, "note": note}


def derived_2018():
    """Prevalences implied by the Boer et al. (2022a) pooled odds ratios, equal cells, mean 7.64%."""
    lor = {"girl": math.log(1.208), "13": math.log(1.484), "15": math.log(1.612)}
    cells = [(s, a) for s in (0, 1) for a in ("11", "13", "15")]
    def lin(c, s, a):
        return c + (lor["girl"] if s == 0 else 0) + (lor[a] if a != "11" else 0)
    c0 = brentq(lambda c: np.mean([expit(lin(c, s, a)) for s, a in cells]) - 0.0764, -10, 0)
    p = {(s, a): expit(lin(c0, s, a)) for s, a in cells}
    girls = np.mean([p[(0, a)] for a in ("11", "13", "15")]); boys = np.mean([p[(1, a)] for a in ("11", "13", "15")])
    age11 = np.mean([p[(s, "11")] for s in (0, 1)]); age1315 = np.mean([p[(s, a)] for s in (0, 1) for a in ("13", "15")])
    return dict(ratio_sex=girls / boys, ratio_age=age11 / age1315, girls=girls, boys=boys, age11=age11, age1315=age1315)


def fit_shift(base, st0, targets, z0=None):
    """Unknowns b0, b1, phi0 (to 2018), phi1 (from 2018), psi. Targets: 2018 total and girls/boys ratio;
    2022 girls, boys, 11-year-olds."""
    t18, rs18, tg, tb, ta = targets

    def resid(z):
        b0, b1, phi0, phi1, psi = np.exp(z)
        st_a = st0.with_(phi=phi0, psi=psi); st_b = st0.with_(phi=phi1, psi=psi)
        p0 = base.with_(b=b0)
        x18 = AS.equilibrium(p0, st_a)
        o18 = AS.hbsc_obs(x18)
        _, y = AS.simulate(p0.with_(b=b1), st_b, x18, 4, t_eval=[4])
        o22 = AS.hbsc_obs(y[:, -1])
        return np.array([o18["total"] - t18, (o18["girls"] / o18["boys"] - rs18) * 0.1,
                         o22["girls"] - tg, o22["boys"] - tb, o22["age11"] - ta]) * 100
    if z0 is None:
        z0 = np.log([1.0, 1.6, 1.3, 2.4, 1.3])
    res = least_squares(resid, z0, method="lm", xtol=1e-12, ftol=1e-12, max_nfev=600)
    if not res.success or np.max(np.abs(res.fun)) > 1e-4:
        return None
    b0, b1, phi0, phi1, psi = np.exp(res.x)
    return dict(b0=b0, b1=b1, phi0=phi0, phi1=phi1, psi=psi)


def project_shift(base, st0, sol, years=24):
    st_a = st0.with_(phi=sol["phi0"], psi=sol["psi"]); st_b = st0.with_(phi=sol["phi1"], psi=sol["psi"])
    p0 = base.with_(b=sol["b0"]); x18 = AS.equilibrium(p0, st_a)
    t, y = AS.simulate(p0.with_(b=sol["b1"]), st_b, x18, years, t_eval=np.arange(0, years + 1, 1.0))
    return x18, t, y


def main():
    t0 = time.time()
    d18 = derived_2018()
    for k, v in d18.items():
        put(f"derived18_{k}", v if "ratio" in k else 100 * v)
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)
    base = pc1.with_(b=0.0); st = AS.Struct()
    # out-of-sample check of the constant-vulnerability model (central fit from ext_agesex)
    asx = json.load(open(OUT / "results_agesex.json"))
    for sc in ["P", "T"]:
        phi, psi, b0 = (asx[f"asx_c_{sc}_{k}"]["value"] for k in ("phi", "psi", "b0"))
        x = AS.equilibrium(base.with_(b=b0), st.with_(phi=phi, psi=psi))
        o = AS.hbsc_obs(x)
        put(f"const_{sc}_2018_ratio_sex", o["girls"] / o["boys"])
        put(f"const_{sc}_2018_ratio_age", o["age11"] / o["age13_15"])
        put(f"const_{sc}_2018_girls", 100 * o["girls"]); put(f"const_{sc}_2018_boys", 100 * o["boys"])
    # variant: girls' vulnerability changes in 2018 (Scenario P structure)
    targets = (0.07, d18["ratio_sex"], 0.13, 0.09, 0.09)
    sol = fit_shift(base, st, targets)
    for k, v in sol.items():
        put(f"c_{k}", v)
    x18s, t, y = project_shift(base, st, sol)
    o18 = AS.hbsc_obs(x18s)
    put("c_2018_ratio_age", o18["age11"] / o18["age13_15"])
    put("c_2018_girls", 100 * o18["girls"]); put("c_2018_boys", 100 * o18["boys"])
    for yr in [2022, 2026, 2032, 2042]:
        o = AS.hbsc_obs(y[:, int(yr - 2018)])
        put(f"c_total_{yr}", 100 * o["total"]); put(f"c_girls_{yr}", 100 * o["girls"]); put(f"c_boys_{yr}", 100 * o["boys"])
        put(f"c_whole_{yr}", 100 * AS.whole_prev(y[:, int(yr - 2018)]))
    # across draws
    g = np.random.default_rng(SEED)
    ok = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]][:60]
    rows = []
    for r in ok:
        d = {k: r[k] for k in cal.PARAM_ORDER}
        b = Params(b=0.0, **d)
        stx = st.with_(eps_a=g.uniform(0.5, 0.9), eps_s=g.uniform(0.3, 0.7))
        tg = (g.uniform(0.065, 0.075), g.uniform(d18["ratio_sex"] - 0.08, d18["ratio_sex"] + 0.08),
              g.uniform(0.125, 0.135), g.uniform(0.085, 0.095), g.uniform(0.085, 0.095))
        s = fit_shift(b, stx, tg)
        if s is None:
            continue
        _, _, yy = project_shift(b, stx, s)
        o = AS.hbsc_obs(yy[:, -1])
        rows.append(dict(phi0=s["phi0"], phi1=s["phi1"], total42=o["total"], girls42=o["girls"], boys42=o["boys"]))
    put("d_n", len(ok)); put("d_n_ok", len(rows))
    for k in ["phi0", "phi1", "total42", "girls42", "boys42"]:
        v = np.array([x[k] for x in rows]); sc = 1.0 if k.startswith("phi") else 100.0
        put(f"d_{k}_median", sc * float(np.median(v))); put(f"d_{k}_p05", sc * float(np.percentile(v, 5)))
        put(f"d_{k}_p95", sc * float(np.percentile(v, 95)))
    put("runtime_minutes", (time.time() - t0) / 60)
    json.dump(R, open(OUT / "results_asxs.json", "w"), indent=1, ensure_ascii=False)
    print(f"done, {len(R)} results, {(time.time() - t0) / 60:.1f} min")
    print({k: round(v["value"], 3) for k, v in R.items() if isinstance(v["value"], float)})


if __name__ == "__main__":
    main()
