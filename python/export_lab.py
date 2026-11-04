"""Export the thesis results that the Interactive Simulation Lab displays (app/).

Writes app/js/model/lab-data.js, an ES module `export const LAB = {...}`. Everything the lab shows as a thesis
result comes from this file: the parameter ranges, the calibrated central values, every kept calibration draw,
the projection percentiles, the likelihood ratios for the next HBSC wave, intervention, network, control,
parameter-recovery and age-sex results, and reference values the lab's JavaScript model is tested against.
Run after run_all.py, derived.py and the extension scripts:

    python export_lab.py
"""
from __future__ import annotations

import csv, json, pathlib, re, subprocess, datetime
import numpy as np

import run_all as ra
from dpm import calibrate as cal
from dpm.model import Params, equilibrium, simulate, R0, relapse_multiplier
from dpm import network as nw
from dpm import agesex as AS

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"
APP = HERE.parent / "app" / "js" / "model" / "lab-data.js"
RESEARCH = HERE.parents[1] / "research"          # thesis/research (verified references), when present


def results():
    R = {}
    for f in sorted(OUT.glob("results*.json")):
        R.update(json.load(open(f)))
    return {k: v["value"] for k, v in R.items()}


KEEP = re.compile(r"^(R0_b[01]$|I_star_b_zero$|eq_b1_I$|elast_I_|ident_n_|cal_n_|cal_kept_median_|cal_rej_median_|"
                  r"proj_[PT]_20(22|26|32|42)_|proj_gap_2042$|p2_(implicit|flat|log)_(lr|pP|pT)_\d+$|"
                  r"int_.*_(median|p05|p95)$|int_start_year$|ctl_mv_ratio|ctl_ref_(const_share_pct|gain_pct|averted_opt_pct|"
                  r"I_2042_(none|opt|const)_pct)$|net_(ER|WS|BA)_|net_meanfield_|net_mean_degree$|net_nodes$|net_reps$|"
                  r"netx_tgt_(cal|r1p5)_(hub|random|acq)_c(5|10|20|30)_(rel_pct|se)$|netx_tgt_(cal|r1p5)_Rnet$|netx_tgt_(cal|r1p5)_base$|"
                  r"asx_c_P_(phi|psi|b0|k2|girls_20(22|42)|boys_20(22|42)|total_20(22|42))$|asx_d_P_(phi|psi)_(median|p05|p95)$|"
                  r"asx_d_P_(girls42|boys42|gap42)_(median|p05|p95)$|asxs_derived18_|asx_ma_c(30|50|70)_|asx_ma_d_|"
                  r"p5_[PT]_ratio18_(sex|age)_(median|p05|p95)$|mathx_ageing_R0_|rec_waves$|seed$|asx_start_year$|env_)")


def r6(x):
    return float(f"{x:.6g}")


def main():
    R = results()
    ranges = json.load(open(HERE / "data" / "parameter_ranges.json"))
    pc0, pc1, da, x18 = ra.central(ranges)

    # ---- calibration draws (kept) and projection percentiles
    draws = [r for r in json.load(open(OUT / "calibration_draws.json")) if r["ok"]]
    keys = cal.PARAM_ORDER + ["b0", "b1", "d_alpha", "t18", "t22"]
    D = [[r6(r[k]) for k in keys] + [r6(v) for v in r["x18"]] for r in draws]
    proj = {}
    for sc in "PT":
        A = np.loadtxt(OUT / f"projection_{sc}.csv", delimiter=",")
        proj[sc] = {q: [r6(v) for v in np.percentile(A, q, axis=0)] for q in (5, 50, 95)}

    # ---- reference values for the JavaScript self-test (central parameters)
    t = np.arange(0, 25.0)
    _, yP = simulate(pc1, x18, 24, t_eval=t)
    _, yT = simulate(pc0, x18, 24, t_eval=t, alpha_t=cal.shock(pc0.alpha, da))
    beta = nw.beta_from_b(pc1, 10)
    ref = {
        "R0_b0": R0(pc0.with_(alpha=0.0)), "R0_b1": R0(pc1.with_(alpha=0.0)),
        "relapse_multiplier": relapse_multiplier(pc1),
        "eq_b1": list(equilibrium(pc1)), "eq_b0": list(equilibrium(pc0)),
        "eq_b1_alpha0": list(equilibrium(pc1.with_(alpha=0.0))),
        "floor": float(equilibrium(pc1.with_(b=0.0))[2]),
        "I_P": list(yP[2]), "I_T": list(yT[2]),
        "beta_k10": beta, "T_career_k10": nw.career_transmissibility(pc1.with_(alpha=0.0), beta),
        "T_single_k10": nw.transmissibility(pc1, beta),
    }
    ref = {k: ([r6(v) for v in x] if isinstance(x, list) else r6(x)) for k, x in ref.items()}

    # ---- age and sex: central constant-vulnerability fit and the variant with a change in 2018
    import ext_paper5 as e5
    st = AS.Struct()
    asx = json.load(open(OUT / "results_agesex.json")); asxs = json.load(open(OUT / "results_asxs.json"))
    v = lambda d, k: d[k]["value"]
    base = pc1.with_(b=0.0)
    stc = st.with_(phi=v(asx, "asx_c_P_phi"), psi=v(asx, "asx_c_P_psi"))
    p0 = base.with_(b=v(asx, "asx_c_P_b0"))
    tt = np.arange(0, 24.001, 0.5)
    _, yc = AS.simulate(p0.with_(b=v(asx, "asx_c_P_k2")), stc, AS.equilibrium(p0, stc), 24, t_eval=tt)
    gc, bc = e5.series(yc)
    sol = {k: v(asxs, f"asxs_c_{k}") for k in ("b0", "b1", "phi0", "phi1", "psi")}
    q0 = base.with_(b=sol["b0"])
    _, ys = AS.simulate(q0.with_(b=sol["b1"]), st.with_(phi=sol["phi1"], psi=sol["psi"]),
                        AS.equilibrium(q0, st.with_(phi=sol["phi0"], psi=sol["psi"])), 24, t_eval=tt)
    gs, bs = e5.series(ys)
    agesex = {"years": [r6(2018 + x) for x in tt], "const": {"girls": [r6(x) for x in gc], "boys": [r6(x) for x in bc]},
              "shift": {"girls": [r6(x) for x in gs], "boys": [r6(x) for x in bs]}}

    # ---- parameter-recovery design table
    rec = [dict(noise=r["noise"], n=int(r["n"]), transition=r["transition"], true=r6(float(r["true"])),
                bias=r6(float(r["bias"])), rmse=r6(float(r["rmse"]))) for r in csv.DictReader(open(OUT / "recovery.csv"))]

    # ---- verified references used on the lab's screens
    want = ["BonielNissim2024", "vandenEijnden2016", "Boer2020", "Boer2022", "Hill2010", "Manski1993", "Shalizi2011",
            "Orben2024", "Ko2007", "Ko2009", "Wartberg2019", "Wartberg2020", "Xin2026", "Chiang2019", "Malygin2025",
            "Jo2025", "vandenDriessche2002", "CastilloChavez2004", "Newman2002", "Bliss1939", "Lenhart2007",
            "Zucchini2016", "StatisticsSouthAfrica2025", "Hethcote2000", "Kermack1927", "Poole2000", "Mossong2008",
            "Cohen2003", "Centola2007", "MarinoHBSCItalia2020"]
    refs = {}
    if RESEARCH.exists():
        for f in [RESEARCH / "references.json"] + sorted(RESEARCH.glob("references_*.json")):
            if f.name == "references_library.json":
                continue
            for e in json.load(open(f)):
                if e["key"] in want:
                    refs[e["key"]] = e["apa"]
    missing = [k for k in want if k not in refs]
    if missing:
        raise SystemExit(f"references not found: {missing}")

    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, capture_output=True, text=True).stdout.strip()
    LAB = {
        "meta": {"generated": datetime.date.today().isoformat(), "commit": commit,
                 "note": "Generated by python/export_lab.py from the thesis analysis outputs. Do not edit by hand."},
        "hbsc": {"2018": 0.07, "2022": 0.11, "girls2022": 0.13, "boys2022": 0.09, "age11_2022": 0.09, "age13_2022": 0.12},
        "ranges": {k: {kk: vv for kk, vv in ranges[k].items()} for k in cal.PARAM_ORDER},
        "central": {"params": {k: r6(getattr(pc1, k)) for k in ["eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"]},
                    "b0": r6(pc0.b), "b1": r6(pc1.b), "d_alpha": r6(da), "x18": [r6(x) for x in x18]},
        "interventions": {   # mirrors run_all.INTERVENTIONS (multipliers on rates; efficacy drawn uniformly)
            "platform": {"label": "Platform design regulation", "acts": "b", "range": [0.10, 0.40],
                         "desc": "Exposure rate b reduced by the efficacy."},
            "school": {"label": "School-based prevention", "acts": "gamma1", "range": [0.10, 0.50],
                       "desc": "Remission from At-Risk (γ₁) increased by the efficacy."},
            "minimum_age": {"label": "Minimum age of 16", "acts": "b, alpha", "range": [0.30, 0.70],
                            "desc": "Compliance c; exposure b and non-social acquisition α each reduced by 0.4c (the share of 10-to-24-year-olds under 16 times c)."},
            "treatment": {"label": "Treatment access", "acts": "gamma2, delta", "range": [0.20, 1.00], "range2": [0.0, 0.30],
                          "desc": "Recovery γ₂ increased by the first efficacy and relapse δ reduced by the second."},
            "start": ra.START},
        "draw_keys": keys + ["x18_S", "x18_A", "x18_I", "x18_R"],
        "draws": D,
        "projections": {"years": list(range(2018, 2043)), **proj},
        "results": {k: (r6(x) if isinstance(x, (int, float)) else x) for k, x in R.items() if KEEP.search(k)},
        "recovery": rec,
        "agesex": agesex,
        "ref": ref,
        "references": refs,
    }
    APP.write_text("// lab-data.js: generated by python/export_lab.py from the thesis analysis outputs. Do not edit.\n"
                   "export const LAB = " + json.dumps(LAB, ensure_ascii=False, separators=(",", ":")) + ";\n")
    print(f"wrote {APP} ({APP.stat().st_size / 1024:.0f} kB), {len(D)} draws, {len(LAB['results'])} result keys")


if __name__ == "__main__":
    main()
