"""Numbers and figures specific to Paper 1 (mathematical analysis of the S-A-I-R model).

Everything here is evaluated at the calibrated central parameter set after 2018 (Chapter 7 of the thesis),
with no re-calibration, so that the comparisons isolate the mathematical effect of each change.

1. Erlang spell durations at fixed parameters: R0, equilibrium prevalence, floor and the slowest decay rate,
   for shapes 1 to 20 and the deterministic limit (Proposition on general durations).
2. The state-transition diagram of the model (Graphviz), in the black-and-white journal style.

Writes outputs/results_p1.json (keys p1_...) and figures/fig_p1_model.png (figures_bw/fig_p1_model.pdf).
"""
from __future__ import annotations

import json, math, os, pathlib, subprocess
import numpy as np

import run_all as ra
from dpm.model import Params, R0, equilibrium
import ext_durations as ed

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"; FIG = HERE / "figures"; FIGBW = HERE / "figures_bw"
R: dict = {}


def put(k, v, note=""):
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    R["p1_" + k] = {"value": v, "note": note}


def slowest_rate(p, k):
    ch = ed.Chain(*k)
    J = ed.jac_num(ch.rhs(p), ch.eq_state(p))
    T = np.vstack([-np.ones((1, ch.n - 1)), np.eye(ch.n - 1)])
    ev = np.linalg.eigvals((J @ T)[1:, :])
    return float(-np.max(ev.real))


def erlang_fixed(pc1):
    shapes = [1, 2, 3, 5, 10, 20]
    for s in shapes:
        k = (s, s, s)
        put(f"erl{s}_R0", ed.R0_dur(pc1.with_(alpha=0.0), *k))
        put(f"erl{s}_mult", 1 / (1 - ed.stage_quantities(pc1, *k)["q"]))
        put(f"erl{s}_Ieq_pct", 100 * ed.I_eq_dur(pc1, *k))
        put(f"erl{s}_floor_pct", 100 * ed.I_eq_dur(pc1.with_(b=0.0), *k))
        r = slowest_rate(pc1, k)
        put(f"erl{s}_halflife", math.log(2) / r)
    big = (5000, 5000, 5000)
    put("erlinf_R0", ed.R0_dur(pc1.with_(alpha=0.0), *big))
    put("erlinf_mult", 1 / (1 - ed.stage_quantities(pc1, *big)["q"]))
    put("erlinf_Ieq_pct", 100 * ed.I_eq_dur(pc1, *big))
    put("erlinf_floor_pct", 100 * ed.I_eq_dur(pc1.with_(b=0.0), *big))
    for nm, k in [("A", (20, 1, 1)), ("I", (1, 20, 1)), ("R", (1, 1, 20))]:
        put(f"only{nm}20_R0", ed.R0_dur(pc1.with_(alpha=0.0), *k))
        put(f"only{nm}20_Ieq_pct", 100 * ed.I_eq_dur(pc1, *k))
    put("shapes", shapes)
    # checks: shape 1 equals the main model
    put("chk_R0_k1", abs(ed.R0_dur(pc1.with_(alpha=0.0), 1, 1, 1) - R0(pc1.with_(alpha=0.0))))
    put("chk_I_k1", abs(ed.I_eq_dur(pc1, 1, 1, 1) - equilibrium(pc1)[2]))


def trajectory_and_noise(pc0, pc1, x18):
    """Overshoot of the illustrative path after 2018, damping of the complex pair, binomial comparison for the LNA."""
    from dpm.model import simulate
    tt = np.linspace(0, 40, 40001)
    _, X = simulate(pc1, x18, 40, t_eval=tt)
    I = X[2]; S = X[0]; k = int(np.argmax(I))
    put("peak_pct", 100 * I[k]); put("peak_t", tt[k]); put("peak_year", 2018 + tt[k])
    put("I2022_pct", 100 * np.interp(4.0, tt, I))
    put("S2018_pct", 100 * S[0]); put("Speak_pct", 100 * S[k]); put("Sinf_pct", 100 * S[-1])
    put("A2018_pct", 100 * X[1][0]); put("Apeak_pct", 100 * X[1][k])
    put("R2018_pct", 100 * X[3][0]); put("Rpeak_pct", 100 * X[3][k]); put("Rinf_pct", 100 * X[3][-1])
    put("Ieq_pct", 100 * I[-1])
    lam = complex(ra_json("mathx_eig_b1_re")[0], ra_json("mathx_eig_b1_im")[0])
    put("damping_ratio", -lam.real / abs(lam)); put("period_years", 2 * math.pi / abs(lam.imag))
    pI = equilibrium(pc1)[2]
    for N in (1000, 10000):
        put(f"binom_sd_N{N}_pp", 100 * math.sqrt(pI * (1 - pI) / N))
    put("lna_over_binom_N1000", ra_json("mathx_lna_N1000_sd_pp") / (100 * math.sqrt(pI * (1 - pI) / 1000)))
    put("lna_excess_N1000_pct", 100 * (ra_json("mathx_lna_N1000_sd_pp") / (100 * math.sqrt(pI * (1 - pI) / 1000)) - 1))
    put("diff_sd_N1000_pp", math.sqrt(2) * ra_json("mathx_lna_N1000_sd_pp"))
    put("diff_mad_N1000_pp", math.sqrt(2) * ra_json("mathx_lna_N1000_sd_pp") * math.sqrt(2 / math.pi))


def ra_json(key):
    for f in sorted(OUT.glob("results*.json")):
        d = json.load(open(f))
        if key in d:
            return d[key]["value"]
    raise KeyError(key)


DOT = r"""
digraph SAIR {
  graph [rankdir=LR, nodesep=0.6, ranksep=1.0, pad=0.15];
  node  [shape=box, style="rounded", width=0.75, height=0.55, fontsize=18, fontname="Latin Modern Math,", color="black"];
  edge  [fontsize=15, fontname="Latin Modern Math,", color="black", arrowsize=0.7];
  in [shape=point, width=0.01, style=invis];
  S [label="𝑆"]; A [label="𝐴"]; I [label="𝐼", style="rounded,filled", fillcolor="#d9d9d9"]; R [label="𝑅"];
  o1 [shape=point, width=0.01, style=invis]; o2 [shape=point, width=0.01, style=invis];
  o3 [shape=point, width=0.01, style=invis]; o4 [shape=point, width=0.01, style=invis];
  in -> S [label="𝜈"];
  S -> A [label="𝜆 + 𝛼", penwidth=1.6];
  A -> S [label=<𝛾<SUB>1</SUB>>, constraint=false];
  A -> I [label="𝜎"];
  I -> R [label=<𝛾<SUB>2</SUB>>];
  R -> A [label="𝛿", style=dashed, constraint=false];
  S -> o1 [label="𝜈", color="#808080", fontcolor="#404040"];
  A -> o2 [label="𝜈", color="#808080", fontcolor="#404040"];
  I -> o3 [label="𝜈", color="#808080", fontcolor="#404040"];
  R -> o4 [label="𝜈", color="#808080", fontcolor="#404040"];
  {rank=same; S; o1;} {rank=same; A; o2;} {rank=same; I; o3;} {rank=same; R; o4;}
}
"""


def model_diagram():
    FIG.mkdir(exist_ok=True); FIGBW.mkdir(exist_ok=True)
    dot = FIGBW / "fig_p1_model.dot"; dot.write_text(DOT)
    subprocess.run(["dot", "-Tpdf", str(dot), "-o", str(FIGBW / "fig_p1_model.pdf")], check=True)
    for d in (FIG, FIGBW):
        subprocess.run(["dot", "-Tpng", "-Gdpi=200", str(dot), "-o", str(d / "fig_p1_model.png")], check=True)
    dot.unlink()


def main():
    pc0, pc1, da, x18 = ra.central(json.load(open(HERE / "data" / "parameter_ranges.json")))
    erlang_fixed(pc1)
    trajectory_and_noise(pc0, pc1, x18)
    model_diagram()
    json.dump(R, open(OUT / "results_p1.json", "w"), indent=1, ensure_ascii=False)
    print("done", len(R))


if __name__ == "__main__":
    main()
