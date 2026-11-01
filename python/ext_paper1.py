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
    model_diagram()
    json.dump(R, open(OUT / "results_p1.json", "w"), indent=1, ensure_ascii=False)
    print("done", len(R))


if __name__ == "__main__":
    main()
