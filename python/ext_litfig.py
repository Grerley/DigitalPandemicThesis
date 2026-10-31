"""Conceptual causal diagram (DAG) for Chapter 3 (literature workstream).

Draws the hypothesised mechanisms behind the two acquisition routes of the
S-A-I-R model with Graphviz (`dot` must be on the PATH) and writes
figures/fig_litx_dag.png. This is a conceptual figure, not model output: it
uses no data and no results file.

Node conventions
  filled box      observed in the data the thesis uses (HBSC published aggregates)
  plain box       observable with other designs (panels, network surveys, logs)
  dashed ellipse  latent in every design considered
  double border   a friendship tie: analyses of peer effects condition on it
Edge conventions (colour, line style and label all carry the meaning)
  blue, solid     influence: the only route that b is meant to represent
  orange, dashed  shared-environment and homophily paths that also load onto b
  green, solid    non-social routes that load onto alpha
  red, dotted     reverse causation (absorbed into alpha or sigma)
  grey            other hypothesised dependencies
"""

import os
import pathlib
import subprocess

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "figures" / "fig_litx_dag.png"

OBS_DATA = 'shape=box, style="filled,rounded", fillcolor="#d9d9d9", color="#333333"'
OBS_OTHER = 'shape=box, style="rounded", color="#333333"'
LATENT = 'shape=ellipse, style="dashed", color="#333333"'
TIE = 'shape=box, peripheries=2, style="rounded", color="#333333"'

INFL = 'color="#1f4e9c", penwidth=2.2, fontcolor="#1f4e9c"'
BUNDLED = 'color="#c26a00", penwidth=1.8, style=dashed, fontcolor="#c26a00"'
NONSOC = 'color="#2e7d32", penwidth=2.0, fontcolor="#2e7d32"'
REVERSE = 'color="#b00020", penwidth=1.8, style=dotted, fontcolor="#b00020"'
OTHER = 'color="#7a7a7a", penwidth=1.0, fontcolor="#555555"'

BW = os.environ.get("DPM_BW") == "1"   # black-and-white edition for the LaTeX thesis
TIE_INFL = 'penwidth=1.2, style=dashed'
KEY_EDGES = """      <TR><TD ALIGN="LEFT"><FONT COLOR="#1f4e9c">blue: influence, what b is meant to represent</FONT></TD></TR>
      <TR><TD ALIGN="LEFT"><FONT COLOR="#c26a00">orange dashed: shared environment, homophily (also load on b)</FONT></TD></TR>
      <TR><TD ALIGN="LEFT"><FONT COLOR="#2e7d32">green: non-social routes (load on α)</FONT></TD></TR>
      <TR><TD ALIGN="LEFT"><FONT COLOR="#b00020">red dotted: reverse causation</FONT></TD></TR>"""
FONT = "Helvetica"
if BW:
    FONT = "Latin Modern Roman,"
    OBS_DATA = 'shape=box, style="filled,rounded", fillcolor="#d9d9d9", color="black"'
    OBS_OTHER = 'shape=box, style="rounded", color="black"'
    LATENT = 'shape=ellipse, style="dashed", color="black"'
    TIE = 'shape=box, peripheries=2, style="rounded", color="black"'
    INFL = 'color="black", penwidth=2.4, fontcolor="black"'
    BUNDLED = 'color="black", penwidth=1.3, style=dashed, fontcolor="black"'
    NONSOC = 'color="black", penwidth=1.3, arrowhead=empty, fontcolor="black"'
    REVERSE = 'color="black", penwidth=1.6, style=dotted, fontcolor="black"'
    OTHER = 'color="#808080", penwidth=0.8, fontcolor="#404040"'
    TIE_INFL = 'penwidth=1.0'
    KEY_EDGES = """      <TR><TD ALIGN="LEFT">thick line: influence, what b is meant to represent</TD></TR>
      <TR><TD ALIGN="LEFT">dashed line: shared environment, homophily (also load on b)</TD></TR>
      <TR><TD ALIGN="LEFT">open arrowhead: non-social routes (load on α)</TD></TR>
      <TR><TD ALIGN="LEFT">dotted line: reverse causation</TD></TR>
      <TR><TD ALIGN="LEFT">thin grey line: other dependencies and model transitions</TD></TR>"""

DOT = f"""
digraph DAG {{
  graph [rankdir=LR, splines=true, nodesep=0.35, ranksep=0.75, fontname="{FONT}",
         fontsize=11, pad=0.3, dpi=200];
  node  [fontname="{FONT}", fontsize=11, margin="0.12,0.06"];
  edge  [fontname="{FONT}", fontsize=9, arrowsize=0.7];

  // exogenous characteristics
  age  [label="Age", {OBS_DATA}];
  sex  [label="Sex", {OBS_DATA}];

  // latent common causes
  env  [label="Shared environment\\n(school, cohort, period,\\nplatform diffusion)", {LATENT}];
  vi   [label="Vulnerability of i\\n(stable traits,\\nearlier symptoms)", {LATENT}];
  vj   [label="Vulnerability of\\nfriend j", {LATENT}];

  // observable with other designs
  plat [label="Access and platform\\nexposure of i\\n(features, feeds)", {OBS_OTHER}];
  tie  [label="Friendship tie\\ni - j", {TIE}];
  ij   [label="Friend j's\\nproblematic use", {OBS_OTHER}];
  m0   [label="Mental health of i\\n(earlier)", {OBS_OTHER}];
  m1   [label="Mental health of i\\n(later)", {OBS_OTHER}];

  // model states of i
  ai   [label="Sub-threshold\\nsymptoms of i\\n(state A)", {LATENT}];
  ii   [label="Problematic use of i\\n(state I, SMDS 6+)", {OBS_DATA}];
  ri   [label="Recovered history\\nof i (state R)", {LATENT}];

  {{ rank=same; env; vi; vj; }}
  {{ rank=same; plat; tie; m0; }}
  {{ rank=same; ij; }}

  // exogenous characteristics modify several routes
  age -> vi   [{OTHER}];
  sex -> vi   [{OTHER}];
  age -> plat [{OTHER}];
  sex -> plat [{OTHER}];
  age -> tie  [label="assortative\\nmixing", {OTHER}];
  sex -> tie  [{OTHER}];

  // the route b is meant to represent
  ij -> ai    [label="influence\\n(b)", {INFL}];

  tie -> ai   [label="tie carries\\ninfluence", {INFL}, {TIE_INFL}];

  // paths that b absorbs in a mean-field model
  env -> ij   [label="shared\\nenvironment", {BUNDLED}];
  env -> ai   [label="shared\\nenvironment\\n(loads on b)", {BUNDLED}];
  vi -> tie   [label="homophily", {BUNDLED}];
  vj -> tie   [{BUNDLED}];
  vj -> ij    [{BUNDLED}];

  // non-social routes (α)
  env -> plat [{OTHER}];
  plat -> ai  [label="platform\\nroute (α)", {NONSOC}];
  vi -> ai    [label="individual\\nroute (α)", {NONSOC}];

  // reverse causation
  vi -> m0    [{OTHER}];
  m0 -> ai    [label="reverse\\ncausation\\n(α, σ)", {REVERSE}];

  // model transitions and the outcome association
  ai -> ii    [label="σ", {OTHER}];
  ii -> ri    [label="γ₂", {OTHER}];
  ri -> ai    [label="δ", {OTHER}, constraint=false];
  ii -> m1    [label="association\\n(not modelled)", {OTHER}];

  // compact key, placed with the exogenous characteristics
  key [shape=plaintext, fontsize=9, label=<
    <TABLE BORDER="1" COLOR="#999999" CELLBORDER="0" CELLSPACING="1" CELLPADDING="2">
      <TR><TD ALIGN="LEFT"><B>Key</B></TD></TR>
      <TR><TD ALIGN="LEFT" BGCOLOR="#d9d9d9">grey box: observed in the data used (HBSC aggregates)</TD></TR>
      <TR><TD ALIGN="LEFT">white box: observable with other designs</TD></TR>
      <TR><TD ALIGN="LEFT">dashed ellipse: latent</TD></TR>
      <TR><TD ALIGN="LEFT">double box: tie, a collider that peer-effect analyses condition on</TD></TR>
{KEY_EDGES}
    </TABLE>>];
  {{ rank=same; age; sex; key; }}
}}
"""


def main() -> None:
    if BW:
        d = pathlib.Path(os.environ.get("DPM_BW_DIR", HERE / "figures_bw")); d.mkdir(parents=True, exist_ok=True)
        dot_path = d / "fig_litx_dag.dot"; dot_path.write_text(DOT)
        for fmt in ("pdf", "png"):
            subprocess.run(["dot", f"-T{fmt}", str(dot_path), "-o", str(d / f"fig_litx_dag.{fmt}")], check=True)
        dot_path.unlink(); print("wrote", d); return
    OUT.parent.mkdir(parents=True, exist_ok=True)
    dot_path = OUT.with_suffix(".dot")
    dot_path.write_text(DOT)
    subprocess.run(["dot", "-Tpng", str(dot_path), "-o", str(OUT)], check=True)
    dot_path.unlink()
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
