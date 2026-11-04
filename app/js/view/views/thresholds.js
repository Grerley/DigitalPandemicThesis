// thresholds.js: Screen 4, the threshold, the floor and why the non-social route matters (thesis Chapter 5).

import { h } from "../dom.js";
import { panel, chartCard, sliderRow } from "../ui.js";
import { xyLines } from "../charts.js";
import { viewHead, pc, source } from "../common.js";
import { centralParams } from "../../model/params.js";
import { equilibrium, floor, blocks, elasticities } from "../../model/sair.js";

export function ThresholdsView() {
  const el = h("div", { class: "view" });
  el.appendChild(viewHead("Screen 4 · Thresholds and floor", "With a non-social route, there is no threshold to cross",
    "Without non-social acquisition (<var>α</var> = 0), problematic use dies out when R₀ ≤ 1 and persists when R₀ &gt; 1. " +
    "With any <var>α</var> &gt; 0 the problem-free state is no longer an equilibrium: there is exactly one equilibrium for every exposure rate, it rises with exposure and with <var>α</var>, " +
    "and as exposure falls to zero prevalence falls only to a positive <strong>floor</strong>. Under the model, no policy that acts on social exposure alone can go below it."));

  let alpha = centralParams().alpha;
  const bs = Array.from({ length: 161 }, (_, i) => (i * 6) / 160);
  const bif = chartCard({ title: "Equilibrium prevalence against the exposure rate", sub: "central values of the other rates; computed live", filename: "threshold",
    render: () => {
      const p = centralParams(), { c, d, D } = blocks(p), bCrit = D / (d * (p.eta * c + p.sigma));
      return xyLines({ xTitle: "Social exposure rate b (per year)", yTitle: "Equilibrium problematic use I*", xFmt: (v) => v.toFixed(1), yFmt: (v) => `${(v * 100).toFixed(0)}%`,
        series: [
          { name: "α = 0", color: "var(--ink-2)", dash: "5 4", ys: bs.map((b) => equilibrium({ ...p, b, alpha: 0 })[2]), xs: bs },
          { name: `α = ${alpha.toFixed(3)}`, color: "var(--c-i)", ys: bs.map((b) => equilibrium({ ...p, b, alpha })[2]), xs: bs },
        ],
        markers: [{ x: bCrit, label: "R₀ = 1", color: "var(--ink-3)" }, { y: floor({ ...p, alpha }), label: `floor ${pc(floor({ ...p, alpha }), 2)}`, color: "var(--c-i)" }] });
    } });
  const aSlider = sliderRow({ sym: "α", label: "Non-social acquisition", unit: "per year · thesis range 0.002 to 0.05", desc: "Move it to zero to recover the classical threshold.",
    min: 0, max: 0.1, step: 0.001, value: alpha, baseline: alpha, onInput: (v) => { alpha = v; bif.update(); fl.update(); } });

  const as = Array.from({ length: 60 }, (_, i) => 0.001 * Math.pow(10, (i / 59) * 2));
  const fl = chartCard({ title: "The floor against non-social acquisition", sub: "for three recovery rates", filename: "floor",
    render: () => {
      const p = centralParams();
      return xyLines({ w: 560, h: 420, xTitle: "Non-social acquisition α (per year, log scale)", yTitle: "Floor I*min", xLog: true, xFmt: (v) => (v < 0.01 ? v.toFixed(3) : v.toFixed(2)), yFmt: (v) => `${(v * 100).toFixed(1)}%`,
        series: [0.3, 0.55, 1.0].map((g, i) => ({ name: `γ₂ = ${g}`, color: ["var(--c-i)", "var(--accent)", "var(--c-s)"][i], dash: ["", "6 3", "2 3"][i], xs: as, ys: as.map((a) => floor({ ...p, alpha: a, gamma2: g })) })),
        markers: [{ x: alpha > 0.001 ? alpha : 0.0011, label: "current α", color: "var(--ink-3)" }] });
    } });

  const bs2 = Array.from({ length: 80 }, (_, i) => 0.2 + (i * 5.8) / 79);
  const el2 = chartCard({ title: "How responsive prevalence is to exposure", sub: "index ∂ln I* / ∂ln b, central rates", filename: "exposure-elasticity",
    render: () => {
      const p = centralParams();
      return xyLines({ w: 560, h: 420, xTitle: "Social exposure rate b (per year)", yTitle: "Sensitivity index", xFmt: (v) => v.toFixed(1), yFmt: (v) => v.toFixed(1),
        series: [{ name: "∂ln I*/∂ln b", color: "var(--accent)", xs: bs2, ys: bs2.map((b) => elasticities({ ...p, b, alpha }, "I").b) }] });
    } });

  el.appendChild(panel({ title: "The threshold disappears", children: [aSlider, bif,
    h("p", { class: "prose", html: "The dashed curve (<var>α</var> = 0) is zero up to R₀ = 1 and then rises: the forward bifurcation of the classical model. The solid curve never reaches zero. Lowering exposure moves prevalence down the solid curve towards the floor, not to zero." })] }));
  el.appendChild(h("div", { class: "grid cols-2" }, [
    panel({ title: "What sets the floor", children: [fl, h("p", { class: "prose", html: "The floor rises with <var>α</var> and falls with faster recovery. Measures on recovery and relapse lower the floor itself; measures on exposure only move prevalence towards it." })] }),
    panel({ title: "Exposure matters more near the threshold", children: [el2, h("p", { class: "prose", html: "The index of prevalence with respect to exposure rises as exposure falls towards the threshold. This is why exposure measures show accelerating returns, and why combinations of exposure measures can show mild synergy." })] }),
  ]));
  el.appendChild(panel({ title: "The results behind this screen", children: [
    h("ul", { class: "prose" }, [
      h("li", { html: "R₀ = b d (η c + σ)/(a c d − σ γ₂ δ), the product of new At-Risk young people generated in one spell and the expected number of spells 1/(1 − q), q = (σ/a)(γ₂/c)(δ/d)." }),
      h("li", { html: "With <var>α</var> = 0 the problem-free state is globally asymptotically stable for R₀ ≤ 1 (Lyapunov function and LaSalle's principle) and the bifurcation at R₀ = 1 is forward (centre manifold theorem)." }),
      h("li", { html: "With <var>α</var> &gt; 0 the equilibrium is unique for every <var>b</var>, increasing in <var>b</var> and <var>α</var>, with the floor I*<sub>min</sub> = (σ/c)·α/(α m + ω) as <var>b</var> → 0, and it is locally asymptotically stable for every admissible parameter set (Routh–Hurwitz)." }),
    ]),
    source("thesis Chapter 5 and the paper on thresholds, floors and stability (Mutibura and Lubuma, manuscript in preparation)."),
  ] }));
  return { el };
}
