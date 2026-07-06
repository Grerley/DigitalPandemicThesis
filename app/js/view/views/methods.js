// methods.js — Screen 7: Methods & Provenance.
// -----------------------------------------------------------------------------

import { h } from "../dom.js";
import { panel } from "../ui.js";
import { selfTestPanel } from "./selftest-panel.js";
import { PARAM_META, PARAM_KEYS, DEFAULT_PARAMS, INITIAL_STATE } from "../../model/params.js";
import { fx } from "../format.js";

export function MethodsView() {
  const el = h("div", { class: "view" });

  el.appendChild(h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: "Screen 7 · Methods & Provenance" }),
    h("h1", { text: "Methods & Provenance" }),
    h("p", { class: "prose", text: "The model specification, equations, parameters and their sources, and the epistemic caveats carried from the thesis. Every number on the other screens traces back to what is written here." }),
  ]));

  // Self-test panel
  el.appendChild(panel({ title: "Correctness self-test", sub: "assertions run live against the model code", children: [selfTestPanel()] }));

  // Model spec (plain language)
  el.appendChild(panel({ title: "Model specification", children: [
    h("p", { class: "prose", html: "The population is described by a state distribution <var>x</var>(t) = [S, A, I, R]ᵀ of fractions summing to one, evolving on a monthly step under a <strong>state-dependent</strong> transition matrix <var>P</var>(x). The compartments are <strong>Susceptible</strong>, <strong>At-Risk</strong> (heavy/pre-clinical use, partially transmitting), <strong>Addicted</strong> (meets problematic-use thresholds), and <strong>Recovered</strong> (reduced use, <em>not</em> absorbing — subject to relapse)." }),
    eq(`<span class="frac"><span class="num">to</span></span>`, "matrix"),
    h("div", { class: "eq", html: transitionMatrixHTML() }),
    h("p", { class: "prose", html: "Each row is a valid probability vector (sums to 1). The relapse entry P<sub>RA</sub> (from δ) is what makes the epidemic self-sustaining; without it R would be absorbing and I would decay to zero." }),
  ] }));

  // Equations
  el.appendChild(panel({ title: "Equations", children: [
    h("h3", { text: "Mass-action force of infection (Eq 4.3)" }),
    h("div", { class: "eq", html: `P<sub>SA</sub>(x<sub>t</sub>) = min&#8202;(&#8202;1 − exp[ −β·k·(I<sub>t</sub> + η·A<sub>t</sub>) ] + μ,&#8202; 1&#8202;)` }),
    h("p", { class: "prose", html: "Because P<sub>SA</sub> depends on the current prevalence, the map is <strong>nonlinear</strong> — a fixed Markov chain would converge monotonically and could not generate a rise–peak–decline epidemic." }),
    h("h3", { text: "State update (Eq 4.2)" }),
    h("div", { class: "eq", html: `x(t+1) = P(x<sub>t</sub>)<sup>⊤</sup> · x(t)` }),
    h("h3", { text: "Reed–Frost network transmission (Eq 4.4)" }),
    h("div", { class: "eq", html: `P<sub>i</sub>(S → A | m) = min&#8202;(&#8202;1 − (1 − β)<sup>m</sup> + μ,&#8202; 1&#8202;),&nbsp;&nbsp; m = (addicted neighbours) + η·(at-risk neighbours)` }),
    h("h3", { text: "Basic reproduction number (next-generation matrix, Appendix A.4)" }),
    h("div", { class: "eq", html: `R₀ = ρ(F·V⁻¹) = <span class="frac"><span class="num">β·k·S₀·(η·γ₂ + σ)</span><span class="den">(σ + γ₁)·γ₂</span></span>` }),
    h("p", { class: "prose", html: "Evaluated on the {A, I} subsystem this gives the operative R₀ ≈ 2.50 at S₀ = 0.823 and basic R₀ ≈ 3.04 at S₀ = 1. Including the relapse compartment R yields the <em>full</em> R₀ that governs stability of the disease-free equilibrium and crosses 1 exactly at the transcritical bifurcation." }),
    h("h3", { text: "Endemic equilibrium & threshold (Appendix A.5–A.6)" }),
    h("div", { class: "eq", html: `S* = S₀ / R₀,&nbsp;&nbsp; I* = (σ / γ₂)·A*&nbsp;&nbsp; for R₀ > 1;&nbsp;&nbsp; transcritical bifurcation at R₀ = 1` }),
    h("h3", { text: "Continuous-time ODE analogue (Appendix A.1)" }),
    h("div", { class: "eq", html: `dS/dt = Λ − β·S·(I+ηA) + γ₁A − μS<br>dA/dt = β·S·(I+ηA) − (σ+γ₁+μ)A + δR<br>dI/dt = σA − (γ₂+μ)I<br>dR/dt = γ₂I − (δ+μ)R` }),
  ] }));

  // Parameter table
  el.appendChild(panel({ title: "Parameters & sources", sub: "calibrated monthly rates (Chapter 6, Table 6.4)", children: [
    h("table", { class: "data" }, [
      h("thead", {}, h("tr", {}, [h("th", { style: { textAlign: "left" }, text: "Symbol" }), h("th", { style: { textAlign: "left" }, text: "Parameter" }), h("th", { text: "Value" }), h("th", { text: "Unit" }), h("th", { style: { textAlign: "left" }, text: "Source" })])),
      h("tbody", {}, [
        ...PARAM_KEYS.map((k) => {
          const m = PARAM_META[k];
          return h("tr", {}, [
            h("td", { style: { textAlign: "left", fontFamily: "var(--serif)", fontStyle: "italic" }, text: m.sym }),
            h("td", { style: { textAlign: "left" }, text: m.label }),
            h("td", { text: fx(DEFAULT_PARAMS[k], DEFAULT_PARAMS[k] < 0.01 ? 4 : 3) }),
            h("td", { class: "hint", text: m.unit }),
            h("td", { class: "hint", style: { textAlign: "left" }, text: m.source }),
          ]);
        }),
        ...["S", "A", "I", "R"].map((s) => h("tr", {}, [
          h("td", { style: { textAlign: "left", fontFamily: "var(--serif)", fontStyle: "italic" }, text: `x₀(${s})` }),
          h("td", { style: { textAlign: "left" }, text: `Initial ${s}` }),
          h("td", { text: fx(INITIAL_STATE[s], 3) }), h("td", { class: "hint", text: "proportion" }),
          h("td", { class: "hint", style: { textAlign: "left" }, text: s === "R" ? "Assumed" : "HBSC 2022" }),
        ])),
      ]),
    ]),
  ] }));

  // Epistemic caveats
  el.appendChild(panel({ title: "Epistemic caveats", children: [
    h("div", { class: "caveat", html: "<strong>All quantitative outputs are conditional model results, not empirical forecasts.</strong> The robust conclusions are qualitative — a self-sustaining epidemic with a critical threshold and an exploitable network structure; specific numbers illustrate rather than measure." }),
    h("ul", { class: "prose" }, [
      h("li", { html: "Calibration relies partly on a <strong>synthetic cohort</strong> standing in for longitudinal transition data that do not yet exist." }),
      h("li", { html: "The transmission coefficient β captures peer influence, homophily and shared context jointly; it is <strong>not</strong> identified causal peer transmission." }),
      h("li", { html: "Intervention efficacies are <strong>assumed</strong>, informed by the direction of published evaluations, not estimated here — scenario outputs are comparative, not predictive." }),
      h("li", { html: "Risk-factor odds ratios are <strong>associational</strong>, not causal. Cost figures are illustrative unit costs, not a programme budget." }),
      h("li", { html: "The South African arm is extrapolated from OECD priors with adjustment factors, not independently estimated." }),
    ]),
    h("p", { class: "prose", style: { marginTop: "12px" }, html: "Source: <em>The Digital Pandemic: An Epidemic-Modelling Framework for Technology-Induced Mental Health in a Connected World</em> (Mutibura, 2026). This lab reimplements the thesis pipeline (R/Stan) in the browser; the deterministic core is verified against the acceptance targets above." }),
  ] }));

  return { el, destroy: () => {} };
}

function eq() { return h("span"); }

function transitionMatrixHTML() {
  return `<table style="border-collapse:collapse;font-family:var(--serif)">
    <tr><td></td><td style="padding:0 10px;color:var(--ink-3)">to S</td><td style="padding:0 10px;color:var(--ink-3)">to A</td><td style="padding:0 10px;color:var(--ink-3)">to I</td><td style="padding:0 10px;color:var(--ink-3)">to R</td></tr>
    <tr><td style="color:var(--ink-3)">from S</td><td>1−P<sub>SA</sub></td><td>P<sub>SA</sub></td><td>0</td><td>0</td></tr>
    <tr><td style="color:var(--ink-3)">from A</td><td>P<sub>AS</sub></td><td>1−P<sub>AS</sub>−P<sub>AI</sub></td><td>P<sub>AI</sub></td><td>0</td></tr>
    <tr><td style="color:var(--ink-3)">from I</td><td>0</td><td>0</td><td>1−P<sub>IR</sub></td><td>P<sub>IR</sub></td></tr>
    <tr><td style="color:var(--ink-3)">from R</td><td>0</td><td>P<sub>RA</sub></td><td>0</td><td>1−P<sub>RA</sub></td></tr>
  </table>`;
}
