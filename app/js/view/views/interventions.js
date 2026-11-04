// interventions.js: Screen 6, intervention mechanisms under assumed efficacies (thesis Chapter 9).

import { h } from "../dom.js";
import { panel, metricGrid, chartCard, legend, sliderRow, condNote } from "../ui.js";
import { timeSeries } from "../charts.js";
import { viewHead, table, source, pc, COLORS } from "../common.js";
import { R, CENTRAL, DRAWS } from "../../model/params.js";
import { SPEC, KEYS, START, END, multipliers, runDraw, acrossDraws, applyAll } from "../../model/interventions.js";
import { floor } from "../../model/sair.js";

const LABEL = { platform: "Platform design", school: "School prevention", minimum_age: "Minimum age 16", treatment: "Treatment access" };

export function InterventionsView() {
  const el = h("div", { class: "view" });
  el.appendChild(viewHead("Screen 6 · Interventions", "Which mechanisms have most leverage?",
    "Four mechanisms are represented as changes in rates from 2024: platform design regulation lowers the exposure rate; school-based prevention speeds remission from risk; a minimum age of 16 lowers exposure and non-social acquisition for the covered share; treatment access speeds recovery and reduces relapse. " +
    "The efficacies are <strong>assumptions</strong>. The results describe relative leverage inside the model, not what any programme would achieve, and no costs are modelled."));

  const centralDraw = { rates: CENTRAL.rates, cal: { b0: CENTRAL.b0, b1: CENTRAL.b1, dAlpha: CENTRAL.dAlpha, x18: CENTRAL.x18 } };
  const on = { platform: true, school: false, minimum_age: false, treatment: true };
  const eff = { platform: 0.25, school: 0.3, minimum_age: 0.5, treatment: 0.6 };
  let eff2 = 0.15, scenario = "P";
  const mults = () => KEYS.filter((k) => on[k]).map((k) => multipliers(k, eff[k], k === "treatment" ? eff2 : 0));
  const yrs = Array.from({ length: END - START + 1 }, (_, i) => START + i);

  const metrics = h("div");
  const chart = chartCard({ title: "Problematic use from 2024", sub: "central parameters (live)", filename: "interventions",
    render: () => {
      const base = runDraw(centralDraw, scenario, []), w = runDraw(centralDraw, scenario, mults());
      const p = scenario === "P" ? { ...CENTRAL.rates, b: CENTRAL.b1 } : { ...CENTRAL.rates, b: CENTRAL.b0 };
      const fl = floor(applyAll(p, mults()));
      metrics.replaceChildren(metricGrid([
        { label: "2042, no intervention", value: pc(base.at(-1)), foot: scenario === "P" ? "persistent shift" : "temporary shock" },
        { label: "2042, with selected", value: pc(w.at(-1)), foot: `${((1 - w.at(-1) / base.at(-1)) * 100).toFixed(0)}% lower`, accent: true },
        { label: "Floor with these rates", value: pc(fl, 2), foot: "exposure measures cannot pass it" },
      ]));
      return timeSeries({ xs: yrs, yTitle: "Problematic use", series: [
        { key: "b", name: "No intervention", color: "var(--ink-2)", ys: base },
        { key: "w", name: "With selected", color: COLORS.I, ys: w, width: 2.6 }],
        markers: [] });
    } });

  const toggles = h("div", { class: "toggle-list" }, KEYS.map((k) => {
    const sp = SPEC[k];
    const box = h("div", { class: "toggle", dataset: { on: String(on[k]) } }, [h("span", { class: "box" }),
      h("div", {}, [h("div", { class: "tl-title", text: sp.label }), h("div", { class: "tl-blurb", text: sp.desc })])]);
    box.addEventListener("click", () => { on[k] = !on[k]; box.dataset.on = String(on[k]); chart.update(); });
    const sl = sliderRow({ sym: "e", label: k === "minimum_age" ? "Compliance" : "Efficacy", unit: `thesis range ${sp.range[0]} to ${sp.range[1]}`, desc: "", min: sp.range[0], max: sp.range[1], step: 0.01, value: eff[k], baseline: eff[k],
      onInput: (v) => { eff[k] = v; chart.update(); } });
    const extra = k === "treatment" ? sliderRow({ sym: "e₂", label: "Relapse reduction", unit: `thesis range ${sp.range2[0]} to ${sp.range2[1]}`, desc: "", min: sp.range2[0], max: sp.range2[1], step: 0.01, value: eff2, baseline: eff2,
      onInput: (v) => { eff2 = v; chart.update(); } }) : null;
    return h("div", {}, [box, sl, extra]);
  }));
  const scen = h("div", { class: "row", style: { gap: "8px", marginBottom: "10px" } }, ["P", "T"].map((s) =>
    h("button", { class: "sm", text: s === "P" ? "Persistent shift" : "Temporary shock", onclick: () => { scenario = s; chart.update(); } })));

  el.appendChild(h("div", { class: "split" }, [
    panel({ title: "Live comparison", children: [scen, metrics, chart, legend([{ name: "No intervention", color: "var(--ink-2)" }, { name: "With selected mechanisms", color: COLORS.I }])] }),
    panel({ title: "Mechanisms", children: [toggles] }),
  ]));

  // ---- across all draws (live) and thesis results
  const acrossOut = h("div");
  const runAcross = () => {
    const combo = KEYS.filter((k) => on[k]);
    if (!combo.length) { acrossOut.replaceChildren(h("p", { class: "hint", text: "Select at least one mechanism." })); return; }
    acrossOut.replaceChildren(h("p", { class: "hint", text: "computing…" }));
    setTimeout(() => {
      const a = acrossDraws(DRAWS, combo, "P"), b = acrossDraws(DRAWS, combo, "T");
      const name = combo.join("+");
      const thesis = (sc) => { try { const pre = sc === "P" ? "int_" : "int_T_"; return `${R(`${pre}${name}_median`).toFixed(0)}% (${R(`${pre}${name}_p05`).toFixed(0)}% to ${R(`${pre}${name}_p95`).toFixed(0)}%)`; } catch { return "not run in the thesis"; } };
      acrossOut.replaceChildren(table(["Reduction in 2042", "Lab, live", "Thesis"], [
        ["Persistent shift", `${(a.median * 100).toFixed(0)}% (${(a.p05 * 100).toFixed(0)}% to ${(a.p95 * 100).toFixed(0)}%)`, thesis("P")],
        ["Temporary shock", `${(b.median * 100).toFixed(0)}% (${(b.p05 * 100).toFixed(0)}% to ${(b.p95 * 100).toFixed(0)}%)`, thesis("T")],
      ], ["l", "r", "r"]), h("p", { class: "hint", text: `Median and 5th to 95th percentile across ${a.n} parameter draws, with each efficacy drawn uniformly from its assumed range. The lab draws its own efficacies, so its figures differ slightly from the thesis's.` }));
    }, 20);
  };

  const single = KEYS.map((k) => [LABEL[k], `${R(`int_${k}_median`).toFixed(0)}% (${R(`int_${k}_p05`).toFixed(0)}% to ${R(`int_${k}_p95`).toFixed(0)}%)`, `${R(`int_T_${k}_median`).toFixed(0)}%`, ""]);
  const combos = ["platform+minimum_age", "platform+school", "school+treatment", "platform+school+minimum_age+treatment"].map((c) =>
    [c.split("+").map((k) => LABEL[k]).join(" + "), `${R(`int_${c}_median`).toFixed(0)}% (${R(`int_${c}_p05`).toFixed(0)}% to ${R(`int_${c}_p95`).toFixed(0)}%)`, `${R(`int_T_${c}_median`).toFixed(0)}%`,
      `${R(`int_${c}_bliss_index_median`).toFixed(2)} (${R(`int_${c}_bliss_index_p05`).toFixed(2)} to ${R(`int_${c}_bliss_index_p95`).toFixed(2)})`]);
  el.appendChild(panel({ title: "Across all calibrated draws", children: [
    h("button", { class: "primary sm", text: "Run the selected combination across all draws", onclick: runAcross }), acrossOut,
    h("h3", { text: "Thesis results: reduction in 2042 prevalence", style: { marginTop: "18px" } }),
    table(["Mechanisms", "Persistent shift, median (range)", "Temporary shock, median", "Bliss index"], single.concat(combos), ["l", "r", "r", "r"]),
    h("p", { class: "prose", html: "A Bliss index above 1 means the combination does better than independent action would predict. Combinations show at most mild synergy. Treatment access has the largest effect partly because of the assumed efficacy ranges. Measures on exposure cannot pass the floor; measures on recovery are not bounded in the same way but show diminishing returns." }),
  ] }));

  el.appendChild(panel({ title: "Dividing effort over time: optimal control", children: [
    metricGrid([
      { label: "Value of recovery effort", value: `${R("ctl_mv_ratio_per_unit_rate").toFixed(1)}×`, foot: "exposure effort, per unit change in the rate", accent: true },
      { label: "Within the efficacy ceilings", value: `${R("ctl_mv_ratio").toFixed(1)}×`, foot: "marginal value ratio" },
      { label: "Benefit from timing", value: `${R("ctl_ref_gain_pct").toFixed(1)}%`, foot: "optimal path against constant effort" },
    ]),
    h("p", { class: "prose", html: "An optimal-control problem with one control on exposure and one on recovery has an optimal solution (Fleming and Rishel) whose necessary conditions follow from the maximum principle; it was solved by a forward-backward sweep. With costs entering only as assumed relative weights, how effort is divided between mechanisms matters far more than when it is applied." }),
    condNote("All intervention effects are conditional on assumed efficacies and on the model. They are statements about relative leverage, not predictions for any programme."),
    source("thesis Chapter 9 and the interventions paper (Mutibura and Lubuma, manuscript in preparation); python/run_all.py, python/ext_interventions.py and python/ext_control.py."),
  ] }));
  return { el };
}
