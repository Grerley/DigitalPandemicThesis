// explanations.js: Screen 2, two explanations of the 2018 to 2022 rise and what the next HBSC wave can decide.

import { h } from "../dom.js";
import { panel, metricGrid, chartCard, legend, sliderRow } from "../ui.js";
import { timeSeries } from "../charts.js";
import { viewHead, table, source, pc, COLORS } from "../common.js";
import { LAB, R, DRAWS, PARAM_META } from "../../model/params.js";
import { scenarioPath } from "../../model/calibrate.js";
import { R0 } from "../../model/sair.js";

export function ExplanationsView() {
  const el = h("div", { class: "view" });
  el.appendChild(viewHead("Screen 2 · Two explanations", "Two explanations fit the data exactly",
    "Two published figures can identify at most two quantities. For each draw of the other rates from literature-informed ranges, two quantities are solved so that the model reproduces " +
    "7% in 2018 and 11% in 2022 exactly. A <strong>persistent shift</strong> raises the exposure rate from 2018 onwards. A <strong>temporary shock</strong> raises non-social acquisition during 2020 and 2021 only. Both fit; they disagree about the future."));

  const P = LAB.projections, yrs = P.years;
  const bandChart = chartCard({
    title: "Projected problematic use, 2018 to 2042", sub: `median and 5th to 95th percentile across ${R("cal_n_ok")} kept draws (thesis results)`,
    filename: "two-explanations",
    csv: () => yrs.map((y, i) => ({ year: y, P_p05: P.P[5][i], P_median: P.P[50][i], P_p95: P.P[95][i], T_p05: P.T[5][i], T_median: P.T[50][i], T_p95: P.T[95][i] })),
    render: () => timeSeries({ xs: yrs, yTitle: "Problematic use", yMax: 0.24,
      series: [
        { key: "P", name: "Persistent shift", color: COLORS.P, ys: P.P[50], band: { lo: P.P[5], hi: P.P[95] } },
        { key: "T", name: "Temporary shock", color: COLORS.T, dash: "6 4", ys: P.T[50], band: { lo: P.T[5], hi: P.T[95] } },
      ],
      points: [{ x: 2018, y: 0.07, label: "HBSC 7%", color: "var(--ink)" }, { x: 2022, y: 0.11, label: "HBSC 11%", color: "var(--ink)" }],
      markers: [{ x: 2026, label: "next HBSC wave", color: "var(--ink-3)" }] }),
  });
  el.appendChild(panel({ title: "The projections", children: [
    legend([{ name: "Persistent shift (median, band)", color: COLORS.P }, { name: "Temporary shock (median, band)", color: COLORS.T, dash: "6 4" }]),
    bandChart,
    metricGrid([
      { label: "Persistent shift, 2042", value: `${R("proj_P_2042_median").toFixed(1)}%`, foot: `${R("proj_P_2042_p05").toFixed(1)}% to ${R("proj_P_2042_p95").toFixed(1)}%` },
      { label: "Temporary shock, 2042", value: `${R("proj_T_2042_median").toFixed(1)}%`, foot: `${R("proj_T_2042_p05").toFixed(1)}% to ${R("proj_T_2042_p95").toFixed(1)}%` },
      { label: "Gap between medians", value: `${R("proj_gap_2042").toFixed(1)} pp`, foot: "2042", accent: true },
      { label: "Draws rejected", value: `${R("cal_n_draws") - R("cal_n_ok")}`, foot: `${R("cal_n_rejected_2018")} at 2018, ${R("cal_n_draws") - R("cal_n_ok") - R("cal_n_rejected_2018")} at 2022` },
    ]),
    h("p", { class: "prose", html: `A draw is rejected when 7% cannot be reached in 2018 even at very high exposure (${R("cal_n_rejected_2018")} draws), or when 11% cannot be reached by 2022 in one of the scenarios. The rejected draws have much faster recovery than the kept ones (median ${R("cal_rej_median_gamma2").toFixed(2)} against ${R("cal_kept_median_gamma2").toFixed(2)} per year): with recovery that fast, prevalence cannot rise from 7% to 11% in four years. The speed of the rise therefore rules out the fastest recovery rates.` }),
    source("thesis Chapter 7; python/run_all.py (calibration and projections)."),
  ] }));

  // ---- one draw at a time (live)
  const drawBox = h("div");
  let idx = 0;
  const drawChart = chartCard({ title: "One parameter draw", sub: "computed live in your browser", filename: "one-draw",
    render: () => {
      const d = DRAWS[idx];
      return timeSeries({ xs: yrs, yTitle: "Problematic use", yMax: 0.24, series: [
        { key: "P", name: "Persistent shift", color: COLORS.P, ys: scenarioPath(d.rates, d.cal, "P") },
        { key: "T", name: "Temporary shock", color: COLORS.T, dash: "6 4", ys: scenarioPath(d.rates, d.cal, "T") },
      ], points: [{ x: 2018, y: d.t18, color: "var(--ink)" }, { x: 2022, y: d.t22, color: "var(--ink)" }] });
    } });
  const info = h("div");
  const showInfo = () => {
    const d = DRAWS[idx];
    info.replaceChildren(table(["Quantity", "Value"], [
      ...["eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"].map((k) => [`${PARAM_META[k].sym} ${PARAM_META[k].label}`, d[k].toFixed(4)]),
      ["b₀ exposure before 2018 (solved)", d.b0.toFixed(3)], ["b₁ exposure from 2018, persistent shift (solved)", d.b1.toFixed(3)],
      ["Δα shock in 2020 and 2021, temporary shock (solved)", d.d_alpha.toFixed(3)],
      ["R₀ before / after 2018 (α = 0)", `${R0({ ...d.rates, b: d.b0, alpha: 0 }).toFixed(2)} / ${R0({ ...d.rates, b: d.b1, alpha: 0 }).toFixed(2)}`],
      ["Targets drawn within the rounding of 7% and 11%", `${pc(d.t18, 2)} and ${pc(d.t22, 2)}`],
    ], ["l", "r"]));
  };
  const slider = sliderRow({ sym: "#", label: "Draw", unit: "", desc: "Pick any of the kept draws", min: 0, max: DRAWS.length - 1, step: 1, value: 0, baseline: 0,
    onInput: (v) => { idx = v; drawChart.update(); showInfo(); } });
  showInfo();
  drawBox.append(slider, drawChart, info);
  el.appendChild(panel({ title: "Look inside one draw", sub: "each draw reproduces both figures exactly, yet the two scenarios part after 2022", children: [drawBox] }));

  // ---- next wave
  const vals = [9, 10, 11, 12, 13, 14];
  const fmt = (v) => (v >= 10 ? v.toFixed(0) : v >= 1 ? v.toFixed(1) : v.toFixed(2));
  el.appendChild(panel({ title: "What the 2025/2026 HBSC wave can decide", sub: "likelihood ratio of the persistent shift against the temporary shock, for each possible published figure", children: [
    table(["Published figure", "Implicit prior, kept draws", "Implicit prior, larger set", "Log-uniform prior", "Flat prior"],
      vals.map((v) => [`${v}%`, fmt(R(`calx_lr_${v}`)), fmt(R(`p2_implicit_lr_${v}`)), fmt(R(`p2_log_lr_${v}`)), fmt(R(`p2_flat_lr_${v}`))]), ["l", "r", "r", "r", "r"]),
    h("p", { class: "hint", html: `The first column uses the ${R("cal_n_ok")} kept draws (thesis Chapter 7). The other columns use the larger set of draws of the calibration paper, which also re-weights to explicit priors; the two implicit-prior columns differ by Monte Carlo and sampling variation.` }),
    h("p", { class: "prose", html:
      "A ratio above 1 favours the persistent shift and below 1 the temporary shock. A figure of 14% favours the persistent shift under the implicit and log-uniform priors. " +
      "For 10% or 11% the direction of the inference reverses between the implicit prior and the explicit priors, and the flat-prior ratios depend on the bounds placed on the solved rates and rest on the tails of the predictive distributions. " +
      "The calibration design is itself a Bayesian calibration with an implicit prior, so it should be reported alongside explicit alternatives." }),
    source("thesis Chapter 7 (python/ext_calibration.py) and the calibration paper (Mutibura and Lubuma, manuscript in preparation; python/ext_paper2.py)."),
  ] }));
  return { el };
}
