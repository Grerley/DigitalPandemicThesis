// explorer.js: Screen 3, the model explorer (live), with calibration to the HBSC figures for any set of rates.

import { h, debounce } from "../dom.js";
import { panel, metricGrid, chartCard, legend, sliderRow, condNote } from "../ui.js";
import { timeSeries, bars } from "../charts.js";
import { stateDiagram } from "../diagram.js";
import { viewHead, pc, COLORS, source } from "../common.js";
import { PARAM_META, RANGES, LAB, centralParams } from "../../model/params.js";
import { simulate, R0, R0Factors, equilibrium, floor, isStable, elasticities, PARAM_KEYS } from "../../model/sair.js";
import { calibrateRates, scenarioPath } from "../../model/calibrate.js";
import { store } from "../../state.js";

const BOUNDS = {
  b: [0, 8, 0.01], eta: [0.2, 0.8, 0.01], alpha: [0, 0.06, 0.0005], sigma: [0.1, 0.3, 0.005], gamma1: [0.15, 0.6, 0.005],
  gamma2: [0.15, 1.25, 0.005], delta: [0.05, 0.25, 0.005], nu: [0.0625, 0.0714, 0.0001],
};

export function ExplorerView() {
  const el = h("div", { class: "view" });
  el.appendChild(viewHead("Screen 3 · Model explorer", "Change the rates and see the model respond",
    "Sliders run over the ranges the thesis carries for each rate (the exposure rate <var>b</var> is solved from data in the thesis, so its slider is wide). " +
    "Everything on this screen is computed live by the JavaScript port of the thesis model, starting from the central parameters after 2018."));

  const P = () => store.get().params;
  const x22 = simulate(centralParams(), LAB.central.x18, [4])[0];
  const central = centralParams();

  const metrics = h("div");
  const traj = chartCard({ title: "Trajectory from 2022", sub: "all four states, 2022 to 2042", filename: "explorer-trajectory",
    csv: () => { const t = Array.from({ length: 21 }, (_, i) => i); return simulate(P(), x22, t).map((x, i) => ({ year: 2022 + i, S: x[0], A: x[1], I: x[2], R: x[3] })); },
    render: () => {
      const t = Array.from({ length: 81 }, (_, i) => i / 4), xs = simulate(P(), x22, t);
      const eq = equilibrium(P());
      return timeSeries({ xs: t.map((v) => 2022 + v), yTitle: "Share of young people aged 10 to 24",
        series: ["S", "A", "I", "R"].map((k, j) => ({ key: k, name: { S: "Susceptible", A: "At-Risk", I: "Problematic use (I)", R: "Recovered" }[k], color: COLORS[k], ys: xs.map((x) => x[j]), width: k === "I" ? 2.6 : 1.8 })),
        markers: [], points: [{ x: 2042, y: eq[2], color: COLORS.I, label: `I* ${pc(eq[2])}` }] });
    } });
  const diagWrap = h("div");
  const elastCard = chartCard({ title: "Sensitivity of equilibrium prevalence", sub: "normalised index ∂ln I* / ∂ln θ", filename: "sensitivity",
    render: () => {
      const e = elasticities(P(), "I");
      return bars({ rows: PARAM_KEYS.map((k) => ({ label: `${PARAM_META[k].sym}  ${PARAM_META[k].label}`, value: e[k], color: e[k] < 0 ? "var(--c-s)" : "var(--c-i)" })),
        xTitle: "1% rise in the rate changes I* by this many per cent" });
    } });

  const calOut = h("div");
  const calChart = h("div");
  function runCalibration() {
    const rates = Object.fromEntries(["eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"].map((k) => [k, P()[k]]));
    const cal = calibrateRates(rates);
    calChart.replaceChildren();
    if (!cal.ok) {
      calOut.replaceChildren(h("div", { class: "caveat", html: `<strong>No exact fit.</strong> ${cal.reason}` }));
      return;
    }
    const yrs = Array.from({ length: 25 }, (_, i) => 2018 + i);
    calChart.appendChild(chartCard({ title: "Both explanations at these rates", filename: "explorer-calibration",
      render: () => timeSeries({ xs: yrs, yTitle: "Problematic use", series: [
        { key: "P", name: "Persistent shift", color: COLORS.P, ys: scenarioPath(rates, cal, "P") },
        { key: "T", name: "Temporary shock", color: COLORS.T, dash: "6 4", ys: scenarioPath(rates, cal, "T") }],
        points: [{ x: 2018, y: 0.07, label: "7%", color: "var(--ink)" }, { x: 2022, y: 0.11, label: "11%", color: "var(--ink)" }] }) }));
    calOut.replaceChildren(metricGrid([
      { label: "b₀ (solved)", value: cal.b0.toFixed(3), foot: "exposure before 2018" },
      { label: "b₁ (solved)", value: cal.b1.toFixed(3), foot: "from 2018, persistent shift" },
      { label: "Δα (solved)", value: cal.dAlpha.toFixed(3), foot: "2020 and 2021, temporary shock" },
      { label: "R₀ after 2018", value: R0({ ...rates, b: cal.b1, alpha: 0 }).toFixed(2), foot: "persistent shift, α = 0" },
    ]), h("div", { class: "row", style: { marginTop: "10px" } }, [
      h("button", { class: "sm", onclick: () => { store.set({ params: { b: cal.b1 } }); refresh(); syncSliders(); }, text: "Use b₁ in the explorer" })]));
  }

  function refresh() {
    const p = P(), f = R0Factors({ ...p, alpha: 0 }), eq = equilibrium(p);
    metrics.replaceChildren(metricGrid([
      { label: "R₀ (α = 0)", value: R0({ ...p, alpha: 0 }).toFixed(2), foot: `${f.single.toFixed(2)} per spell × ${f.multiplier.toFixed(3)} spells`, accent: true },
      { label: "Equilibrium I*", value: pc(eq[2]), foot: isStable(p, eq) ? "locally stable" : "not stable" },
      { label: "Floor at b = 0", value: pc(floor(p), 2), foot: p.alpha > 0 ? "no exposure policy goes below" : "zero when α = 0" },
      { label: "Relapse multiplier", value: f.multiplier.toFixed(3), foot: "expected spells of risk" },
    ]));
    traj.update(); elastCard.update();
    diagWrap.replaceChildren(stateDiagram({ state: { S: eq[0], A: eq[1], I: eq[2], R: eq[3] }, params: p }));
  }
  const refreshSoon = debounce(refresh, 30);

  const rows = {};
  const sliders = h("div", {}, PARAM_KEYS.map((k) => {
    const [mn, mx, st] = BOUNDS[k], m = PARAM_META[k];
    const range = RANGES[k] ? ` · thesis range ${RANGES[k].low} to ${RANGES[k].high}` : " · solved from data in the thesis";
    rows[k] = sliderRow({ sym: m.sym, label: m.label, unit: `${m.unit}${range}`, desc: m.desc, min: mn, max: mx, step: st, value: P()[k], baseline: central[k],
      onInput: (v) => { store.set({ params: { [k]: v } }, { silent: true }); refreshSoon(); } });
    return rows[k];
  }));
  const syncSliders = () => { for (const k of PARAM_KEYS) rows[k].setValue(P()[k]); };

  el.appendChild(h("div", { class: "split" }, [
    h("div", {}, [panel({ title: "Live model", children: [metrics, traj] }), panel({ title: "Where young people go", sub: "flows per year at the equilibrium", children: [diagWrap] })]),
    panel({ title: "Rates", actions: h("button", { class: "sm ghost", onclick: () => { store.resetParams(); syncSliders(); refresh(); }, text: "Reset to central" }), children: [sliders] }),
  ]));
  el.appendChild(panel({ title: "Sensitivity", children: [elastCard,
    h("p", { class: "prose", html: "At the central parameters equilibrium prevalence is about three times as sensitive to recovery as to exposure. Remission and recovery pull prevalence down; exposure, progression, relapse and the influence of At-Risk contacts push it up; at these values non-social acquisition and turnover have little effect on the level, although non-social acquisition sets the floor." })] }));
  el.appendChild(panel({ title: "Calibrate these rates to the HBSC figures", sub: "solve b₀, b₁ and Δα exactly for the current values of the other rates (the exposure slider is ignored)", children: [
    h("button", { class: "primary sm", onclick: runCalibration, text: "Calibrate" }), calOut, calChart,
    condNote("The two scenarios always agree on 2018 and 2022 by construction. Where they part after 2022 is the identifiability point of the thesis."),
    source("thesis Chapters 4, 5 and 7; JavaScript port checked against python/dpm/model.py and python/dpm/calibrate.py by the self-test."),
  ] }));
  refresh();
  return { el };
}
