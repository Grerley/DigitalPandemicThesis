// explorer.js — Screen 2: Parameter Explorer.
// -----------------------------------------------------------------------------

import { h, clear, debounce } from "../dom.js";
import { panel, metricGrid, chartCard, legend, sliderRow, condNote } from "../ui.js";
import { timeSeries, tornado } from "../charts.js";
import { computeScenario } from "../scenario.js";
import { simulate, summarise } from "../../model/sair.js";
import { r0Trio } from "../../model/r0.js";
import { COMPARTMENTS, PARAM_META, PARAM_KEYS, DEFAULT_PARAMS, baselineInit } from "../../model/params.js";
import { pct, fx, yr } from "../format.js";
import { store } from "../../state.js";

export function ExplorerView() {
  const el = h("div", { class: "view" });
  const metricsWrap = h("div");
  const sliderRows = {};
  // The Explorer manipulates the base parameters directly, so it deliberately
  // ignores any active intervention preset (those live on Screen 6).
  const baseState = () => ({ ...store.get(), interventions: [], customMods: {} });

  const chart = chartCard({
    title: "Live trajectory", filename: "explorer-trajectory",
    csv: () => {
      const sc = computeScenario(baseState());
      return sc.traj.month.map((m, i) => ({ month: m, year: fx(sc.traj.year[i], 3),
        S: fx(sc.traj.S[i], 5), A: fx(sc.traj.A[i], 5), I: fx(sc.traj.I[i], 5), R: fx(sc.traj.R[i], 5) }));
    },
    render: () => {
      const sc = computeScenario(baseState());
      const base = simulate(DEFAULT_PARAMS, baselineInit(), store.get().horizonYears || 20);
      const series = COMPARTMENTS.map((c) => ({ key: c.key, name: c.name, color: c.color, dash: c.dash, ys: sc.traj[c.key] }));
      // faint baseline I for reference
      series.push({ key: "baseI", name: "baseline I", color: "var(--ink-3)", dash: "1 3", ys: base.I, width: 1.4, opacity: 0.7, label: false });
      return timeSeries({ xs: sc.traj.year, series, yTitle: "Population fraction",
        markers: [{ x: sc.summary.peakYear, label: "peak I", color: "var(--c-i)" }] });
    },
  });

  const tornadoCard = chartCard({
    title: "One-way sensitivity", sub: "peak prevalence across each parameter's range",
    filename: "tornado-sensitivity",
    render: () => tornado(buildTornado()),
  });

  function buildTornado() {
    const st = store.get();
    const p0 = { ...st.params };
    const init = baselineInit();
    const basePeak = summarise(simulate(p0, init, 20)).peakPrevalence;
    const rows = PARAM_KEYS.map((key) => {
      const meta = PARAM_META[key];
      const lowP = { ...p0, [key]: meta.min };
      const highP = { ...p0, [key]: meta.max };
      return {
        label: `${meta.sym}  ${meta.label}`,
        low: summarise(simulate(lowP, init, 20)).peakPrevalence,
        high: summarise(simulate(highP, init, 20)).peakPrevalence,
        span: 0,
      };
    });
    // sort by influence (span)
    for (const r of rows) r.span = Math.abs(r.high - r.low);
    rows.sort((a, b) => b.span - a.span);
    return { rows, baseline: basePeak, valueFmt: (v) => pct(v, 0), metricLabel: "Peak addicted prevalence" };
  }

  function renderMetrics() {
    const sc = computeScenario(baseState());
    const s = sc.summary, r = sc.r0;
    clear(metricsWrap);
    metricsWrap.appendChild(metricGrid([
      { label: "R₀ operative", value: fx(r.operative, 2), accent: true, foot: "vs 2.50 baseline" },
      { label: "R₀ basic", value: fx(r.basic, 2), foot: "S₀ = 1" },
      { label: "Peak addicted", value: pct(s.peakPrevalence), foot: `at ${yr(s.peakYear)}`, accent: true },
      { label: "Endemic addicted", value: pct(s.endemicPrevalence), foot: "t = 20 yr" },
    ]));
  }

  const update = debounce(() => { chart.update(); renderMetrics(); tornadoCard.update(); }, 40);

  function buildSliders() {
    const st = store.get();
    return PARAM_KEYS.map((key) => {
      const meta = PARAM_META[key];
      const row = sliderRow({
        sym: meta.sym, label: meta.label, unit: meta.unit, desc: meta.desc,
        min: meta.min, max: meta.max, step: meta.step,
        value: st.params[key], baseline: DEFAULT_PARAMS[key],
        onInput: (v) => { store.setParam(key, v, { silent: true }); update(); },
      });
      sliderRows[key] = row;
      return row;
    });
  }

  const resetBtn = h("button", { class: "primary sm", onclick: () => {
    store.resetParams({ silent: true });
    for (const key of PARAM_KEYS) sliderRows[key].setValue(DEFAULT_PARAMS[key]);
    update();
  } }, ["Reset to calibrated baseline"]);

  el.appendChild(h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: "Screen 2 · Parameter Explorer" }),
    h("h1", { text: "Parameter Explorer" }),
    h("p", { class: "prose", text: "Move any of the eight monthly rates and watch the trajectory, R₀, peak and endemic level recompute in real time. The faint dotted line marks the calibrated baseline for I(t). The tornado ranks each parameter's leverage over the peak." }),
  ]));

  el.appendChild(metricsWrap);
  el.appendChild(h("div", { class: "split" }, [
    panel({ title: "Trajectory", children: [chart, legend(COMPARTMENTS.map((c) => ({ name: c.name, color: c.color, dash: c.dash })))] }),
    panel({ title: "Parameters", actions: resetBtn, children: buildSliders() }),
  ]));
  el.appendChild(panel({ title: "Sensitivity (tornado)", sub: "each bar spans the peak prevalence as one parameter sweeps its full range, others held at current values",
    badge: "conditional", children: [tornadoCard] }));
  el.appendChild(condNote("Sweeping parameters explores model behaviour under counterfactual assumptions; only the calibrated baseline is fitted to data."));

  renderMetrics();

  const unsub = store.subscribe(() => {
    for (const key of PARAM_KEYS) sliderRows[key]?.setValue(store.get().params[key]);
    chart.update(); renderMetrics(); tornadoCard.update();
  });
  return { el, destroy: unsub };
}
