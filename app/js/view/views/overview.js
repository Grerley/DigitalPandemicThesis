// overview.js — Screen 1: Overview dashboard (baseline S–A–I–R run).
// -----------------------------------------------------------------------------

import { h, clear } from "../dom.js";
import { panel, metricGrid, chartCard, legend, timeline, condNote } from "../ui.js";
import { timeSeries } from "../charts.js";
import { stateDiagram } from "../diagram.js";
import { computeScenario, stateAt } from "../scenario.js";
import { COMPARTMENTS } from "../../model/params.js";
import { pct, fx, yr } from "../format.js";
import { store } from "../../state.js";

export function OverviewView() {
  const el = h("div", { class: "view" });
  let scrubIdx = null; // null = show full series; number = scrub month
  let tl;

  function seriesFor(traj) {
    return COMPARTMENTS.map((c) => ({
      key: c.key, name: c.name, color: c.color, dash: c.dash, ys: traj[c.key],
    }));
  }

  const modeBtn = h("button", { class: "sm", dataset: { mode: "line" }, onclick: () => {
    modeBtn.dataset.mode = modeBtn.dataset.mode === "line" ? "stacked" : "line";
    modeBtn.textContent = modeBtn.dataset.mode === "line" ? "Stacked view" : "Line view";
    chart.update();
  } }, ["Stacked view"]);

  const chart = chartCard({
    title: "Population dynamics", sub: "0–20 year horizon",
    filename: "sair-baseline",
    csv: () => {
      const sc = computeScenario(store.get());
      return sc.traj.month.map((m, i) => ({
        month: m, year: fx(sc.traj.year[i], 3),
        S: fx(sc.traj.S[i], 5), A: fx(sc.traj.A[i], 5), I: fx(sc.traj.I[i], 5), R: fx(sc.traj.R[i], 5),
        incidence: fx(sc.traj.incidence[i], 6),
      }));
    },
    render: () => {
      const st = store.get();
      const sc = computeScenario(st);
      const s = sc.summary;
      const series = seriesFor(sc.traj);
      // Compare overlay: pinned scenario's Addicted curve.
      if (st.compare && modeBtn.dataset.mode === "line") {
        const cmp = computeScenario({ ...st, params: st.compare.params, interventions: st.compare.interventions || [], customMods: {} });
        series.push({ key: "cmpI", name: `${st.compare.label} · I`, color: "var(--ink-2)", dash: "5 3", ys: cmp.traj.I, width: 1.6, opacity: 0.85 });
      }
      return timeSeries({
        xs: sc.traj.year, series,
        mode: modeBtn.dataset.mode, yTitle: "Population fraction", xTitle: "Year",
        markers: [{ x: s.peakYear, label: "peak I", color: "var(--c-i)" }],
        points: modeBtn.dataset.mode === "line" ? [{ x: s.peakYear, y: s.peakPrevalence, color: "var(--c-i)", label: pct(s.peakPrevalence) }] : [],
        scrubX: scrubIdx == null ? null : sc.traj.year[scrubIdx],
      });
    },
  });

  chart.querySelector(".chart-toolbar").insertBefore(modeBtn, chart.querySelector(".chart-toolbar .spacer").nextSibling);

  const diagWrap = h("div");
  const metricsWrap = h("div");

  function renderDiagram() {
    const st = store.get();
    const sc = computeScenario(st);
    const state = scrubIdx == null ? sc.summary.endemicState : stateAt(sc.traj, scrubIdx);
    // For the diagram we show flows at the scrubbed (or final) state.
    clear(diagWrap);
    diagWrap.appendChild(stateDiagram({ state, params: sc.params }));
  }

  function renderMetrics() {
    const sc = computeScenario(store.get());
    const s = sc.summary, r = sc.r0;
    clear(metricsWrap);
    metricsWrap.appendChild(metricGrid([
      { label: "R₀ operative", value: fx(r.operative, 2), foot: "NGM at S₀ = 0.823", accent: true },
      { label: "R₀ basic", value: fx(r.basic, 2), foot: "NGM at S₀ = 1" },
      { label: "R₀ crude", value: fx(r.crude, 2), foot: "β·k·D, D = 1/γ₂" },
      { label: "Peak addicted", value: pct(s.peakPrevalence), foot: `at ${yr(s.peakYear)}`, accent: true },
      { label: "Endemic addicted", value: pct(s.endemicPrevalence), foot: "t = 20 yr" },
      { label: "Peak incidence", value: pct(s.peakIncidence, 2), foot: `at ${yr(s.peakIncidenceYear)}` },
    ]));
  }

  function rebuild() {
    chart.update();
    renderDiagram();
    renderMetrics();
  }

  // timeline over months
  function buildTimeline() {
    const sc = computeScenario(store.get());
    const max = sc.traj.month.length - 1;
    tl = timeline({
      max, value: scrubIdx ?? max,
      labelFor: (v) => `${(v / 12).toFixed(1)} yr · month ${v}`,
      onScrub: (v) => { scrubIdx = v; chart.update(); renderDiagram(); },
    });
    return tl.el;
  }

  el.appendChild(h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: "Screen 1 · Deterministic core" }),
    h("h1", { text: "Overview dashboard" }),
    h("p", { class: "prose", text: "The calibrated S–A–I–R baseline: a nonlinear difference-equation system on a monthly step. Watch the epidemic rise, peak, and settle to an endemic plateau sustained by relapse. Scrub the timeline to inspect the state and live transition rates at any month." }),
  ]));

  el.appendChild(metricsWrap);
  const tlHolder = h("div");
  el.appendChild(panel({
    title: "Four-compartment trajectory",
    sub: "Susceptible → At-Risk → Addicted → Recovered, with relapse R→A",
    children: [chart, legend(COMPARTMENTS.map((c) => ({ name: `${c.key} — ${c.name}`, color: c.color, dash: c.dash }))), tlHolder],
  }));

  el.appendChild(h("div", { class: "split" }, [
    panel({ title: "State-transition diagram", sub: "Edge thickness ∝ live monthly transition probability", badge: "conditional",
      children: [diagWrap, h("p", { class: "chart-cap", text: "Flows animate in the direction of transition; the S→A edge is the mass-action force of infection and grows as the epidemic spreads." })] }),
    panel({ title: "Reading this model", children: [
      h("p", { class: "prose", html: "Because the S→A probability depends on the current prevalence, the map is <em>nonlinear</em> — a fixed Markov chain could not produce a rise–peak–decline curve. The relapse term δ keeps Recovered non-absorbing, which is what sustains the endemic plateau near 8%." }),
      condNote("All values are outputs of the calibrated model, not empirical forecasts. See Methods & Provenance for the equations and parameter sources."),
    ] }),
  ]));

  rebuild();
  tlHolder.appendChild(buildTimeline());

  const unsub = store.subscribe(() => { scrubIdx = null; clear(tlHolder); rebuild(); tlHolder.appendChild(buildTimeline()); });
  return { el, destroy: () => { unsub(); tl?.stop(); } };
}
