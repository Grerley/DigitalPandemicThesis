// interventions.js — Screen 6: Intervention Simulator.
// -----------------------------------------------------------------------------

import { h, clear } from "../dom.js";
import { panel, metricGrid, chartCard, legend, condNote } from "../ui.js";
import { timeSeries } from "../charts.js";
import { simulate, summarise } from "../../model/sair.js";
import { INTERVENTIONS, COMPREHENSIVE, COST_ASSUMPTIONS, applyMods } from "../../model/interventions.js";
import { activeMods, effectiveParams } from "../scenario.js";
import { baselineInit, DEFAULT_PARAMS } from "../../model/params.js";
import { pct, fx, signedPct, zar } from "../format.js";
import { store } from "../../state.js";

export function InterventionsView() {
  const el = h("div", { class: "view" });
  const init = baselineInit();
  const horizon = () => store.get().horizonYears || 20;

  function baseTraj() { return simulate(store.get().params, init, horizon()); }
  function scenTraj() { return simulate(effectiveParams(store.get()), init, horizon()); }

  const chart = chartCard({
    title: "Intervention vs baseline", sub: "Addicted prevalence I(t)",
    filename: "intervention-overlay",
    csv: () => { const b = baseTraj(), s = scenTraj(); return b.year.map((y, i) => ({ year: fx(y, 3), baseline_I: fx(b.I[i], 5), scenario_I: fx(s.I[i], 5) })); },
    render: () => {
      const b = baseTraj(), sc = scenTraj();
      return timeSeries({ xs: b.year, yTitle: "Addicted fraction", xTitle: "Year", series: [
        { key: "base", name: "baseline", color: "var(--ink-3)", dash: "3 3", ys: b.I },
        { key: "scen", name: "your scenario", color: "var(--c-i)", ys: sc.I },
      ], markers: [
        { x: summarise(b).peakYear, color: "var(--ink-3)" },
        { x: summarise(sc).peakYear, label: "peak", color: "var(--c-i)" },
      ] });
    },
  });

  const metricsWrap = h("div");
  const tableWrap = h("div");
  const costWrap = h("div");
  const modsWrap = h("div");

  function renderAll() {
    chart.update();
    renderMetrics();
    renderTable();
    renderCost();
    renderMods();
  }

  function renderMetrics() {
    const b = summarise(baseTraj()), s = summarise(scenTraj());
    const red = (b.peakPrevalence - s.peakPrevalence) / b.peakPrevalence;
    const delay = s.peakYear - b.peakYear;
    clear(metricsWrap);
    metricsWrap.appendChild(metricGrid([
      { label: "Baseline peak", value: pct(b.peakPrevalence), foot: `at ${fx(b.peakYear, 1)} yr` },
      { label: "Scenario peak", value: pct(s.peakPrevalence), accent: true, foot: `at ${fx(s.peakYear, 1)} yr` },
      { label: "Peak reduction", value: signedPct(red), accent: true, foot: "vs baseline" },
      { label: "Peak delay", value: `${delay >= 0 ? "+" : ""}${fx(delay, 1)} yr`, foot: "later peak" },
      { label: "Endemic", value: pct(s.endemicPrevalence), foot: `from ${pct(b.endemicPrevalence)}` },
    ]));
  }

  function renderTable() {
    const b = summarise(baseTraj());
    const rows = [];
    for (const [key, spec] of Object.entries(INTERVENTIONS)) {
      const s = summarise(simulate(applyMods(store.get().params, spec.mods), init, horizon()));
      rows.push({ key, label: spec.label, s, mods: spec.mods });
    }
    const comp = summarise(simulate(applyMods(store.get().params, COMPREHENSIVE.mods), init, horizon()));
    rows.push({ key: "comprehensive", label: COMPREHENSIVE.label, s: comp, mods: COMPREHENSIVE.mods });

    const table = h("table", { class: "data" }, [
      h("thead", {}, h("tr", {}, [
        h("th", { text: "Scenario" }), h("th", { text: "Peak I" }), h("th", { text: "Peak yr" }),
        h("th", { text: "Endemic I" }), h("th", { text: "Peak reduction" }), h("th", {}),
      ])),
      h("tbody", {}, [
        h("tr", {}, [
          h("td", {}, [h("span", { class: "scen-dot", style: { background: "var(--ink-3)" } }), "Baseline"]),
          h("td", { text: pct(b.peakPrevalence) }), h("td", { text: fx(b.peakYear, 1) }),
          h("td", { text: pct(b.endemicPrevalence) }), h("td", { text: "—" }), h("td", {}),
        ]),
        ...rows.map((r) => {
          const red = (b.peakPrevalence - r.s.peakPrevalence) / b.peakPrevalence;
          const applyBtn = h("button", { class: "sm ghost", onclick: () => {
            if (r.key === "comprehensive") store.set({ interventions: ["comprehensive"], customMods: {} });
            else store.set({ interventions: [r.key], customMods: {} });
          } }, ["apply"]);
          return h("tr", {}, [
            h("td", {}, [h("span", { class: "scen-dot", style: { background: r.key === "comprehensive" ? "var(--c-i)" : "var(--accent)" } }), r.label]),
            h("td", { text: pct(r.s.peakPrevalence) }), h("td", { text: fx(r.s.peakYear, 1) }),
            h("td", { text: pct(r.s.endemicPrevalence) }),
            h("td", { class: "pos", text: signedPct(red) }), h("td", {}, applyBtn),
          ]);
        }),
      ]),
    ]);
    clear(tableWrap); tableWrap.appendChild(table);
  }

  function renderCost() {
    const b = baseTraj(), s = scenTraj();
    // Addicted person-years (illustrative): monthly I summed × (1/12) × population.
    const pop = COST_ASSUMPTIONS.populationYouth;
    const py = (traj) => traj.I.reduce((a, v) => a + v, 0) / 12 * pop;
    const averted = py(b) - py(s);
    const mods = activeMods(store.get());
    // Unit cost: sum of active named interventions' per-youth-year cost.
    let unit = 0;
    const ivs = store.get().interventions;
    if (ivs.includes("comprehensive")) unit = COMPREHENSIVE.cost;
    else for (const k of ivs) unit += INTERVENTIONS[k]?.cost || 0;
    const totalCost = unit * pop * horizon();
    const perAverted = averted > 1 ? totalCost / averted : null;
    clear(costWrap);
    if (!Object.keys(mods).length) {
      costWrap.appendChild(h("p", { class: "hint", text: "Select one or more interventions to see the illustrative cost framing." }));
      return;
    }
    costWrap.appendChild(metricGrid([
      { label: "Addicted-person-years averted", value: Math.round(averted).toLocaleString(), foot: `over ${horizon()} yr, ${(pop / 1e6).toFixed(1)}M youth` },
      { label: "Programme cost (illustrative)", value: zar(totalCost), foot: `${zar(unit)}/youth/yr` },
      { label: "Cost per person-year averted", value: perAverted ? zar(perAverted) : "n/a", accent: true, foot: "lower is better" },
    ]));
  }

  function renderMods() {
    const mods = activeMods(store.get());
    clear(modsWrap);
    if (!Object.keys(mods).length) { modsWrap.appendChild(h("span", { class: "hint", text: "No levers active — baseline parameters." })); return; }
    modsWrap.appendChild(h("div", { class: "row" }, Object.entries(mods).map(([k, v]) =>
      h("span", { class: "pill", text: `${k} × ${fx(v, 2)}` }))));
  }

  // Toggle list
  const toggleList = h("div", { class: "toggle-list" });
  function renderToggles() {
    const ivs = store.get().interventions;
    const comp = ivs.includes("comprehensive");
    clear(toggleList);
    for (const [key, spec] of Object.entries(INTERVENTIONS)) {
      const on = !comp && ivs.includes(key);
      const t = h("div", { class: "toggle", dataset: { on: String(on) }, role: "checkbox", tabindex: "0",
        "aria-checked": String(on), onclick: () => toggle(key), onkeydown: (e) => { if (e.key === " " || e.key === "Enter") { e.preventDefault(); toggle(key); } } }, [
        h("div", { class: "box" }),
        h("div", {}, [h("div", { class: "tl-title", text: spec.label }), h("div", { class: "tl-blurb", text: spec.blurb })]),
      ]);
      toggleList.appendChild(t);
    }
  }
  function toggle(key) {
    let ivs = store.get().interventions.filter((k) => k !== "comprehensive");
    if (ivs.includes(key)) ivs = ivs.filter((k) => k !== key); else ivs.push(key);
    store.set({ interventions: ivs, customMods: {} });
  }

  const compBtn = h("button", { class: "primary sm", onclick: () => store.set({ interventions: ["comprehensive"], customMods: {} }) }, ["Comprehensive bundle"]);
  const clearBtn = h("button", { class: "sm", onclick: () => store.set({ interventions: [], customMods: {} }) }, ["Clear all"]);

  el.appendChild(h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: "Screen 6 · Policy levers" }),
    h("h1", { text: "Intervention Simulator" }),
    h("p", { class: "prose", text: "Interventions modify model parameters multiplicatively — efficacies are assumed, not estimated, so results are comparative. Compose your own bundle and watch the peak fall and shift against the baseline curve. The comprehensive package lands near a 63% peak reduction." }),
  ]));
  el.appendChild(metricsWrap);
  el.appendChild(h("div", { class: "split" }, [
    panel({ title: "Trajectory overlay", children: [chart, legend([{ name: "baseline", color: "var(--ink-3)", dash: "3 3" }, { name: "your scenario", color: "var(--c-i)" }])] }),
    panel({ title: "Compose levers", actions: h("div", { class: "row" }, [compBtn, clearBtn]), children: [toggleList, h("hr", { class: "sep" }), h("div", { style: { marginBottom: "6px", fontSize: "12px", color: "var(--ink-3)" } }, ["Active modifiers"]), modsWrap] }),
  ]));
  el.appendChild(panel({ title: "Peak-reduction table", sub: "each named scenario applied to current base parameters", badge: "conditional", children: [tableWrap] }));
  el.appendChild(panel({ title: "Cost-effectiveness framing (ZAR)", sub: COST_ASSUMPTIONS.note, badge: "illustrative", children: [costWrap] }));
  el.appendChild(condNote("Intervention efficacies are assumed values informed by the direction of published evaluations, not estimated here. Cost figures are illustrative unit costs for comparison, not a programme budget."));

  renderToggles();
  renderAll();
  const unsub = store.subscribe(() => { renderToggles(); renderAll(); });
  return { el, destroy: unsub };
}
