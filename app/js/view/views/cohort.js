// cohort.js — Screen 5: Cohort & Risk Factors.
// -----------------------------------------------------------------------------

import { h, clear, toast } from "../dom.js";
import { panel, metricGrid, chartCard, legend, condNote } from "../ui.js";
import { timeSeries, forest, sankey, heatmap } from "../charts.js";
import { generateCohort } from "../../model/cohort.js";
import { COMPARTMENTS } from "../../model/params.js";
import { computeScenario } from "../scenario.js";
import { pct, fx } from "../format.js";
import { store } from "../../state.js";

const COMP_COLORS = ["var(--c-s)", "var(--c-a)", "var(--c-i)", "var(--c-r)"];

export function CohortView() {
  const el = h("div", { class: "view" });
  let seed = 42;
  let data = null;

  function regen() {
    const p = computeScenario(store.get()).params;
    data = generateCohort(p, { n: 6000, years: 10, seed });
  }
  regen();

  const modeBtn = h("button", { class: "sm", dataset: { mode: "stacked" }, onclick: () => {
    modeBtn.dataset.mode = modeBtn.dataset.mode === "line" ? "stacked" : "line";
    modeBtn.textContent = modeBtn.dataset.mode === "line" ? "Stacked view" : "Line view";
    prevCard.update();
  } }, ["Line view"]);

  const prevCard = chartCard({
    title: "Cohort prevalence", sub: "synthetic 6,000-youth cohort over 10 years",
    filename: "cohort-prevalence",
    csv: () => data.byYear.map((r) => ({ year: r.year, S: fx(r.S, 4), A: fx(r.A, 4), I: fx(r.I, 4), R: fx(r.R, 4) })),
    render: () => {
      const xs = data.byYear.map((r) => r.year);
      const series = COMPARTMENTS.map((c, i) => ({ key: c.key, name: c.name, color: COMP_COLORS[i], dash: c.dash, ys: data.byYear.map((r) => r[c.key]) }));
      return timeSeries({ xs, series, xTitle: "Year", yTitle: "Cohort share", mode: modeBtn.dataset.mode });
    },
  });
  prevCard.querySelector(".chart-toolbar").insertBefore(modeBtn, prevCard.querySelector(".chart-toolbar .spacer").nextSibling);

  const forestCard = chartCard({
    title: "Risk-factor odds ratios", sub: "logistic regression recovers the published ORs",
    filename: "odds-ratio-forest",
    csv: () => data.oddsRatios.map((o) => ({ factor: o.label, known_OR: o.knownOR, est_OR: fx(o.estOR, 3), lo: fx(o.lo, 3), hi: fx(o.hi, 3) })),
    render: () => forest({ rows: data.oddsRatios }),
  });

  const sankeyCard = chartCard({
    title: "State flows", sub: "baseline → final state over 10 years",
    filename: "cohort-sankey",
    render: () => sankey({ flows: data.flows, labels: data.states, colors: COMP_COLORS }),
  });

  const heatKnown = chartCard({ title: "", filename: "transition-known", render: () => heatmap({ matrix: data.knownMatrix, labels: data.states, title: "Known (generating)", accent: "--c-s" }) });
  const heatEst = chartCard({ title: "", filename: "transition-estimated", render: () => heatmap({ matrix: data.estMatrix, labels: data.states, title: "Estimated (empirical)", accent: "--c-i" }) });

  const metricsWrap = h("div");
  function renderMetrics() {
    const last = data.byYear[data.byYear.length - 1];
    clear(metricsWrap);
    metricsWrap.appendChild(metricGrid([
      { label: "Cohort size", value: data.n.toLocaleString(), foot: `${data.years}-year follow-up` },
      { label: "Final addicted", value: pct(last.I), accent: true, foot: "synthetic prevalence" },
      { label: "Matrix recovery MAE", value: fx(data.matrixMAE, 3), foot: "known vs estimated", accent: true },
      { label: "Strongest OR", value: fx(Math.max(...data.oddsRatios.map((o) => o.estOR)), 2), foot: "high depression" },
    ]));
  }

  function rebuild() { prevCard.update(); forestCard.update(); sankeyCard.update(); heatKnown.update(); heatEst.update(); renderMetrics(); }

  const regenBtn = h("button", { class: "sm", onclick: () => { seed = (seed * 1103515245 + 12345) % 2147483647; regen(); rebuild(); toast(`Regenerated cohort (seed ${seed})`); } }, ["New random seed"]);

  el.appendChild(h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: "Screen 5 · Statistical layer" }),
    h("h1", { text: "Cohort & Risk Factors" }),
    h("p", { class: "prose", text: "A synthetic longitudinal cohort — a modelling instrument, not observed data — whose individual onset probability is tilted by covariates through published odds ratios. It lets us recover the risk-factor forest plot by logistic regression and validate transition-matrix estimation against a known generating matrix." }),
  ]));
  el.appendChild(metricsWrap);
  el.appendChild(h("div", { class: "split" }, [
    panel({ title: "Cohort dynamics", actions: regenBtn, children: [prevCard, legend(COMPARTMENTS.map((c, i) => ({ name: c.name, color: COMP_COLORS[i], dash: c.dash })))] }),
    panel({ title: "State-flow Sankey", children: [sankeyCard] }),
  ]));
  el.appendChild(panel({ title: "Risk-factor forest plot", sub: "estimated OR (●) with 95% CI vs published OR (◇); reference line at 1", badge: "conditional", children: [
    forestCard,
    legend([{ name: "estimated (95% CI)", color: "var(--accent)" }, { name: "published OR", color: "var(--c-i)" }]),
  ] }));
  el.appendChild(panel({ title: "Transition-matrix recovery", sub: "the empirically estimated matrix reproduces the known generating matrix", children: [
    h("div", { class: "grid cols-2" }, [heatKnown, heatEst]),
    h("p", { class: "chart-cap", text: `Mean absolute error between known and estimated matrices: ${fx(data.matrixMAE, 3)}. Rows are annual transition probabilities from each state.` }),
  ] }));
  el.appendChild(condNote("The cohort is fully regenerable from a fixed seed and stands in for the longitudinal transition data that do not yet exist; odds ratios are associational, not causal."));

  renderMetrics();
  const unsub = store.subscribe(() => { regen(); rebuild(); });
  return { el, destroy: unsub };
}
