// network-lab.js — Screen 4: Network Lab.
// -----------------------------------------------------------------------------
// Rendered contagion animation on an explicit graph + Monte-Carlo ensembles
// (mean + 95% envelope) comparing topologies, plus hub-vs-random immunisation.
// Ensembles run in a Web Worker (fallback to main thread) so the UI stays live.

import { h, clear, toast } from "../dom.js";
import { panel, metricGrid, chartCard, legend } from "../ui.js";
import { timeline } from "../ui.js";
import { timeSeries } from "../charts.js";
import { makeNetwork, simulateNetworkOnce, STATE } from "../../model/network.js";
import { RNG } from "../../model/rng.js";
import { layout, drawGraph } from "../netgraph.js";
import { runEnsemble, workerAvailable } from "../ensemble-runner.js";
import { computeScenario } from "../scenario.js";
import { pct, fx, yr } from "../format.js";
import { store } from "../../state.js";

const TOPO = { ER: "Erdős–Rényi (random)", WS: "Watts–Strogatz (small-world)", BA: "Barabási–Albert (scale-free)" };
const TOPO_COLOR = { ER: "var(--c-s)", WS: "var(--c-a)", BA: "var(--c-i)" };

export function NetworkLabView() {
  const el = h("div", { class: "view" });
  const st = () => store.get();
  const cfg = { topology: st().network.type, n: st().network.n, meanDegree: st().network.meanDegree,
    runs: st().network.runs, immunizeFrac: st().network.immunizeFrac, immunizeStrategy: st().network.immunizeStrategy,
    horizonMonths: 180, seed: st().network.seed || 1, seedFrac: 0.004 };

  // Calibration-consistent Reed–Frost parameters: the network uses the full
  // node degree, so the per-contact β is scaled to keep the mean transmission
  // pressure β·⟨k⟩ matched to the mean-field calibration (β·k ≈ 0.122/mo). μ is
  // set to 0 so pure social contagion drives the topology comparison. This puts
  // the model near the epidemic threshold, where network structure matters.
  function netParams() {
    const p = computeScenario(st()).params;
    return { ...p, beta: (p.beta * p.k) / cfg.meanDegree, mu: 0 };
  }

  // ---------- live animation ----------
  const VIS_N = () => Math.min(cfg.n, 240);
  let anim = null; // { g, pos, frames }
  const canvas = h("canvas", { role: "img", "aria-label": "Contagion spreading across the network" });
  const ctx = canvas.getContext("2d");
  const animStatus = h("div", { class: "hint" });
  let tlHolder = h("div");

  function sizeCanvas() {
    const cssW = canvas.clientWidth || 640, cssH = Math.round(cssW * 0.6);
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    canvas.width = cssW * dpr; canvas.height = cssH * dpr;
    canvas.style.height = cssH + "px";
    return { w: cssW, h: cssH, dpr };
  }

  function buildAnimation() {
    const dims = sizeCanvas();
    const p = netParams();
    const rng = new RNG(1234);
    const g = makeNetwork(cfg.topology, VIS_N(), cfg.meanDegree, rng);
    const runRng = new RNG(99);
    const run = simulateNetworkOnce(g, p, {
      horizonMonths: cfg.horizonMonths, seedFrac: 0.02, seedHub: false,
      immunizeFrac: cfg.immunizeFrac, immunizeStrategy: cfg.immunizeStrategy,
      recordStates: true, frameStride: 1,
    }, runRng);
    animStatus.textContent = "Laying out graph…";
    // Defer layout so the status paints.
    setTimeout(() => {
      const pos = layout(g.adj, dims.w, dims.h, 7, VIS_N() > 180 ? 130 : 170);
      anim = { g, pos, frames: run.frames, dims };
      drawFrame(0);
      animStatus.textContent = `${VIS_N()} nodes · ⟨k⟩≈${fx(mean(g.degree), 1)} · max degree ${Math.max(...g.degree)}`;
      clear(tlHolder); tlHolder.appendChild(buildTimeline());
    }, 20);
  }

  function drawFrame(fi) {
    if (!anim) return;
    const f = anim.frames[Math.max(0, Math.min(fi, anim.frames.length - 1))];
    drawGraph(ctx, anim.g, anim.pos, f, anim.dims);
  }
  function buildTimeline() {
    const max = (anim?.frames.length || 1) - 1;
    const tl = timeline({ max, value: 0, labelFor: (v) => `${(v / 12).toFixed(1)} yr · month ${v}`,
      onScrub: drawFrame, fps: 18 });
    return tl.el;
  }

  // ---------- ensembles ----------
  const topoStatus = h("div", { class: "hint" });
  const topoMetrics = h("div");
  let topoResults = {};
  const topoCard = chartCard({
    title: "Ensemble trajectories", sub: "mean + 95% envelope, Addicted fraction",
    filename: "network-ensembles",
    csv: () => {
      const keys = Object.keys(topoResults);
      if (!keys.length) return [];
      const yrs = topoResults[keys[0]].year;
      return yrs.map((y, i) => { const row = { year: fx(y, 3) }; for (const k of keys) row[k + "_mean"] = fx(topoResults[k].mean[i], 5); return row; });
    },
    render: () => {
      const keys = Object.keys(topoResults);
      if (!keys.length) return placeholderSVG("Run ensembles to compare topologies");
      const xs = topoResults[keys[0]].year;
      const series = keys.map((k) => ({ key: k, name: TOPO[k].split(" ")[0], color: TOPO_COLOR[k],
        ys: topoResults[k].mean, band: { lo: topoResults[k].lo, hi: topoResults[k].hi } }));
      return timeSeries({ xs, series, yTitle: "Addicted fraction", xTitle: "Year" });
    },
  });

  const immStatus = h("div", { class: "hint" });
  const immMetrics = h("div");
  let immResults = null;
  const immCard = chartCard({
    title: "Immunisation strategy", sub: "protect the same coverage as hubs vs at random (scale-free graph)",
    filename: "immunisation-compare",
    render: () => {
      if (!immResults) return placeholderSVG("Run the immunisation comparison");
      const xs = immResults.hub.year;
      return timeSeries({ xs, yTitle: "Addicted fraction", xTitle: "Year", series: [
        { key: "none", name: "no immunisation", color: "var(--ink-3)", dash: "2 3", ys: immResults.none.mean },
        { key: "random", name: "random coverage", color: "var(--c-a)", ys: immResults.random.mean, band: { lo: immResults.random.lo, hi: immResults.random.hi } },
        { key: "hub", name: "hub-targeted", color: "var(--c-i)", ys: immResults.hub.mean, band: { lo: immResults.hub.lo, hi: immResults.hub.hi } },
      ] });
    },
  });

  let running = false;
  async function runTopologies() {
    if (running) { toast("A run is already in progress"); return; }
    running = true;
    topoResults = {};
    const params = netParams();
    const topos = ["ER", "WS", "BA"];
    for (const t of topos) {
      topoStatus.textContent = `Running ${TOPO[t]} — ${cfg.runs} runs of ${cfg.n} nodes…`;
      const res = await runEnsemble({ topology: t, n: cfg.n, meanDegree: cfg.meanDegree, runs: cfg.runs,
        horizonMonths: cfg.horizonMonths, seed: 42, params, seedHub: false, immunizeFrac: 0, seedFrac: cfg.seedFrac },
        (frac) => { topoStatus.textContent = `Running ${TOPO[t]} — ${Math.round(frac * 100)}%`; });
      topoResults[t] = res;
      topoCard.update();
    }
    renderTopoMetrics();
    topoStatus.textContent = `Done · ${cfg.runs} runs each · ${workerAvailable() ? "web worker" : "main thread"}`;
    running = false;
  }

  function renderTopoMetrics() {
    const items = ["ER", "WS", "BA"].filter((k) => topoResults[k]).map((k) => ({
      label: TOPO[k].split(" (")[0], value: pct(topoResults[k].meanPeak), foot: `peak ~${yr(topoResults[k].meanPeakYear)}`,
      accent: k === "BA",
    }));
    clear(topoMetrics);
    if (items.length) topoMetrics.appendChild(metricGrid(items));
  }

  async function runImmunisation() {
    if (running) { toast("A run is already in progress"); return; }
    running = true;
    const params = netParams();
    const cov = cfg.immunizeFrac > 0 ? cfg.immunizeFrac : 0.3;
    const base = { topology: "BA", n: cfg.n, meanDegree: cfg.meanDegree, runs: cfg.runs, horizonMonths: cfg.horizonMonths, seed: 11, params, seedHub: false, seedFrac: cfg.seedFrac };
    immStatus.textContent = "No immunisation…";
    const none = await runEnsemble({ ...base, immunizeFrac: 0 }, (f) => immStatus.textContent = `No immunisation — ${Math.round(f * 100)}%`);
    immStatus.textContent = "Random coverage…";
    const random = await runEnsemble({ ...base, immunizeFrac: cov, immunizeStrategy: "random" }, (f) => immStatus.textContent = `Random ${Math.round(cov * 100)}% — ${Math.round(f * 100)}%`);
    immStatus.textContent = "Hub-targeted…";
    const hub = await runEnsemble({ ...base, immunizeFrac: cov, immunizeStrategy: "hub" }, (f) => immStatus.textContent = `Hub ${Math.round(cov * 100)}% — ${Math.round(f * 100)}%`);
    immResults = { none, random, hub };
    immCard.update();
    const redHub = (none.meanPeak - hub.meanPeak) / none.meanPeak;
    const redRand = (none.meanPeak - random.meanPeak) / none.meanPeak;
    clear(immMetrics);
    immMetrics.appendChild(metricGrid([
      { label: `Coverage`, value: pct(cov, 0), foot: "of nodes protected" },
      { label: "Hub-targeted peak", value: pct(hub.meanPeak), accent: true, foot: `${pct(redHub, 0)} reduction` },
      { label: "Random peak", value: pct(random.meanPeak), foot: `${pct(redRand, 0)} reduction` },
      { label: "Targeting advantage", value: `${fx((redHub - redRand) * 100, 0)} pts`, accent: true, foot: "extra peak reduction" },
    ]));
    immStatus.textContent = `Done · hub-targeting removes ${pct(redHub, 0)} of peak vs ${pct(redRand, 0)} at random.`;
    running = false;
  }

  // ---------- controls ----------
  const topoSel = selectControl("Topology", TOPO, cfg.topology, (v) => { cfg.topology = v; store.set({ network: { type: v } }, { silent: true }); buildAnimation(); });
  const nSel = selectControl("Nodes", { 200: "200", 400: "400", 600: "600", 800: "800", 1000: "1000" }, String(cfg.n), (v) => { cfg.n = +v; store.set({ network: { n: +v } }, { silent: true }); buildAnimation(); });
  const degSel = selectControl("Mean degree ⟨k⟩", { 4: "4", 6: "6", 8: "8", 10: "10", 12: "12", 16: "16" }, String(cfg.meanDegree), (v) => { cfg.meanDegree = +v; store.set({ network: { meanDegree: +v } }, { silent: true }); buildAnimation(); });
  const runsSel = selectControl("Runs / ensemble", { 20: "20", 40: "40", 60: "60", 80: "80" }, String(cfg.runs), (v) => { cfg.runs = +v; store.set({ network: { runs: +v } }, { silent: true }); });
  const covSel = selectControl("Immunise coverage", { 0: "off", 0.1: "10%", 0.2: "20%", 0.3: "30%", 0.4: "40%", 0.5: "50%" }, String(cfg.immunizeFrac), (v) => { cfg.immunizeFrac = +v; store.set({ network: { immunizeFrac: +v } }, { silent: true }); buildAnimation(); });

  el.appendChild(h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: "Screen 4 · Network Monte-Carlo" }),
    h("h1", { text: "Network Lab" }),
    h("p", { class: "prose", text: "Stochastic Reed–Frost contagion on explicit social graphs. Watch it spread across a rendered network, compare topologies over Monte-Carlo ensembles, and test whether protecting high-degree hubs beats random coverage at equal reach." }),
  ]));

  el.appendChild(panel({ title: "Configuration", children: [
    h("div", { class: "row", style: { gap: "18px", flexWrap: "wrap" } }, [topoSel, nSel, degSel, runsSel, covSel]),
    h("p", { class: "hint", style: { marginTop: "10px" }, html: `Per-contact β is scaled to the node degree so β·⟨k⟩ matches the mean-field calibration (≈ 0.12/mo); exogenous μ is set to 0 so pure social contagion drives the comparison. ${workerAvailable() ? "Ensembles run in a background worker; the UI stays responsive." : "Running on the main thread (open via a local server to enable the background worker)."}` }),
  ] }));

  el.appendChild(panel({ title: "Live contagion", sub: `one stochastic run on a ${VIS_N()}-node ${TOPO[cfg.topology].split(" (")[0]} graph`,
    badge: "conditional", children: [
      h("div", { class: "net-wrap" }, [canvas]),
      animStatus,
      tlHolder,
      legend([{ name: "S Susceptible", color: "var(--c-s)" }, { name: "A At-Risk", color: "var(--c-a)" }, { name: "I Addicted", color: "var(--c-i)" }, { name: "R Recovered", color: "var(--c-r)" }, { name: "protected", color: "var(--ink-3)" }]),
    ] }));

  const runTopoBtn = h("button", { class: "primary sm", onclick: runTopologies }, ["Run topology comparison"]);
  el.appendChild(panel({ title: "Topology comparison", sub: "scale-free peaks higher and earlier", actions: runTopoBtn, children: [
    topoMetrics, topoCard, topoStatus,
    legend(["ER", "WS", "BA"].map((k) => ({ name: TOPO[k], color: TOPO_COLOR[k] }))),
  ] }));

  const runImmBtn = h("button", { class: "primary sm", onclick: runImmunisation }, ["Run immunisation comparison"]);
  el.appendChild(panel({ title: "Hub vs random immunisation", sub: "equal coverage, scale-free graph", actions: runImmBtn, children: [
    immMetrics, immCard, immStatus,
  ] }));

  // init
  requestAnimationFrame(buildAnimation);
  const onResize = () => { if (anim) { sizeCanvas(); drawFrame(0); } };
  window.addEventListener("resize", onResize);

  return { el, destroy: () => window.removeEventListener("resize", onResize) };
}

function selectControl(label, opts, value, onchange) {
  const sel = h("select", { "aria-label": label, onchange: (e) => onchange(e.target.value) },
    Object.entries(opts).map(([v, t]) => h("option", { value: v, selected: v === value, text: t })));
  return h("label", { class: "grp", style: { flexDirection: "column", alignItems: "flex-start", gap: "3px", fontSize: "12px", color: "var(--ink-3)" } }, [label, sel]);
}
function mean(a) { return a.reduce((x, y) => x + y, 0) / a.length; }
function placeholderSVG(text) {
  const svg = document.createElementNS("http://www.w3.org/2000/svg", "svg");
  svg.setAttribute("viewBox", "0 0 720 260");
  const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
  t.setAttribute("x", 360); t.setAttribute("y", 130); t.setAttribute("text-anchor", "middle");
  t.setAttribute("fill", "var(--ink-3)"); t.setAttribute("font-size", 14); t.textContent = text;
  svg.appendChild(t); return svg;
}
