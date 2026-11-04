// networks.js: Screen 5, the model on explicit contact networks (thesis Chapter 8).

import { h } from "../dom.js";
import { panel, metricGrid, chartCard, legend, timeline } from "../ui.js";
import { timeSeries } from "../charts.js";
import { layout, drawGraph } from "../netgraph.js";
import { viewHead, table, source, pc, COLORS } from "../common.js";
import { R, LAB, centralParams } from "../../model/params.js";
import { simulate } from "../../model/sair.js";
import { makeGraph, simulateNetwork, betaFromB, Rnet, GRAPH_LABELS } from "../../model/network.js";
import { runJob } from "../ensemble-runner.js";

const STRATS = { none: "No protection", random: "Random young people", degree: "Best-connected young people", acquaintance: "Nominated friends (acquaintance)" };

export function NetworksView() {
  const el = h("div", { class: "view" });
  el.appendChild(viewHead("Screen 5 · Networks", "Does who is connected to whom matter?",
    "The calibrated model is placed on three families of random networks with the same mean number of contacts (10): Erdős–Rényi (random), Watts–Strogatz (clustered) and Barabási–Albert (a few highly connected young people). " +
    "Time runs in monthly steps; a Susceptible young person is exposed through At-Risk and problematic-use contacts and through the non-social route. Protection reduces a young person's susceptibility to social exposure to 20%."));

  const pc1 = centralParams();
  const x22 = simulate(pc1, LAB.central.x18, [4])[0];

  // ---- controls
  const cfg = { kind: "BA", n: 1000, regime: "calibrated", strategy: "degree", cover: 0.1, reps: 12 };
  const select = (label, key, options) => h("label", { class: "hint", style: { display: "flex", flexDirection: "column", gap: "4px" } }, [label,
    h("select", { onchange: (e) => { cfg[key] = isNaN(+e.target.value) ? e.target.value : +e.target.value; } },
      Object.entries(options).map(([v, t]) => h("option", { value: v, text: t, selected: String(cfg[key]) === v })))]);
  const status = h("span", { class: "hint" });
  const result = h("div");
  const runBtn = h("button", { class: "primary sm", text: "Run experiment", onclick: run });

  async function run() {
    runBtn.disabled = true; status.textContent = "running…";
    const res = await runJob({ ...cfg, k: 10, years: 20, params: pc1, x22, damp: 0.2, seed: 4242 },
      (f) => { status.textContent = `running… ${Math.round(f * 100)}%`; });
    runBtn.disabled = false; status.textContent = "";
    const yrs = Array.from({ length: res.base.mean.length }, (_, i) => 2022 + i);
    const series = [{ key: "b", name: "No protection", color: "var(--ink-2)", ys: res.base.mean, band: { lo: res.base.lo, hi: res.base.hi } }];
    if (res.prot) series.push({ key: "p", name: "With protection", color: COLORS.I, dash: "6 4", ys: res.prot.mean, band: { lo: res.prot.lo, hi: res.prot.hi } });
    const b42 = res.base.mean.at(-1), p42 = res.prot ? res.prot.mean.at(-1) : null;
    result.replaceChildren(
      metricGrid([
        { label: "Network reproduction number", value: res.rnet.toFixed(2), foot: "heuristic R_net, mean over graphs", accent: true },
        { label: "Exposure rate used", value: res.bUsed.toFixed(2), foot: cfg.regime === "near" ? "set so that R_net = 1.5" : "calibrated, after 2018" },
        { label: "2042, no protection", value: pc(b42), foot: `mean of ${cfg.reps} runs` },
        p42 != null ? { label: "2042, with protection", value: pc(p42), foot: `${((1 - p42 / b42) * 100).toFixed(0)}% lower` } : { label: "Protection", value: "none", foot: "" },
      ]),
      chartCard({ title: `${GRAPH_LABELS[cfg.kind]} network, ${cfg.n.toLocaleString("en-GB")} young people`, sub: "mean and 5th to 95th percentile across runs (live)", filename: "network-experiment",
        render: () => timeSeries({ xs: yrs, series, yTitle: "Problematic use", yMax: Math.max(...series.flatMap((s) => s.band.hi)) * 1.15 }) }));
  }

  el.appendChild(panel({ title: "Run a network experiment", sub: "each run uses the same graph and random numbers with and without protection, as in the thesis", children: [
    h("div", { class: "row", style: { flexWrap: "wrap", gap: "14px", alignItems: "flex-end", marginBottom: "12px" } }, [
      select("Network", "kind", { ER: "Erdős–Rényi", WS: "Watts–Strogatz", BA: "Barabási–Albert" }),
      select("Young people", "n", { 500: "500", 1000: "1,000", 2000: "2,000" }),
      select("Exposure", "regime", { calibrated: "Calibrated (far above threshold)", near: "Near threshold (R_net = 1.5)" }),
      select("Protect 10% of", "strategy", STRATS),
      select("Runs", "reps", { 6: "6", 12: "12", 24: "24" }),
      runBtn, status,
    ]),
    result,
  ] }));

  // ---- small animated network
  const W = 640, H = 420, N = 220;
  const canvas = h("canvas", { width: W * 2, height: H * 2, style: { width: "100%", maxWidth: `${W}px` } });
  const anim = { kind: "BA" };
  const animBox = h("div");
  function buildAnim() {
    const adj = makeGraph(anim.kind, N, 6, 7);
    const g = { adj, degree: adj.map((a) => a.length) };
    const pos = layout(adj, W, H, 3, 140);
    const p = { ...pc1, b: pc1.b * 2 }, beta = betaFromB(p, 6);
    const { snapshots, frac } = simulateNetwork(adj, p, beta, { months: 120, x0: [0.97, 0.0, 0.03, 0], seed: 11, recordEvery: 1, keepStates: true });
    const ctx = canvas.getContext("2d");
    const draw = (t) => {
      const st = Int8Array.from(snapshots[t]);
      drawGraph(ctx, g, pos, st, { w: W, h: H, dpr: 2 });
    };
    const tl = timeline({ max: snapshots.length - 1, value: 0, onScrub: draw, labelFor: (t) => `month ${t} · problematic use ${pc(frac[t][2])}`, fps: 12 });
    animBox.replaceChildren(tl.el);
    draw(0);
  }
  el.appendChild(panel({ title: "Watch spread on a small network", sub: `${N} young people, mean degree 6, exposure doubled and 3% starting with problematic use so the spread is visible (illustration only)`, children: [
    h("div", { class: "row", style: { gap: "10px", marginBottom: "8px" } }, Object.entries(GRAPH_LABELS).map(([k, t]) =>
      h("button", { class: "sm", text: t, onclick: () => { anim.kind = k; buildAnim(); } }))),
    h("div", { class: "net-wrap" }, canvas), animBox,
    legend([{ name: "Susceptible", color: COLORS.S }, { name: "At-Risk", color: COLORS.A }, { name: "Problematic use", color: COLORS.I }, { name: "Recovered", color: COLORS.R }]),
  ] }));
  setTimeout(buildAnim, 30);

  // ---- thesis results
  const kinds = ["ER", "WS", "BA"];
  el.appendChild(panel({ title: "Thesis results", children: [
    table(["Network", "⟨k²⟩/⟨k⟩", "R_net", "Problematic use 2042, mean (range)"],
      kinds.map((k) => [GRAPH_LABELS[k], R(`net_${k}_k2_over_k`).toFixed(1), R(`net_${k}_R0net`).toFixed(2),
        `${R(`net_${k}_I2042_mean`).toFixed(1)}% (${R(`net_${k}_I2042_p05`).toFixed(1)}% to ${R(`net_${k}_I2042_p95`).toFixed(1)}%)`])
        .concat([["Population model", "", R("net_meanfield_R0").toFixed(2), `${R("net_meanfield_I2042").toFixed(1)}%`]]), ["l", "r", "r", "r"]),
    h("p", { class: "prose", html: "Degree heterogeneity changes the threshold strongly but the projected level little. Far above the threshold, networks lower prevalence slightly through an exposure deficit." }),
    table(["Protect 10% of young people", "Calibrated exposure", "Near threshold (R_net = 1.5)"],
      [["Best-connected", "hub"], ["Random", "random"], ["Nominated friends", "acq"]].map(([t, k]) =>
        [t, `${R(`netx_tgt_cal_${k}_c10_rel_pct`).toFixed(0)}% lower`, `${R(`netx_tgt_r1p5_${k}_c10_rel_pct`).toFixed(0)}% lower`]), ["l", "r", "r"]),
    h("p", { class: "prose", html: "Targeting the best-connected adds little far above the threshold and a great deal near it. Protecting friends nominated by randomly chosen young people recovers about half of that advantage near the threshold without knowing the network." }),
    h("div", { class: "caveat", html: "<strong>Friendship data cannot separate contagion from homophily.</strong> In simulated two-wave surveys, a world with no contagion but a modest latent trait shared by friends gives the same naive peer-effect estimate as contagion. Targeting helps only if the exposure rate reflects influence." }),
    source("thesis Chapter 8 and the networks paper (Mutibura and Lubuma, manuscript in preparation); python/run_all.py and python/ext_networks.py."),
  ] }));
  return { el };
}
