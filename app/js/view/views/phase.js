// phase.js — Screen 3: Phase & Equilibrium.
// -----------------------------------------------------------------------------

import { h, clear } from "../dom.js";
import { panel, metricGrid, chartCard, legend, condNote } from "../ui.js";
import { timeSeries, bifurcation, phasePortrait } from "../charts.js";
import { computeScenario } from "../scenario.js";
import { bifurcationCurve, endemicEquilibrium, r0Full } from "../../model/r0.js";
import { transitionMatrix } from "../../model/sair.js";
import { baselineInit } from "../../model/params.js";
import { pct, fx } from "../format.js";
import { store } from "../../state.js";

export function PhaseView() {
  const el = h("div", { class: "view" });
  const metricsWrap = h("div");

  // Vector field for the (A, I) phase plane from the discrete map.
  function vectorField(p) {
    const field = [];
    const nA = 8, nI = 8;
    const Amax = 0.32, Imax = 0.20;
    for (let ia = 1; ia <= nA; ia++) {
      for (let ii = 1; ii <= nI; ii++) {
        const A = (Amax * ia) / (nA + 1), I = (Imax * ii) / (nI + 1);
        const R = Math.min(0.8, (p.gamma2 * I) / (p.delta + 1e-9));
        const S = Math.max(0, 1 - A - I - R);
        const sum = S + A + I + R || 1;
        const st = { S: S / sum, A: A / sum, I: I / sum, R: R / sum };
        const { M } = transitionMatrix(st, p);
        const nA_ = st.S * M[0][1] + st.A * M[1][1] + st.I * M[2][1] + st.R * M[3][1];
        const nI_ = st.A * M[1][2] + st.I * M[2][2];
        field.push({ x: A, y: I, dx: nA_ - st.A, dy: nI_ - st.I });
      }
    }
    return field;
  }

  const bifCard = chartCard({
    title: "Transcritical bifurcation", sub: "endemic prevalence I* vs basic R₀",
    filename: "bifurcation",
    csv: () => { const bc = bifurcationCurve(store.get().params); return bc.r0.map((r, i) => ({ R0: fx(r, 3), I_star: fx(bc.iStar[i], 5) })); },
    render: () => bifurcation(bifurcationCurve(store.get().params)),
  });

  const phaseCard = chartCard({
    title: "Phase portrait", sub: "At-Risk (A) vs Addicted (I)",
    filename: "phase-portrait",
    render: () => {
      const sc = computeScenario(store.get());
      const eq = endemicEquilibrium(sc.params, baselineInit(), 400);
      return phasePortrait({ traj: sc.traj, eq, xKey: "A", yKey: "I", field: vectorField(sc.params) });
    },
  });

  const lagCard = chartCard({
    title: "Incidence leads prevalence", sub: "the prevention window closes early",
    filename: "incidence-prevalence-lag",
    csv: () => { const sc = computeScenario(store.get()); return sc.traj.year.map((y, i) => ({ year: fx(y, 3), incidence: fx(sc.traj.incidence[i], 6), prevalence_I: fx(sc.traj.I[i], 5) })); },
    render: () => {
      const sc = computeScenario(store.get());
      const inc = sc.traj.incidence, I = sc.traj.I;
      const incMax = Math.max(...inc) || 1, IMax = Math.max(...I) || 1;
      const incN = inc.map((v) => v / incMax), IN = I.map((v) => v / IMax);
      const incPeakY = sc.traj.year[argmax(inc)], IPeakY = sc.traj.year[argmax(I)];
      return timeSeries({
        xs: sc.traj.year,
        series: [
          { key: "inc", name: "incidence (new S→A)", color: "var(--c-a)", ys: incN, dash: "6 3" },
          { key: "prev", name: "prevalence (I)", color: "var(--c-i)", ys: IN },
        ],
        yTitle: "share of own peak", yMax: 1.1, yFmt: (v) => `${Math.round(v * 100)}%`,
        markers: [
          { x: incPeakY, label: `incidence peak`, color: "var(--c-a)" },
          { x: IPeakY, label: `prevalence peak`, color: "var(--c-i)" },
        ],
      });
    },
  });

  function renderMetrics() {
    const sc = computeScenario(store.get());
    const p = sc.params;
    const eq = endemicEquilibrium(p, baselineInit(), 400);
    const incPeakY = sc.traj.year[argmax(sc.traj.incidence)];
    const IPeakY = sc.traj.year[argmax(sc.traj.I)];
    clear(metricsWrap);
    metricsWrap.appendChild(metricGrid([
      { label: "R₀ full (threshold)", value: fx(r0Full(p, 1), 2), foot: "governs DFE stability", accent: true },
      { label: "Endemic S*", value: pct(eq.S) },
      { label: "Endemic A*", value: pct(eq.A) },
      { label: "Endemic I*", value: pct(eq.I), accent: true },
      { label: "Endemic R*", value: pct(eq.R) },
      { label: "Incidence lead", value: `${fx(IPeakY - incPeakY, 1)} yr`, foot: "prevalence peak − incidence peak" },
    ]));
  }

  el.appendChild(h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: "Screen 3 · Dynamical structure" }),
    h("h1", { text: "Phase & Equilibrium" }),
    h("p", { class: "prose", text: "The threshold structure of the model. Below R₀ = 1 the disease-free state is stable; above it a positive endemic equilibrium emerges through a transcritical bifurcation. The phase portrait traces the state's approach to that equilibrium, and the lag view shows incidence peaking years before prevalence." }),
  ]));
  el.appendChild(metricsWrap);
  el.appendChild(panel({ title: "Bifurcation diagram", sub: "swept by scaling β so R₀ ranges over [0.2, 4]; μ set to 0 to isolate the endogenous threshold",
    badge: "conditional", children: [bifCard] }));
  el.appendChild(h("div", { class: "split" }, [
    panel({ title: "Phase portrait", children: [phaseCard,
      h("p", { class: "chart-cap", text: "Arrows show the one-step direction of the map; the trajectory spirals into the endemic equilibrium (●)." })] }),
    panel({ title: "Incidence vs prevalence", children: [lagCard,
      legend([{ name: "incidence (new S→A)", color: "var(--c-a)", dash: "6 3" }, { name: "prevalence (I)", color: "var(--c-i)" }])] }),
  ]));
  el.appendChild(condNote("The ≈3-year lead of incidence over prevalence implies interventions timed to observed addiction rates arrive after the transmission wave has largely passed."));

  renderMetrics();
  const unsub = store.subscribe(() => { bifCard.update(); phaseCard.update(); lagCard.update(); renderMetrics(); });
  return { el, destroy: unsub };
}

function argmax(arr) { let bi = 0; for (let i = 1; i < arr.length; i++) if (arr[i] > arr[bi]) bi = i; return bi; }
