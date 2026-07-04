// diagram.js — Animated S–A–I–R state-transition diagram.
// -----------------------------------------------------------------------------
// Four compartment nodes with directed edges whose thickness tracks the live
// per-step transition probabilities. Animated dash offset conveys flow
// direction (disabled under prefers-reduced-motion). Node fill encodes current
// occupancy. Pure of app state — driven by (state, params) passed in.

import { s } from "./dom.js";
import { transitionMatrix } from "../model/sair.js";
import { reduceMotion } from "./dom.js";

const NODES = {
  S: { x: 90, y: 80, label: "S", full: "Susceptible", color: "var(--c-s)" },
  A: { x: 300, y: 80, label: "A", full: "At-Risk", color: "var(--c-a)" },
  I: { x: 300, y: 230, label: "I", full: "Addicted", color: "var(--c-i)" },
  R: { x: 90, y: 230, label: "R", full: "Recovered", color: "var(--c-r)" },
};

// edges: from,to, probability key, curvature, label position
const EDGES = [
  { from: "S", to: "A", key: "P_SA", name: "β·k·(I+ηA)", curve: -28 },
  { from: "A", to: "I", key: "P_AI", name: "σ", curve: 0 },
  { from: "I", to: "R", key: "P_IR", name: "γ₂", curve: -28 },
  { from: "R", to: "A", key: "P_RA", name: "δ (relapse)", curve: 0, cross: true },
  { from: "A", to: "S", key: "P_AS", name: "γ₁", curve: 28 },
];

const R = 34;

export function stateDiagram({ state, params, w = 400, h = 310 }) {
  const svg = s("svg", { viewBox: `0 0 ${w} ${h}`, role: "img" });
  svg.appendChild(s("title", { text: "S–A–I–R state-transition diagram with live flow rates" }));
  const { probs } = transitionMatrix(state, params);
  const maxP = Math.max(probs.P_SA, probs.P_AI, probs.P_AS, probs.P_IR, probs.P_RA, 1e-6);
  const defs = s("defs");
  const anim = !reduceMotion();

  const edgeG = s("g");
  for (const e of EDGES) {
    const a = NODES[e.from], b = NODES[e.to];
    const p = probs[e.key];
    const width = 1.2 + 7 * Math.sqrt(p / maxP);
    const { path, mid, ang } = edgePath(a, b, e.curve, e.cross);
    const mkId = `arw-${e.from}${e.to}`;
    const mk = s("marker", { id: mkId, viewBox: "0 0 8 8", refX: 6, refY: 4, markerWidth: 6, markerHeight: 6, orient: "auto-start-reverse" });
    mk.appendChild(s("path", { d: "M0,0L8,4L0,8Z", fill: a.color, opacity: 0.9 }));
    defs.appendChild(mk);
    // base line
    edgeG.appendChild(s("path", { d: path, fill: "none", stroke: a.color, "stroke-width": width, opacity: 0.28, "stroke-linecap": "round", "marker-end": `url(#${mkId})` }));
    // animated flow dashes
    const flow = s("path", { d: path, fill: "none", stroke: a.color, "stroke-width": Math.max(1, width - 2), "stroke-linecap": "round",
      "stroke-dasharray": "1.5 10", opacity: 0.95 });
    if (anim) {
      const speed = Math.max(0.6, 3.2 * (p / maxP));
      const an = s("animate", { attributeName: "stroke-dashoffset", from: "0", to: "-11.5", dur: `${(1 / speed).toFixed(2)}s`, repeatCount: "indefinite" });
      flow.appendChild(an);
    }
    edgeG.appendChild(flow);
    // label
    edgeG.appendChild(s("text", { x: mid.x, y: mid.y - 3, "text-anchor": "middle", class: "serieslabel", fill: a.color,
      "font-size": 10.5, text: `${e.name}` }));
    edgeG.appendChild(s("text", { x: mid.x, y: mid.y + 9, "text-anchor": "middle", class: "axis", fill: "var(--ink-3)",
      "font-size": 9.5, text: `${(p * 100).toFixed(p < 0.01 ? 2 : 1)}%/mo` }));
  }
  svg.appendChild(defs);
  svg.appendChild(edgeG);

  // nodes
  for (const key of ["S", "A", "I", "R"]) {
    const n = NODES[key];
    const occ = state[key] ?? 0;
    const g = s("g");
    g.appendChild(s("circle", { cx: n.x, cy: n.y, r: R, fill: "var(--panel)", stroke: n.color, "stroke-width": 2 }));
    // occupancy ring fill
    g.appendChild(s("circle", { cx: n.x, cy: n.y, r: R, fill: n.color, "fill-opacity": 0.1 + 0.55 * occ }));
    g.appendChild(s("text", { x: n.x, y: n.y - 2, "text-anchor": "middle", "font-size": 22, "font-weight": 700, fill: n.color, "font-family": "var(--serif)", "font-style": "italic", text: n.label }));
    g.appendChild(s("text", { x: n.x, y: n.y + 14, "text-anchor": "middle", "font-size": 10.5, "font-variant-numeric": "tabular-nums", fill: "var(--ink-2)", text: `${(occ * 100).toFixed(1)}%` }));
    g.appendChild(s("text", { x: n.x, y: n.y + R + 15, "text-anchor": "middle", class: "axis", "font-size": 10.5, fill: "var(--ink-3)", text: n.full }));
    svg.appendChild(g);
  }
  return svg;
}

function edgePath(a, b, curve, cross) {
  const dx = b.x - a.x, dy = b.y - a.y;
  const len = Math.hypot(dx, dy);
  const ux = dx / len, uy = dy / len;
  // start/end on circle boundaries
  const sx = a.x + ux * R, sy = a.y + uy * R;
  const ex = b.x - ux * R, ey = b.y - uy * R;
  // control point offset perpendicular
  const nx = -uy, ny = ux;
  const cx = (sx + ex) / 2 + nx * curve, cy = (sy + ey) / 2 + ny * curve;
  const path = `M${sx},${sy}Q${cx},${cy} ${ex},${ey}`;
  const mid = { x: 0.25 * sx + 0.5 * cx + 0.25 * ex, y: 0.25 * sy + 0.5 * cy + 0.25 * ey };
  return { path, mid, ang: Math.atan2(dy, dx) };
}
