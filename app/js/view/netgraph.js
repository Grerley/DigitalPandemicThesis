// netgraph.js — Force-directed layout + canvas rendering of a contagion graph.
// -----------------------------------------------------------------------------
// A compact, seeded force layout (repulsion + spring attraction + centring)
// computed once, then fast canvas redraws for each animation frame. Node colour
// encodes state (S/A/I/R/protected). Kept small (≤ ~260 nodes) for legibility
// and smooth animation; ensembles run at full size elsewhere.

import { RNG } from "../model/rng.js";
import { STATE } from "../model/network.js";

const STATE_COLORS = {
  css: ["--c-s", "--c-a", "--c-i", "--c-r", "--ink-3"],
};

/** Compute a 2D layout for adjacency `adj` within [0,w]×[0,h]. Deterministic. */
export function layout(adj, w, h, seed = 1, iters = 160) {
  const n = adj.length;
  const rng = new RNG(seed);
  const pos = new Array(n);
  const cx = w / 2, cy = h / 2, rad = Math.min(w, h) * 0.42;
  for (let i = 0; i < n; i++) {
    const a = rng.next() * Math.PI * 2, r = Math.sqrt(rng.next()) * rad;
    pos[i] = { x: cx + Math.cos(a) * r, y: cy + Math.sin(a) * r, dx: 0, dy: 0 };
  }
  const area = w * h;
  const kk = 0.9 * Math.sqrt(area / Math.max(1, n));
  let temp = w * 0.12;
  for (let it = 0; it < iters; it++) {
    for (let i = 0; i < n; i++) { pos[i].dx = 0; pos[i].dy = 0; }
    // repulsion (O(n^2), fine for small n)
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        let dx = pos[i].x - pos[j].x, dy = pos[i].y - pos[j].y;
        let d2 = dx * dx + dy * dy;
        if (d2 < 0.01) { dx = (rng.next() - 0.5); dy = (rng.next() - 0.5); d2 = 0.01; }
        const f = (kk * kk) / d2;
        const fx = dx * f, fy = dy * f;
        pos[i].dx += fx; pos[i].dy += fy;
        pos[j].dx -= fx; pos[j].dy -= fy;
      }
    }
    // attraction along edges
    for (let i = 0; i < n; i++) {
      for (const j of adj[i]) {
        if (j <= i) continue;
        const dx = pos[i].x - pos[j].x, dy = pos[i].y - pos[j].y;
        const d = Math.hypot(dx, dy) || 0.01;
        const f = (d * d) / kk / 12;
        const fx = (dx / d) * f, fy = (dy / d) * f;
        pos[i].dx -= fx; pos[i].dy -= fy;
        pos[j].dx += fx; pos[j].dy += fy;
      }
    }
    // centring + integrate with cooling
    for (let i = 0; i < n; i++) {
      pos[i].dx += (cx - pos[i].x) * 0.012;
      pos[i].dy += (cy - pos[i].y) * 0.012;
      const d = Math.hypot(pos[i].dx, pos[i].dy) || 1;
      const step = Math.min(d, temp);
      pos[i].x += (pos[i].dx / d) * step;
      pos[i].y += (pos[i].dy / d) * step;
      pos[i].x = Math.max(6, Math.min(w - 6, pos[i].x));
      pos[i].y = Math.max(6, Math.min(h - 6, pos[i].y));
    }
    temp *= 0.975;
  }
  return pos;
}

/** Resolve a CSS custom property to a concrete colour for canvas. */
function resolveColors() {
  const cs = getComputedStyle(document.documentElement);
  return STATE_COLORS.css.map((v) => cs.getPropertyValue(v).trim() || "#888");
}

/**
 * Draw the graph on a canvas.
 * @param {CanvasRenderingContext2D} ctx
 * @param {{adj:number[][],degree:number[]}} g
 * @param {Array} pos layout positions
 * @param {Int8Array} states node states for this frame
 * @param {object} o { w, h, dpr }
 */
export function drawGraph(ctx, g, pos, states, o) {
  const { w, h, dpr = 1 } = o;
  const colors = resolveColors();
  ctx.save();
  ctx.scale(dpr, dpr);
  ctx.clearRect(0, 0, w, h);
  // edges
  ctx.lineWidth = 0.5;
  ctx.strokeStyle = "rgba(130,130,130,0.22)";
  ctx.beginPath();
  const n = g.adj.length;
  for (let i = 0; i < n; i++) {
    for (const j of g.adj[i]) {
      if (j <= i) continue;
      ctx.moveTo(pos[i].x, pos[i].y);
      ctx.lineTo(pos[j].x, pos[j].y);
    }
  }
  ctx.stroke();
  // nodes sized by degree
  const maxDeg = Math.max(1, ...g.degree);
  for (let i = 0; i < n; i++) {
    const st = states[i];
    const r = 2.2 + 4.5 * Math.sqrt(g.degree[i] / maxDeg);
    ctx.beginPath();
    ctx.arc(pos[i].x, pos[i].y, r, 0, Math.PI * 2);
    ctx.fillStyle = colors[st];
    ctx.globalAlpha = st === STATE.V ? 0.45 : 1;
    ctx.fill();
    if (st === STATE.I) { ctx.globalAlpha = 0.9; ctx.lineWidth = 1; ctx.strokeStyle = colors[2]; ctx.stroke(); }
    ctx.globalAlpha = 1;
  }
  ctx.restore();
}
