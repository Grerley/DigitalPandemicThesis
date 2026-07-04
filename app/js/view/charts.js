// charts.js — Reusable, export-ready SVG chart primitives.
// -----------------------------------------------------------------------------
// Every chart returns a self-contained <svg> element with a fixed viewBox
// (scaled to its container by CSS). Charts favour direct labelling over
// legends, colour-blind-safe encodings paired with line style, crisp axes and
// restrained motion. All are pure functions of their data + options.

import { s } from "./dom.js";
import { fx } from "./format.js";

// ---- scales & ticks ---------------------------------------------------------
export function linScale([d0, d1], [r0, r1]) {
  const m = d1 === d0 ? 0 : (r1 - r0) / (d1 - d0);
  const f = (v) => r0 + (v - d0) * m;
  f.invert = (p) => d0 + (p - r0) / m;
  f.domain = [d0, d1]; f.range = [r0, r1];
  return f;
}
export function logScale([d0, d1], [r0, r1]) {
  const l0 = Math.log10(d0), l1 = Math.log10(d1);
  const f = (v) => r0 + ((Math.log10(v) - l0) / (l1 - l0)) * (r1 - r0);
  f.domain = [d0, d1]; f.range = [r0, r1];
  return f;
}
export function niceTicks(min, max, count = 5) {
  const span = max - min || 1;
  const step0 = span / count;
  const mag = Math.pow(10, Math.floor(Math.log10(step0)));
  const norm = step0 / mag;
  const step = (norm < 1.5 ? 1 : norm < 3 ? 2 : norm < 7 ? 5 : 10) * mag;
  const start = Math.ceil(min / step) * step;
  const ticks = [];
  for (let v = start; v <= max + step * 1e-6; v += step) ticks.push(Math.round(v / step) * step);
  return ticks;
}

// ---- frame ------------------------------------------------------------------
function frame({ w = 720, h = 420, m = { t: 22, r: 74, b: 44, l: 56 }, title } = {}) {
  const svg = s("svg", { viewBox: `0 0 ${w} ${h}`, role: "img", preserveAspectRatio: "xMidYMid meet" });
  const g = s("g", { transform: `translate(${m.l},${m.t})` });
  svg.appendChild(g);
  const iw = w - m.l - m.r, ih = h - m.t - m.b;
  if (title) svg.appendChild(s("title", { text: title }));
  return { svg, g, iw, ih, m, w, h };
}

function axes(g, x, y, iw, ih, { xTitle, yTitle, xTicks, yTicks, yFmt = (v) => fx(v, 0), xFmt = (v) => fx(v, 0), yGrid = true } = {}) {
  const gridG = s("g", { class: "grid" });
  const axG = s("g", { class: "axis" });
  // y ticks + grid
  for (const t of yTicks) {
    const yy = y(t);
    if (yGrid) gridG.appendChild(s("line", { x1: 0, x2: iw, y1: yy, y2: yy }));
    axG.appendChild(s("text", { x: -9, y: yy + 3.5, "text-anchor": "end", text: yFmt(t) }));
  }
  // x ticks
  for (const t of xTicks) {
    const xx = x(t);
    axG.appendChild(s("line", { class: "dom", x1: xx, x2: xx, y1: ih, y2: ih + 5, stroke: "currentColor" }));
    axG.appendChild(s("text", { x: xx, y: ih + 18, "text-anchor": "middle", text: xFmt(t) }));
  }
  // domain line
  axG.appendChild(s("path", { d: `M0,${ih}H${iw}`, fill: "none", stroke: "currentColor" }));
  if (xTitle) axG.appendChild(s("text", { class: "axis-title", x: iw / 2, y: ih + 38, "text-anchor": "middle", text: xTitle }));
  if (yTitle) axG.appendChild(s("text", { class: "axis-title", transform: `translate(${-42},${ih / 2}) rotate(-90)`, "text-anchor": "middle", text: yTitle }));
  g.appendChild(gridG); g.appendChild(axG);
  return axG;
}

function linePath(xs, ys, x, y) {
  let d = "";
  for (let i = 0; i < xs.length; i++) d += `${i ? "L" : "M"}${x(xs[i]).toFixed(2)},${y(ys[i]).toFixed(2)}`;
  return d;
}
function areaPath(xs, lo, hi, x, y) {
  let d = "";
  for (let i = 0; i < xs.length; i++) d += `${i ? "L" : "M"}${x(xs[i]).toFixed(2)},${y(hi[i]).toFixed(2)}`;
  for (let i = xs.length - 1; i >= 0; i--) d += `L${x(xs[i]).toFixed(2)},${y(lo[i]).toFixed(2)}`;
  return d + "Z";
}

/**
 * Time-series / trajectory chart. The workhorse for most screens.
 * @param {object} cfg
 *   xs: number[] (shared x, e.g. years)
 *   series: [{ key,name,color,dash,ys, width?, labelSide?, band?:{lo,hi}, opacity? }]
 *   xTitle,yTitle, yMax?, yFmt?, xFmt?, mode: 'line'|'stacked'
 *   markers: [{ x, label, color }]  vertical reference lines
 *   scrubX: number|null  vertical time cursor
 *   points: [{x,y,color,label}]  emphasised points
 *   directLabels: boolean (default true)
 */
export function timeSeries(cfg) {
  const {
    xs, series, xTitle = "Year", yTitle = "Proportion", yMax = null, yMin = 0,
    yFmt = (v) => `${Math.round(v * 100)}%`, xFmt = (v) => fx(v, 0), mode = "line",
    markers = [], scrubX = null, points = [], directLabels = true, w = 720, h = 420, xTickCount = 6, yTickCount = 5,
  } = cfg;
  const { svg, g, iw, ih } = frame({ w, h });
  const xExtent = [xs[0], xs[xs.length - 1]];
  let ymax = yMax;
  if (ymax == null) {
    ymax = 0;
    if (mode === "stacked") ymax = 1;
    else for (const se of series) {
      const arr = se.band ? se.band.hi : se.ys;
      for (const v of arr) if (v > ymax) ymax = v;
    }
    ymax *= 1.08;
  }
  const x = linScale(xExtent, [0, iw]);
  const y = linScale([yMin, ymax], [ih, 0]);
  const yTicks = niceTicks(yMin, ymax, yTickCount);
  const xTicks = niceTicks(xExtent[0], xExtent[1], xTickCount);
  axes(g, x, y, iw, ih, { xTitle, yTitle, xTicks, yTicks, yFmt, xFmt });

  if (mode === "stacked") {
    // stacked areas bottom→top in series order
    const cum = new Array(xs.length).fill(0);
    for (const se of series) {
      const lo = cum.slice();
      for (let i = 0; i < xs.length; i++) cum[i] += se.ys[i];
      g.appendChild(s("path", { d: areaPath(xs, lo, cum, x, y), fill: se.color, "fill-opacity": 0.85, stroke: "none" }));
    }
    // thin separators
    const cum2 = new Array(xs.length).fill(0);
    for (const se of series) {
      for (let i = 0; i < xs.length; i++) cum2[i] += se.ys[i];
      g.appendChild(s("path", { d: linePath(xs, cum2.slice(), x, y), fill: "none", stroke: "var(--panel)", "stroke-width": 0.8, opacity: 0.6 }));
    }
  } else {
    // envelope bands first
    for (const se of series) {
      if (se.band) g.appendChild(s("path", { d: areaPath(xs, se.band.lo, se.band.hi, x, y), fill: se.color, "fill-opacity": 0.14, stroke: "none" }));
    }
    for (const se of series) {
      g.appendChild(s("path", {
        d: linePath(xs, se.ys, x, y), fill: "none", stroke: se.color,
        "stroke-width": se.width || 2, "stroke-dasharray": se.dash && se.dash !== "none" ? se.dash : null,
        "stroke-linejoin": "round", "stroke-linecap": "round", opacity: se.opacity ?? 1,
      }));
    }
  }

  // markers (vertical reference lines)
  for (const mk of markers) {
    const xx = x(mk.x);
    g.appendChild(s("line", { x1: xx, x2: xx, y1: 0, y2: ih, stroke: mk.color || "var(--ink-3)", "stroke-width": 1, "stroke-dasharray": "3 3", opacity: 0.8 }));
    if (mk.label) g.appendChild(s("text", { x: xx + 4, y: 11, class: "serieslabel", fill: mk.color || "var(--ink-3)", text: mk.label }));
  }
  // emphasised points
  for (const pt of points) {
    g.appendChild(s("circle", { cx: x(pt.x), cy: y(pt.y), r: 4, fill: pt.color || "var(--ink)", stroke: "var(--panel)", "stroke-width": 1.5 }));
    if (pt.label) g.appendChild(s("text", { x: x(pt.x), y: y(pt.y) - 9, "text-anchor": "middle", class: "serieslabel", fill: pt.color || "var(--ink)", text: pt.label }));
  }
  // scrub cursor
  if (scrubX != null) {
    const xx = x(scrubX);
    g.appendChild(s("line", { x1: xx, x2: xx, y1: 0, y2: ih, stroke: "var(--accent)", "stroke-width": 1.5, opacity: 0.9 }));
    for (const se of series) {
      const i = nearestIndex(xs, scrubX);
      let yv = se.ys[i];
      if (mode === "stacked") { /* dots omitted in stacked */ }
      else g.appendChild(s("circle", { cx: xx, cy: y(yv), r: 3.5, fill: se.color, stroke: "var(--panel)", "stroke-width": 1.2 }));
    }
  }
  // direct labels at right end
  if (directLabels && mode !== "stacked") {
    const placed = [];
    for (const se of series.filter((z) => z.label !== false)) {
      let yy = y(se.ys[se.ys.length - 1]);
      for (const p of placed) if (Math.abs(p - yy) < 12) yy = p + 12;
      placed.push(yy);
      g.appendChild(s("text", { x: iw + 6, y: yy + 3.5, class: "serieslabel", fill: se.color, text: se.name }));
    }
  }
  return svg;
}

function nearestIndex(xs, v) {
  let best = 0, bd = Infinity;
  for (let i = 0; i < xs.length; i++) { const d = Math.abs(xs[i] - v); if (d < bd) { bd = d; best = i; } }
  return best;
}

// ---- bifurcation ------------------------------------------------------------
export function bifurcation({ r0, iStar, threshold = 1, w = 720, h = 380 }) {
  const { svg, g, iw, ih } = frame({ w, h });
  const x = linScale([r0[0], r0[r0.length - 1]], [0, iw]);
  const ymax = Math.max(...iStar) * 1.12 || 0.1;
  const y = linScale([0, ymax], [ih, 0]);
  axes(g, x, y, iw, ih, {
    xTitle: "Basic reproduction number  R₀  (full system)", yTitle: "Endemic prevalence  I*",
    xTicks: niceTicks(r0[0], r0[r0.length - 1], 6), yTicks: niceTicks(0, ymax, 5),
    yFmt: (v) => `${Math.round(v * 100)}%`, xFmt: (v) => fx(v, 1),
  });
  // shade sub/supercritical
  const xt = x(threshold);
  g.insertBefore(s("rect", { x: 0, y: 0, width: xt, height: ih, fill: "var(--ink-3)", opacity: 0.06 }), g.firstChild);
  g.appendChild(s("line", { x1: xt, x2: xt, y1: 0, y2: ih, stroke: "var(--c-i)", "stroke-width": 1.4, "stroke-dasharray": "4 3" }));
  g.appendChild(s("text", { x: xt + 5, y: 13, class: "serieslabel", fill: "var(--c-i)", text: "R₀ = 1  (transcritical)" }));
  g.appendChild(s("text", { x: 6, y: ih - 8, class: "axis", fill: "var(--ink-3)", "font-size": 10.5, text: "disease-free stable" }));
  // disease-free branch (I*=0) up to threshold
  g.appendChild(s("path", { d: `M0,${ih}L${xt},${ih}`, fill: "none", stroke: "var(--ink-3)", "stroke-width": 2, "stroke-dasharray": "2 3" }));
  // endemic branch
  const xs2 = [], ys2 = [];
  for (let i = 0; i < r0.length; i++) if (r0[i] >= threshold) { xs2.push(r0[i]); ys2.push(iStar[i]); }
  g.appendChild(s("path", { d: linePath(xs2, ys2, x, y), fill: "none", stroke: "var(--accent)", "stroke-width": 2.4, "stroke-linecap": "round" }));
  g.appendChild(s("text", { x: iw - 4, y: y(ys2[ys2.length - 1]) - 8, "text-anchor": "end", class: "serieslabel", fill: "var(--accent)", text: "stable endemic branch" }));
  return svg;
}

// ---- phase portrait ---------------------------------------------------------
export function phasePortrait({ traj, eq, xKey = "A", yKey = "I", field = null, w = 480, h = 420 }) {
  const { svg, g, iw, ih } = frame({ w, h, m: { t: 20, r: 24, b: 46, l: 54 } });
  const xmax = Math.max(...traj[xKey]) * 1.15 || 0.3;
  const ymax = Math.max(...traj[yKey]) * 1.15 || 0.2;
  const x = linScale([0, xmax], [0, iw]);
  const y = linScale([0, ymax], [ih, 0]);
  axes(g, x, y, iw, ih, {
    xTitle: `${xKey}  (At-Risk fraction)`, yTitle: `${yKey}  (Addicted fraction)`,
    xTicks: niceTicks(0, xmax, 5), yTicks: niceTicks(0, ymax, 5),
    yFmt: (v) => `${Math.round(v * 100)}%`, xFmt: (v) => `${Math.round(v * 100)}%`,
  });
  // vector field
  if (field) {
    for (const v of field) {
      const x0 = x(v.x), y0 = y(v.y);
      const len = 11, ang = Math.atan2(-(v.dy) * (ih / ymax), v.dx * (iw / xmax));
      const x1 = x0 + Math.cos(ang) * len, y1 = y0 + Math.sin(ang) * len;
      g.appendChild(s("line", { x1: x0, y1: y0, x2: x1, y2: y1, stroke: "var(--ink-3)", "stroke-width": 0.8, opacity: 0.45, "marker-end": "url(#vfarrow)" }));
    }
    const defs = s("defs");
    const mk = s("marker", { id: "vfarrow", viewBox: "0 0 6 6", refX: 5, refY: 3, markerWidth: 5, markerHeight: 5, orient: "auto" });
    mk.appendChild(s("path", { d: "M0,0L6,3L0,6Z", fill: "var(--ink-3)", opacity: 0.5 }));
    defs.appendChild(mk); svg.insertBefore(defs, svg.firstChild);
  }
  // trajectory
  g.appendChild(s("path", { d: linePath(traj[xKey], traj[yKey], x, y), fill: "none", stroke: "var(--accent)", "stroke-width": 2, "stroke-linecap": "round" }));
  // start & equilibrium
  g.appendChild(s("circle", { cx: x(traj[xKey][0]), cy: y(traj[yKey][0]), r: 4, fill: "var(--panel)", stroke: "var(--ink)", "stroke-width": 1.5 }));
  g.appendChild(s("text", { x: x(traj[xKey][0]) + 7, y: y(traj[yKey][0]) + 3, class: "serieslabel", fill: "var(--ink-2)", text: "t = 0" }));
  if (eq) {
    g.appendChild(s("circle", { cx: x(eq[xKey]), cy: y(eq[yKey]), r: 5, fill: "var(--c-i)", stroke: "var(--panel)", "stroke-width": 1.5 }));
    g.appendChild(s("text", { x: x(eq[xKey]) + 7, y: y(eq[yKey]) + 3, class: "serieslabel", fill: "var(--c-i)", text: "endemic equilibrium" }));
  }
  return svg;
}

// ---- tornado (one-way sensitivity) -----------------------------------------
export function tornado({ rows, w = 720, h = null, valueFmt = (v) => fx(v, 3), baseline, metricLabel = "Peak prevalence" }) {
  const rowH = 30, top = 30, bot = 40;
  const H = h || top + rows.length * rowH + bot;
  const { svg, g, iw, ih } = frame({ w, h: H, m: { t: top, r: 30, b: bot, l: 128 } });
  let lo = Infinity, hi = -Infinity;
  for (const r of rows) { lo = Math.min(lo, r.low, r.high); hi = Math.max(hi, r.low, r.high); }
  lo = Math.min(lo, baseline); hi = Math.max(hi, baseline);
  const pad = (hi - lo) * 0.08 || 0.01;
  const x = linScale([lo - pad, hi + pad], [0, iw]);
  for (const t of niceTicks(lo - pad, hi + pad, 5)) {
    g.appendChild(s("line", { class: "grid", x1: x(t), x2: x(t), y1: 0, y2: rows.length * rowH, stroke: "var(--grid)" }));
    g.appendChild(s("text", { class: "axis", x: x(t), y: rows.length * rowH + 16, "text-anchor": "middle", fill: "var(--ink-3)", "font-size": 10, text: valueFmt(t) }));
  }
  // baseline line
  g.appendChild(s("line", { x1: x(baseline), x2: x(baseline), y1: -6, y2: rows.length * rowH, stroke: "var(--ink)", "stroke-width": 1.2 }));
  g.appendChild(s("text", { class: "serieslabel", x: x(baseline), y: -12, "text-anchor": "middle", fill: "var(--ink-2)", text: "baseline" }));
  rows.forEach((r, i) => {
    const yy = i * rowH + rowH / 2;
    const x0 = x(Math.min(r.low, r.high)), x1 = x(Math.max(r.low, r.high));
    g.appendChild(s("rect", { x: x0, y: yy - 8, width: Math.max(1, x1 - x0), height: 16, rx: 3, fill: "var(--accent)", opacity: 0.75 }));
    g.appendChild(s("text", { x: -10, y: yy + 3.5, "text-anchor": "end", class: "serieslabel", fill: "var(--ink-2)", text: r.label }));
  });
  g.appendChild(s("text", { class: "axis-title", x: iw / 2, y: rows.length * rowH + 34, "text-anchor": "middle", fill: "var(--ink-2)", text: `${metricLabel} across ±parameter range` }));
  return svg;
}

// ---- forest plot ------------------------------------------------------------
export function forest({ rows, w = 720, h = null }) {
  const rowH = 34, top = 24, bot = 44;
  const H = h || top + rows.length * rowH + bot;
  const { svg, g, iw, ih } = frame({ w, h: H, m: { t: top, r: 92, b: bot, l: 168 } });
  let lo = 0.5, hi = 3.5;
  for (const r of rows) { lo = Math.min(lo, r.lo, r.knownOR); hi = Math.max(hi, r.hi, r.knownOR); }
  const x = logScale([Math.max(0.4, lo * 0.9), hi * 1.1], [0, iw]);
  const ticks = [0.5, 0.7, 1, 1.5, 2, 3].filter((t) => t >= x.domain[0] && t <= x.domain[1]);
  for (const t of ticks) {
    g.appendChild(s("line", { class: "grid", x1: x(t), x2: x(t), y1: 0, y2: rows.length * rowH, stroke: t === 1 ? "var(--ink-3)" : "var(--grid)", "stroke-dasharray": t === 1 ? "4 3" : null }));
    g.appendChild(s("text", { class: "axis", x: x(t), y: rows.length * rowH + 16, "text-anchor": "middle", fill: "var(--ink-3)", "font-size": 10, text: fx(t, t < 1 ? 1 : t % 1 ? 1 : 0) }));
  }
  g.appendChild(s("text", { class: "axis", x: x(1), y: -8, "text-anchor": "middle", fill: "var(--ink-3)", "font-size": 10, text: "no effect" }));
  rows.forEach((r, i) => {
    const yy = i * rowH + rowH / 2;
    // CI whisker
    g.appendChild(s("line", { x1: x(r.lo), x2: x(r.hi), y1: yy, y2: yy, stroke: "var(--accent)", "stroke-width": 1.6 }));
    for (const cap of [r.lo, r.hi]) g.appendChild(s("line", { x1: x(cap), x2: x(cap), y1: yy - 4, y2: yy + 4, stroke: "var(--accent)", "stroke-width": 1.4 }));
    // estimated point
    g.appendChild(s("circle", { cx: x(r.estOR), cy: yy, r: 4.5, fill: "var(--accent)" }));
    // known OR marker (open diamond)
    const kx = x(r.knownOR);
    g.appendChild(s("path", { d: `M${kx},${yy - 5}L${kx + 5},${yy}L${kx},${yy + 5}L${kx - 5},${yy}Z`, fill: "none", stroke: "var(--c-i)", "stroke-width": 1.6 }));
    g.appendChild(s("text", { x: -12, y: yy + 3.5, "text-anchor": "end", class: "serieslabel", fill: "var(--ink-2)", text: r.label }));
    g.appendChild(s("text", { x: iw + 8, y: yy + 3.5, class: "axis", "font-size": 11, fill: "var(--ink-3)", text: `${fx(r.estOR, 2)} (${fx(r.lo, 2)}–${fx(r.hi, 2)})` }));
  });
  g.appendChild(s("text", { class: "axis-title", x: iw / 2, y: rows.length * rowH + 34, "text-anchor": "middle", fill: "var(--ink-2)", text: "Odds ratio for S→A onset (log scale)" }));
  return svg;
}

// ---- heatmap (transition matrix) -------------------------------------------
export function heatmap({ matrix, labels, title, w = 340, h = 340, accent = "--accent" }) {
  const n = matrix.length;
  const m = { t: 34, r: 12, b: 16, l: 34 };
  const { svg, g, iw, ih } = frame({ w, h, m });
  const cell = Math.min(iw, ih) / n;
  const size = cell * n;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const v = matrix[i][j];
      g.appendChild(s("rect", { x: j * cell, y: i * cell, width: cell - 1.5, height: cell - 1.5, rx: 3,
        fill: `var(${accent})`, "fill-opacity": 0.08 + 0.9 * Math.min(1, v) }));
      g.appendChild(s("text", { x: j * cell + cell / 2, y: i * cell + cell / 2 + 3.5, "text-anchor": "middle",
        class: "axis", "font-size": 11, fill: v > 0.5 ? "var(--panel)" : "var(--ink-2)", text: fx(v, 2) }));
    }
    g.appendChild(s("text", { x: -8, y: i * cell + cell / 2 + 3.5, "text-anchor": "end", class: "serieslabel", fill: "var(--ink-3)", text: labels[i] }));
    g.appendChild(s("text", { x: i * cell + cell / 2, y: -8, "text-anchor": "middle", class: "serieslabel", fill: "var(--ink-3)", text: labels[i] }));
  }
  if (title) svg.querySelector("g").appendChild(s("text", { x: size / 2, y: -22, "text-anchor": "middle", class: "axis-title", fill: "var(--ink-2)", text: title }));
  return svg;
}

// ---- sankey (2-column start→end flows) -------------------------------------
export function sankey({ flows, labels, colors, w = 560, h = 360 }) {
  const n = flows.length;
  const m = { t: 20, r: 90, b: 20, l: 90 };
  const { svg, g, iw, ih } = frame({ w, h, m });
  const total = flows.flat().reduce((a, b) => a + b, 0);
  const gap = 10;
  const avail = ih - gap * (n - 1);
  const outSum = flows.map((r) => r.reduce((a, b) => a + b, 0));
  const inSum = labels.map((_, j) => flows.reduce((a, r) => a + r[j], 0));
  // node y positions
  const leftY = []; let cy = 0;
  for (let i = 0; i < n; i++) { const hgt = (outSum[i] / total) * avail; leftY.push([cy, cy + hgt]); cy += hgt + gap; }
  const rightY = []; cy = 0;
  for (let j = 0; j < n; j++) { const hgt = (inSum[j] / total) * avail; rightY.push([cy, cy + hgt]); cy += hgt + gap; }
  const nodeW = 12;
  const leftCursor = leftY.map((r) => r[0]);
  const rightCursor = rightY.map((r) => r[0]);
  // ribbons
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const v = flows[i][j];
      if (v <= 0) continue;
      const th = (v / total) * avail;
      const y0 = leftCursor[i], y1 = leftCursor[i] + th; leftCursor[i] += th;
      const y2 = rightCursor[j], y3 = rightCursor[j] + th; rightCursor[j] += th;
      const x0 = nodeW, x1 = iw - nodeW;
      const cxm = (x0 + x1) / 2;
      const d = `M${x0},${y0}C${cxm},${y0} ${cxm},${y2} ${x1},${y2}L${x1},${y3}C${cxm},${y3} ${cxm},${y1} ${x0},${y1}Z`;
      g.appendChild(s("path", { d, fill: colors[i], "fill-opacity": i === j ? 0.16 : 0.34, stroke: "none" }));
    }
  }
  // nodes
  for (let i = 0; i < n; i++) {
    g.appendChild(s("rect", { x: 0, y: leftY[i][0], width: nodeW, height: Math.max(2, leftY[i][1] - leftY[i][0]), rx: 2, fill: colors[i] }));
    g.appendChild(s("text", { x: -8, y: (leftY[i][0] + leftY[i][1]) / 2 + 3.5, "text-anchor": "end", class: "serieslabel", fill: "var(--ink-2)", text: labels[i] }));
    g.appendChild(s("rect", { x: iw - nodeW, y: rightY[i][0], width: nodeW, height: Math.max(2, rightY[i][1] - rightY[i][0]), rx: 2, fill: colors[i] }));
    g.appendChild(s("text", { x: iw + 8, y: (rightY[i][0] + rightY[i][1]) / 2 + 3.5, class: "serieslabel", fill: "var(--ink-2)", text: labels[i] }));
  }
  g.appendChild(s("text", { x: nodeW, y: ih + 14, class: "axis", "font-size": 10, fill: "var(--ink-3)", text: "baseline state" }));
  g.appendChild(s("text", { x: iw - nodeW, y: ih + 14, "text-anchor": "end", class: "axis", "font-size": 10, fill: "var(--ink-3)", text: "final state" }));
  return svg;
}
