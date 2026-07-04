// export.js — CSV / SVG / PNG export helpers.
// -----------------------------------------------------------------------------

import { toast } from "./dom.js";

/** Trigger a browser download of a Blob. */
export function download(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

/** Export an array of row-objects as CSV. */
export function exportCSV(rows, filename = "trajectory.csv") {
  if (!rows.length) return;
  const cols = Object.keys(rows[0]);
  const esc = (v) => {
    const str = v == null ? "" : String(v);
    return /[",\n]/.test(str) ? `"${str.replace(/"/g, '""')}"` : str;
  };
  const lines = [cols.join(",")];
  for (const r of rows) lines.push(cols.map((c) => esc(r[c])).join(","));
  download(new Blob([lines.join("\n")], { type: "text/csv;charset=utf-8" }), filename);
  toast(`Exported ${rows.length} rows → ${filename}`);
}

/** Serialize an <svg> element (inlining computed CSS colours) to an SVG string. */
export function svgToString(svg) {
  const clone = svg.cloneNode(true);
  inlineStyles(svg, clone);
  clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
  clone.setAttribute("xmlns:xlink", "http://www.w3.org/1999/xlink");
  // Paint a solid background matching the page.
  const bg = getComputedStyle(document.body).backgroundColor;
  const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
  rect.setAttribute("x", 0); rect.setAttribute("y", 0);
  rect.setAttribute("width", "100%"); rect.setAttribute("height", "100%");
  rect.setAttribute("fill", bg);
  clone.insertBefore(rect, clone.firstChild);
  return new XMLSerializer().serializeToString(clone);
}

// Copy resolved presentation styles from live nodes onto the clone so exports
// don't depend on the external stylesheet or CSS variables.
const STYLE_PROPS = ["fill", "stroke", "stroke-width", "stroke-dasharray", "stroke-linecap",
  "stroke-linejoin", "opacity", "fill-opacity", "stroke-opacity", "font-size", "font-family",
  "font-weight", "font-style", "text-anchor", "dominant-baseline"];
function inlineStyles(src, dst) {
  const s = getComputedStyle(src);
  let style = "";
  for (const p of STYLE_PROPS) {
    const v = s.getPropertyValue(p);
    if (v) style += `${p}:${v};`;
  }
  if (style) dst.setAttribute("style", style);
  const sc = src.children, dc = dst.children;
  for (let i = 0; i < sc.length; i++) inlineStyles(sc[i], dc[i]);
}

export function exportSVG(svg, filename = "chart.svg") {
  const str = svgToString(svg);
  download(new Blob([str], { type: "image/svg+xml;charset=utf-8" }), filename);
  toast(`Exported ${filename}`);
}

/** Rasterize an <svg> element to a PNG at `scale`× resolution. */
export function exportPNG(svg, filename = "chart.png", scale = 2) {
  const str = svgToString(svg);
  const vb = svg.viewBox.baseVal;
  const w = (vb && vb.width) || svg.clientWidth || 800;
  const hgt = (vb && vb.height) || svg.clientHeight || 500;
  const img = new Image();
  const blob = new Blob([str], { type: "image/svg+xml;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  img.onload = () => {
    const canvas = document.createElement("canvas");
    canvas.width = w * scale; canvas.height = hgt * scale;
    const ctx = canvas.getContext("2d");
    ctx.scale(scale, scale);
    ctx.drawImage(img, 0, 0, w, hgt);
    URL.revokeObjectURL(url);
    canvas.toBlob((b) => { download(b, filename); toast(`Exported ${filename}`); }, "image/png");
  };
  img.onerror = () => { URL.revokeObjectURL(url); toast("PNG export failed"); };
  img.src = url;
}
