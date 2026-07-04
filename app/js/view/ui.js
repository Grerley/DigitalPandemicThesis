// ui.js — Shared, composable view components.
// -----------------------------------------------------------------------------

import { h, s, clear, toast } from "./dom.js";
import { exportPNG, exportSVG, exportCSV } from "./export.js";

/** A titled panel. */
export function panel({ title, sub, badge, children = [], actions } = {}) {
  const head = title
    ? h("div", { class: "row between", style: { alignItems: "baseline", marginBottom: sub ? "2px" : "12px" } }, [
        h("div", {}, [
          h("h2", { text: title }),
          sub ? h("div", { class: "panel-sub", text: sub }) : null,
        ]),
        actions || (badge ? h("span", { class: "badge-cond", text: badge }) : null),
      ])
    : null;
  return h("section", { class: "panel" }, [head, ...(Array.isArray(children) ? children : [children])].filter(Boolean));
}

/** Metric strip. items: [{label, value, unit?, foot?, accent?}]. */
export function metricGrid(items) {
  return h("div", { class: "metrics" }, items.map((it) =>
    h("div", { class: "metric" + (it.accent ? " accent" : "") }, [
      h("div", { class: "label", text: it.label }),
      h("div", { class: "metric-value" }, [
        String(it.value),
        it.unit ? h("small", { text: " " + it.unit }) : null,
      ].filter(Boolean)),
      it.foot ? h("div", { class: "foot", text: it.foot }) : null,
    ])
  ));
}

/**
 * A chart card: holds an <svg> plus a toolbar with SVG/PNG (and optional CSV)
 * export. `render()` returns an svg element; re-render by calling update().
 */
export function chartCard({ title, sub, caption, filename = "chart", render, csv, controls }) {
  const holder = h("figure", { class: "chart" });
  const svgWrap = h("div");
  const build = () => {
    clear(svgWrap);
    const svg = render();
    svg.__filename = filename;
    svgWrap.appendChild(svg);
    return svg;
  };
  let svg = build();
  const btns = [
    h("button", { class: "sm ghost", title: "Download SVG", onclick: () => exportSVG(svgWrap.querySelector("svg"), `${filename}.svg`) }, ["SVG"]),
    h("button", { class: "sm ghost", title: "Download PNG", onclick: () => exportPNG(svgWrap.querySelector("svg"), `${filename}.png`) }, ["PNG"]),
  ];
  if (csv) btns.unshift(h("button", { class: "sm ghost", title: "Download CSV", onclick: () => exportCSV(csv(), `${filename}.csv`) }, ["CSV"]));
  const toolbar = h("div", { class: "chart-toolbar" }, [
    title ? h("strong", { style: { fontSize: "14px" }, text: title }) : null,
    sub ? h("span", { class: "hint", text: sub }) : null,
    h("span", { class: "spacer" }),
    ...(controls || []), ...btns,
  ].filter(Boolean));
  holder.appendChild(toolbar);
  holder.appendChild(svgWrap);
  if (caption) holder.appendChild(h("figcaption", { class: "chart-cap", text: caption }));
  holder.update = () => { svg = build(); };
  holder.svgWrap = svgWrap;
  return holder;
}

/** Legend row for compartments/series. items: [{name,color,dash?}]. */
export function legend(items) {
  return h("div", { class: "legend" }, items.map((it) =>
    h("span", { class: "item" }, [
      h("span", { class: "swatch", style: { borderTopColor: it.color, borderTopStyle: it.dash && it.dash !== "none" ? "dashed" : "solid" } }),
      it.name,
    ])
  ));
}

/**
 * Slider row bound to a parameter.
 * @param {object} o { sym,label,unit,desc,min,max,step,value,baseline,onInput }
 */
export function sliderRow(o) {
  const valEl = h("span", { class: "val" });
  const fmt = (v) => (o.step < 0.001 ? v.toFixed(4) : o.step < 0.01 ? v.toFixed(3) : v.toFixed(o.step < 0.1 ? 2 : 1));
  const range = h("input", {
    type: "range", min: o.min, max: o.max, step: o.step, value: o.value,
    "aria-label": `${o.label} (${o.sym})`,
    oninput: (e) => { const v = parseFloat(e.target.value); valEl.textContent = fmt(v); mark(v); o.onInput(v); },
  });
  const mark = (v) => valEl.classList.toggle("changed", Math.abs(v - o.baseline) > o.step / 2);
  valEl.textContent = fmt(o.value); mark(o.value);
  const row = h("div", { class: "slider-row" }, [
    h("div", { class: "slider-top" }, [
      h("span", { class: "name" }, [h("span", { class: "sym", text: o.sym }), o.label]),
      valEl,
    ]),
    range,
    h("div", { class: "slider-desc" }, [o.desc, o.unit && o.unit !== "—" ? ` · ${o.unit}` : ""]),
  ]);
  row.setValue = (v) => { range.value = v; valEl.textContent = fmt(v); mark(v); };
  return row;
}

/**
 * Play / pause / scrub timeline controller.
 * @param {object} o { max, value, onScrub, labelFor, fps }
 * returns { el, setValue, stop }
 */
export function timeline(o) {
  let playing = false, raf = null, last = 0;
  const label = h("span", { class: "tl-label", text: o.labelFor(o.value) });
  const range = h("input", { type: "range", min: 0, max: o.max, step: 1, value: o.value, "aria-label": "Timeline scrubber",
    oninput: (e) => { stop(); update(parseInt(e.target.value)); } });
  const playBtn = h("button", { class: "tl-btn", "aria-label": "Play", title: "Play / pause", onclick: () => (playing ? stop() : play()) });
  playBtn.innerHTML = playIcon();
  function update(v) { range.value = v; label.textContent = o.labelFor(v); o.onScrub(v); }
  function play() {
    if (parseInt(range.value) >= o.max) update(0);
    playing = true; playBtn.innerHTML = pauseIcon(); playBtn.setAttribute("aria-label", "Pause");
    last = performance.now(); loop();
  }
  function loop() {
    if (!playing) return;
    raf = requestAnimationFrame((now) => {
      const dt = now - last;
      const fps = o.fps || 24;
      if (dt > 1000 / fps) {
        last = now;
        let v = parseInt(range.value) + 1;
        if (v > o.max) { stop(); return; }
        update(v);
      }
      loop();
    });
  }
  function stop() { playing = false; cancelAnimationFrame(raf); playBtn.innerHTML = playIcon(); playBtn.setAttribute("aria-label", "Play"); }
  const el = h("div", { class: "timeline" }, [playBtn, range, label]);
  return { el, setValue: (v) => update(v), stop, max: o.max };
}

function playIcon() { return '<svg width="14" height="14" viewBox="0 0 14 14" aria-hidden="true"><path d="M3 2l9 5-9 5z" fill="currentColor"/></svg>'; }
function pauseIcon() { return '<svg width="14" height="14" viewBox="0 0 14 14" aria-hidden="true"><rect x="3" y="2" width="3" height="10" fill="currentColor"/><rect x="8" y="2" width="3" height="10" fill="currentColor"/></svg>'; }

/** Small conditional-results note. */
export function condNote(text) {
  return h("div", { class: "caveat", html: `<strong>Conditional model result.</strong> ${text}` });
}
