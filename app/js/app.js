// app.js: application shell (routing, navigation, theme, self-test, sharing).
// -----------------------------------------------------------------------------

import { h, $, clear, toast } from "./view/dom.js";
import { store } from "./state.js";
import { ContextView } from "./view/views/context.js";
import { ExplanationsView } from "./view/views/explanations.js";
import { ExplorerView } from "./view/views/explorer.js";
import { ThresholdsView } from "./view/views/thresholds.js";
import { NetworksView } from "./view/views/networks.js";
import { InterventionsView } from "./view/views/interventions.js";
import { PanelView } from "./view/views/panel.js";
import { AgeSexView } from "./view/views/agesex.js";
import { MethodsView } from "./view/views/methods.js";
import { selfTestPanel } from "./view/views/selftest-panel.js";

const VIEWS = [
  { id: "context", title: "Research context", desc: "What the research asks", make: ContextView },
  { id: "explanations", title: "Two explanations", desc: "The rise, projections, next wave", make: ExplanationsView },
  { id: "explorer", title: "Model explorer", desc: "Rates, R₀, equilibrium, calibration", make: ExplorerView },
  { id: "thresholds", title: "Thresholds and floor", desc: "Why the threshold disappears", make: ThresholdsView },
  { id: "networks", title: "Networks", desc: "Structure and targeting", make: NetworksView },
  { id: "interventions", title: "Interventions", desc: "Mechanisms and leverage", make: InterventionsView },
  { id: "panel", title: "Panel design", desc: "Data that would identify the rates", make: PanelView },
  { id: "agesex", title: "Age and sex", desc: "Structured model, 2018 test", make: AgeSexView },
  { id: "methods", title: "Methods and sources", desc: "Equations, parameters, references", make: MethodsView },
];

let current = null;
const contentEl = () => $("#view-holder");

function mountView(id) {
  const v = VIEWS.find((x) => x.id === id) || VIEWS[0];
  if (current) current.destroy?.();
  const holder = contentEl();
  clear(holder);
  current = v.make();
  holder.appendChild(current.el);
  holder.scrollTop = 0; window.scrollTo({ top: 0, behavior: "instant" in window ? "instant" : "auto" });
  // nav highlight
  for (const a of document.querySelectorAll("nav.sidebar a.navlink"))
    a.setAttribute("aria-current", a.dataset.view === v.id ? "page" : "false");
  document.title = `${v.title} · Digital Pandemic Lab`;
}

function navigate(id) {
  $("nav.sidebar")?.classList.remove("open");
  if (store.get().view === id && current) return;
  if (location.hash !== `#${id}`) history.pushState(null, "", `#${id}`);
  store.set({ view: id }, { silent: true });
  $("nav.sidebar")?.classList.remove("open");
  mountView(id);
}

// ---- theme ----
function initTheme() {
  const saved = localStorage.getItem("dpl.theme");
  if (saved) document.documentElement.setAttribute("data-theme", saved);
  updateThemeBtn();
}
function cycleTheme() {
  const cur = document.documentElement.getAttribute("data-theme");
  const next = cur === "dark" ? "light" : cur === "light" ? "" : "dark";
  if (next) { document.documentElement.setAttribute("data-theme", next); localStorage.setItem("dpl.theme", next); }
  else { document.documentElement.removeAttribute("data-theme"); localStorage.removeItem("dpl.theme"); }
  updateThemeBtn();
}
function updateThemeBtn() {
  const btn = $("#themeBtn"); if (!btn) return;
  const cur = document.documentElement.getAttribute("data-theme") || "auto";
  btn.textContent = cur === "dark" ? "◑ Dark" : cur === "light" ? "◐ Light" : "◒ Auto";
}

// ---- modal ----
function openModal(title, node) {
  const overlay = h("div", { class: "modal-overlay", role: "dialog", "aria-modal": "true", "aria-label": title,
    style: { position: "fixed", inset: 0, background: "rgba(0,0,0,.45)", zIndex: 90, display: "grid", placeItems: "center", padding: "20px" },
    onclick: (e) => { if (e.target === overlay) close(); } });
  const box = h("div", { class: "panel", style: { maxWidth: "820px", width: "100%", maxHeight: "86vh", overflow: "auto", margin: 0 } }, [
    h("div", { class: "row between", style: { marginBottom: "12px" } }, [h("h2", { text: title }), h("button", { class: "ghost sm", onclick: () => close(), "aria-label": "Close", text: "✕ Close" })]),
    node,
  ]);
  overlay.appendChild(box);
  document.body.appendChild(overlay);
  const onKey = (e) => { if (e.key === "Escape") close(); };
  document.addEventListener("keydown", onKey);
  function close() { overlay.remove(); document.removeEventListener("keydown", onKey); }
  return close;
}

async function share() {
  const url = store.shareUrl();
  try { await navigator.clipboard.writeText(url); toast("Shareable link copied to clipboard"); }
  catch { openModal("Shareable link", h("div", {}, [h("p", { class: "prose", text: "Copy this URL to share the current scenario:" }), h("input", { type: "text", value: url, readonly: true, style: { width: "100%", padding: "8px", border: "1px solid var(--rule-2)", borderRadius: "6px", background: "var(--panel)", color: "var(--ink)" }, onclick: (e) => e.target.select() })])); }
}

// ---- build shell ----
function build() {
  const app = h("div", { class: "app" }, [
    h("a", { class: "skip-link", href: "#content", text: "Skip to content" }),
    h("header", { class: "topbar" }, [
      h("button", { class: "navtoggle ghost sm", "aria-label": "Toggle navigation", onclick: () => $("nav.sidebar")?.classList.toggle("open") }, ["☰"]),
      h("div", { class: "brand" }, [
        h("span", { class: "mark", text: "Digital Pandemic Lab" }),
        h("span", { class: "sub", text: "Problematic social media use among young people · Interactive model lab" }),
      ]),
      h("span", { class: "spacer" }),
      h("div", { class: "tools" }, [
        h("button", { class: "sm ghost", onclick: () => openModal("Self-test against the thesis model", selfTestPanel()), title: "Check the lab model against the thesis model" }, ["✓ Self-test"]),
        h("button", { class: "sm ghost", onclick: share, title: "Copy shareable URL" }, ["⇪ Share"]),
        h("button", { id: "themeBtn", class: "sm ghost", onclick: cycleTheme }, ["◒ Auto"]),
      ]),
    ]),
    h("div", { class: "layout" }, [
      h("nav", { class: "sidebar", "aria-label": "Screens" }, [
        h("ol", {}, VIEWS.map((v, i) => h("li", {}, h("a", {
          class: "navlink", href: `#${v.id}`, dataset: { view: v.id },
          onclick: (e) => { e.preventDefault(); navigate(v.id); },
        }, [
          h("span", { class: "n", text: String(i + 1) }),
          h("span", {}, [h("span", { class: "t", text: v.title }), h("span", { class: "d", text: v.desc })]),
        ])))),
        h("div", { style: { padding: "12px 11px", marginTop: "8px", borderTop: "1px solid var(--rule)" } }, [
          h("p", { class: "hint", html: "All outputs are <strong>model results conditional on stated assumptions</strong>, not forecasts." }),
        ]),
      ]),
      h("main", { class: "content", id: "content", tabindex: "-1" }, [
        h("div", { id: "view-holder" }),
      ]),
    ]),
  ]);
  document.body.appendChild(app);
}

// ---- boot ----
function boot() {
  store.loadFromUrl();
  build();
  initTheme();
  mountView(store.get().view);
  window.addEventListener("hashchange", () => {
    const id = location.hash.replace("#", "");
    if (VIEWS.find((v) => v.id === id)) navigate(id);
  });
}

if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", boot);
else boot();
