// app.js — Application shell: routing, navigation, scenario bar, theme.
// -----------------------------------------------------------------------------

import { h, $, clear, toast } from "./view/dom.js";
import { store, PRESETS } from "./state.js";
import { OverviewView } from "./view/views/overview.js";
import { ExplorerView } from "./view/views/explorer.js";
import { PhaseView } from "./view/views/phase.js";
import { NetworkLabView } from "./view/views/network-lab.js";
import { CohortView } from "./view/views/cohort.js";
import { InterventionsView } from "./view/views/interventions.js";
import { MethodsView } from "./view/views/methods.js";
import { selfTestPanel } from "./view/views/selftest-panel.js";

const VIEWS = [
  { id: "overview", title: "Overview", desc: "Baseline S–A–I–R run", make: OverviewView },
  { id: "explorer", title: "Parameter Explorer", desc: "Sliders, R₀, sensitivity", make: ExplorerView },
  { id: "phase", title: "Phase & Equilibrium", desc: "Bifurcation, phase, lag", make: PhaseView },
  { id: "network", title: "Network Lab", desc: "Topologies, immunisation", make: NetworkLabView },
  { id: "cohort", title: "Cohort & Risk", desc: "Forest, Sankey, HMM", make: CohortView },
  { id: "interventions", title: "Intervention Simulator", desc: "Compose levers, cost", make: InterventionsView },
  { id: "methods", title: "Methods & Provenance", desc: "Equations, sources, caveats", make: MethodsView },
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
  if (store.get().view === id && current) return;
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

// ---- scenario bar ----
function buildScenarioBar() {
  const bar = h("div", { class: "scenariobar" });
  const rebuild = () => { clear(bar); fill(); };

  function fill() {
    bar.appendChild(h("span", { class: "grp", style: { color: "var(--ink-3)", fontWeight: 600 } }, ["Presets"]));
    for (const [key, p] of Object.entries(PRESETS)) {
      bar.appendChild(h("button", { class: "sm", title: p.note, onclick: () => { store.applyPreset(key); toast(`Applied: ${p.label}`); } }, [p.label]));
    }
    bar.appendChild(h("span", { style: { width: "1px", height: "20px", background: "var(--rule-2)" } }));
    // compare
    const cmp = store.get().compare;
    if (cmp) {
      bar.appendChild(h("span", { class: "pill", text: `compare: ${cmp.label}` }));
      bar.appendChild(h("button", { class: "sm ghost", onclick: () => { store.set({ compare: null }); toast("Compare cleared"); } }, ["clear compare"]));
    } else {
      bar.appendChild(h("button", { class: "sm", title: "Pin the current scenario to overlay against edits", onclick: () => {
        const s = store.get();
        store.set({ compare: { label: describe(s), params: { ...s.params }, interventions: [...s.interventions] } });
        toast("Pinned current scenario for comparison");
      } }, ["⎘ Pin for compare"]));
    }
    bar.appendChild(h("span", { class: "grp spacer", style: { flex: 1 } }));
    // save/load
    bar.appendChild(h("button", { class: "sm ghost", onclick: openSaveMenu }, ["✔ Save scenario"]));
    bar.appendChild(h("button", { class: "sm ghost", onclick: openLoadMenu }, ["↺ Saved…"]));
  }
  fill();
  store.subscribe(rebuild);
  return bar;
}

function describe(s) {
  if (s.interventions.includes("comprehensive")) return "comprehensive";
  if (s.interventions.length) return s.interventions.join("+");
  const changed = Object.keys(s.params).filter((k) => s.params[k] !== PRESETS.baseline.params[k]);
  return changed.length ? "custom" : "baseline";
}

function openSaveMenu() {
  const input = h("input", { type: "text", placeholder: "Scenario name", style: { padding: "8px 10px", border: "1px solid var(--rule-2)", borderRadius: "6px", width: "100%", background: "var(--panel)", color: "var(--ink)" } });
  const body = h("div", {}, [
    h("p", { class: "prose", text: "Save the current parameters and interventions to this browser (localStorage)." }),
    input,
    h("div", { class: "row", style: { marginTop: "12px" } }, [
      h("button", { class: "primary sm", onclick: () => { const n = input.value.trim() || `scenario ${new Date().toISOString().slice(0, 16).replace("T", " ")}`; store.saveScenario(n); toast(`Saved “${n}”`); close(); } }, ["Save"]),
    ]),
  ]);
  const close = openModal("Save scenario", body);
  setTimeout(() => input.focus(), 30);
}

function openLoadMenu() {
  const list = store.savedScenarios();
  const body = h("div", {});
  if (!list.length) body.appendChild(h("p", { class: "prose", text: "No saved scenarios yet. Use “Save scenario” to store the current setup." }));
  else {
    const table = h("table", { class: "data" }, [
      h("thead", {}, h("tr", {}, [h("th", { style: { textAlign: "left" }, text: "Name" }), h("th", { text: "Saved" }), h("th", {})])),
      h("tbody", {}, list.slice().reverse().map((e) => h("tr", {}, [
        h("td", { style: { textAlign: "left" }, text: e.name }),
        h("td", { class: "hint", text: new Date(e.ts).toLocaleDateString() }),
        h("td", {}, h("div", { class: "row", style: { justifyContent: "flex-end" } }, [
          h("button", { class: "sm", onclick: () => { store.loadScenario(e); toast(`Loaded “${e.name}”`); close(); } }, ["load"]),
          h("button", { class: "sm ghost", onclick: () => { store.deleteScenario(e.name); close(); openLoadMenu(); } }, ["delete"]),
        ])),
      ]))),
    ]);
    body.appendChild(table);
  }
  const close = openModal("Saved scenarios", body);
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
        h("span", { class: "mark", text: "The Digital Pandemic" }),
        h("span", { class: "sub", text: "Interactive Simulation Lab · Mutibura 2026" }),
      ]),
      h("span", { class: "spacer" }),
      h("div", { class: "tools" }, [
        h("button", { class: "sm ghost", onclick: () => openModal("Correctness self-test", selfTestPanel()), title: "Run acceptance self-tests" }, ["✓ Self-test"]),
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
          h("p", { class: "hint", html: "All outputs are <strong>conditional model results</strong>, not empirical forecasts." }),
        ]),
      ]),
      h("main", { class: "content", id: "content", tabindex: "-1" }, [
        buildScenarioBar(),
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
