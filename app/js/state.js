// state.js — Central scenario state: params, presets, URL sharing, localStorage.
// -----------------------------------------------------------------------------
// A single observable store. Views subscribe and re-render on change. The full
// scenario is URL-encoded (query string) so any state is shareable, and user
// scenarios persist to localStorage.

import { DEFAULT_PARAMS, PARAM_KEYS } from "./model/params.js";

const LS_KEY = "dpl.scenarios.v1";

/** Built-in presets. Each carries params + optional network/intervention hints. */
export const PRESETS = {
  baseline: {
    label: "Calibrated baseline",
    params: { ...DEFAULT_PARAMS },
    interventions: [],
    network: { type: "BA", n: 800, meanDegree: 10, runs: 60, immunizeFrac: 0, immunizeStrategy: "hub" },
    note: "Chapter-6 calibration reproducing the thesis targets.",
  },
  scalefree: {
    label: "Scale-free outbreak",
    params: { ...DEFAULT_PARAMS },
    interventions: [],
    network: { type: "BA", n: 1000, meanDegree: 10, runs: 80, immunizeFrac: 0, immunizeStrategy: "hub" },
    note: "Influencer-hub topology: higher, earlier peak via degree heterogeneity.",
    view: "network",
  },
  comprehensive: {
    label: "Comprehensive intervention",
    params: { ...DEFAULT_PARAMS },
    interventions: ["comprehensive"],
    network: { type: "BA", n: 800, meanDegree: 10, runs: 60, immunizeFrac: 0, immunizeStrategy: "hub" },
    note: "All levers combined: ≈ −63% peak (conditional on assumed efficacies).",
    view: "interventions",
  },
};

// Short keys for compact URLs.
const URL_KEYS = { beta: "b", k: "k", eta: "e", mu: "m", sigma: "s", gamma1: "g", gamma2: "G", delta: "d" };

class Store {
  constructor() {
    this.state = {
      params: { ...DEFAULT_PARAMS },
      horizonYears: 20,
      interventions: [], // keys into INTERVENTIONS or ['comprehensive']
      customMods: {}, // user-built multipliers { param: multiplier }
      compare: null, // { label, params } overlay
      network: { type: "BA", n: 800, meanDegree: 10, runs: 60, seed: 1, immunizeFrac: 0, immunizeStrategy: "hub" },
      view: "overview",
    };
    this.subs = new Set();
    this._silent = false;
  }
  get() { return this.state; }
  subscribe(fn) { this.subs.add(fn); return () => this.subs.delete(fn); }
  emit() { if (!this._silent) for (const fn of this.subs) fn(this.state); }

  /** Patch top-level fields (shallow merge on nested params/network). */
  set(patch, { silent = false, pushUrl = true } = {}) {
    if (patch.params) patch.params = { ...this.state.params, ...patch.params };
    if (patch.network) patch.network = { ...this.state.network, ...patch.network };
    this.state = { ...this.state, ...patch };
    if (pushUrl) this.syncUrl();
    if (!silent) this.emit();
  }
  setParam(key, val, opts) { this.set({ params: { [key]: val } }, opts); }
  resetParams(opts) { this.set({ params: { ...DEFAULT_PARAMS }, interventions: [], customMods: {} }, opts); }

  applyPreset(name, opts = {}) {
    const p = PRESETS[name];
    if (!p) return;
    this.set({
      params: { ...p.params },
      interventions: [...(p.interventions || [])],
      customMods: {},
      network: { ...this.state.network, ...(p.network || {}) },
      ...(p.view ? { view: p.view } : {}),
    }, opts);
  }

  // ---- URL sharing ----------------------------------------------------------
  syncUrl() {
    const q = new URLSearchParams();
    for (const key of PARAM_KEYS) {
      const v = this.state.params[key];
      if (v !== DEFAULT_PARAMS[key]) q.set(URL_KEYS[key], round(v));
    }
    if (this.state.horizonYears !== 20) q.set("H", this.state.horizonYears);
    if (this.state.interventions.length) q.set("iv", this.state.interventions.join("."));
    const custom = Object.entries(this.state.customMods).filter(([, v]) => v !== 1);
    if (custom.length) q.set("cm", custom.map(([k, v]) => `${URL_KEYS[k]}~${round(v)}`).join("."));
    const nw = this.state.network;
    if (nw.type !== "BA") q.set("nt", nw.type);
    if (nw.immunizeFrac) q.set("im", `${round(nw.immunizeFrac)}${nw.immunizeStrategy[0]}`);
    if (this.state.view !== "overview") q.set("v", this.state.view);
    const qs = q.toString();
    const url = `${location.pathname}${qs ? "?" + qs : ""}`;
    history.replaceState(null, "", url);
  }

  loadFromUrl() {
    const q = new URLSearchParams(location.search);
    const params = { ...DEFAULT_PARAMS };
    for (const [key, sk] of Object.entries(URL_KEYS)) {
      if (q.has(sk)) { const v = parseFloat(q.get(sk)); if (isFinite(v)) params[key] = v; }
    }
    const patch = { params };
    if (q.has("H")) patch.horizonYears = clampInt(parseInt(q.get("H")), 5, 40);
    if (q.has("iv")) patch.interventions = q.get("iv").split(".").filter(Boolean);
    if (q.has("cm")) {
      const cm = {};
      for (const tok of q.get("cm").split(".")) {
        const [sk, v] = tok.split("~");
        const key = Object.keys(URL_KEYS).find((k) => URL_KEYS[k] === sk);
        if (key) cm[key] = parseFloat(v);
      }
      patch.customMods = cm;
    }
    const nw = {};
    if (q.has("nt")) nw.type = q.get("nt");
    if (q.has("im")) {
      const raw = q.get("im"); const frac = parseFloat(raw);
      nw.immunizeFrac = isFinite(frac) ? frac : 0;
      nw.immunizeStrategy = /r$/.test(raw) ? "random" : "hub";
    }
    if (Object.keys(nw).length) patch.network = nw;
    if (q.has("v")) patch.view = q.get("v");
    this.set(patch, { pushUrl: false, silent: true });
  }

  shareUrl() { this.syncUrl(); return location.href; }

  // ---- localStorage scenarios ----------------------------------------------
  savedScenarios() {
    try { return JSON.parse(localStorage.getItem(LS_KEY) || "[]"); } catch { return []; }
  }
  saveScenario(name) {
    const list = this.savedScenarios();
    const entry = {
      name, ts: Date.now(),
      params: { ...this.state.params },
      interventions: [...this.state.interventions],
      network: { ...this.state.network },
    };
    const idx = list.findIndex((e) => e.name === name);
    if (idx >= 0) list[idx] = entry; else list.push(entry);
    localStorage.setItem(LS_KEY, JSON.stringify(list));
    return list;
  }
  deleteScenario(name) {
    const list = this.savedScenarios().filter((e) => e.name !== name);
    localStorage.setItem(LS_KEY, JSON.stringify(list));
    return list;
  }
  loadScenario(entry, opts) {
    this.set({ params: { ...entry.params }, interventions: [...(entry.interventions || [])],
      network: { ...this.state.network, ...(entry.network || {}) }, customMods: {} }, opts);
  }
}

function round(v) { return Math.round(v * 1e5) / 1e5; }
function clampInt(v, lo, hi) { return Math.max(lo, Math.min(hi, v || lo)); }

export const store = new Store();
