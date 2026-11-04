// state.js: the lab's small observable store (current screen and the explorer's parameters), with URL sharing.

import { centralParams } from "./model/params.js";
import { PARAM_KEYS } from "./model/sair.js";

class Store {
  constructor() {
    this.state = { view: "context", params: centralParams() };
    this.subs = new Set();
  }
  get() { return this.state; }
  subscribe(fn) { this.subs.add(fn); return () => this.subs.delete(fn); }
  set(patch, { silent = false } = {}) {
    const params = patch.params ? { ...this.state.params, ...patch.params } : this.state.params;
    this.state = { ...this.state, ...patch, params };
    if (!silent) for (const fn of this.subs) fn(this.state);
  }
  resetParams() { this.set({ params: centralParams() }); }
  loadFromUrl() {
    const view = location.hash.replace("#", "");
    if (view) this.state.view = view;
    const q = new URLSearchParams(location.search);
    for (const k of PARAM_KEYS) if (q.has(k) && Number.isFinite(+q.get(k))) this.state.params[k] = +q.get(k);
  }
  shareUrl() {
    const q = new URLSearchParams();
    for (const k of PARAM_KEYS) q.set(k, String(+this.state.params[k].toPrecision(6)));
    return `${location.origin}${location.pathname}?${q}#${this.state.view}`;
  }
}

export const store = new Store();
