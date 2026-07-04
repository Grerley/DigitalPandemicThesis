// scenario.js — Derive effective parameters & trajectories from store state.
// -----------------------------------------------------------------------------
// Bridges the store and the model layer: applies selected interventions and any
// user-built custom modifiers to the base parameters, then runs the model.

import { simulate, summarise, transitionMatrix } from "../model/sair.js";
import { r0Trio, r0Full } from "../model/r0.js";
import { baselineInit } from "../model/params.js";
import { INTERVENTIONS, COMPREHENSIVE, composeMods, applyMods } from "../model/interventions.js";

/** Combine intervention presets + custom modifiers into one multiplier map. */
export function activeMods(state) {
  let mods = {};
  const ivs = state.interventions || [];
  if (ivs.includes("comprehensive")) {
    mods = { ...COMPREHENSIVE.mods };
  } else {
    mods = composeMods(ivs.filter((k) => k in INTERVENTIONS));
  }
  for (const [k, v] of Object.entries(state.customMods || {})) {
    mods[k] = (mods[k] ?? 1) * v;
  }
  return mods;
}

/** Effective parameters after applying active interventions/custom modifiers. */
export function effectiveParams(state) {
  const mods = activeMods(state);
  return Object.keys(mods).length ? applyMods(state.params, mods) : { ...state.params };
}

/** Full derived scenario: effective params, trajectory, summary, R0 values. */
export function computeScenario(state) {
  const p = effectiveParams(state);
  const init = baselineInit();
  const traj = simulate(p, init, state.horizonYears || 20);
  const summary = summarise(traj);
  const trio = r0Trio(p);
  return {
    params: p, init, traj, summary,
    r0: { ...trio, full: r0Full(p, 1) },
    hasIntervention: Object.keys(activeMods(state)).length > 0,
  };
}

/** State vector at a given month index of a trajectory. */
export function stateAt(traj, i) {
  const idx = Math.max(0, Math.min(i, traj.I.length - 1));
  return { S: traj.S[idx], A: traj.A[idx], I: traj.I[idx], R: traj.R[idx] };
}
