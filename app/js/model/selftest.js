// selftest.js — Automated acceptance tests, surfaced in the UI.
// -----------------------------------------------------------------------------
// Makes correctness VISIBLE: asserts the thesis acceptance targets against the
// live model code and returns PASS/FAIL rows. Run at startup and on demand.

import { simulate, summarise, transitionMatrix, rowsSumToOne, foiSA } from "./sair.js";
import { r0Trio } from "./r0.js";
import { baselineParams, baselineInit } from "./params.js";
import { runScenarioTable } from "./interventions.js";

/** Run all acceptance checks; returns { pass, results:[{name, got, want, pass, tol}] }. */
export function runSelfTests() {
  const p = baselineParams();
  const init = baselineInit();
  const traj = simulate(p, init, 20);
  const s = summarise(traj);
  const r0 = r0Trio(p);
  const scen = runScenarioTable(p, init, 20);
  const comp = scen.find((r) => r.key === "comprehensive");

  const results = [];
  const add = (name, got, want, tol, unit = "", detail = "") =>
    results.push({ name, got, want, tol, pass: Math.abs(got - want) <= tol, unit, detail });

  add("R₀ — NGM operative (S₀ = 0.823)", r0.operative, 2.50, 0.05, "", "spectral radius of F·V⁻¹ on {A,I}");
  add("R₀ — NGM basic (S₀ = 1)", r0.basic, 3.04, 0.05, "", "next-generation matrix at fully susceptible");
  add("R₀ — crude β·k·D (D = 1/γ₂)", r0.crude, 2.71, 0.05, "", "well-mixed upper estimate");
  add("Peak addicted prevalence", s.peakPrevalence, 0.14, 0.01, "%", "baseline S–A–I–R run");
  add("Peak timing", s.peakYear, 7.0, 1.0, " yr", "year of maximum I(t)");
  add("Endemic prevalence (t = 20 yr)", s.endemicPrevalence, 0.08, 0.01, "%", "sustained by relapse δ");
  add("Comprehensive intervention peak reduction", comp.pctReduction, 0.63, 0.05, "%", "combined levers vs baseline");

  // Structural invariants (boolean checks encoded as 1/1).
  const { M } = transitionMatrix(init, p);
  results.push({ name: "Transition-matrix rows sum to 1", got: rowsSumToOne(M) ? 1 : 0, want: 1, tol: 0, pass: rowsSumToOne(M), unit: "", detail: "each row is a valid probability vector", bool: true });
  const foi = foiSA({ S: 0.5, A: 0.9, I: 0.9 }, { ...p, beta: 5, k: 5 });
  results.push({ name: "Force of infection clamped ≤ 1", got: foi <= 1 ? 1 : 0, want: 1, tol: 0, pass: foi <= 1, unit: "", detail: "probabilities constrained to [0,1]", bool: true });

  const allStateSums = traj.S.every((_, i) => Math.abs(traj.S[i] + traj.A[i] + traj.I[i] + traj.R[i] - 1) < 1e-6);
  results.push({ name: "State vector normalised each step", got: allStateSums ? 1 : 0, want: 1, tol: 0, pass: allStateSums, unit: "", detail: "S+A+I+R = 1 for all t", bool: true });

  return { pass: results.every((r) => r.pass), results };
}
