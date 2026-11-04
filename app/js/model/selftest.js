// selftest.js: checks that the lab's JavaScript model reproduces the thesis's Python model.
// -----------------------------------------------------------------------------
// Reference values come from lab-data.js (python/export_lab.py). Each check compares the JavaScript
// computation at the thesis's central parameters with the Python value, to a stated tolerance.

import { LAB } from "./lab-data.js";
import { R0, R0Factors, equilibrium, floor, isStable, elasticities } from "./sair.js";
import { calibrateRates, scenarioPath } from "./calibrate.js";
import { careerTransmissibility, betaFromB } from "./network.js";
import { CENTRAL, DRAWS } from "./params.js";
import { runDraw } from "./interventions.js";

const rel = (a, b) => Math.abs(a - b) / Math.max(1e-12, Math.abs(b));

export function runSelfTests() {
  const ref = LAB.ref, out = [];
  const add = (name, got, want, tol, detail, fmt = (v) => v.toFixed(4)) =>
    out.push({ name, got: fmt(got), want: fmt(want), tol: `${tol.toExponential(0)} relative`, pass: rel(got, want) <= tol, detail });
  const addBool = (name, ok, detail) => out.push({ name, got: ok ? "yes" : "no", want: "yes", tol: "", pass: ok, detail });

  const p0 = { ...CENTRAL.rates, b: CENTRAL.b0 }, p1 = { ...CENTRAL.rates, b: CENTRAL.b1 };
  add("R₀ before 2018 (α = 0)", R0({ ...p0, alpha: 0 }), ref.R0_b0, 1e-5, "Closed form (5.3)");
  add("R₀ after 2018 (α = 0)", R0({ ...p1, alpha: 0 }), ref.R0_b1, 1e-5, "Closed form (5.3)");
  add("Relapse multiplier", R0Factors(p1).multiplier, ref.relapse_multiplier, 1e-5, "1/(1 − q)");
  add("Equilibrium prevalence after 2018", equilibrium(p1)[2], ref.eq_b1[2], 1e-5, "Quadratic (5.9)");
  add("Equilibrium prevalence with α = 0", equilibrium({ ...p1, alpha: 0 })[2], ref.eq_b1_alpha0[2], 1e-5, "S* = 1/R₀");
  add("Floor of prevalence at b = 0", floor(p1), ref.floor, 1e-5, "Equation (5.10)");
  addBool("Interior equilibrium is locally stable", isStable(p1, equilibrium(p1)), "Routh–Hurwitz (Theorem 5.9)");

  const el = elasticities(p1, "I");
  add("Sensitivity of I* to recovery γ₂", el.gamma2, LAB.results.elast_I_gamma2, 1e-4, "Normalised index");
  add("Sensitivity of I* to exposure b", el.b, LAB.results.elast_I_b, 1e-4, "Normalised index");
  const cal = calibrateRates(CENTRAL.rates);
  add("Solved b₀ at central rates", cal.b0, CENTRAL.b0, 1e-5, "Exact solution, 7% in 2018");
  add("Solved b₁ at central rates", cal.b1, CENTRAL.b1, 1e-4, "Exact solution, 11% in 2022");
  add("Solved Δα at central rates", cal.dAlpha, CENTRAL.dAlpha, 1e-4, "Exact solution, 11% in 2022");
  const pP = scenarioPath(CENTRAL.rates, { b0: CENTRAL.b0, b1: CENTRAL.b1, dAlpha: CENTRAL.dAlpha, x18: CENTRAL.x18 }, "P");
  const pT = scenarioPath(CENTRAL.rates, { b0: CENTRAL.b0, b1: CENTRAL.b1, dAlpha: CENTRAL.dAlpha, x18: CENTRAL.x18 }, "T");
  add("Persistent shift: prevalence in 2042", pP[24], ref.I_P[24], 1e-5, "ODE integration", (v) => `${(v * 100).toFixed(3)}%`);
  add("Temporary shock: prevalence in 2026", pT[8], ref.I_T[8], 1e-5, "ODE with step change in α", (v) => `${(v * 100).toFixed(3)}%`);
  add("Temporary shock: prevalence in 2042", pT[24], ref.I_T[24], 1e-5, "ODE with step change in α", (v) => `${(v * 100).toFixed(3)}%`);

  // projection medians across all kept draws
  const yrs = [8, 24], P = yrs.map(() => []), T = yrs.map(() => []);
  for (const d of DRAWS) {
    const a = scenarioPath(d.rates, d.cal, "P"), b = scenarioPath(d.rates, d.cal, "T");
    yrs.forEach((y, i) => { P[i].push(a[y]); T[i].push(b[y]); });
  }
  const med = (v) => { const s = v.slice().sort((x, y) => x - y); const m = s.length >> 1; return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; };
  yrs.forEach((y, i) => {
    add(`Median across ${DRAWS.length} draws, persistent shift, ${2018 + y}`, med(P[i]), LAB.projections.P[50][y], 1e-3, "Projection percentiles", (v) => `${(v * 100).toFixed(2)}%`);
    add(`Median across ${DRAWS.length} draws, temporary shock, ${2018 + y}`, med(T[i]), LAB.projections.T[50][y], 1e-3, "Projection percentiles", (v) => `${(v * 100).toFixed(2)}%`);
  });

  const centralDraw = { rates: CENTRAL.rates, cal: { b0: CENTRAL.b0, b1: CENTRAL.b1, dAlpha: CENTRAL.dAlpha, x18: CENTRAL.x18 } };
  const mults = [{ b: 0.75 }, { b: 0.8, alpha: 0.8 }, { gamma2: 1.6, delta: 0.85 }];
  add("Three mechanisms together, persistent shift, 2042", runDraw(centralDraw, "P", mults).at(-1), ref.interv_P_2042, 1e-5, "Interventions from 2024", (v) => `${(v * 100).toFixed(3)}%`);
  add("Three mechanisms together, temporary shock, 2042", runDraw(centralDraw, "T", mults).at(-1), ref.interv_T_2042, 1e-5, "Interventions from 2024", (v) => `${(v * 100).toFixed(3)}%`);
  const beta = betaFromB(p1, 10);
  add("Per-edge monthly transmission at k = 10", beta, ref.beta_k10, 1e-5, "β from b");
  add("Career transmissibility at k = 10", careerTransmissibility({ ...p1, alpha: 0 }, beta), ref.T_career_k10, 1e-5, "Chapter 8");

  return { pass: out.every((r) => r.pass), results: out };
}
