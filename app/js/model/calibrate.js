// calibrate.js: exact-solution calibration to the two HBSC figures (thesis Chapter 7; python/dpm/calibrate.py).
// -----------------------------------------------------------------------------
// Two published figures (7% in 2018, 11% in 2022) identify at most two quantities. For given values of the
// other rates, two quantities are solved exactly:
//   Scenario P (persistent shift): equilibrium in 2018 at exposure b0; exposure b1 from 2018 so that I(2022) = 11%.
//   Scenario T (temporary shock):  exposure stays b0; non-social acquisition rises by Δα during 2020 and 2021.
// Time t is measured in years from 2018.

import { equilibrium, simulate, prevalencePath } from "./sair.js";

export const TARGETS = { y2018: 0.07, y2022: 0.11 };
export const BMAX = 200, DMAX = 20;
export const shock = (dAlpha) => ({ dAlpha, start: 2, end: 4 });

/** Brent's method on [a, b] with f(a), f(b) of opposite sign. */
export function brent(f, a, b, tol = 1e-10, maxIter = 200) {
  let fa = f(a), fb = f(b);
  if (fa * fb > 0) return null;
  if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a]; [fa, fb] = [fb, fa]; }
  let c = a, fc = fa, d = b - a, mflag = true;
  for (let i = 0; i < maxIter; i++) {
    if (fb === 0 || Math.abs(b - a) < tol) return b;
    let s;
    if (fa !== fc && fb !== fc) {
      s = (a * fb * fc) / ((fa - fb) * (fa - fc)) + (b * fa * fc) / ((fb - fa) * (fb - fc)) + (c * fa * fb) / ((fc - fa) * (fc - fb));
    } else s = b - (fb * (b - a)) / (fb - fa);
    const cond = (s - (3 * a + b) / 4) * (s - b) >= 0 || (mflag && Math.abs(s - b) >= Math.abs(b - c) / 2) ||
      (!mflag && Math.abs(s - b) >= Math.abs(c - d) / 2) || (mflag && Math.abs(b - c) < tol) || (!mflag && Math.abs(c - d) < tol);
    if (cond) { s = (a + b) / 2; mflag = true; } else mflag = false;
    const fs = f(s);
    d = c; c = b; fc = fb;
    if (fa * fs < 0) { b = s; fb = fs; } else { a = s; fa = fs; }
    if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a]; [fa, fb] = [fb, fa]; }
  }
  return b;
}

const I4 = (p, x0, opts) => simulate(p, x0, [4], opts)[0][2];

export function solveB0(base, target = TARGETS.y2018) {
  return brent((b) => equilibrium({ ...base, b })[2] - target, 0, BMAX);
}
export function solveB1(p0, x0, target = TARGETS.y2022) {
  return brent((b) => I4({ ...p0, b }, x0) - target, p0.b, BMAX, 1e-9);
}
export function solveDalpha(p0, x0, target = TARGETS.y2022) {
  return brent((d) => I4(p0, x0, { shock: shock(d) }) - target, 0, DMAX, 1e-10);
}

/**
 * Calibrate both scenarios for the rates in `rates` (η, α, σ, γ₁, γ₂, δ, ν).
 * Returns { ok, reason, b0, b1, dAlpha, x18 } with the solved quantities.
 */
export function calibrateRates(rates, targets = TARGETS) {
  const base = { ...rates, b: 0 };
  const b0 = solveB0(base, targets.y2018);
  if (b0 === null) return { ok: false, reason: "The 2018 level cannot be reached with these rates." };
  const p0 = { ...base, b: b0 };
  const x18 = equilibrium(p0);
  const b1 = solveB1(p0, x18, targets.y2022);
  const dAlpha = solveDalpha(p0, x18, targets.y2022);
  let reason = null;
  if (b1 === null && dAlpha === null) reason = "The 2022 level cannot be reached in either scenario: recovery is too fast for prevalence to rise from 7% to 11% in four years.";
  else if (b1 === null) reason = "The 2022 level cannot be reached by a persistent shift in exposure.";
  else if (dAlpha === null) reason = "The 2022 level cannot be reached by a temporary shock.";
  return { ok: reason === null, reason, b0, b1, dAlpha, x18 };
}

/** Prevalence 2018..2018+years under scenario "P" or "T" for a calibrated set. */
export function scenarioPath(rates, cal, scenario, years = 24) {
  const p0 = { ...rates, b: cal.b0 };
  return scenario === "P" ? prevalencePath({ ...p0, b: cal.b1 }, cal.x18, years)
    : prevalencePath(p0, cal.x18, years, { shock: shock(cal.dAlpha) });
}
