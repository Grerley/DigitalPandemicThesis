// interventions.js: intervention mechanisms under assumed efficacies (thesis Chapter 9; python/run_all.py).
// -----------------------------------------------------------------------------
// Each mechanism multiplies one or more rates from 2024. Efficacy ranges are assumptions, taken from the thesis.
// Effects are statements about the relative leverage of mechanisms inside the model, not predictions.

import { LAB } from "./lab-data.js";
import { simulate } from "./sair.js";
import { shock } from "./calibrate.js";
import { RNG } from "./rng.js";

export const SPEC = LAB.interventions;
export const KEYS = ["platform", "school", "minimum_age", "treatment"];
export const START = SPEC.start;          // 2024
export const END = 2042;

/** Multipliers for one mechanism at efficacy e (and e2 for the second effect of treatment). */
export function multipliers(key, e, e2 = 0) {
  if (key === "platform") return { b: 1 - e };
  if (key === "school") return { gamma1: 1 + e };
  if (key === "minimum_age") return { b: 1 - 0.4 * e, alpha: 1 - 0.4 * e };
  if (key === "treatment") return { gamma2: 1 + e, delta: 1 - e2 };
  throw new Error(key);
}

export function applyAll(p, mults) {
  const q = { ...p };
  for (const m of mults) for (const [k, v] of Object.entries(m)) q[k] *= v;
  return q;
}

/**
 * Prevalence 2024..2042 for a calibrated draw under scenario P or T, with the given multipliers applied from 2024.
 * The state in 2024 comes from the calibrated path (in T the 2020 to 2021 shock has ended by then).
 */
export function runDraw(draw, scenario, mults) {
  const p0 = { ...draw.rates, b: draw.cal.b0 };
  const pre = scenario === "P" ? { ...p0, b: draw.cal.b1 } : p0;
  const opts = scenario === "T" ? { shock: shock(draw.cal.dAlpha) } : {};
  const x24 = simulate(pre, draw.cal.x18, [START - 2018], opts)[0];
  const q = applyAll(pre, mults);
  const t = Array.from({ length: END - START + 1 }, (_, i) => i);
  return simulate(q, x24, t).map((x) => x[2]);
}

/** Random efficacy draw for a mechanism, as in the thesis. */
export function drawEfficacy(key, rng) {
  const r = SPEC[key].range;
  const e = r[0] + rng.next() * (r[1] - r[0]);
  if (key === "treatment") { const r2 = SPEC.treatment.range2; return multipliers(key, e, r2[0] + rng.next() * (r2[1] - r2[0])); }
  return multipliers(key, e);
}

/**
 * Relative reduction in 2042 prevalence across draws for a combination of mechanisms, with efficacies drawn
 * uniformly from the assumed ranges (one draw per mechanism per parameter set).
 */
export function acrossDraws(draws, combo, scenario, seed = 20261031) {
  const rng = new RNG(seed), red = [];
  for (const d of draws) {
    const mults = combo.map((k) => drawEfficacy(k, rng));
    const base = runDraw(d, scenario, []);
    const with_ = runDraw(d, scenario, mults);
    red.push(1 - with_[with_.length - 1] / base[base.length - 1]);
  }
  red.sort((a, b) => a - b);
  const q = (f) => red[Math.min(red.length - 1, Math.floor(f * red.length))];
  return { median: q(0.5), p05: q(0.05), p95: q(0.95), n: red.length };
}
