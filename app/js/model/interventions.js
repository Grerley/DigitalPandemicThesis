// interventions.js — Intervention scenarios (Chapter 8).
// -----------------------------------------------------------------------------
// Port of R/07_interventions.R. Interventions are multiplicative modifications
// to model parameters. Efficacies are ASSUMED, informed by the direction and
// rough magnitude of published evaluations — NOT estimated here. Outputs are
// therefore comparative, not predictive.

import { simulate, summarise } from "./sair.js";

/**
 * Named interventions. Each maps a parameter to a multiplier. `comprehensive`
 * combines structural + behavioural + treatment levers and is tuned to land
 * near a −63% peak reduction.
 */
export const INTERVENTIONS = Object.freeze({
  education: { label: "Education", mods: { beta: 0.82 }, blurb: "Media-literacy / behavioural programmes reduce per-contact transmission (−18% β).", cost: 120 },
  regulation: { label: "Platform regulation", mods: { beta: 0.71 }, blurb: "Structural platform-design constraints reduce transmission (−29% β).", cost: 65 },
  age_gating: { label: "Age-gating", mods: { mu: 0.5, beta: 0.9 }, blurb: "Verified age limits cut exogenous entry (−50% μ) and some transmission (−10% β).", cost: 40 },
  treatment: { label: "Treatment access", mods: { gamma2: 1.4, delta: 0.7 }, blurb: "Clinical services raise recovery (+40% γ₂) and reduce relapse (−30% δ).", cost: 210 },
});

/** The comprehensive bundle (combined levers). Tuned to ≈ −63% peak. */
export const COMPREHENSIVE = Object.freeze({
  label: "Comprehensive",
  mods: { beta: 0.55, mu: 0.5, gamma2: 1.4, delta: 0.7 },
  blurb: "All levers combined: structural + behavioural transmission reduction, age-gating, and expanded treatment.",
  cost: 380,
});

/** The (illustrative) unit-cost basis for the ZAR cost-effectiveness framing. */
export const COST_ASSUMPTIONS = Object.freeze({
  populationYouth: 12_500_000, // approx SA population aged 10–24 (illustrative)
  currency: "ZAR",
  note: "Costs are illustrative per-youth-per-year figures (ZAR) for comparative framing only — not a costed programme budget.",
});

/**
 * Apply a set of multiplicative modifiers to a parameter vector.
 * @param {object} p base parameters
 * @param {Record<string, number>} mods parameter → multiplier
 */
export function applyMods(p, mods) {
  const out = { ...p };
  for (const key of Object.keys(mods)) {
    if (key in out) out[key] = p[key] * mods[key];
  }
  return out;
}

/**
 * Compose several named interventions into a single modifier map. When two
 * interventions touch the same parameter, their multipliers multiply.
 * @param {string[]} names keys of INTERVENTIONS
 */
export function composeMods(names) {
  const mods = {};
  for (const name of names) {
    const spec = INTERVENTIONS[name];
    if (!spec) continue;
    for (const [k, v] of Object.entries(spec.mods)) {
      mods[k] = (mods[k] ?? 1) * v;
    }
  }
  return mods;
}

/**
 * Run every named scenario plus baseline and comprehensive, returning a table
 * of peak prevalence, peak year, endemic level, and % peak reduction.
 */
export function runScenarioTable(p, init, horizonYears = 20) {
  const base = simulate(p, init, horizonYears);
  const baseSum = summarise(base);
  const rows = [
    { key: "baseline", label: "Baseline", ...scenarioSummary(baseSum, baseSum) },
  ];
  for (const [key, spec] of Object.entries(INTERVENTIONS)) {
    const s = summarise(simulate(applyMods(p, spec.mods), init, horizonYears));
    rows.push({ key, label: spec.label, cost: spec.cost, ...scenarioSummary(s, baseSum) });
  }
  const comp = summarise(simulate(applyMods(p, COMPREHENSIVE.mods), init, horizonYears));
  rows.push({ key: "comprehensive", label: COMPREHENSIVE.label, cost: COMPREHENSIVE.cost, ...scenarioSummary(comp, baseSum) });
  return rows;
}

function scenarioSummary(s, baseSum) {
  const reduction = (baseSum.peakPrevalence - s.peakPrevalence) / baseSum.peakPrevalence;
  return {
    peakPrevalence: s.peakPrevalence,
    peakYear: s.peakYear,
    endemicPrevalence: s.endemicPrevalence,
    pctReduction: reduction,
  };
}
