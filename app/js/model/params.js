// params.js — Calibrated parameters, metadata, initial state, presets.
// -----------------------------------------------------------------------------
// Single source of truth for the S–A–I–R model parameters, mirroring
// R/utils.R (default_params / default_init) and data/parameters.csv. All rates
// are per MONTH. These values are calibrated to reproduce the thesis targets:
//   R0 (NGM, operative) ≈ 2.50, peak addicted prevalence ≈ 14% at ~year 7,
//   endemic prevalence ≈ 8%.

/** Calibrated baseline parameters (monthly rates). */
export const DEFAULT_PARAMS = Object.freeze({
  beta: 0.033, // per-contact transmission probability / month
  k: 3.7, // effective mean-field transmitting contacts
  eta: 0.5, // reduced infectivity weight of At-Risk vs Addicted
  mu: 0.0015, // exogenous (non-social) acquisition rate / month
  sigma: 0.02, // At-Risk -> Addicted progression / month
  gamma1: 0.018, // At-Risk -> Susceptible remission / month
  gamma2: 0.045, // Addicted -> Recovered recovery / month
  delta: 0.006, // Recovered -> At-Risk relapse / month
});

/** Initial state distribution [S, A, I, R] (HBSC 2022 baseline). */
export const INITIAL_STATE = Object.freeze({
  S: 0.823,
  A: 0.065,
  I: 0.047,
  R: 0.065,
});

/**
 * Slider metadata for the Parameter Explorer. `min/max` bound the slider;
 * `step` is the granularity; `label` is display text; `sym` is the symbol;
 * `unit` a short unit; `desc` a one-line description; `source` provenance.
 */
export const PARAM_META = Object.freeze({
  beta: { sym: "β", label: "Transmission", min: 0.0, max: 0.08, step: 0.001, unit: "/contact·mo", desc: "Per-contact social-influence transmission probability.", source: "HBSC Bayesian calibration (Ch 6.5)" },
  k: { sym: "k", label: "Effective contacts", min: 0.5, max: 10, step: 0.1, unit: "contacts", desc: "Effective mean-field transmitting contacts.", source: "Mean-field calibration (Ch 6.5.3)" },
  eta: { sym: "η", label: "At-Risk infectivity", min: 0.0, max: 1.0, step: 0.01, unit: "—", desc: "Reduced infectivity of At-Risk relative to Addicted.", source: "Assumed (Ch 4)" },
  mu: { sym: "μ", label: "Exogenous acquisition", min: 0.0, max: 0.01, step: 0.0001, unit: "/mo", desc: "Non-social (exogenous) acquisition rate.", source: "Calibration" },
  sigma: { sym: "σ", label: "Progression A→I", min: 0.0, max: 0.1, step: 0.001, unit: "/mo", desc: "At-Risk → Addicted progression rate.", source: "HMM / survival (Ch 6.4)" },
  gamma1: { sym: "γ₁", label: "Remission A→S", min: 0.0, max: 0.1, step: 0.001, unit: "/mo", desc: "At-Risk → Susceptible remission rate.", source: "HMM (Ch 6.2)" },
  gamma2: { sym: "γ₂", label: "Recovery I→R", min: 0.005, max: 0.15, step: 0.001, unit: "/mo", desc: "Addicted → Recovered recovery rate (mean ~1.9 yr).", source: "HMM + survival + Winkler 2013" },
  delta: { sym: "δ", label: "Relapse R→A", min: 0.0, max: 0.05, step: 0.0005, unit: "/mo", desc: "Recovered → At-Risk relapse rate.", source: "Multi-state model + clinical lit" },
});

export const PARAM_KEYS = Object.keys(DEFAULT_PARAMS);

/** Return a fresh mutable copy of the calibrated baseline parameters. */
export function baselineParams() {
  return { ...DEFAULT_PARAMS };
}

/** Return a fresh mutable copy of the initial state. */
export function baselineInit() {
  return { ...INITIAL_STATE };
}

/** Compartment display metadata. Colour-blind-safe, print-legible. */
export const COMPARTMENTS = Object.freeze([
  { key: "S", name: "Susceptible", color: "var(--c-s)", dash: "none", desc: "Not problematically engaged; at risk of exposure." },
  { key: "A", name: "At-Risk", color: "var(--c-a)", dash: "6 3", desc: "Heavy/pre-clinical use; elevated risk; partially transmitting." },
  { key: "I", name: "Addicted", color: "var(--c-i)", dash: "none", desc: "Meets problematic-use thresholds (operational construct)." },
  { key: "R", name: "Recovered", color: "var(--c-r)", dash: "2 3", desc: "Reduced use; NOT absorbing — subject to relapse." },
]);
