// params.js: parameter ranges, central values and calibration draws, read from the thesis outputs.
// -----------------------------------------------------------------------------
// Every value here comes from lab-data.js, which python/export_lab.py generates from the thesis analysis.

import { LAB } from "./lab-data.js";

export const RATE_KEYS = ["eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"];

export const PARAM_META = {
  b: { sym: "b", label: "Social exposure rate", unit: "per year", desc: "Effective rate of acquisition through exposure to At-Risk and problematic users. Bundles influence, homophily and shared environment.", source: "Solved exactly from the HBSC figures" },
  eta: { sym: "η", label: "Influence of At-Risk contacts", unit: "relative", desc: "Exposure from an At-Risk young person relative to one with problematic use.", source: "Assumption" },
  alpha: { sym: "α", label: "Non-social acquisition", unit: "per year", desc: "Entry into risk that does not depend on contacts, for example through platform design.", source: "Assumption" },
  sigma: { sym: "σ", label: "Progression, At-Risk to problematic use", unit: "per year", desc: "Rate at which At-Risk young people reach problematic use.", source: "Xin et al. (2026); Chiang et al. (2019)" },
  gamma1: { sym: "γ₁", label: "Remission from At-Risk", unit: "per year", desc: "Rate at which At-Risk young people return to no risk.", source: "Xin et al. (2026); Chiang et al. (2019)" },
  gamma2: { sym: "γ₂", label: "Recovery from problematic use", unit: "per year", desc: "Rate of leaving problematic use.", source: "Ko et al. (2007); Wartberg et al. (2019, 2020); Xin et al. (2026); Malygin et al. (2025)" },
  delta: { sym: "δ", label: "Relapse into risk", unit: "per year", desc: "Rate at which recovered young people return to the At-Risk state.", source: "Jo et al. (2025)" },
  nu: { sym: "ν", label: "Turnover", unit: "per year", desc: "Entry at age 10 and exit at 25.", source: "Age band 10 to 24" },
};

export const RANGES = LAB.ranges;

/** Central rates and the two solved scenarios at those rates. */
export const CENTRAL = {
  rates: { ...LAB.central.params },
  b0: LAB.central.b0, b1: LAB.central.b1, dAlpha: LAB.central.d_alpha, x18: LAB.central.x18,
};

/** Central parameter set after 2018 under the persistent shift (b = b1). */
export const centralParams = () => ({ ...CENTRAL.rates, b: CENTRAL.b1 });

/** Kept calibration draws as objects. */
export const DRAWS = LAB.draws.map((row) => {
  const o = {};
  LAB.draw_keys.forEach((k, i) => { o[k] = row[i]; });
  o.rates = Object.fromEntries(RATE_KEYS.map((k) => [k, o[k]]));
  o.cal = { b0: o.b0, b1: o.b1, dAlpha: o.d_alpha, x18: [o.x18_S, o.x18_A, o.x18_I, o.x18_R] };
  return o;
});

/** A thesis result by key (values as written by the analysis code). */
export function R(key) {
  const v = LAB.results[key];
  if (v === undefined) throw new Error(`missing thesis result: ${key}`);
  return v;
}

export { LAB };
