// cohort.js — Synthetic longitudinal cohort + statistical layer.
// -----------------------------------------------------------------------------
// Port of R/03_synthetic_data.R + R/05_regression.R. A MODELLING INSTRUMENT,
// not observed data: it produces individual trajectories that reproduce known
// prevalence and risk-factor odds ratios, enabling (a) an odds-ratio forest
// plot recovered by logistic regression, (b) a start→end state Sankey, and
// (c) a "known vs estimated" transition-matrix heat-map. Fully seeded.

import { RNG } from "./rng.js";
import { rateToProb } from "./sair.js";
import { INITIAL_STATE } from "./params.js";

/** Published odds ratios that tilt the individual S→A probability. */
export const RISK_FACTORS = Object.freeze([
  { key: "age", label: "Age (per year)", or: 0.94, kind: "cont" },
  { key: "female", label: "Female", or: 1.38, kind: "bin" },
  { key: "urban", label: "Urban residence", or: 2.18, kind: "bin" },
  { key: "low_ses", label: "Low SES", or: 1.45, kind: "bin" },
  { key: "high_depression", label: "High depression", or: 2.84, kind: "bin" },
  { key: "high_literacy", label: "High digital literacy", or: 0.72, kind: "bin" },
]);

const STATES = ["S", "A", "I", "R"];

/** annual transition probability from a monthly rate. */
function annualProb(monthlyRate) {
  return 1 - Math.exp(-12 * monthlyRate);
}

function logit(x) { return Math.log(x / (1 - x)); }
function expit(x) { return 1 / (1 + Math.exp(-x)); }

/**
 * Generate a synthetic cohort and derive statistical summaries.
 * @param {object} p parameters
 * @param {object} opts { n, years, seed }
 */
export function generateCohort(p, opts = {}) {
  const { n = 6000, years = 10, seed = 42 } = opts;
  const rng = new RNG(seed);

  // Covariates.
  const age = new Float64Array(n), female = new Int8Array(n), urban = new Int8Array(n),
        lowSes = new Int8Array(n), highDep = new Int8Array(n), highLit = new Int8Array(n);
  for (let i = 0; i < n; i++) {
    age[i] = 10 + rng.int(11); // 10..20
    female[i] = rng.bernoulli(0.51);
    urban[i] = rng.bernoulli(0.68);
    lowSes[i] = rng.bernoulli(0.33);
    highDep[i] = rng.bernoulli(0.20);
    highLit[i] = rng.bernoulli(0.40);
  }

  // Initial states drawn from baseline distribution.
  const init = INITIAL_STATE;
  const state = new Array(n);
  const startState = new Int8Array(n);
  for (let i = 0; i < n; i++) {
    const s = rng.categorical([init.S, init.A, init.I, init.R]);
    state[i] = s; startState[i] = s;
  }

  const pAI = annualProb(p.sigma), pAS = annualProb(p.gamma1),
        pIR = annualProb(p.gamma2), pRA = annualProb(p.delta);

  const byYear = [];
  // Onset records for the risk-factor regression: for each person-year that
  // begins in Susceptible, covariates + whether they transitioned S→A. This
  // directly matches the generative logit, so the published ORs are recoverable.
  const reg = { rows: [], y: [] };
  // Empirical annual transition counts (estimated matrix).
  const transCount = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
  let sumPSA = 0, nPSAsteps = 0;

  const covOddsShift = (i) =>
    Math.log(0.94) * (age[i] - 14) +
    Math.log(1.38) * female[i] +
    Math.log(2.18) * urban[i] +
    Math.log(1.45) * lowSes[i] +
    Math.log(2.84) * highDep[i] +
    Math.log(0.72) * highLit[i];

  for (let yr = 0; yr <= years; yr++) {
    // Record prevalence this year.
    const cnt = [0, 0, 0, 0];
    for (let i = 0; i < n; i++) cnt[state[i]]++;
    byYear.push({ year: yr, S: cnt[0] / n, A: cnt[1] / n, I: cnt[2] / n, R: cnt[3] / n });
    if (yr === years) break;

    // Cohort-average base S→A probability this step (mass-action, annual).
    const prev = byYear[byYear.length - 1];
    const foiRate = p.beta * p.k * (prev.I + p.eta * prev.A);
    let basePSA = 1 - Math.exp(-12 * foiRate) + p.mu;
    basePSA = Math.min(Math.max(basePSA, 1e-4), 0.9);
    sumPSA += basePSA; nPSAsteps++;

    const nextState = new Array(n);
    for (let i = 0; i < n; i++) {
      const s = state[i];
      let ns = s;
      if (s === 0) {
        const pSA_i = expit(logit(basePSA) + covOddsShift(i));
        const onset = rng.next() < pSA_i ? 1 : 0;
        if (onset) ns = 1;
        // Record the onset event for the risk-factor regression.
        reg.rows.push([age[i] - 14, female[i], urban[i], lowSes[i], highDep[i], highLit[i]]);
        reg.y.push(onset);
      } else if (s === 1) {
        const u = rng.next();
        if (u < pAI) ns = 2; else if (u < pAI + pAS) ns = 0;
      } else if (s === 2) {
        if (rng.next() < pIR) ns = 3;
      } else if (s === 3) {
        if (rng.next() < pRA) ns = 1;
      }
      transCount[s][ns]++;
      nextState[i] = ns;
    }
    for (let i = 0; i < n; i++) state[i] = nextState[i];
  }

  // Sankey: start (year 0) → end (final year) flows.
  const flows = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
  for (let i = 0; i < n; i++) flows[startState[i]][state[i]]++;

  // Estimated (empirical) annual transition matrix from counts.
  const estMatrix = transCount.map((row) => {
    const tot = row.reduce((a, b) => a + b, 0) || 1;
    return row.map((v) => v / tot);
  });
  // Known generating annual transition matrix (S-row uses cohort-avg base P_SA).
  const avgPSA = sumPSA / nPSAsteps;
  const knownMatrix = [
    [1 - avgPSA, avgPSA, 0, 0],
    [pAS, 1 - pAS - pAI, pAI, 0],
    [0, 0, 1 - pIR, pIR],
    [0, pRA, 0, 1 - pRA],
  ];

  // Logistic regression to recover odds ratios (forest plot).
  const fit = logisticRegression(reg.rows, reg.y);

  return {
    byYear, flows, estMatrix, knownMatrix, states: STATES,
    oddsRatios: buildOddsTable(fit),
    matrixMAE: matMAE(knownMatrix, estMatrix),
    n, years,
  };
}

function buildOddsTable(fit) {
  // fit.beta[0] is intercept; predictors align with RISK_FACTORS order.
  return RISK_FACTORS.map((rf, i) => {
    const b = fit.beta[i + 1], se = fit.se[i + 1];
    return {
      key: rf.key, label: rf.label,
      knownOR: rf.or,
      estOR: Math.exp(b),
      lo: Math.exp(b - 1.96 * se),
      hi: Math.exp(b + 1.96 * se),
    };
  });
}

function matMAE(a, b) {
  let s = 0, c = 0;
  for (let i = 0; i < a.length; i++)
    for (let j = 0; j < a[i].length; j++) { s += Math.abs(a[i][j] - b[i][j]); c++; }
  return s / c;
}

/**
 * Binary logistic regression via iteratively reweighted least squares (IRLS).
 * Returns coefficients and standard errors (sqrt of diagonal of (XᵀWX)⁻¹).
 * @param {number[][]} X rows of predictors (without intercept)
 * @param {number[]} y binary outcomes
 */
export function logisticRegression(X, y, maxIter = 30, tol = 1e-8) {
  const n = X.length, k = X[0].length + 1;
  // Design with intercept.
  const D = X.map((row) => [1, ...row]);
  let beta = new Array(k).fill(0);
  let cov = identity(k);
  for (let iter = 0; iter < maxIter; iter++) {
    const XtWX = zeros(k, k);
    const XtWz = new Array(k).fill(0);
    for (let i = 0; i < n; i++) {
      let eta = 0;
      for (let j = 0; j < k; j++) eta += D[i][j] * beta[j];
      const mu = expit(eta);
      const w = Math.max(mu * (1 - mu), 1e-8);
      const z = eta + (y[i] - mu) / w;
      for (let a = 0; a < k; a++) {
        XtWz[a] += D[i][a] * w * z;
        for (let b = 0; b < k; b++) XtWX[a][b] += D[i][a] * w * D[i][b];
      }
    }
    const inv = invert(XtWX);
    const newBeta = matVec(inv, XtWz);
    let delta = 0;
    for (let j = 0; j < k; j++) delta += Math.abs(newBeta[j] - beta[j]);
    beta = newBeta; cov = inv;
    if (delta < tol) break;
  }
  const se = cov.map((row, i) => Math.sqrt(Math.max(row[i], 0)));
  return { beta, se };
}

// ---- tiny linear algebra ----------------------------------------------------
function zeros(r, c) { return Array.from({ length: r }, () => new Array(c).fill(0)); }
function identity(k) { const m = zeros(k, k); for (let i = 0; i < k; i++) m[i][i] = 1; return m; }
function matVec(M, v) { return M.map((row) => row.reduce((s, x, j) => s + x * v[j], 0)); }
function invert(A) {
  const n = A.length;
  const M = A.map((row, i) => [...row, ...identity(n)[i]]);
  for (let col = 0; col < n; col++) {
    let piv = col;
    for (let r = col + 1; r < n; r++) if (Math.abs(M[r][col]) > Math.abs(M[piv][col])) piv = r;
    [M[col], M[piv]] = [M[piv], M[col]];
    const d = M[col][col] || 1e-12;
    for (let j = 0; j < 2 * n; j++) M[col][j] /= d;
    for (let r = 0; r < n; r++) {
      if (r === col) continue;
      const f = M[r][col];
      for (let j = 0; j < 2 * n; j++) M[r][j] -= f * M[col][j];
    }
  }
  return M.map((row) => row.slice(n));
}
