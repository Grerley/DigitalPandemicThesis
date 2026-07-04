// r0.js — Reproduction number, equilibrium, and bifurcation.
// -----------------------------------------------------------------------------
// Exact port of R/06_r0_analysis.R. Provides the three R0 values reported in the
// thesis, the endemic equilibrium (as the fixed point of the deterministic
// map), and the transcritical bifurcation curve (equilibrium I* vs R0).

import { simulate } from "./sair.js";
import { baselineInit } from "./params.js";
import { INITIAL_STATE } from "./params.js";

/** Crude well-mixed R0 = β·k·D, with D = 1/γ₂ (mean months Addicted). */
export function r0Crude(p) {
  const D = 1 / p.gamma2;
  return p.beta * p.k * D;
}

/**
 * Next-generation-matrix R0 (Appendix A.4) on the infected subsystem {A, I},
 * evaluated at susceptible fraction S0.
 *   F = [[b·η·S0, b·S0],[0,0]],  V = [[σ+γ₁, 0],[−σ, γ₂]],  b = β·k
 *   R0 = ρ(F·V⁻¹) = b·S0·(η·γ₂ + σ) / ((σ+γ₁)·γ₂)
 * The matrix F·V⁻¹ is rank-1 with a single non-zero eigenvalue, so the
 * spectral radius has the closed form above.
 */
export function r0NGM(p, S0 = INITIAL_STATE.S) {
  const b = p.beta * p.k;
  return (b * S0 * (p.eta * p.gamma2 + p.sigma)) / ((p.sigma + p.gamma1) * p.gamma2);
}

/**
 * Full-system basic reproduction number on the {A, I, R} infected subsystem,
 * i.e. WITHOUT eliminating the relapse loop. In the closed model (no background
 * exit) every Recovered eventually relapses, so R recirculates individuals
 * through the transmitting compartments; the full NGM therefore governs the
 * true stability of the disease-free equilibrium and crosses 1 exactly at the
 * transcritical bifurcation. It is larger than the reduced (operative) R0 that
 * the thesis reports as its headline figure. Evaluated at susceptible S0.
 *   F = [[b·η·S0, b·S0, 0],[0,0,0],[0,0,0]],  b = β·k
 *   V = [[σ+γ₁, 0, −δ],[−σ, γ₂, 0],[0, −γ₂, δ]]
 */
export function r0Full(p, S0 = 1) {
  const b = p.beta * p.k;
  const V = [
    [p.sigma + p.gamma1, 0, -p.delta],
    [-p.sigma, p.gamma2, 0],
    [0, -p.gamma2, p.delta],
  ];
  const det =
    V[0][0] * (V[1][1] * V[2][2] - V[1][2] * V[2][1]) -
    V[0][1] * (V[1][0] * V[2][2] - V[1][2] * V[2][0]) +
    V[0][2] * (V[1][0] * V[2][1] - V[1][1] * V[2][0]);
  // Only the first row of F is non-zero, so ρ(F·V⁻¹) = (F·V⁻¹)[0][0].
  const inv00 = (V[1][1] * V[2][2] - V[1][2] * V[2][1]) / det;
  const inv10 = (V[1][2] * V[2][0] - V[1][0] * V[2][2]) / det;
  return b * p.eta * S0 * inv00 + b * S0 * inv10;
}

/** The three labelled R0 values reported in the thesis. */
export function r0Trio(p) {
  return {
    operative: r0NGM(p, INITIAL_STATE.S), // NGM at baseline S0 = 0.823 ≈ 2.50
    basic: r0NGM(p, 1), // NGM at S0 = 1 ≈ 3.04
    crude: r0Crude(p), // β·k·D ≈ 2.71
  };
}

/**
 * Networked R0 with degree heterogeneity (Appendix A.7):
 *   R0_net = R0_base · (⟨k²⟩ − ⟨k⟩) / ⟨k⟩,  R0_base = β / γ₂.
 */
export function r0Network(p, degrees) {
  const n = degrees.length;
  let k1 = 0, k2 = 0;
  for (const d of degrees) { k1 += d; k2 += d * d; }
  k1 /= n; k2 /= n;
  const base = p.beta / p.gamma2;
  return (base * (k2 - k1)) / k1;
}

/**
 * Endemic equilibrium as the fixed point of the deterministic map: iterate the
 * simulation to convergence and read off the stationary state. This is robust
 * and always consistent with the trajectory the user sees. Returns S*,A*,I*,R*
 * plus the operative R0 and a converged flag.
 */
export function endemicEquilibrium(p, init = baselineInit(), maxYears = 400, tol = 1e-10) {
  // Long-run integration; the closed system relaxes to its attractor.
  const traj = simulate(p, init, maxYears);
  const n = traj.I.length;
  const last = { S: traj.S[n - 1], A: traj.A[n - 1], I: traj.I[n - 1], R: traj.R[n - 1] };
  const prev = { S: traj.S[n - 2], A: traj.A[n - 2], I: traj.I[n - 2], R: traj.R[n - 2] };
  const drift = Math.max(
    Math.abs(last.S - prev.S), Math.abs(last.A - prev.A),
    Math.abs(last.I - prev.I), Math.abs(last.R - prev.R)
  );
  return { ...last, R0: r0NGM(p, INITIAL_STATE.S), converged: drift < tol * 1e3 };
}

/**
 * Transcritical bifurcation curve: equilibrium prevalence I* as a function of
 * R0. We sweep β (holding all else fixed) so that R0 ranges over [rMin, rMax],
 * and record I* at each. Below R0 = 1 the disease-free equilibrium is stable
 * (I* = 0); above it the positive endemic branch emerges.
 * @returns {{r0:number[], iStar:number[], threshold:number}}
 */
export function bifurcationCurve(p, rMin = 0.2, rMax = 4.0, steps = 90) {
  // The disease-free equilibrium sits at S = 1, so its stability is governed by
  // the FULL basic reproduction number (including the relapse recirculation
  // R→A). The transcritical bifurcation is at full R0 = 1, so we sweep and
  // label the x-axis by full R0. R0 scales linearly in β, so invert for β.
  const r0AtBase = r0Full(p, 1);
  const betaBase = p.beta;
  const r0 = [], iStar = [];
  for (let s = 0; s <= steps; s++) {
    const targetR0 = rMin + ((rMax - rMin) * s) / steps;
    const beta = betaBase * (targetR0 / r0AtBase);
    // Set μ = 0 so the disease-free equilibrium is genuinely disease-free: the
    // exogenous seeding term μ would otherwise sustain a small I even below the
    // threshold and blur the transcritical bifurcation.
    const pp = { ...p, beta, mu: 0 };
    const seed = targetR0 <= 1
      ? { S: 0.999, A: 0.0005, I: 0.0005, R: 0 }
      : { S: 0.9, A: 0.05, I: 0.03, R: 0.02 };
    const eq = endemicEquilibrium(pp, seed, 600);
    r0.push(targetR0);
    iStar.push(eq.I < 1e-4 ? 0 : eq.I);
  }
  return { r0, iStar, threshold: 1 };
}
