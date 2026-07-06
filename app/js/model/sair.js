// sair.js — Deterministic S–A–I–R difference-equation engine.
// -----------------------------------------------------------------------------
// Exact port of R/01_sair_model.R. The population state x = [S, A, I, R]
// (fractions summing to 1) evolves under a STATE-DEPENDENT transition matrix
// P(x_t): x(t+1) = P(x_t)^T x(t), with a monthly time step. The S→A entry is a
// mass-action force of infection, which makes the map nonlinear (this is what
// produces a rise–peak–decline epidemic curve, not a linear Markov chain).
//
// This module is PURE: no DOM, no globals. It can be audited independently of
// the UI and is exercised directly by the self-test panel.

/** Convert a per-step rate to a transition probability p = 1 − exp(−rate). */
export function rateToProb(rate) {
  return 1 - Math.exp(-rate);
}

/** Clamp x to [lo, hi]. */
export function clamp(x, lo = 0, hi = 1) {
  return x < lo ? lo : x > hi ? hi : x;
}

/**
 * Mass-action force of infection → S→A transition probability (Eq 4.3).
 * P_SA = min(1 − exp(−β·k·(I + η·A)) + μ, 1)
 */
export function foiSA(state, p) {
  return clamp(1 - Math.exp(-p.beta * p.k * (state.I + p.eta * state.A)) + p.mu, 0, 1);
}

/**
 * Build the 4×4 state-dependent transition matrix P(x_t) (Eq 4.1).
 * Rows/cols ordered [S, A, I, R]; every row sums to exactly 1.
 * Returns a 4×4 array plus the individual probabilities for the flow diagram.
 */
export function transitionMatrix(state, p) {
  const P_SA = foiSA(state, p);
  const P_AI = rateToProb(p.sigma);
  const P_AS = rateToProb(p.gamma1);
  const P_IR = rateToProb(p.gamma2);
  const P_RA = rateToProb(p.delta);

  // [from][to] with order S,A,I,R
  const M = [
    [1 - P_SA, P_SA, 0, 0],
    [P_AS, 1 - P_AS - P_AI, P_AI, 0],
    [0, 0, 1 - P_IR, P_IR],
    [0, P_RA, 0, 1 - P_RA],
  ];
  return { M, probs: { P_SA, P_AI, P_AS, P_IR, P_RA } };
}

/** Verify every row of a transition matrix sums to 1 (within tol). */
export function rowsSumToOne(M, tol = 1e-9) {
  return M.every((row) => Math.abs(row.reduce((a, b) => a + b, 0) - 1) < tol);
}

function renorm(x) {
  const s = x.S + x.A + x.I + x.R;
  return { S: x.S / s, A: x.A / s, I: x.I / s, R: x.R / s };
}

/**
 * Simulate the deterministic trajectory.
 * @param {object} p parameters
 * @param {object} init initial state {S,A,I,R}
 * @param {number} horizonYears projection length (dt = 1 month)
 * @returns {{month:number[],year:number[],S:number[],A:number[],I:number[],R:number[],incidence:number[]}}
 */
export function simulate(p, init, horizonYears = 20) {
  const steps = Math.round(horizonYears * 12);
  let x = renorm({ ...init });
  const out = {
    month: [0],
    year: [0],
    S: [x.S], A: [x.A], I: [x.I], R: [x.R],
    incidence: [0],
  };
  for (let t = 1; t <= steps; t++) {
    const { M } = transitionMatrix(x, p);
    const inc = x.S * M[0][1]; // new S→A entries this step
    // x_{t+1} = P^T x  (column update)
    const nx = {
      S: x.S * M[0][0] + x.A * M[1][0] + x.I * M[2][0] + x.R * M[3][0],
      A: x.S * M[0][1] + x.A * M[1][1] + x.I * M[2][1] + x.R * M[3][1],
      I: x.S * M[0][2] + x.A * M[1][2] + x.I * M[2][2] + x.R * M[3][2],
      R: x.S * M[0][3] + x.A * M[1][3] + x.I * M[2][3] + x.R * M[3][3],
    };
    x = renorm(nx);
    out.month.push(t);
    out.year.push(t / 12);
    out.S.push(x.S); out.A.push(x.A); out.I.push(x.I); out.R.push(x.R);
    out.incidence.push(inc);
  }
  return out;
}

/** Summarise peak and endemic behaviour of a trajectory. */
export function summarise(traj) {
  let peakI = -Infinity, peakIdx = 0;
  let peakInc = -Infinity, peakIncIdx = 0;
  for (let i = 0; i < traj.I.length; i++) {
    if (traj.I[i] > peakI) { peakI = traj.I[i]; peakIdx = i; }
    if (traj.incidence[i] > peakInc) { peakInc = traj.incidence[i]; peakIncIdx = i; }
  }
  const n = traj.I.length;
  return {
    peakPrevalence: peakI,
    peakYear: traj.year[peakIdx],
    endemicPrevalence: traj.I[n - 1],
    endemicState: { S: traj.S[n - 1], A: traj.A[n - 1], I: traj.I[n - 1], R: traj.R[n - 1] },
    peakIncidence: peakInc,
    peakIncidenceYear: traj.year[peakIncIdx],
  };
}
