// sair.js: the S-A-I-R model of the thesis (Chapter 4), its threshold, equilibria and sensitivity (Chapter 5).
// -----------------------------------------------------------------------------
// A faithful port of python/dpm/model.py. All rates are per YEAR. State x = [S, A, I, R], fractions of the
// population aged 10 to 24.
//
//   dS/dt = ν − (λ + α) S + γ₁ A − ν S
//   dA/dt = (λ + α) S − (σ + γ₁ + ν) A + δ R
//   dI/dt = σ A − (γ₂ + ν) I
//   dR/dt = γ₂ I − (δ + ν) R,            λ = b (I + η A)
//
// Integration uses an adaptive Dormand–Prince 5(4) method (as scipy's RK45 in the Python code), integrated
// piecewise so that a step change in α (the temporary-shock scenario) falls on a step boundary.

export const STATES = ["S", "A", "I", "R"];
export const PARAM_KEYS = ["b", "eta", "alpha", "sigma", "gamma1", "gamma2", "delta", "nu"];

export function rhs(x, p, alpha = p.alpha) {
  const [S, A, I, R] = x;
  const lam = p.b * (I + p.eta * A);
  const inc = (lam + alpha) * S;
  return [
    p.nu - inc + p.gamma1 * A - p.nu * S,
    inc - (p.sigma + p.gamma1 + p.nu) * A + p.delta * R,
    p.sigma * A - (p.gamma2 + p.nu) * I,
    p.gamma2 * I - (p.delta + p.nu) * R,
  ];
}

// Dormand–Prince 5(4) tableau
const C = [0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1, 1];
const A = [
  [],
  [1 / 5],
  [3 / 40, 9 / 40],
  [44 / 45, -56 / 15, 32 / 9],
  [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729],
  [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656],
  [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
];
const B5 = [35 / 384, 0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0];
const B4 = [5179 / 57600, 0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40];

/** Integrate dx/dt = f(x) from t0 to t1 adaptively, returning x(t1). */
function dopri(f, x0, t0, t1, rtol = 1e-9, atol = 1e-11) {
  let t = t0, x = x0.slice(), hstep = Math.min(0.05, t1 - t0);
  const n = x.length;
  if (t1 <= t0) return x;
  let guard = 0;
  while (t < t1 - 1e-14) {
    if (++guard > 200000) throw new Error("integrator did not converge");
    if (t + hstep > t1) hstep = t1 - t;
    const k = [];
    for (let s = 0; s < 7; s++) {
      const xs = x.slice();
      for (let j = 0; j < s; j++) for (let i = 0; i < n; i++) xs[i] += hstep * A[s][j] * k[j][i];
      k.push(f(xs, t + C[s] * hstep));
    }
    const x5 = x.slice(), x4 = x.slice();
    for (let i = 0; i < n; i++) for (let s = 0; s < 7; s++) { x5[i] += hstep * B5[s] * k[s][i]; x4[i] += hstep * B4[s] * k[s][i]; }
    let err = 0;
    for (let i = 0; i < n; i++) {
      const sc = atol + rtol * Math.max(Math.abs(x[i]), Math.abs(x5[i]));
      err = Math.max(err, Math.abs(x5[i] - x4[i]) / sc);
    }
    if (err <= 1) { t += hstep; x = x5; }
    const fac = err === 0 ? 5 : Math.min(5, Math.max(0.2, 0.9 * Math.pow(err, -0.2)));
    hstep *= fac;
  }
  return x;
}

/**
 * Simulate the model and return the state at each time in tEval (years since start).
 * opts.shock = { dAlpha, start, end } adds dAlpha to α on [start, end) (the temporary-shock scenario).
 */
export function simulate(p, x0, tEval, opts = {}) {
  const sh = opts.shock;
  const brk = sh ? [sh.start, sh.end] : [];
  const alphaAt = (t) => (sh && t >= sh.start && t < sh.end ? p.alpha + sh.dAlpha : p.alpha);
  const out = [];
  let x = x0.slice(), t = 0;
  for (const te of tEval) {
    // advance from t to te, splitting at breakpoints
    const stops = brk.filter((b) => b > t && b < te).concat([te]);
    for (const s of stops) {
      const a = alphaAt((t + s) / 2);
      x = dopri((y) => rhs(y, p, a), x, t, s);
      t = s;
    }
    out.push(x.slice());
  }
  return out;
}

/** Prevalence of problematic use (I) at each year 0..years. */
export function prevalencePath(p, x0, years, opts = {}) {
  const t = Array.from({ length: years + 1 }, (_, i) => i);
  return simulate(p, x0, t, opts).map((x) => x[2]);
}

// ---------------------------------------------------------------- threshold
export function blocks(p) {
  const a = p.sigma + p.gamma1 + p.nu, c = p.gamma2 + p.nu, d = p.delta + p.nu;
  return { a, c, d, D: a * c * d - p.sigma * p.gamma2 * p.delta };
}

/** Basic reproduction number with α = 0: R₀ = b d (η c + σ) / D, D = a c d − σ γ₂ δ. */
export function R0(p) {
  const { c, d, D } = blocks(p);
  return (p.b * d * (p.eta * c + p.sigma)) / D;
}

/** Single-spell factor b(ηc + σ)/(ac) and relapse multiplier 1/(1 − q), q = (σ/a)(γ₂/c)(δ/d). */
export function R0Factors(p) {
  const { a, c, d } = blocks(p);
  const q = (p.sigma / a) * (p.gamma2 / c) * (p.delta / d);
  return { single: (p.b * (p.eta * c + p.sigma)) / (a * c), multiplier: 1 / (1 - q), q };
}

// ---------------------------------------------------------------- equilibria
/** The unique equilibrium in the simplex (Theorems 5.4 and 5.8). */
export function equilibrium(p) {
  const { c, d, D } = blocks(p);
  const m = 1 + p.sigma / c + (p.sigma * p.gamma2) / (c * d);
  const kappa = (p.b * (p.sigma + p.eta * c)) / c;
  const omega = D / (c * d);
  let Astar;
  if (p.alpha === 0) {
    const r0 = R0(p);
    if (r0 <= 1) return [1, 0, 0, 0];
    Astar = (1 - 1 / r0) / m;
  } else if (kappa === 0) {
    Astar = p.alpha / (p.alpha * m + omega);
  } else {
    const B = kappa - p.alpha * m - omega;
    Astar = (B + Math.sqrt(B * B + 4 * kappa * m * p.alpha)) / (2 * kappa * m);
  }
  const I = (p.sigma * Astar) / c;
  const R = (p.sigma * p.gamma2 * Astar) / (c * d);
  return [1 - Astar - I - R, Astar, I, R];
}

/** Floor of prevalence as social exposure falls to zero (equation 5.10). */
export function floor(p) {
  return equilibrium({ ...p, b: 0 })[2];
}

/** Local stability of an equilibrium on the simplex: Routh–Hurwitz on the 3×3 reduced Jacobian. */
export function isStable(p, x) {
  const [S, A, I] = x;
  const lam = p.b * (I + p.eta * A);
  const J = [
    [-(lam + p.alpha) - p.nu, -p.b * p.eta * S + p.gamma1, -p.b * S, 0],
    [lam + p.alpha, p.b * p.eta * S - (p.sigma + p.gamma1 + p.nu), p.b * S, p.delta],
    [0, p.sigma, -(p.gamma2 + p.nu), 0],
    [0, 0, p.gamma2, -(p.delta + p.nu)],
  ];
  // eliminate S = 1 − A − I − R: reduced matrix for (A, I, R)
  const T = [[-1, -1, -1], [1, 0, 0], [0, 1, 0], [0, 0, 1]];
  const Jr = [1, 2, 3].map((i) => [0, 1, 2].map((j) => J[i].reduce((acc, v, k) => acc + v * T[k][j], 0)));
  const tr = Jr[0][0] + Jr[1][1] + Jr[2][2];
  const m2 = Jr[0][0] * Jr[1][1] - Jr[0][1] * Jr[1][0] + Jr[0][0] * Jr[2][2] - Jr[0][2] * Jr[2][0] + Jr[1][1] * Jr[2][2] - Jr[1][2] * Jr[2][1];
  const det = Jr[0][0] * (Jr[1][1] * Jr[2][2] - Jr[1][2] * Jr[2][1]) - Jr[0][1] * (Jr[1][0] * Jr[2][2] - Jr[1][2] * Jr[2][0]) + Jr[0][2] * (Jr[1][0] * Jr[2][1] - Jr[1][1] * Jr[2][0]);
  const a1 = -tr, a2 = m2, a3 = -det; // characteristic polynomial λ³ + a1λ² + a2λ + a3
  return a1 > 0 && a3 > 0 && a1 * a2 > a3;
}

/** Normalised sensitivity indices ∂ln Q / ∂ln θ by central differences, for Q = R₀ (α = 0) or I*. */
export function elasticities(p, which = "I") {
  const f = which === "R0" ? (q) => R0({ ...q, alpha: 0 }) : (q) => equilibrium(q)[2];
  const keys = which === "R0" ? ["b", "eta", "sigma", "gamma1", "gamma2", "delta", "nu"] : PARAM_KEYS;
  const base = f(p), out = {};
  for (const k of keys) {
    if (p[k] === 0) { out[k] = 0; continue; }
    const hh = 1e-4;
    const up = f({ ...p, [k]: p[k] * (1 + hh) }), dn = f({ ...p, [k]: p[k] * (1 - hh) });
    out[k] = (up - dn) / (2 * hh * base);
  }
  return out;
}
