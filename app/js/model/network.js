// network.js — Stochastic Reed–Frost contagion on explicit graphs.
// -----------------------------------------------------------------------------
// Port of R/02_network_simulation.R with added node-immunisation. Three
// topologies (Erdős–Rényi, Watts–Strogatz, Barabási–Albert) are generated from
// a seeded RNG so ensembles are exactly reproducible. States are encoded as
// integers: 0=S, 1=A(At-Risk), 2=I(Addicted), 3=R(Recovered), 4=V(protected).
//
// PURE module: no DOM. Imported by the web worker and by a main-thread fallback.

import { RNG } from "./rng.js";
import { rateToProb } from "./sair.js";

export const STATE = { S: 0, A: 1, I: 2, R: 3, V: 4 };

// ---- Graph generators -------------------------------------------------------

/** Erdős–Rényi G(n, p) with p = meanDegree/(n−1). */
function erdosRenyi(n, meanDegree, rng) {
  const adj = Array.from({ length: n }, () => []);
  const p = meanDegree / (n - 1);
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      if (rng.next() < p) { adj[i].push(j); adj[j].push(i); }
    }
  }
  return adj;
}

/** Watts–Strogatz small-world: ring lattice (nei each side) rewired w.p. β. */
function wattsStrogatz(n, meanDegree, rng, beta = 0.1) {
  const nei = Math.max(1, Math.floor(meanDegree / 2));
  const set = Array.from({ length: n }, () => new Set());
  for (let i = 0; i < n; i++) {
    for (let d = 1; d <= nei; d++) {
      const j = (i + d) % n;
      set[i].add(j); set[j].add(i);
    }
  }
  // Rewire each forward edge with probability beta.
  for (let i = 0; i < n; i++) {
    for (let d = 1; d <= nei; d++) {
      const j = (i + d) % n;
      if (rng.next() < beta) {
        set[i].delete(j); set[j].delete(i);
        let t = rng.int(n), guard = 0;
        while ((t === i || set[i].has(t)) && guard++ < 50) t = rng.int(n);
        set[i].add(t); set[t].add(i);
      }
    }
  }
  return set.map((s) => [...s]);
}

/** Barabási–Albert scale-free by preferential attachment (m new edges/node). */
function barabasiAlbert(n, meanDegree, rng) {
  const m = Math.max(1, Math.floor(meanDegree / 2));
  const adj = Array.from({ length: n }, () => []);
  const targets = []; // multiset of node ids, degree-weighted
  // Seed: m0 = m fully connected nodes.
  const m0 = m + 1;
  for (let i = 0; i < m0; i++) {
    for (let j = i + 1; j < m0; j++) { adj[i].push(j); adj[j].push(i); targets.push(i, j); }
  }
  for (let v = m0; v < n; v++) {
    const chosen = new Set();
    let guard = 0;
    while (chosen.size < m && guard++ < 200) {
      const t = targets.length ? targets[rng.int(targets.length)] : rng.int(v);
      if (t !== v) chosen.add(t);
    }
    for (const t of chosen) {
      adj[v].push(t); adj[t].push(v);
      targets.push(v, t);
    }
  }
  return adj;
}

/**
 * Build a network of the requested topology.
 * @param {"ER"|"WS"|"BA"} type
 * @returns {{adj:number[][], degree:number[]}}
 */
export function makeNetwork(type, n, meanDegree, rng) {
  let adj;
  if (type === "ER") adj = erdosRenyi(n, meanDegree, rng);
  else if (type === "WS") adj = wattsStrogatz(n, meanDegree, rng);
  else adj = barabasiAlbert(n, meanDegree, rng);
  const degree = adj.map((a) => a.length);
  return { adj, degree };
}

// ---- Contagion --------------------------------------------------------------

/**
 * One stochastic Reed–Frost run on a fixed graph.
 * @param {{adj:number[][],degree:number[]}} g
 * @param {object} p parameters
 * @param {object} opts { horizonMonths, seedFrac, seedHub, immunizeFrac,
 *                        immunizeStrategy ('hub'|'random'), recordStates }
 * @param {RNG} rng
 * @returns trajectory arrays; if recordStates, also `frames` (per-step int8 states)
 */
export function simulateNetworkOnce(g, p, opts, rng) {
  const {
    horizonMonths = 240, seedFrac = 0.01, seedHub = false,
    immunizeFrac = 0, immunizeStrategy = "hub", recordStates = false, frameStride = 1,
  } = opts;
  const n = g.adj.length;
  const state = new Int8Array(n); // all S = 0

  // Immunise a coverage fraction (protected forever).
  const nImm = Math.round(immunizeFrac * n);
  if (nImm > 0) {
    let idx;
    if (immunizeStrategy === "hub") {
      idx = [...g.degree.keys()].sort((a, b) => g.degree[b] - g.degree[a]);
    } else {
      idx = rng.shuffle([...Array(n).keys()]);
    }
    for (let i = 0; i < nImm; i++) state[idx[i]] = STATE.V;
  }

  // Seed Addicted among non-immune nodes.
  const nSeed = Math.max(1, Math.round(seedFrac * n));
  let seedPool;
  if (seedHub) seedPool = [...g.degree.keys()].sort((a, b) => g.degree[b] - g.degree[a]);
  else seedPool = rng.shuffle([...Array(n).keys()]);
  let placed = 0;
  for (const node of seedPool) {
    if (placed >= nSeed) break;
    if (state[node] === STATE.S) { state[node] = STATE.I; placed++; }
  }

  const pAI = rateToProb(p.sigma), pAS = rateToProb(p.gamma1);
  const pIR = rateToProb(p.gamma2), pRA = rateToProb(p.delta);

  const T = horizonMonths;
  const S = new Float64Array(T + 1), A = new Float64Array(T + 1),
        I = new Float64Array(T + 1), R = new Float64Array(T + 1);
  const frames = recordStates ? [] : null;
  const tally = (t) => {
    let s = 0, a = 0, i = 0, r = 0;
    for (let k = 0; k < n; k++) {
      const st = state[k];
      if (st === STATE.S) s++; else if (st === STATE.A) a++;
      else if (st === STATE.I) i++; else if (st === STATE.R) r++;
    }
    S[t] = s / n; A[t] = a / n; I[t] = i / n; R[t] = r / n;
  };
  tally(0);
  if (recordStates) frames.push(Int8Array.from(state));

  const next = new Int8Array(n);
  for (let t = 1; t <= T; t++) {
    next.set(state);
    for (let i = 0; i < n; i++) {
      const s = state[i];
      if (s === STATE.S) {
        const nb = g.adj[i];
        let mInf = 0;
        for (let q = 0; q < nb.length; q++) {
          const ns = state[nb[q]];
          if (ns === STATE.I) mInf += 1;
          else if (ns === STATE.A) mInf += p.eta;
        }
        const pinf = Math.min(1 - Math.pow(1 - p.beta, mInf) + p.mu, 1);
        if (rng.next() < pinf) next[i] = STATE.A;
      } else if (s === STATE.A) {
        const u = rng.next();
        if (u < pAI) next[i] = STATE.I;
        else if (u < pAI + pAS) next[i] = STATE.S;
      } else if (s === STATE.I) {
        if (rng.next() < pIR) next[i] = STATE.R;
      } else if (s === STATE.R) {
        if (rng.next() < pRA) next[i] = STATE.A;
      }
      // V stays V.
    }
    state.set(next);
    tally(t);
    if (recordStates && t % frameStride === 0) frames.push(Int8Array.from(state));
  }
  return { S, A, I, R, frames, months: T };
}

/**
 * Monte-Carlo ensemble: many graph draws / runs. Returns mean trajectory and a
 * 95% envelope (2.5–97.5 percentiles) for the Addicted compartment.
 */
export function simulateEnsemble(type, opts, onProgress) {
  const {
    n = 800, meanDegree = 10, runs = 60, horizonMonths = 240, seed = 1,
    params, seedHub = false, immunizeFrac = 0, immunizeStrategy = "hub", seedFrac = 0.01,
  } = opts;
  const master = new RNG(seed);
  const T = horizonMonths;
  const all = []; // runs × (T+1) Addicted fraction
  let peakSum = 0, peakYearSum = 0;
  for (let r = 0; r < runs; r++) {
    const rng = new RNG((master.int(2 ** 30) ^ (r * 2654435761)) >>> 0);
    const g = makeNetwork(type, n, meanDegree, rng);
    const run = simulateNetworkOnce(g, params, {
      horizonMonths, seedHub, immunizeFrac, immunizeStrategy, seedFrac,
    }, rng);
    all.push(run.I);
    let pk = 0, pkT = 0;
    for (let t = 0; t <= T; t++) if (run.I[t] > pk) { pk = run.I[t]; pkT = t; }
    peakSum += pk; peakYearSum += pkT / 12;
    if (onProgress && r % 5 === 0) onProgress((r + 1) / runs);
  }
  const months = [], year = [], mean = [], lo = [], hi = [];
  const col = new Float64Array(runs);
  for (let t = 0; t <= T; t++) {
    let sum = 0;
    for (let r = 0; r < runs; r++) { col[r] = all[r][t]; sum += all[r][t]; }
    const sorted = Float64Array.from(col).sort();
    months.push(t); year.push(t / 12);
    mean.push(sum / runs);
    lo.push(percentile(sorted, 0.025));
    hi.push(percentile(sorted, 0.975));
  }
  if (onProgress) onProgress(1);
  return {
    type, months, year, mean, lo, hi, runs,
    meanPeak: peakSum / runs, meanPeakYear: peakYearSum / runs,
  };
}

function percentile(sorted, q) {
  const idx = q * (sorted.length - 1);
  const lo = Math.floor(idx), hi = Math.ceil(idx);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (sorted[hi] - sorted[lo]) * (idx - lo);
}
