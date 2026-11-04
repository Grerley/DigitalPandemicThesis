// network.js: stochastic S-A-I-R dynamics on explicit contact networks (thesis Chapter 8; python/dpm/network.py).
// -----------------------------------------------------------------------------
// Discrete monthly steps, one node per young person. Each month:
//   S → A with probability 1 − (1 − α_m)(1 − β)^(w_i (m_I + η m_A)), m_I and m_A the Addicted and At-Risk neighbours,
//          w_i the node's susceptibility to social exposure (1, or a reduced value if protected);
//   A, I and R leave by competing hazards converted with P = 1 − exp(−h Δt);
//   turnover: with probability 1 − exp(−ν Δt) a node ages out and is replaced by a new Susceptible.
// Graphs are generated in JavaScript with the same families and settings as the thesis; they are not the
// same random graphs, so individual runs differ from the thesis runs while their distributions agree.

import { RNG } from "./rng.js";

export const DT = 1 / 12;

// ---------------------------------------------------------------- graphs (adjacency lists)
function edgesToAdj(n, edges) {
  const adj = Array.from({ length: n }, () => []);
  for (const [u, v] of edges) { adj[u].push(v); adj[v].push(u); }
  return adj;
}

/** Erdős–Rényi G(n, p) with p = k/(n − 1), by geometric skipping. */
function erdosRenyi(n, k, rng) {
  const p = k / (n - 1), edges = [];
  const lp = Math.log(1 - p);
  let v = 1, w = -1;
  while (v < n) {
    w += 1 + Math.floor(Math.log(1 - rng.next()) / lp);
    while (w >= v && v < n) { w -= v; v++; }
    if (v < n) edges.push([v, w]);
  }
  return edgesToAdj(n, edges);
}

/** Watts–Strogatz: ring lattice with k neighbours, each edge rewired with probability pr. */
function wattsStrogatz(n, k, rng, pr = 0.1) {
  const half = Math.round(k) / 2, set = new Set(), key = (u, v) => (u < v ? u * n + v : v * n + u);
  for (let u = 0; u < n; u++) for (let j = 1; j <= half; j++) set.add(key(u, (u + j) % n));
  for (let j = 1; j <= half; j++) {
    for (let u = 0; u < n; u++) {
      const v = (u + j) % n;
      if (rng.next() < pr) {
        let w = rng.int(n), guard = 0;
        while ((w === u || set.has(key(u, w))) && guard++ < 50) w = rng.int(n);
        if (guard < 50) { set.delete(key(u, v)); set.add(key(u, w)); }
      }
    }
  }
  return edgesToAdj(n, [...set].map((e) => [Math.floor(e / n), e % n]));
}

/** Barabási–Albert preferential attachment with m = k/2 edges per new node, started from a star on m + 1 nodes. */
function barabasiAlbert(n, k, rng) {
  const m = Math.round(k / 2), edges = [], repeated = [];
  for (let i = 1; i <= m; i++) { edges.push([0, i]); repeated.push(0, i); }
  for (let v = m + 1; v < n; v++) {
    const targets = new Set();
    while (targets.size < m) targets.add(repeated[rng.int(repeated.length)]);
    for (const t of targets) { edges.push([v, t]); repeated.push(v, t); }
  }
  return edgesToAdj(n, edges);
}

export const GRAPH_LABELS = { ER: "Erdős–Rényi", WS: "Watts–Strogatz", BA: "Barabási–Albert" };

export function makeGraph(kind, n, k, seed) {
  const rng = new RNG(seed);
  if (kind === "ER") return erdosRenyi(n, k, rng);
  if (kind === "WS") return wattsStrogatz(n, k, rng);
  if (kind === "BA") return barabasiAlbert(n, k, rng);
  throw new Error(kind);
}

export function degreeMoments(adj) {
  let k1 = 0, k2 = 0;
  for (const nb of adj) { k1 += nb.length; k2 += nb.length * nb.length; }
  return { k1: k1 / adj.length, k2: k2 / adj.length };
}

// ---------------------------------------------------------------- transmissibility and R_net
/** Per-edge monthly transmission probability matching the mean-field exposure rate b = rate × k. */
export const betaFromB = (p, k) => 1 - Math.exp(-p.b / k / 12);

function rates(p, beta) {
  return { theta: -Math.log(1 - beta) / DT, a: p.sigma + p.gamma1 + p.nu, c: p.gamma2 + p.nu, d: p.delta + p.nu };
}

/** Career transmissibility: probability of exposing a given neighbour over all spells, relapses included. */
export function careerTransmissibility(p, beta) {
  const { theta, a, c, d } = rates(p, beta);
  const hA = a + p.eta * theta, hI = c + theta;
  const X = (p.gamma1 + p.nu) / hA + (p.sigma / hA) * (p.nu / hI + (p.gamma2 / hI) * (p.nu / d));
  const Y = (p.sigma / hA) * (p.gamma2 / hI) * (p.delta / d);
  return 1 - X / (1 - Y);
}

/** Heuristic network reproduction number T_career (⟨k²⟩ − ⟨k⟩)/⟨k⟩ (α set to zero). */
export function Rnet(p, beta, adj) {
  const { k1, k2 } = degreeMoments(adj);
  return careerTransmissibility({ ...p, alpha: 0 }, beta) * (k2 - k1) / k1;
}

// ---------------------------------------------------------------- simulation
/**
 * Run the monthly process. Returns { frac: [[S,A,I,R] per recorded month], states } .
 * opts: { months, x0, seed, recordEvery = 12, weights (Float64Array of susceptibility multipliers) }
 */
export function simulateNetwork(adj, p, beta, opts) {
  const { months, x0, seed, recordEvery = 12, weights = null, keepStates = false } = opts;
  const rng = new RNG(seed), n = adj.length;
  let st = new Uint8Array(n);
  const tot = x0[0] + x0[1] + x0[2] + x0[3];
  for (let i = 0; i < n; i++) {
    let r = rng.next() * tot, s = 0;
    while (s < 3 && r >= x0[s]) { r -= x0[s]; s++; }
    st[i] = s;
  }
  const am = 1 - Math.exp(-p.alpha * DT);
  const pA = 1 - Math.exp(-(p.sigma + p.gamma1) * DT), shareI = p.sigma / (p.sigma + p.gamma1);
  const pI = 1 - Math.exp(-p.gamma2 * DT), pR = 1 - Math.exp(-p.delta * DT), pexit = 1 - Math.exp(-p.nu * DT);
  const l1mb = Math.log(1 - beta);
  const frac = [], snapshots = [];
  const record = () => {
    const c = [0, 0, 0, 0];
    for (let i = 0; i < n; i++) c[st[i]]++;
    frac.push(c.map((v) => v / n));
    if (keepStates) snapshots.push(st.slice());
  };
  const next = new Uint8Array(n);
  for (let t = 0; t <= months; t++) {
    if (t % recordEvery === 0) record();
    if (t === months) break;
    for (let i = 0; i < n; i++) {
      const s = st[i], u = rng.next();
      let ns = s;
      if (s === 0) {
        let e = 0;
        for (const j of adj[i]) { const sj = st[j]; if (sj === 2) e += 1; else if (sj === 1) e += p.eta; }
        if (weights) e *= weights[i];
        if (u < 1 - (1 - am) * Math.exp(l1mb * e)) ns = 1;
      } else if (s === 1) {
        if (u < pA) ns = u < pA * shareI ? 2 : 0;
      } else if (s === 2) {
        if (u < pI) ns = 3;
      } else if (u < pR) ns = 1;
      next[i] = ns;
    }
    for (let i = 0; i < n; i++) st[i] = rng.next() < pexit ? 0 : next[i];
  }
  return { frac, snapshots };
}

/** Indices of the m highest-degree nodes. */
export function topDegree(adj, m) {
  return adj.map((nb, i) => [nb.length, i]).sort((a, b) => b[0] - a[0]).slice(0, m).map((x) => x[1]);
}

/** Node-state codes for drawing (V marks a protected Susceptible). */
export const STATE = { S: 0, A: 1, I: 2, R: 3, V: 4 };
