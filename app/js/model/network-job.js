// network-job.js: one network experiment, as run on the Networks screen (worker or main thread).
// For each replicate the same graph and random seed are used without protection and with the chosen
// protection, so the difference measures the protection itself (paired design, as in the thesis). Near the
// threshold, the exposure rate is solved separately for each graph so that its R_net is 1.5; the thesis solved it
// once on one graph. Watts–Strogatz graphs are not forced to be connected, unlike networkx's connected variant.

import { makeGraph, simulateNetwork, betaFromB, Rnet, topDegree } from "./network.js";
import { calibrateRates } from "./calibrate.js";
import { brent } from "./calibrate.js";
import { RNG } from "./rng.js";
import { simulate, equilibrium } from "./sair.js";

/**
 * opts: { kind, n, k, reps, years, params (b = post-2018 exposure), x22, regime: "calibrated"|"near",
 *         strategy: "none"|"random"|"degree"|"acquaintance", cover, damp, seed }
 */
export function runNetworkJob(opts, onProgress = () => {}) {
  const { kind, n, k, reps, years, regime, strategy, cover, damp, seed } = opts;
  let p = { ...opts.params };
  const months = years * 12;
  const base = [], prot = [], rnets = [];
  const bs = [];
  for (let r = 0; r < reps; r++) {
    const adj = makeGraph(kind, n, k, seed + r);
    let q = p, x0 = opts.x22;
    if (regime === "near") {
      // exposure rate at which the heuristic R_net of this graph equals 1.5
      const b = brent((bb) => Rnet({ ...p, b: bb }, betaFromB({ ...p, b: bb }, k), adj) - 1.5, 1e-4, 200, 1e-8);
      q = { ...p, b }; x0 = equilibrium(q);
    }
    bs.push(q.b);
    const beta = betaFromB(q, k);
    rnets.push(Rnet(q, beta, adj));
    const run = (w) => simulateNetwork(adj, q, beta, { months, x0, seed: seed + 1000 + r, recordEvery: 12, weights: w }).frac.map((f) => f[2]);
    base.push(run(null));
    if (strategy !== "none") {
      const m = Math.round(cover * n), w = new Float64Array(n).fill(1), rng = new RNG(seed + 2000 + r);
      let chosen;
      if (strategy === "random") chosen = rng.sampleIndices(n, m);
      else if (strategy === "degree") chosen = topDegree(adj, m);
      else { // acquaintance: a random neighbour of a randomly chosen young person, until m distinct are protected
        const set = new Set(); let guard = 0;
        while (set.size < m && guard++ < 50 * m) { const i = rng.int(n); if (adj[i].length) set.add(adj[i][rng.int(adj[i].length)]); }
        chosen = [...set];
      }
      for (const i of chosen) w[i] = damp;
      prot.push(run(w));
    }
    onProgress((r + 1) / reps);
  }
  const summarise = (runs) => {
    if (!runs.length) return null;
    const T = runs[0].length, mean = [], lo = [], hi = [];
    for (let t = 0; t < T; t++) {
      const v = runs.map((x) => x[t]).sort((a, b) => a - b);
      mean.push(v.reduce((a, b) => a + b, 0) / v.length);
      lo.push(v[Math.floor(0.05 * (v.length - 1))]); hi.push(v[Math.ceil(0.95 * (v.length - 1))]);
    }
    return { mean, lo, hi, final: runs.map((x) => x[x.length - 1]) };
  };
  return { base: summarise(base), prot: summarise(prot), rnet: rnets.reduce((a, b) => a + b, 0) / rnets.length, bUsed: bs.reduce((a, b) => a + b, 0) / bs.length };
}
