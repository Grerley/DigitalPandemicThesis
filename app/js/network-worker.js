// network-worker.js — Runs network Monte-Carlo ensembles off the main thread.
// -----------------------------------------------------------------------------
// Module worker. Receives a job spec, streams progress, posts the final
// ensemble (mean + 95% envelope). Keeps the UI responsive during heavy runs.

import { simulateEnsemble } from "./model/network.js";

self.onmessage = (ev) => {
  const { id, type, opts } = ev.data;
  if (type !== "ensemble") return;
  try {
    const result = simulateEnsemble(opts.topology, {
      n: opts.n, meanDegree: opts.meanDegree, runs: opts.runs,
      horizonMonths: opts.horizonMonths, seed: opts.seed, params: opts.params,
      seedHub: opts.seedHub, immunizeFrac: opts.immunizeFrac, immunizeStrategy: opts.immunizeStrategy,
      seedFrac: opts.seedFrac,
    }, (frac) => self.postMessage({ id, type: "progress", frac }));
    self.postMessage({ id, type: "done", result });
  } catch (err) {
    self.postMessage({ id, type: "error", message: String(err && err.message || err) });
  }
};
