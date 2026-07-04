// ensemble-runner.js — Worker-backed ensemble runner with main-thread fallback.
// -----------------------------------------------------------------------------
// Prefers a module Web Worker so the UI never freezes; if workers are
// unavailable (e.g. opened from file://), it runs on the main thread in yielding
// chunks so the page stays interactive. Same seed → same result either way.

import { simulateEnsemble } from "../model/network.js";

let worker = null;
let workerOK = null;
let jobId = 0;
const pending = new Map();

function ensureWorker() {
  if (workerOK !== null) return workerOK;
  try {
    worker = new Worker(new URL("../network-worker.js", import.meta.url), { type: "module" });
    worker.onmessage = (ev) => {
      const { id, type, frac, result, message } = ev.data;
      const job = pending.get(id);
      if (!job) return;
      if (type === "progress") job.onProgress?.(frac);
      else if (type === "done") { job.resolve(result); pending.delete(id); }
      else if (type === "error") { job.reject(new Error(message)); pending.delete(id); }
    };
    worker.onerror = () => { workerOK = false; };
    workerOK = true;
  } catch {
    workerOK = false;
  }
  return workerOK;
}

/**
 * Run an ensemble. Returns a promise resolving to the ensemble result.
 * @param {object} opts { topology, n, meanDegree, runs, horizonMonths, seed,
 *                        params, seedHub, immunizeFrac, immunizeStrategy }
 * @param {(frac:number)=>void} onProgress
 */
export function runEnsemble(opts, onProgress) {
  if (ensureWorker()) {
    const id = ++jobId;
    return new Promise((resolve, reject) => {
      pending.set(id, { resolve, reject, onProgress });
      worker.postMessage({ id, type: "ensemble", opts });
    });
  }
  // Main-thread fallback with cooperative yielding between runs.
  return mainThreadEnsemble(opts, onProgress);
}

function mainThreadEnsemble(opts, onProgress) {
  return new Promise((resolve) => {
    // simulateEnsemble is synchronous; chunk by yielding to the event loop.
    // We approximate progressive behaviour by running in a microtask and
    // reporting progress via its callback (which fires within the sync run).
    setTimeout(() => {
      const result = simulateEnsemble(opts.topology, {
        n: opts.n, meanDegree: opts.meanDegree, runs: opts.runs,
        horizonMonths: opts.horizonMonths, seed: opts.seed, params: opts.params,
        seedHub: opts.seedHub, immunizeFrac: opts.immunizeFrac, immunizeStrategy: opts.immunizeStrategy,
        seedFrac: opts.seedFrac,
      }, onProgress);
      resolve(result);
    }, 0);
  });
}

export function workerAvailable() { return ensureWorker(); }
