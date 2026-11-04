// ensemble-runner.js: run a network job in a Web Worker, or on the main thread when workers are unavailable.
import { runNetworkJob } from "../model/network-job.js";

let worker = null, workerOK = null, jobId = 0;
const pending = new Map();

function ensureWorker() {
  if (workerOK !== null) return workerOK;
  if (globalThis.__DPL_FORCE_MAIN_THREAD__) return (workerOK = false);
  try {
    worker = new Worker(new URL("../network-worker.js", import.meta.url), { type: "module" });
    worker.onmessage = (ev) => {
      const { id, type, frac, result, message } = ev.data, job = pending.get(id);
      if (!job) return;
      if (type === "progress") job.onProgress?.(frac);
      else if (type === "done") { job.resolve(result); pending.delete(id); }
      else { job.reject(new Error(message)); pending.delete(id); }
    };
    worker.onerror = () => { workerOK = false; };
    workerOK = true;
  } catch { workerOK = false; }
  return workerOK;
}

export function runJob(opts, onProgress) {
  if (ensureWorker()) {
    const id = ++jobId;
    return new Promise((resolve, reject) => { pending.set(id, { resolve, reject, onProgress }); worker.postMessage({ id, opts }); });
  }
  return new Promise((resolve) => setTimeout(() => resolve(runNetworkJob(opts, onProgress)), 20));
}
