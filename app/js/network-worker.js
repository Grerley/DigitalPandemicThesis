// network-worker.js: runs network ensembles off the main thread.
import { runNetworkJob } from "./model/network-job.js";

self.onmessage = (ev) => {
  const { id, opts } = ev.data;
  try {
    const result = runNetworkJob(opts, (frac) => self.postMessage({ id, type: "progress", frac }));
    self.postMessage({ id, type: "done", result });
  } catch (err) {
    self.postMessage({ id, type: "error", message: String((err && err.message) || err) });
  }
};
