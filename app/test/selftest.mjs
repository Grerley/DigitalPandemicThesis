// selftest.mjs: headless check that the lab's JavaScript model reproduces the thesis's Python model.
//   node app/test/selftest.mjs        (exit 0 on success, 1 on any failure)

import { runSelfTests } from "../js/model/selftest.js";
import { makeGraph, simulateNetwork, betaFromB, Rnet } from "../js/model/network.js";
import { centralParams, R, LAB } from "../js/model/params.js";
import { simulate } from "../js/model/sair.js";
import { runNetworkJob } from "../js/model/network-job.js";

let failures = 0;
const p = (ok, name, detail = "") => {
  console.log(`${ok ? "  PASS" : "  FAIL"}  ${name}${detail ? "  ·  " + detail : ""}`);
  if (!ok) failures++;
};

console.log("\nDigital Pandemic Lab · self-test against the thesis model\n" + "-".repeat(60));
const t0 = Date.now();
const { results } = runSelfTests();
for (const r of results) p(r.pass, r.name, `lab ${r.got}, thesis ${r.want}${r.tol ? ` (tolerance ${r.tol})` : ""}`);

// Network: mean 2042 prevalence over a few runs should fall inside the thesis's 5th to 95th percentile range.
console.log("\n[network] stochastic runs compared with the thesis distribution");
const pc = centralParams();
const x22 = simulate(pc, LAB.central.x18, [4])[0];
const n = R("net_nodes"), k = R("net_mean_degree"), beta = betaFromB(pc, k);
for (const kind of ["ER", "WS", "BA"]) {
  const vals = [], rn = [];
  for (let s = 0; s < 4; s++) {
    const adj = makeGraph(kind, n, k, 100 + s);
    rn.push(Rnet(pc, beta, adj));
    const { frac } = simulateNetwork(adj, pc, beta, { months: 240, x0: x22, seed: 200 + s });
    vals.push(frac[frac.length - 1][2]);
  }
  const m = vals.reduce((a, b) => a + b) / vals.length, rm = rn.reduce((a, b) => a + b) / rn.length;
  const lo = R(`net_${kind}_I2042_p05`) / 100, hi = R(`net_${kind}_I2042_p95`) / 100;
  p(m >= lo && m <= hi, `${kind}: mean 2042 prevalence inside thesis 5th to 95th percentile range`,
    `lab ${(m * 100).toFixed(1)}%, thesis ${(lo * 100).toFixed(1)}% to ${(hi * 100).toFixed(1)}%`);
  const rt = R(`net_${kind}_R0net`);
  p(Math.abs(rm - rt) / rt < 0.1, `${kind}: network reproduction number within 10% of thesis`, `lab ${rm.toFixed(2)}, thesis ${rt.toFixed(2)}`);
}
// Targeting near the threshold: protecting the best-connected 10% beats random protection (thesis Chapter 8).
console.log("\n[network] targeting experiment near the threshold (Barabási–Albert)");
const job = (strategy) => runNetworkJob({ kind: "BA", n: 1000, k: 10, reps: 8, years: 20, params: pc, x22, regime: "near", strategy, cover: 0.1, damp: 0.2, seed: 4242 });
const rel = (r) => 1 - r.prot.mean.at(-1) / r.base.mean.at(-1);
const hub = rel(job("degree")), rnd = rel(job("random"));
p(Math.abs(job("none").rnet - 1.5) < 1e-6, "Near-threshold exposure solved so that R_net = 1.5");
p(hub > rnd + 0.15, "Degree-targeted protection beats random protection near the threshold",
  `lab ${(hub * 100).toFixed(0)}% against ${(rnd * 100).toFixed(0)}%; thesis ${R("netx_tgt_r1p5_hub_c10_rel_pct").toFixed(0)}% against ${R("netx_tgt_r1p5_random_c10_rel_pct").toFixed(0)}%`);

console.log(`\n${failures ? failures + " check(s) failed" : "All checks passed"} (${((Date.now() - t0) / 1000).toFixed(1)} s)\n`);
process.exit(failures ? 1 : 0);
