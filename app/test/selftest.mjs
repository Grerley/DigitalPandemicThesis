// selftest.mjs: headless check that the lab's JavaScript model reproduces the thesis's Python model.
//   node app/test/selftest.mjs        (exit 0 on success, 1 on any failure)

import { runSelfTests } from "../js/model/selftest.js";
import { makeGraph, simulateNetwork, betaFromB, Rnet } from "../js/model/network.js";
import { centralParams, R, LAB } from "../js/model/params.js";
import { simulate } from "../js/model/sair.js";

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
console.log(`\n${failures ? failures + " check(s) failed" : "All checks passed"} (${((Date.now() - t0) / 1000).toFixed(1)} s)\n`);
process.exit(failures ? 1 : 0);
