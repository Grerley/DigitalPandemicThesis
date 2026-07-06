// selftest.mjs — Headless acceptance self-test for the Digital Pandemic Lab.
// -----------------------------------------------------------------------------
// Runs the model layer under Node (no browser) and asserts the thesis
// acceptance targets plus the qualitative network/cohort claims. Exits 0 on
// success, 1 on any failure — suitable for CI and the SessionStart hook.
//
//   node app/test/selftest.mjs

import { runSelfTests } from "../js/model/selftest.js";
import { simulateEnsemble } from "../js/model/network.js";
import { generateCohort } from "../js/model/cohort.js";
import { baselineParams } from "../js/model/params.js";

let failures = 0;
const p = (ok, name, detail = "") => {
  console.log(`${ok ? "  ✅ PASS" : "  ‼️  FAIL"}  ${name}${detail ? "  —  " + detail : ""}`);
  if (!ok) failures++;
};

console.log("\nDigital Pandemic Lab · acceptance self-test\n" + "-".repeat(52));

// 1) Deterministic acceptance targets (R₀ trio, peak, endemic, comprehensive, invariants)
console.log("\n[1] Deterministic core & acceptance targets");
const { pass, results } = runSelfTests();
for (const r of results) {
  const got = r.bool ? (r.pass ? "yes" : "no")
    : r.unit === "%" ? `${(r.got * 100).toFixed(1)}%`
    : r.unit === " yr" ? `${r.got.toFixed(1)} yr` : r.got.toFixed(2);
  p(r.pass, r.name, `model ${got}`);
}

// 2) Network: scale-free peaks higher AND earlier than random; hubs beat random.
console.log("\n[2] Network Monte-Carlo (reduced runs)");
const base = baselineParams();
const netP = { ...base, beta: (base.beta * base.k) / 10, mu: 0 };
const opt = { n: 400, meanDegree: 10, runs: 20, horizonMonths: 180, params: netP, seed: 42, seedFrac: 0.004 };
const ER = simulateEnsemble("ER", opt);
const WS = simulateEnsemble("WS", opt);
const BA = simulateEnsemble("BA", opt);
p(BA.meanPeak > ER.meanPeak && BA.meanPeak > WS.meanPeak,
  "Scale-free peaks higher than random & small-world",
  `BA ${(BA.meanPeak * 100).toFixed(1)}% vs ER ${(ER.meanPeak * 100).toFixed(1)}% / WS ${(WS.meanPeak * 100).toFixed(1)}%`);
p(BA.meanPeakYear < ER.meanPeakYear && BA.meanPeakYear < WS.meanPeakYear,
  "Scale-free peaks earlier than random & small-world",
  `BA ${BA.meanPeakYear.toFixed(1)}yr vs ER ${ER.meanPeakYear.toFixed(1)}yr / WS ${WS.meanPeakYear.toFixed(1)}yr`);

const immBase = { ...opt, topology: "BA", seed: 11 };
const none = simulateEnsemble("BA", { ...opt, seed: 11 });
const hub = simulateEnsemble("BA", { ...opt, seed: 11, immunizeFrac: 0.3, immunizeStrategy: "hub" });
const rnd = simulateEnsemble("BA", { ...opt, seed: 11, immunizeFrac: 0.3, immunizeStrategy: "random" });
p(hub.meanPeak < rnd.meanPeak,
  "Hub-targeted immunisation beats random at equal coverage",
  `hub ${(hub.meanPeak * 100).toFixed(1)}% < random ${(rnd.meanPeak * 100).toFixed(1)}% (from ${(none.meanPeak * 100).toFixed(1)}%)`);

// 3) Cohort: transition-matrix recovery is accurate.
console.log("\n[3] Synthetic cohort & statistics");
const cohort = generateCohort(base, { n: 3000, years: 10, seed: 42 });
p(cohort.matrixMAE < 0.05, "Estimated transition matrix recovers the known matrix",
  `MAE ${cohort.matrixMAE.toFixed(3)} < 0.05`);
p(cohort.oddsRatios.every((o) => (o.knownOR > 1) === (o.estOR > 1)),
  "Recovered odds ratios match the direction of published ORs");

// Determinism: same seed ⇒ identical result.
const c2 = generateCohort(base, { n: 3000, years: 10, seed: 42 });
p(Math.abs(c2.matrixMAE - cohort.matrixMAE) < 1e-12, "Seeded RNG is reproducible (same seed ⇒ same result)");

console.log("\n" + "-".repeat(52));
if (failures === 0) {
  console.log(`✅ ALL CHECKS PASSED (${results.length + 6} assertions)\n`);
  process.exit(0);
} else {
  console.log(`‼️  ${failures} CHECK(S) FAILED\n`);
  process.exit(1);
}
