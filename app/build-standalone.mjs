// build-standalone.mjs — Produce a single, self-contained HTML build.
// -----------------------------------------------------------------------------
// Bundles the ES-module app into one inline <script> and inlines the CSS, so the
// whole lab is a single portable .html file that runs by double-clicking or from
// any static host (GitHub Pages, Netlify, an email attachment). Network Monte-
// Carlo ensembles run on the main thread in this build (the separate Web Worker
// only applies to the multi-file served version).
//
//   node app/build-standalone.mjs            → writes app/dist/digital-pandemic-lab.html
//
// Requires network access the first time to fetch esbuild via npx.

import { execFileSync } from "node:child_process";
import { readFileSync, writeFileSync, mkdirSync, rmSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const root = dirname(fileURLToPath(import.meta.url));
const dist = join(root, "dist");
mkdirSync(dist, { recursive: true });
const bundlePath = join(dist, "_bundle.js");

console.log("• bundling with esbuild…");
execFileSync("npx", ["-y", "esbuild@0.24.0", join(root, "js/app.js"),
  "--bundle", "--format=iife", `--outfile=${bundlePath}`], { stdio: "inherit" });

const css = readFileSync(join(root, "css/main.css"), "utf8");
const js = readFileSync(bundlePath, "utf8");
rmSync(bundlePath, { force: true });

const html = `<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<meta name="color-scheme" content="light dark" />
<title>The Digital Pandemic · Interactive Simulation Lab</title>
<meta name="description" content="Interactive lab for the S–A–I–R epidemic model of technology-induced mental health (Mutibura, 2026). All outputs are conditional model results." />
<style>
${css}
</style>
</head>
<body>
<script>window.__DPL_FORCE_MAIN_THREAD__ = true;</script>
<script>
${js}
</script>
</body>
</html>
`;

const out = join(dist, "digital-pandemic-lab.html");
writeFileSync(out, html);
console.log(`✓ wrote ${out}  (${(html.length / 1024).toFixed(0)} kB, self-contained)`);
