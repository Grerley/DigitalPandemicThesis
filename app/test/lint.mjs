// lint.mjs — Dependency-free syntax check for every source file.
// -----------------------------------------------------------------------------
// Parses each .js/.mjs file with `node --check`. Catches syntax errors without
// requiring any external linter or install step. Exits 1 if any file fails.
//
//   node app/test/lint.mjs

import { readdirSync, statSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath } from "node:url";
import { execFileSync } from "node:child_process";

const root = join(dirname(fileURLToPath(import.meta.url)), "..");

function walk(dir) {
  const out = [];
  for (const name of readdirSync(dir)) {
    if (name === "node_modules") continue;
    const full = join(dir, name);
    if (statSync(full).isDirectory()) out.push(...walk(full));
    else if (/\.(mjs|js)$/.test(name)) out.push(full);
  }
  return out;
}

const files = walk(root);
let failed = 0;
for (const f of files) {
  try {
    execFileSync(process.execPath, ["--check", f], { stdio: "pipe" });
    console.log(`  ✅ ${f.replace(root + "/", "")}`);
  } catch (e) {
    failed++;
    console.log(`  ‼️  ${f.replace(root + "/", "")}\n${e.stderr?.toString() || e.message}`);
  }
}
console.log(`\n${failed ? "‼️  " + failed + " file(s) failed" : "✅ " + files.length + " files OK"}`);
process.exit(failed ? 1 : 0);
