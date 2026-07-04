// selftest-panel.js — Renders acceptance self-tests as PASS/FAIL rows.
// -----------------------------------------------------------------------------

import { h } from "../dom.js";
import { runSelfTests } from "../../model/selftest.js";
import { fx } from "../format.js";

export function selfTestPanel() {
  const { pass, results } = runSelfTests();
  const rows = results.map((r) => {
    const got = r.bool ? (r.pass ? "yes" : "no")
      : r.unit === "%" ? `${fx(r.got * 100, 1)}%`
      : r.unit === " yr" ? `${fx(r.got, 1)} yr`
      : fx(r.got, 2);
    const want = r.bool ? "true"
      : r.unit === "%" ? `${fx(r.want * 100, 0)}% ±${fx(r.tol * 100, 0)}`
      : r.unit === " yr" ? `${fx(r.want, 0)} ±${fx(r.tol, 0)} yr`
      : `${fx(r.want, 2)} ±${fx(r.tol, 2)}`;
    return h("tr", {}, [
      h("td", {}, [h("span", { class: `pill ${r.pass ? "pass" : "fail"}`, text: r.pass ? "PASS" : "FAIL" })]),
      h("td", { style: { textAlign: "left" }, text: r.name }),
      h("td", { text: got }),
      h("td", { class: "hint", text: want }),
      h("td", { class: "hint", style: { textAlign: "left" }, text: r.detail }),
    ]);
  });
  return h("div", {}, [
    h("div", { class: "row between", style: { marginBottom: "10px" } }, [
      h("div", {}, [
        h("strong", { text: "Acceptance self-tests " }),
        h("span", { class: `pill ${pass ? "pass" : "fail"}`, text: pass ? "ALL PASS" : "FAILURES" }),
      ]),
      h("span", { class: "hint", text: `${results.filter((r) => r.pass).length}/${results.length} checks` }),
    ]),
    h("table", { class: "data" }, [
      h("thead", {}, h("tr", {}, [h("th", { text: "" }), h("th", { style: { textAlign: "left" }, text: "Assertion" }), h("th", { text: "Model" }), h("th", { text: "Target" }), h("th", { style: { textAlign: "left" }, text: "Basis" })])),
      h("tbody", {}, rows),
    ]),
  ]);
}
