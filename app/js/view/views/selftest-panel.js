// selftest-panel.js: the self-test as PASS/FAIL rows.
import { h } from "../dom.js";
import { runSelfTests } from "../../model/selftest.js";

export function selfTestPanel() {
  const { pass, results } = runSelfTests();
  return h("div", {}, [
    h("div", { class: "row between", style: { marginBottom: "10px" } }, [
      h("div", {}, [h("strong", { text: "Lab model against the thesis model " }), h("span", { class: `pill ${pass ? "pass" : "fail"}`, text: pass ? "ALL PASS" : "FAILURES" })]),
      h("span", { class: "hint", text: `${results.filter((r) => r.pass).length}/${results.length} checks` }),
    ]),
    h("table", { class: "data" }, [
      h("thead", {}, h("tr", {}, ["", "Check", "Lab", "Thesis", "Tolerance", "Basis"].map((t, i) => h("th", { style: { textAlign: i === 1 || i === 5 ? "left" : "right" }, text: t })))),
      h("tbody", {}, results.map((r) => h("tr", {}, [
        h("td", {}, [h("span", { class: `pill ${r.pass ? "pass" : "fail"}`, text: r.pass ? "PASS" : "FAIL" })]),
        h("td", { style: { textAlign: "left" }, text: r.name }), h("td", { text: r.got }), h("td", { text: r.want }),
        h("td", { class: "hint", text: r.tol }), h("td", { class: "hint", style: { textAlign: "left" }, text: r.detail }),
      ]))),
    ]),
  ]);
}
