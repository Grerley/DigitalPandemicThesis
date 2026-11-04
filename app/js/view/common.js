// common.js: small helpers shared by the screens.
import { h } from "./dom.js";
import { LAB } from "../model/params.js";

export function viewHead(kicker, title, intro) {
  return h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: kicker }),
    h("h1", { text: title }),
    intro ? h("p", { class: "prose", html: intro }) : null,
  ]);
}

/** Simple data table. rows: array of arrays of cell contents; align: "l" or "r" per column. */
export function table(head, rows, align = []) {
  return h("div", { class: "tablewrap" }, h("table", { class: "data" }, [
    h("thead", {}, h("tr", {}, head.map((c, i) => h("th", { style: { textAlign: align[i] === "l" ? "left" : "right" }, html: c })))),
    h("tbody", {}, rows.map((r) => h("tr", {}, r.map((c, i) => h("td", { style: { textAlign: align[i] === "l" ? "left" : "right" }, html: String(c) }))))),
  ]));
}

/** Short source note under a panel. */
export const source = (html) => h("p", { class: "hint", style: { marginTop: "10px" }, html: `Source: ${html}` });

export const pc = (x, dp = 1) => `${(x * 100).toFixed(dp)}%`;
export const COLORS = { P: "var(--c-i)", T: "var(--c-s)", S: "var(--c-s)", A: "var(--c-a)", I: "var(--c-i)", R: "var(--c-r)" };
export const HBSC = LAB.hbsc;
