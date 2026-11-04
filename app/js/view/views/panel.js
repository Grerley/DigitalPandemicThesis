// panel.js: Screen 7, designing the panel study that could estimate the model's rates (thesis Chapter 6).

import { h } from "../dom.js";
import { panel, chartCard } from "../ui.js";
import { xyLines } from "../charts.js";
import { viewHead, table, source } from "../common.js";
import { LAB, R } from "../../model/params.js";

const TR = ["S→A", "A→S", "A→I", "I→R", "R→A"];
const NAME = { "S→A": "Incidence (S→A)", "A→S": "Remission (A→S)", "A→I": "Progression (A→I)", "I→R": "Recovery (I→R)", "R→A": "Relapse (R→A)" };
const COL = ["var(--c-s)", "var(--accent)", "var(--c-a)", "var(--c-i)", "var(--c-r)"];

const fmtPct = (v) => { const r = Math.round(100 * v); return `${r === 0 ? 0 : r}%`; };

export function PanelView() {
  const el = h("div", { class: "view" });
  el.appendChild(viewHead("Screen 7 · Panel design", "What data would identify the rates?",
    "Two of the four states are never observed directly: At-Risk lies below the survey cut-off, and Recovered looks the same as Susceptible in a single survey. " +
    "The thesis simulated panels of young people with known rates, observed each year through a noisy symptom count, and re-estimated the annual transition probabilities with a hidden Markov model fitted by the EM algorithm (Zucchini et al., 2016). " +
    "The generating rates were chosen for the design study and are not estimates for any population."));

  let noise = "low";
  const rows = (nz) => LAB.recovery.filter((r) => r.noise === nz);
  const chart = chartCard({ title: "Relative error of the recovered annual transition probabilities", sub: `${R("rec_waves")} annual waves (thesis results)`, filename: "panel-design",
    csv: () => LAB.recovery,
    render: () => xyLines({ xLog: true, xTitle: "Young people followed (log scale)", yTitle: "Root mean squared error ÷ true value", xFmt: (v) => v.toLocaleString("en-GB"), yFmt: (v) => `${(v * 100).toFixed(0)}%`,
      series: TR.map((t, i) => {
        const r = rows(noise).filter((x) => x.transition === t).sort((a, b) => a.n - b.n);
        return { name: NAME[t], color: COL[i], dash: i % 2 ? "6 3" : "", markers: true, xs: r.map((x) => x.n), ys: r.map((x) => x.rmse / x.true) };
      }) }) });
  const btns = h("div", { class: "row", style: { gap: "8px", marginBottom: "8px" } }, ["low", "high"].map((nz) =>
    h("button", { class: "sm", text: `${nz === "low" ? "Low" : "High"} measurement noise`, onclick: () => { noise = nz; chart.update(); tbl.replaceChildren(makeTable()); } })));
  const makeTable = () => table(["Transition", "True annual probability", "500", "2,000", "8,000"],
    TR.map((t) => {
      const r = rows(noise).filter((x) => x.transition === t).sort((a, b) => a.n - b.n);
      return [NAME[t], r[0].true.toFixed(3), ...r.map((x) => `bias ${fmtPct(x.bias / x.true)}, error ${fmtPct(x.rmse / x.true)}`)];
    }), ["l", "r", "r", "r", "r"]);
  const tbl = h("div", {}, makeTable());

  el.appendChild(panel({ title: "Recovery of the rates by panel size", children: [btns, chart, tbl] }));
  el.appendChild(panel({ title: "What a panel would need", children: [
    h("ul", { class: "prose" }, [
      h("li", { html: "<strong>The full symptom count, not the classification.</strong> Recording only whether a young person is above the HBSC cut-off leaves the latent transitions unrecoverable at any realistic size." }),
      h("li", { html: "<strong>Waves more than size.</strong> For the latent rates, ten annual waves of 500 young people recover incidence better than three waves of 8,000." }),
      h("li", { html: "<strong>A baseline question on past problematic use.</strong> When incidence is low, a first wave cannot tell never-affected from recovered young people; a baseline item with error rates known from validation removes the resulting bias in incidence." }),
      h("li", { html: "<strong>A continuous-time model.</strong> An annual panel observes a process that runs in continuous time; fitting one-step annual transitions biases the hazards, and a continuous-time hidden Markov model removes most of that bias." }),
      h("li", { html: "<strong>Relapse is the hardest rate.</strong> It is identified only through the order of states within a young person's history." }),
    ]),
    h("p", { class: "prose", html: "These requirements define the South African panel study that the thesis proposes: translation and validation of the Social Media Disorder Scale first, then a panel across urban and rural settings sized from this design table." }),
    source("thesis Chapter 6; python/run_all.py and python/ext_recovery.py."),
  ] }));
  return { el };
}
