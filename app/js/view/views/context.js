// context.js: Screen 1, what the research asks and how to read the lab.

import { h } from "../dom.js";
import { panel, metricGrid } from "../ui.js";
import { viewHead, pc } from "../common.js";
import { R, LAB } from "../../model/params.js";

export function ContextView() {
  const el = h("div", { class: "view" });
  el.appendChild(viewHead("Screen 1 · Research context", "What the research asks",
    "This lab is the interactive companion to the doctoral research <em>An epidemic-modelling framework for problematic social media use among young people</em> " +
    "(Grerley Mutibura, School of Computer Science and Applied Mathematics, University of the Witwatersrand). Every number on its screens is produced by the same model and analysis code as the thesis."));

  el.appendChild(panel({ title: "The starting point", children: [
    h("p", { class: "prose", html:
      "In the Health Behaviour in School-aged Children (HBSC) study, the share of 11-, 13- and 15-year-olds with problematic social media use rose from <strong>7% in 2017/2018 to 11% in 2021/2022</strong> (Boniel-Nissim et al., 2024). " +
      "Problematic use means six or more of the nine features on the Social Media Disorder Scale (van den Eijnden et al., 2016), such as preoccupation, failed attempts to cut down and conflict with others. " +
      "It is a measured behavioural category, not a clinical diagnosis." }),
    h("p", { class: "prose", html:
      "A prevalence changes only through flows: young people starting problematic use, stopping, relapsing, and entering or leaving the age group. Two survey figures cannot say which flow changed. " +
      "South Africa does not take part in HBSC and has no measure of problematic social media use among adolescents, so its results here are scenarios and a research design, not estimates." }),
    metricGrid([
      { label: "HBSC 2018", value: "7%", foot: "11-, 13- and 15-year-olds" },
      { label: "HBSC 2022", value: "11%", foot: "girls 13%, boys 9%" },
      { label: "Kept calibration draws", value: R("cal_n_ok").toLocaleString("en-GB"), foot: `of ${R("cal_n_draws").toLocaleString("en-GB")} drawn` },
      { label: "Equilibrium floor at b = 0", value: `${R("I_star_b_zero").toFixed(2)}%`, foot: "central parameters", accent: true },
    ]),
  ] }));

  el.appendChild(panel({ title: "The model in one paragraph", children: [
    h("p", { class: "prose", html:
      "Young people aged 10 to 24 are Susceptible (S), At-Risk (A), Addicted (I, the state the survey observes as problematic use) or Recovered (R). " +
      "Susceptible young people move into risk through <strong>social exposure</strong> to At-Risk and problematic users, at rate <var>b</var>, and through a <strong>non-social route</strong>, at rate <var>α</var>, which stands for influences such as platform design. " +
      "At-Risk young people progress (<var>σ</var>) or remit (<var>γ</var>₁); problematic use ends in recovery (<var>γ</var>₂); recovery can end in relapse (<var>δ</var>). Young people enter at 10 and leave at 25 (<var>ν</var>)." }),
    h("div", { class: "caveat", html:
      "<strong>The epidemic language is a modelling lens.</strong> Similarity among friends can come from influence, from choosing similar friends or from a shared environment, and observational data cannot in general separate these (Manski, 1993; Shalizi &amp; Thomas, 2011). " +
      "The exposure rate <var>b</var> bundles all three. It is not a measured peer effect, and nothing in the lab claims that problematic use is transmitted or that social media causes poor mental health." }),
  ] }));

  const rq = [
    ["Thresholds and the floor", "What thresholds, equilibria and stability does the model have, and how does a non-social route change them?", "Model explorer; Thresholds and floor"],
    ["What the data can decide", "What can two published figures identify, and what do competing explanations of the rise predict for the next survey?", "Two explanations"],
    ["Structure", "How do network structure, age and sex change the dynamics and the value of targeting?", "Networks; Age and sex"],
    ["Interventions and data", "Under stated efficacy assumptions, which mechanisms have most leverage, and what panel study would estimate the rates?", "Interventions; Panel design"],
  ];
  el.appendChild(panel({ title: "Research questions and where to explore them", children: [
    h("table", { class: "data" }, [
      h("thead", {}, h("tr", {}, ["Theme", "Question", "Screen"].map((t) => h("th", { style: { textAlign: "left" }, text: t })))),
      h("tbody", {}, rq.map((r) => h("tr", {}, r.map((c) => h("td", { style: { textAlign: "left" }, text: c }))))),
    ]),
  ] }));

  el.appendChild(panel({ title: "How to read the numbers", children: [
    h("p", { class: "prose", html:
      "Two kinds of output appear. <strong>Thesis results</strong> are read from the analysis outputs (for example the projection bands across all " + R("cal_n_ok") + " kept parameter draws) and are labelled with their source. " +
      "<strong>Live results</strong> are computed in your browser by a JavaScript port of the thesis model, which the self-test (top right) checks against the thesis values. " +
      "Both are conditional on stated assumptions. Projections are model output, not forecasts, and intervention effects describe relative leverage inside the model under assumed efficacies, not what any programme would achieve." }),
    h("p", { class: "hint", text: `Data generated from the thesis analysis, commit ${LAB.meta.commit}, on ${LAB.meta.generated}.` }),
  ] }));
  return { el };
}
