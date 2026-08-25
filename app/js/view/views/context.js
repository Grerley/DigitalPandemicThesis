// context.js — Screen 1: Research Context.
// -----------------------------------------------------------------------------
// Orientation for the whole lab: the framing, the research questions, the method
// each question is answered with, the data behind the calibration, a reading
// guide to the other screens, and the epistemic status. Content is drawn from
// the thesis and the repository docs (methods_summary.md, model_specification.md,
// README) — no numbers are invented here.

import { h } from "../dom.js";
import { panel, metricGrid, condNote } from "../ui.js";

export function ContextView() {
  const el = h("div", { class: "view" });

  el.appendChild(h("div", { class: "view-head" }, [
    h("div", { class: "kicker", text: "Screen 1 · Research Context" }),
    h("h1", { text: "Research Context" }),
    h("p", { class: "prose", html: "This lab is the interactive companion to the doctoral thesis <em>The Digital Pandemic: An Epidemic-Modelling Framework for Technology-Induced Mental Health in a Connected World</em> (Mutibura, University of KwaZulu-Natal, 2026). It treats youth digital addiction not as a fixed trait but as a <strong>behavioural contagion</strong> — something that spreads through social contact — and asks what happens if we model it with the mathematics built for epidemics." }),
  ]));

  // At-a-glance orientation strip (all illustrative model outputs)
  el.appendChild(panel({
    title: "The headline, at a glance",
    sub: "illustrative baseline outputs — conditional model results, not empirical forecasts",
    badge: "conditional",
    children: [
      metricGrid([
        { label: "Basic reproduction no. R₀", value: "≈ 2.5", foot: "> 1 ⇒ self-sustaining", accent: true },
        { label: "Peak addicted prevalence", value: "≈ 14%", foot: "baseline, around year 7" },
        { label: "Endemic prevalence", value: "≈ 8%", foot: "sustained by relapse" },
        { label: "Incidence–prevalence lag", value: "≈ 3 yr", foot: "prevention window closes early" },
      ]),
    ],
  }));

  // Why this framing
  el.appendChild(panel({ title: "Why model addiction as an epidemic?", children: [
    h("p", { class: "prose", html: "Adolescent problematic technology use rose steeply and in step across many countries — the pattern of a spreading process, not of isolated individual choices. Peers influence one another; platforms are engineered to maximise engagement; norms and behaviours diffuse through social networks. Those are exactly the ingredients epidemiology was built to reason about: a <strong>force of infection</strong> that depends on how many others are already affected, a <strong>threshold</strong> that separates a fizzle from an outbreak, and a <strong>network structure</strong> that decides who is exposed to whom." }),
    h("p", { class: "prose", html: "The pay-off of the framing is leverage. If technology-induced mental-health harm behaves like a contagion, then the epidemiological toolkit — R₀, herd effects, targeted immunisation, endemic equilibria — becomes available for <em>understanding and controlling</em> it. The thesis makes that analogy precise and then tests where it holds and where it breaks." }),
    h("div", { class: "caveat", html: "<strong>The analogy is a modelling lens, not a medical claim.</strong> “Addiction” here is an operational construct (meeting problematic-use thresholds), not a settled clinical diagnosis, and “transmission” bundles peer influence, homophily and shared context together — it is not identified causal peer transmission." }),
  ] }));

  // Research questions
  el.appendChild(panel({ title: "Research questions", sub: "each is answered by a specific method — and a screen in this lab", children: [
    rqList([
      { q: "Can technology-induced mental-health harm be represented as a self-sustaining epidemic — and if so, does it cross a critical threshold?",
        method: "Deterministic S–A–I–R compartmental model with a mass-action force of infection and a relapse pathway.",
        screens: [["overview", "Overview"], ["phase", "Phase & Equilibrium"]] },
      { q: "What determines whether an outbreak takes off, peaks, or settles into a stable endemic level?",
        method: "Reproduction-number (next-generation-matrix) analysis, endemic-equilibrium and bifurcation structure; parameter sensitivity.",
        screens: [["explorer", "Parameter Explorer"], ["phase", "Phase & Equilibrium"]] },
      { q: "How does the structure of the social network shape how far and how fast the behaviour spreads?",
        method: "Stochastic Reed–Frost Monte-Carlo simulation over random, small-world and scale-free graphs.",
        screens: [["network", "Network Lab"]] },
      { q: "Who is most at risk, and can the underlying transition process be recovered from observable data?",
        method: "Synthetic longitudinal cohort, logistic regression of risk factors, and Hidden-Markov recovery of the transition matrix.",
        screens: [["cohort", "Cohort & Risk"]] },
      { q: "Which interventions reduce harm most — and does targeting the network beat treating individuals?",
        method: "Scenario analysis mapping education, regulation, age-gating and treatment to parameter changes; hub-targeted vs random immunisation.",
        screens: [["interventions", "Intervention Simulator"], ["network", "Network Lab"]] },
    ]),
  ] }));

  // Approach — three methods
  el.appendChild(panel({ title: "The approach — three complementary methods", children: [
    h("p", { class: "prose", html: "No single technique answers all of the questions above, so the thesis triangulates with three:" }),
    h("table", { class: "data" }, [
      h("thead", {}, h("tr", {}, [
        h("th", { style: { textAlign: "left" }, text: "Method" }),
        h("th", { style: { textAlign: "left" }, text: "What it contributes" }),
        h("th", { style: { textAlign: "left" }, text: "Explore on" }),
      ])),
      h("tbody", {}, [
        methodRow("Compartmental dynamics", "S–A–I–R nonlinear difference equations with a state-dependent transition matrix. Generates the rise–peak–decline curve, the endemic plateau, R₀ and the threshold.", [["overview", "Overview"], ["phase", "Phase"]]),
        methodRow("Network Monte-Carlo", "Reed–Frost transmission over Erdős–Rényi, Watts–Strogatz and Barabási–Albert graphs. Shows how topology — especially hubs — changes the outbreak and where intervention bites.", [["network", "Network Lab"]]),
        methodRow("Statistical calibration", "Hierarchical Bayesian estimation of the transmission rate (Stan, 44 HBSC countries), Hidden-Markov transition estimation, and logistic risk-factor regression.", [["cohort", "Cohort & Risk"], ["methods", "Methods"]]),
      ]),
    ]),
  ] }));

  // Data & calibration
  el.appendChild(panel({ title: "Data & calibration", children: [
    h("ul", { class: "prose" }, [
      h("li", { html: "<strong>WHO HBSC 2022</strong> (Health Behaviour in School-aged Children) — cross-national adolescent prevalence and the initial state distribution, across 44 countries." }),
      h("li", { html: "<strong>OECD surveillance priors</strong> — used for the higher-exposure comparison arm and for the extrapolated South African arm (adjusted from OECD priors, not independently estimated)." }),
      h("li", { html: "<strong>A documented synthetic cohort</strong> — a seeded, fully regenerable modelling instrument that stands in for the individual-level <em>longitudinal transition data that do not yet exist</em>. It reproduces the empirical prevalence and risk-factor odds ratios, enabling direct transition-probability computation. It is not real data." }),
    ]),
    h("p", { class: "prose", html: "Rates are calibrated as monthly probabilities and converted with <span class=\"eq-inline\">P = 1 − exp(−rate·Δt)</span>. The full parameter table, equations and sources live on the <a href=\"#methods\">Methods &amp; Provenance</a> screen; every number elsewhere in the lab traces back to it." }),
  ] }));

  // How to read this lab
  el.appendChild(panel({ title: "How to read this lab", sub: "eight screens, front to back", children: [
    h("ol", { class: "prose readguide" }, [
      guide("overview", "Overview", "the baseline epidemic — one run, all four compartments over time."),
      guide("explorer", "Parameter Explorer", "move the parameters yourself and watch R₀ and the curve respond."),
      guide("phase", "Phase & Equilibrium", "the dynamical skeleton: threshold, bifurcation, endemic point, incidence–prevalence lag."),
      guide("network", "Network Lab", "the same contagion on real network structures, with targeted vs random immunisation."),
      guide("cohort", "Cohort & Risk", "who gets addicted (risk factors) and how well the hidden transition process can be recovered."),
      guide("interventions", "Intervention Simulator", "compose policy levers and read off the change in peak harm."),
      guide("methods", "Methods & Provenance", "equations, parameters, sources, the live self-test, and every caveat."),
    ]),
    h("p", { class: "hint", html: "The bar at the top of each screen lets you apply <strong>presets</strong>, <strong>pin a scenario to compare</strong>, and <strong>save or share</strong> a setup by URL. Nothing here needs a server or an account." }),
  ] }));

  // Glossary
  el.appendChild(panel({ title: "Key terms", children: [
    h("dl", { class: "glossary" }, [].concat(...[
      ["S · Susceptible", "Not problematically engaged, but exposed and at risk."],
      ["A · At-Risk", "Heavy or pre-clinical use; elevated risk and partially transmitting."],
      ["I · Addicted", "Meets problematic-use thresholds (operational construct, not a diagnosis)."],
      ["R · Recovered", "Reduced use — but not absorbing: subject to relapse back toward At-Risk."],
      ["R₀ (basic reproduction number)", "Expected new cases generated by one case in a fully susceptible population. R₀ > 1 means the behaviour can sustain itself."],
      ["Force of infection", "The per-step probability of onset, rising with how many contacts are already At-Risk or Addicted — the source of the nonlinearity."],
      ["Endemic equilibrium", "The non-zero level at which the epidemic stabilises when R₀ > 1, held up here by relapse."],
      ["Relapse", "The R → A pathway. Removing it makes Recovery absorbing and the epidemic dies out — it is what keeps the plateau alive."],
    ].map(([t, d]) => [h("dt", { text: t }), h("dd", { class: "prose", text: d })]))),
  ] }));

  // Epistemic status + citation
  el.appendChild(panel({ title: "Epistemic status & how to cite", children: [
    condNote("Every quantitative output in this lab is a result <em>conditional on the model and its assumed parameters</em>. The robust findings are qualitative — a self-sustaining epidemic with a critical threshold and an exploitable network structure. Specific numbers illustrate the mechanism; they do not measure the world."),
    h("p", { class: "prose", style: { marginTop: "12px" }, html: "Intervention efficacies are <strong>assumed</strong> (informed by the direction of published evaluations, not estimated here), risk-factor odds ratios are <strong>associational</strong>, and cost figures are illustrative unit costs. See <a href=\"#methods\">Methods &amp; Provenance</a> for the complete list of caveats." }),
    h("p", { class: "prose", html: "<strong>Cite as:</strong> Mutibura, G. (2026). <em>The Digital Pandemic: An Epidemic-Modelling Framework for Technology-Induced Mental Health in a Connected World</em> (PhD thesis, University of KwaZulu-Natal). This lab reimplements the thesis pipeline (R/Stan) in the browser; the deterministic core is verified against the acceptance targets on the Methods screen." }),
  ] }));

  return { el, destroy: () => {} };
}

// A research-question card with the method that answers it and screen links.
function rqList(items) {
  return h("ol", { class: "rq-list" }, items.map((it, i) =>
    h("li", { class: "rq" }, [
      h("div", { class: "rq-q" }, [h("span", { class: "rq-n", text: `RQ${i + 1}` }), h("span", { class: "prose", text: it.q })]),
      h("div", { class: "rq-m prose", html: `<span class="rq-lab">Method</span> ${it.method}` }),
      h("div", { class: "rq-links" }, [h("span", { class: "rq-lab", text: "Explore" })].concat(
        it.screens.map(([id, label]) => h("a", { class: "chip", href: `#${id}`, text: label + " →" })))),
    ])));
}

function methodRow(name, contributes, screens) {
  return h("tr", {}, [
    h("td", { style: { textAlign: "left", fontWeight: "600" }, text: name }),
    h("td", { style: { textAlign: "left" }, class: "prose", text: contributes }),
    h("td", { style: { textAlign: "left" } }, h("div", { class: "rq-links" },
      screens.map(([id, label]) => h("a", { class: "chip", href: `#${id}`, text: label + " →" })))),
  ]);
}

function guide(id, title, desc) {
  return h("li", {}, [
    h("a", { class: "rg-link", href: `#${id}` }, [h("strong", { text: title })]),
    h("span", { class: "prose", html: ` — ${desc}` }),
  ]);
}
