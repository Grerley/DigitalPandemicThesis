# The Digital Pandemic — Interactive Simulation Lab

A research-grade, browser-based "flight simulator" for the dynamical models in the
thesis *The Digital Pandemic: An Epidemic-Modelling Framework for Technology-Induced
Mental Health in a Connected World* (Mutibura, 2026). It reimplements the R/Stan
pipeline (`../R/`) as a live, interactive web application — no backend, no build step.

> **Epistemic status.** Every quantity on screen is a **conditional model result**,
> not an empirical forecast. See the *Methods & Provenance* screen for the equations,
> parameter sources, and caveats.

## Run it

The app uses native ES modules and a Web Worker, so it must be served over HTTP
(opening `index.html` directly from `file://` disables the worker and module loading).
Any static server works:

```bash
cd app
python3 -m http.server 8000
# then open http://localhost:8000/
```

No dependencies, no `npm install`. Runs fully client-side in any modern browser.

## What's inside

Seven screens, each a view over the shared, audited model layer:

1. **Overview** — the calibrated S–A–I–R baseline: four-compartment trajectory, the
   R₀ trio, peak/endemic readouts, and an animated state-transition diagram whose edge
   thickness tracks live monthly transition rates. Scrub the timeline.
2. **Parameter Explorer** — sliders for all eight monthly rates with live R₀ / peak /
   endemic, a calibrated baseline overlay, one-way sensitivity (tornado), and reset.
3. **Phase & Equilibrium** — transcritical bifurcation (I\* vs basic R₀), phase portrait
   with vector field, and the incidence-vs-prevalence lead.
4. **Network Lab** — stochastic Reed–Frost contagion on Erdős–Rényi, Watts–Strogatz and
   Barabási–Albert graphs. Animated spread on a rendered network, Monte-Carlo ensembles
   (mean + 95% envelope) run in a Web Worker, and hub-vs-random node immunisation.
5. **Cohort & Risk Factors** — a synthetic longitudinal cohort, a risk-factor odds-ratio
   forest plot recovered by logistic regression, a state-flow Sankey, and a
   known-vs-estimated transition-matrix comparison.
6. **Intervention Simulator** — compose interventions (multiplicative parameter changes),
   overlay the baseline, a peak-reduction table, and an illustrative ZAR cost framing.
7. **Methods & Provenance** — the model spec, rendered equations, parameter table with
   sources, epistemic caveats, and a live acceptance self-test panel (PASS/FAIL).

## Architecture

The **model layer** (`js/model/`) is pure, framework-free, and decoupled from the UI so
the science can be audited independently. It is a faithful port of the R reference and is
verified against the thesis acceptance targets by `js/model/selftest.js`.

```
js/
  model/            ← scientific core (pure, no DOM)
    params.js         calibrated parameters, metadata, presets
    rng.js            seeded RNG (mulberry32) — same seed ⇒ same result
    sair.js           deterministic S–A–I–R difference-equation engine
    r0.js             R₀ (reduced/full NGM, crude), equilibrium, bifurcation
    network.js        graph generators + Reed–Frost contagion + ensembles
    cohort.js         synthetic cohort, logistic regression (IRLS), matrices
    interventions.js  intervention presets & multiplicative modifiers
    selftest.js       acceptance-target assertions
  view/             ← rendering, charts, screens
    charts.js         reusable SVG charts (time-series, bifurcation, phase,
                      tornado, forest, heatmap, Sankey) — export-ready
    diagram.js        animated state-transition diagram
    netgraph.js       force layout + canvas network renderer
    views/            one module per screen
    export.js         CSV / SVG / PNG export
  state.js          observable store: scenario state, URL sharing, localStorage
  network-worker.js Web Worker entry for ensembles (main-thread fallback provided)
  app.js            shell: routing, scenario bar, theme
```

## Reproducibility & correctness

- The deterministic core reproduces the thesis targets exactly: **R₀ ≈ 2.50 / 3.04 / 2.71**,
  **peak ≈ 14% at ~year 7**, **endemic ≈ 8%**, **comprehensive intervention ≈ −63% peak**.
  These are asserted in the self-test panel (Methods screen or the "Self-test" button).
- All stochastic components draw from a seeded RNG — identical seed ⇒ identical result.
- Transition-matrix rows sum to 1, all probabilities are clamped to [0, 1], and the state
  vector is renormalised every step.

## Sharing & export

- Every scenario is encoded in the URL query string — copy the link (Share button) to
  reproduce an exact configuration.
- Trajectories export to **CSV**; every chart exports to **SVG** and **PNG**.
- Presets: *Calibrated baseline*, *Scale-free outbreak*, *Comprehensive intervention*,
  plus user "save scenario" slots (localStorage) and a two-scenario compare overlay.

## Accessibility

Keyboard-navigable, colour-blind-safe encodings paired with line style (never colour
alone), light/dark themes, and full support for `prefers-reduced-motion`.
