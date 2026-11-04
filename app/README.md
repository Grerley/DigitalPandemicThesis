# Digital Pandemic Lab

An interactive, browser-based companion to the doctoral research *An epidemic-modelling framework for problematic
social media use among young people* (Grerley Mutibura, School of Computer Science and Applied Mathematics,
University of the Witwatersrand). It runs entirely in the browser, with no back end and no build step.

> Every output is a model result conditional on stated assumptions. Projections are not forecasts, intervention
> effects describe relative leverage inside the model under assumed efficacies, and the exposure rate is not a
> causal peer effect.

## Where the numbers come from

- **Thesis results** (projection bands, likelihood ratios for the next HBSC wave, intervention, network, control,
  panel-design and age-sex results) are read from `js/model/lab-data.js`, which `python/export_lab.py` generates
  from the outputs of the thesis analysis in `python/`. The file records the commit it was generated from.
- **Live results** are computed by a JavaScript port of the thesis model (`js/model/`): the S-A-I-R equations
  integrated by an adaptive Dormand–Prince method, the closed-form threshold, equilibrium and floor, the exact
  calibration to the two HBSC figures, the monthly network process and the intervention scenarios.
- The **self-test** (button in the app, or `node app/test/selftest.mjs`) checks the port against the Python values:
  thresholds, equilibria, the floor, sensitivity indices, the exactly solved exposure rates, trajectories under both
  scenarios, projection medians across all kept draws, network transmissibility, and network runs against the
  thesis's distributions. The deployment workflow runs it before publishing.

To refresh the data after the analysis changes: run the thesis scripts, then `python python/export_lab.py`.

## Screens

1. **Research context**: the question, the measure, the model and how to read the numbers.
2. **Two explanations**: the persistent shift and the temporary shock, projection bands across all kept draws, any
   single draw computed live, and what each possible 2025/2026 HBSC figure would imply under three priors.
3. **Model explorer**: sliders over the thesis's parameter ranges, live R₀ (with its relapse multiplier),
   equilibrium, floor, stability, flows and sensitivity indices, and exact calibration of any set of rates.
4. **Thresholds and floor**: why a non-social route removes the threshold, and what sets the floor.
5. **Networks**: paired experiments on Erdős–Rényi, Watts–Strogatz and Barabási–Albert networks with random,
   degree-targeted or acquaintance protection, an animated small network, and the thesis results.
6. **Interventions**: four mechanisms with efficacies inside the assumed ranges, combinations across all draws,
   Bliss indices and the optimal-control result.
7. **Panel design**: how well a hidden Markov model recovers the rates from simulated panels, and what a panel needs.
8. **Age and sex**: the structured model, the 2018 out-of-sample test and a minimum age of 16 by age band.
9. **Methods and sources**: equations, parameter ranges and sources, scope of the claims, references.

## Run locally

```bash
cd app && python3 -m http.server 8000     # then open http://localhost:8000/
node test/selftest.mjs                    # model checks
node test/lint.mjs                        # syntax check of every file
```

## Earlier version

The first version of the lab, built on the first version of the thesis, is kept for reference in
`archive/simulation-lab-v1/`. Its parameters, starting state, intervention
effects and cost framing are not consistent with the current thesis.
