# The Digital Pandemic

**Reproducible modelling pipeline for the doctoral thesis**
*The Digital Pandemic: An Epidemic-Modelling Framework for Technology-Induced Mental Health in a Connected World*
Grerley Mutibura · University of KwaZulu-Natal · 2026

[![License: MIT](https://img.shields.io/badge/License-MIT-black.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/Data-CC%20BY%204.0-black.svg)](data/README.md)

---

## Overview

This repository contains the full, reproducible modelling pipeline behind the thesis. It models youth
digital addiction as a process of **behavioural contagion** and analyses its dynamics, drivers, and
control using three complementary methods:

1. **An S–A–I–R compartmental model** (Susceptible → At-Risk → Addicted → Recovered) with a **relapse**
   pathway, formulated as a system of nonlinear difference equations with a mass-action force of infection.
2. **Network Monte Carlo simulation** over Erdős–Rényi (random), Watts–Strogatz (small-world), and
   Barabási–Albert (scale-free) topologies, to quantify how social-network structure shapes spread.
3. **Statistical calibration** — hierarchical Bayesian estimation of the transmission rate (Stan),
   Hidden Markov Model estimation of transition probabilities, and logistic regression for risk factors.

## Interactive Simulation Lab

A browser-based, research-grade **interactive lab** ([`app/`](app/)) reimplements this entire
pipeline as a live "flight simulator" for the thesis: explore every model in real time, sweep
parameters, animate network contagion, and compose interventions — with all five model families
interactive, URL-shareable scenarios, CSV/PNG/SVG export, and a built-in acceptance self-test that
verifies R₀ ≈ 2.50 / 3.04 / 2.71, peak ≈ 14% @ ~yr 7, endemic ≈ 8%, and a ≈ −63% comprehensive
intervention.

**▶ Live demo: https://grerley.github.io/DigitalPandemicThesis/** — deployed automatically from
`app/` on every push to `main` by [`.github/workflows/deploy-pages.yml`](.github/workflows/deploy-pages.yml).
(One-time: in *Settings → Pages → Build and deployment*, set the source to **GitHub Actions**.)

Or run it locally — no backend or build step:

```bash
cd app && python3 -m http.server 8000   # then open http://localhost:8000/
```

See [`app/README.md`](app/README.md) for the full description and architecture.

> **Epistemic status.** All quantitative outputs are *conditional model results*, not empirical
> forecasts. Parameters are calibrated against aggregate WHO HBSC and OECD surveillance data together
> with a documented **synthetic cohort** that stands in for the longitudinal transition data that do
> not yet exist. The robust conclusions are qualitative (a self-sustaining epidemic with a critical
> threshold and an exploitable network structure); specific numbers illustrate rather than measure.

## Headline (illustrative) results

| Quantity | Model output | Notes |
|---|---|---|
| Basic reproduction number R₀ | ≈ 2.4 | next-generation-matrix derivation (`06_r0_analysis.R`) |
| Peak addicted prevalence | ≈ 14% (baseline) → higher under scale-free / OECD | `01`, `02` |
| Endemic prevalence | ≈ 8% | positive equilibrium, sustained by relapse |
| Incidence–prevalence lag | ≈ 3 years | prevention window closes early |
| Comprehensive intervention | up to −68% peak | conditional on assumed efficacies (`07`) |

## Quick start

```bash
git clone https://github.com/GrerleyMutibura/DigitalPandemicThesis.git
cd DigitalPandemicThesis
Rscript -e 'install.packages(c("deSolve","igraph","depmixS4","dplyr","tidyr","ggplot2","survival"))'
# (optional, for Bayesian calibration) install cmdstanr / rstan
Rscript analysis/run_all.R      # reproduces core results and figures
```

Outputs (trajectories, calibrated parameters, and black-and-white journal-style figures) are written to
`outputs/` and `figures/`.

## Repository layout

```
.
├── R/                     Model + analysis source
│   ├── utils.R                shared helpers, default parameters
│   ├── 01_sair_model.R        deterministic S–A–I–R difference-equation model
│   ├── 02_network_simulation.R  ER / WS / BA Monte Carlo contagion
│   ├── 03_synthetic_data.R    synthetic longitudinal cohort generator
│   ├── 04_hmm_estimation.R    Hidden Markov transition estimation
│   ├── 05_regression.R        logistic regression of risk factors
│   ├── 06_r0_analysis.R       R₀, next-generation matrix, endemic equilibrium, stability
│   └── 07_interventions.R     intervention scenario simulation
├── stan/
│   └── hierarchical_beta.stan hierarchical Bayesian calibration of β
├── data/
│   ├── parameters.csv         calibrated parameters with uncertainty ranges
│   └── README.md              data provenance and access (HBSC/OECD are external)
├── figures/make_figures.R     regenerates all thesis figures (B&W)
├── analysis/run_all.R         master pipeline
├── docs/
│   ├── model_specification.md  full equations (4.1–4.4, A.1–A.7)
│   └── methods_summary.md      methods & calibration overview
├── CITATION.cff
├── LICENSE
└── .gitignore
```

## Model summary

The population state vector **x**(t) = [S, A, I, R]ᵀ evolves under a state-dependent transition matrix
**P**(**x**ₜ) whose S→A entry is a mass-action force of infection

```
P_SA(x_t) = 1 − exp[ −β·k·(I_t + η·A_t) ] + μ
```

with per-contact transmission β, mean active contacts k, reduced At-Risk infectivity η ∈ [0,1], and
exogenous acquisition μ. The Recovered→At-Risk relapse term δ keeps the epidemic self-sustaining and
produces a positive endemic equilibrium; the disease-free equilibrium loses stability through a
transcritical bifurcation at R₀ = 1. Full derivations are in [`docs/model_specification.md`](docs/model_specification.md).

## Data availability

No identifiable human-subjects data are held in this repository. Country-level prevalence targets derive
from the publicly reported **WHO Health Behaviour in School-aged Children (HBSC)** study and **OECD**
sources; access instructions are in [`data/README.md`](data/README.md). The synthetic cohort is fully
regenerable from `R/03_synthetic_data.R` with a fixed seed.

## Citation

If you use this code, please cite the thesis (see [`CITATION.cff`](CITATION.cff)).

## Licence

Code is released under the MIT Licence ([`LICENSE`](LICENSE)). Documentation and generated data tables
are under CC BY 4.0.
