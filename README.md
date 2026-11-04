# An epidemic-modelling framework for problematic social media use among young people

**Code, analysis outputs and interactive lab for the doctoral thesis of Grerley Mutibura**
School of Computer Science and Applied Mathematics, University of the Witwatersrand

[![License: MIT](https://img.shields.io/badge/License-MIT-black.svg)](LICENSE)

---

## Overview

In the Health Behaviour in School-aged Children (HBSC) study, problematic social media use among 11-, 13- and
15-year-olds rose from 7% in 2017/2018 to 11% in 2021/2022 (Boniel-Nissim et al., 2024). Problematic use means six
or more of the nine items of the Social Media Disorder Scale; it is a measured behavioural category, not a diagnosis.

The thesis builds a four-state compartmental model of young people aged 10 to 24 (Susceptible, At-Risk, Addicted,
Recovered) with population turnover, a social and a non-social route into risk, remission, recovery and relapse. It

- proves the model's threshold, stability and bifurcation properties, and shows that any non-social acquisition
  removes the threshold and sets a floor that measures acting only on social exposure cannot pass;
- calibrates the model to the two HBSC figures by solving exactly the two quantities they identify, makes the prior
  implied by that design explicit, and states in advance what the next HBSC wave would imply;
- designs, by a hidden Markov parameter-recovery study, the panel data that would identify the transition rates;
- examines network structure and targeting, intervention mechanisms under assumed efficacies, optimal control,
  and an age- and sex-structured version tested against the 2018 survey.

The social exposure rate bundles influence, homophily and shared environment and is not a causal peer effect.
The thesis makes no causal claim about social media and mental health. South Africa, which does not take part in
HBSC, is treated through scenarios and a research design, not estimates.

## Selected results

All are model outputs conditional on stated assumptions, produced by the code in `python/`.

| Quantity | Value |
|---|---|
| Parameter draws reproducing both HBSC figures | 805 of 1,000 |
| Projected problematic use in 2042, persistent shift in exposure (median, 90% range) | 13.2% (9.8% to 21.3%) |
| Projected problematic use in 2042, temporary shock in 2020 and 2021 (median, 90% range) | 7.3% (6.6% to 8.6%) |
| Floor of prevalence with no social exposure, central parameters | 0.59% |
| R₀ after 2018 without non-social acquisition, central parameters | 4.02 |
| Reduction in 2042 prevalence, treatment access (median, assumed efficacy) | 37% |
| Reduction in 2042 prevalence, platform design regulation (median, assumed efficacy) | 13% |

## Interactive lab

A browser-based lab ([`app/`](app/)) lets readers explore the model: the two explanations of the rise and what the
next wave can decide, the threshold and the floor, networks and targeting, interventions, panel design, and age and
sex. Thesis results are generated from the Python analysis by `python/export_lab.py`; live results come from a
JavaScript port of the model that a built-in self-test checks against the Python values.

**Live: https://grerley.github.io/DigitalPandemicThesis/**, deployed from `app/` on every push to `main` by
[`.github/workflows/deploy-pages.yml`](.github/workflows/deploy-pages.yml). See [`app/README.md`](app/README.md).

## Reproduce

```bash
cd python
pip install -r requirements.txt
python run_all.py          # calibration, projections, mathematics, networks, interventions, recovery study
python derived.py
python ext_*.py            # extension analyses, each writing outputs/results_<name>.json
python export_lab.py       # data file for the interactive lab
python -m pytest tests     # numerical checks of the theorems
```

See [`python/README.md`](python/README.md). Running any script with `DPM_BW=1` renders black-and-white figures.

## Repository layout

```
.
├── python/                 the thesis pipeline (model, calibration, networks, recovery study, interventions)
│   ├── dpm/                    model package
│   ├── data/                   literature-based parameter ranges and scenario inputs
│   ├── outputs/                every number reported in the thesis (results*.json)
│   ├── figures/, figures_bw/   figures
│   └── tests/                  numerical checks of the mathematical results
├── app/                    interactive lab (static site, no build step)
├── archive/
│   └── simulation-lab-v1/  first version of the lab, kept for reference; not consistent with the thesis
├── R/, stan/, analysis/, figures/, docs/, data/
│                           code and documents from the first version of the thesis; not used for any result
└── .github/workflows/      Python tests and lab deployment
```

## Data

No individual-level or identifiable data are held here. The calibration uses only published HBSC aggregates;
South African figures are published household and child-user statistics used for context and scenarios. The
parameter-recovery study uses synthetic data generated from known values, which are not estimates for any
population.

## Citation and licence

Please cite the thesis (see [`CITATION.cff`](CITATION.cff)). Code is released under the MIT Licence ([`LICENSE`](LICENSE)).
