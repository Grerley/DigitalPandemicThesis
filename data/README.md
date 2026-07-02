# Data

This repository contains **no identifiable human-subjects data**. It holds only (a) a small table of
calibrated model parameters and (b) code to regenerate a synthetic cohort.

## `parameters.csv`
Calibrated S–A–I–R parameters with 95% uncertainty ranges and provenance (Chapter 6, Table 6.4).
Rates are given in their natural units (per month or per year as noted); the code converts to
monthly transition probabilities via `rate_to_prob()`.

## Synthetic cohort
The individual-level longitudinal data used for Hidden Markov Model validation are **synthetic**,
produced by `R/03_synthetic_data.R` with a fixed seed (`generate_cohort(seed = 42)`). They are a
modelling instrument that reproduces empirical prevalence and risk-factor patterns while allowing
direct computation of transition probabilities; they are not observations of real individuals.

## External data (not redistributed here)
Calibration targets derive from publicly reported aggregate statistics. Obtain them directly from the
custodians:

| Source | Use | Access |
|---|---|---|
| WHO **HBSC** (Health Behaviour in School-aged Children) | country-level prevalence, trends | https://hbsc.org (Data Use Agreement) |
| **OECD** child well-being / digital indicators | cross-national access & context | https://data.oecd.org |
| CDC **YRBS** | US adolescent behaviour | https://www.cdc.gov/yrbs |
| **Statistics South Africa** | SA context | https://www.statssa.gov.za |

Only aggregate, published figures are used; no data-use agreement is required to reproduce the
calibration targets from the cited reports.
