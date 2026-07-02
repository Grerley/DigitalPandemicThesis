# Methods Summary

A short guide to the estimation and simulation methods, mapping each thesis method to its code.

## 1. Deterministic dynamics — `R/01_sair_model.R`
The population state distribution evolves under a state-dependent transition matrix with a mass-action
force of infection. Produces the baseline trajectory (peak ≈13% around year 8, endemic ≈8%).

## 2. Network contagion — `R/02_network_simulation.R`
Stochastic Reed–Frost transmission on Erdős–Rényi, Watts–Strogatz and Barabási–Albert graphs
(`igraph`). Monte Carlo ensembles give mean trajectories with 95% envelopes. Scale-free networks peak
higher and earlier; hub-targeted intervention beats random at equal coverage.

## 3. Synthetic cohort — `R/03_synthetic_data.R`
Generates individual-level longitudinal trajectories reproducing empirical prevalence and risk-factor
odds ratios, enabling direct transition-probability computation. A modelling instrument, seeded and
fully regenerable; not real data.

## 4. Hidden Markov Models — `R/04_hmm_estimation.R`
Estimates the 4-state transition matrix from longitudinal observations (screen time + symptom score)
using `depmixS4` (Baum–Welch / EM; Viterbi decoding). Validated by recovering the known generating
matrix from synthetic data (mean absolute error < 0.05).

## 5. Regression — `R/05_regression.R`
Logistic regression of the Addicted state on age, gender, urbanicity, SES, baseline depression and
digital literacy, reproducing Table 6.2 odds ratios and mapping them to heterogeneous S→A
probabilities. Associational, not causal.

## 6. R₀, equilibria, stability — `R/06_r0_analysis.R`
Crude and next-generation-matrix R₀ (≈2.5), networked R₀ with degree heterogeneity, and the endemic
equilibrium with a transcritical bifurcation at R₀ = 1.

## 7. Bayesian calibration — `stan/hierarchical_beta.stan`
Hierarchical partial-pooling estimation of country-specific transmission rates across 44 HBSC
countries: log(β_c) ~ Normal(μ_β, τ_β), with weakly informative hyperpriors. Fit via Stan
(cmdstanr/rstan).

## 8. Interventions — `R/07_interventions.R`
Scenario analysis mapping education, platform regulation, age-gating and treatment to parameter
changes. **Efficacies are assumed, not estimated** — outputs are comparative, not predictive.

## Reproducibility
`Rscript analysis/run_all.R` runs the full pipeline. Random components use fixed seeds. Package
versions are pinned in `DESCRIPTION` (and, if used, `renv.lock`).

## Epistemic caveats (carried from the thesis)
- All numeric outputs are **conditional model results**, not empirical forecasts.
- Calibration relies partly on a synthetic cohort standing in for absent longitudinal data.
- The transmission coefficient captures peer influence, homophily and shared context jointly; it is
  not identified causal peer transmission.
- The South African arm is extrapolated from OECD priors with adjustment factors, not independently
  estimated.
