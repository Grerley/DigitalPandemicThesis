# Changelog

## 2.1.0 (thesis v2, expanded analyses)
Python extension scripts in `python/`, each writing its own `outputs/results_<name>.json`:

- `ext_model.py`: parameter interpretation (durations, exit probabilities, life-course risk)
- `ext_math.py`: closed-form sensitivity indices, Routh–Hurwitz stability of the interior equilibrium, time scales, branching process and Gillespie checks, linear noise approximation, profile likelihoods, simpler structures, ageing against exponential turnover
- `ext_durations.py`: Erlang spell durations and re-calibration
- `ext_agesex.py`, `ext_agesex_shift.py`: age- and sex-structured model, calibration and out-of-sample 2018 check
- `ext_calibration.py`, `ext_bayes.py`, `ext_misclass.py`: further scenarios, plausibility weighting, sampling sensitivity, implicit-prior re-weighting, misclassification
- `ext_recovery.py`: parameter-recovery design grid, baseline item, measurement and continuous-time HMM
- `ext_networks.py`: network extensions
- `ext_interventions.py`, `ext_control.py`: intervention response curves, timing, uncertainty, optimal control
- Tests for the new results in `python/tests/test_model.py`


## 1.0.0 (2026)
Initial public release accompanying the doctoral thesis
*The Digital Pandemic: An Epidemic-Modelling Framework for Technology-Induced Mental Health in a Connected World*.

- S–A–I–R compartmental model with relapse (nonlinear difference-equation form, mass-action force of infection)
- Network Monte Carlo contagion (Erdős–Rényi, Watts–Strogatz, Barabási–Albert)
- Synthetic longitudinal cohort generator (seeded, reproducible)
- Hidden Markov Model transition estimation (validated on synthetic data)
- Logistic regression risk-factor model
- Next-generation-matrix R₀, endemic equilibrium, transcritical-bifurcation analysis
- Hierarchical Bayesian calibration of transmission rate (Stan)
- Intervention scenario analysis (assumed efficacies; comparative outputs)
- Black-and-white journal-style figure generation
- Parameters calibrated so the model reproduces the thesis targets (R₀ ≈ 2.4, peak ≈ 14%, endemic ≈ 8%)
