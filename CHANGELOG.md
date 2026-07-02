# Changelog

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
