# Thesis v2 pipeline (Python)

This folder reproduces every number, table and figure in the corrected thesis
*An Epidemic-Modelling Framework for Problematic Social Media Use among Young People*
(G. Mutibura, University of the Witwatersrand, 2026).

```bash
pip install numpy scipy networkx matplotlib
python run_all.py        # about 12 minutes; writes outputs/results.json and figures/
python derived.py        # adds recovery-study bias and paired targeting runs to results.json
python -m pytest tests   # numerical checks of every theorem in Chapter 5
```

* `dpm/model.py` S-A-I-R ODE with turnover, R0, equilibria, sensitivity
* `dpm/network.py` stochastic network version and network reproduction number
* `dpm/calibrate.py` calibration to the published HBSC 2018 and 2022 figures under two scenarios
* `dpm/recovery.py` hidden Markov parameter-recovery study (synthetic data with known truth)
* `data/parameter_ranges.json` literature-based parameter ranges (sources in the thesis, Table 7.1)
* `data/sa_access.json` South African connected-share scenarios (sources in the thesis, Table 7.4)

All outputs are model outputs. The synthetic data in the recovery study are generated from known
values and are not estimates for any population. The R code in `../R` belongs to the first version
of the thesis and is not used for any result in this version.
