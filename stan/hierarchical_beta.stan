// hierarchical_beta.stan --------------------------------------------------
// Hierarchical Bayesian calibration of the country-specific transmission rate
// beta_c across the 44 HBSC countries (Chapter 6.5.4, Eq 6.4 / D.4).
// Partial pooling: log(beta_c) ~ Normal(mu_beta, tau_beta). The likelihood
// links each country's transmission rate to its observed prevalence through a
// deterministic map f() supplied as data (precomputed steady-state or
// short-horizon model output per candidate beta), keeping the Stan model light.

data {
  int<lower=1> C;                 // number of countries
  vector[C] prevalence_obs;       // observed addicted prevalence (proportion)
  vector<lower=0>[C] prevalence_se; // measurement SE per country
  // Precomputed lookup: model-predicted prevalence is approx a*log(beta)+b
  // (linearised f); supply slope/intercept from the deterministic model.
  real map_slope;
  real map_intercept;
}
parameters {
  vector[C] log_beta;             // log transmission rate per country
  real mu_beta;                   // population mean (log scale)
  real<lower=0> tau_beta;         // between-country SD (log scale)
}
transformed parameters {
  vector[C] pred_prev;
  for (c in 1:C)
    pred_prev[c] = map_intercept + map_slope * log_beta[c];
}
model {
  // Priors (weakly informative)
  mu_beta   ~ normal(-3, 1);          // beta ~ 0.05/contact on natural scale
  tau_beta  ~ cauchy(0, 0.5);         // half-Cauchy via <lower=0>
  log_beta  ~ normal(mu_beta, tau_beta);
  // Likelihood
  prevalence_obs ~ normal(pred_prev, prevalence_se);
}
generated quantities {
  vector[C] beta_c = exp(log_beta);   // natural-scale transmission rates
  real beta_pop = exp(mu_beta);       // population-level transmission rate
}
