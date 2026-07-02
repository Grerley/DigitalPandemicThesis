# 05_regression.R ---------------------------------------------------------
# Logistic regression for risk factors and heterogeneous transition
# probabilities (Chapter 6.3). Reproduces the odds ratios in Table 6.2 and
# maps them onto individual S->A probabilities. Associational, not causal.

suppressWarnings(source(file.path("R", "utils.R")))

fit_risk_model <- function(cohort) {
  # Outcome: currently Addicted (I) vs not
  cohort$addicted <- as.integer(cohort$state == "I")
  glm(addicted ~ age + female + urban + low_ses + high_depression + high_literacy,
      data = cohort, family = binomial())
}

odds_ratios <- function(model) {
  est <- exp(cbind(OR = coef(model), confint.default(model)))
  round(est, 3)
}

# Map fitted coefficients to a heterogeneous S->A probability for a covariate
# profile (Eq 5.1 / 6.1).
predict_pSA <- function(model, newdata, base_pSA = 0.21) {
  lo <- qlogis(base_pSA) + predict(model, newdata) - coef(model)[1]
  plogis(lo)
}

if (sys.nframe() == 0) {
  message("Run with: source('R/03_synthetic_data.R'); d <- generate_cohort();",
          " m <- fit_risk_model(d); odds_ratios(m)")
}
