# 05_regression.R ---------------------------------------------------------
# Logistic regression for risk factors and heterogeneous transition
# probabilities (Chapter 6.3). Two complementary models:
#
#   * fit_transition_model() estimates the S->A *transition* odds ratios. This
#     is the validation target for the synthetic generator: it recovers the ORs
#     that 03_synthetic_data.R applies to the S->A hazard (Table 6.2).
#   * fit_risk_model() is an associational cross-sectional model of *current*
#     addiction (I) prevalence. It does NOT recover the transition ORs (a
#     cross-sectional snapshot confounds incidence with duration) and is kept
#     for descriptive use only.
#
# Associational, not causal.

suppressWarnings(source(file.path("R", "utils.R")))

# Build the person-year S->A transition frame from a longitudinal cohort. -----
# One row per person-year spent Susceptible, with outcome `transitioned` = 1 if
# the individual is At-Risk (A) at the next observation, else 0. Covariates are
# the time-invariant risk factors. Fully vectorised over the ordered cohort.
build_transition_frame <- function(cohort) {
  covars <- c("age","female","urban","low_ses","high_depression","high_literacy")
  cohort <- cohort[order(cohort$id, cohort$year), ]
  n <- nrow(cohort)
  next_state <- c(cohort$state[-1], NA)          # state at the next observation
  same_id    <- c(cohort$id[-1], NA) == cohort$id # next row is the same person
  keep <- !is.na(same_id) & same_id & cohort$state == "S"
  data.frame(cohort[keep, covars, drop = FALSE],
             year = cohort$year[keep],
             transitioned = as.integer(next_state[keep] == "A"))
}

# Estimate S->A transition odds ratios (validation target, Table 6.2). --------
# `factor(year)` absorbs the year-varying baseline hazard (the force of
# infection changes as the epidemic evolves and the susceptible pool depletes);
# without it the pooled slopes are attenuated toward 1. The year terms are
# nuisance baseline-hazard parameters -- the risk-factor coefficients are the
# quantities of interest.
fit_transition_model <- function(cohort) {
  tf <- build_transition_frame(cohort)
  glm(transitioned ~ age + female + urban + low_ses + high_depression +
        high_literacy + factor(year),
      data = tf, family = binomial())
}

# Convenience: just the risk-factor odds ratios (drops the year nuisance terms).
transition_odds_ratios <- function(model) {
  or <- odds_ratios(model)
  or[!grepl("^factor\\(year\\)|Intercept", rownames(or)), , drop = FALSE]
}

# Cross-sectional prevalence model (associational; NOT the transition ORs). ----
fit_risk_model <- function(cohort) {
  # Outcome: currently Addicted (I) vs not.
  cohort$addicted <- as.integer(cohort$state == "I")
  glm(addicted ~ age + female + urban + low_ses + high_depression + high_literacy,
      data = cohort, family = binomial())
}

odds_ratios <- function(model) {
  est <- exp(cbind(OR = coef(model), confint.default(model)))
  round(est, 3)
}

# Map fitted transition coefficients to a heterogeneous S->A probability for a
# covariate profile (Eq 5.1 / 6.1). Applies only the risk-factor log-odds to a
# chosen baseline (dropping the intercept and any year nuisance terms), so it
# works with either model and needs no `year` column. Pass the transition model
# for calibrated ORs.
predict_pSA <- function(model, newdata, base_pSA = 0.21) {
  cf   <- coef(model)
  risk <- intersect(c("age","female","urban","low_ses",
                       "high_depression","high_literacy"), names(cf))
  contrib <- as.matrix(newdata[, risk, drop = FALSE]) %*% cf[risk]
  plogis(qlogis(base_pSA) + as.numeric(contrib))
}

if (sys.nframe() == 0) {
  message("Run with: source('R/03_synthetic_data.R'); d <- generate_cohort();",
          " m <- fit_transition_model(d); odds_ratios(m)")
}
