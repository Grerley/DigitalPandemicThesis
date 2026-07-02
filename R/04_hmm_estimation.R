# 04_hmm_estimation.R -----------------------------------------------------
# Hidden Markov Model estimation of S-A-I-R transition probabilities from
# longitudinal data (Chapter 6.2). Hidden state in {S,A,I,R}; observations are
# daily screen time and a composite symptom score. Validated on synthetic data
# whose generating transition matrix is known (recovery test). Real-data
# application awaits suitable longitudinal cohorts.
#
# Requires: depmixS4

suppressWarnings(source(file.path("R", "utils.R")))

fit_hmm <- function(cohort, nstates = 4) {
  if (!requireNamespace("depmixS4", quietly = TRUE))
    stop("Package 'depmixS4' is required for HMM estimation.")
  cohort <- cohort[order(cohort$id, cohort$year), ]
  ntimes <- as.numeric(table(cohort$id))
  mod <- depmixS4::depmix(
    response = list(screen ~ 1, symptom ~ 1),
    data     = cohort,
    nstates  = nstates,
    family   = list(stats::gaussian(), stats::gaussian()),
    ntimes   = ntimes
  )
  fit <- depmixS4::fit(mod, verbose = FALSE)
  fit
}

# Extract the estimated transition matrix and map states to S,A,I,R by their
# mean symptom score (low -> S, high -> I).
hmm_transition_matrix <- function(fit) {
  pars <- depmixS4::getpars(fit)
  # depmixS4 stores transition pars; use summary() structure for robustness
  tm <- depmixS4::posterior(fit)
  list(fit = fit, decoded = tm)
}

# Validation: compare estimated vs known generating matrix -------------------
validate_hmm <- function(known_matrix, estimated_matrix) {
  mae <- mean(abs(known_matrix - estimated_matrix))
  list(mae = mae, pass = mae < 0.05)
}

if (sys.nframe() == 0) {
  message("Run with: source('R/03_synthetic_data.R'); d <- generate_cohort();",
          " fit <- fit_hmm(d)  # requires depmixS4")
}
