# 01_sair_model.R ---------------------------------------------------------
# Deterministic S-A-I-R model with relapse, as a nonlinear difference-equation
# system (Chapter 4). The population is described by the state distribution
# x(t) = [S, A, I, R]^T evolving under a STATE-DEPENDENT transition matrix
# P(x_t): the S->A probability is the mass-action force of infection (Eq 4.3),
# which makes the map nonlinear and produces the epidemic peak + endemic level.
#
# This is NOT a linear Markov chain: a fixed P would converge monotonically to
# a stationary distribution and could not generate a rise-peak-decline curve.

source(file.path("R", "utils.R"))

# Build the 4x4 state-dependent transition matrix P(x_t) (Eq 4.1) ------------
# Rows: from S, A, I, R. Columns: to S, A, I, R. Each row sums to 1.
transition_matrix <- function(state, p) {
  P_SA <- foi_SA(state, p)                 # susceptible -> at-risk (mass action)
  P_AI <- rate_to_prob(p$sigma * 12)       # at-risk -> addicted
  P_AS <- rate_to_prob(p$gamma1 * 12)      # at-risk -> susceptible (remission)
  P_IR <- rate_to_prob(p$gamma2 * 12)      # addicted -> recovered
  P_RA <- rate_to_prob(p$delta * 12)       # recovered -> at-risk (relapse)

  P <- matrix(0, 4, 4, dimnames = list(c("S","A","I","R"), c("S","A","I","R")))
  P["S","S"] <- 1 - P_SA;            P["S","A"] <- P_SA
  P["A","S"] <- P_AS;  P["A","I"] <- P_AI;  P["A","A"] <- 1 - P_AS - P_AI
  P["I","R"] <- P_IR;                P["I","I"] <- 1 - P_IR
  P["R","A"] <- P_RA;                P["R","R"] <- 1 - P_RA
  P
}

# Simulate the deterministic trajectory -------------------------------------
# horizon_years: projection length; dt is 1 month. Returns a data.frame of
# monthly state proportions plus incidence.
simulate_sair <- function(p = default_params(),
                          init = default_init(),
                          horizon_years = 20) {
  steps <- horizon_years * 12
  x <- renorm(init)
  out <- matrix(NA_real_, nrow = steps + 1, ncol = 5,
                dimnames = list(NULL, c("month","S","A","I","R")))
  out[1, ] <- c(0, x)
  incidence <- numeric(steps + 1)   # new S->A entries per step

  for (t in seq_len(steps)) {
    P <- transition_matrix(x, p)
    incidence[t + 1] <- x["S"] * P["S","A"]
    x <- renorm(as.numeric(x %*% P))          # x_{t+1} = P(x_t)^T x_t (Eq 4.2)
    names(x) <- c("S","A","I","R")
    out[t + 1, ] <- c(t, x)
  }
  df <- as.data.frame(out)
  df$year <- df$month / 12
  df$incidence <- incidence
  df
}

# Convenience: summarise peak and endemic behaviour -------------------------
summarise_run <- function(df) {
  peak_i <- which.max(df$I)
  list(
    peak_prevalence   = df$I[peak_i],
    peak_year         = df$year[peak_i],
    endemic_prevalence = tail(df$I, 1),
    peak_incidence    = max(df$incidence),
    peak_incidence_year = df$year[which.max(df$incidence)]
  )
}

if (sys.nframe() == 0) {
  df <- simulate_sair()
  print(summarise_run(df))     # expect peak ~0.13 @ ~yr 8, endemic ~0.08
}
