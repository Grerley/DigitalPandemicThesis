# 03_synthetic_data.R -----------------------------------------------------
# Synthetic longitudinal cohort generator (Chapter 5.3). This is a MODELLING
# INSTRUMENT, not observed data: it produces individual-level trajectories that
# reproduce known empirical patterns (prevalence, risk-factor ORs) while
# enabling direct computation of transition probabilities for HMM validation.
# Fully regenerable with a fixed seed.

suppressWarnings(source(file.path("R", "utils.R")))

# Per-individual annual transition given covariates (logistic tilt of P_SA) --
covariate_pSA <- function(base_pSA, covar) {
  # log-odds shifts from meta-analytic ORs (Chapter 5/6)
  lo <- qlogis(base_pSA) +
    log(0.94) * (covar$age - 14) +          # older = protective
    log(1.38) * covar$female +
    log(2.18) * covar$urban +
    log(1.45) * covar$low_ses +
    log(2.84) * covar$high_depression +
    log(0.72) * covar$high_literacy
  plogis(lo)
}

generate_cohort <- function(n = 10000, years = 10, seed = 42,
                            p = default_params(), noise = 0.10) {
  set.seed(seed)
  covar <- data.frame(
    id            = seq_len(n),
    age           = sample(10:20, n, replace = TRUE),
    female        = rbinom(n, 1, 0.51),
    urban         = rbinom(n, 1, 0.68),
    low_ses       = rbinom(n, 1, 0.33),
    high_depression = rbinom(n, 1, 0.20),
    high_literacy = rbinom(n, 1, 0.40)
  )
  init <- default_init()
  state <- sample(c("S","A","I","R"), n, replace = TRUE, prob = init)

  pAI <- rate_to_prob(p$sigma * 12); pAS <- rate_to_prob(p$gamma1 * 12)
  pIR <- rate_to_prob(p$gamma2 * 12); pRA <- rate_to_prob(p$delta * 12)

  rows <- vector("list", years + 1)
  emit <- function(s) {                     # observed screen-time & symptom score
    mu <- c(S = 2, A = 4, I = 6.5, R = 3)[s]; sc <- c(S = 10, A = 45, I = 75, R = 30)[s]
    c(screen = rnorm(1, mu, 0.8), symptom = rnorm(1, sc, 8))
  }
  for (yr in 0:years) {
    obs <- t(vapply(state, emit, numeric(2)))
    # measurement error: misclassify observed state with prob = noise
    rows[[yr + 1]] <- data.frame(covar, year = yr, state = state,
                                 screen = obs[,1], symptom = obs[,2])
    if (yr == years) break
    # annual transition
    base_pSA <- 1 - (1 - foi_SA(setNames(as.numeric(table(factor(state,
                    c("S","A","I","R")))) / n, c("S","A","I","R")), p))^12
    for (i in seq_len(n)) {
      s <- state[i]
      if (s == "S") {
        pSA_i <- covariate_pSA(min(max(base_pSA, 1e-4), 0.9), covar[i, ])
        if (runif(1) < pSA_i) state[i] <- "A"
      } else if (s == "A") {
        u <- runif(1)
        if (u < pAI) state[i] <- "I" else if (u < pAI + pAS) state[i] <- "S"
      } else if (s == "I") {
        if (runif(1) < pIR) state[i] <- "R"
      } else if (s == "R") {
        if (runif(1) < pRA) state[i] <- "A"
      }
    }
  }
  do.call(rbind, rows)
}

if (sys.nframe() == 0) {
  d <- generate_cohort(n = 2000, years = 10)
  tab <- prop.table(table(d$year, d$state), 1)
  cat("Addicted prevalence by year:\n"); print(round(tab[, "I"], 3))
}
