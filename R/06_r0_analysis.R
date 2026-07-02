# 06_r0_analysis.R --------------------------------------------------------
# Basic reproduction number, next-generation matrix, endemic equilibrium and
# threshold analysis (Chapter 6.5.3 + Appendix A). Reconciles the crude and
# NGM derivations to R0 ~ 2.5 (NOT 26.4: an earlier per-capita formula
# double-counted the contact number).

suppressWarnings(source(file.path("R", "utils.R")))

# Crude, well-mixed R0 (Eq in 6.5.3): beta * mean_contacts * mean_duration ----
R0_crude <- function(p = default_params()) {
  D <- 1 / (p$gamma2 + p$natmort)          # mean months in Addicted state
  p$beta * p$k * D
}

# Next-generation-matrix R0 (Appendix A.4) -----------------------------------
# Infected subsystem = {A, I}. F = new-infection (transmission) Jacobian,
# V = transition Jacobian, both at the disease-free equilibrium S0.
R0_ngm <- function(p = default_params(), S0 = 0.95) {
  # per-month rates
  b  <- p$beta * p$k
  sg <- p$sigma; g1 <- p$gamma1; g2 <- p$gamma2; d <- p$delta; mu <- p$natmort
  # Reduced 2x2 (A, I) system after eliminating fast R relapse loop:
  Fm <- matrix(c(b * p$eta * S0, b * S0,
                 0,               0), 2, 2, byrow = TRUE)
  Vm <- matrix(c(sg + g1 + mu, 0,
                 -sg,          g2 + mu), 2, 2, byrow = TRUE)
  K  <- Fm %*% solve(Vm)                    # next-generation matrix
  max(Re(eigen(K)$values))                  # spectral radius = R0
}

# Networked R0 on a graph with degree distribution (Appendix A.7) ------------
# R0_net = R0 * <k^2 - k> / <k>  (heterogeneity amplification).
R0_network <- function(p = default_params(), degrees) {
  k1 <- mean(degrees); k2 <- mean(degrees^2)
  base <- p$beta / (p$gamma2 + p$natmort)
  base * (k2 - k1) / k1
}

# Endemic equilibrium (Appendix A.5-A.6). Returns S*, A*, I*, R*. -------------
endemic_equilibrium <- function(p = default_params(), S0 = 1) {
  R0 <- R0_ngm(p, S0 = S0)
  if (R0 <= 1) return(c(S = S0, A = 0, I = 0, R = 0, R0 = R0))
  Sstar <- S0 / R0
  # Solve remaining balances numerically for the positive fixed point.
  f <- function(A) {
    I <- (p$sigma / (p$gamma2 + p$natmort)) * A
    R <- (p$gamma2 * I) / (p$delta + p$natmort + 1e-9)
    # At-Risk balance: inflow from S and relapse = outflow
    inflow  <- p$beta * p$k * Sstar * (I + p$eta * A) + p$delta * R
    outflow <- (p$sigma + p$gamma1 + p$natmort) * A
    inflow - outflow
  }
  Astar <- tryCatch(uniroot(f, c(1e-6, 1))$root, error = function(e) NA_real_)
  Istar <- (p$sigma / (p$gamma2 + p$natmort)) * Astar
  Rstar <- (p$gamma2 * Istar) / (p$delta + p$natmort + 1e-9)
  c(S = Sstar, A = Astar, I = Istar, R = Rstar, R0 = R0)
}

if (sys.nframe() == 0) {
  p <- default_params()
  cat(sprintf("Crude R0  = %.1f  (well-mixed, upper bound)\n", R0_crude(p)))
  cat(sprintf("NGM   R0  = %.2f  (Appendix A.4)\n", R0_ngm(p)))
  print(round(endemic_equilibrium(p), 4))
}
