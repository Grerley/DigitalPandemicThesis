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
  # Compartment ratios at the fixed point, from the I- and R-balance equations:
  #   sigma*A = (gamma2 + natmort)*I   ->  I = kI * A
  #   gamma2*I = (delta + natmort)*R   ->  R = kR * A
  kI <- p$sigma / (p$gamma2 + p$natmort)
  kR <- (p$gamma2 * kI) / (p$delta + p$natmort + 1e-9)
  # The At-Risk balance alone is homogeneous in A (inflow and outflow both scale
  # linearly with A), so it cannot pin down A*'s magnitude. Close the system with
  # the conservation constraint S* + A* + I* + R* = 1.
  Astar <- (1 - Sstar) / (1 + kI + kR)
  Istar <- kI * Astar
  Rstar <- kR * Astar
  c(S = Sstar, A = Astar, I = Istar, R = Rstar, R0 = R0)
}

if (sys.nframe() == 0) {
  p <- default_params()
  cat(sprintf("Crude R0  = %.1f  (well-mixed, upper bound)\n", R0_crude(p)))
  cat(sprintf("NGM   R0  = %.2f  (Appendix A.4)\n", R0_ngm(p)))
  print(round(endemic_equilibrium(p), 4))
}
