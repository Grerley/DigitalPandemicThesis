# 02_network_simulation.R -------------------------------------------------
# Stochastic S-A-I-R contagion on explicit social networks (Chapter 4 / 7.4).
# Compares Erdos-Renyi (random), Watts-Strogatz (small-world) and
# Barabasi-Albert (scale-free) topologies via Monte Carlo. Scale-free networks
# (influencer hubs) produce higher, earlier peaks and near-zero epidemic
# threshold; targeting high-degree nodes is disproportionately effective.
#
# Requires: igraph

suppressWarnings(source(file.path("R", "utils.R")))

make_network <- function(n = 1000, type = c("BA","ER","WS"), mean_degree = 10) {
  type <- match.arg(type)
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required for network simulation.")
  switch(type,
    ER = igraph::sample_gnp(n, p = mean_degree / (n - 1)),
    WS = igraph::sample_smallworld(1, n, nei = mean_degree %/% 2, p = 0.1),
    BA = igraph::sample_pa(n, power = 1, m = mean_degree %/% 2, directed = FALSE)
  )
}

# One stochastic run on a fixed graph. States: 1=S 2=A 3=I 4=R ---------------
simulate_network_once <- function(g, p = default_params(),
                                   horizon_months = 240, seed_frac = 0.01,
                                   seed_hub = FALSE) {
  if (!requireNamespace("igraph", quietly = TRUE))
    stop("Package 'igraph' is required.")
  n   <- igraph::vcount(g)
  adj <- igraph::as_adj_list(g)
  state <- rep(1L, n)                                   # all Susceptible
  n_seed <- max(1, round(seed_frac * n))
  seeds <- if (seed_hub) order(igraph::degree(g), decreasing = TRUE)[seq_len(n_seed)]
           else sample.int(n, n_seed)
  state[seeds] <- 3L                                    # seed Addicted

  pAI <- rate_to_prob(p$sigma * 12); pAS <- rate_to_prob(p$gamma1 * 12)
  pIR <- rate_to_prob(p$gamma2 * 12); pRA <- rate_to_prob(p$delta * 12)
  prev <- numeric(horizon_months + 1)
  prev[1] <- mean(state == 3L)

  for (t in seq_len(horizon_months)) {
    new <- state
    for (i in seq_len(n)) {
      s <- state[i]
      if (s == 1L) {                                    # Susceptible
        nb <- adj[[i]]
        m_inf <- sum(state[nb] == 3L) + p$eta * sum(state[nb] == 2L)
        pinf <- 1 - (1 - p$beta)^m_inf + p$mu           # Reed-Frost (Eq 4.4)
        if (runif(1) < pinf) new[i] <- 2L
      } else if (s == 2L) {                             # At-Risk
        u <- runif(1)
        if (u < pAI) new[i] <- 3L else if (u < pAI + pAS) new[i] <- 1L
      } else if (s == 3L) {                             # Addicted
        if (runif(1) < pIR) new[i] <- 4L
      } else {                                          # Recovered
        if (runif(1) < pRA) new[i] <- 2L
      }
    }
    state <- new
    prev[t + 1] <- mean(state == 3L)
  }
  data.frame(month = 0:horizon_months, year = (0:horizon_months)/12, I = prev)
}

# Monte Carlo ensemble over many graph draws / runs --------------------------
simulate_network_ensemble <- function(type = "BA", n = 1000, runs = 100,
                                       p = default_params(), horizon_months = 240,
                                       seed_hub = FALSE, seed = 1) {
  set.seed(seed)
  mat <- replicate(runs, {
    g <- make_network(n, type)
    simulate_network_once(g, p, horizon_months, seed_hub = seed_hub)$I
  })
  months <- 0:horizon_months
  data.frame(
    month = months, year = months / 12,
    mean  = rowMeans(mat),
    lo    = apply(mat, 1, quantile, 0.025),
    hi    = apply(mat, 1, quantile, 0.975)
  )
}

if (sys.nframe() == 0) {
  for (tp in c("ER","WS","BA")) {
    e <- simulate_network_ensemble(tp, n = 500, runs = 20, horizon_months = 180)
    cat(sprintf("%s: peak mean I = %.3f at year %.1f\n",
                tp, max(e$mean), e$year[which.max(e$mean)]))
  }
}
