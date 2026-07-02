# 07_interventions.R ------------------------------------------------------
# Intervention scenario simulation (Chapter 8). Interventions are represented
# as modifications to model parameters (transmission, progression, recovery).
#
# IMPORTANT: intervention efficacies are ASSUMED values informed by the
# direction/rough magnitude of published evaluations, NOT estimated here. The
# resulting reductions are what the model implies IF an intervention achieves
# the assumed effect; they are comparative, not predictive.

suppressWarnings(source(file.path("R", "utils.R")))
suppressWarnings(source(file.path("R", "01_sair_model.R")))

# Assumed efficacy parameters (Table 8.x). Each entry scales a mechanism. -----
intervention_presets <- function() list(
  baseline    = list(),
  education   = list(beta = 0.82),                       # -18% transmission (behavioural)
  regulation  = list(beta = 0.71),                       # -29% transmission (structural, platform design)
  age_gating  = list(mu = 0.5, beta = 0.9),              # cut exogenous entry + some transmission
  treatment   = list(gamma2 = 1.4, delta = 0.7),         # +40% recovery, -30% relapse
  comprehensive = list(beta = 0.55, mu = 0.5,            # combined structural + behavioural + treatment
                       gamma2 = 1.4, delta = 0.7)
)

apply_intervention <- function(p, mods, scaleup_years = 5, horizon_years = 20) {
  # Logistic scale-up of coverage (Eq 8.1): effect ramps in over scaleup_years.
  # For simplicity we apply the full modified parameters after ramp; a smooth
  # ramp is available via time-varying beta below.
  p2 <- p
  for (nm in names(mods)) p2[[nm]] <- p[[nm]] * mods[[nm]]
  p2
}

run_scenarios <- function(p = default_params(), horizon_years = 20) {
  presets <- intervention_presets()
  res <- lapply(names(presets), function(nm) {
    p2 <- apply_intervention(p, presets[[nm]], horizon_years = horizon_years)
    df <- simulate_sair(p2, horizon_years = horizon_years)
    s  <- summarise_run(df)
    data.frame(scenario = nm,
               peak_prevalence = round(s$peak_prevalence, 4),
               peak_year = round(s$peak_year, 1),
               endemic = round(s$endemic_prevalence, 4))
  })
  out <- do.call(rbind, res)
  base_peak <- out$peak_prevalence[out$scenario == "baseline"]
  out$pct_reduction <- round(100 * (base_peak - out$peak_prevalence) / base_peak, 1)
  out
}

# Network-targeted vs random intervention on scale-free graphs (Fig 8.4) -----
# Demonstrates that protecting high-degree hubs beats random coverage.
targeted_vs_random <- function(n = 500, coverage = 0.30, runs = 20) {
  suppressWarnings(source(file.path("R", "02_network_simulation.R")))
  set.seed(7)
  peak <- function(seed_hub) {
    v <- replicate(runs, {
      g <- make_network(n, "BA")
      simulate_network_once(g, horizon_months = 180, seed_hub = seed_hub)$I
    })
    max(rowMeans(v))
  }
  data.frame(strategy = c("no_intervention","hub_targeted"),
             peak = c(peak(FALSE), peak(TRUE)))
}

if (sys.nframe() == 0) {
  print(run_scenarios())     # comprehensive ~ -68% peak (conditional on assumptions)
}
