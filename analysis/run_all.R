# run_all.R ---------------------------------------------------------------
# Master pipeline: reproduces the core results and figures of the thesis.
# Run from the repository root:  Rscript analysis/run_all.R

options(warn = 1)
source(file.path("R", "utils.R"))
source(file.path("R", "01_sair_model.R"))
source(file.path("R", "06_r0_analysis.R"))
source(file.path("R", "07_interventions.R"))
ensure_dir("outputs")

cat("== 1. Reproduction number ==\n")
p <- default_params()
cat(sprintf("  Crude (well-mixed) R0 = %.1f\n", R0_crude(p)))
cat(sprintf("  NGM R0                = %.2f\n", R0_ngm(p)))
eq <- endemic_equilibrium(p); print(round(eq, 4))

cat("\n== 2. Baseline S-A-I-R trajectory ==\n")
base <- simulate_sair(p, horizon_years = 20)
print(summarise_run(base))
write.csv(base, "outputs/baseline_trajectory.csv", row.names = FALSE)

cat("\n== 3. Intervention scenarios (conditional on assumed efficacies) ==\n")
sc <- run_scenarios(p); print(sc)
write.csv(sc, "outputs/intervention_scenarios.csv", row.names = FALSE)

cat("\n== 4. Network topology comparison (small MC for speed) ==\n")
if (requireNamespace("igraph", quietly = TRUE)) {
  source(file.path("R", "02_network_simulation.R"))
  for (tp in c("ER","WS","BA")) {
    e <- simulate_network_ensemble(tp, n = 400, runs = 10, horizon_months = 150)
    cat(sprintf("  %s: peak mean I = %.3f at year %.1f\n",
                tp, max(e$mean), e$year[which.max(e$mean)]))
  }
} else cat("  (install 'igraph' to run network simulations)\n")

cat("\n== 5. Figures ==\n")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  source(file.path("figures", "make_figures.R")); make_all_figures()
  cat("  figures written to figures/output/\n")
} else cat("  (install 'ggplot2' to render figures)\n")

cat("\nDone. See outputs/ and figures/output/.\n")
