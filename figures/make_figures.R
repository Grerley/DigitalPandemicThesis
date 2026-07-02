# make_figures.R ----------------------------------------------------------
# Regenerates the core thesis figures as black-and-white, journal-style plots.
# Requires ggplot2. Figures are written to figures/output/.

suppressWarnings(source(file.path("R", "utils.R")))
suppressWarnings(source(file.path("R", "01_sair_model.R")))
suppressWarnings(source(file.path("R", "07_interventions.R")))

make_all_figures <- function(outdir = file.path("figures", "output")) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 required."); library(ggplot2)
  ensure_dir(outdir)

  # Fig 7.1 baseline compartments -------------------------------------------
  base <- simulate_sair(horizon_years = 20)
  long <- data.frame(
    year = rep(base$year, 4),
    prop = c(base$S, base$A, base$I, base$R),
    state = rep(c("S Susceptible","A At-Risk","I Addicted","R Recovered"),
                each = nrow(base))
  )
  p1 <- ggplot(long, aes(year, prop*100, linetype = state)) +
    geom_line(colour = "black", linewidth = 0.6) +
    scale_linetype_manual(values = c("solid","longdash","dotdash","dotted")) +
    labs(x = "Year", y = "Prevalence (%)", linetype = NULL,
         title = "Figure 7.1  Baseline S-A-I-R trajectories") +
    theme_bw_journal()
  ggsave(file.path(outdir, "fig_7_1_baseline.png"), p1, width = 7, height = 4, dpi = 200)

  # Fig 8.x intervention comparison -----------------------------------------
  sc <- run_scenarios()
  p2 <- ggplot(sc, aes(reorder(scenario, -peak_prevalence), peak_prevalence*100)) +
    geom_col(fill = "grey20", width = 0.65) +
    geom_text(aes(label = sprintf("%.1f%%", peak_prevalence*100)),
              vjust = -0.4, size = 3) +
    labs(x = NULL, y = "Peak addicted prevalence (%)",
         title = "Intervention scenarios (conditional on assumed efficacies)") +
    theme_bw_journal() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))
  ggsave(file.path(outdir, "fig_8_interventions.png"), p2, width = 7, height = 4, dpi = 200)

  invisible(TRUE)
}

if (sys.nframe() == 0 && requireNamespace("ggplot2", quietly = TRUE)) make_all_figures()
