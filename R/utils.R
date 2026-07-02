# utils.R -----------------------------------------------------------------
# Shared helpers and default (calibrated) parameters for the S-A-I-R model of
# youth digital addiction. All rates are expressed per MONTH unless noted.
#
# Thesis: "The Digital Pandemic" (Mutibura, 2026). See docs/model_specification.md.

# Default calibrated parameter set (Chapter 6, Table 6.4) --------------------
# Rates per MONTH. This set is calibrated so the model reproduces the thesis
# target dynamics: R0 (NGM) ~ 2.4, peak addicted prevalence ~ 14% around year 7,
# and endemic prevalence ~ 8%. `beta` is the EFFECTIVE per-contact transmission
# probability (so beta * k is the per-month transmission scale); an earlier
# per-capita mis-specification that produced R0 ~ 26 is not used (see
# docs/model_specification.md, "Basic reproduction number").
default_params <- function() {
  list(
    beta   = 0.033,   # per-contact transmission probability per month (matches thesis)
    k      = 3.7,     # effective mean-field transmitting contacts (beta*k = 0.122/mo).
                      # NB: the NETWORK model (02) uses the full friendship degree (~10)
                      # in the Reed-Frost kernel; k here is the smaller mean-field
                      # effective contact number, since homogeneous mixing over the full
                      # degree overstates spread (an earlier such mis-use gave R0 ~ 26).
    eta    = 0.5,     # reduced infectivity weight of At-Risk vs Addicted, in [0,1]
    mu     = 0.0015,  # exogenous (non-social) acquisition rate per month
    sigma  = 0.02,    # At-Risk -> Addicted progression per month (0.24/yr)
    gamma1 = 0.018,   # At-Risk -> Susceptible remission per month (0.22/yr)
    gamma2 = 0.045,   # Addicted -> Recovered recovery per month (0.54/yr; mean ~1.9 yr)
    delta  = 0.006,   # Recovered -> At-Risk relapse per month (0.07/yr)
    Lambda = 0.0,     # inflow (set >0 for open population / demography)
    natmort = 0.0     # background exit; 0 for closed cohort over the horizon
  )
}

# Initial state distribution (Chapter 5 baseline): S, A, I, R -----------------
default_init <- function() c(S = 0.823, A = 0.065, I = 0.047, R = 0.065)

# Mass-action force of infection -> S->A transition probability (Eq 4.3) ------
foi_SA <- function(state, p) {
  1 - exp(-p$beta * p$k * (state["I"] + p$eta * state["A"])) + p$mu
}

# Normalise a state vector to sum to 1 (guards against numerical drift) -------
renorm <- function(x) x / sum(x)

# Convert an annual rate to a monthly transition probability (Eq 6.2) ---------
rate_to_prob <- function(rate_per_year, dt_months = 1) {
  1 - exp(-(rate_per_year / 12) * dt_months)
}

# Simple palette-free (black & white) ggplot theme for journal figures -------
theme_bw_journal <- function(base_size = 11) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(NULL)
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey88", linewidth = 0.3),
      axis.line  = ggplot2::element_line(colour = "black", linewidth = 0.4),
      legend.position = "right",
      text = ggplot2::element_text(colour = "black")
    )
}

ensure_dir <- function(path) if (!dir.exists(path)) dir.create(path, recursive = TRUE)
