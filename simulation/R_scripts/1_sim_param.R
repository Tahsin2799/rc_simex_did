# Simulation parameters, design v2.0 (2026-09-12).
#   source("R_scripts/1_sim_param.R")
#
# p = 20 layered covariates, moments calibrated to BIHS Round 1 (2011).
#
# CHANGE FROM v1.3: the outcome equation now has covariate-dependent TRENDS
# (delta). In v1.3 covariates entered only the level, so they cancelled in
# Delta Y, unconditional parallel trends held, OLS-DiD without covariates was
# unbiased, and covariate measurement error could not bias any Delta-Y
# estimator. With delta != 0, conditional parallel trends holds but
# unconditional does not, so covariate adjustment matters and covariate
# error has a bias channel (PreReg v2.0 Proposition 2).

default_params <- function(panel_strength = c("low", "high")) {
  panel_strength <- rlang::arg_match(panel_strength)

  beta <- list(
    b1 =  0.30, b2 =  0.10, b3 =  0.40, b4 = -0.20, b5 =  0.30,
    b6 =  0.15, b7 = -0.10, b8 =  0.10,
    b9_d1 = 0.10, b9_d2 = 0.20, b9_d3 = 0.05, b9_d4 = -0.10, b9_d5 = -0.15, b9_d6 = -0.05,
    b10_d1 = 0.20, b10_d2 = 0.40
  )
  # Trend loadings: covariates that predict levels also predict growth, at
  # half strength. Calibrated so OLS-DiD (no covariates) carries a bias of
  # roughly 0.5-0.8 theta; see smoke_v2.R output recorded in PreReg v2 s10.
  delta <- lapply(beta, \(b) 0.5 * b)
  names(delta) <- sub("^b", "d", names(beta))

  list(
    cov_moments = list(
      c1_mean = 28.015,  c1_sd = 16.973,
      c2_p    = 0.500,
      c3_mean = 58.249,  c3_sd = 18.577,
      c4_mean = 4.950,   c4_sd = 1.824,
      c5_p    = 0.79,
      c6_mean = 38.241,  c6_sd = 12.861,
      c7_p    = 0.157,
      c8_mean = 0.119,   c8_sd = 0.201,
      c9_probs = c(0.0683, 0.1915, 0.3177, 0.0699, 0.0953, 0.0966, 0.1607),
      c10_probs = c(0.1978, 0.3006, 0.5016),
      c11_p   = 0.746,
      c12_mean = 0.619,  c12_sd = 0.189,
      c13_p   = 0.521,
      c14_mean = 2.268,  c14_sd = 1.236,
      c15_mean = 0.741,  c15_sd = 0.121,
      c16_mean = 14.568, c16_sd = 10.486,
      c17_mean = 0.539,  c17_sd = 0.095,
      c18_mean = 25.480, c18_sd = 0.425,
      c19_mean = 3864.78, c19_sd = 1852.87,
      c20_probs = c(0.907, 0.090, 0.003)
    ),
    kappa = list(
      a0 = -0.55,
      a1 = 0.20, a2 = 0.00, a3 = -0.35, a4 = 0.20, a5 = -0.40,
      a6 = -0.10, a7 = 0.30, a8 = 0.10,
      a9_d1 = 0.20, a9_d2 = 0.00, a9_d3 = 0.10, a9_d4 = 0.30, a9_d5 = 0.40, a9_d6 = 0.20,
      a10_d1 = -0.20, a10_d2 = -0.40
    ),
    lambda_0 = 5.00,
    lambda_t = 0.20,
    lambda_a = 0.00,
    theta    = 0.20,
    beta  = beta,
    delta = delta,
    gamma = list(g11 = 0.30, g12 = 0.20, g13 = -0.10, g14 = 0.10, g15 = 0.10),
    sigma_alpha = if (panel_strength == "low") 0.7 else 1.0,
    sigma_eps   = if (panel_strength == "low") 0.7 else 0.5,
    panel_strength = panel_strength
  )
}
