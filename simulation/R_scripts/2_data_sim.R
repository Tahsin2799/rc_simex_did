# Data-generating process, design v2.0 (2026-09-12).
#   source("R_scripts/2_data_sim.R")
#
# Layer A c1-c5 confounders; Layer B c6-c10 confounders with non-linear
# loading; Layer C c11-c15 effect modifiers (nonlinear family only);
# Layer D c16-c20 pure noise.
#
# Outcome (v2.0):
#   Y_i0 = lambda_0 + lambda_a A + beta'X + nl(X) + alpha_i + eps_i0
#   Y_i1 = Y_i0 - eps_i0 + eps_i1 + lambda_t + delta'X + nlt(X) + theta_i A
# so   dY_i = lambda_t + delta'X + nlt(X) + theta_i A + (eps_i1 - eps_i0).
# Conditional PT holds (given X); unconditional PT fails because A depends
# on X and delta != 0.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})
source("R_scripts/1_sim_param.R")

COV_NAMES <- paste0("c", 1:20)
COV_FORMULA_STR <- paste(COV_NAMES, collapse = " + ")
# Continuous confounders that receive covariate measurement error.
ME_X_COVS <- c("c1", "c3", "c4", "c6", "c8")

.dummies <- function(x, levels, prefix) {
  rest <- levels[-1]
  out  <- lapply(rest, \(lv) as.integer(x == lv))
  names(out) <- paste0(prefix, "_d", seq_along(rest))
  as.data.frame(out)
}

# Numeric design matrix from a data frame holding c1..c20: 27 columns, less
# any column with no variation in this sample. c20 level 2 has prevalence
# 0.003, so its dummy is empty in about 3.5% of N = 1000 draws; lm and glm
# alias such a column, but DRDID solves its own normal equations and dies on
# the Cholesky (std_ipw_did_panel). Dropping a constant column is an exact
# reparameterisation, so no estimate changes, and it is the predicate
# .dict_transform (30_core.R) already applies per fold. PreReg v2.0 s5.
design_matrix <- function(dat) {
  num <- \(x) if (is.factor(x)) as.integer(as.character(x)) else as.numeric(x)
  M <- cbind(
    c1 = num(dat$c1), c2 = num(dat$c2), c3 = num(dat$c3), c4 = num(dat$c4),
    c5 = num(dat$c5), c6 = num(dat$c6), c7 = num(dat$c7), c8 = num(dat$c8),
    .dummies(num(dat$c9), 0:6, "c9"),
    .dummies(num(dat$c10), c(1, 2, 4), "c10"),
    c11 = num(dat$c11), c12 = num(dat$c12), c13 = num(dat$c13),
    c14 = num(dat$c14), c15 = num(dat$c15), c16 = num(dat$c16),
    c17 = num(dat$c17), c18 = num(dat$c18), c19 = num(dat$c19),
    .dummies(num(dat$c20), 0:2, "c20")
  )
  M <- as.matrix(M)
  M[, apply(M, 2, \(v) any(v != v[1])), drop = FALSE]
}

# Wide one-row-per-unit frame used by every estimator:
#   id, A, y0, y1, dY, c1..c20 (as in the long frame).
wide_frame <- function(dat_long) {
  d0 <- dat_long |> filter(time == 0) |> select(id, A, all_of(COV_NAMES), y0 = Y)
  d1 <- dat_long |> filter(time == 1) |> select(id, y1 = Y)
  inner_join(d0, d1, by = "id") |> mutate(dY = y1 - y0)
}

# Standardised non-linear trend surface used by the nonlinear family and
# by the redesigned covariate-dependent outcome-ME wrapper. Returns a
# vector with mean ~0, sd ~1 over the population.
.h_nonlin <- function(z1, z3, z4, z6, c2, c7) {
  h <- 0.10 * z1 * z6 + 0.15 * sin(z3) + 0.10 * (z4^2 - 1) +
       0.10 * c2 * z3 + 0.08 * c7 * z6
  h / 0.19   # population sd of h under the calibrated moments (approx.)
}

data_gen_v2 <- function(N,
                        params     = default_params("low"),
                        dgp_family = c("linear", "nonlinear"),
                        seed       = 123,
                        return_cf  = FALSE) {
  dgp_family <- rlang::arg_match(dgp_family)
  set.seed(seed)
  m <- params$cov_moments; k <- params$kappa; b <- params$beta
  dl <- params$delta; g <- params$gamma

  c1  <- rnorm(N, m$c1_mean, m$c1_sd)
  c2  <- rbinom(N, 1, m$c2_p)
  c3  <- rnorm(N, m$c3_mean, m$c3_sd)
  c4  <- pmax(2L, as.integer(round(rnorm(N, m$c4_mean, m$c4_sd))))
  c5  <- rbinom(N, 1, m$c5_p)
  c6  <- rnorm(N, m$c6_mean, m$c6_sd)
  c7  <- rbinom(N, 1, m$c7_p)
  c8  <- rnorm(N, m$c8_mean, m$c8_sd)
  c9  <- sample(0:6, N, replace = TRUE, prob = m$c9_probs)
  c10 <- sample(c(1, 2, 4), N, replace = TRUE, prob = m$c10_probs)
  c11 <- rbinom(N, 1, m$c11_p)
  c12 <- pmin(1, pmax(0, rnorm(N, m$c12_mean, m$c12_sd)))
  c13 <- rbinom(N, 1, m$c13_p)
  c14 <- pmax(0L, as.integer(round(rnorm(N, m$c14_mean, m$c14_sd))))
  c15 <- pmin(1, pmax(0, rnorm(N, m$c15_mean, m$c15_sd)))
  c16 <- pmax(0, rnorm(N, m$c16_mean, m$c16_sd))
  c17 <- pmin(1, pmax(0, rnorm(N, m$c17_mean, m$c17_sd)))
  c18 <- rnorm(N, m$c18_mean, m$c18_sd)
  c19 <- pmax(0, rnorm(N, m$c19_mean, m$c19_sd))
  c20 <- sample(0:2, N, replace = TRUE, prob = m$c20_probs)

  z1 <- (c1 - m$c1_mean) / m$c1_sd;  z3 <- (c3 - m$c3_mean) / m$c3_sd
  z4 <- (c4 - m$c4_mean) / m$c4_sd;  z6 <- (c6 - m$c6_mean) / m$c6_sd
  z8 <- (c8 - m$c8_mean) / m$c8_sd
  z12 <- (c12 - m$c12_mean) / m$c12_sd; z14 <- (c14 - m$c14_mean) / m$c14_sd
  z15 <- (c15 - m$c15_mean) / m$c15_sd

  c9d  <- .dummies(c9, 0:6, "c9");  c10d <- .dummies(c10, c(1, 2, 4), "c10")

  lin <- function(coef, pre) {
    coef[[paste0(pre, "1")]] * z1 + coef[[paste0(pre, "2")]] * c2 +
    coef[[paste0(pre, "3")]] * z3 + coef[[paste0(pre, "4")]] * z4 +
    coef[[paste0(pre, "5")]] * c5 + coef[[paste0(pre, "6")]] * z6 +
    coef[[paste0(pre, "7")]] * c7 + coef[[paste0(pre, "8")]] * z8 +
    coef[[paste0(pre, "9_d1")]] * c9d$c9_d1 + coef[[paste0(pre, "9_d2")]] * c9d$c9_d2 +
    coef[[paste0(pre, "9_d3")]] * c9d$c9_d3 + coef[[paste0(pre, "9_d4")]] * c9d$c9_d4 +
    coef[[paste0(pre, "9_d5")]] * c9d$c9_d5 + coef[[paste0(pre, "9_d6")]] * c9d$c9_d6 +
    coef[[paste0(pre, "10_d1")]] * c10d$c10_d1 + coef[[paste0(pre, "10_d2")]] * c10d$c10_d2
  }

  e_X <- plogis(lin(k, "a") + k$a0)
  A   <- rbinom(N, 1, prob = e_X)

  lvl_X <- lin(b, "b")
  trd_X <- lin(dl, "d")

  if (dgp_family == "nonlinear") {
    h <- .h_nonlin(z1, z3, z4, z6, c2, c7)
    nl_lvl <- 0.19 * h                 # same surface as v1.3 in the level
    nl_trd <- 0.10 * h                 # non-linear trend heterogeneity
    z11 <- (c11 - m$c11_p) / sqrt(m$c11_p * (1 - m$c11_p))
    z13 <- (c13 - m$c13_p) / sqrt(m$c13_p * (1 - m$c13_p))
    theta_i <- params$theta * (1 + g$g11 * z11 + g$g12 * z12 + g$g13 * z13 +
                               g$g14 * z14 + g$g15 * z15)
  } else {
    nl_lvl <- 0; nl_trd <- 0
    theta_i <- rep(params$theta, N)
  }

  alpha_i <- rnorm(N, 0, params$sigma_alpha)
  eps_0   <- rnorm(N, 0, params$sigma_eps)
  eps_1   <- rnorm(N, 0, params$sigma_eps)

  base_0 <- params$lambda_0 + params$lambda_a * A + lvl_X + nl_lvl + alpha_i
  trend  <- params$lambda_t + trd_X + nl_trd
  Y_0 <- base_0 + eps_0
  Y_1 <- base_0 + trend + theta_i * A + eps_1

  dat_wide <- tibble(
    id = seq_len(N), A,
    c1, c2 = factor(c2), c3, c4 = as.integer(c4), c5 = factor(c5),
    c6, c7 = factor(c7), c8, c9 = factor(c9), c10 = factor(c10),
    c11 = factor(c11), c12, c13 = factor(c13), c14 = as.integer(c14), c15,
    c16, c17, c18, c19, c20 = factor(c20),
    Y_0, Y_1
  )
  dat_long <- dat_wide |>
    pivot_longer(cols = c(Y_0, Y_1), names_to = "time", values_to = "Y") |>
    mutate(time = as.integer(str_remove(time, "Y_"))) |>
    arrange(id, time)

  if (!return_cf) return(dat_long)

  # Realised ATET: mean of theta_i over treated units (sampling-error-free
  # counterfactual contrast; eps cancels exactly in the DiD contrast).
  true_atet <- mean(theta_i[A == 1])
  # Bias of the covariate-blind contrast under this draw (diagnostic).
  ols_pt_gap <- mean((trd_X + nl_trd)[A == 1]) - mean((trd_X + nl_trd)[A == 0])

  list(observed = dat_long, e_X = e_X, theta_i = theta_i,
       true_atet = true_atet, ols_pt_gap = ols_pt_gap)
}

# ---- Measurement-error wrappers ---------------------------------------------
# All take and return the LONG frame.

# Classical additive outcome ME.
me_y_classical <- function(dat, sigma_eps_y) {
  dat |> mutate(Y = Y + rnorm(n(), 0, sigma_eps_y))
}

# Differential (by A) outcome ME, mean zero.
me_y_diff <- function(dat, sigma0, sigma1) {
  sd_vec <- if_else(dat$A == 1, sigma1, sigma0)
  dat |> mutate(Y = Y + rnorm(n(), 0, sd_vec))
}

# Contraction toward the pooled mean: dY* = (1+gamma) dY + noise.
me_y_nc_contraction <- function(dat, gamma, sigma) {
  ybar <- mean(dat$Y)
  dat |> mutate(Y = Y + gamma * (Y - ybar) + rnorm(n(), 0, sigma))
}

# v2.0: covariate-dependent, time-varying, NON-LINEAR outcome ME.
#   Y*_1 = Y_1 + rho * h(X) + noise,  h standardised non-linear in X.
# rho is expressed as a multiple of sigma_eps. Applied at t=1 only so it
# does not cancel in dY. Discriminates parametric from ML nuisances because
# h(X) is not in the span of the linear OR/PS models.
me_y_nc_covariate <- function(dat, rho, sigma, params = default_params("low")) {
  m <- params$cov_moments
  num <- \(x) if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
  h <- .h_nonlin((dat$c1 - m$c1_mean) / m$c1_sd, (dat$c3 - m$c3_mean) / m$c3_sd,
                 (num(dat$c4) - m$c4_mean) / m$c4_sd, (dat$c6 - m$c6_mean) / m$c6_sd,
                 num(dat$c2), num(dat$c7))
  dat |> mutate(Y = Y + if_else(time == 1L, rho * params$sigma_eps * h, 0) +
                        rnorm(n(), 0, sigma))
}

# v2.0: classical covariate ME parameterised by RELIABILITY r in (0,1]:
#   X*_k = X_k + U_k, U_k ~ N(0, sd(X_k)^2 (1/r - 1)),  so Var(X_k)/Var(X*_k) = r.
# Same draw in both periods (time-invariant per unit).
me_x_classical <- function(dat, reliability = 0.75, covs = ME_X_COVS) {
  unit <- dat |> distinct(id, across(all_of(covs)))
  for (v in covs) {
    sd_v <- sd(as.numeric(unit[[v]]))
    unit[[paste0("u_", v)]] <- rnorm(nrow(unit), 0, sd_v * sqrt(1 / reliability - 1))
  }
  out <- dat |> left_join(unit |> select(id, starts_with("u_")), by = "id")
  for (v in covs) out[[v]] <- as.numeric(out[[v]]) + out[[paste0("u_", v)]]
  out |> select(-starts_with("u_"))
}

# v2.0: treatment-correlated covariate ME, scaled PER VARIABLE:
#   X*_k = X_k + sd(X_k) * (delta_a * A + N(0, sigma^2)).
me_x_nc <- function(dat, delta_a = 0.30, sigma = 0.10, covs = ME_X_COVS) {
  unit <- dat |> distinct(id, A, across(all_of(covs))) |>
    mutate(shock = delta_a * A + rnorm(n(), 0, sigma))
  out <- dat |> left_join(unit |> select(id, shock), by = "id")
  for (v in covs) {
    sd_v <- sd(as.numeric(unit[[v]]))
    out[[v]] <- as.numeric(out[[v]]) + sd_v * out$shock
  }
  out |> select(-shock)
}

# Treatment misclassification (X-independent), flips A per unit.
me_treat_misclass <- function(dat, alpha0, alpha1) {
  unit <- dat |> distinct(id, A) |>
    mutate(flip = if_else(A == 1, rbinom(n(), 1, alpha1), rbinom(n(), 1, alpha0)),
           A_star = if_else(flip == 1, 1L - as.integer(A), as.integer(A))) |>
    select(id, A_star)
  dat |> left_join(unit, by = "id") |> mutate(A = A_star) |> select(-A_star)
}
