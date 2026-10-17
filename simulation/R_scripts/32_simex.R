# SIMEX-DR-DiD and reliability bounds, design v2.0 (paper Section 3.4-3.5,
# Theorems 1-2; PreReg v2.0 s4.3).
#   source("R_scripts/32_simex.R")
#
# Setting: X = (W, V); W (the ME_X_COVS columns) observed as W* = W + U with
# U ~ N(0, diag(sd(W*_k)^2 (1 - r_k))), reliability r_k = Var(W_k)/Var(W*_k).
#
# SIMEX (Cook & Stefanski 1994): for lambda in LAMBDAS and b = 1..B, add
# noise with variance lambda * sigma_u^2 to W*, recompute the cross-fit
# DR-DiD, average over b, fit a quadratic in lambda through the (J+1) grid
# points (lambda_0 = 0 is the naive estimate), extrapolate to lambda = -1.
#
# The fitted quadratic evaluated at lambda is a FIXED linear combination of
# the grid means, q_hat(lambda) = c(lambda)' tau_bar with
#     c(lambda) = Z (Z'Z)^{-1} z(lambda),  z(lambda) = (1, lambda, lambda^2)'.
# Standard error (Theorem 1(iii)): the same linear combination of the
# per-unit influence values, averaged over b within each lambda,
#     psi_S,i(lambda) = sum_j c_j(lambda) * mean_b psi_{i,j,b},
# se = sd(psi_S) / sqrt(N). This is exact for the linear combination and
# accounts for the covariance across lambda through the common units. The
# Stefanski-Cook (1995) per-lambda variance extrapolation is kept only as an
# attribute for comparison (it was ~20 % too small in the 2026-09-12 check).
#
# One fitted curve gives the estimate for ANY assumed reliability r' at no
# extra cost: the extrapolation point for r' is
#     lambda(r') = -(1 - r') / (1 - r_ref),
# where r_ref is the reliability used to generate the noise.
#
# Bounds (unknown common r in [r_L, 1]; Theorem 2): identified set = range of
# q(lambda) over [-1, 0]; estimated by the min/max of the fitted quadratic
# over a grid PLUS its exact vertex when the vertex is interior; endpoint
# SEs from the same influence-function construction at c(lambda_L),
# c(lambda_U); Imbens-Manski (2004) / Stoye (2009) interval.

LAMBDAS_DEFAULT <- c(0.5, 1.0, 2.0)
B_DEFAULT       <- 20L

# Extrapolation weights c(lambda) for a quadratic fit on grid `lam`.
.extrap_weights <- function(lam, at) {
  Z <- cbind(1, lam, lam^2)
  H <- Z %*% solve(crossprod(Z))            # (J+1) x 3
  # vapply pins the (J+1) x length(at) shape for both callers (a scalar `at` from
  # SIMEX, a grid from the bounds). sapply agrees on every `at` used here but
  # returns a list for an empty one, which would fail further down in crossprod.
  vapply(at, \(l) as.numeric(H %*% c(1, l, l^2)), numeric(nrow(H)))
}

# Fit the SIMEX curve. Returns the grid means tau_bar, the per-unit
# influence matrix psi_bar (N x (J+1), averaged over b), and the naive fit.
.simex_curve <- function(d, base_fn, reliability, lambdas, B, covs = ME_X_COVS) {
  N <- nrow(d)
  sd_star  <- vapply(covs, \(v) sd(as.numeric(d[[v]])), numeric(1))
  sigma_u2 <- sd_star^2 * (1 - reliability)          # implied ME variance
  J <- length(lambdas)
  tau_bar <- numeric(J + 1); psi_bar <- matrix(0, N, J + 1)
  se_l <- matrix(NA_real_, J, B); est_l <- matrix(NA_real_, J, B)
  r0 <- base_fn(d, return_if = TRUE)
  tau_bar[1] <- r0[["est"]]; psi_bar[, 1] <- attr(r0, "psi")
  for (j in seq_len(J)) {
    for (b in seq_len(B)) {
      dd <- d
      for (k in seq_along(covs)) {
        v <- covs[k]
        dd[[v]] <- as.numeric(dd[[v]]) + rnorm(N, 0, sqrt(lambdas[j] * sigma_u2[k]))
      }
      r <- base_fn(dd, return_if = TRUE)
      est_l[j, b] <- r[["est"]]; se_l[j, b] <- r[["se"]]
      psi_bar[, j + 1] <- psi_bar[, j + 1] + attr(r, "psi") / B
    }
    tau_bar[j + 1] <- mean(est_l[j, ])
  }
  lam <- c(0, lambdas)
  # Stefanski-Cook variance function on the grid (comparison only).
  tau2_sc <- c(r0[["se"]]^2, rowMeans(se_l^2) - apply(est_l, 1, var))
  list(lam = lam, tau_bar = tau_bar, psi_bar = psi_bar, naive = r0,
       coef_est = coef(lm(tau_bar ~ lam + I(lam^2))),
       coef_var_sc = coef(lm(pmax(tau2_sc, 1e-12) ~ lam + I(lam^2))))
}

# Evaluate the fitted curve and its IF-based SE at one or more lambda.
.simex_eval <- function(cv, at) {
  Cw  <- .extrap_weights(cv$lam, at)                 # (J+1) x length(at)
  est <- as.numeric(crossprod(Cw, cv$tau_bar))
  psi <- cv$psi_bar %*% Cw                            # N x length(at)
  se  <- unname(apply(psi, 2, sd) / sqrt(nrow(psi)))
  list(est = est, se = se, psi = psi)
}

.eval_quad <- function(cf, x) cf[1] + cf[2] * x + cf[3] * x^2

# SIMEX-corrected estimate assuming the reliability is KNOWN (= r_ref).
simex_drdid <- function(d, base_fn, reliability, lambdas = LAMBDAS_DEFAULT, B = B_DEFAULT) {
  cv <- .simex_curve(d, base_fn, reliability, lambdas, B)
  ev <- .simex_eval(cv, -1)
  out <- .out(ev$est, ev$se)
  attr(out, "se_sc") <- sqrt(max(unname(.eval_quad(cv$coef_var_sc, -1)), 1e-12))
  out
}

# Imbens-Manski / Stoye interval for theta in [theta_L, theta_U].
# se_lo, se_hi are on the estimate scale (already divided by sqrt N).
imbens_manski <- function(lo, hi, se_lo, se_hi, alpha = 0.05) {
  lo <- unname(lo); hi <- unname(hi); se_lo <- unname(se_lo); se_hi <- unname(se_hi)
  delta <- max(hi - lo, 0); smax <- max(se_lo, se_hi, 1e-12)
  f <- \(c) pnorm(c + delta / smax) - pnorm(-c) - (1 - alpha)
  c_a <- if (f(qnorm(1 - alpha)) >= 0) qnorm(1 - alpha) else
           uniroot(f, c(qnorm(1 - alpha), qnorm(1 - alpha / 2)))$root
  c(ci_lo = lo - c_a * se_lo, ci_hi = hi + c_a * se_hi)
}

# Partial identification when reliability is only known to lie in
# [r_L, 1]. r_ref (noise-generation reliability) is set to r_L so the
# extrapolation never goes beyond lambda = -1.
bounds_drdid <- function(d, base_fn, r_L = 0.6, lambdas = LAMBDAS_DEFAULT,
                         B = B_DEFAULT, grid = 41L) {
  cv <- .simex_curve(d, base_fn, reliability = r_L, lambdas, B)
  lam_grid <- seq(-1, 0, length.out = grid)
  a <- cv$coef_est
  if (is.finite(a[3]) && a[3] != 0) {                 # exact vertex if interior
    # unname(): a = coef(lm(...)) carries names, so vtx would be named "lam",
    # which propagates through .simex_eval() to the imbens_manski() result and
    # turns the ci_hi lookup below into NA whenever the upper endpoint is the
    # interior vertex (found by the 2026-09-12 dry run).
    vtx <- unname(-a[2] / (2 * a[3]))
    if (vtx > -1 && vtx < 0) lam_grid <- sort(c(lam_grid, vtx))
  }
  ev <- .simex_eval(cv, lam_grid)
  i_lo <- which.min(ev$est); i_hi <- which.max(ev$est)
  lo <- ev$est[i_lo]; hi <- ev$est[i_hi]
  ci <- imbens_manski(lo, hi, ev$se[i_lo], ev$se[i_hi])
  out <- .out(est = (lo + hi) / 2, se = NA_real_,
              ci_lo = unname(ci["ci_lo"]), ci_hi = unname(ci["ci_hi"]),
              lo = lo, hi = hi)
  attr(out, "lambda_LU") <- c(lam_grid[i_lo], lam_grid[i_hi])
  attr(out, "se_LU")     <- c(ev$se[i_lo], ev$se[i_hi])
  out
}
