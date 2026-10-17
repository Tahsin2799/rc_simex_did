# Closed-form Gaussian design for the fold-calibrated SIMEX (amendment 2026-10-06).
#   source("R_scripts/36_gauss_design.R")
#
# A ~ Bern(pi); latent W_c | A ~ N(A dmu_c, v_c - s0_c); observed W*_c = W_c + U_c,
# U_c ~ N(0, s0_c) with s0 = (1 - r0) Var(W*); V ~ N(0, 1);
# dY = beta'W + gamma V + tau A + N(0, 1).
# W* | A then has variance v_c, so after adding lambda s_L of Gaussian noise
# (s_L = (1 - r_L) Var(W*), the reference covariance) every pseudo-law has a
# logit-linear propensity and a linear control outcome mean in (R, V). The
# parametric learner is correctly specified and the SIMEX curve is closed form:
#   q(lambda) = tau + sum_c beta_c dmu_c [1 - (v_c - s0_c) / (v_c + lambda s_L,c)],
# with q(t0) = tau at t0 = -s0 / s_L = -(1 - r0) / (1 - r_L).

gauss_design <- function(r0 = 0.6, r_L = r0) {
  gd <- list(pi = 0.4, dmu = c(0.6, 0.4), v = c(1, 1), beta = c(1, 0.8),
             gamma = 0.5, tau = 0.5, r0 = r0, r_L = r_L)
  gd$varWs <- gd$v + gd$pi * (1 - gd$pi) * gd$dmu^2
  gd$s0 <- (1 - r0) * gd$varWs
  gd$sL <- (1 - r_L) * gd$varWs
  gd$t0 <- -(1 - r0) / (1 - r_L)
  gd
}

gauss_draw <- function(gd, n) {
  A  <- rbinom(n, 1, gd$pi)
  W  <- vapply(seq_len(2), \(c) A * gd$dmu[c] + rnorm(n, 0, sqrt(gd$v[c] - gd$s0[c])), numeric(n))
  Ws <- W + vapply(seq_len(2), \(c) rnorm(n, 0, sqrt(gd$s0[c])), numeric(n))
  V  <- rnorm(n)
  dY <- drop(W %*% gd$beta) + gd$gamma * V + gd$tau * A + rnorm(n)
  list(X = cbind(w1 = Ws[, 1], w2 = Ws[, 2], v = V), A = A, dY = dY)
}

# Population curve on the reference scale (s = reference covariance vector).
gauss_q <- function(gd, lam, s = gd$sL)
  vapply(lam, \(l) gd$tau + sum(gd$beta * gd$dmu * (1 - (gd$v - gd$s0) / (gd$v + l * s))), 0)
# Gradient of q(lambda) in the reference covariance s, at s = sL.
gauss_D <- function(gd, lam)
  gd$beta * gd$dmu * (gd$v - gd$s0) * lam / (gd$v + lam * gd$sL)^2
# Third derivative; positive and decreasing in lambda.
gauss_q3 <- function(gd, lam)
  vapply(lam, \(l) sum(6 * gd$beta * gd$dmu * (gd$v - gd$s0) * gd$sL^3 / (gd$v + l * gd$sL)^4), 0)
# Population quadratic coefficients for grid `LAM` and the true M on [t0, max LAM].
gauss_truth <- function(gd, LAM = c(0, 0.5, 1, 2)) {
  Zg <- cbind(1, LAM, LAM^2)
  coef <- drop(solve(crossprod(Zg), crossprod(Zg, gauss_q(gd, LAM))))
  tg <- seq(-1, 0, length.out = 20001); q2 <- coef[1] + coef[2] * tg + coef[3] * tg^2
  vtx <- -coef[2] / (2 * coef[3])
  ext <- c(-1, 0, if (vtx > -1 && vtx < 0) vtx)
  q2e <- coef[1] + coef[2] * ext + coef[3] * ext^2
  list(coef = coef, q2 = \(t) coef[1] + coef[2] * t + coef[3] * t^2,
       range = c(min(q2e), max(q2e)), M = gauss_q3(gd, gd$t0), Q = gauss_q(gd, LAM))
}
