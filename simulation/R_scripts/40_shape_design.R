# Covariate-shape x trend design (third amendment, amendment_2026-10-11/AMENDMENT_3.md).
#   source("R_scripts/36_gauss_design.R"); source("R_scripts/40_shape_design.R")
#
# Extends the closed-form Gaussian design of 36_gauss_design.R while keeping its
# first two moments: A ~ Bern(0.4); latent W_c = A dmu_c + sqrt(v_c - s0_c) L_c
# with L_c standardised (mean 0, variance 1) and independent across c; observed
# W*_c = W_c + N(0, s0_c), s0 = (1 - r0) Var(W*); V ~ N(0, 1);
#   dY = beta'W + kappa_q (W_1^2 - E W_1^2) + gamma V + tau A + N(0, 1).
# shape = "gauss": L Gaussian (the closed-form design when kappa_q = 0).
# shape = "lnorm": L a standardised log-normal, exp(omega G) centred and scaled,
#   omega = 0.8 (skewness about 3.7, excess kurtosis about 31).
# trend = "linear" (kappa_q = 0) or "quad" (kappa_q = 0.5). omega (default 0.8) is the log-scale SD
# (fourth amendment: omega = 1.0, skewness 6.2).
# Means, variances and covariances of (W, W*, V, A) are the same for both shapes,
# so the best linear predictor of W given (W*, V) within each group, and hence
# regression calibration's map, has the same population form in all four cells.
# The effect is constant, so the realised ATT equals tau = 0.5 in every replicate.

SHAPE_OMEGA <- 0.8
SHAPE_KAPPA_Q <- 0.5

shape_design <- function(shape = c("gauss", "lnorm"), trend = c("linear", "quad"), r0 = 0.6, omega = SHAPE_OMEGA) {
  shape <- rlang::arg_match(shape); trend <- rlang::arg_match(trend)
  gd <- gauss_design(r0, r0)
  gd$shape <- shape; gd$trend <- trend; gd$omega <- omega
  gd$kappa_q <- if (trend == "quad") SHAPE_KAPPA_Q else 0
  gd$EW1sq <- (gd$v[1] - gd$s0[1]) + gd$pi * gd$dmu[1]^2   # E[W_1^2] = Var(L part) + E[(A dmu)^2]
  gd
}

.std_draw <- function(n, shape, w = SHAPE_OMEGA) {
  if (shape == "gauss") return(rnorm(n))
  (exp(w * rnorm(n)) - exp(w^2 / 2)) / sqrt((exp(w^2) - 1) * exp(w^2))
}

# Returns observed X = (w1, w2, v) and latent XL; RNG order: A, L1, L2, U1, U2, V, noise.
shape_draw <- function(gd, n) {
  A  <- rbinom(n, 1, gd$pi)
  W  <- vapply(seq_len(2), \(c) A * gd$dmu[c] + sqrt(gd$v[c] - gd$s0[c]) * .std_draw(n, gd$shape, gd$omega), numeric(n))
  Ws <- W + vapply(seq_len(2), \(c) rnorm(n, 0, sqrt(gd$s0[c])), numeric(n))
  V  <- rnorm(n)
  dY <- drop(W %*% gd$beta) + gd$kappa_q * (W[, 1]^2 - gd$EW1sq) + gd$gamma * V + gd$tau * A + rnorm(n)
  list(X = cbind(w1 = Ws[, 1], w2 = Ws[, 2], v = V), XL = cbind(w1 = W[, 1], w2 = W[, 2], v = V),
       A = A, dY = dY)
}
