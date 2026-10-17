# Population targets of regression calibration and quadratic SIMEX in the skewed
# (log-normal) covariate design, by numerical integration. No estimator is fitted.
#   cd simulation; Rscript runs/amend4_population_lnorm.R
# Writes amendment_2026-10-11b/population_map_lnorm.csv.
#
# Design (R_scripts/40_shape_design.R, generalised): A ~ Bern(0.4); W_c = dmu_c A + sw_c L_c,
# L_c a standardised log-normal with log-scale SD sl (sl -> 0 gives the Gaussian design);
# W*_c = W_c + N(0, s0_c), s0 = (1 - r) Var(W*), Var(W* | A) = 1; V ~ N(0, 1);
#   D = b1 W1 + b2 W2 + bq W1^2 + 0.5 V + tau A + N(0, 1).
# Coordinates are independent given A, so each target is a sum of one-dimensional terms.
# With nonparametric nuisances:
#  * naive / SIMEX at total error variance e: target - tau = sum_c E_1[f_c(W)] - E_{y ~ W*_e | A=1} E_0[f_c(W) | W*_e = y];
#  * regression calibration: h_a(y) = y - s0 (y - dmu a) within group a (Var(W* | A) = 1), so a treated
#    unit's calibrated point has control preimage y + s0 dmu / (1 - s0), and
#    target - tau = sum_c E_1[f_c(W)] - E_{y ~ W* | A=1} E_0[f_c(W) | W* = y + s0 dmu / (1 - s0)].
suppressPackageStartupMessages(library(dplyr))
dir.create("amendment_2026-10-11b", showWarnings = FALSE)

gq <- local({ g <- seq(-8, 8, length.out = 4001); w <- dnorm(g); list(g = g, w = w / sum(w)) })
latent <- function(a, dmu, sw, sl) {
  L <- if (sl == 0) gq$g else (exp(sl * gq$g) - exp(sl^2 / 2)) / sqrt((exp(sl^2) - 1) * exp(sl^2))
  a * dmu + sw * L
}
# E_0[f(W) | W*_e = y] for a vector y (control group), by quadrature over the latent law.
post_mean <- function(y, f, dmu, sw, sl, e) {
  w <- latent(0, dmu, sw, sl); K <- outer(y, w, \(yy, ww) dnorm(yy, ww, sqrt(e)))
  drop(K %*% (f(w) * gq$w)) / drop(K %*% gq$w)
}
# E over y ~ W*_e | A = 1 of g(y): integrate over latent W_1 and the error by Gauss-Hermite.
gh <- local({ x <- seq(-7, 7, length.out = 141); w <- dnorm(x); list(x = x, w = w / sum(w)) })
# The posterior mean is evaluated on a fine grid of y and interpolated (it is smooth in y).
treated_mean <- function(f, dmu, sw, sl, e, shift = 0) {
  w1 <- latent(1, dmu, sw, sl)
  y <- outer(w1, sqrt(e) * gh$x, "+") + shift                  # n_latent x n_error
  yg <- seq(min(y), max(y), length.out = 3001)
  pm <- post_mean(yg, f, dmu, sw, sl, e)
  vals <- matrix(approx(yg, pm, xout = as.vector(y), rule = 2)$y, nrow(y))
  sum(gq$w * drop(vals %*% gh$w))
}
term <- function(f, dmu, sw, sl, e, shift = 0) {
  E1 <- sum(gq$w * f(latent(1, dmu, sw, sl)))
  E1 - treated_mean(f, dmu, sw, sl, e, shift)
}

pop_lnorm <- function(sl, r, bq, b1 = 1, b2 = 0.8, dmu = c(0.6, 0.4), pi = 0.4) {
  varWs <- 1 + pi * (1 - pi) * dmu^2; s0 <- (1 - r) * varWs; sw <- sqrt(1 - s0)
  f1 <- \(w) b1 * w + bq * w^2; f2 <- \(w) b2 * w
  q <- function(lam) { e <- s0 * (1 + lam)
    term(f1, dmu[1], sw[1], sl, e[1]) + term(f2, dmu[2], sw[2], sl, e[2]) }
  rc <- term(f1, dmu[1], sw[1], sl, s0[1], s0[1] * dmu[1] / (1 - s0[1])) +
        term(f2, dmu[2], sw[2], sl, s0[2], s0[2] * dmu[2] / (1 - s0[2]))
  lam <- c(0, 0.5, 1, 2); Q <- vapply(lam, q, 0); Zg <- cbind(1, lam, lam^2)
  cf <- solve(crossprod(Zg), crossprod(Zg, Q))
  tibble(sl, r, bq, naive_bias = Q[1], rc_bias = rc, simex_bias = cf[1] - cf[2] + cf[3])
}
if (sys.nframe() == 0) {
grid <- expand.grid(sl = c(0, 0.8, 1.0), r = c(0.6, 0.75, 0.9), bq = c(0, 0.5))
map <- bind_rows(lapply(seq_len(nrow(grid)), \(i) do.call(pop_lnorm, as.list(grid[i, ]))))
write.csv(map, "amendment_2026-10-11b/population_map_lnorm.csv", row.names = FALSE)
print(as.data.frame(map |> mutate(rc_over_simex = abs(rc_bias) / abs(simex_bias))), digits = 3)
}
