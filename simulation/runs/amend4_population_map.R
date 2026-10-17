# Population map for the fourth amendment: closed-form targets of regression
# calibration and quadratic SIMEX in a Gaussian design with group-specific
# covariate dispersion. No estimator is fitted.
#   cd simulation; Rscript runs/amend4_population_map.R
# Writes amendment_2026-10-11b/population_map.csv.
#
# Design: A ~ Bern(pi); W1 | A = a ~ N(d1 a, sig_a^2); W2 | A ~ N(0.4 A, sig2^2); V ~ N(0, 1),
# all independent given A; W*_c = W_c + N(0, s0_c), s0 = (1 - r) Var(W*); reliability r known;
#   D = b1 W1 + b2 W2 + bq (W1^2 - E W1^2) + 0.5 V + tau A + N(0, 1).
# Regression calibration (group-wise BLP, Appendix C.8): its target is tau + bq (Om1 - Om0), with
#   Om_a = sig_a^2 s0 / (sig_a^2 + s0), the within-group posterior variance of W1.
# SIMEX: with total error variance e = s0 (1 + lambda), the DR-DiD target is
#   q(lambda) = tau + b1 d1 (1 - k0) + b2 d2 (1 - k2)
#              + bq [d1^2 + sig1^2 - k0^2 (d1^2 + sig1^2 + e) - sig0^2 e / (sig0^2 + e)],
#   k0 = sig0^2 / (sig0^2 + e), k2 = sig2^2 / (sig2^2 + e2), and q(-1) = tau.
suppressPackageStartupMessages(library(dplyr))
dir.create("amendment_2026-10-11b", showWarnings = FALSE)

pop <- function(ratio, d1, bq, r = 0.6, pi = 0.4, vw = 0.5654, b1 = 1, b2 = 0.8, d2 = 0.4, sig2 = sqrt(0.5846), tau = 0.5) {
  s0sq <- vw / (pi * ratio + 1 - pi); s1sq <- ratio * s0sq           # pi sig1^2 + (1 - pi) sig0^2 = vw
  varW1s <- pi * s1sq + (1 - pi) * s0sq + pi * (1 - pi) * d1^2       # Var(W1) before error
  varW1s <- varW1s / r                                               # Var(W1*) = Var(W1) / r
  varW2s <- (sig2^2 + pi * (1 - pi) * d2^2) / r
  s01 <- (1 - r) * varW1s; s02 <- (1 - r) * varW2s
  q <- function(lam) {
    e1 <- s01 * (1 + lam); e2 <- s02 * (1 + lam)
    k0 <- s0sq / (s0sq + e1); k2 <- sig2^2 / (sig2^2 + e2)
    tau + b1 * d1 * (1 - k0) + b2 * d2 * (1 - k2) +
      bq * (d1^2 + s1sq - k0^2 * (d1^2 + s1sq + e1) - s0sq * e1 / (s0sq + e1))
  }
  lam <- c(0, 0.5, 1, 2); Zg <- cbind(1, lam, lam^2)
  cf <- solve(crossprod(Zg), crossprod(Zg, q(lam)))
  q2m1 <- cf[1] - cf[2] + cf[3]
  # Third derivative of q on [-1, 2] by central differences on a fine grid.
  tg <- seq(-0.99, 2, by = 0.01); h <- 1e-3
  q3 <- (q(tg + 2 * h) - 2 * q(tg + h) + 2 * q(tg - h) - q(tg - 2 * h)) / (2 * h^3)
  om <- \(s) s * s01 / (s + s01)
  tibble(r, ratio, d1, bq, sig0sq = s0sq, sig1sq = s1sq, naive_bias = q(0) - tau,
         rc_bias = bq * (om(s1sq) - om(s0sq)), simex_bias = q2m1 - tau, M_true = max(abs(q3)),
         q_check = q(-1) - tau)
}
grid <- expand.grid(ratio = c(1, 1.5, 2, 3, 4, 6), d1 = c(0, 0.3, 0.6), bq = c(0.5, 1), r = c(0.6, 0.75, 0.9))
map <- bind_rows(lapply(seq_len(nrow(grid)), \(i) do.call(pop, as.list(grid[i, ]))))
if (max(abs(map$q_check)) >= 1e-12) cli::cli_abort("Closed-form check of the population map failed: largest |q_check| = {signif(max(abs(map$q_check)), 3)}.")
write.csv(map, "amendment_2026-10-11b/population_map.csv", row.names = FALSE)
print(as.data.frame(map |> select(-q_check)), digits = 3)
