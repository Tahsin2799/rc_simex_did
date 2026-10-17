# True third-derivative bound M of the population SIMEX curve q(lambda) on [-1, 2] for the
# fourth amendment's cells (log-normal, omega = 1.0; Gaussian reference). q is evaluated by
# quadrature (amend4_population_lnorm.R) at Chebyshev nodes on [-1, 2], interpolated by a
# degree-16 Chebyshev polynomial, and differentiated analytically; also returns q2(-1) - tau and
# the regression-calibration population bias. No estimator is fitted.
#   cd simulation; Rscript runs/amend4_population_M.R
suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
})
cells <- data.frame(sl = c(1, 1, 1, 1, 0), r = c(0.9, 0.9, 0.75, 0.75, 0.9), bq = c(0, 0.5, 0, 0.5, 0))
one <- function(i) {
  source("runs/amend4_population_lnorm.R")
  sl <- cells$sl[i]; r <- cells$r[i]; bq <- cells$bq[i]
  b1 <- 1; b2 <- 0.8; dmu <- c(0.6, 0.4); pA <- 0.4
  varWs <- 1 + pA * (1 - pA) * dmu^2; s0 <- (1 - r) * varWs; sw <- sqrt(1 - s0)
  f1 <- \(w) b1 * w + bq * w^2; f2 <- \(w) b2 * w
  q <- function(lam) { e <- s0 * (1 + lam)
    term(f1, dmu[1], sw[1], sl, e[1]) + term(f2, dmu[2], sw[2], sl, e[2]) }
  deg <- 16; k <- 0:deg; xk <- cos((2 * k + 1) * pi / (2 * (deg + 1)))     # Chebyshev nodes on [-1, 1]
  a <- -0.995; b <- 2; lam_k <- (xk + 1) / 2 * (b - a) + a
  qk <- vapply(lam_k, q, 0)
  Tm <- outer(xk, 0:deg, \(x, j) cos(j * acos(x)))
  cf <- solve(Tm, qk)                                                    # interpolating coefficients
  p <- polynom_fun <- function(x, d) {                                   # d-th derivative in lambda
    co <- cf; for (it in seq_len(d)) { # Chebyshev derivative recurrence
      n <- length(co) - 1; dco <- numeric(max(n, 1)); if (n >= 1) { dco[n] <- 2 * n * co[n + 1]
        if (n >= 2) dco[n - 1] <- 2 * (n - 1) * co[n]
        if (n >= 3) for (j in (n - 2):1) dco[j] <- dco[j + 2] + 2 * j * co[j + 1] }
      dco[1] <- dco[1] / 2; co <- dco }
    drop(outer(x, 0:(length(co) - 1), \(xx, j) cos(j * acos(pmin(pmax(xx, -1), 1)))) %*% co) * (2 / (b - a))^d
  }
  xg <- seq(-1, 1, length.out = 2001)
  M <- max(abs(p(xg, 3)))
  lam <- c(0, 0.5, 1, 2); Q <- vapply(lam, q, 0); Zg <- cbind(1, lam, lam^2)
  c2 <- solve(crossprod(Zg), crossprod(Zg, Q))
  rc <- term(f1, dmu[1], sw[1], sl, s0[1], s0[1] * dmu[1] / (1 - s0[1])) +
        term(f2, dmu[2], sw[2], sl, s0[2], s0[2] * dmu[2] / (1 - s0[2]))
  interp_err <- abs(p(-1 + 2 * (0.3 - a) / (b - a) - 0, 0) - q(0.3))
  data.frame(cells[i, ], naive_bias = Q[1], rc_bias = rc, simex_bias = c2[1] - c2[2] + c2[3], M_true = M,
             q_at_minus1 = q(-0.995), interp_check = interp_err)
}
cl <- makeCluster(5); clusterExport(cl, "cells"); out <- do.call(rbind, parLapply(cl, seq_len(nrow(cells)), one)); stopCluster(cl)
write.csv(out, "amendment_2026-10-11b/population_M.csv", row.names = FALSE)
print(out, digits = 4)
