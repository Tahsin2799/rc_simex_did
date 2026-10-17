# Independent numerical check of (I19), using Gaussian quadrature on a dense
# bounded grid. Does not call the manuscript's SIMEX or estimation code.
# Run from the repository root with base R (cli formats the failure message).
w <- seq(-9, 9, length.out = 1201L)
z <- seq(-14, 14, length.out = 1801L)
dw <- w[2L] - w[1L]
dz <- z[2L] - z[1L]
trap <- function(x, dx) dx * (sum(x) - (x[1L] + x[length(x)]) / 2)
p <- dnorm(w)
e <- 0.2 + 0.6 * plogis(0.8 * w)
mu <- 0.5 * sin(w) + 0.2 * w
f0 <- (1 - e) * p
f1 <- e * p
h0 <- mu * f0
pi1 <- trap(f1, dw)
treated_mean <- 0.2 + trap(mu * f1, dw) / pi1
weights <- rep(dw, length(w))
weights[c(1L, length(w))] <- dw / 2
dist <- outer(z, w, "-")
curve <- function(t, gradient = FALSE) {
  K <- dnorm(dist, sd = sqrt(t))
  a <- as.vector(K %*% (weights * f1))
  b <- as.vector(K %*% (weights * f0))
  h <- as.vector(K %*% (weights * h0))
  m <- h / b
  q <- treated_mean - trap(a * m, dz) / pi1
  if (!gradient) return(q)
  Kz <- -dist * K / t
  ap <- as.vector(Kz %*% (weights * f1))
  bp <- as.vector(Kz %*% (weights * f0))
  hp <- as.vector(Kz %*% (weights * h0))
  mp <- (hp - m * bp) / b
  lp <- ap / a - bp / b
  rhs <- trap(a * mp * lp, dz) / pi1
  c(q = q, heatflow_rhs = rhs)
}
result <- do.call(rbind, lapply(c(0.1, 0.6, 1.5), \(t) {
  analytic <- curve(t, TRUE)
  step <- 1e-4
  numeric <- (curve(t + step) - curve(t - step)) / (2 * step)
  data.frame(noise_variance = t, q = unname(analytic["q"]),
             derivative_finite_difference = numeric,
             derivative_heatflow = unname(analytic["heatflow_rhs"]),
             absolute_error = abs(numeric - analytic["heatflow_rhs"]))
}))
if (!all(is.finite(as.matrix(result))) || max(result$absolute_error) >= 1e-7) {
  cli::cli_abort("Heat-flow derivative check failed: largest error {signif(max(result$absolute_error), 3)}.")
}
write.csv(result, "simulation/audit_2026-10-06/identification_heatflow_checks.csv",
          row.names = FALSE)
print(result, row.names = FALSE, digits = 12)
cat("PASS: independent finite-difference and heat-flow integrals agree.\n")
