# Independent population-quadrature checks of the SIMEX repair derivations.
# No fitted learners, production estimates, survey data, or asymptotic claims.
# Run at repo root: Rscript --vanilla simulation/runs/math_simex_repair_checks.R
out <- "simulation/audit_2026-10-06"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
record <- function(name, actual, expected = 0, tolerance = 1e-8) {
  discrepancy <- max(abs(actual - expected))
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, discrepancy = discrepancy, tolerance = tolerance,
    pass = is.finite(discrepancy) && discrepancy < tolerance)
}

# Derive weights independently of production .extrap_weights().
lambda <- c(0, .5, 1, 2)
design <- cbind(1, lambda, lambda^2)
weights <- as.vector(design %*% solve(crossprod(design), c(1, -1, 1)))
record("Production-grid weights equal exact rational derivation",
       weights, c(36, -8, -27, 10) / 11)
record("Extrapolation reproduces powers zero through two",
       as.vector(crossprod(design, weights)), c(1, -1, 1))

kernel_sum <- function(u) vapply(u, \(v) {
  sum(weights * pmax(lambda - v, 0)^2) / 2
}, numeric(1))
kernel_piece <- function(u) {
  dplyr::case_when(
    u <= 0  ~ (u + 1)^2 / 2,
    u <= .5 ~ -25 * u^2 / 22 + u + .5,
    u <= 1  ~ -17 * u^2 / 22 + 7 * u / 11 + 13 / 22,
    .default = 5 * (u - 2)^2 / 11
  )
}
knots <- c(-1, 0, .5, 1, 2)
record("Piecewise Peano kernel agrees at all knots", kernel_piece(knots), kernel_sum(knots))
locations <- seq(-1, 2, length.out = 1001)
record("Piecewise Peano kernel agrees between knots", kernel_piece(locations), kernel_sum(locations))
record("Peano kernel is nonnegative on checked domain", min(kernel_piece(locations)), 0)
integrals <- vapply(seq_len(4), \(j) {
  integrate(kernel_piece, knots[j], knots[j + 1], rel.tol = 1e-12)$value
}, numeric(1))
record("Four exact Peano kernel integrals", integrals, c(1/6, 173/528, 163/528, 5/33))
record("Sharp smooth-function bias constant", sum(integrals), 21/22)
M <- .7
cubic <- function(x) .4 - .2 * x + .1 * x^2 + M * x^3 / 6
record("Cubic with constant third derivative attains bound",
       sum(weights * cubic(lambda)) - cubic(-1), M * 21/22)
for (curve in c("exp", "sin")) {
  fn <- get(curve)
  third <- if (curve == "exp") exp else function(x) -cos(x)
  remainder <- sum(vapply(seq_len(4), \(j) {
    integrate(\(u) kernel_piece(u) * third(u),
              knots[j], knots[j+1], rel.tol = 1e-12)$value
  }, numeric(1)))
  record(paste("Peano remainder identity for", curve),
         sum(weights * fn(lambda)) - fn(-1), remainder)
}

# An observed law with W* in {-1,1}, P(A=1|W*)=(.2,.8),
# and Y=.4 W*+.2 A. Further Gaussian noise R=W*+sqrt(lambda*s) Z
# creates nontrivial conditional nuisances that are available exactly.
# Integrating over each original atom and Z preserves the common-row law.
w_values <- c(-1, 1)
pw <- c(.5, .5)
prop <- c(.2, .8)
pi_a <- sum(pw * prop)
atoms <- expand.grid(w_index = 1:2, a = 0:1)
atoms$prob <- pw[atoms$w_index] * dplyr::if_else(atoms$a == 1,
                                               prop[atoms$w_index], 1 - prop[atoms$w_index])
atoms$y <- .4 * w_values[atoms$w_index] + .2 * atoms$a
normal_expectation <- function(fn) integrate(
  \(z) fn(z) * dnorm(z), -9, 9, rel.tol = 1e-10, subdivisions = 500L)$value
nuisances <- function(r, lam, s) {
  d1 <- dnorm(r, w_values[1], sqrt(lam*s))
  d2 <- dnorm(r, w_values[2], sqrt(lam*s))
  marg <- pw[1]*d1 + pw[2]*d2
  untreated_density <- pw[1]*(1-prop[1])*d1 + pw[2]*(1-prop[2])*d2
  ps <- (pw[1]*prop[1]*d1 + pw[2]*prop[2]*d2) / marg
  m0 <- (pw[1]*(1-prop[1])*.4*w_values[1]*d1 +
         pw[2]*(1-prop[2])*.4*w_values[2]*d2) / untreated_density
  list(ps = ps, m0 = m0)
}
target <- function(lam, s) {
  sum(vapply(1:2, \(k) {
    pw[k]*prop[k]/pi_a * normal_expectation(\(z) {
      r <- w_values[k] + sqrt(lam*s)*z
      .4*w_values[k] + .2 - nuisances(r, lam, s)$m0
    })
  }, numeric(1)))
}
score_moments <- function(lam, s) {
  q <- target(lam, s)
  ans <- lapply(seq_len(nrow(atoms)), \(k) {
    aa <- atoms$a[k]; yy <- atoms$y[k]; ww <- w_values[atoms$w_index[k]]
    influence <- function(z) {
      nu <- nuisances(ww + sqrt(lam*s)*z, lam, s)
      aa/pi_a*(yy-nu$m0-q) - (1-aa)*nu$ps/(pi_a*(1-nu$ps))*(yy-nu$m0)
    }
    c(mean = normal_expectation(influence),
      gradient = normal_expectation(\(z) influence(z)*(z^2-1)/(2*s)))
  })
  colSums(do.call(rbind, ans) * atoms$prob)
}
gradients <- list()
s0 <- .6
for (lam in c(.5, 1, 2)) {
  step <- 1e-4
  finite_diff <- (target(lam, s0+step)-target(lam, s0-step))/(2*step)
  finite_diff_half <- (target(lam, s0+step/2)-target(lam, s0-step/2))/step
  moments <- score_moments(lam, s0)
  record(paste("Pseudo-law EIF centered at lambda", lam), moments[["mean"]], 0, 1e-9)
  record(paste("Gaussian score gradient matches population difference at lambda", lam),
         moments[["gradient"]], finite_diff_half, 2e-7)
  record(paste("Finite difference stable to halving step at lambda", lam),
         finite_diff, finite_diff_half, 2e-7)
  gradients[[length(gradients)+1L]] <- data.frame(
    lambda = lam, covariance = s0, population_target = target(lam,s0),
    gaussian_score_gradient = moments[["gradient"]],
    central_difference = finite_diff_half)
}
gradient_table <- do.call(rbind, gradients)
if (any(abs(gradient_table$gaussian_score_gradient) <= 1e-3)) {
  cli::cli_abort("A score-based gradient is too close to zero to test: smallest |value| {signif(min(abs(gradient_table$gaussian_score_gradient)), 3)}.")
}
tab <- do.call(rbind, checks)
write.csv(tab, file.path(out, "math_simex_repair_checks.csv"), row.names = FALSE)
write.csv(gradient_table, file.path(out, "math_simex_calibration_gradient.csv"), row.names = FALSE)
report <- capture.output({
  cat("Independent base-R diagnostics for the SIMEX derivation repair\n\n")
  print(tab, row.names = FALSE, digits = 10)
  cat("\nPopulation Gaussian-kernel score derivative diagnostics:\n")
  print(gradient_table, row.names = FALSE, digits = 10)
  cat("\nNormal integrals truncated at |Z|=9; omitted tail probability:",
      2*pnorm(-9), "\n")
  cat("\nThese checks test algebra and population derivative identities.\n")
  cat("They do not prove learner rates, uniform asymptotics, causal identification, or novelty.\n")
})
writeLines(report, file.path(out, "math_simex_repair_report.txt"))
cat(paste(report, collapse = "\n"), "\n")
if (!all(tab$pass)) cli::cli_abort("Failed checks: {.val {tab$check[!tab$pass]}}.")
