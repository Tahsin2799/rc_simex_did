# Finite-support, algebraic and finite-difference checks for the full derivation
# in simulation/audit_2026-10-06/derivation_misclass_review.md.
# Run from repository root. Base R plus dplyr::if_else and cli error messages; no simulation production files changed.
source("simulation/R_scripts/33_misclass.R")
out_dir <- "simulation/audit_2026-10-06"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
checks <- list()
record <- function(name, value, expected = 0, tolerance = 1e-10) {
  discrepancy <- max(abs(value - expected))
  checks[[length(checks) + 1L]] <<- data.frame(
    check = name, value = if (length(value) == 1L) value else max(abs(value)),
    expected = if (length(expected) == 1L) expected else NA_real_,
    discrepancy = discrepancy, tolerance = tolerance,
    passed = is.finite(discrepancy) && discrepancy <= tolerance)
}
make_law <- function(a0 = .1, a1 = .2) {
  px <- c(.2, .3, .1, .4)
  e <- c(.12, .31, .58, .82)
  mu <- c(-.7, .2, 1.3, .4)
  theta <- c(.8, .2, .6, 1.1)
  kap <- 1 - a0 - a1
  q <- a0 + kap * e
  p1 <- (1 - a1) * e / q
  p0 <- a1 * e / (1 - q)
  m1 <- mu + theta * p1
  m0 <- mu + theta * p0
  dat <- expand.grid(x = 1:4, A = 0:1, B = 0:1, eps = c(-.4, .4))
  pb1 <- a0 + kap * dat$A
  dat$p <- px[dat$x] * dplyr::if_else(dat$A == 1, e[dat$x], 1 - e[dat$x]) *
    dplyr::if_else(dat$B == 1, pb1, 1 - pb1) / 2
  dat$Y <- mu[dat$x] + dat$A * theta[dat$x] + dat$eps
  dat$S <- (dat$B - a0) / kap
  list(px = px, e = e, mu = mu, theta = theta, kap = kap, q = q,
       m1 = m1, m0 = m0, c = m1 - m0, dat = dat, a0 = a0, a1 = a1,
       pi = sum(px * e), tau = sum(px * e * theta) / sum(px * e))
}
score <- function(law, h, m1, m0, target = "H") {
  dat <- law$dat; id <- dat$x; B <- dat$B; Y <- dat$Y
  f <- .mc_funs(law$a0, law$a1)[[target]]
  cbar <- m1[id] - m0[id]; hi <- h[id]
  N <- cbar * f$K(hi) + f$K(hi) *
    (B * (Y - m1[id]) / hi - (1 - B) * (Y - m0[id]) / (1 - hi)) +
    cbar * f$Kp(hi) * (B - hi)
  D <- f$L(hi) + f$Lp(hi) * (B - hi)
  list(N = N, D = D, value = sum(dat$p * N) / sum(dat$p * D))
}
law <- make_law()
with(law, record("law probabilities sum to one", sum(dat$p), 1))
record("Bayes attenuation identity", law$c - law$theta * law$kap *
         law$e * (1 - law$e) / (law$q * (1 - law$q)))
oracle <- score(law, law$q, law$m1, law$m0)
record("oracle correction identifies ATT", oracle$value, law$tau)
h <- law$a0 + law$kap * c(.24, .49, .39, .69)
m1bar <- law$m1 + c(.17, -.31, .23, .09)
m0bar <- law$m0 + c(-.12, .21, -.18, .14)
fit <- score(law, h, m1bar, m0bar)
mu_bar <- ((1 - law$a1) * (1 - h) * m0bar - law$a1 * h * m1bar) /
  (1 - law$a1 - h)
g <- (h - law$a0) / law$kap
id <- law$dat$x
signed_num <- (law$dat$S - g[id]) / (1 - g[id]) * (law$dat$Y - mu_bar[id])
record("pointwise signed score equals current H score", fit$N - signed_num)
record("pointwise denominator is corrected treatment", fit$D - law$dat$S)
for (target in c("H", "R")) {
  f <- .mc_funs(law$a0, law$a1)[[target]]
  cur <- score(law, h, m1bar, m0bar, target)
  r <- h - law$q; u1 <- m1bar - law$m1; u0 <- m0bar - law$m0
  RN <- sum(law$px * (law$c * (f$K(h) - f$K(law$q) - f$Kp(h) * r) +
      r * u1 * (f$K(h) / h - f$Kp(h)) +
      r * u0 * (f$K(h) / (1 - h) + f$Kp(h))))
  RD <- sum(law$px * (f$L(h) - f$L(law$q) - f$Lp(h) * r))
  N0 <- sum(law$px * law$c * f$K(law$q))
  D0 <- sum(law$px * f$L(law$q))
  record(paste(target, "exact numerator remainder"), sum(law$dat$p * cur$N) - N0, RN)
  record(paste(target, "exact denominator remainder"), sum(law$dat$p * cur$D) - D0, RD)
  record(paste(target, "exact ratio remainder"), cur$value - N0 / D0,
         (RN - N0 / D0 * RD) / (D0 + RD))
}
r <- h - law$q; b <- 1 - law$a1
RN_H <- sum(law$px * (-law$c * b * law$a1 * r^2 / ((b - h)^2 * (b - law$q)) +
                       r * (-law$a1 * h * (m1bar - law$m1) +
                         b * (1 - h) * (m0bar - law$m0)) / (b - h)^2))
record("H simplified exact numerator remainder", sum(law$dat$p * fit$N) -
         law$pi * law$tau, RN_H)
latent_value <- function(g, m) sum(law$dat$p *
  (law$dat$S - g[id]) / (1 - g[id]) * (law$dat$Y - m[id])) / law$pi
latent_remainder <- sum(law$px * (law$e - g) * (law$mu - mu_bar) / (1 - g)) / law$pi
record("latent exact global DR remainder", latent_value(g, mu_bar) - law$tau,
       latent_remainder)
record("latent DR correct trend arbitrary propensity", latent_value(g, law$mu), law$tau)
record("latent DR correct propensity arbitrary trend", latent_value(law$e, mu_bar), law$tau)
record("observed DR correct propensity wrong means", score(law, law$q, m1bar, m0bar)$value,
       law$tau)
observed_mean_bias <- score(law, h, law$m1, law$m0)$value - law$tau
record("correct observed means wrong q retains nonzero bias", as.numeric(observed_mean_bias < -1e-4), 1)
one_side <- make_law(a1 = 0)
h_one <- one_side$a0 + one_side$kap * c(.24, .49, .39, .69)
record("false-positive-only observed outcome DR exception",
       score(one_side, h_one, one_side$m1, one_side$m0)$value, one_side$tau)
fval <- law$kap * law$e * (1 - law$e) / (law$q * (1 - law$q))
w1 <- law$px * law$e / law$pi
tauR <- sum(w1 * law$theta * fval) / sum(w1 * fval)
cov1 <- sum(w1 * (law$theta - law$tau) * (fval - sum(w1 * fval)))
record("R target exact covariance displacement", tauR - law$tau, cov1 / sum(w1 * fval))
params <- expand.grid(a0 = c(0, .02, .15, .35), a1 = c(0, .05, .2, .4))
largest_r <- 0
for (ii in seq_len(nrow(params))) {
  a0 <- params$a0[ii]; a1 <- params$a1[ii]; kap <- 1 - a0 - a1
  if (kap <= 0) next
  for (ep in seq(.05, .95, .05)) {
    pp <- c(.3, .7); es <- c(ep, .97); pi <- sum(pp * es)
    qs <- a0 + kap * es; pB <- sum(pp * qs)
    rr <- pi * kap * (1 - es) / (pB * (1 - qs))
    largest_r <- max(largest_r, rr)
    if (any(rr <= 0) || any(rr > 1 + 1e-12)) cli::cli_abort("Attenuation weights outside (0, 1]: range {signif(range(rr), 4)}.")
  }
}
record("common-sign attenuation multiplier upper bound", largest_r, 1, 1e-12)
es <- c(.2, .8); th <- c(-2.8, .8); qs <- .1 + .7 * es
fs <- .7 * es * (1 - es) / (qs * (1 - qs))
tau_mixed <- sum(es * th) / sum(es)
naive_mixed <- sum(qs * fs * th) / sum(qs)
record("mixed-sign heterogeneous ATT sign reversal", as.numeric(tau_mixed > 0 && naive_mixed < 0), 1)
target_rates <- function(a0, a1) {
  kap <- 1 - a0 - a1
  sum(law$px * law$c * law$q * (1 - law$q) / (1 - a1 - law$q)) /
    sum(law$px * (law$q - a0) / kap)
}
gamma0 <- law$tau * (1 - law$pi) / (law$kap * law$pi)
gamma1 <- sum(law$px * law$theta * law$e^2 / (1 - law$e)) / (law$kap * law$pi)
eps <- 1e-6
der0 <- (target_rates(law$a0 + eps, law$a1) - target_rates(law$a0 - eps, law$a1)) / (2 * eps)
der1 <- (target_rates(law$a0, law$a1 + eps) - target_rates(law$a0, law$a1 - eps)) / (2 * eps)
record("external-rate alpha0 derivative", der0, gamma0, 1e-7)
record("external-rate alpha1 derivative", der1, gamma1, 1e-7)
psi <- (oracle$N - law$tau * law$dat$S) / law$pi
record("known-rate influence function centered", sum(law$dat$p * psi))
nu0 <- (1 - law$dat$A) * (law$dat$B - law$a0) / (1 - law$pi)
nu1 <- law$dat$A * ((1 - law$dat$B) - law$a1) / law$pi
record("validation alpha0 IF centered", sum(law$dat$p * nu0))
record("validation alpha1 IF centered", sum(law$dat$p * nu1))
record("validation IF covariance zero", sum(law$dat$p * nu0 * nu1))
record("validation alpha0 IF variance", sum(law$dat$p * nu0^2), law$a0 * (1 - law$a0) / (1 - law$pi))
record("validation alpha1 IF variance", sum(law$dat$p * nu1^2), law$a1 * (1 - law$a1) / law$pi)
all_checks <- do.call(rbind, checks)
write.csv(all_checks, file.path(out_dir, "math_misclass_repair_checks.csv"), row.names = FALSE)
report <- c(
  "Treatment-misclassification repair checks, 2026-10-06",
  sprintf("Passed %d of %d exact/finite-difference checks.", sum(all_checks$passed), nrow(all_checks)),
  sprintf("Correct observed m0,m1 + wrong q: asymptotic bias = %.12f", observed_mean_bias),
  sprintf("Same wrong q + correct latent mu: bias = %.3g", latent_value(g, law$mu) - law$tau),
  sprintf("Mixed-sign counterexample: ATT = %.12f, naive target = %.12f", tau_mixed, naive_mixed),
  sprintf("Validation derivatives: Gamma0 = %.12f, Gamma1 = %.12f", gamma0, gamma1),
  sprintf("Internal-validation cross-covariance (before sampling factor) = %.12f",
          sum(law$dat$p * psi * (gamma0 * nu0 + gamma1 * nu1))),
  "These finite-support checks verify algebra/counterexamples, not coverage or novelty."
)
writeLines(report, file.path(out_dir, "math_misclass_repair_report.txt"))
cat(paste(report, collapse = "\n"), "\n")
if (!all(all_checks$passed)) stop("At least one misclassification repair check failed.")
