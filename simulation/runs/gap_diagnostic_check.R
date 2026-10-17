# Can the calibration gap be estimated from the error-prone covariate alone? (14 October 2026, before the
# rewrite.) The gap share G = E_0[(E_0[W | W*] - BLP_0(W*))^2] / Var_0(W*) equals (s / Var W*)^2 J(W*) by
# Brown's identity. Estimator: NPMLE posterior means (R_scripts/43_pm_calibration.R) minus the sample best
# linear predictor, in one control group of n0 draws. Compared with the population value from
# R_scripts/44_population_np.R. 20 samples per design.
#   cd simulation; Rscript runs/gap_diagnostic_check.R
suppressPackageStartupMessages(library(dplyr))
source("R_scripts/44_population_np.R")
source("R_scripts/43_pm_calibration.R")
OUT <- "population_map"; dir.create(OUT, showWarnings = FALSE)
std_ln <- \(n, sl) if (sl == 0) rnorm(n) else (exp(sl * rnorm(n)) - exp(sl^2 / 2)) / sqrt((exp(sl^2) - 1) * exp(sl^2))
res <- bind_rows(lapply(c(0, 0.5, 1), \(sl) bind_rows(lapply(c(0.6, 0.9), \(r) {
  D <- np_design(sl, r, 0); s <- D$s0[1]; P <- D$P1
  truth <- s^2 * np_J(0, P, s)                       # Var(W*_1 | A = 0) = 1
  est <- vapply(1:20, \(b) {
    set.seed(5000 + 100 * sl + 10 * r + b)
    y <- P$sw0 * std_ln(3000, sl) + rnorm(3000, 0, sqrt(s))
    fit <- npmle_fit(y, s); pm <- pm_mean(y, fit, s)
    blp <- mean(y) + (1 - s / var(y)) * (y - mean(y))
    mean((pm - blp)^2) / var(y)
  }, 0)
  data.frame(sl, r, truth, est_mean = mean(est), est_sd = sd(est), est_q10 = quantile(est, 0.1), est_q90 = quantile(est, 0.9))
}))))
write.csv(res, file.path(OUT, "gap_diagnostic_check.csv"), row.names = FALSE)
print(res, digits = 3)
