# Resolution of the gap-share diagnostic (14 October 2026, post hoc, for Sections 6.3 and 8).
# The gap share of a Gaussian covariate is zero; its estimate is not, because the NPMLE overfits.
# This script gives the distribution of the estimate for Gaussian samples of the sizes used in the
# manuscript (the BIHS groups, 220 and 682), so an observed share can be compared with it.
#   cd simulation; Rscript runs/gap_resolution.R
suppressPackageStartupMessages(library(dplyr))
source("R_scripts/43_pm_calibration.R")
gap_share <- function(y, s) {
  fit <- npmle_fit(y, s); pm <- pm_mean(y, fit, s)
  blp <- mean(y) + (1 - s / var(y)) * (y - mean(y))
  mean((pm - blp)^2) / var(y)
}
res <- bind_rows(lapply(c(220L, 682L), \(n) bind_rows(lapply(c(0.6, 0.7, 0.8, 0.9), \(r) {
  set.seed(9000 + n + round(100 * r))
  est <- vapply(1:200, \(b) { y <- rnorm(n); gap_share(y, (1 - r) * var(y)) }, 0)
  data.frame(n, reliability = r, mean = mean(est), q50 = median(est), q95 = quantile(est, 0.95), q99 = quantile(est, 0.99))
}))))
write.csv(res, "population_map/gap_resolution.csv", row.names = FALSE)
print(res, digits = 3)
