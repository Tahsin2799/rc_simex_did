# Calibration-gap share in the original design (14 October 2026, post hoc, written for Section 6).
# The diagnostic of Section 4.4 (NPMLE posterior mean minus the best linear predictor, divided by
# Var(W*)) is applied to each error-prone covariate within each treatment group of one large draw of the
# original design (nonlinear family, classical error). A Gaussian draw of the same size and error variance
# gives the resolution of the diagnostic, whose true share is zero.
#   cd simulation; Rscript runs/gap_original_design.R
suppressPackageStartupMessages(library(dplyr))
source("R_scripts/2_data_sim.R")
source("R_scripts/43_pm_calibration.R")
OUT <- "population_map"; dir.create(OUT, showWarnings = FALSE)
N0 <- 20000L                                                     # draws per group used by the diagnostic

gap_share <- function(y, s) {
  fit <- npmle_fit(y, s); pm <- pm_mean(y, fit, s)
  blp <- mean(y) + (1 - s / var(y)) * (y - mean(y))
  mean((pm - blp)^2) / var(y)
}

dat <- data_gen_v2(80000L, dgp_family = "nonlinear", seed = 20261014L, return_cf = TRUE)$observed |>
  filter(time == 0)
res <- bind_rows(lapply(c(0.6, 0.75), \(r) {
  set.seed(round(1e4 * r))
  bind_rows(lapply(ME_X_COVS, \(v) {
    w <- as.numeric(dat[[v]]); s <- var(w) * (1 / r - 1)        # error variance as in me_x_classical()
    ws <- w + rnorm(length(w), 0, sqrt(s))
    bind_rows(lapply(0:1, \(a) {
      y <- ws[dat$A == a]; y <- y[seq_len(min(N0, length(y)))]
      g <- rnorm(length(y), mean(y), sqrt(var(y)))               # Gaussian with the same size and variance
      data.frame(reliability = r, covariate = v, group = a, n = length(y), skew_w = mean(((w[dat$A == a] - mean(w[dat$A == a])) / sd(w[dat$A == a]))^3),
                 gap_share = gap_share(y, s), gaussian_resolution = gap_share(g, s))
    }))
  }))
}))
write.csv(res, file.path(OUT, "gap_original_design.csv"), row.names = FALSE)
print(res, digits = 3)
