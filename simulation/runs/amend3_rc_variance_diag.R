# Diagnostics for the regression-calibration variance (computed after the third amendment's study;
# explanatory, not part of its frozen design). Writes amendment_2026-10-11/checks/rc_var_components.csv
# and rc_var_tails.csv.
#   cd simulation; Rscript runs/amend3_rc_variance_diag.R
# 1. Which part of the calibration drives the calibration term? The contribution c_j of
#    39_regcal_fold.R is split into its group-mean, group-covariance and error-variance parts, in
#    the closed-form design (Gaussian, linear and quadratic trend) at n = 200,000, parametric
#    learner for the linear trend and lasso for the quadratic one.
# 2. Is the under-coverage at n = 1,000 driven by heavy tails? For rc_fold in the third amendment's
#    Gaussian cells, the mean SE is compared with the empirical SD and with a robust SD (IQR / 1.349).
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})
source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
source("R_scripts/34_simex_fold.R"); source("R_scripts/36_gauss_design.R")
source("R_scripts/39_regcal_fold.R"); source("R_scripts/40_shape_design.R")
OUT <- "amendment_2026-10-11/checks"

# Split .rc_contrib into parts: mean (lin), covariance (qf) and error variance (trS).
contrib_parts <- function(Z, A, parts, reliability) {
  cp <- parts$cp; me <- cp$me; n <- nrow(Z); p <- ncol(Z)
  Wm <- Z[, me, drop = FALSE]; ds <- (1 - reliability) * sweep(sweep(Wm, 2, colMeans(Wm))^2, 2, apply(Wm, 2, var))
  out <- list(mean = numeric(n), cov = numeric(n), s = numeric(n))
  for (a in 0:1) {
    ga <- cp$g[[a + 1]]; Ua <- parts$U[[a + 1]]; pia <- mean(A == a); sg <- if (a == 1) 1 else -1
    Sp <- matrix(0, p, p); Sp[cbind(me, me)] <- ga$s_eff
    SiU <- ga$Si %*% Ua; K <- SiU %*% Sp %*% ga$Si; Ks <- (K + t(K)) / 2
    ua <- sweep(Z, 2, ga$mu); ind <- (A == a) / pia
    out$s    <- out$s    + sg * (-drop(ds %*% diag(SiU)[me]))
    out$cov  <- out$cov  + sg * (ind * (rowSums((ua %*% Ks) * ua) - sum(ga$Sig * Ks)))
    out$mean <- out$mean + sg * (ind * drop(ua %*% t(Sp %*% ga$Si) %*% parts$gsum))
  }
  lapply(out, \(x) -x / parts$n1)
}
comp <- bind_rows(lapply(list(c("linear", "glm"), c("quad", "glmnet")), \(tl) {
  gd <- shape_design("gauss", tl[1]); set.seed(501); g <- shape_draw(gd, 200000); set.seed(502)
  f <- regcal_fold_drdid(g$X, g$A, g$dY, 1:2, 0.6, tl[2], .strat_folds(g$A, 5L), "fold")
  pc <- contrib_parts(g$X, g$A, f$parts, 0.6)
  parts_gap <- max(abs(pc$mean + pc$cov + pc$s - f$c))
  if (!isTRUE(parts_gap < 1e-10)) cli::cli_abort("The calibration-term components do not add up to the total: largest gap {signif(parts_gap, 3)}.")
  tibble(trend = tl[1], learner = tl[2], var_phi = var(f$phi), var_chi = var(f$chi),
         share_calib = var(f$c) / var(f$chi),
         var_mean = var(pc$mean), var_cov = var(pc$cov), var_s = var(pc$s),
         cov_phi_c = cov(f$phi, f$c))
}))
write.csv(comp, file.path(OUT, "rc_var_components.csv"), row.names = FALSE)
print(as.data.frame(comp), digits = 3)

fs <- list.files("amendment_2026-10-11/study/raw/shape", "\\.rds$", recursive = TRUE, full.names = TRUE)
rows <- bind_rows(lapply(fs, \(f) readRDS(f)$rows)) |> filter(is.na(error), method == "rc_fold", grepl("^gauss", cell))
tails <- rows |> group_by(cell, learner) |>
  summarise(reps = n(), se_over_sd = mean(se) / sd(est), se_over_robust_sd = mean(se) / (IQR(est) / 1.349),
            excess_kurtosis = mean(((est - mean(est)) / sd(est))^4) - 3, .groups = "drop")
write.csv(tails, file.path(OUT, "rc_var_tails.csv"), row.names = FALSE)
print(as.data.frame(tails), digits = 3)
