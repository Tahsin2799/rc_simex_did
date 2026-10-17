# BIHS illustration, additions for the comparison of the two corrections (14 October 2026, written after
# all simulation results and the SIMEX analysis of 02_analysis.R had been seen; reported as post hoc).
#   Rscript application/R_v2/04_rc_gap.R          (from repo root; ~10 min)
# Needs application/derived/analysis_child_panel.rds from 01_build_panel.R. Analytic cohort only.
#
# 1. Regression calibration with cross-fitted calibration (simulation/R_scripts/39_regcal_fold.R), lasso
#    learner, at assumed reliabilities 0.9-0.6 of log assets and log land, with the SE that includes the
#    calibration term (Proposition 2) and the SE that treats the calibration as fixed. Village clusters.
# 2. Calibration by the posterior mean: the latent law of each error-prone covariate is estimated within
#    each treatment group on the training folds by the NPMLE (R_scripts/43_pm_calibration.R) and every unit
#    is replaced by its posterior mean; then the same cross-fitted DR-DiD. The SE treats the calibration as
#    fixed. This is the check suggested in Section 5 of the manuscript, not an estimator with known inference.
# 3. The calibration-gap share of each error-prone covariate within each group (Section 4.4), with a
#    village bootstrap of its sampling variability.
# Folds are the same as in 02_analysis.R for the analytic cohort (same seed, first draw).

suppressPackageStartupMessages({ library(dplyr) })
source("simulation/R_scripts/30_core.R")
source("simulation/R_scripts/34_simex_fold.R")
source("simulation/R_scripts/39_regcal_fold.R")
source("simulation/R_scripts/43_pm_calibration.R")
out_t <- "application/output"; out_d <- "application/derived"
an <- readRDS(file.path(out_d, "analysis_child_panel.rds"))

ME <- c("log_assets", "log_land")
COVS <- c("age11", "female", "hh_size", "n_child_5_15", "max_adult_yrs", "head_female",
          "head_age", "head_yrs", "log_assets", "log_land")
design <- function(d) {
  X <- cbind(as.matrix(d[, COVS]), model.matrix(~ factor(division), d)[, -1, drop = FALSE])
  colnames(X) <- make.names(colnames(X))
  X[, apply(X, 2, \(v) any(v != v[1])), drop = FALSE]
}
R_GRID <- c(0.9, 0.8, 0.7, 0.6)
csum_se <- \(psi, cl) sqrt(drop(.if_vcov(psi, cl)))

set.seed(20261006)
X <- design(an); A <- an$A; dY <- an$dY; cl <- an$cluster
fold <- .cluster_folds(cl, 5L)                                   # identical to 02_analysis.R's first draw
me <- match(ME, colnames(X)); K <- max(fold); nk <- tabulate(fold, K); n <- nrow(X)

# Posterior-mean calibration on the training folds, then cross-fitted DR-DiD with the lasso.
pm_rc <- function(r) {
  L <- FOLD_LEARNERS[["glmnet"]]; tau <- numeric(K); psi <- numeric(n)
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k); Xc <- X
    for (j in me) {
      s <- (1 - r) * var(X[tr, j])
      for (a in 0:1) {
        fit <- npmle_fit(X[tr[A[tr] == a], j], s)
        rows <- which(A == a); Xc[rows, j] <- pm_mean(X[rows, j], fit, s)
      }
    }
    ctrl <- tr[A[tr] == 0]
    e <- L(Xc[tr, , drop = FALSE], A[tr], "ps")(Xc[te, , drop = FALSE])
    m0 <- L(Xc[ctrl, , drop = FALSE], dY[ctrl], "or")(Xc[te, , drop = FALSE])
    fs <- .fold_score(A[te], dY[te], e, m0, PS_TRIM)
    tau[k] <- fs$tau; psi[te] <- fs$psi * n / (K * nk[k])
  }
  list(est = mean(tau), psi = psi)
}

t0 <- Sys.time()
rc <- bind_rows(lapply(R_GRID, function(r) {
  f <- regcal_fold_drdid(X, A, dY, ME, r, learner = "glmnet", fold = fold, calibration = "fold")
  p <- pm_rc(r)
  bind_rows(
    data.frame(reliability = r, method = "Regression calibration", estimate = f$est,
               se = csum_se(f$chi, cl), se_fixed = csum_se(f$phi, cl), calib_share = f$calib_share,
               clip_frac = f$clip_frac),
    data.frame(reliability = r, method = "Posterior-mean calibration", estimate = p$est,
               se = NA_real_, se_fixed = csum_se(p$psi, cl), calib_share = NA_real_, clip_frac = NA_real_))
})) |> mutate(ci_lo = estimate - Z975 * coalesce(se, se_fixed), ci_hi = estimate + Z975 * coalesce(se, se_fixed))
write.csv(rc, file.path(out_t, "rc_estimates.csv"), row.names = FALSE)

# Gap share (Section 4.4): NPMLE posterior mean minus the within-group best linear predictor, over Var(W*).
gap_share <- function(y, s) {
  fit <- npmle_fit(y, s); pm <- pm_mean(y, fit, s)
  blp <- mean(y) + (1 - s / var(y)) * (y - mean(y))
  mean((pm - blp)^2) / var(y)
}
B_BOOT <- 100L
uv <- unique(cl)
gap <- bind_rows(lapply(R_GRID, function(r) bind_rows(lapply(ME, function(v) bind_rows(lapply(0:1, function(a) {
  w <- an[[v]]; s <- (1 - r) * var(w)
  est <- gap_share(w[A == a], s)
  set.seed(4000 + 100 * r + 10 * a + match(v, ME))
  bs <- vapply(seq_len(B_BOOT), function(b) {
    idx <- unlist(split(seq_along(cl), cl)[as.character(sample(uv, replace = TRUE))], use.names = FALSE)
    ww <- w[idx]; aa <- A[idx]; gap_share(ww[aa == a], (1 - r) * var(ww))
  }, 0)
  data.frame(reliability = r, covariate = v, group = a, n = sum(A == a), gap_share = est,
             boot_sd = sd(bs), boot_q05 = quantile(bs, 0.05), boot_q95 = quantile(bs, 0.95),
             skewness = mean(((w[A == a] - mean(w[A == a])) / sd(w[A == a]))^3),
             share_zero = mean(w[A == a] == min(w)))
}))))))
write.csv(gap, file.path(out_t, "gap_share.csv"), row.names = FALSE)

cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(rc, digits = 3, row.names = FALSE)
print(gap, digits = 3, row.names = FALSE)
