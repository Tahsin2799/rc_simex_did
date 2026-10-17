# BIHS illustration: estimates, diagnostics and measurement sensitivity.
#   Rscript application/R_v2/02_analysis.R          (from repo root; ~5 min)
# Needs application/derived/analysis_child_panel.rds from 01_build_panel.R.
#
# Interpretation is fixed in advance (paper1_rc_simex_did/notes/post_prod_sim_steps.md, Gate A): stipend
# receipt in 2015 is observed only among children attending in 2015, so the
# contrast below is an ASSUMPTION-DEPENDENT ILLUSTRATION of the measurement
# workflow, not an estimate of the policy effect.
#
# Target: unweighted ATT of 2015 stipend receipt on 2011-2018 gain in completed
# school years, analytic cohort. Inference: villages as independent clusters,
# whole villages assigned to cross-fitting folds. Survey weights are not used.

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(DRDID) })
source("simulation/R_scripts/30_core.R")
source("simulation/R_scripts/34_simex_fold.R")
out_t <- "application/output"; out_d <- "application/derived"
an <- readRDS(file.path(out_d, "analysis_child_panel.rds"))
set.seed(20261006)

ME <- c("log_assets", "log_land")
COVS <- c("age11", "female", "hh_size", "n_child_5_15", "max_adult_yrs", "head_female",
          "head_age", "head_yrs", "log_assets", "log_land")
design <- function(d) {
  X <- cbind(as.matrix(d[, COVS]), model.matrix(~ factor(division), d)[, -1, drop = FALSE])
  colnames(X) <- make.names(colnames(X))
  X[, apply(X, 2, \(v) any(v != v[1])), drop = FALSE]
}
R_L <- 0.6; R_GRID <- c(1, 0.9, 0.8, 0.7, 0.6); M_GRID <- c(0, 0.02, 0.05, 0.1, 0.2)
t_of <- \(r) -(1 - r) / (1 - R_L)
csum_se <- \(psi, cl) sqrt(drop(.if_vcov(psi, cl)))

# fold_drdid() (cross-fit DR-DiD, fold-normalised) comes from 34_simex_fold.R.

analyse <- function(d, label, with_xgb = TRUE) {
  X <- design(d); A <- d$A; dY <- d$dY; cl <- d$cluster
  fold <- .cluster_folds(cl, 5L)
  rows <- list()
  add <- \(method, est, se, note = "") rows[[length(rows) + 1L]] <<- data.frame(sample = label,
    method, estimate = est, se, ci_lo = est - Z975 * se, ci_hi = est + Z975 * se, note)
  # Unadjusted DiD
  p1 <- mean(A); mu1 <- mean(dY[A == 1]); mu0 <- mean(dY[A == 0])
  psi <- A / p1 * (dY - mu1) - (1 - A) / (1 - p1) * (dY - mu0)
  add("Unadjusted DiD", mu1 - mu0, csum_se(psi, cl))
  # Parametric improved DR-DiD (Sant'Anna-Zhao), logit/OLS on the raw columns
  r <- drdid_imp_panel(y1 = d$yrs_now18, y0 = d$yrs11, D = A, covariates = cbind(1, X),
                       boot = FALSE, inffunc = TRUE)
  add("Parametric DR-DiD", r$ATT, csum_se(r$att.inf.func, cl))
  fl <- fold_drdid(X, A, dY, "glmnet", fold)
  add("Lasso DR-DiD (cross-fit)", fl$est, csum_se(fl$psi, cl))
  if (with_xgb) { fx <- fold_drdid(X, A, dY, "xgb", fold)
    add("Boosted DR-DiD (cross-fit)", fx$est, csum_se(fx$psi, cl)) }
  # Fold-calibrated SIMEX on the reference reliability r_L = 0.6
  cv <- simex_fold_curve(X, A, dY, ME, R_L, B = 20L, K = 5L, learner = "glmnet",
                         calibration = "fold", fold = fold, cluster = cl)
  for (r in R_GRID[-1]) { p <- simex_fold_point(cv, t_of(r), 0)
    add(sprintf("SIMEX-lasso, reliability %.1f (M = 0)", r), p[["est"]], p[["se"]],
        sprintf("se without calibration term %.4f", p[["se_fixed"]])) }
  grid <- expand.grid(reliability = R_GRID, M = M_GRID) |> rowwise() |>
    mutate(p = list(simex_fold_point(cv, t_of(reliability), M)),
           estimate = p[["est"]], se = p[["se"]], kappa = p[["kappa"]],
           ci_lo = p[["ci_lo"]], ci_hi = p[["ci_hi"]]) |> select(-p) |> ungroup() |>
    mutate(sample = label, .before = 1)
  rng <- bind_rows(lapply(M_GRID, \(M) { z <- simex_fold_range(cv, -1, 0, M)
    data.frame(sample = label, M, range_lo = z[["range_lo"]], range_hi = z[["range_hi"]],
               outer_lo = z[["ci_lo"]], outer_hi = z[["ci_hi"]],
               outer_lo_conservative = z[["ci_lo_conservative"]],
               outer_hi_conservative = z[["ci_hi_conservative"]]) }))
  curve <- data.frame(sample = label, lambda = cv$lam, Q = cv$Q,
                      se = sqrt(diag(cv$V_eta)), sd_over_b = c(NA, apply(cv$Qjb[-1, ], 1, sd)))
  # Overlap and balance
  e <- pmin(pmax(fl$e_hat, PS_TRIM[1]), PS_TRIM[2]); w0 <- (1 - A) * e / (1 - e)
  smd <- \(x, w) (mean(x[A == 1]) - weighted.mean(x[A == 0], w[A == 0])) /
                 sqrt((var(x[A == 1]) + var(x[A == 0])) / 2)
  bal <- data.frame(sample = label, covariate = colnames(X),
                    mean_treated = colMeans(X[A == 1, , drop = FALSE]),
                    mean_control = colMeans(X[A == 0, , drop = FALSE]),
                    smd_raw = apply(X, 2, smd, w = rep(1, length(A))),
                    smd_weighted = apply(X, 2, smd, w = w0), row.names = NULL)
  ov <- data.frame(sample = label, group = c("treated", "control"),
                   n = c(sum(A), sum(1 - A)), clusters = c(n_distinct(cl[A == 1]), n_distinct(cl[A == 0])),
                   ps_p01 = c(quantile(fl$e_hat[A == 1], .01), quantile(fl$e_hat[A == 0], .01)),
                   ps_median = c(median(fl$e_hat[A == 1]), median(fl$e_hat[A == 0])),
                   ps_p99 = c(quantile(fl$e_hat[A == 1], .99), quantile(fl$e_hat[A == 0], .99)),
                   frac_clipped = c(mean(fl$e_hat[A == 1] < .02 | fl$e_hat[A == 1] > .98),
                                    mean(fl$e_hat[A == 0] < .02 | fl$e_hat[A == 0] > .98)),
                   control_weight_ess = c(NA, sum(w0)^2 / sum(w0^2)))
  list(est = do.call(rbind, rows), grid = grid, range = rng, curve = curve, balance = bal,
       overlap = ov, cv = cv, D = cv$D)
}

t0 <- Sys.time()
main <- analyse(an, "analytic cohort")
sec  <- analyse(filter(an, enrolled_now15 %in% 1L), "attending in 2015 (conditions on a post-baseline variable)",
                with_xgb = FALSE)
res <- list(main, sec)
cat_res <- \(k) do.call(rbind, lapply(res, `[[`, k))
est <- cat_res("est"); write.csv(est, file.path(out_t, "estimates.csv"), row.names = FALSE)
write.csv(cat_res("grid"), file.path(out_t, "reliability_M_grid.csv"), row.names = FALSE)
write.csv(cat_res("range"), file.path(out_t, "range_intervals.csv"), row.names = FALSE)
write.csv(cat_res("curve"), file.path(out_t, "simex_curve.csv"), row.names = FALSE)
write.csv(cat_res("balance"), file.path(out_t, "balance.csv"), row.names = FALSE)
write.csv(cat_res("overlap"), file.path(out_t, "overlap.csv"), row.names = FALSE)
saveRDS(lapply(res, `[[`, "cv"), file.path(out_d, "simex_curves.rds"))

# Parallel-trends sensitivity: an assumed untreated-trend gap delta (years)
# shifts the target to ATT - delta; delta* is where the 95% interval reaches 0.
pt <- est |> filter(method %in% c("Parametric DR-DiD", "Lasso DR-DiD (cross-fit)")) |>
  crossing(delta = seq(-0.5, 0.5, by = 0.1)) |>
  mutate(estimate_delta = estimate - delta, ci_lo_delta = ci_lo - delta, ci_hi_delta = ci_hi - delta)
breakdown <- est |> mutate(delta_star_lower = ci_lo, delta_star_upper = ci_hi) |>
  select(sample, method, estimate, delta_star_lower, delta_star_upper)
write.csv(pt, file.path(out_t, "pt_sensitivity.csv"), row.names = FALSE)
write.csv(breakdown, file.path(out_t, "pt_breakdown.csv"), row.names = FALSE)

cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
print(est, digits = 3, row.names = FALSE)
print(main$range, digits = 3, row.names = FALSE)
print(main$overlap, digits = 3, row.names = FALSE)
print(main$curve, digits = 4, row.names = FALSE)
cat("Calibration gradient D (rows lambda, cols log_assets, log_land):\n"); print(signif(main$D, 3))
