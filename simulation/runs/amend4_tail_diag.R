# Diagnostics for the fourth amendment's study (computed after its results were seen; explanatory,
# not part of its frozen design). Writes amendment_2026-10-11b/checks/tail_robust.csv,
# tail_extrapolation.csv, tail_rep389.csv, task_runtime.csv and tail_diag.txt.
#   cd simulation; Rscript runs/amend4_tail_diag.R
# 1. Robust summaries (median error, MAD, extreme replicates) for every cell x learner x method,
#    because the mean-based metrics of runs/amend3_aggregate.R are dominated by a few replicates in
#    the log-normal cells.
# 2. Extrapolation: for each n = 5,000 log-normal replicate, the largest amount by which a test-fold
#    W_1 exceeds the maximum W_1 of its training folds (data and folds only, no fitting), set against
#    the replicate's error.
# 3. One replicate (spline lasso, quadratic trend, rep 389, latent) refitted and decomposed by fold
#    and unit, to show how a single extreme control produces the error.
# 4. Task wall time by cell and learner, and by method.
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(digest)
})
source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
source("R_scripts/34_simex_fold.R"); source("R_scripts/36_gauss_design.R"); source("R_scripts/38_regcal.R")
source("R_scripts/39_regcal_fold.R"); source("R_scripts/40_shape_design.R"); source("R_scripts/41_spline_learner.R")
ROOT <- "amendment_2026-10-11b/study"; OUT <- "amendment_2026-10-11b/checks"
dir.create(OUT, showWarnings = FALSE)
seed_of <- \(...) as.integer(digest::digest2int(paste(..., sep = "|")) %% .Machine$integer.max)   # as in amend4_study.R

fs <- list.files(file.path(ROOT, "raw"), "\\.rds$", recursive = TRUE, full.names = TRUE)
obj <- lapply(fs, readRDS)
rows <- bind_rows(lapply(obj, `[[`, "rows")) |> mutate(err = est - true_att)

# 1. Robust summaries.
robust <- rows |> group_by(cell, learner, method) |> summarise(
  mean_err = mean(err), median_err = median(err), sd = sd(est), mad = mad(est),
  se_over_sd = mean(se) / sd(est), median_se_over_mad = median(se) / mad(est),
  min_err = min(err), max_err = max(err), n_abs_err_gt1 = sum(abs(err) > 1),
  median_clip_frac = median(clip_frac), .groups = "drop")
write.csv(robust, file.path(OUT, "tail_robust.csv"), row.names = FALSE)

# 2. Extrapolation beyond the training folds' range of W_1.
excess <- bind_rows(lapply(c("lnorm_linear_r090_n5000", "lnorm_quad_r090_n5000"), \(cell) {
  p <- strsplit(cell, "_")[[1]]; gd <- shape_design(p[1], p[2], 0.9, omega = 1.0)
  bind_rows(lapply(1:500, \(rep) {
    set.seed(seed_of(cell, rep, "data")); g <- shape_draw(gd, 5000L)
    set.seed(seed_of(cell, rep, "folds")); fold <- .strat_folds(g$A, 5L)
    ex <- \(w) max(vapply(1:5, \(k) max(w[fold == k]) - max(w[fold != k]), 0))
    tibble(cell, rep, excess_latent = ex(g$XL[, 1]), excess_obs = ex(g$X[, 1]))
  }))
}))
extrap <- rows |> inner_join(excess, by = c("cell", "rep")) |>
  mutate(excess = if_else(method == "latent", excess_latent, excess_obs),
         excess_band = cut(excess, c(-Inf, 0, 2, 5, 10, Inf))) |>
  group_by(cell, learner, method, excess_band) |>
  summarise(n = n(), n_abs_err_gt1 = sum(abs(err) > 1), mean_err = mean(err), median_err = median(err), .groups = "drop")
write.csv(extrap, file.path(OUT, "tail_extrapolation.csv"), row.names = FALSE)

# 3. Replicate 389 of the spline-lasso quadratic cell, latent covariates, refitted with the study's seeds.
cell <- "lnorm_quad_r090_n5000"; learner <- "glmnet_ns"; REP <- 389L
gd <- shape_design("lnorm", "quad", 0.9, omega = 1.0)
set.seed(seed_of(cell, REP, "data")); g <- shape_draw(gd, 5000L)
set.seed(seed_of(cell, REP, "folds")); fold <- .strat_folds(g$A, 5L)
set.seed(seed_of(cell, REP, learner, "latent"))
L <- FOLD_LEARNERS[[learner]]; X <- g$XL; A <- g$A; dY <- g$dY
rep389 <- bind_rows(lapply(1:5, \(k) {
  tr <- which(fold != k); te <- which(fold == k); ctrl <- tr[A[tr] == 0]
  e_raw <- L(X[tr, , drop = FALSE], A[tr], "ps")(X[te, , drop = FALSE])
  m0 <- L(X[ctrl, , drop = FALSE], dY[ctrl], "or")(X[te, , drop = FALSE])
  e <- pmin(pmax(e_raw, PS_TRIM[1]), PS_TRIM[2]); a <- A[te]; r <- dY[te] - m0
  w0 <- (1 - a) * e / (1 - e); w0 <- w0 / mean(w0); w1 <- a / mean(a)
  j <- which.max(abs(w0 * r))
  tibble(fold = k, tau_fold = mean(w1 * r) - mean(w0 * r), top_unit = te[j], top_w1 = X[te[j], 1],
         train_max_w1 = max(X[tr, 1]), top_e_raw = e_raw[j], top_w0_normalised = w0[j], top_dY = dY[te[j]],
         top_m0 = m0[j], top_contrib_to_tau0 = w0[j] * r[j] / length(te))
}))
write.csv(rep389, file.path(OUT, "tail_rep389.csv"), row.names = FALSE)

# 4. Task wall time.
rt <- tibble(cl = basename(dirname(fs)), task_s = vapply(obj, `[[`, 0, "runtime_s")) |> group_by(cl) |>
  summarise(n = n(), median_s = median(task_s), mean_s = mean(task_s), p90_s = quantile(task_s, 0.9),
            max_s = max(task_s), hours = sum(task_s) / 3600)
rt_method <- rows |> group_by(cell, learner, method) |>
  summarise(median_s = median(runtime_s), mean_s = mean(runtime_s), p95_s = quantile(runtime_s, 0.95), .groups = "drop")
write.csv(bind_rows(mutate(rt, level = "task"), mutate(rt_method, level = "method")),
          file.path(OUT, "task_runtime.csv"), row.names = FALSE)

options(width = 200)
txt <- capture.output({
  cat("Fourth amendment tail diagnostics,", format(Sys.time()), "-", length(fs), "task files\n\n")
  cat("1. Robust summaries\n"); print(as.data.frame(robust), digits = 3)
  cat("\n2. Error by extrapolation band (n = 5,000 log-normal cells)\n"); print(as.data.frame(extrap), digits = 3)
  cat(sprintf("\n3. Spline lasso, quadratic trend, rep %d, latent: estimate %.3f\n", REP, mean(rep389$tau_fold)))
  print(as.data.frame(rep389), digits = 4)
  cat("\n4. Task wall time\n"); print(as.data.frame(rt), digits = 4); print(as.data.frame(rt_method), digits = 4)
})
writeLines(txt, file.path(OUT, "tail_diag.txt")); cat(txt, sep = "\n")
