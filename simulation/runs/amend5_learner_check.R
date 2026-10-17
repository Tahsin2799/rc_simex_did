# Code check for the fifth amendment's clamped spline learner (42_spline_clamped.R), run before the
# study. Fits "glmnet_ns" and "glmnet_nsc" with the same RNG state on fold 1 of replicate 1 of each
# spline cell (latent and observed covariates, both nuisances) and checks that (a) predictions at
# test points inside the training range are identical, (b) predictions beyond it equal the
# prediction at the clamped point. No estimator is computed.
#   cd simulation; Rscript runs/amend5_learner_check.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(digest)
})
source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
source("R_scripts/34_simex_fold.R"); source("R_scripts/36_gauss_design.R"); source("R_scripts/40_shape_design.R")
source("R_scripts/41_spline_learner.R"); source("R_scripts/42_spline_clamped.R")
seed_of <- \(...) as.integer(digest::digest2int(paste(..., sep = "|")) %% .Machine$integer.max)

res <- bind_rows(lapply(c("lnorm_linear_r090_n5000", "lnorm_quad_r090_n5000", "lnorm_linear_r090_n1000"), \(cell) {
  p <- strsplit(cell, "_")[[1]]; gd <- shape_design(p[1], p[2], 0.9, omega = 1.0)
  set.seed(seed_of(cell, 1L, "data")); g <- shape_draw(gd, as.integer(sub("n", "", p[4])))
  set.seed(seed_of(cell, 1L, "folds")); fold <- .strat_folds(g$A, 5L)
  tr <- which(fold != 1); te <- which(fold == 1); ctrl <- tr[g$A[tr] == 0]
  bind_rows(lapply(c("X", "XL"), \(xn) bind_rows(lapply(c("ps", "or"), \(type) {
    X <- g[[xn]]; idx <- if (type == "ps") tr else ctrl; y <- if (type == "ps") g$A[idx] else g$dY[idx]
    set.seed(1L); f0 <- FOLD_LEARNERS$glmnet_ns(X[idx, ], y, type)
    set.seed(1L); f1 <- FOLD_LEARNERS$glmnet_nsc(X[idx, ], y, type)
    lo <- apply(X[idx, ], 2, min); hi <- apply(X[idx, ], 2, max)
    Xte <- X[te, ]; inside <- apply(sweep(Xte, 2, lo, ">=") & sweep(Xte, 2, hi, "<="), 1, all)
    Xcl <- sweep(sweep(Xte, 2, lo, pmax), 2, hi, pmin)
    tibble(cell, covariates = xn, nuisance = type, n_test = length(te), n_outside = sum(!inside),
           max_diff_inside = max(abs(f0(Xte[inside, ]) - f1(Xte[inside, ]))),
           max_diff_clamped = max(abs(f1(Xte) - f0(Xcl))),
           max_change_outside = if (any(!inside)) max(abs(f0(Xte[!inside, , drop = FALSE]) - f1(Xte[!inside, , drop = FALSE]))) else 0)
  }))))
}))
options(width = 200); print(as.data.frame(res), digits = 3)
if (any(res$max_diff_inside != 0) || any(res$max_diff_clamped >= 1e-12)) {
  cli::cli_abort("Clamped learner check failed: largest difference inside the range {max(res$max_diff_inside)}, at the clamp {signif(max(res$max_diff_clamped), 3)}.")
}
cat("PASS: identical inside the training range; constant (clamped) beyond it\n")
