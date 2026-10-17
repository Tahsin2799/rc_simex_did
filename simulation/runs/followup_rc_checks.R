# Checks of the regression-calibration comparator before its run (second amendment,
# Section 2). Writes amendment_2026-10-08/checks/rc_checks.csv.
#   cd simulation; Rscript runs/followup_rc_checks.R
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(digest)
})
source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
source("R_scripts/34_simex_fold.R"); source("R_scripts/36_gauss_design.R"); source("R_scripts/38_regcal.R")
dir.create("amendment_2026-10-08/checks", recursive = TRUE, showWarnings = FALSE)
gauss_draw_latent <- function(gd, n) {     # copy of the runner's function
  A  <- rbinom(n, 1, gd$pi)
  W  <- vapply(seq_len(2), \(c) A * gd$dmu[c] + rnorm(n, 0, sqrt(gd$v[c] - gd$s0[c])), numeric(n))
  Ws <- W + vapply(seq_len(2), \(c) rnorm(n, 0, sqrt(gd$s0[c])), numeric(n))
  V  <- rnorm(n)
  dY <- drop(W %*% gd$beta) + gd$gamma * V + gd$tau * A + rnorm(n)
  list(X = cbind(w1 = Ws[, 1], w2 = Ws[, 2], v = V), XL = cbind(w1 = W[, 1], w2 = W[, 2], v = V), A = A, dY = dY)
}
gd <- gauss_design(0.6, 0.6); res <- list()
set.seed(11); a <- gauss_draw(gd, 2000); set.seed(11); b <- gauss_draw_latent(gd, 2000)
res$draw_identical <- identical(a$X, b$X) && identical(a$A, b$A) && identical(a$dY, b$dY)

set.seed(12); g <- gauss_draw_latent(gd, 200000); fold <- .strat_folds(g$A, 5L)
est <- \(X) fold_drdid(X, g$A, g$dY, "glm", fold)$est
naive <- est(g$X); rc_g <- est(regcal(g$X, g$A, c("w1", "w2"), 0.6, TRUE))
rc_p <- est(regcal(g$X, g$A, c("w1", "w2"), 0.6, FALSE)); lat <- est(g$XL)
q0 <- gauss_q(gd, 0)
res$naive_vs_q0 <- naive - q0                 # naive targets q(0)
res$rc_group_minus_tau <- rc_g - gd$tau       # group-wise calibration is exact here
res$rc_pooled_minus_naive <- rc_p - naive     # affine map: same linear span
res$latent_minus_tau <- lat - gd$tau
# A finite-sample tolerance of 0.01 at n = 200,000 (SD of each estimate about 0.005).
pass <- c(draw_identical = res$draw_identical, naive_vs_q0 = abs(res$naive_vs_q0) < 0.01,
          rc_group = abs(res$rc_group_minus_tau) < 0.01, rc_pooled = abs(res$rc_pooled_minus_naive) < 1e-6,
          latent = abs(res$latent_minus_tau) < 0.01)
out <- data.frame(check = names(pass), value = c(NA, res$naive_vs_q0, res$rc_group_minus_tau,
                                                 res$rc_pooled_minus_naive, res$latent_minus_tau), pass = pass)
write.csv(out, "amendment_2026-10-08/checks/rc_checks.csv", row.names = FALSE)
print(out, digits = 4); cat(sprintf("q(0) = %.4f, tau = %.2f\n", q0, gd$tau))
cat(if (all(pass)) "ALL PASS\n" else "FAILURE\n")
