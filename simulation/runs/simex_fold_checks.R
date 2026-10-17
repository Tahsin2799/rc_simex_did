# Implementation checks for the fold-calibrated SIMEX (R_scripts/34_simex_fold.R).
#   cd simulation; Rscript --vanilla runs/simex_fold_checks.R [--reps=400] [--workers=10]
#
# Three layers, written to amendment_2026-10-06/checks/:
#   A. deterministic algebra: weights, exact Peano constants, fold scores,
#      internal consistency of the returned curve, draw reproducibility;
#   B. one large sample from a closed-form Gaussian design: estimated grid
#      means and calibration gradient D against their exact population values;
#   C. Monte Carlo in that design: coverage of q2(-1), of the true ATT with and
#      without the declared third-derivative budget, and of the simultaneous
#      band, for fixed and fold-estimated calibration.
#
# The design is theorem-compatible: latent W | A Gaussian, classical Gaussian
# error, Y linear in W, so every pseudo-law has a logit-linear propensity and
# a linear outcome mean and the parametric "glm" learner is correctly
# specified. It checks the code against the derivation; it says nothing about
# lasso/boosting rates or about the production DGP.

suppressPackageStartupMessages({
  library(parallel)
})
source("R_scripts/30_core.R")
source("R_scripts/34_simex_fold.R")
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
R_MC <- as.integer(kv("reps", "400")); WORKERS <- as.integer(kv("workers", "10"))
out <- "amendment_2026-10-06/checks"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

checks <- list()
record <- function(layer, name, actual, expected = 0, tolerance = 1e-10) {
  disc <- max(abs(actual - expected))
  checks[[length(checks) + 1L]] <<- data.frame(layer = layer, check = name,
    discrepancy = disc, tolerance = tolerance, pass = is.finite(disc) && disc <= tolerance)
}

# ---- closed-form Gaussian design -----------------------------------------------
source("R_scripts/36_gauss_design.R")
GD <- gauss_design(r0 = 0.6)            # reference = true reliability
GD$r <- GD$r0
gen_gd <- \(n) gauss_draw(GD, n)
# Population curve at noise level lambda (reference covariance s) and its
# s-gradient at s0; q(-1) at s0 is the true ATT.
Q_pop  <- \(lam, s = GD$s0) gauss_q(GD, lam, s)
D_pop  <- \(lam) gauss_D(GD, lam)
q3_pop <- \(lam) gauss_q3(GD, lam)
LAM <- c(0, 0.5, 1, 2)
M_TRUE <- q3_pop(-1)                                     # q''' is decreasing in lambda
Qg <- vapply(LAM, Q_pop, 0)
coef_true <- drop(solve(crossprod(cbind(1, LAM, LAM^2)), crossprod(cbind(1, LAM, LAM^2), Qg)))
q2_true <- \(t) coef_true[1] + coef_true[2] * t + coef_true[3] * t^2

# ---- A. deterministic algebra ---------------------------------------------------
record("A", "c(-1) on production grid is (36,-8,-27,10)/11",
       drop(.c_weights(LAM, -1)), c(36, -8, -27, 10) / 11, 1e-12)
for (tt in c(-1, -0.4, 0.3)) record("A", sprintf("c(t) reproduces 1, t, t^2 at t = %.1f", tt),
  drop(crossprod(cbind(1, LAM, LAM^2), .c_weights(LAM, tt))), c(1, tt, tt^2), 1e-12)
record("A", "kappa(-1) = 21/22 on production grid (exact roots)", peano_kappa(LAM, -1), 21 / 22, 1e-12)
num_kappa <- function(lam, tt) {
  cw <- drop(.c_weights(lam, tt)); kn <- sort(unique(c(tt, lam[lam > tt])))
  sum(vapply(seq_len(length(kn) - 1), \(i) integrate(\(u) abs(vapply(u, \(x)
    sum(cw * pmax(lam - x, 0)^2) / 2, 0)), kn[i], kn[i + 1], rel.tol = 1e-11,
    subdivisions = 500L)$value, 0))
}
ALT <- c(0, 0.25, 0.5, 1, 1.5, 2)
for (g in list(LAM, ALT)) for (tt in c(-1, -0.7, -0.4, -0.1, 0))
  record("A", sprintf("kappa(%.1f) exact = numerical, grid (%s)", tt, paste(g, collapse = ",")),
         peano_kappa(g, tt), num_kappa(g, tt), 1e-7)
peano_rem <- \(tt) { cw <- drop(.c_weights(LAM, tt)); kn <- sort(unique(c(tt, LAM[LAM > tt])))
  sum(vapply(seq_len(length(kn) - 1), \(i) integrate(\(u) vapply(u, \(x)
    sum(cw * pmax(LAM - x, 0)^2) / 2 * q3_pop(x), 0), kn[i], kn[i + 1], rel.tol = 1e-11)$value, 0)) }
for (tt in c(-1, -0.5)) record("A", sprintf("Peano identity q2(t) - q(t) for closed-form curve, t = %.1f", tt),
  q2_true(tt) - Q_pop(tt), peano_rem(tt), 1e-9)
record("A", "closed-form extrapolation bias within 21M/22 (excess over bound, clipped at 0)",
       max(abs(q2_true(-1) - GD$tau) - 21 / 22 * M_TRUE, 0), 0, 0)

set.seed(11)
dd <- gen_gd(3000); e0 <- plogis(rnorm(3000)); m00 <- rnorm(3000)
one <- .fold_score(dd$A, dd$dY, e0, m00, PS_TRIM)
prod <- drdid_score(dd$A, dd$dY, e0, m00, return_if = TRUE)
record("A", "K = 1 fold score equals production pooled score (estimate)", one$tau, prod[["est"]], 1e-12)
record("A", "K = 1 fold score equals production pooled score (IF)", one$psi, attr(prod, "psi"), 1e-12)
record("A", "held-out score sums to zero within the fold", sum(one$psi), 0, 1e-9)

set.seed(12)
cvF <- simex_fold_curve(dd$X, dd$A, dd$dY, 1:2, GD$r, B = 3L, K = 5L, learner = "glm",
                        calibration = "fold")
record("A", "rowwise grid IF averages to zero (fold-normalised)", colMeans(cvF$zeta), 0, 1e-10)
record("A", "coef = L Q", cvF$coef,
       drop(solve(crossprod(cbind(1, cvF$lam, cvF$lam^2)), crossprod(cbind(1, cvF$lam, cvF$lam^2), cvF$Q))), 1e-12)
record("A", "Q is the mean of its (j, b) fits", cvF$Q[-1], rowMeans(cvF$Qjb[-1, ]), 1e-12)
record("A", "fold calibration = (1 - r) Var(W*) on training rows",
       cvF$S_folds[2, ], (1 - GD$r) * apply(dd$X[cvF$fold != 2, 1:2], 2, var), 1e-12)
record("A", "average fold calibration minus s_hat equals mean(h_s) up to O(1/n)",
       colMeans(cvF$S_folds) - (1 - GD$r) * apply(dd$X[, 1:2], 2, var) - colMeans(cvF$h_s), 0, 5e-3)
pt <- simex_fold_point(cvF, -1)
cw <- drop(.c_weights(cvF$lam, -1))
record("A", "point SE = sqrt(c' Var(eta) c / n)", pt[["se"]],
       sqrt(drop(t(cw) %*% cov(cvF$eta) %*% cw) / cvF$n), 1e-12)
record("A", "eta = zeta + h_s D'", cvF$eta, cvF$zeta + cvF$h_s %*% t(cvF$D), 1e-12)
record("A", "singleton clusters reproduce the iid variance",
       .if_vcov(cvF$eta, seq_len(cvF$n)), cov(cvF$eta) / cvF$n, 1e-15)
cl_id <- rep(seq_len(cvF$n / 3), each = 3)
Sg <- rowsum(cvF$eta %*% cw, cl_id)
record("A", "cluster variance = G/(G-1) sum (S_g - mean S)^2 / n^2",
       drop(t(cw) %*% .if_vcov(cvF$eta, cl_id) %*% cw),
       length(Sg) / (length(Sg) - 1) * sum((Sg - mean(Sg))^2) / cvF$n^2, 1e-15)
set.seed(12)
cvF2 <- simex_fold_curve(dd$X, dd$A, dd$dY, 1:2, GD$r, B = 3L, K = 5L, learner = "glm",
                         calibration = "fold", save_draws = TRUE)
record("A", "same seed reproduces the curve exactly", cvF2$Q, cvF$Q, 0)
set.seed(13)
cvX <- simex_fold_curve(dd$X, dd$A, dd$dY, 1:2, GD$r, B = 3L, K = 5L, learner = "glm",
                        calibration = "fixed", s_fixed = GD$s0)
record("A", "fixed calibration: eta equals zeta", cvX$eta, cvX$zeta, 0)
rg <- simex_fold_range(cvF, M = 0)
record("A", "range interval contains the estimated range",
       max(rg[["ci_lo"]] - rg[["range_lo"]], rg[["range_hi"]] - rg[["ci_hi"]], 0), 0, 0)
record("A", "numerical range interval inside closed-form conservative one",
       max(rg[["ci_lo_conservative"]] - rg[["ci_lo"]], rg[["ci_hi"]] - rg[["ci_hi_conservative"]], 0), 0, 1e-10)

# ---- B. one large sample, closed form -------------------------------------------
set.seed(2026)
big <- gen_gd(200000)
t0 <- Sys.time()
cvB <- simex_fold_curve(big$X, big$A, big$dY, 1:2, GD$r, B = 2L, K = 2L,
                        learner = "glm", calibration = "fold")
tB <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
seQ <- apply(cvB$eta, 2, sd) / sqrt(cvB$n)
large <- data.frame(j = 0:3, lambda = LAM, Q_hat = cvB$Q, Q_pop = Qg, se = seQ,
                    z = (cvB$Q - Qg) / seQ)
Dtab <- expand.grid(j = 1:3, c = 1:2)
Dtab$D_hat <- cvB$D[cbind(Dtab$j + 1, Dtab$c)]
Dtab$D_pop <- mapply(\(j, c) D_pop(LAM[j + 1])[c], Dtab$j, Dtab$c)
Dtab$mcse  <- cvB$D_mcse[cbind(Dtab$j + 1, Dtab$c)]
Dtab$z     <- (Dtab$D_hat - Dtab$D_pop) / Dtab$mcse
for (j in 0:3) record("B", sprintf("Q_hat(lambda=%.1f) within 4 SE of closed form", LAM[j + 1]),
                      abs(large$z[j + 1]), 0, 4)
for (i in seq_len(nrow(Dtab))) record("B", sprintf("D_hat[j=%d, c=%d] within 4 MCSE of closed form",
                                                   Dtab$j[i], Dtab$c[i]), abs(Dtab$z[i]), 0, 4)
write.csv(large, file.path(out, "large_n_grid.csv"), row.names = FALSE)
write.csv(Dtab, file.path(out, "large_n_gradient.csv"), row.names = FALSE)
nB <- cvB$n; rm(big, cvB, cvF, cvF2, cvX)   # keep the cluster export small

# ---- C. Monte Carlo coverage ------------------------------------------------------
mc_rep <- function(rep, n = 2000L, B = 10L) {
  set.seed(20261006L + rep)
  g <- gen_gd(n)
  fold <- .strat_folds(g$A, 5L)
  res <- list()
  for (cal in c("fixed", "fold")) {
    cv <- simex_fold_curve(g$X, g$A, g$dY, 1:2, GD$r, B = B, K = 5L, learner = "glm",
                           calibration = cal, s_fixed = GD$s0, fold = fold)
    p0 <- simex_fold_point(cv, -1, M = 0); pM <- simex_fold_point(cv, -1, M = M_TRUE)
    zq <- qnorm(0.975)
    res[[cal]] <- data.frame(rep = rep, calibration = cal, est = p0[["est"]],
      se_fixed = p0[["se_fixed"]], se = p0[["se"]],
      cover_q2_fixedse = abs(p0[["est"]] - q2_true(-1)) <= zq * p0[["se_fixed"]],
      cover_q2 = abs(p0[["est"]] - q2_true(-1)) <= zq * p0[["se"]],
      cover_tau_M0 = p0[["ci_lo"]] <= GD$tau & GD$tau <= p0[["ci_hi"]],
      cover_tau_Mtrue = pM[["ci_lo"]] <= GD$tau & GD$tau <= pM[["ci_hi"]],
      width_Mtrue = pM[["ci_hi"]] - pM[["ci_lo"]],
      band_covers_q2 = simex_fold_band_covers(cv, coef_true),
      clip_frac = max(cv$clip_frac, na.rm = TRUE))
  }
  do.call(rbind, res)
}
cl <- makeCluster(WORKERS)
clusterExport(cl, setdiff(ls(), c("cl", "checks", "record")))
invisible(clusterEvalQ(cl, { source("R_scripts/30_core.R"); source("R_scripts/34_simex_fold.R"); source("R_scripts/36_gauss_design.R"); NULL }))
t0 <- Sys.time()
mc <- do.call(rbind, parLapply(cl, seq_len(R_MC), mc_rep))
tC <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
stopCluster(cl)
write.csv(mc, file.path(out, "mc_closed_form_reps.csv"), row.names = FALSE)

wilson <- function(x) { n <- length(x); p <- mean(x); z <- qnorm(0.975)
  ctr <- (p + z^2 / (2 * n)) / (1 + z^2 / n); hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(ctr - hw, ctr + hw) }
summ <- do.call(rbind, lapply(split(mc, mc$calibration), \(z) {
  rows <- c("cover_q2_fixedse", "cover_q2", "cover_tau_M0", "cover_tau_Mtrue", "band_covers_q2")
  data.frame(calibration = z$calibration[1], metric = rows, n = nrow(z),
    value = vapply(rows, \(r) mean(z[[r]]), numeric(1)),
    mcse = vapply(rows, \(r) sqrt(mean(z[[r]]) * (1 - mean(z[[r]])) / nrow(z)), numeric(1)),
    wilson_lo = vapply(rows, \(r) wilson(z[[r]])[1], numeric(1)), wilson_hi = vapply(rows, \(r) wilson(z[[r]])[2], numeric(1)),
    row.names = NULL)
}))
se_tab <- do.call(rbind, lapply(split(mc, mc$calibration), \(z) data.frame(
  calibration = z$calibration[1], mean_est = mean(z$est), q2_target = q2_true(-1), tau = GD$tau,
  emp_sd = sd(z$est), mean_se_fixed = mean(z$se_fixed), mean_se_cal = mean(z$se),
  ratio_fixed = mean(z$se_fixed) / sd(z$est), ratio_cal = mean(z$se) / sd(z$est),
  max_clip_frac = max(z$clip_frac))))
write.csv(summ, file.path(out, "mc_closed_form_coverage.csv"), row.names = FALSE)
write.csv(se_tab, file.path(out, "mc_closed_form_se.csv"), row.names = FALSE)
# Pass rule for the intervals the theory says are valid: Wilson upper limit
# reaches 0.95 (no detectable undercoverage at this Monte Carlo size).
for (cal in c("fixed", "fold")) for (m in c("cover_q2", "cover_tau_Mtrue", "band_covers_q2")) {
  w <- summ[summ$calibration == cal & summ$metric == m, ]
  record("C", sprintf("%s calibration: %s, Wilson upper >= 0.95", cal, m), max(0.95 - w$wilson_hi, 0), 0, 0)
}

chk <- do.call(rbind, checks)
write.csv(chk, file.path(out, "checks.csv"), row.names = FALSE)
rep_txt <- capture.output({
  cat("Fold-calibrated SIMEX implementation checks,", SIMEX_FOLD_VERSION, "\n")
  cat("Run:", format(Sys.time()), "| R", R.version$major, ".", R.version$minor, "\n\n", sep = "")
  cat("Closed-form design: pi =", GD$pi, " dmu =", GD$dmu, " r =", GD$r, " s0 =", round(GD$s0, 5),
      " tau =", GD$tau, "\n")
  cat("q2(-1) =", q2_true(-1), " extrapolation bias q2(-1) - tau =", q2_true(-1) - GD$tau,
      " M_true =", M_TRUE, " 21M/22 =", 21 / 22 * M_TRUE, "\n\n")
  cat(sprintf("%d of %d checks pass\n\n", sum(chk$pass), nrow(chk)))
  print(chk, row.names = FALSE, digits = 4)
  cat(sprintf("\nB. Large sample n = %d, K = 2, B = 2 (%.0f s):\n", nB, tB))
  print(large, row.names = FALSE, digits = 5); print(Dtab, row.names = FALSE, digits = 5)
  cat(sprintf("\nC. Monte Carlo: %d reps, n = 2000, K = 5, B = 10, glm learner (%.1f min)\n", R_MC, tC))
  print(summ, row.names = FALSE, digits = 4); print(se_tab, row.names = FALSE, digits = 4)
  cat("\nThese check the implementation against its derivation in a correctly specified\n",
      "design. They are not evidence about lasso/boosting rates or the production DGP.\n", sep = "")
})
writeLines(rep_txt, file.path(out, "checks_report.txt"))
cat(rep_txt, sep = "\n")
