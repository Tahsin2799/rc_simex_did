# Checks of the fold-calibrated regression calibration and its calibration-aware
# influence function (third amendment). Writes amendment_2026-10-11/checks/rc_fold_checks.csv
# and rc_fold_mc.csv.
#   cd simulation; Rscript runs/rc_fold_checks.R [--mc_reps=400] [--workers=10]
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
MC_REPS <- as.integer(kv("mc_reps", "400")); WORKERS <- as.integer(kv("workers", "10"))
src <- function() {
  source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
  source("R_scripts/34_simex_fold.R"); source("R_scripts/36_gauss_design.R")
  source("R_scripts/38_regcal.R"); source("R_scripts/39_regcal_fold.R"); source("R_scripts/40_shape_design.R")
}
src()
OUT <- "amendment_2026-10-11/checks"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
res <- list(); add <- function(check, value, pass, note = "") res[[length(res) + 1]] <<- tibble(check, value, pass, note)

# 1. The map: the full-sample calibration equals 38_regcal.R; the preimage inverts the control map.
gd <- shape_design("gauss", "linear"); set.seed(101); g <- shape_draw(gd, 5000)
cp <- .rc_params(g$X, g$A, 1:2, 0.6)
x_new <- .rc_apply(g$X, g$A, cp); x_old <- regcal(g$X, g$A, c("w1", "w2"), 0.6, TRUE)
add("map_equals_38_regcal", max(abs(x_new - x_old)), max(abs(x_new - x_old)) < 1e-10)
X1 <- x_new[g$A == 1, ]; zp <- .rc_preimage(X1, cp)
back <- .rc_apply(zp, rep(0L, nrow(zp)), cp)
add("preimage_inverts_control_map", max(abs(back - X1)), max(abs(back - X1)) < 1e-10)

# 2. Vectorised contributions equal the generic directional derivative at IF directions.
set.seed(102); fold <- .strat_folds(g$A, 5L)
fit <- regcal_fold_drdid(g$X, g$A, g$dY, 1:2, 0.6, "glm", fold, "full")
pr <- fit$parts; Z <- g$X; me <- 1:2; n <- nrow(Z)
wbar <- colMeans(Z[, me]); v <- apply(Z[, me], 2, var)
if_dir <- function(j) {
  dmu <- dSig <- list()
  for (a in 0:1) {
    ga <- pr$cp$g[[a + 1]]; pia <- mean(g$A == a); u <- Z[j, ] - ga$mu
    dmu[[a + 1]] <- (g$A[j] == a) / pia * u
    dSig[[a + 1]] <- (g$A[j] == a) / pia * (tcrossprod(u) - ga$Sig)
  }
  list(dmu = dmu, dSig = dSig, ds = 0.4 * ((Z[j, me] - wbar)^2 - v))
}
js <- c(1, 2, 17, 450, 3333)
gen <- vapply(js, \(j) rc_dQ(pr, if_dir(j)), 0)
add("contrib_equals_generic_dQ", max(abs(gen - fit$c[js])), max(abs(gen - fit$c[js])) < 1e-10)

# 3. dQ against central differences, n = 200,000, parametric learner, full-sample calibration.
# (a) The outcome-regression part of the estimator, mean over treated of D - m(Xhat) with m
#     refitted at each calibration: its derivative is exactly the analytic dQ (calculus check).
# (b) The full DR-DiD estimator: it adds the propensity-weight path, whose population value is
#     zero when the outcome regression is correct (Gaussian design only); reported, not tested.
or_at <- function(X, A, dY, cp) {
  Xh <- .rc_apply(X, A, cp); mf <- .learner_glm(Xh[A == 0, , drop = FALSE], dY[A == 0], "or")
  mean(dY[A == 1] - mf(Xh[A == 1, , drop = FALSE]))
}
num_check <- function(shape, trend, seed, ndir = 4, eps = 1e-3) {
  gd <- shape_design(shape, trend); set.seed(seed); g <- shape_draw(gd, 200000)
  set.seed(seed + 1); fold <- .strat_folds(g$A, 5L)
  fit <- regcal_fold_drdid(g$X, g$A, g$dY, 1:2, 0.6, "glm", fold, "full")
  cp <- fit$parts$cp; p <- ncol(g$X)
  set.seed(seed + 2)
  bind_rows(lapply(seq_len(ndir), \(d) {
    rs <- \(k) { M <- matrix(rnorm(k * k), k); (M + t(M)) / 2 }
    delta <- list(dmu = list(rnorm(p), rnorm(p)), dSig = list(rs(p) * 0.2, rs(p) * 0.2), ds = rnorm(2) * 0.2)
    shift <- function(sgn) {
      c2 <- cp; c2$s <- cp$s + sgn * eps * delta$ds
      for (a in 0:1) {
        ga <- cp$g[[a + 1]]; Sig <- ga$Sig + sgn * eps * delta$dSig[[a + 1]]
        c2$g[[a + 1]] <- list(mu = ga$mu + sgn * eps * delta$dmu[[a + 1]], Sig = Sig, Si = solve(Sig),
                              s_eff = c2$s, n_flag = 0L)
      }
      c2
    }
    tibble(direction = d, analytic = rc_dQ(fit$parts, delta),
           numeric_or = (or_at(g$X, g$A, g$dY, shift(1)) - or_at(g$X, g$A, g$dY, shift(-1))) / (2 * eps),
           numeric_dr = (regcal_at(g$X, g$A, g$dY, shift(1), "glm", fold) -
                           regcal_at(g$X, g$A, g$dY, shift(-1), "glm", fold)) / (2 * eps))
  })) |> mutate(shape, trend)
}
nc <- bind_rows(num_check("gauss", "linear", 201), num_check("lnorm", "linear", 211), num_check("lnorm", "quad", 221))
rel_or <- with(nc, abs(numeric_or - analytic) / pmax(abs(numeric_or), 1e-3))
add("dQ_vs_central_differences_outcome_part", max(rel_or), max(rel_or) < 1e-3,
    sprintf("max relative error over %d directions, three designs; analytic G uses the fold fits, the numeric side one fit on all controls", nrow(nc)))
gl <- nc[nc$shape == "gauss", ]
add("dQ_vs_central_differences_full_DR_gauss", max(abs(gl$numeric_dr - gl$analytic)), NA,
    "diagnostic: difference is the propensity-weight path, mean zero when the outcome model is correct")
write.csv(nc, file.path(OUT, "rc_fold_dQ_numeric.csv"), row.names = FALSE)

# 4. Population check: fold-calibrated RC recovers tau in the closed-form design (n = 200,000).
gd <- shape_design("gauss", "linear"); set.seed(301); g <- shape_draw(gd, 200000); set.seed(302)
fit <- regcal_fold_drdid(g$X, g$A, g$dY, 1:2, 0.6, "glm", .strat_folds(g$A, 5L), "fold")
add("closed_form_large_n_bias", fit$est - gd$tau, abs(fit$est - gd$tau) < 4 * fit$se,
    sprintf("SE %.4f", fit$se))

# 5. Monte Carlo: SE calibration in the closed-form design, n = 2,000, parametric learner.
mc_one <- function(r) {
  gd <- shape_design("gauss", "linear"); set.seed(1e5 + r); g <- shape_draw(gd, 2000)
  set.seed(2e5 + r); fold <- .strat_folds(g$A, 5L)
  f <- regcal_fold_drdid(g$X, g$A, g$dY, 1:2, 0.6, "glm", fold, "fold")
  c(est = f$est, se = f$se, se_fixed = f$se_fixed, share = f$calib_share)
}
cl <- makeCluster(WORKERS); clusterCall(cl, src); clusterExport(cl, "mc_one")
mc <- do.call(rbind, parLapply(cl, seq_len(MC_REPS), mc_one)); stopCluster(cl)
mc <- as_tibble(mc); z <- qnorm(0.975)
mcs <- tibble(reps = nrow(mc), bias = mean(mc$est) - 0.5, bias_mcse = sd(mc$est) / sqrt(nrow(mc)), emp_sd = sd(mc$est),
              se_ratio = mean(mc$se) / sd(mc$est), se_ratio_fixed = mean(mc$se_fixed) / sd(mc$est),
              cover = mean(abs(mc$est - 0.5) <= z * mc$se), cover_fixed = mean(abs(mc$est - 0.5) <= z * mc$se_fixed),
              share = mean(mc$share))
write.csv(mcs, file.path(OUT, "rc_fold_mc.csv"), row.names = FALSE)
add("mc_se_ratio_calibration_aware", mcs$se_ratio, abs(mcs$se_ratio - 1) < 0.08,
    sprintf("%d reps, n = 2,000; fixed-calibration SE ratio %.3f", mcs$reps, mcs$se_ratio_fixed))
add("mc_coverage_calibration_aware", mcs$cover, abs(mcs$cover - 0.95) < 3 * sqrt(0.95 * 0.05 / mcs$reps),
    sprintf("fixed-calibration coverage %.3f", mcs$cover_fixed))

out <- bind_rows(res); write.csv(out, file.path(OUT, "rc_fold_checks.csv"), row.names = FALSE)
print(out, n = Inf, width = 200); print(mcs)
cat(if (all(out$pass, na.rm = TRUE)) "ALL PASS\n" else "FAILURE\n")
