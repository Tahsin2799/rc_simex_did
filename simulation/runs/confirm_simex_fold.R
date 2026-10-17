# CONFIRMATORY supplemental study, post-production amendment 2026-10-06
# (amendment_2026-10-06/AMENDMENT.md, Section 4, frozen before this runs).
#   cd simulation
#   Rscript runs/confirm_simex_fold.R [--workers=11] [--reps=500] [--only=<study>]
#
# Resumable: one file per (study, cell, rep) under amendment_2026-10-06/confirm/raw/,
# skipped if present, so an interrupted run restarts where it stopped.
# Progress is appended to amendment_2026-10-06/confirm/_progress.log.
#
# Studies
#   gauss_known   closed-form design, r0 = r_L = 0.6; glm n = 1000, 5000; lasso n = 1000
#   gauss_range   closed-form design, r0 = 0.75, reference r_L = 0.6 (unknown
#                 reliability in [0.6, 1]); same learners and sizes
#   dgp_simex     production DGP, x_classical nonlinear: r0.60 N1000, r0.60 N5000,
#                 r0.75 N1000; lasso; production seeds, so samples equal production's
#   oracle        production DGP, x_classical r0.60 and treat_misclass large, linear
#                 and nonlinear, N5000: dml_lasso / dml_xgb on observed and latent data
#   bihs_semi     BIHS-calibrated semi-synthetic design (application/R_v2/
#                 03_semisynth_population.R): villages resampled, known ATT 0.25,
#                 reliability 0.7 on log assets and log land; cluster folds and
#                 village-clustered variance. Added 2026-10-06 before its first rep.
# Arms for the SIMEX studies: calibration "fixed" (population reference
# covariance) and "fold" (estimated per training fold), sharing folds and data.

suppressPackageStartupMessages({
  library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
WORKERS <- as.integer(kv("workers", "11")); REPS <- as.integer(kv("reps", "500"))
ONLY <- kv("only", NA_character_)
ROOT <- kv("root", "amendment_2026-10-06/confirm")      # --root= only for smoke tests
dir.create(file.path(ROOT, "raw"), recursive = TRUE, showWarnings = FALSE)
M_GRID <- c(0, 0.005, 0.01, 0.02, 0.05)
R_L <- 0.6
SAVE_IF_REPS <- 20L        # rowwise influence values kept for the first reps only

cells <- rbind(
  data.frame(study = "dgp_simex", cell = c("x_classical__r0.60__nonlinear__N5000",
             "x_classical__r0.60__nonlinear__N1000", "x_classical__r0.75__nonlinear__N1000"),
             learner = "glmnet", cost = c(290, 126, 126)),   # measured s per task, 2 arms
  data.frame(study = "gauss_known", cell = c("gauss_r0.60_n5000", "gauss_r0.60_n1000", "gauss_r0.60_n1000"),
             learner = c("glm", "glm", "glmnet"), cost = c(6, 2, 48)),
  data.frame(study = "gauss_range", cell = c("gauss_r0.75_n5000", "gauss_r0.75_n1000", "gauss_r0.75_n1000"),
             learner = c("glm", "glm", "glmnet"), cost = c(6, 2, 48)),
  data.frame(study = "oracle", cell = c("x_classical__r0.60__linear__N5000", "x_classical__r0.60__nonlinear__N5000",
             "treat_misclass__large__linear__N5000", "treat_misclass__large__nonlinear__N5000"),
             learner = "both", cost = 10),
  data.frame(study = "bihs_semi", cell = "bihs_semi_r0.70", learner = "glmnet", cost = 90))
if (!is.na(ONLY)) cells <- cells[cells$study == ONLY, ]

worker_init <- function(root) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(digest)
  })
  source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R")
  source("R_scripts/3_registry.R"); source("R_scripts/34_simex_fold.R")
  source("R_scripts/35_population_targets.R"); source("R_scripts/36_gauss_design.R")
  pin_threads()
  ROOT <<- root
  pop_file <- "../application/derived/semisynth_population.rds"
  if (file.exists(pop_file)) SS_POP <<- readRDS(pop_file)
  NULL
}

seed_of <- \(...) as.integer(digest::digest2int(paste(..., sep = "|")) %% .Machine$integer.max)
num <- \(x) as.numeric(x)

# Point and range summaries of one fitted curve.
summarise_curve <- function(cv, t_point, M_vals) {
  pts <- lapply(M_vals, \(M) simex_fold_point(cv, t_point, M))
  rgs <- lapply(M_vals, \(M) simex_fold_range(cv, -1, 0, M))
  p0 <- pts[[1]]; r0 <- rgs[[1]]
  c(est = p0[["est"]], se_fixed = p0[["se_fixed"]], se = p0[["se"]],
    range_lo = r0[["range_lo"]], range_hi = r0[["range_hi"]],
    setNames(vapply(pts, \(p) p[["ci_lo"]], 0), paste0("pt_lo_M", M_vals)),
    setNames(vapply(pts, \(p) p[["ci_hi"]], 0), paste0("pt_hi_M", M_vals)),
    setNames(vapply(rgs, \(r) r[["ci_lo"]], 0), paste0("rg_lo_M", M_vals)),
    setNames(vapply(rgs, \(r) r[["ci_hi"]], 0), paste0("rg_hi_M", M_vals)),
    rg_lo_cons_M0 = r0[["ci_lo_conservative"]], rg_hi_cons_M0 = r0[["ci_hi_conservative"]],
    max_clip_frac = max(cv$clip_frac, na.rm = TRUE),
    cal_var_share = var(drop(cv$h_s %*% t(cv$D) %*% .c_weights(cv$lam, t_point))) /
                    var(drop(cv$eta %*% .c_weights(cv$lam, t_point))))
}
slim_curve <- \(cv, keep_if) if (keep_if) cv[setdiff(names(cv), "Z")] else
  cv[setdiff(names(cv), c("Z", "zeta", "eta", "h_s", "fold"))]

run_task <- function(study, cell, learner, rep) {
  f <- file.path(ROOT, "raw", study, paste0(cell, "__", learner), sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) return(invisible("skip"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  out <- tryCatch(switch(study,
    gauss_known = , gauss_range = {
      r0 <- if (study == "gauss_known") 0.6 else 0.75
      n  <- as.integer(sub(".*_n", "", cell))
      gd <- gauss_design(r0, R_L); tr <- gauss_truth(gd)
      set.seed(seed_of(cell, rep, "data")); g <- gauss_draw(gd, n)
      set.seed(seed_of(cell, rep, "folds")); fold <- .strat_folds(g$A, 5L)
      rows <- list(); curves <- list()
      for (arm in c("fixed", "fold")) {
        set.seed(seed_of(cell, rep, learner, arm))
        cv <- simex_fold_curve(g$X, g$A, g$dY, 1:2, R_L, B = 20L, K = 5L, learner = learner,
                               calibration = arm, s_fixed = gd$sL, fold = fold)
        s <- summarise_curve(cv, gd$t0, c(M_GRID, tr$M))
        names(s) <- sub(paste0("M", tr$M, "$"), "Mtrue", names(s))
        rows[[arm]] <- data.frame(study, cell, learner, rep, arm, n, t_point = gd$t0,
          true_att = gd$tau, pop_att = gd$tau, q2_target = tr$q2(gd$t0),
          range_target_lo = tr$range[1], range_target_hi = tr$range[2], M_true = tr$M,
          band_covers = simex_fold_band_covers(cv, tr$coef), t(s))
        curves[[arm]] <- slim_curve(cv, rep <= SAVE_IF_REPS)
      }
      list(rows = do.call(rbind, rows), curves = curves)
    },
    dgp_simex = {
      r0 <- as.numeric(sub(".*__r([0-9.]+)__.*", "\\1", cell))
      fam <- strsplit(cell, "__")[[1]][3]; N <- as.integer(sub(".*__N", "", cell))
      seed_r <- seed_of_prod(cell, rep)
      set.seed(seed_r)
      gen <- data_gen_v2(N, default_params("low"), fam, seed = seed_r, return_cf = TRUE)
      d <- wide_frame(me_x_classical(gen$observed, reliability = r0))
      X <- design_matrix(d)
      set.seed(seed_of(cell, rep, "folds")); fold <- .strat_folds(d$A, 5L)
      t_point <- -(1 - r0) / (1 - R_L)
      s_ref <- population_me_variance(0.5) / r0 * (1 - R_L)   # Var(W) / r0 = Var(W*); times (1 - r_L)
      rows <- list(); curves <- list()
      for (arm in c("fixed", "fold")) {
        set.seed(seed_of(cell, rep, learner, arm))
        cv <- simex_fold_curve(X, d$A, d$dY, ME_X_COVS, R_L, B = 20L, K = 5L, learner = learner,
                               calibration = arm, s_fixed = s_ref, fold = fold)
        s <- summarise_curve(cv, t_point, M_GRID)
        rows[[arm]] <- data.frame(study, cell, learner, rep, arm, n = N, t_point,
          true_att = gen$true_atet, pop_att = population_att(fam), q2_target = NA_real_,
          range_target_lo = NA_real_, range_target_hi = NA_real_, M_true = NA_real_,
          band_covers = NA, t(s))
        curves[[arm]] <- slim_curve(cv, rep <= SAVE_IF_REPS)
      }
      list(rows = do.call(rbind, rows), curves = curves)
    },
    oracle = {
      parts <- strsplit(cell, "__")[[1]]
      fam <- parts[3]; N <- as.integer(sub("N", "", parts[4]))
      seed_r <- seed_of_prod(cell, rep)
      set.seed(seed_r)
      gen <- data_gen_v2(N, default_params("low"), fam, seed = seed_r, return_cf = TRUE)
      obs <- if (parts[1] == "x_classical") me_x_classical(gen$observed, reliability = 0.6)
             else me_treat_misclass(gen$observed, alpha0 = 0.10, alpha1 = 0.20)
      dl <- list(observed = wide_frame(obs), latent = wide_frame(gen$observed))
      rows <- list()
      for (est in c("dml_lasso", "dml_xgb")) for (ver in names(dl)) {
        set.seed(seed_of(cell, rep, est, ver))
        v <- get(est)(dl[[ver]])
        rows[[paste(est, ver)]] <- data.frame(study, cell, learner = est, rep, arm = ver, n = N,
          true_att = gen$true_atet, pop_att = population_att(fam),
          est = v[["est"]], se = v[["se"]], ci_lo = v[["ci_lo"]], ci_hi = v[["ci_hi"]])
      }
      list(rows = do.call(rbind, rows), curves = NULL)
    },
    bihs_semi = {
      pop <- SS_POP
      set.seed(seed_of(cell, rep, "data"))
      ug <- unique(pop$cluster); draw <- sample(ug, length(ug), replace = TRUE)
      idx <- unlist(lapply(draw, \(g) which(pop$cluster == g)))
      cl <- rep(seq_along(draw), times = vapply(draw, \(g) sum(pop$cluster == g), 0L))
      X <- pop$X[idx, , drop = FALSE]; n <- nrow(X)
      A <- rbinom(n, 1, plogis(drop(cbind(1, X) %*% pop$bA)))
      za <- (X[, pop$me[1]] - pop$mu_me[1]) / pop$sd_me[1]
      zl <- (X[, pop$me[2]] - pop$mu_me[2]) / pop$sd_me[2]
      dY <- drop(cbind(1, X) %*% pop$bY) + pop$nl * za * zl + pop$tau * A +
            rnorm(length(draw), 0, sqrt(pop$icc) * pop$sigma)[cl] +
            rnorm(n, 0, sqrt(1 - pop$icc) * pop$sigma)
      Xs <- X
      for (k in seq_along(pop$me)) Xs[, pop$me[k]] <- X[, pop$me[k]] + rnorm(n, 0, sqrt(pop$s0[k]))
      set.seed(seed_of(cell, rep, "folds")); fold <- .cluster_folds(cl, 5L)
      set.seed(seed_of(cell, rep, "oracle"))
      orc <- fold_drdid(X, A, dY, learner, fold)
      rows <- list(); curves <- list()
      for (arm in c("fixed", "fold")) {
        set.seed(seed_of(cell, rep, learner, arm))
        cv <- simex_fold_curve(Xs, A, dY, pop$me, pop$r0, B = 20L, K = 5L, learner = learner,
                               calibration = arm, s_fixed = pop$s0, fold = fold, cluster = cl)
        s <- summarise_curve(cv, -1, M_GRID)
        cw <- drop(.c_weights(cv$lam, -1))
        rows[[arm]] <- data.frame(study, cell, learner, rep, arm, n, t_point = -1,
          true_att = pop$tau, pop_att = pop$tau, q2_target = NA_real_,
          range_target_lo = NA_real_, range_target_hi = NA_real_, M_true = NA_real_,
          band_covers = NA, t(s),
          se_iid = sqrt(drop(t(cw) %*% .if_vcov(cv$eta, NULL) %*% cw)),
          naive_est = cv$Q[1], naive_se = sqrt(cv$V_zeta[1, 1]),
          naive_se_iid = sqrt(.if_vcov(cv$zeta[, 1, drop = FALSE], NULL)[1, 1]),
          oracle_est = orc$est, oracle_se = sqrt(drop(.if_vcov(orc$psi, cl))),
          clusters = length(draw))
        curves[[arm]] <- slim_curve(cv, rep <= SAVE_IF_REPS)
      }
      list(rows = do.call(rbind, rows), curves = curves)
    }),
    error = \(e) list(rows = data.frame(study, cell, learner, rep, error = conditionMessage(e)), curves = NULL))
  out$runtime_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out$method_version <- SIMEX_FOLD_VERSION
  saveRDS(out, f)
  cat(sprintf("%s %s %s %s rep %d %.0fs%s\n", format(Sys.time(), "%H:%M:%S"), study, cell, learner,
              rep, out$runtime_s, if ("error" %in% names(out$rows)) " ERROR" else ""),
      file = file.path(ROOT, "_progress.log"), append = TRUE)
  invisible("done")
}
# run_batch.R's replayable production seed (PreReg s7).
seed_of_prod <- \(cell, rep) as.integer((digest::digest2int(cell) + rep) %% .Machine$integer.max)

tasks <- do.call(rbind, lapply(seq_len(nrow(cells)), \(i)
  data.frame(cells[i, ], rep = seq_len(REPS), row.names = NULL)))
tasks <- tasks[order(-tasks$cost, tasks$rep), ]
done <- file.exists(file.path(ROOT, "raw", tasks$study, paste0(tasks$cell, "__", tasks$learner),
                              sprintf("rep_%03d.rds", tasks$rep)))
cat(sprintf("%d tasks, %d already done, %d to run on %d workers. Estimated %.1f h wall.\n",
            nrow(tasks), sum(done), sum(!done), WORKERS, sum(tasks$cost[!done]) / WORKERS / 3600))
tasks <- tasks[!done, ]
if (nrow(tasks)) {
  cl <- makeCluster(WORKERS)
  clusterExport(cl, c("run_task", "seed_of", "seed_of_prod", "summarise_curve", "slim_curve",
                      "M_GRID", "R_L", "SAVE_IF_REPS", "num"))
  invisible(clusterCall(cl, worker_init, ROOT))
  t0 <- Sys.time()
  # chunk.size = 1: tasks differ in cost by a factor of 150, and the default
  # (2 chunks per worker) hands one worker hundreds of the slowest tasks.
  invisible(parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))),
                        \(t) run_task(t$study, t$cell, t$learner, t$rep), chunk.size = 1L))
  stopCluster(cl)
  cat(sprintf("Finished in %.1f h. Next: Rscript runs/confirm_aggregate.R\n",
              as.numeric(difftime(Sys.time(), t0, units = "hours"))))
}
