# Fourth amendment (amendment_2026-10-11b/AMENDMENT_4.md): a design where regression
# calibration should fail -- skewed covariates (log-normal, omega = 1.0) at reliability 0.9 --
# with the quadratic-dictionary lasso and a natural-spline lasso. Same methods as the third
# amendment (runs/amend3_study.R), from which this runner is copied.
#   cd simulation
#   Rscript runs/amend4_study.R [--workers=10] [--reps=500] [--root=amendment_2026-10-11b/study]
# One task = (study, cell, learner, rep); one file per task under <root>/raw/, skipped if present.
#   shape: naive, latent, rc_group (full-sample calibration), rc_fold (fold calibration,
#          calibration-aware SE) and simex_fold (fold calibration, M in {0, 0.05}), same folds.
#   rc_dgp: naive, latent, rc_group, rc_fold on the production design, production seeds and
#          the second amendment's folds and learner seeds, so rc_group reproduces its rows.

suppressPackageStartupMessages({
  library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
WORKERS <- as.integer(kv("workers", "10")); REPS <- as.integer(kv("reps", "500"))
ROOT <- kv("root", "amendment_2026-10-11b/study"); ONLY <- kv("only", NA_character_)
dir.create(file.path(ROOT, "raw"), recursive = TRUE, showWarnings = FALSE)

cells <- data.frame(
  study = "shape4",
  cell = c("lnorm_linear_r090_n5000", "lnorm_quad_r090_n5000", "lnorm_linear_r090_n5000", "lnorm_quad_r090_n5000",
           "gauss_linear_r090_n5000", "lnorm_linear_r090_n1000", "lnorm_linear_r090_n1000"),
  learner = c("glmnet", "glmnet", "glmnet_ns", "glmnet_ns", "glmnet", "glmnet", "glmnet_ns"),
  cost = c(120, 120, 120, 120, 120, 40, 40))
if (!is.na(ONLY)) cells <- cells[cells$study == ONLY, ]

worker_init <- function(root) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(digest)
  })
  source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
  source("R_scripts/34_simex_fold.R"); source("R_scripts/35_population_targets.R")
  source("R_scripts/36_gauss_design.R"); source("R_scripts/38_regcal.R")
  source("R_scripts/39_regcal_fold.R"); source("R_scripts/40_shape_design.R"); source("R_scripts/41_spline_learner.R")
  source("runs/grid.R")
  pin_threads()
  ROOT <<- root; GRID <<- build_grid()
  NULL
}
seed_for <- function(cell_id, rep) as.integer((digest::digest2int(cell_id) + rep) %% .Machine$integer.max)
seed_of  <- \(...) as.integer(digest::digest2int(paste(..., sep = "|")) %% .Machine$integer.max)
M_GRID <- c(0, 0.05)

row_of <- function(method, est, se, se_fixed, true_att, extra = list()) {
  z <- qnorm(0.975)
  data.frame(method, est, se, se_fixed, ci_lo = est - z * se, ci_hi = est + z * se, true_att,
             calib_share = if (is.null(extra$calib_share)) NA_real_ else extra$calib_share,
             clip_frac = if (is.null(extra$clip_frac)) NA_real_ else extra$clip_frac,
             n_flag = if (is.null(extra$n_flag)) 0L else extra$n_flag,
             ci_lo_M005 = if (is.null(extra$ci_lo_M005)) NA_real_ else extra$ci_lo_M005,
             ci_hi_M005 = if (is.null(extra$ci_hi_M005)) NA_real_ else extra$ci_hi_M005,
             runtime_s = if (is.null(extra$runtime_s)) NA_real_ else extra$runtime_s)
}
timed <- function(expr) { t0 <- Sys.time(); v <- force(expr); attr(v, "rt") <- as.numeric(difftime(Sys.time(), t0, units = "secs")); v }

run_task <- function(study, cell, learner, rep) {
  f <- file.path(ROOT, "raw", study, paste0(cell, "__", learner), sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) return(invisible("skip"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  rows <- tryCatch({
    if (study == "shape4") {
      p <- strsplit(cell, "_")[[1]]; r0 <- as.numeric(sub("r", "", p[3])) / 100
      gd <- shape_design(p[1], p[2], r0, omega = 1.0); n <- as.integer(sub("n", "", p[4]))
      set.seed(seed_of(cell, rep, "data")); g <- shape_draw(gd, n)
      set.seed(seed_of(cell, rep, "folds")); fold <- .strat_folds(g$A, 5L)
      X <- g$X; XL <- g$XL; A <- g$A; dY <- g$dY; r <- gd$r0; me <- c("w1", "w2"); true_att <- gd$tau
    } else {
      row <- GRID[GRID$cell_id == cell, ]
      seed_r <- seed_for(cell, rep); set.seed(seed_r)
      gen <- data_gen_v2(N = row$N, params = default_params(row$panel_strength),
                         dgp_family = row$dgp_family, seed = seed_r, return_cf = TRUE)
      d  <- wide_frame(do.call(get(row$fn_name, mode = "function"), c(list(gen$observed), row$args[[1]])))
      dl <- wide_frame(gen$observed)
      if (!identical(d$id, dl$id) || !identical(d$A, dl$A)) {
        cli::cli_abort("Cell {.val {cell}} rep {rep}: the measurement wrapper reordered units or changed A, so the error-prone and true covariates no longer align.")
      }
      X <- design_matrix(d); XL <- design_matrix(dl)
      common <- intersect(colnames(X), colnames(XL)); X <- X[, common]; XL <- XL[, common]
      A <- d$A; dY <- d$dY; me <- ME_X_COVS; r <- row$args[[1]]$reliability
      set.seed(seed_of(cell, rep, "folds")); fold <- .strat_folds(A, 5L)       # = second amendment's folds
      true_att <- gen$true_atet
    }
    se1 <- \(psi) sqrt(drop(.if_vcov(psi)))
    set.seed(seed_of(cell, rep, learner, "naive"));  nv <- timed(fold_drdid(X, A, dY, learner, fold))
    set.seed(seed_of(cell, rep, learner, "latent")); lt <- timed(fold_drdid(XL, A, dY, learner, fold))
    set.seed(seed_of(cell, rep, learner, "rc_group")); rg <- timed(regcal_fold_drdid(X, A, dY, me, r, learner, fold, "full"))
    set.seed(seed_of(cell, rep, learner, "rc_fold"));  rf <- timed(regcal_fold_drdid(X, A, dY, me, r, learner, fold, "fold"))
    out <- list(
      row_of("naive", nv$est, se1(nv$psi), se1(nv$psi), true_att,
             list(clip_frac = mean(nv$e_hat < PS_TRIM[1] | nv$e_hat > PS_TRIM[2]), runtime_s = attr(nv, "rt"))),
      row_of("latent", lt$est, se1(lt$psi), se1(lt$psi), true_att,
             list(clip_frac = mean(lt$e_hat < PS_TRIM[1] | lt$e_hat > PS_TRIM[2]), runtime_s = attr(lt, "rt"))),
      row_of("rc_group", rg$est, rg$se, rg$se_fixed, true_att,
             list(calib_share = rg$calib_share, clip_frac = rg$clip_frac, n_flag = rg$n_flag, runtime_s = attr(rg, "rt"))),
      row_of("rc_fold", rf$est, rf$se, rf$se_fixed, true_att,
             list(calib_share = rf$calib_share, clip_frac = rf$clip_frac, n_flag = rf$n_flag, runtime_s = attr(rf, "rt"))))
    if (study == "shape4") {
      set.seed(seed_of(cell, rep, learner, "simex_fold"))
      cv <- timed(simex_fold_curve(X, A, dY, me, r, learner = learner, calibration = "fold", fold = fold))
      p0 <- simex_fold_point(cv, -1, 0); p5 <- simex_fold_point(cv, -1, 0.05)
      cw <- drop(.c_weights(cv$lam, -1))
      share <- var(drop(cv$h_s %*% t(cv$D) %*% cw)) / var(drop(cv$eta %*% cw))   # calibration term's variance share
      out <- c(out, list(row_of("simex_fold", p0[["est"]], p0[["se"]], p0[["se_fixed"]], true_att,
                                list(calib_share = share, clip_frac = mean(cv$clip_frac),
                                     ci_lo_M005 = p5[["ci_lo"]], ci_hi_M005 = p5[["ci_hi"]],
                                     runtime_s = attr(cv, "rt")))))
    }
    cbind(data.frame(study, cell, learner, rep, n = length(A), reliability = r), do.call(rbind, out), error = NA_character_)
  }, error = \(e) data.frame(study, cell, learner, rep, error = conditionMessage(e)))
  out <- list(rows = rows, runtime_s = as.numeric(difftime(Sys.time(), t0, units = "secs")))
  saveRDS(out, f)
  cat(sprintf("%s %s %s %s rep %d %.0fs%s\n", format(Sys.time(), "%m-%d %H:%M:%S"), study, cell, learner, rep,
              out$runtime_s, if (!is.na(rows$error[1])) " ERROR" else ""),
      file = file.path(ROOT, "_progress.log"), append = TRUE)
  invisible("done")
}

tasks <- do.call(rbind, lapply(seq_len(nrow(cells)), \(i) data.frame(cells[i, ], rep = seq_len(REPS), row.names = NULL)))
tasks <- tasks[order(-tasks$cost, tasks$rep), ]
done <- file.exists(file.path(ROOT, "raw", tasks$study, paste0(tasks$cell, "__", tasks$learner),
                              sprintf("rep_%03d.rds", tasks$rep)))
cat(sprintf("%d tasks, %d already done, %d to run on %d workers. Estimated %.1f h wall.\n",
            nrow(tasks), sum(done), sum(!done), WORKERS, sum(tasks$cost[!done]) / WORKERS / 3600))
tasks <- tasks[!done, ]
if (nrow(tasks)) {
  cl <- makeCluster(WORKERS)
  clusterExport(cl, c("run_task", "row_of", "timed", "seed_for", "seed_of", "M_GRID"))
  invisible(clusterCall(cl, worker_init, ROOT))
  t0 <- Sys.time()
  invisible(parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))),
                        \(t) run_task(t$study, t$cell, t$learner, t$rep), chunk.size = 1L))
  stopCluster(cl)
  cat(sprintf("Finished in %.2f h.\n", as.numeric(difftime(Sys.time(), t0, units = "hours"))))
}
