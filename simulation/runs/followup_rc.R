# Regression-calibration comparator (second amendment, Section 2;
# amendment_2026-10-08/AMENDMENT_2.md).
#   cd simulation
#   Rscript runs/followup_rc.R [--workers=5] [--reps=500] [--root=amendment_2026-10-08/rc] [--only=<study>]
# One task = (study, cell, learner, rep); one file per task under <root>/raw/, skipped
# if present. Every task fits, on the same folds: naive DR-DiD, DR-DiD on
# group-wise regression-calibrated covariates, and DR-DiD on the latent covariates;
# closed-form tasks add pooled calibration as a check.

suppressPackageStartupMessages({
  library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
WORKERS <- as.integer(kv("workers", "5")); REPS <- as.integer(kv("reps", "500"))
ROOT <- kv("root", "amendment_2026-10-08/rc"); ONLY <- kv("only", NA_character_)
dir.create(file.path(ROOT, "raw"), recursive = TRUE, showWarnings = FALSE)

xc <- expand.grid(r = c("0.60", "0.75", "0.90"), fam = c("linear", "nonlinear"), N = c(1000, 5000),
                  stringsAsFactors = FALSE)
cells <- rbind(
  data.frame(study = "rc_gauss", cell = c("gauss_r0.60_n5000", "gauss_r0.60_n1000", "gauss_r0.60_n1000"),
             learner = c("glm", "glm", "glmnet"), cost = c(3, 1, 10)),
  data.frame(study = "rc_dgp", cell = sprintf("x_classical__r%s__%s__N%d", xc$r, xc$fam, xc$N),
             learner = "glmnet", cost = dplyr::if_else(xc$N == 5000, 8, 3)),
  data.frame(study = "rc_xnc", cell = c("x_nc__large__linear__N5000", "x_nc__large__nonlinear__N5000"),
             learner = "glmnet", cost = 8))
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
  source("R_scripts/36_gauss_design.R"); source("R_scripts/38_regcal.R"); source("runs/grid.R")
  pin_threads()
  ROOT <<- root; GRID <<- build_grid()
  NULL
}
seed_for <- function(cell_id, rep) as.integer((digest::digest2int(cell_id) + rep) %% .Machine$integer.max)
seed_of  <- \(...) as.integer(digest::digest2int(paste(..., sep = "|")) %% .Machine$integer.max)

# gauss_draw() with the latent W returned as well; the RNG calls are identical,
# so X, A and dY equal gauss_draw()'s for the same seed (checked in the smoke run).
gauss_draw_latent <- function(gd, n) {
  A  <- rbinom(n, 1, gd$pi)
  W  <- vapply(seq_len(2), \(c) A * gd$dmu[c] + rnorm(n, 0, sqrt(gd$v[c] - gd$s0[c])), numeric(n))
  Ws <- W + vapply(seq_len(2), \(c) rnorm(n, 0, sqrt(gd$s0[c])), numeric(n))
  V  <- rnorm(n)
  dY <- drop(W %*% gd$beta) + gd$gamma * V + gd$tau * A + rnorm(n)
  list(X = cbind(w1 = Ws[, 1], w2 = Ws[, 2], v = V), XL = cbind(w1 = W[, 1], w2 = W[, 2], v = V), A = A, dY = dY)
}

fit_row <- function(method, Xm, A, dY, learner, fold, seed) {
  set.seed(seed)
  f <- fold_drdid(Xm, A, dY, learner, fold)
  se <- sqrt(drop(.if_vcov(f$psi)))
  data.frame(method, est = f$est, se, ci_lo = f$est - qnorm(0.975) * se, ci_hi = f$est + qnorm(0.975) * se,
             clip_frac = mean(f$e_hat < PS_TRIM[1] | f$e_hat > PS_TRIM[2]),
             n_flag = if (is.null(attr(Xm, "n_flag"))) 0L else attr(Xm, "n_flag"))
}

run_task <- function(study, cell, learner, rep) {
  f <- file.path(ROOT, "raw", study, paste0(cell, "__", learner), sprintf("rep_%03d.rds", rep))
  if (file.exists(f)) return(invisible("skip"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  rows <- tryCatch({
    if (study == "rc_gauss") {
      gd <- gauss_design(0.6, 0.6); n <- as.integer(sub(".*_n", "", cell))
      set.seed(seed_of(cell, rep, "data")); g <- gauss_draw_latent(gd, n)       # = first amendment's data
      set.seed(seed_of(cell, rep, "folds")); fold <- .strat_folds(g$A, 5L)
      X <- g$X; XL <- g$XL; A <- g$A; dY <- g$dY; r <- 0.6; me <- c("w1", "w2")
      true_att <- gd$tau; pop_att <- gd$tau
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
      A <- d$A; dY <- d$dY; me <- ME_X_COVS
      r <- if (row$me_design == "x_classical") row$args[[1]]$reliability else 0.75   # production SIMEX's value
      set.seed(seed_of(cell, rep, "folds")); fold <- .strat_folds(A, 5L)          # = first amendment's folds
      true_att <- gen$true_atet; pop_att <- population_att(row$dgp_family)
    }
    out <- list(fit_row("naive", X, A, dY, learner, fold, seed_of(cell, rep, learner, "naive")),
                fit_row("rc_group", regcal(X, A, me, r, TRUE), A, dY, learner, fold, seed_of(cell, rep, learner, "rc_group")),
                fit_row("latent", XL, A, dY, learner, fold, seed_of(cell, rep, learner, "latent")))
    if (study == "rc_gauss")
      out <- c(out, list(fit_row("rc_pooled", regcal(X, A, me, r, FALSE), A, dY, learner, fold,
                                 seed_of(cell, rep, learner, "rc_pooled"))))
    cbind(data.frame(study, cell, learner, rep, n = length(A), reliability = r, true_att, pop_att),
          do.call(rbind, out), error = NA_character_)
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
  clusterExport(cl, c("run_task", "fit_row", "gauss_draw_latent", "seed_for", "seed_of"))
  invisible(clusterCall(cl, worker_init, ROOT))
  t0 <- Sys.time()
  invisible(parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))),
                        \(t) run_task(t$study, t$cell, t$learner, t$rep), chunk.size = 1L))
  stopCluster(cl)
  cat(sprintf("Finished in %.2f h.\n", as.numeric(difftime(Sys.time(), t0, units = "hours"))))
}
