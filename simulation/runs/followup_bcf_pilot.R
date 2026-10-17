# BCF diagnosis pilot (second amendment, Section 1; amendment_2026-10-08/AMENDMENT_2.md).
#   cd simulation
#   Rscript runs/followup_bcf_pilot.R [--workers=6] [--reps=20] [--root=amendment_2026-10-08/bcf]
# One task = (cell, rep, config, chain); one file per task under <root>/raw/, skipped
# if present, so the run is resumable. Each file holds the chain's draws of the
# treated-average tau. Data are production's: seed_for(cell, rep), data_gen_v2, then
# the cell's measurement-error function, exactly as runs/run_batch.R.

suppressPackageStartupMessages({
  library(parallel)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
WORKERS <- as.integer(kv("workers", "6")); REPS <- as.integer(kv("reps", "20"))
ROOT <- kv("root", "amendment_2026-10-08/bcf"); NSIM <- as.integer(kv("nsim", "1000"))
dir.create(file.path(ROOT, "raw"), recursive = TRUE, showWarnings = FALSE)

CELLS <- c("baseline__none__linear__N5000", "y_diff__medium__linear__N5000", "y_diff__medium__nonlinear__N5000")
CONFIGS <- list(prod = 250L, long = 1000L)          # nburn; nsim = 1000 in both
CHAINS <- 4L

worker_init <- function(root, nsim) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(digest)
    library(bcf)
  })
  source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
  source("runs/grid.R")
  pin_threads()
  ROOT <<- root; NSIM <<- nsim; GRID <<- build_grid()
  NULL
}
seed_for <- function(cell_id, rep) as.integer((digest::digest2int(cell_id) + rep) %% .Machine$integer.max)
seed_of  <- \(...) as.integer(digest::digest2int(paste(..., sep = "|")) %% .Machine$integer.max)

prod_data <- function(cell, rep) {
  row <- GRID[GRID$cell_id == cell, ]
  seed_r <- seed_for(cell, rep); set.seed(seed_r)
  gen <- data_gen_v2(N = row$N, params = default_params(row$panel_strength),
                     dgp_family = row$dgp_family, seed = seed_r, return_cf = TRUE)
  dat <- gen$observed
  if (!is.na(row$fn_name)) dat <- do.call(get(row$fn_name, mode = "function"), c(list(dat), row$args[[1]]))
  list(d = wide_frame(dat), true_atet = gen$true_atet)
}

run_task <- function(cell, rep, config, chain) {
  f <- file.path(ROOT, "raw", cell, config, sprintf("rep_%03d_chain_%d.rds", rep, chain))
  if (file.exists(f)) return(invisible("skip"))
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  t0 <- Sys.time()
  out <- tryCatch({
    pd <- prod_data(cell, rep); d <- pd$d
    X  <- design_matrix(d)
    ps <- glm(d$A ~ X, family = binomial)
    pihat <- pmin(pmax(as.numeric(fitted(ps)), PS_TRIM[1]), PS_TRIM[2])
    invisible(capture.output(
      fit <- bcf(y = d$dY, z = d$A, x_control = X, x_moderate = X, pihat = pihat,
                 random_seed = seed_of(cell, rep, config, chain), n_chains = 1,
                 nburn = CONFIGS[[config]], nsim = NSIM, ntree_control = 100, ntree_moderate = 50,
                 update_interval = NSIM, no_output = TRUE, verbose = FALSE, n_threads = 1),
      type = "output"))
    list(draws = rowMeans(fit$tau[, d$A == 1, drop = FALSE]), true_atet = pd$true_atet, error = NA_character_)
  }, error = \(e) list(draws = NULL, true_atet = NA_real_, error = conditionMessage(e)))
  out <- c(list(cell = cell, rep = rep, config = config, chain = chain, nburn = CONFIGS[[config]],
                nsim = NSIM), out, list(runtime_s = as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  saveRDS(out, f)
  cat(sprintf("%s %s %s chain %d rep %d %.0fs%s\n", format(Sys.time(), "%m-%d %H:%M:%S"), cell, config, chain,
              rep, out$runtime_s, if (!is.na(out$error)) " ERROR" else ""),
      file = file.path(ROOT, "_progress.log"), append = TRUE)
  invisible("done")
}

tasks <- expand.grid(chain = seq_len(CHAINS), config = names(CONFIGS), rep = seq_len(REPS), cell = CELLS,
                     stringsAsFactors = FALSE)
tasks <- tasks[order(tasks$config != "long", tasks$rep), ]     # long chains first for load balance
# --reverse=1 runs the list from the other end, so a second runner started later on
# free cores meets the first in the middle; each task re-checks its file before starting.
if (kv("reverse", "0") == "1") tasks <- tasks[rev(seq_len(nrow(tasks))), ]
done <- file.exists(file.path(ROOT, "raw", tasks$cell, tasks$config,
                              sprintf("rep_%03d_chain_%d.rds", tasks$rep, tasks$chain)))
cat(sprintf("%d tasks, %d already done, %d to run on %d workers.\n", nrow(tasks), sum(done), sum(!done), WORKERS))
tasks <- tasks[!done, ]
if (nrow(tasks)) {
  cl <- makeCluster(WORKERS)
  clusterExport(cl, c("run_task", "prod_data", "seed_for", "seed_of", "CONFIGS"))
  invisible(clusterCall(cl, worker_init, ROOT, NSIM))
  t0 <- Sys.time()
  invisible(parLapplyLB(cl, split(tasks, seq_len(nrow(tasks))),
                        \(t) run_task(t$cell, t$rep, t$config, t$chain), chunk.size = 1L))
  stopCluster(cl)
  cat(sprintf("Finished in %.2f h.\n", as.numeric(difftime(Sys.time(), t0, units = "hours"))))
}
