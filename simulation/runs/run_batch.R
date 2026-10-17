# One job = one (cell, batch of reps) in one R process. Design v2.0.
#
#   Rscript runs/run_batch.R --cell=<cell_id> --batch=<k> [--n_reps=500]
#          [--batch_size=25] [--out=sim_results_v2] [--estimators=a,b,c]
#
# Writes <out>/<cell_id>/batch_<k>.rds (a tibble, one row per estimator per
# rep) and is skipped by launch.sh if that file already exists. No nested
# parallelism: this process runs single-threaded (pin_threads()) and the
# shell driver decides how many jobs run at once.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(digest)
})
source("R_scripts/1_sim_param.R")
source("R_scripts/2_data_sim.R")
source("R_scripts/3_registry.R")
source("runs/grid.R")
pin_threads()

.kv <- function(args, key, default = NULL) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", key, "="), "", hit[[1]])
}

# Seed: replayable from (cell_id, rep) alone (PreReg s7).
seed_for <- function(cell_id, rep) {
  as.integer((digest::digest2int(cell_id) + rep) %% .Machine$integer.max)
}

# An estimator that ran but is no longer the estimator on its label says so
# in a "status" attribute (drdid_imp's glm fallback). A convergence warning
# nobody claimed is recorded generically. aggregate.R holds both out of the
# primary metrics and reports them as rates (PreReg v2.0 s5).
status_of <- function(err, warn, val, ok) {
  if (!is.null(err)) "error"
  else if (!ok) "NA_estimate"
  else if (!is.null(attr(val, "status"))) attr(val, "status")
  else if (!is.null(warn) &&
           grepl("converge", conditionMessage(warn), ignore.case = TRUE)) "non_convergence"
  else "ok"
}

run_one_rep <- function(rep, row, estimators) {
  seed_r <- seed_for(row$cell_id, rep)
  set.seed(seed_r)
  gen <- data_gen_v2(N = row$N, params = default_params(row$panel_strength),
                     dgp_family = row$dgp_family, seed = seed_r, return_cf = TRUE)
  dat <- gen$observed
  if (!is.na(row$fn_name)) {
    dat <- do.call(get(row$fn_name, mode = "function"), c(list(dat), row$args[[1]]))
  }
  d <- wide_frame(dat)
  out <- vector("list", length(estimators))
  for (i in seq_along(estimators)) {
    nm <- names(estimators)[i]
    set.seed(seed_r + i)
    tic <- Sys.time(); err <- NULL; warn <- NULL
    val <- tryCatch(
      withCallingHandlers(estimators[[i]](d),
                          warning = \(w) { warn <<- w; invokeRestart("muffleWarning") }),
      error = \(e) { err <<- e; NULL })
    rt <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
    ok <- is.null(err) && length(val) >= 4 && all(is.finite(val[c("est", "ci_lo", "ci_hi")]))
    status <- status_of(err, warn, val, ok)
    g <- \(k) if (ok || (is.null(err) && length(val) >= 4)) unname(val[[k]]) else NA_real_
    out[[i]] <- tibble(
      cell_id = row$cell_id, rep = rep, estimator = nm,
      est = g("est"), se = g("se"), ci_lo = g("ci_lo"), ci_hi = g("ci_hi"),
      lo = g("lo"), hi = g("hi"),
      true_atet = gen$true_atet, ols_pt_gap = gen$ols_pt_gap,
      status = status, runtime_s = rt,
      error_msg = if (!is.null(err)) substr(conditionMessage(err), 1, 200) else NA_character_)
  }
  bind_rows(out)
}

if (!interactive() && identical(sys.nframe(), 0L)) {
  args <- commandArgs(trailingOnly = TRUE)
  cell   <- .kv(args, "cell"); batch <- as.integer(.kv(args, "batch"))
  n_reps <- as.integer(.kv(args, "n_reps", "500"))
  bs     <- as.integer(.kv(args, "batch_size", "25"))
  out_dir <- .kv(args, "out", "sim_results_v2")
  only   <- .kv(args, "estimators")
  if (is.null(cell) || is.na(batch)) cli::cli_abort(c(
    "Both {.code --cell=<cell_id>} and {.code --batch=<k>} are required.",
    i = "Got --cell={.val {cell}}, --batch={.val {batch}}."
  ))

  grid <- build_grid()
  row  <- grid[grid$cell_id == cell, ]
  if (nrow(row) != 1L) {
    near <- grid$cell_id[order(utils::adist(cell, grid$cell_id))]
    cli::cli_abort(c(
      "Cell {.val {cell}} matched {nrow(row)} of the {nrow(grid)} grid cells, expected 1.",
      i = "Closest ids: {.val {near[seq_len(min(3L, length(near)))]}}."
    ))
  }
  ests <- estimator_set(row$me_design, row$args[[1]])
  if (!is.null(only)) ests <- ests[strsplit(only, ",")[[1]]]

  # seq() counts backwards when from > to, so an out-of-range --batch would
  # silently rerun earlier reps and save them under the wrong batch number.
  rep_from <- (batch - 1L) * bs + 1L
  rep_to   <- min(batch * bs, n_reps)
  if (rep_from > rep_to) cli::cli_abort(c(
    "Batch {batch} is empty: it covers reps {rep_from}-{rep_to} of n_reps={n_reps}.",
    i = "With batch_size={bs} the last non-empty batch is {ceiling(n_reps / bs)}."
  ))
  reps <- seq.int(rep_from, rep_to)

  dir.create(file.path(out_dir, cell), showWarnings = FALSE, recursive = TRUE)
  f <- file.path(out_dir, cell, sprintf("batch_%03d.rds", batch))
  t0 <- Sys.time()
  # One line per replicate, flushed. Without it this file writes nothing at
  # all until the batch ends, which at N = 5000 is up to three hours of an
  # empty log with no way to tell a slow job from a hung one. stdout is a
  # redirected file here, so it block-buffers unless flushed by hand.
  res <- vector("list", length(reps))
  for (k in seq_along(reps)) {
    tk <- Sys.time()
    res[[k]] <- run_one_rep(reps[k], row, ests)
    slow <- which.max(res[[k]]$runtime_s)
    cat(sprintf("  rep %d/%d (id %d) %.0fs, %d non-ok, slowest %s %.0fs\n",
                k, length(reps), reps[k],
                as.numeric(difftime(Sys.time(), tk, units = "secs")),
                sum(res[[k]]$status != "ok"),
                res[[k]]$estimator[slow], res[[k]]$runtime_s[slow]))
    flush(stdout())
  }
  res <- bind_rows(res)
  attr(res, "meta") <- list(cell_id = cell, batch = batch, reps = range(reps),
                            n_estimators = length(ests), started = t0, ended = Sys.time(),
                            r_version = R.version.string, platform = R.version$platform,
                            host = Sys.info()[["nodename"]])
  saveRDS(res, f)
  cat(sprintf("%s batch %d: %d reps, %d estimators, %.0fs, %d non-ok\n",
              cell, batch, length(reps), length(ests),
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              sum(res$status != "ok")))
}
