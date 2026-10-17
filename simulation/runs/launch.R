# Portable flat-job driver, design v2.0. Run from simulation/.
#   Rscript runs/launch.R [--jobs=N] [--n_reps=500] [--batch_size=25]
#     [--out=sim_results_v2] [--timeout=21600] [--heartbeat=600] [--cell=REGEX]
#
# Same contract as runs/launch.sh: one single-threaded R process per
# (cell, batch), a batch whose .rds already exists is skipped, and re-running
# the same command resumes. This version replaces bash, `xargs -P` and
# `timeout` with processx, so it also runs on Windows.
#
# Two things it does that the shell driver does not, both from the
# 2026-05-20 post-mortem (SIM_ISSUES_2026-09-11.md):
#   I1  Workers start with --no-init-file and an explicit R_LIBS, so renv
#       never activates inside them. Putting six R processes through renv
#       activation at once is what corrupted the R heap on Windows, and
#       simulation/.Rprofile would do exactly that once per job.
#   I2  Every worker is held in a handle, killed when it exceeds --timeout,
#       and killed again on exit, so no job survives the run to contend with
#       the next relaunch.

suppressPackageStartupMessages({
  library(processx)
})
source("runs/grid.R")

.kv <- function(args, key, default) {
  hit <- grep(paste0("^--", key, "="), args, value = TRUE)
  if (!length(hit)) default else sub(paste0("^--", key, "="), "", hit[[1]])
}

args     <- commandArgs(trailingOnly = TRUE)
n_jobs   <- as.integer(.kv(args, "jobs", max(1L, parallel::detectCores(logical = FALSE) - 1L)))
n_reps   <- as.integer(.kv(args, "n_reps", "500"))
bs       <- as.integer(.kv(args, "batch_size", "25"))
out_dir  <- .kv(args, "out", "sim_results_v2")
timeout_s <- as.numeric(.kv(args, "timeout", as.character(6 * 3600)))
beat_s   <- as.numeric(.kv(args, "heartbeat", "600"))
cell_re  <- .kv(args, "cell", ".")

# --- Pending work ------------------------------------------------------------
jobs <- build_jobs(n_reps = n_reps, batch_size = bs)
f    <- file.path(out_dir, jobs$cell_id, sprintf("batch_%03d.rds", jobs$batch))
jobs <- jobs[!file.exists(f) & grepl(cell_re, jobs$cell_id), ]
dir.create(file.path(out_dir, "_logs"), showWarnings = FALSE, recursive = TRUE)
writeLines(paste(jobs$cell_id, jobs$batch), "runs/_pending.txt")
cat(nrow(jobs), "pending jobs\n")
if (!nrow(jobs)) quit(save = "no")

# --- Worker environment ------------------------------------------------------
# The renv library is resolved once, here, and handed to the workers as a
# plain R_LIBS. A worker therefore never reads renv.lock or runs the
# autoloader. Thread pinning must be in the environment before the child
# starts: a BLAS reads its thread count when it initialises, so setting it
# from inside a running R process is too late.
lib <- tryCatch(renv::paths$library(), error = \(e) .libPaths()[1])
if (!dir.exists(lib)) cli::cli_abort(c(
  "renv library not found at {.path {lib}}.",
  i = "Run {.code renv::restore()} from the repository root first."
))
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
worker_env <- c(Sys.getenv(),
                R_LIBS = lib, R_LIBS_USER = lib,
                RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE",
                OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

worker_args <- function(i) c(
  "--no-init-file", "runs/run_batch.R",
  paste0("--cell=", jobs$cell_id[i]), paste0("--batch=", jobs$batch[i]),
  paste0("--n_reps=", n_reps), paste0("--batch_size=", bs),
  paste0("--out=", out_dir))

log_of <- function(i) file.path(out_dir, "_logs",
                               sprintf("%s__b%d.log", jobs$cell_id[i], jobs$batch[i]))

# --- Pool --------------------------------------------------------------------
# In a function, so that on.exit() has a context to fire in: at script top
# level it silently never runs. processx is asked for cleanup = TRUE as well,
# which kills a worker if its handle is garbage-collected or R exits
# abnormally. Between them no worker outlives this driver (I2: the May run
# accumulated 21 stalled Rscript processes across three relaunches).
run_pool <- function(n_jobs, timeout_s, beat_s) {
  n        <- nrow(jobs)
  running  <- list()
  next_job <- 1L
  done     <- 0L
  progress <- file(file.path(out_dir, "_progress.log"), open = "at")
  on.exit({
    for (r in running) if (r$proc$is_alive()) r$proc$kill()
    close(progress)
  }, add = TRUE)

  say <- function(txt) { cat(txt, "\n", sep = ""); writeLines(txt, progress); flush(progress) }
  say(sprintf("== launch.R %s: %d jobs, %d at a time, timeout %.0fs, out=%s",
              format(Sys.time(), "%Y-%m-%d %H:%M:%S"), n, n_jobs, timeout_s, out_dir))
  last_beat <- Sys.time()

  # A finished job is the only thing that would otherwise reach this log, and
  # the heaviest batches run for hours. The heartbeat names the longest-running
  # workers so a stall is visible against their expected cost.
  beat <- function(done, running) {
    if (!length(running)) return(invisible(NULL))
    mins <- vapply(running, \(r) as.numeric(difftime(Sys.time(), r$started, units = "mins")), numeric(1))
    ord  <- order(mins, decreasing = TRUE)[seq_len(min(3L, length(running)))]
    say(sprintf("-- %s  %d/%d done, %d running; oldest: %s",
                format(Sys.time(), "%H:%M:%S"), done, n, length(running),
                paste(sprintf("%s b%d %.0fm", jobs$cell_id[vapply(running[ord], \(r) r$idx, integer(1))],
                              jobs$batch[vapply(running[ord], \(r) r$idx, integer(1))], mins[ord]),
                      collapse = ", ")))
  }

  while (done < n) {
    while (length(running) < n_jobs && next_job <= n) {
      i <- next_job
      running[[length(running) + 1L]] <- list(
        proc = process$new(rscript, worker_args(i), env = worker_env,
                           stdout = log_of(i), stderr = "2>&1", cleanup = TRUE),
        idx = i, started = Sys.time())
      next_job <- next_job + 1L
    }
    Sys.sleep(1)
    if (as.numeric(difftime(Sys.time(), last_beat, units = "secs")) >= beat_s) {
      beat(done, running); last_beat <- Sys.time()
    }
    keep <- logical(length(running))
    for (k in seq_along(running)) {
      r <- running[[k]]
      over <- as.numeric(difftime(Sys.time(), r$started, units = "secs")) > timeout_s
      if (over && r$proc$is_alive()) {
        r$proc$kill()
        say(sprintf("TIMEOUT %s %d after %.0fs (see %s)",
                    jobs$cell_id[r$idx], jobs$batch[r$idx], timeout_s, log_of(r$idx)))
      }
      if (r$proc$is_alive()) { keep[k] <- TRUE; next }
      done <- done + 1L
      st <- r$proc$get_exit_status()
      last <- tryCatch(tail(readLines(log_of(r$idx), warn = FALSE), 1), error = \(e) "")
      say(if (isTRUE(st == 0L)) sprintf("[%d/%d] %s", done, n, last)
          else sprintf("[%d/%d] FAILED %s %d (exit %s, see %s)", done, n,
                       jobs$cell_id[r$idx], jobs$batch[r$idx], st, log_of(r$idx)))
    }
    running <- running[keep]
  }
  say(sprintf("== done %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
}

run_pool(n_jobs, timeout_s, beat_s)
