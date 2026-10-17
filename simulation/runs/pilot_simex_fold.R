# DIAGNOSTIC PILOT (not confirmatory) for the fold-calibrated SIMEX, lasso.
#   cd simulation; Rscript runs/pilot_simex_fold.R [--reps=40] [--workers=11]
#
# Post-production amendment 2026-10-06, P0 "revised SIMEX implementation and
# calibration": find implementation failures and time the method before a
# confirmatory design is frozen. Results here must not be reported as
# confirmatory coverage.
#
# Cell x_classical__r0.60__nonlinear. Data are regenerated with the production
# seeds, so each replicate is the same sample production analysed:
#   original  production simex_lasso, read from _per_rep.rds (replayed for the
#             first reps to confirm the replay is exact);
#   fixed     simex_fold, population error covariance supplied;
#   fold      simex_fold, calibration estimated on each training fold.
# New arms use stable name-based seeds, not list positions (Section 6).

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
REPS <- as.integer(kv("reps", "40")); WORKERS <- as.integer(kv("workers", "11"))
N_REPLAY <- 3L; N5000_REPS <- as.integer(kv("n5000", "2"))
FAMILY <- "nonlinear"; R_TRUE <- 0.60
M_GRID <- c(0, 0.005, 0.01, 0.02, 0.05)
out <- "amendment_2026-10-06/pilot"
dir.create(file.path(out, "curves"), recursive = TRUE, showWarnings = FALSE)

worker_init <- function() {
  suppressPackageStartupMessages({
    library(dplyr)
    library(tibble)
    library(tidyr)
    library(digest)
  })
  source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R")
  source("R_scripts/3_registry.R"); source("R_scripts/34_simex_fold.R")
  source("R_scripts/35_population_targets.R")
  pin_threads()
  NULL
}

task <- function(N, rep, replay) {
  cell <- sprintf("x_classical__r0.60__%s__N%d", FAMILY, N)
  seed_r <- as.integer((digest::digest2int(cell) + rep) %% .Machine$integer.max)  # run_batch.R
  set.seed(seed_r)
  gen <- data_gen_v2(N = N, params = default_params("low"), dgp_family = FAMILY,
                     seed = seed_r, return_cf = TRUE)
  d <- wide_frame(me_x_classical(gen$observed, reliability = R_TRUE))
  rows <- list()
  if (replay) {                       # simex_lasso is estimator 13 in x_classical cells
    set.seed(seed_r + 13L); tic <- Sys.time()
    v <- simex_drdid(d, dml_lasso, reliability = R_TRUE)
    rows$original <- tibble(arm = "original_replay", est = v[["est"]], se = v[["se"]], se_fixed = NA,
                            ci_lo = v[["ci_lo"]], ci_hi = v[["ci_hi"]],
                            runtime_s = as.numeric(difftime(Sys.time(), tic, units = "secs")))
  }
  for (cal in c("fixed", "fold")) {
    set.seed(as.integer(digest::digest2int(paste(cell, rep, "simex_fold_lasso", cal)) %% .Machine$integer.max))
    tic <- Sys.time()
    v <- simex_fold_drdid(d, R_TRUE, "glmnet", calibration = cal,
                          s_fixed = population_me_variance(R_TRUE))
    rt <- as.numeric(difftime(Sys.time(), tic, units = "secs"))
    cv <- attr(v, "curve")
    saveRDS(cv[setdiff(names(cv), "Z")], file.path(out, "curves", sprintf("N%d_rep%03d_%s.rds", N, rep, cal)))
    pts <- lapply(M_GRID, \(M) simex_fold_point(cv, -1, M))
    rows[[cal]] <- tibble(arm = cal, est = v[["est"]], se = v[["se"]], se_fixed = attr(v, "se_fixed"),
                          ci_lo = v[["ci_lo"]], ci_hi = v[["ci_hi"]], runtime_s = rt,
                          !!!setNames(lapply(pts, \(p) p[["ci_lo"]]), paste0("ci_lo_M", M_GRID)),
                          !!!setNames(lapply(pts, \(p) p[["ci_hi"]]), paste0("ci_hi_M", M_GRID)),
                          max_clip_frac = max(cv$clip_frac, na.rm = TRUE),
                          cD_h_share = var(drop(cv$h_s %*% t(cv$D) %*% .c_weights(cv$lam, -1))) /
                                       var(drop(cv$eta %*% .c_weights(cv$lam, -1))))
  }
  bind_rows(rows) |> mutate(N = N, rep = rep, true_atet = gen$true_atet,
                            pop_att = population_att(FAMILY), .before = 1)
}

tasks <- c(lapply(seq_len(REPS), \(r) list(N = 1000L, rep = r, replay = r <= N_REPLAY)),
           lapply(seq_len(N5000_REPS), \(r) list(N = 5000L, rep = r, replay = FALSE)))
tasks <- tasks[order(-vapply(tasks, \(t) t$N, 0L))]          # slow ones first
cl <- makeCluster(WORKERS)
clusterExport(cl, c("task", "worker_init", "FAMILY", "R_TRUE", "M_GRID", "out"))
invisible(clusterCall(cl, worker_init))
t0 <- Sys.time()
res <- bind_rows(parLapplyLB(cl, tasks, \(t) task(t$N, t$rep, t$replay)))
wall <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
stopCluster(cl)

source("R_scripts/1_sim_param.R"); source("R_scripts/35_population_targets.R")
# Production rows for the same replicates.
prod <- readRDS("sim_results_v2/_per_rep.rds") |>
  filter(cell_id == sprintf("x_classical__r0.60__%s__N1000", FAMILY),
         estimator == "simex_lasso", rep %in% seq_len(REPS)) |>
  transmute(N = 1000L, rep, arm = "original", est, se, ci_lo, ci_hi, runtime_s, true_atet,
            status)
pop <- population_att(FAMILY)
res <- bind_rows(res, prod |> mutate(pop_att = pop))
write.csv(res, file.path(out, "pilot_reps.csv"), row.names = FALSE)

replay <- res |> filter(arm == "original_replay") |>
  inner_join(prod |> select(rep, est_saved = est, se_saved = se), by = "rep") |>
  transmute(rep, est_diff = est - est_saved, se_diff = se - se_saved)

wil <- \(x) { n <- length(x); p <- mean(x); z <- qnorm(.975)
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)) }
summ <- res |> filter(N == 1000L, arm %in% c("original", "fixed", "fold")) |>
  group_by(arm) |>
  summarise(n = n(), bias_realised = mean(est - true_atet),
            bias_realised_mcse = sd(est - true_atet) / sqrt(n()),
            bias_population = mean(est - pop_att), emp_sd = sd(est), mean_se = mean(se),
            se_ratio = mean(se) / sd(est), mean_se_fixed = mean(se_fixed),
            coverage = mean(ci_lo <= true_atet & true_atet <= ci_hi),
            cov_wilson_lo = wil(ci_lo <= true_atet & true_atet <= ci_hi)[1],
            cov_wilson_hi = wil(ci_lo <= true_atet & true_atet <= ci_hi)[2],
            runtime_median_s = median(runtime_s), .groups = "drop")
mcov <- do.call(rbind, lapply(M_GRID, \(M) res |> filter(N == 1000L, arm %in% c("fixed", "fold")) |>
  group_by(arm) |>
  summarise(M = M, coverage = mean(.data[[paste0("ci_lo_M", M)]] <= true_atet &
                                     true_atet <= .data[[paste0("ci_hi_M", M)]]),
            mean_width = mean(.data[[paste0("ci_hi_M", M)]] - .data[[paste0("ci_lo_M", M)]]),
            .groups = "drop")))
timing <- res |> filter(arm %in% c("fixed", "fold", "original_replay")) |>
  group_by(N, arm) |> summarise(reps = n(), median_s = median(runtime_s), max_s = max(runtime_s), .groups = "drop")
diag <- res |> filter(arm %in% c("fixed", "fold")) |> group_by(N, arm) |>
  summarise(max_clip_frac = max(max_clip_frac), median_calibration_var_share = median(cD_h_share), .groups = "drop")
write.csv(summ, file.path(out, "pilot_summary.csv"), row.names = FALSE)
write.csv(mcov, file.path(out, "pilot_M_grid.csv"), row.names = FALSE)
write.csv(timing, file.path(out, "pilot_timing.csv"), row.names = FALSE)

txt <- capture.output({
  cat("DIAGNOSTIC PILOT, not confirmatory.", SIMEX_FOLD_VERSION <- "simex_fold_3.0", "\n")
  cat(sprintf("Cell x_classical__r0.60__%s; %d reps at N = 1000, %d timing reps at N = 5000; %.1f min wall on %d workers\n",
              FAMILY, REPS, N5000_REPS, wall, WORKERS))
  cat("Population ATT:", pop, "\n\nReplay of production simex_lasso (est/se differences):\n")
  print(as.data.frame(replay), digits = 4)
  cat("\nN = 1000 summary (coverage of the realised ATT; Wald, M = 0):\n"); print(as.data.frame(summ), digits = 4)
  cat("\nCoverage of the realised ATT against the exploratory M grid:\n"); print(as.data.frame(mcov), digits = 4)
  cat("\nTiming:\n"); print(as.data.frame(timing), digits = 4)
  cat("\nClipping and calibration-variance share of the q2(-1) influence:\n"); print(as.data.frame(diag), digits = 4)
  cat("\n", REPS, " replicates resolve coverage only to about +/- ", round(2 * sqrt(.95 * .05 / REPS), 3),
      ". Use this pilot to find failures and choose a budget; do not report it as validation.\n", sep = "")
})
writeLines(txt, file.path(out, "pilot_report.txt"))
cat(txt, sep = "\n")
