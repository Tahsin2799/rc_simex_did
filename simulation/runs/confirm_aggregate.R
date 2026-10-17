# Aggregate the confirmatory supplemental study (runs/confirm_simex_fold.R).
#   cd simulation; Rscript runs/confirm_aggregate.R [--root=amendment_2026-10-06/confirm]
# Fast (reads small files). Writes to <root>/:
#   confirm_metrics_long.csv  one row per study x cell x learner x arm x metric
#   confirm_paired.csv        fold - fixed, and new arms - production, same reps
#   confirm_oracle.csv        observed versus latent (oracle) fits, paired gap
#   confirm_report.txt        human-readable summary

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
ROOT <- kv("root", "amendment_2026-10-06/confirm")
files <- list.files(file.path(ROOT, "raw"), "^rep_\\d+\\.rds$", recursive = TRUE, full.names = TRUE)
if (!length(files)) cli::cli_abort("No replicate files under {.path {file.path(ROOT, 'raw')}}.")
raw <- lapply(files, readRDS)
rows <- bind_rows(lapply(raw, \(x) mutate(x$rows, runtime_task_s = x$runtime_s)))
errs <- if ("error" %in% names(rows)) filter(rows, !is.na(error)) else rows[0, ]
rows <- if ("error" %in% names(rows)) filter(rows, is.na(error)) |> select(-error) else rows
z <- qnorm(0.975)
wilson <- function(x) { n <- length(x); p <- mean(x); den <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den) }
prop_row <- \(x, metric) { w <- wilson(x)
  tibble(metric = metric, value = mean(x), mcse = sqrt(mean(x) * (1 - mean(x)) / length(x)),
         wilson_lo = w[1], wilson_hi = w[2], n = length(x)) }
mean_row <- \(x, metric) tibble(metric = metric, value = mean(x), mcse = sd(x) / sqrt(length(x)),
                                wilson_lo = NA_real_, wilson_hi = NA_real_, n = length(x))

# ---- SIMEX studies ------------------------------------------------------------
sx <- rows |> filter(study != "oracle")
m_tags <- sub("^pt_lo_M", "", grep("^pt_lo_M", names(sx), value = TRUE))
simex_metrics <- function(g) {
  out <- list(
    mean_row(g$est - g$true_att, "bias_att"),
    mean_row(g$est - g$pop_att, "bias_pop_att"),
    tibble(metric = "emp_sd", value = sd(g$est), mcse = sd(g$est) / sqrt(2 * (nrow(g) - 1)),
           wilson_lo = NA, wilson_hi = NA, n = nrow(g)),
    mean_row(g$se, "mean_se"), mean_row(g$se_fixed, "mean_se_fixed"),
    tibble(metric = "se_ratio", value = mean(g$se) / sd(g$est), mcse = NA, wilson_lo = NA,
           wilson_hi = NA, n = nrow(g)),
    tibble(metric = "se_fixed_ratio", value = mean(g$se_fixed) / sd(g$est), mcse = NA,
           wilson_lo = NA, wilson_hi = NA, n = nrow(g)),
    mean_row(g$cal_var_share, "calibration_var_share"),
    tibble(metric = "max_clip_frac", value = max(g$max_clip_frac), mcse = NA, wilson_lo = NA,
           wilson_hi = NA, n = nrow(g)),
    mean_row(g$runtime_task_s / 2, "runtime_s_per_curve"),
    # added 2026-10-07 after the run (second amendment, Section 4; reporting only)
    tibble(metric = "runtime_task_median_s", value = median(g$runtime_task_s), mcse = NA, wilson_lo = NA,
           wilson_hi = NA, n = nrow(g)),
    tibble(metric = "runtime_task_p90_s", value = unname(quantile(g$runtime_task_s, 0.9)), mcse = NA,
           wilson_lo = NA, wilson_hi = NA, n = nrow(g)),
    tibble(metric = "median_clip_frac", value = median(g$max_clip_frac), mcse = NA, wilson_lo = NA,
           wilson_hi = NA, n = nrow(g)))
  for (m in m_tags) {
    lo <- g[[paste0("pt_lo_M", m)]]; hi <- g[[paste0("pt_hi_M", m)]]
    out <- c(out, list(prop_row(lo <= g$true_att & g$true_att <= hi, paste0("cover_att_point_M", m)),
                       mean_row(hi - lo, paste0("width_point_M", m))))
    lo <- g[[paste0("rg_lo_M", m)]]; hi <- g[[paste0("rg_hi_M", m)]]
    out <- c(out, list(prop_row(lo <= g$true_att & g$true_att <= hi, paste0("cover_att_outer_M", m)),
                       mean_row(hi - lo, paste0("width_outer_M", m))))
  }
  if (all(!is.na(g$q2_target))) {
    out <- c(out, list(
      mean_row(g$est - g$q2_target, "bias_q2"),
      prop_row(abs(g$est - g$q2_target) <= z * g$se, "cover_q2_calaware"),
      prop_row(abs(g$est - g$q2_target) <= z * g$se_fixed, "cover_q2_fixedse"),
      prop_row(g$band_covers, "cover_band_q2"),
      prop_row(g$rg_lo_M0 <= g$range_target_lo & g$range_target_hi <= g$rg_hi_M0, "cover_range_q2"),
      prop_row(g$rg_lo_cons_M0 <= g$range_target_lo & g$range_target_hi <= g$rg_hi_cons_M0,
               "cover_range_q2_conservative"),
      prop_row(g$range_lo <= g$true_att & g$true_att <= g$range_hi, "att_in_estimated_range")))
  } else out <- c(out, list(prop_row(g$range_lo <= g$true_att & g$true_att <= g$range_hi,
                                     "att_in_estimated_range")))
  if ("naive_est" %in% names(g) && all(!is.na(g$naive_est))) {      # bihs_semi only
    cv_ <- \(e, s) abs(e - g$true_att) <= z * s
    out <- c(out, list(
      prop_row(cv_(g$est, g$se_iid), "cover_att_point_M0_iidse"),
      mean_row(g$naive_est - g$true_att, "naive_bias_att"),
      prop_row(cv_(g$naive_est, g$naive_se), "naive_cover_att"),
      prop_row(cv_(g$naive_est, g$naive_se_iid), "naive_cover_att_iidse"),
      mean_row(g$oracle_est - g$true_att, "oracle_bias_att"),
      prop_row(cv_(g$oracle_est, g$oracle_se), "oracle_cover_att"),
      tibble(metric = "naive_emp_sd", value = sd(g$naive_est), mcse = NA, wilson_lo = NA,
             wilson_hi = NA, n = nrow(g)),
      mean_row(g$naive_se, "naive_mean_se")))
  }
  bind_rows(out)
}
long <- sx |> group_by(study, cell, learner, arm) |> group_modify(\(g, k) simex_metrics(g)) |> ungroup()

# Production comparators on the same replicates (dgp_simex cells share production seeds).
prod_file <- "sim_results_v2/_per_rep.rds"
prod <- if (file.exists(prod_file)) readRDS(prod_file) |>
  filter(cell_id %in% unique(sx$cell[sx$study == "dgp_simex"]),
         estimator %in% c("simex_lasso", "bounds_lasso", "dml_lasso"), status == "ok") |>
  select(cell = cell_id, rep, estimator, est, se, ci_lo, ci_hi, lo, hi, true_atet) else NULL
if (!is.null(prod) && nrow(prod)) {
  reps_used <- distinct(filter(sx, study == "dgp_simex"), cell, rep)
  pl <- prod |> semi_join(reps_used, by = c("cell", "rep")) |> group_by(cell, estimator) |>
    group_modify(\(g, k) bind_rows(
      mean_row(g$est - g$true_atet, "bias_att"),
      tibble(metric = "emp_sd", value = sd(g$est), mcse = NA, wilson_lo = NA, wilson_hi = NA, n = nrow(g)),
      mean_row(g$se, "mean_se"),
      prop_row(g$ci_lo <= g$true_atet & g$true_atet <= g$ci_hi, "cover_att_production_interval"),
      mean_row(g$ci_hi - g$ci_lo, "width_production_interval"))) |>
    ungroup() |> mutate(study = "dgp_simex", learner = "glmnet", arm = paste0("production_", estimator)) |>
    select(-estimator)
  long <- bind_rows(long, pl)
}
write.csv(long, file.path(ROOT, "confirm_metrics_long.csv"), row.names = FALSE)

# Compute actually used, per study and overall (second amendment, Section 4). Task
# files are written on completion, so their times bracket the run, restart included.
cmp <- tibble(study = basename(dirname(dirname(files))), runtime_s = vapply(raw, \(x) x$runtime_s, 0),
              finished = file.mtime(files))
cmp_sum <- \(d) summarise(d, tasks = n(), task_hours = sum(runtime_s) / 3600, median_s = median(runtime_s),
                          p90_s = unname(quantile(runtime_s, 0.9)), first_finished = min(finished),
                          last_finished = max(finished),
                          span_hours = as.numeric(difftime(max(finished), min(finished), units = "hours")))
write.csv(bind_rows(cmp |> group_by(study) |> cmp_sum() |> ungroup(), cmp_sum(cmp) |> mutate(study = "all")),
          file.path(ROOT, "confirm_compute.csv"), row.names = FALSE)

# ---- paired contrasts -----------------------------------------------------------
cov_pt0 <- \(d) as.numeric(d$pt_lo_M0 <= d$true_att & d$true_att <= d$pt_hi_M0)
pair_stats <- \(d, a, b, label) {
  d <- d |> filter(!is.na(est_a), !is.na(est_b))
  tibble(contrast = label, n = nrow(d),
         est_diff = mean(d$est_a - d$est_b), est_diff_mcse = sd(d$est_a - d$est_b) / sqrt(nrow(d)),
         abs_err_diff = mean(abs(d$est_a - d$true_att) - abs(d$est_b - d$true_att)),
         abs_err_diff_mcse = sd(abs(d$est_a - d$true_att) - abs(d$est_b - d$true_att)) / sqrt(nrow(d)),
         cover_diff = mean(d$cov_a - d$cov_b), cover_diff_mcse = sd(d$cov_a - d$cov_b) / sqrt(nrow(d)))
}
wide <- sx |> mutate(cov = cov_pt0(sx)) |> select(study, cell, learner, rep, arm, est, cov, true_att)
paired <- wide |> filter(arm == "fold") |>
  inner_join(wide |> filter(arm == "fixed") |> select(study, cell, learner, rep, est_b = est, cov_b = cov),
             by = c("study", "cell", "learner", "rep")) |>
  rename(est_a = est, cov_a = cov) |> group_by(study, cell, learner) |>
  group_modify(\(d, k) pair_stats(d, "fold", "fixed", "fold - fixed")) |> ungroup()
if (!is.null(prod) && nrow(prod)) {
  pp <- prod |> filter(estimator == "simex_lasso") |>
    transmute(cell, rep, est_b = est, cov_b = as.numeric(ci_lo <= true_atet & true_atet <= ci_hi))
  paired <- bind_rows(paired, wide |> filter(study == "dgp_simex") |> rename(est_a = est, cov_a = cov) |>
    inner_join(pp, by = c("cell", "rep")) |> group_by(study, cell, learner, arm) |>
    group_modify(\(d, k) pair_stats(d, k$arm, "production", paste(k$arm, "- production simex_lasso"))) |>
    ungroup() |> select(-arm))
}
# bihs_semi: SIMEX, naive and latent-data fits share each replicate (added 2026-10-07,
# after the run; reporting only, see the amendment's deviations log).
bs <- sx |> filter(study == "bihs_semi")
if (nrow(bs) && "naive_est" %in% names(bs)) {
  cov_se <- \(e, s, t) as.numeric(abs(e - t) <= z * s)
  fits <- list(SIMEX = transmute(bs, e = est, c = cov_pt0(bs)),
               naive = transmute(bs, e = naive_est, c = cov_se(naive_est, naive_se, true_att)),
               latent = transmute(bs, e = oracle_est, c = cov_se(oracle_est, oracle_se, true_att)))
  key_bs <- select(bs, study, cell, learner, arm, rep, true_att)
  paired <- bind_rows(paired, bind_rows(lapply(list(c("SIMEX", "naive"), c("SIMEX", "latent"),
                                                    c("naive", "latent")), \(p)
    bind_cols(key_bs, est_a = fits[[p[1]]]$e, cov_a = fits[[p[1]]]$c,
              est_b = fits[[p[2]]]$e, cov_b = fits[[p[2]]]$c) |>
      group_by(study, cell, learner, arm) |>
      group_modify(\(d, k) pair_stats(d, p[1], p[2], paste(p[1], "-", p[2]))) |> ungroup())))
}
write.csv(paired, file.path(ROOT, "confirm_paired.csv"), row.names = FALSE)

# ---- oracle -------------------------------------------------------------------------
orc <- rows |> filter(study == "oracle")
oracle <- NULL
if (nrow(orc)) {
  per <- orc |> group_by(cell, learner, arm) |>
    group_modify(\(g, k) bind_rows(mean_row(g$est - g$true_att, "bias_att"),
      tibble(metric = "rmse", value = sqrt(mean((g$est - g$true_att)^2)), mcse = NA,
             wilson_lo = NA, wilson_hi = NA, n = nrow(g)),
      prop_row(g$ci_lo <= g$true_att & g$true_att <= g$ci_hi, "cover_att"))) |> ungroup()
  gap <- orc |> select(cell, learner, rep, arm, est) |> pivot_wider(names_from = arm, values_from = est) |>
    group_by(cell, learner) |>
    summarise(metric = "observed_minus_latent", value = mean(observed - latent),
              mcse = sd(observed - latent) / sqrt(n()), wilson_lo = NA, wilson_hi = NA, n = n(),
              .groups = "drop") |> mutate(arm = "paired")
  oracle <- bind_rows(per, gap)
  write.csv(oracle, file.path(ROOT, "confirm_oracle.csv"), row.names = FALSE)
}

key <- c("bias_att", "bias_q2", "emp_sd", "mean_se", "se_ratio", "se_fixed_ratio", "cover_q2_calaware",
         "cover_q2_fixedse", "cover_band_q2", "cover_range_q2", "cover_att_point_M0",
         "cover_att_point_Mtrue", "width_point_Mtrue", "cover_att_outer_M0", "cover_att_outer_Mtrue",
         "cover_att_outer_M0.02", "att_in_estimated_range", "calibration_var_share",
         "cover_att_production_interval", "cover_att_point_M0_iidse", "naive_bias_att", "naive_cover_att",
         "oracle_bias_att", "oracle_cover_att")
txt <- capture.output({
  cat("Confirmatory supplemental study, aggregated", format(Sys.time()), "\n")
  cat(length(files), "task files;", nrow(errs), "errored tasks\n\n")
  if (nrow(errs)) print(as.data.frame(errs))
  w <- long |> filter(metric %in% key) |> mutate(v = sprintf("%.4f", value)) |>
    select(study, cell, learner, arm, metric, v) |> pivot_wider(names_from = metric, values_from = v)
  print(as.data.frame(w), right = FALSE)
  cat("\nPaired contrasts:\n"); print(as.data.frame(paired), digits = 3)
  if (!is.null(oracle)) { cat("\nOracle diagnostic:\n"); print(as.data.frame(oracle), digits = 3) }
})
writeLines(txt, file.path(ROOT, "confirm_report.txt"))
cat(txt, sep = "\n")
