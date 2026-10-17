# Aggregate batch outputs into per-cell summaries with MCSEs. Design v2.0.
#   Rscript runs/aggregate.R [--out=sim_results_v2]
# Writes <out>/_per_rep.rds (all rows) and <out>/_summary.rds (one row per
# cell x estimator: bias, empSE, modSE, RMSE, coverage, CI width, failure
# rate, mean runtime, bounds width; each with MCSE per Morris et al. 2019).
# Coverage is scored on (ci_lo, ci_hi), so Wald, posterior-quantile and
# Imbens-Manski intervals are all handled uniformly.
#
# Every rep falls in exactly one of three outcomes (PreReg v2.0 s5):
#   ok       clean run. The primary metrics are computed on these alone.
#   flagged  a usable estimate, but the solver complained or fell back to
#            something that is no longer the estimator on the label
#            (drdid_imp -> ipt_fallback). Reported as flag_rate and carried
#            in the _all sensitivity columns, never in the primary ones.
#   failed   no usable estimate. Reported as fail_rate.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})
args <- commandArgs(trailingOnly = TRUE)
out_dir <- sub("^--out=", "", grep("^--out=", args, value = TRUE)[1])
if (is.na(out_dir)) out_dir <- "sim_results_v2"

files <- list.files(out_dir, pattern = "^batch_\\d+\\.rds$", recursive = TRUE, full.names = TRUE)
if (!length(files)) cli::cli_abort(c(
  "No {.file batch_*.rds} found under {.path {out_dir}}.",
  i = "Pass {.code --out=<dir>}; the default is {.path sim_results_v2}."
))

# A batch truncated by a killed job fails in readRDS with a message that does not
# name the file. Say which one, so it can be deleted and that job rerun.
read_batch <- \(p) tryCatch(readRDS(p), error = \(e) cli::cli_abort(c(
  "Unreadable batch {.path {p}}.",
  i = "Most likely truncated by a killed job; delete it and rerun that batch.",
  x = "{conditionMessage(e)}"
), call = NULL))
per_rep <- bind_rows(lapply(files, read_batch))
saveRDS(per_rep, file.path(out_dir, "_per_rep.rds"))

per_rep <- per_rep |>
  mutate(failed  = status %in% c("error", "NA_estimate"),
         flagged = !failed & status != "ok")

# Post-production reporting repairs (2026-10-06; post_prod_sim_steps.md s3.2):
#   * bias MCSE from sd(est - true_atet): the nonlinear truth varies by rep;
#   * coverage carries a Wilson interval, so 0 % / 100 % coverage does not
#     report an MCSE of zero as if it were certain; denominators are the
#     intervals actually scored;
#   * CI width has an MCSE; runtime has median and 90th percentile;
#   * range estimators (bounds_*): est is only the midpoint of [lo, hi], so it
#     is kept out of bias/RMSE (midpoint_* instead) and the range is reported
#     by endpoint bias, width and set inclusion;
#   * bias against the population ATT (35_population_targets.R) beside the
#     registered realised-sample target.
source("R_scripts/1_sim_param.R"); source("R_scripts/35_population_targets.R")
RANGE_ESTIMATORS <- c("bounds_lasso", "bounds_xgb")
.wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n; ctr <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  cbind(lo = ctr - hw, hi = ctr + hw)
}
per_rep <- per_rep |>
  mutate(pop_att = if_else(grepl("__nonlinear__", cell_id),
                          population_att("nonlinear"), population_att("linear")),
         is_range = estimator %in% RANGE_ESTIMATORS)

# The same metric block on two bases: clean reps (primary) and every rep that
# returned a usable estimate (sensitivity, suffix _all). Morris et al. (2019)
# s5.4 ask for the failure mechanism and for the metrics computed both ways.
summarise_perf <- function(d, suffix = "") {
  out <- d |>
    mutate(err = est - true_atet,
           covered = ci_lo <= true_atet & true_atet <= ci_hi,
           in_set = lo <= true_atet & true_atet <= hi) |>
    group_by(cell_id, estimator) |>
    summarise(
      is_range    = first(is_range),
      n_used      = sum(!is.na(est)),
      bias        = mean(err, na.rm = TRUE),
      bias_mcse   = sd(err, na.rm = TRUE) / sqrt(n_used),
      bias_pop    = mean(est - pop_att, na.rm = TRUE),
      bias_pop_mcse = sd(est - pop_att, na.rm = TRUE) / sqrt(n_used),
      emp_se      = sd(est, na.rm = TRUE),
      emp_se_mcse = emp_se / sqrt(2 * (n_used - 1)),
      mod_se      = sqrt(mean(se^2, na.rm = TRUE)),
      rmse        = sqrt(mean(err^2, na.rm = TRUE)),
      rmse_mcse   = sqrt(var(err^2, na.rm = TRUE) / n_used) / (2 * rmse),
      n_cov       = sum(!is.na(covered)),
      coverage    = mean(covered, na.rm = TRUE),
      cov_mcse    = sqrt(coverage * (1 - coverage) / n_cov),
      cov_wilson_lo = .wilson(sum(covered, na.rm = TRUE), n_cov)[, "lo"],
      cov_wilson_hi = .wilson(sum(covered, na.rm = TRUE), n_cov)[, "hi"],
      ci_width    = mean(ci_hi - ci_lo, na.rm = TRUE),
      ci_width_mcse = sd(ci_hi - ci_lo, na.rm = TRUE) / sqrt(n_cov),
      set_width   = mean(hi - lo, na.rm = TRUE),
      set_width_mcse = sd(hi - lo, na.rm = TRUE) / sqrt(sum(!is.na(hi - lo))),
      set_inclusion = mean(in_set, na.rm = TRUE),
      lo_bias     = mean(lo - true_atet, na.rm = TRUE),
      hi_bias     = mean(hi - true_atet, na.rm = TRUE),
      runtime_s   = mean(runtime_s),
      runtime_median_s = median(runtime_s),
      runtime_p90_s = unname(quantile(runtime_s, 0.9)),
      .groups = "drop") |>
    mutate(midpoint_bias = if_else(is_range, bias, NA_real_),
           midpoint_rmse = if_else(is_range, rmse, NA_real_),
           across(c(bias, bias_mcse, bias_pop, bias_pop_mcse, rmse, rmse_mcse),
                  \(x) if_else(is_range, NA_real_, x)),
           across(c(set_width, set_width_mcse, set_inclusion, lo_bias, hi_bias),
                  \(x) if_else(is_range, x, NA_real_)))
  if (!nzchar(suffix)) out
  else rename_with(out, \(x) paste0(x, suffix), -c(cell_id, estimator))
}

counts <- per_rep |>
  group_by(cell_id, estimator) |>
  summarise(
    n_sim      = n(),
    n_ok       = sum(status == "ok"),
    n_flag     = sum(flagged),
    n_fail     = sum(failed),
    fail_rate  = mean(failed),
    fail_mcse  = sqrt(fail_rate * (1 - fail_rate) / n_sim),
    flag_rate  = mean(flagged),
    flag_mcse  = sqrt(flag_rate * (1 - flag_rate) / n_sim),
    ols_pt_gap = mean(ols_pt_gap),
    .groups = "drop")

summ <- counts |>
  left_join(summarise_perf(filter(per_rep, status == "ok")),
            by = c("cell_id", "estimator"), relationship = "one-to-one") |>
  left_join(summarise_perf(filter(per_rep, !failed), "_all") |>
              select(cell_id, estimator, n_used_all, bias_all, bias_mcse_all, rmse_all,
                     coverage_all, cov_wilson_lo_all, cov_wilson_hi_all, set_inclusion_all),
            by = c("cell_id", "estimator"), relationship = "one-to-one") |>
  separate(cell_id, c("me_design", "mag", "dgp_family", "N"), sep = "__", remove = FALSE) |>
  mutate(N = as.integer(sub("^N", "", N)))
saveRDS(summ, file.path(out_dir, "_summary.rds"))
cat(sprintf("%d rows per-rep, %d cell x estimator summaries, %d cells\n",
            nrow(per_rep), nrow(summ), n_distinct(summ$cell_id)))
rng <- summ |> count(cell_id) |> summarise(min_est = min(n), max_est = max(n))
cat(sprintf("estimators per cell: %d-%d\n", rng$min_est, rng$max_est))
st <- per_rep |> count(status) |> arrange(desc(n))
cat(sprintf("status: %s\n", paste(sprintf("%s=%d", st$status, st$n), collapse = "  ")))
