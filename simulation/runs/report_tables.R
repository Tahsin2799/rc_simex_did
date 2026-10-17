# Reporting tables for the frozen v2.0 production run (no estimator is refit).
#   cd simulation; Rscript runs/aggregate.R; Rscript runs/report_tables.R
# Reads sim_results_v2/_per_rep.rds and _summary.rds (from aggregate.R) and
# writes small CSVs to tables_v2/, which is tracked in git so the manuscript
# can cite them without the git-ignored batches.
#
#   metrics_all_cells.csv     every cell x estimator, clean and all-usable
#   paired_comparisons.csv    correction minus its naive counterpart on the SAME
#                             replicates, paired MCSE, clean and all-usable bases
#   range_estimators.csv      bounds_*: endpoints, width, set inclusion, interval
#   flags_failures.csv        non-ok outcomes by cell x estimator x status
#   population_targets.csv    realised versus population ATT by cell
#   headline_N5000.csv        post_prod_sim_steps.md Table 3.1 rows

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
src <- "sim_results_v2"; out <- "tables_v2"
dir.create(out, showWarnings = FALSE)
per_rep <- readRDS(file.path(src, "_per_rep.rds"))
summ    <- readRDS(file.path(src, "_summary.rds"))
if (!"bias_pop" %in% names(summ)) {      # aggregate.R must be the 2026-10-06 version
  cli::cli_abort("The summary has no {.field bias_pop} column; rerun {.file runs/aggregate.R} (2026-10-06 version or later).")
}
source("R_scripts/1_sim_param.R"); source("R_scripts/35_population_targets.R")

write.csv(summ, file.path(out, "metrics_all_cells.csv"), row.names = FALSE)

# Paired comparisons. Two estimators in one cell analyse the same generated
# sample, so their difference is summarised per replicate.
pairings <- tribble(~corrected, ~naive,
  "simex_lasso", "dml_lasso", "simex_xgb", "dml_xgb",
  "misclass_lasso", "dml_lasso", "misclass_xgb", "dml_xgb",
  "misclass_lasso_R", "dml_lasso", "misclass_xgb_R", "dml_xgb")
usable <- \(s) !(s %in% c("error", "NA_estimate"))
paired_one <- function(basis) {
  keep <- if (basis == "clean") \(s) s == "ok" else usable
  x <- per_rep |> filter(keep(status)) |>
    select(cell_id, rep, estimator, est, ci_lo, ci_hi, true_atet)
  pairings |>
    inner_join(x, by = c(corrected = "estimator"), relationship = "many-to-many") |>
    inner_join(x |> select(cell_id, rep, naive = estimator, est_n = est, lo_n = ci_lo, hi_n = ci_hi),
               by = c("cell_id", "rep", "naive")) |>
    mutate(d_est = est - est_n,
           d_abs = abs(est - true_atet) - abs(est_n - true_atet),
           d_sq  = (est - true_atet)^2 - (est_n - true_atet)^2,
           d_cov = (ci_lo <= true_atet & true_atet <= ci_hi) - (lo_n <= true_atet & true_atet <= hi_n)) |>
    group_by(cell_id, corrected, naive) |>
    summarise(basis = basis, n_pairs = n(),
              across(c(d_est, d_abs, d_sq, d_cov), list(mean = mean, mcse = \(v) sd(v) / sqrt(length(v)))),
              .groups = "drop")
}
paired <- bind_rows(paired_one("clean"), paired_one("all_usable")) |>
  rename(estimate_difference = d_est_mean, estimate_difference_mcse = d_est_mcse,
         abs_error_difference = d_abs_mean, abs_error_difference_mcse = d_abs_mcse,
         mse_difference = d_sq_mean, mse_difference_mcse = d_sq_mcse,
         coverage_difference = d_cov_mean, coverage_difference_mcse = d_cov_mcse)
write.csv(paired, file.path(out, "paired_comparisons.csv"), row.names = FALSE)

range_tab <- summ |> filter(is_range) |>
  select(cell_id, me_design, mag, dgp_family, N, estimator, n_used, lo_bias, hi_bias,
         set_width, set_width_mcse, set_inclusion, set_inclusion_all,
         interval_coverage = coverage, cov_wilson_lo, cov_wilson_hi, ci_width, ci_width_mcse,
         midpoint_bias, midpoint_rmse)
write.csv(range_tab, file.path(out, "range_estimators.csv"), row.names = FALSE)

flags <- per_rep |> filter(status != "ok") |> count(cell_id, estimator, status, name = "n") |>
  left_join(per_rep |> count(cell_id, estimator, name = "n_sim"), by = c("cell_id", "estimator")) |>
  mutate(rate = n / n_sim)
write.csv(flags, file.path(out, "flags_failures.csv"), row.names = FALSE)

targets <- per_rep |> distinct(cell_id, rep, true_atet) |>
  group_by(cell_id) |>
  summarise(realised_mean = mean(true_atet), realised_sd = sd(true_atet),
            realised_min = min(true_atet), realised_max = max(true_atet), .groups = "drop") |>
  mutate(population_att = if_else(grepl("__nonlinear__", cell_id),
                                 population_att("nonlinear"), population_att("linear")),
         realised_minus_population = realised_mean - population_att)
write.csv(targets, file.path(out, "population_targets.csv"), row.names = FALSE)

headline <- summ |>
  filter(N == 5000L, dgp_family == "nonlinear",
         (me_design == "baseline" & estimator %in% c("drdid_imp", "dml_lasso", "dml_xgb", "cf_did")) |
         (me_design == "x_classical" & mag == "r0.60" & estimator %in% c("dml_lasso", "simex_lasso", "dml_xgb", "simex_xgb")) |
         (me_design == "treat_misclass" & mag == "large" & estimator %in% c("dml_lasso", "misclass_lasso", "misclass_xgb")) |
         (me_design == "x_nc" & mag == "large" & estimator %in% c("dml_lasso", "simex_lasso"))) |>
  select(cell_id, estimator, n_used, bias, bias_mcse, bias_pop, rmse, rmse_mcse,
         coverage, cov_wilson_lo, cov_wilson_hi, ci_width, n_flag, n_fail)
write.csv(headline, file.path(out, "headline_N5000.csv"), row.names = FALSE)

cat(sprintf("tables_v2/: %d metric rows, %d paired rows, %d range rows, %d flag rows\n",
            nrow(summ), nrow(paired), nrow(range_tab), nrow(flags)))
cat("Population ATT: linear", population_att("linear"), " nonlinear", population_att("nonlinear"), "\n")
print(as.data.frame(headline |> mutate(across(where(is.numeric), \(x) signif(x, 4)))), row.names = FALSE)
