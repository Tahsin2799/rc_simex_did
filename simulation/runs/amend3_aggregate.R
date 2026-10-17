# Aggregate the third amendment's study (amendment_2026-10-11/AMENDMENT_3.md).
#   cd simulation; Rscript runs/amend3_aggregate.R [--root=amendment_2026-10-11/study] [--quiet=1]
# Writes <root>/amend3_metrics_long.csv, amend3_paired.csv and amend3_report.txt.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
ROOT <- kv("root", "amendment_2026-10-11/study"); QUIET <- kv("quiet", "0") == "1"
z <- qnorm(0.975)
wilson <- function(x) { n <- length(x); p <- mean(x); den <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den) }
r_mean <- \(x, metric) tibble(metric, value = mean(x), mcse = sd(x) / sqrt(length(x)), wilson_lo = NA_real_, wilson_hi = NA_real_, n = length(x))
r_prop <- \(x, metric) { w <- wilson(x); tibble(metric, value = mean(x), mcse = sqrt(mean(x) * (1 - mean(x)) / length(x)),
                                                 wilson_lo = w[1], wilson_hi = w[2], n = length(x)) }
r_val  <- \(v, metric, n) tibble(metric, value = v, mcse = NA_real_, wilson_lo = NA_real_, wilson_hi = NA_real_, n = n)

fs <- list.files(file.path(ROOT, "raw"), "\\.rds$", recursive = TRUE, full.names = TRUE)
obj <- lapply(fs, readRDS)
rows <- bind_rows(lapply(obj, `[[`, "rows"))
n_err <- sum(!is.na(rows$error)); rows <- filter(rows, is.na(error))
task_rt <- tibble(file = fs, task_s = vapply(obj, `[[`, 0, "runtime_s"))

rows <- rows |> mutate(err = est - true_att, covers = ci_lo <= true_att & true_att <= ci_hi,
                       covers_fixed = abs(err) <= z * se_fixed,
                       covers_M005 = if_else(method == "simex_fold", ci_lo_M005 <= true_att & true_att <= ci_hi_M005, NA))
key <- c("study", "cell", "learner", "method")
long <- rows |> group_by(across(all_of(key))) |> group_modify(\(d, k) {
  out <- bind_rows(
    r_mean(d$err, "bias_att"), r_val(sd(d$est), "emp_sd", nrow(d)),
    r_val(mean(d$se) / sd(d$est), "se_ratio", nrow(d)), r_val(mean(d$se_fixed) / sd(d$est), "se_ratio_fixed", nrow(d)),
    r_val(sqrt(mean(d$err^2)), "rmse", nrow(d)), r_prop(d$covers, "cover_att"), r_prop(d$covers_fixed, "cover_att_fixedse"),
    r_mean(d$ci_hi - d$ci_lo, "width"), r_val(median(d$clip_frac, na.rm = TRUE), "median_clip_frac", nrow(d)),
    r_val(max(d$clip_frac, na.rm = TRUE), "max_clip_frac", nrow(d)), r_val(sum(d$n_flag), "n_flag", nrow(d)),
    r_val(median(d$runtime_s), "runtime_median_s", nrow(d)))
  if (!all(is.na(d$calib_share))) out <- bind_rows(out, r_mean(d$calib_share, "calibration_var_share"))
  if (k$method == "simex_fold") out <- bind_rows(out, r_prop(d$covers_M005, "cover_att_M0.05"),
                                                 r_mean(d$ci_hi_M005 - d$ci_lo_M005, "width_M0.05"))
  out
}) |> ungroup()

contrasts <- list(c("rc_fold", "simex_fold"), c("rc_fold", "naive"), c("rc_fold", "latent"), c("rc_fold", "rc_group"),
                  c("simex_fold", "latent"), c("naive", "latent"), c("rc_group", "latent"))
paired <- bind_rows(lapply(contrasts, \(cc) {
  a <- rows |> filter(method == cc[1]) |> select(study, cell, learner, rep, est_a = est, err_a = err, cov_a = covers)
  b <- rows |> filter(method == cc[2]) |> select(study, cell, learner, rep, est_b = est, err_b = err, cov_b = covers)
  inner_join(a, b, by = c("study", "cell", "learner", "rep")) |> group_by(study, cell, learner) |>
    summarise(contrast = paste(cc[1], "-", cc[2]), n = n(),
              est_diff = mean(est_a - est_b), est_diff_mcse = sd(est_a - est_b) / sqrt(n()),
              abs_err_diff = mean(abs(err_a) - abs(err_b)), abs_err_diff_mcse = sd(abs(err_a) - abs(err_b)) / sqrt(n()),
              cover_diff = mean(cov_a - cov_b), cover_diff_mcse = sd(cov_a - cov_b) / sqrt(n()),
              rmse_a = sqrt(mean(err_a^2)), rmse_b = sqrt(mean(err_b^2)), .groups = "drop")
}))
write.csv(long, file.path(ROOT, "amend3_metrics_long.csv"), row.names = FALSE)
write.csv(paired, file.path(ROOT, "amend3_paired.csv"), row.names = FALSE)
compute <- tibble(tasks = length(fs), task_hours = sum(task_rt$task_s) / 3600, median_task_s = median(task_rt$task_s),
                  p90_task_s = unname(quantile(task_rt$task_s, 0.9)), errored = n_err)
write.csv(compute, file.path(ROOT, "amend3_compute.csv"), row.names = FALSE)
if (QUIET) {
  cat(sprintf("aggregated %d task files (%d errored): %d metric rows, %d paired rows, columns ok: %s\n", length(fs), n_err,
              nrow(long), nrow(paired), all(c("value", "mcse", "wilson_lo") %in% names(long))))
} else {
  wide <- long |> filter(metric %in% c("bias_att", "rmse", "se_ratio", "se_ratio_fixed", "cover_att", "cover_att_fixedse",
                                       "cover_att_M0.05", "calibration_var_share")) |>
    select(all_of(key), metric, value) |> pivot_wider(names_from = metric, values_from = value)
  txt <- capture.output({
    cat("Third amendment study,", format(Sys.time()), "-", length(fs), "task files,", n_err, "errored\n\n")
    print(as.data.frame(wide), digits = 3)
    cat("\nPaired:\n"); print(as.data.frame(paired), digits = 3)
    cat("\nCompute:\n"); print(as.data.frame(compute))
  })
  writeLines(txt, file.path(ROOT, "amend3_report.txt")); cat(txt, sep = "\n")
}
