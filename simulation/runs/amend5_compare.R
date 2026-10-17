# Fifth amendment (amendment_2026-10-14/AMENDMENT_5.md, Section 3): paired comparison of the clamped
# spline learner (glmnet_nsc) with the fourth amendment's glmnet_ns rows, replicate by replicate
# (same data, folds and method seeds), plus robust summaries.
#   cd simulation; Rscript runs/amend5_compare.R [--root=amendment_2026-10-14/study]
# Writes <root>/amend5_compare.csv and amend5_compare.txt.
suppressPackageStartupMessages({
  library(dplyr)
})
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
ROOT <- kv("root", "amendment_2026-10-14/study"); OLD <- "amendment_2026-10-11b/study/raw/shape4"
rd <- \(dir) bind_rows(lapply(list.files(dir, "[.]rds$", recursive = TRUE, full.names = TRUE), \(p) readRDS(p)$rows))
new <- rd(file.path(ROOT, "raw")) |> mutate(err = est - true_att)
cells <- unique(new$cell)
old <- bind_rows(lapply(cells, \(ce) rd(file.path(OLD, paste0(ce, "__glmnet_ns"))))) |> mutate(err = est - true_att)
if (any(!is.na(c(new$error, old$error)))) {
  cli::cli_abort("{sum(!is.na(new$error))} new and {sum(!is.na(old$error))} old task{?s} ended in an error; the comparison needs complete pairs.")
}
cov <- \(d) d$ci_lo <= d$true_att & d$true_att <= d$ci_hi
pr <- inner_join(old |> transmute(cell, rep, method, err_old = err, cov_old = cov(old), se_old = se),
                 new |> transmute(cell, rep, method, err_new = err, cov_new = cov(new), se_new = se),
                 by = c("cell", "rep", "method"))
if (nrow(pr) != nrow(new)) cli::cli_abort("Pairing matched {nrow(pr)} of {nrow(new)} new rows to the fourth amendment's rows.")
out <- pr |> group_by(cell, method) |> summarise(
  n = n(),
  mean_err_old = mean(err_old), mean_err_new = mean(err_new),
  median_err_old = median(err_old), median_err_new = median(err_new),
  mad_new = mad(err_new), rmse_old = sqrt(mean(err_old^2)), rmse_new = sqrt(mean(err_new^2)),
  mean_change = mean(err_new - err_old), mcse_change = sd(err_new - err_old) / sqrt(n()),
  median_change = median(err_new - err_old), identical_reps = sum(abs(err_new - err_old) < 1e-12),
  n_big_old = sum(abs(err_old) > 1), n_big_new = sum(abs(err_new) > 1),
  big_fixed = sum(abs(err_old) > 1 & abs(err_new) <= 1), big_created = sum(abs(err_old) <= 1 & abs(err_new) > 1),
  cover_old = mean(cov_old), cover_new = mean(cov_new),
  se_sd_old = mean(se_old) / sd(err_old), se_sd_new = mean(se_new) / sd(err_new), .groups = "drop")
write.csv(out, file.path(ROOT, "amend5_compare.csv"), row.names = FALSE)
txt <- capture.output({
  cat("Fifth amendment: clamped spline (new) vs fourth amendment spline (old), paired by replicate,", format(Sys.time()), "\n\n")
  options(width = 220); print(as.data.frame(out), digits = 3)
})
writeLines(txt, file.path(ROOT, "amend5_compare.txt")); cat(txt, sep = "\n")
