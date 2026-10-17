# Aggregate the second amendment's runs (amendment_2026-10-08/AMENDMENT_2.md).
#   cd simulation; Rscript runs/followup_aggregate.R [--only=bcf|rc]
# Writes amendment_2026-10-08/bcf/{bcf_fits,bcf_summary,bcf_paired,bcf_trace_rep1}.csv
# and amendment_2026-10-08/rc/{rc_metrics_long,rc_paired}.csv, plus a report each.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
source("R_scripts/37_mcmc_diag.R")
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
ONLY <- kv("only", "all"); FU <- "amendment_2026-10-08"
z <- qnorm(0.975)
wilson <- function(x) { n <- length(x); p <- mean(x); den <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den) }
row_mean <- \(x, metric) tibble(metric, value = mean(x), mcse = sd(x) / sqrt(length(x)),
                                wilson_lo = NA_real_, wilson_hi = NA_real_, n = length(x))
row_prop <- \(x, metric) { w <- wilson(x); tibble(metric, value = mean(x), mcse = sqrt(mean(x) * (1 - mean(x)) / length(x)),
                                                 wilson_lo = w[1], wilson_hi = w[2], n = length(x)) }
row_val  <- \(v, metric, n) tibble(metric, value = v, mcse = NA_real_, wilson_lo = NA_real_, wilson_hi = NA_real_, n = n)
prod_rows <- readRDS("sim_results_v2/_per_rep.rds")

# ---- BCF pilot ---------------------------------------------------------------
if (ONLY %in% c("all", "bcf")) {
  fs <- list.files(file.path(FU, "bcf/raw"), "\\.rds$", recursive = TRUE, full.names = TRUE)
  ch <- lapply(fs, readRDS)
  meta <- bind_rows(lapply(ch, \(x) tibble(cell = x$cell, rep = x$rep, config = x$config, chain = x$chain,
                                           true_atet = x$true_atet, runtime_s = x$runtime_s, error = x$error)))
  if (any(!is.na(meta$error))) cli::cli_abort("{sum(!is.na(meta$error))} task{?s} ended in an error, e.g. {.val {meta$error[!is.na(meta$error)][1]}}.")
  key <- paste(meta$cell, meta$rep, meta$config)
  fits <- bind_rows(lapply(split(seq_along(ch), key), \(ix) {
    ix <- ix[order(meta$chain[ix])]
    if (length(ix) < 4) return(NULL)                                  # incomplete fit
    D <- vapply(ch[ix], `[[`, numeric(length(ch[[ix[1]]]$draws)), "draws")   # iterations x chains
    cm <- colMeans(D); m <- meta[ix[1], ]
    tibble(cell = m$cell, rep = m$rep, config = m$config, true_atet = m$true_atet, chains = ncol(D), draws = nrow(D),
           rhat = rhat_rank(D), ess_bulk = ess_bulk(D), ess_tail = ess_tail(D),
           post_mean = mean(D), post_sd = sd(as.vector(D)),
           chain1_mean = cm[1], chain1_sd = sd(D[, 1]), sd_chain_means = sd(cm),
           ci_lo = unname(quantile(D, 0.025)), ci_hi = unname(quantile(D, 0.975)),
           runtime_s = sum(meta$runtime_s[ix]))
  }))
  fits <- mutate(fits, covers = ci_lo <= true_atet & true_atet <= ci_hi,
                 chain1_covers = abs(chain1_mean - true_atet) <= z * chain1_sd)
  emp <- read.csv("tables_v2/metrics_all_cells.csv") |> filter(estimator == "bcf_did", cell_id %in% unique(fits$cell)) |>
    transmute(cell = cell_id, prod_emp_sd = emp_se, prod_mod_se = mod_se, prod_coverage = coverage)
  summ <- fits |> group_by(cell, config) |>
    summarise(fits = n(), rhat_median = median(rhat), rhat_max = max(rhat), share_rhat_lt_1.01 = mean(rhat < 1.01),
              ess_bulk_median = median(ess_bulk), ess_bulk_min = min(ess_bulk), share_bulk_gt_400 = mean(ess_bulk > 400),
              share_both = mean(rhat < 1.01 & ess_bulk > 400), ess_tail_median = median(ess_tail), ess_tail_min = min(ess_tail),
              mean_post_sd = mean(post_sd), mean_chain1_sd = mean(chain1_sd), mean_sd_chain_means = mean(sd_chain_means),
              bias = mean(post_mean - true_atet), bias_mcse = sd(post_mean - true_atet) / sqrt(n()),
              coverage = mean(covers), chain1_coverage = mean(chain1_covers), runtime_h = sum(runtime_s) / 3600,
              .groups = "drop") |>
    left_join(emp, by = "cell") |> mutate(post_sd_over_prod_emp_sd = mean_post_sd / prod_emp_sd)
  pw <- fits |> select(cell, rep, config, post_sd, post_mean) |>
    pivot_wider(names_from = config, values_from = c(post_sd, post_mean))
  for (v in c("post_sd_long", "post_sd_prod", "post_mean_long", "post_mean_prod")) if (!v %in% names(pw)) pw[[v]] <- NA_real_
  pw <- filter(pw, !is.na(post_sd_long), !is.na(post_sd_prod))
  paired <- pw |> group_by(cell) |>
    summarise(n = n(), sd_diff = mean(post_sd_long - post_sd_prod), sd_diff_mcse = sd(post_sd_long - post_sd_prod) / sqrt(n()),
              sd_pct_change = 100 * (mean(post_sd_long) / mean(post_sd_prod) - 1),
              mean_diff = mean(post_mean_long - post_mean_prod),
              mean_diff_mcse = sd(post_mean_long - post_mean_prod) / sqrt(n()), .groups = "drop") |>
    mutate(burnin_matters_sd = abs(sd_pct_change) >= 5 & abs(sd_diff) >= 2 * sd_diff_mcse)
  rule <- summ |> filter(config == "prod") |> select(cell, share_both) |> left_join(paired, by = "cell") |>
    mutate(mcmc_length_explains = !(share_both >= 0.9 & !burnin_matters_sd))
  write.csv(fits, file.path(FU, "bcf/bcf_fits.csv"), row.names = FALSE)
  write.csv(summ, file.path(FU, "bcf/bcf_summary.csv"), row.names = FALSE)
  write.csv(left_join(paired, select(rule, cell, share_both, mcmc_length_explains), by = "cell"),
            file.path(FU, "bcf/bcf_paired.csv"), row.names = FALSE)
  tr <- bind_rows(lapply(which(meta$rep == 1), \(i) tibble(cell = meta$cell[i], config = meta$config[i],
                                                           chain = meta$chain[i], iter = seq_along(ch[[i]]$draws), draw = ch[[i]]$draws)))
  write.csv(tr, file.path(FU, "bcf/bcf_trace_rep1.csv"), row.names = FALSE)
  txt <- capture.output({
    cat("BCF pilot,", format(Sys.time()), "-", nrow(meta), "chain files,", nrow(fits), "complete fits\n\n")
    print(as.data.frame(summ), digits = 3); cat("\nPaired long - prod and the reading rule:\n"); print(as.data.frame(rule), digits = 3)
  })
  writeLines(txt, file.path(FU, "bcf/bcf_report.txt")); cat(txt, sep = "\n")
}

# ---- regression calibration --------------------------------------------------
if (ONLY %in% c("all", "rc")) {
  fs <- list.files(file.path(FU, "rc/raw"), "^rep_\\d+\\.rds$", recursive = TRUE, full.names = TRUE)
  raw <- lapply(fs, readRDS)
  rows <- bind_rows(lapply(raw, \(x) mutate(x$rows, runtime_s = x$runtime_s)))
  errs <- filter(rows, !is.na(error)); rows <- filter(rows, is.na(error))
  mets <- \(g) bind_rows(
      row_mean(g$est - g$true_att, "bias_att"), row_mean(g$est - g$pop_att, "bias_pop_att"),
      row_val(sd(g$est), "emp_sd", nrow(g)), row_mean(g$se, "mean_se"),
      row_val(mean(g$se) / sd(g$est), "se_ratio", nrow(g)),
      row_val(sqrt(mean((g$est - g$true_att)^2)), "rmse", nrow(g)),
      row_prop(g$ci_lo <= g$true_att & g$true_att <= g$ci_hi, "cover_att"),
      row_mean(g$ci_hi - g$ci_lo, "width"))
  long <- rows |> group_by(study, cell, learner, method) |>
    group_modify(\(g, k) bind_rows(mets(g),
      row_val(median(g$clip_frac), "median_clip_frac", nrow(g)),
      row_val(max(g$clip_frac), "max_clip_frac", nrow(g)), row_val(sum(g$n_flag), "n_flag", nrow(g)),
      row_val(median(g$runtime_s), "runtime_task_median_s", nrow(g)))) |> ungroup()
  # SIMEX on the same replicates, scored the same way: production SIMEX-lasso
  # (production cells) and revised SIMEX, fold arm, M = 0 (first-amendment cells).
  simex_rows <- rows |> distinct(study, cell, learner) |> group_by(study, cell, learner) |>
    group_modify(\(g, k) {
      reps <- unique(rows$rep[rows$study == k$study & rows$cell == k$cell & rows$learner == k$learner])
      tru <- rows |> filter(study == k$study, cell == k$cell, learner == k$learner, method == "naive") |>
        select(rep, true_att, pop_att)
      out <- list()
      if (k$study != "rc_gauss")
        out$simex_production <- prod_rows |> filter(cell_id == k$cell, estimator == "simex_lasso", status == "ok", rep %in% reps) |>
          transmute(rep, est, se, ci_lo, ci_hi) |> inner_join(tru, by = "rep")
      src <- if (k$study == "rc_gauss") c("gauss_known", k$learner) else c("dgp_simex", "glmnet")
      f <- list.files(file.path("amendment_2026-10-06/confirm/raw", src[1], paste0(k$cell, "__", src[2])), full.names = TRUE)
      if (length(f))
        out$simex_fold <- bind_rows(lapply(f, \(x) readRDS(x)$rows)) |> filter(arm == "fold", rep %in% reps) |>
          transmute(rep, est, se, ci_lo = pt_lo_M0, ci_hi = pt_hi_M0) |> inner_join(tru, by = "rep")
      bind_rows(lapply(names(out), \(m) mutate(mets(out[[m]]), method = m)))
    }) |> ungroup()
  long <- bind_rows(long, simex_rows)

  # Paired comparators on the same replicates.
  pair <- function(a, b, label) {                       # a, b: rep, est, cov, true_att
    d <- inner_join(a, b, by = "rep", suffix = c("_a", "_b"))
    e <- d$est_a - d$est_b; ae <- abs(d$est_a - d$true_att_a) - abs(d$est_b - d$true_att_a); cdf <- d$cov_a - d$cov_b
    tibble(contrast = label, n = nrow(d), est_diff = mean(e), est_diff_mcse = sd(e) / sqrt(nrow(d)),
           abs_err_diff = mean(ae), abs_err_diff_mcse = sd(ae) / sqrt(nrow(d)),
           cover_diff = mean(cdf), cover_diff_mcse = sd(cdf) / sqrt(nrow(d)),
           rmse_a = sqrt(mean((d$est_a - d$true_att_a)^2)), rmse_b = sqrt(mean((d$est_b - d$true_att_a)^2)))
  }
  own <- \(g, m) g |> filter(method == m) |>
    transmute(rep, est, cov = as.numeric(ci_lo <= true_att & true_att <= ci_hi), true_att)
  cf_files <- "amendment_2026-10-06/confirm/raw"
  confirm_fold <- function(study, cell, learner) {
    f <- list.files(file.path(cf_files, study, paste0(cell, "__", learner)), full.names = TRUE)
    if (!length(f)) return(NULL)
    bind_rows(lapply(f, \(x) readRDS(x)$rows)) |> filter(arm == "fold") |>
      transmute(rep, est, cov = as.numeric(pt_lo_M0 <= true_att & true_att <= pt_hi_M0), true_att)
  }
  paired <- rows |> group_by(study, cell, learner) |> group_modify(\(g, k) {
    rc <- own(g, "rc_group")
    out <- list(pair(rc, own(g, "naive"), "RC - naive"), pair(rc, own(g, "latent"), "RC - latent"),
                pair(own(g, "naive"), own(g, "latent"), "naive - latent"))
    if (k$study == "rc_gauss") {
      sf <- confirm_fold("gauss_known", k$cell, k$learner)
      if (!is.null(sf)) out <- c(out, list(pair(rc, sf, "RC - SIMEX fold (M = 0)")))
    } else {
      ps <- prod_rows |> filter(cell_id == k$cell, estimator == "simex_lasso", status == "ok") |>
        transmute(rep, est, cov = as.numeric(ci_lo <= true_atet & true_atet <= ci_hi), true_att = true_atet)
      out <- c(out, list(pair(rc, ps, "RC - production SIMEX-lasso")))
      sf <- confirm_fold("dgp_simex", k$cell, "glmnet")
      if (!is.null(sf)) out <- c(out, list(pair(rc, sf, "RC - SIMEX fold (M = 0)")))
    }
    bind_rows(out)
  }) |> ungroup()
  write.csv(long, file.path(FU, "rc/rc_metrics_long.csv"), row.names = FALSE)
  write.csv(paired, file.path(FU, "rc/rc_paired.csv"), row.names = FALSE)
  txt <- capture.output({
    cat("Regression-calibration comparator,", format(Sys.time()), "-", length(fs), "task files,", nrow(errs), "errored\n\n")
    w <- long |> filter(metric %in% c("bias_att", "emp_sd", "se_ratio", "rmse", "cover_att", "median_clip_frac", "n_flag")) |>
      select(study, cell, learner, method, metric, value) |> pivot_wider(names_from = metric, values_from = value)
    print(as.data.frame(w), digits = 3); cat("\nPaired:\n"); print(as.data.frame(paired), digits = 3)
  })
  writeLines(txt, file.path(FU, "rc/rc_report.txt")); cat(txt, sep = "\n")
}
