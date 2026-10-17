# Build every table, figure and inline number the manuscript uses.
#   Rscript paper1_rc_simex_did/submission/R/build_exhibits.R    (from the repository root; manuscript/R/... in the code release)
# Reads saved results only; never fits an estimator. Writes exhibits/exhibits.rds and figures/*.pdf in the
# manuscript folder, which is the parent of this script's R/ folder.
# Sections whose inputs are missing (e.g. the confirmatory study before it
# finishes) are stored as NULL and the manuscript says so.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
OUT <- if (length(script)) dirname(dirname(script)) else "paper1_rc_simex_did/submission"
rd  <- \(p) if (file.exists(p)) read.csv(p, stringsAsFactors = FALSE) else NULL
dir.create(file.path(OUT, "exhibits"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT, "figures"), showWarnings = FALSE, recursive = TRUE)
ex <- list(built = format(Sys.time(), "%Y-%m-%d %H:%M"))

# Reference categorical palette (dataviz skill, slots 1-4 in fixed order),
# always paired with shape and linetype so identity survives grayscale print.
PAL <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100")
theme_paper <- function() theme_minimal(base_size = 9) + theme(
  panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
  strip.text = element_text(face = "plain", size = 8.5), legend.position = "bottom",
  legend.title = element_blank(), axis.title = element_text(size = 8.5),
  plot.margin = margin(4, 6, 4, 4))
save_pdf <- \(p, f, w = 6.5, h = 3.6) ggsave(file.path(OUT, "figures", f), p, width = w, height = h,
                                             device = cairo_pdf)

# ---- original production study ----------------------------------------------
m  <- rd("simulation/tables_v2/metrics_all_cells.csv")
pr <- rd("simulation/tables_v2/paired_comparisons.csv")
rg <- rd("simulation/tables_v2/range_estimators.csv")
fl <- rd("simulation/tables_v2/flags_failures.csv")
tg <- rd("simulation/tables_v2/population_targets.csv")
ex$prod <- list(metrics = m, paired = pr, range = rg, flags = fl, targets = tg,
                headline = rd("simulation/tables_v2/headline_N5000.csv"))
ex$prod$status <- c(rows = sum(m$n_sim), ok = sum(m$n_ok), flagged = sum(m$n_flag), failed = sum(m$n_fail))
pp <- "simulation/sim_results_v2/_per_rep.rds"
# The per-replicate file is not part of the code release; its total (computed 2026-10-17) stands in for it.
ex$prod$elapsed_days <- if (file.exists(pp)) sum(readRDS(pp)$runtime_s) / 86400 else 113.6652701198
ex$prod$pop_att <- c(linear = 0.2, nonlinear = unique(tg$population_att[grepl("nonlinear", tg$cell_id)]))
ex$prod$realised_range <- range(c(tg$realised_min[grepl("nonlinear", tg$cell_id)],
                                  tg$realised_max[grepl("nonlinear", tg$cell_id)]))

# Figure: bias under classical covariate error, naive versus SIMEX.
fig_est <- c(dml_lasso = "Lasso DR-DiD", simex_lasso = "SIMEX-lasso (original)",
             dml_xgb = "Boosted DR-DiD", simex_xgb = "SIMEX-boosting (original)")
bias_df <- m |>
  filter(estimator %in% names(fig_est),
         me_design == "x_classical" | (me_design == "baseline" & estimator %in% c("dml_lasso", "dml_xgb"))) |>
  mutate(reliability = if_else(me_design == "baseline", 1,
                               as.numeric(sub("^r", "", if_else(grepl("^r[0-9]", mag), mag, NA_character_)))),
         est = factor(fig_est[estimator], levels = fig_est),
         panel = sprintf("%s trends, n = %s", dgp_family, format(N, big.mark = ",")))
p_bias <- ggplot(bias_df, aes(reliability, bias, colour = est, shape = est, linetype = est)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
  geom_line(linewidth = 0.5, position = position_dodge(0.025)) +
  geom_pointrange(aes(ymin = bias - 2 * bias_mcse, ymax = bias + 2 * bias_mcse),
                  size = 0.25, linewidth = 0.35, position = position_dodge(0.025)) +
  facet_wrap(~panel, nrow = 1) + scale_x_reverse(breaks = c(1, 0.9, 0.75, 0.6)) +
  scale_colour_manual(values = PAL) + scale_shape_manual(values = c(16, 17, 15, 18)) +
  scale_linetype_manual(values = c("solid", "22", "solid", "22")) +
  labs(x = "Covariate reliability (1 = error-free baseline cell)", y = "Bias (sample ATT), ±2 MCSE") +
  theme_paper()
save_pdf(p_bias, "fig_sim_bias.pdf", h = 3.2)

# ---- implementation checks and pilot ----------------------------------------
chk <- rd("simulation/amendment_2026-10-06/checks/checks.csv")
ex$checks <- list(n = nrow(chk), pass = sum(chk$pass),
                  mc = rd("simulation/amendment_2026-10-06/checks/mc_closed_form_coverage.csv"),
                  se = rd("simulation/amendment_2026-10-06/checks/mc_closed_form_se.csv"),
                  large_grid = rd("simulation/amendment_2026-10-06/checks/large_n_grid.csv"),
                  large_D = rd("simulation/amendment_2026-10-06/checks/large_n_gradient.csv"))
ex$pilot <- list(summary = rd("simulation/amendment_2026-10-06/pilot/pilot_summary.csv"),
                 timing = rd("simulation/amendment_2026-10-06/pilot/pilot_timing.csv"))

# ---- confirmatory supplemental study (may be pending) -----------------------
cl <- rd("simulation/amendment_2026-10-06/confirm/confirm_metrics_long.csv")
ex$confirm <- if (is.null(cl)) NULL else list(
  long = cl, paired = rd("simulation/amendment_2026-10-06/confirm/confirm_paired.csv"),
  oracle = rd("simulation/amendment_2026-10-06/confirm/confirm_oracle.csv"),
  compute = rd("simulation/amendment_2026-10-06/confirm/confirm_compute.csv"))

# ---- second amendment: tau_R truth, regression calibration, BCF pilot ----------
fu <- "simulation/amendment_2026-10-08"
ex$followup <- list(
  tau_r = list(scores = rd(file.path(fu, "tau_r/tau_r_scores.csv")),
               population = rd(file.path(fu, "tau_r/tau_r_population.csv"))),
  rc = list(long = rd(file.path(fu, "rc/rc_metrics_long.csv")), paired = rd(file.path(fu, "rc/rc_paired.csv"))),
  bcf = list(fits = rd(file.path(fu, "bcf/bcf_fits.csv")), summary = rd(file.path(fu, "bcf/bcf_summary.csv"))),
  checks = list(mcmc = rd(file.path(fu, "checks/mcmc_diag_checks.csv")), rc = rd(file.path(fu, "checks/rc_checks.csv"))))
ex$followup$bcf$paired <- rd(file.path(fu, "bcf/bcf_paired.csv"))

# ---- third amendment: covariate shape and fold-calibrated regression calibration ----
a3 <- "simulation/amendment_2026-10-11"
ex$amend3 <- list(checks = rd(file.path(a3, "checks/rc_fold_checks.csv")), mc = rd(file.path(a3, "checks/rc_fold_mc.csv")),
                  long = rd(file.path(a3, "study/amend3_metrics_long.csv")), paired = rd(file.path(a3, "study/amend3_paired.csv")),
                  compute = rd(file.path(a3, "study/amend3_compute.csv")),
                  var_components = rd(file.path(a3, "checks/rc_var_components.csv")),
                  var_tails = rd(file.path(a3, "checks/rc_var_tails.csv")))
a4 <- "simulation/amendment_2026-10-11b"
ex$amend4 <- list(long = rd(file.path(a4, "study/amend3_metrics_long.csv")), paired = rd(file.path(a4, "study/amend3_paired.csv")),
                  compute = rd(file.path(a4, "study/amend3_compute.csv")), map = rd(file.path(a4, "population_map.csv")),
                  map_lnorm = rd(file.path(a4, "population_map_lnorm.csv")), M = rd(file.path(a4, "population_M.csv")))
ex$amend4$tail <- rd(file.path(a4, "checks/tail_robust.csv"))
ex$amend4$rep389 <- rd(file.path(a4, "checks/tail_rep389.csv"))
ex$amend4$extrap <- rd(file.path(a4, "checks/tail_extrapolation.csv"))
# ---- fifth amendment (clamped spline, sensitivity) ----
a5 <- "simulation/amendment_2026-10-14"
ex$amend5 <- list(long = rd(file.path(a5, "study/amend3_metrics_long.csv")), compare = rd(file.path(a5, "study/amend5_compare.csv")),
                  compute = rd(file.path(a5, "study/amend3_compute.csv")))
# ---- population map (Section 5), calibration-gap checks, bias-aware intervals (post hoc) ----
pm <- "simulation/population_map"
ex$popmap <- list(map = rd(file.path(pm, "population_map.csv")), gapchk = rd(file.path(pm, "gap_diagnostic_check.csv")),
                  biascheck = rd(file.path(pm, "gap_bias_check.csv")), original = rd(file.path(pm, "gap_original_design.csv")))
ba <- "simulation/bias_aware"
ex$ak <- list(confirm = rd(file.path(ba, "ak_confirm.csv")), trueM = rd(file.path(ba, "ak_confirm_trueM.csv")),
              amend = rd(file.path(ba, "ak_amend.csv")), app = rd(file.path(ba, "ak_application.csv")),
              breakdown = rd(file.path(ba, "ak_breakdown.csv")))
ex$popmap$convergence <- rd(file.path(pm, "population_map_convergence.csv"))
ex$popmap$mc <- rd(file.path(pm, "population_map_mc_check.csv"))
ex$popmap$resolution <- rd(file.path(pm, "gap_resolution.csv"))

# Figure: population biases of regression calibration and SIMEX (nonparametric nuisances, no sampling),
# equal group dispersion. Columns are the skewness of the latent covariates, rows the trend.
if (!is.null(ex$popmap$map)) {
  sk_lab <- \(sl, sk) if_else(sl == 0, "Gaussian", sprintf("Skewness %.2f", sk))
  pmap <- ex$popmap$map |> filter(ratio == 1) |>
    select(sl, skewness, r, bq, rc, simex) |>
    pivot_longer(c(rc, simex), names_to = "method", values_to = "bias") |>
    mutate(method = factor(recode(method, rc = "Regression calibration", simex = "SIMEX"),
                           levels = c("Regression calibration", "SIMEX")),
           shape_lab = factor(sk_lab(sl, skewness), levels = unique(sk_lab(sort(unique(sl)), sort(unique(skewness))))),
           trend = factor(if_else(bq == 0, "Linear trend", "Quadratic trend"), levels = c("Linear trend", "Quadratic trend")))
  p_map <- ggplot(pmap, aes(r, bias, colour = method, shape = method, linetype = method)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_line(linewidth = 0.5) + geom_point(size = 1.3) +
    facet_grid(trend ~ shape_lab) +
    scale_colour_manual(values = PAL[c(1, 2)]) + scale_shape_manual(values = c(16, 17)) +
    scale_linetype_manual(values = c("solid", "22")) +
    scale_x_continuous(breaks = c(0.5, 0.7, 0.9), labels = c("0.5", "0.7", "0.9")) +
    labs(x = "Reliability", y = "Population bias (target minus ATT)") +
    theme_paper() + theme(panel.border = element_rect(fill = NA, colour = "grey80", linewidth = 0.3),
                          panel.spacing.x = unit(6, "pt"))
  save_pdf(p_map, "fig_pop_map.pdf", w = 6.5, h = 3.6)
}

# Figure: coverage of fold-calibrated regression calibration with the fixed-calibration and the
# calibration-aware SE, Wilson 95% intervals; colour plus shape so identity survives grayscale.
rc_cov <- bind_rows(ex$amend3$long |> mutate(src = "a3"), if (!is.null(ex$amend4$long)) ex$amend4$long |> mutate(src = "a4")) |>
  filter(method == "rc_fold", metric %in% c("cover_att", "cover_att_fixedse")) |>
  mutate(group = case_when(src == "a3" & study == "rc_dgp" ~ "Original design",
                           grepl("^gauss", cell) & src == "a3" ~ "Gaussian covariates, r = 0.6",
                           grepl("^lnorm", cell) & src == "a3" ~ "Log-normal covariates, r = 0.6",
                           TRUE ~ "Reliability 0.9 (fourth amendment)"),
         n = as.numeric(sub(".*_[nN]", "", cell)),
         trend = case_when(grepl("_quad_", cell) ~ "quadratic", grepl("nonlinear", cell) ~ "nonlinear",
                           TRUE ~ "linear"),
         shape = case_when(grepl("^lnorm", cell) ~ "log-normal, ", grepl("^gauss", cell) & src == "a4" ~ "Gaussian, ",
                           .default = ""),
         rel = case_when(grepl("r0.75", cell) ~ ", r = 0.75", grepl("r0.60", cell) ~ ", r = 0.6", .default = ""),
         lab = sprintf("%s%s%s, n = %s, %s", shape, trend, rel, format(n, big.mark = ","),
                       recode(learner, glm = "parametric", glmnet = "lasso", glmnet_ns = "spline lasso")),
         se = factor(if_else(metric == "cover_att", "SE with the first-step adjustment (Proposition 2)", "SE treating the calibration as fixed"),
                     levels = c("SE treating the calibration as fixed", "SE with the first-step adjustment (Proposition 2)")),
         group = factor(group, levels = c("Gaussian covariates, r = 0.6", "Log-normal covariates, r = 0.6",
                                          "Original design", "Reliability 0.9 (fourth amendment)")))
p_rc <- ggplot(rc_cov, aes(value, reorder(lab, desc(lab)), colour = se, shape = se)) +
  geom_vline(xintercept = 0.95, colour = "grey55", linewidth = 0.3, linetype = "22") +
  geom_pointrange(aes(xmin = wilson_lo, xmax = wilson_hi), size = 0.28, linewidth = 0.45,
                  position = position_dodge(width = 0.55)) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = PAL[c(2, 1)]) + scale_shape_manual(values = c(17, 16)) +
  scale_x_continuous(breaks = seq(0.3, 1, 0.1)) +
  labs(x = "Coverage of the ATT by nominal 95% intervals, with Wilson 95% intervals", y = NULL) +
  theme_paper() + theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5), axis.text.y = element_text(size = 7))
save_pdf(p_rc, "fig_rc_coverage.pdf", w = 6.5, h = 4.6)
trc <- rd(file.path(fu, "bcf/bcf_trace_rep1.csv"))
if (!is.null(trc)) {
  # Trace plots, replicate 1 of each pilot cell: four chains, both burn-in settings.
  trc <- trc |> mutate(cell = factor(recode(cell, baseline__none__linear__N5000 = "Baseline,\nlinear",
                                            y_diff__medium__linear__N5000 = "Diff. error,\nlinear",
                                            y_diff__medium__nonlinear__N5000 = "Diff. error,\nnonlinear"),
                                     levels = c("Baseline,\nlinear", "Diff. error,\nlinear", "Diff. error,\nnonlinear")),
                       config = factor(recode(config, prod = "250 burn-in (original)", long = "1,000 burn-in"),
                                       levels = c("250 burn-in (original)", "1,000 burn-in")),
                       chain = factor(chain))
  p_tr <- ggplot(trc, aes(iter, draw, colour = chain)) + geom_line(linewidth = 0.2, alpha = 0.85) +
    facet_grid(cell ~ config, scales = "free_y") + scale_colour_manual(values = PAL, name = "Chain") +
    labs(x = "Retained draw", y = "Treated-average effect") + theme_paper() +
    theme(legend.title = element_text(size = 8), strip.text.y = element_text(size = 7.5))
  save_pdf(p_tr, "fig_bcf_trace.pdf", w = 6.5, h = 5)
}
if (!is.null(cl)) {
  # Figure: coverage of each target by its own interval, Wilson 95% intervals.
  # Rows are studies, columns are targets, y is the design; x is shared within a
  # column (curve targets near 0.95 on a narrow scale, ATT targets on [0, 1]).
  # The SE-without-calibration-term and outer-true-M columns are left to the text
  # and table (they differ from their neighbours by at most 0.004 and are all 1).
  # Facet labels in plotmath, so q_2 and t* are typeset (the plot font has no subscript-two glyph).
  key <- c(cover_q2_calaware = 'atop(italic(q)[2](italic(t)^"*"), "pointwise")',
           cover_band_q2 = 'atop(italic(q)[2]~"on [−1, 0]", "band")',
           cover_att_point_M0 = 'atop("ATT,"~italic(M)==0, "interval (5)")',
           cover_att_point_Mtrue = 'atop("ATT, true"~italic(M), "interval (5)")',
           cover_att_outer_M0 = 'atop("ATT,"~italic(M)==0, "interval (6)")')
  st_key <- c(gauss_known = "Closed form,\nreliability\nknown", gauss_range = "Closed form,\nreliability\nin [0.6, 1]",
              dgp_simex = "Original\ndesign", bihs_semi = "BIHS semi-\nsynthetic")
  # case_when evaluates every branch, so cnt() must accept labels that are not counts without warning.
  cnt <- \(x) formatC(as.integer(if_else(grepl("^[0-9]+$", x), x, NA_character_)), big.mark = ",", format = "d")
  dlab <- \(st, ce, le) case_when(
    startsWith(ce, "gauss") ~ paste0("n = ", cnt(sub(".*_n", "", ce)), if_else(le == "glm", ", param.", ", lasso")),
    st == "dgp_simex"       ~ sprintf("r = %s, n = %s", sub(".*__r([0-9.]+)__.*", "\\1", ce), cnt(sub(".*__N", "", ce))),
    .default               = "254 villages")
  cd <- cl |> filter(metric %in% names(key), arm %in% c("fixed", "fold"), !is.na(value)) |>   # true M: closed form only
    mutate(target = factor(key[metric], levels = key), study = factor(st_key[study], levels = st_key),
           design = dlab(as.character(names(st_key)[as.integer(study)]), cell, learner))
  cd$design <- factor(cd$design, levels = rev(unique(cd$design[order(cd$study, cd$cell, cd$learner)])))
  # Column-wise x limits through invisible points (facet_grid shares x within a column).
  lims <- expand.grid(target = factor(key, levels = key), study = factor(st_key[1], levels = st_key)) |>
    mutate(lo = if_else(grepl("italic(q)", as.character(target), fixed = TRUE), 0.85, 0), hi = 1) |>
    pivot_longer(c(lo, hi), values_to = "value")
  if (nrow(cd)) {
    p_cov <- ggplot(cd, aes(value, design, colour = arm, shape = arm)) +
      geom_blank(data = lims, inherit.aes = FALSE, aes(x = value)) +
      geom_vline(xintercept = 0.95, colour = "grey55", linewidth = 0.3, linetype = "22") +
      geom_linerange(aes(xmin = wilson_lo, xmax = wilson_hi), linewidth = 0.45,
                     position = position_dodge(0.6)) +
      geom_point(size = 1.5, position = position_dodge(0.6)) +
      facet_grid(study ~ target, scales = "free", space = "free_y", labeller = labeller(target = label_parsed)) +
      scale_colour_manual(values = PAL[c(1, 2)], labels = c(fixed = "known covariance", fold = "cross-fitted calibration")) +
      scale_shape_manual(values = c(16, 17), labels = c(fixed = "known covariance", fold = "cross-fitted calibration")) +
      scale_x_continuous(breaks = \(l) if (l[1] > 0.5) c(0.9, 0.95, 1) else c(0, 0.5, 1),
                         labels = \(x) sub("\\.?0+$", "", formatC(x, format = "f", digits = 2))) +
      labs(x = "Coverage, with Wilson 95% intervals", y = NULL) + theme_paper() +
      theme(strip.text.y = element_text(angle = 0, hjust = 0, size = 7.5), strip.text.x = element_text(size = 7.5),
            axis.text = element_text(size = 7), panel.spacing.x = unit(8, "pt"),
            panel.border = element_rect(fill = NA, colour = "grey80", linewidth = 0.3))
    save_pdf(p_cov, "fig_confirm_coverage.pdf", w = 6.5, h = 4.6)
  }
}

# ---- BIHS application -------------------------------------------------------
ao <- "application/output"
ex$app <- list(flow = rd(file.path(ao, "sample_flow.csv")),
               history = rd(file.path(ao, "exposure_history.csv")),
               receipt = rd(file.path(ao, "receipt_by_attendance_2015.csv")),
               est = rd(file.path(ao, "estimates.csv")), grid = rd(file.path(ao, "reliability_M_grid.csv")),
               range = rd(file.path(ao, "range_intervals.csv")), curve = rd(file.path(ao, "simex_curve.csv")),
               balance = rd(file.path(ao, "balance.csv")), overlap = rd(file.path(ao, "overlap.csv")),
               pt = rd(file.path(ao, "pt_breakdown.csv")), semipop = rd(file.path(ao, "semisynth_population.csv")),
               rc = rd(file.path(ao, "rc_estimates.csv")), gap = rd(file.path(ao, "gap_share.csv")))
gd <- ex$app$grid
if (!is.null(gd)) {
  # SIMEX at three values of M; regression calibration (no M) repeated in every panel, with the lasso
  # DR-DiD estimate at reliability 1, where the calibration is the identity.
  sim <- gd |> filter(sample == "analytic cohort", M %in% c(0, 0.05, 0.2)) |>
    transmute(reliability, M, estimate, ci_lo, ci_hi, method = "SIMEX")
  if (!is.null(ex$app$rc)) {
    l1 <- ex$app$est |> filter(sample == "analytic cohort", method == "Lasso DR-DiD (cross-fit)")
    rcd <- bind_rows(ex$app$rc |> filter(method == "Regression calibration") |> select(reliability, estimate, ci_lo, ci_hi),
                     data.frame(reliability = 1, estimate = l1$estimate, ci_lo = l1$ci_lo, ci_hi = l1$ci_hi)) |>
      mutate(method = "Regression calibration") |> crossing(M = c(0, 0.05, 0.2))
    sim <- bind_rows(rcd, sim)
  }
  sim <- sim |> mutate(Mlab = factor(sprintf("M = %s", M), levels = sprintf("M = %s", c(0, 0.05, 0.2))),
                       method = factor(method, levels = c("Regression calibration", "SIMEX")))
  p_app <- ggplot(sim, aes(reliability, estimate, colour = method, fill = method, shape = method, linetype = method)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.3) +
    geom_linerange(aes(ymin = ci_lo, ymax = ci_hi), linewidth = 0.45, linetype = "solid",
                   position = position_dodge(width = 0.045)) +
    geom_line(linewidth = 0.5, position = position_dodge(width = 0.045)) +
    geom_point(size = 1.5, position = position_dodge(width = 0.045)) +
    facet_wrap(~Mlab, nrow = 1) + scale_x_reverse(breaks = c(1, 0.9, 0.8, 0.7, 0.6)) +
    scale_colour_manual(values = PAL[c(1, 2)]) + scale_fill_manual(values = PAL[c(1, 2)]) +
    scale_shape_manual(values = c(16, 17)) + scale_linetype_manual(values = c("solid", "22")) +
    labs(x = "Assumed reliability of log assets and log land",
         y = "Adjusted contrast (school years), 95% interval") + theme_paper()
  save_pdf(p_app, "fig_app_sensitivity.pdf", h = 3.1)
}

# Figure 1 (causal structures) is drawn in TikZ inside paper_all_in_one.Rmd.

saveRDS(ex, file.path(OUT, "exhibits", "exhibits.rds"))
cat("exhibits written:", paste(names(ex), collapse = ", "), "| confirm available:", !is.null(ex$confirm), "\n")
