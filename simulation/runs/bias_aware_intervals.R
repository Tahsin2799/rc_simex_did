# POST HOC (14 October 2026; added after all results were seen, at the author's request).
# Bias-aware intervals with the Armstrong-Kolesar critical value (Armstrong & Kolesar 2020):
#   est +/- cv_alpha(M kappa(t) / se) * se,  cv_alpha(b) = sqrt of the 1 - alpha quantile of a noncentral
#   chi-square(1) with noncentrality b^2,
# against the pre-declared additive interval est +/- (z se + M kappa(t)). Both are valid when |q'''| <= M;
# the Armstrong-Kolesar interval is never longer. Everything is computed from saved estimates and SEs.
#   cd simulation; Rscript runs/bias_aware_intervals.R
# Writes bias_aware/ak_confirm.csv, ak_amend.csv, ak_application.csv and ak_report.txt.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})
source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
source("R_scripts/34_simex_fold.R")
OUT <- "bias_aware"; dir.create(OUT, showWarnings = FALSE)
z <- qnorm(0.975); LAM <- c(0, 0.5, 1, 2)
# For large b the quantile is b + z_{1-alpha} up to exponentially small terms (the lower tail is negligible),
# and qchisq() stops converging there, so the limit is used for b > 8.
cv_ak <- function(b, alpha = 0.05) if_else(b > 8, b + qnorm(1 - alpha), sqrt(qchisq(1 - alpha, df = 1, ncp = pmin(b, 8)^2)))
rd <- \(dir) bind_rows(lapply(list.files(dir, "[.]rds$", recursive = TRUE, full.names = TRUE), \(p) readRDS(p)$rows))

# A. First supplemental (confirmatory) study: pointwise ATT intervals at t*, M in the pre-declared grid.
cf <- rd("amendment_2026-10-06/confirm/raw") |> filter(study != "oracle")   # oracle rows carry no SIMEX interval
kap <- tibble(t_point = unique(cf$t_point)) |> mutate(kappa = vapply(t_point, \(t) peano_kappa(LAM, t), 0))
cf <- cf |> left_join(kap, by = "t_point")
Ms <- c(0, 0.005, 0.01, 0.02, 0.05)
long <- bind_rows(lapply(Ms, \(M) {
  lo <- cf[[paste0("pt_lo_M", M)]]; hi <- cf[[paste0("pt_hi_M", M)]]
  h_ak <- cv_ak(M * cf$kappa / cf$se) * cf$se
  tibble(study = cf$study, cell = cf$cell, learner = cf$learner, arm = cf$arm, M = M, est = cf$est, se = cf$se,
         kappa = cf$kappa, att = cf$true_att, lo_add = lo, hi_add = hi,
         repro = pmax(abs(lo - (cf$est - z * cf$se - M * cf$kappa)), abs(hi - (cf$est + z * cf$se + M * cf$kappa))),
         lo_ak = cf$est - h_ak, hi_ak = cf$est + h_ak)
}))
if (max(long$repro, na.rm = TRUE) >= 1e-8) {                # the saved additive intervals are reproduced
  cli::cli_abort("First amendment: saved intervals not reproduced, largest gap {signif(max(long$repro, na.rm = TRUE), 3)}.")
}
akA <- long |> group_by(study, cell, learner, arm, M) |> summarise(
  n = n(), cover_add = mean(lo_add <= att & att <= hi_add), cover_ak = mean(lo_ak <= att & att <= hi_ak),
  width_add = mean(hi_add - lo_add), width_ak = mean(hi_ak - lo_ak), .groups = "drop") |>
  mutate(width_ratio = width_ak / width_add)
write.csv(akA, file.path(OUT, "ak_confirm.csv"), row.names = FALSE)
# A2. Where the bound is known to hold: closed-form designs at their true M (|q'''| <= M_true by construction).
akT <- cf |> filter(!is.na(M_true)) |>
  mutate(h_add = z * se + M_true * kappa, h_ak = cv_ak(M_true * kappa / se) * se) |>
  group_by(study, cell, learner, arm) |> summarise(
    n = n(), M_true = first(M_true), cover_add = mean(abs(est - true_att) <= h_add),
    cover_ak = mean(abs(est - true_att) <= h_ak), width_ratio = mean(h_ak) / mean(h_add), .groups = "drop")
write.csv(akT, file.path(OUT, "ak_confirm_trueM.csv"), row.names = FALSE)

# B. Third to fifth amendments: SIMEX intervals at M = 0.05, t = -1.
k1 <- peano_kappa(LAM, -1)
am <- bind_rows(lapply(c("amendment_2026-10-11/study/raw", "amendment_2026-10-11b/study/raw", "amendment_2026-10-14/study/raw"),
                       \(d) rd(d) |> filter(method == "simex_fold") |> mutate(src = d)))
am <- am |> mutate(h_ak = cv_ak(0.05 * k1 / se) * se,
                   repro = pmax(abs(ci_lo_M005 - (est - z * se - 0.05 * k1)), abs(ci_hi_M005 - (est + z * se + 0.05 * k1))))
if (max(am$repro) >= 1e-8) cli::cli_abort("Later amendments: saved intervals not reproduced, largest gap {signif(max(am$repro), 3)}.")
akB <- am |> group_by(src, study, cell, learner) |> summarise(
  n = n(), cover_add = mean(ci_lo_M005 <= true_att & true_att <= ci_hi_M005),
  cover_ak = mean(est - h_ak <= true_att & true_att <= est + h_ak),
  width_add = mean(ci_hi_M005 - ci_lo_M005), width_ak = mean(2 * h_ak), .groups = "drop") |>
  mutate(width_ratio = width_ak / width_add, src = basename(dirname(dirname(src))))
write.csv(akB, file.path(OUT, "ak_amend.csv"), row.names = FALSE)

# C. Application: Armstrong-Kolesar intervals over (reliability, M) and the breakdown value of M,
#    the largest M at which the interval still excludes zero (Rambachan & Roth 2023 report the same quantity).
ag <- read.csv("../application/output/reliability_M_grid.csv")
ag <- ag |> mutate(h_ak = cv_ak(M * kappa / se) * se, lo_ak = estimate - h_ak, hi_ak = estimate + h_ak,
                   repro = pmax(abs(ci_lo - (estimate - z * se - M * kappa)), abs(ci_hi - (estimate + z * se + M * kappa))))
if (max(ag$repro) >= 1e-8) cli::cli_abort("Application: saved intervals not reproduced, largest gap {signif(max(ag$repro), 3)}.")
breakdown <- function(est, se, kappa, ak = TRUE) {
  f <- \(M) abs(est) - (if (ak) cv_ak(M * kappa / se) * se else z * se + M * kappa)
  if (f(0) <= 0) return(0)                                   # zero is inside the interval already at M = 0
  uniroot(f, c(0, 100), tol = 1e-10)$root
}
bd <- ag |> distinct(sample, reliability, estimate, se, kappa) |> rowwise() |>
  mutate(breakdown_M_ak = breakdown(estimate, se, kappa), breakdown_M_add = breakdown(estimate, se, kappa, FALSE)) |>
  ungroup()
write.csv(ag |> select(-repro), file.path(OUT, "ak_application.csv"), row.names = FALSE)
write.csv(bd, file.path(OUT, "ak_breakdown.csv"), row.names = FALSE)

txt <- capture.output({
  options(width = 200)
  cat("POST HOC bias-aware intervals (Armstrong-Kolesar critical value),", format(Sys.time()), "\n")
  cat("Saved additive intervals reproduced to < 1e-8 in all three sources.\n\nA. Confirmatory study, pointwise ATT intervals:\n")
  print(as.data.frame(akA |> filter(M > 0) |> group_by(M) |> summarise(
    cells = n(), mean_width_ratio = mean(width_ratio), min_width_ratio = min(width_ratio),
    max_cover_drop = max(cover_add - cover_ak), mean_cover_add = mean(cover_add), mean_cover_ak = mean(cover_ak))), digits = 4)
  cat("\nA2. Closed-form designs at the true M (bound holds):\n"); print(as.data.frame(akT), digits = 4)
  cat("\nB. Amendments, SIMEX at M = 0.05:\n"); print(as.data.frame(akB), digits = 3)
  cat("\nC. Application, breakdown value of M (interval excludes zero for M below it):\n"); print(as.data.frame(bd), digits = 4)
})
writeLines(txt, file.path(OUT, "ak_report.txt")); cat(txt, sep = "\n")
