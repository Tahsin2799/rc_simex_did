# Helpers shared by paper_all_in_one.Rmd (sourced from the abstract and setup).
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(knitr)
  library(kableExtra)
})
ex <- readRDS("exhibits/exhibits.rds")
m  <- ex$prod$metrics
fmt <- function(x, d = 3) formatC(round(x, d) + 0, format = "f", digits = d)   # + 0 drops a negative zero
nm <- function(x, d = 3) paste0("$", fmt(x, d), "$")   # inline: keep the minus with its number
v <- function(cell, est, col = "bias", d = 3) nm(m[m$cell_id == cell & m$estimator == est, col], d)
big <- function(x) format(x, big.mark = ",", scientific = FALSE)
ae <- ex$app$est
app_est <- function(meth, smp = "analytic cohort", d = 2) {
  r <- ae[ae$method == meth & ae$sample == smp, ]
  sprintf("%s (SE %s)", nm(r$estimate, d), fmt(r$se, d))
}
app_simex <- function(r, smp = "analytic cohort") app_est(sprintf("SIMEX-lasso, reliability %.1f (M = 0)", r), smp)
cf <- ex$confirm
# Confirmatory supplemental study: one value (cv), a min--max range over cells (crng),
# paired contrasts (cp) and oracle rows (co). Default arm is the revised (fold) one.
cv <- function(st, ce, le, me, arm = "fold", col = "value") {
  L <- cf$long; r <- L[L$study == st & L$cell == ce & L$learner == le & L$arm == arm & L$metric == me, col]
  if (length(r) != 1) stop("confirm value not unique: ", paste(st, ce, le, arm, me)); r
}
crng <- function(me, st, arm = "fold", d = 3, scale = 1, col = "value", signed = FALSE) {
  x <- scale * cf$long[[col]][cf$long$study %in% st & cf$long$arm %in% arm & cf$long$metric == me]
  if (signed) paste(nm(min(x), d), "to", nm(max(x), d)) else paste0(fmt(min(x), d), "--", fmt(max(x), d))
}
cp <- function(st, ce, contrast, arm = NA, col = "est_diff") {
  P <- cf$paired; keep <- P$study == st & P$cell == ce & P$contrast == contrast
  if (!is.na(arm)) keep <- keep & P$arm %in% arm
  r <- P[keep, col]; if (length(r) != 1) stop("paired value not unique: ", contrast); r
}
co <- function(ce, le, arm, me = "bias_att", col = "value")
  cf$oracle[cf$oracle$cell == ce & cf$oracle$learner == le & cf$oracle$arm == arm & cf$oracle$metric == me, col]
if (!is.null(cf)) {
  sim_st <- c("gauss_known", "gauss_range", "dgp_simex", "bihs_semi")
  fl <- cf$long[cf$long$arm == "fold" & cf$long$study %in% sim_st, ]
  wide_fl <- function(ms) tidyr::pivot_wider(fl[fl$metric %in% ms, c("study", "cell", "learner", "metric", "value")],
                                             names_from = metric, values_from = value)
  se_w <- wide_fl(c("mean_se", "mean_se_fixed")); q2_w <- wide_fl(c("cover_q2_calaware", "cover_q2_fixedse"))
  ck <- list(
    n_cells = nrow(se_w),
    se_infl_pct = 100 * max(se_w$mean_se / se_w$mean_se_fixed - 1),
    share_pct = 100 * max(fl$value[fl$metric == "calibration_var_share"]),
    q2_cov_diff = max(abs(q2_w$cover_q2_calaware - q2_w$cover_q2_fixedse), na.rm = TRUE),
    fold_fixed_diff = max(abs(cf$paired$est_diff[cf$paired$contrast == "fold - fixed"])),
    fold_fixed_mcse = max(cf$paired$est_diff_mcse[cf$paired$contrast == "fold - fixed"]),
    gk5 = "gauss_r0.60_n5000", gk1 = "gauss_r0.60_n1000", gr5 = "gauss_r0.75_n5000", gr1 = "gauss_r0.75_n1000")
  ck$width_ratio <- cv("gauss_known", ck$gk5, "glm", "width_point_Mtrue") / cv("gauss_known", ck$gk5, "glm", "width_point_M0")
}
# Share of the variance of regression calibration carried by its first-step adjustment (Proposition 2),
# Gaussian cells of the third amendment, lasso learner; used in the abstract.
rc_term_rng <- function(d = 0) {
  L <- ex$amend3$long
  x <- 100 * L$value[L$study == "shape" & grepl("^gauss", L$cell) & L$learner == "glmnet" &
                     L$method == "rc_fold" & L$metric == "calibration_var_share"]
  paste0(fmt(min(x), d), "--", fmt(max(x), d))
}
# Largest share of SIMEX's variance carried by its first-step adjustment, in percent: over the Gaussian cells
# of the third amendment (where regression calibration's share is rc_term_rng()), or over every cell of the
# first, third, fourth and fifth amendments.
simex_share <- function(which = c("gauss", "all")) {
  which <- rlang::arg_match(which)
  pick <- \(L) L$value[L$method == "simex_fold" & L$metric == "calibration_var_share" &
                       (which == "all" | grepl("^gauss", L$cell))]
  cl <- ex$confirm$long
  a1 <- cl$value[cl$arm == "fold" & cl$metric == "calibration_var_share"]
  x <- if (which == "gauss") pick(ex$amend3$long) else c(a1, pick(ex$amend3$long), pick(ex$amend4$long), pick(ex$amend5$long))
  100 * max(x, na.rm = TRUE)
}
# Table cells: a true minus sign (U+2212, mapped to $-$ in the preamble) for every negative number,
# including the bounds of bracketed intervals.
mns <- function(x) gsub("(^|[ \\[(])-(?=\\.?[0-9])", "\\1\u2212", x, perl = TRUE)
mns_df <- function(x) { x[] <- lapply(x, \(col) if (is.character(col)) mns(col) else col); x }
kbl2 <- function(x, cap, ..., font = 8.5) kbl(mns_df(x), booktabs = TRUE, caption = cap, linesep = "", ...) |>
  kable_styling(latex_options = c("hold_position"), font_size = font)
# Long appendix tables. keep_groups() ends every body row with \\* except the last row of a group (the
# row before \addlinespace), so LaTeX breaks pages only between groups and no page opens on rows whose
# group label is on the previous page.
kbl_long <- function(x, cap, font = 6.5, gaps = "", ...)
  kbl(mns_df(x), booktabs = TRUE, longtable = TRUE, escape = FALSE, linesep = gaps, caption = cap, ...) |>
  kable_styling(latex_options = c("repeat_header"), font_size = font)
keep_groups <- function(k) {
  s <- strsplit(as.character(k), "\n", fixed = TRUE)[[1]]
  body <- max(grep("\\endhead", s, fixed = TRUE)) + 1L
  for (i in seq(body, length(s) - 1L))
    if (grepl("\\\\$", s[i]) && !grepl("^\\\\(addlinespace|bottomrule|midrule|end)", s[i + 1L]))
      s[i] <- sub("\\\\$", "\\\\*", s[i])
  k[] <- paste(s, collapse = "\n"); k
}
first_only <- function(x) if_else(duplicated(x), "", x)
group_gaps <- function(x) if_else(c(x[-1] != x[-length(x)], FALSE), "\\addlinespace", "")
# Readable labels for the original study's cells (cell_id = design__magnitude__family__N) and estimators.
# One map, so a change of term is made once.
DESIGN_LAB <- c(baseline = "No measurement error", x_classical = "Classical covariate error",
                x_nc = "Differential covariate error", treat_misclass = "Treatment misclassification",
                y_classical = "Classical outcome error", y_diff = "Differential outcome error",
                y_nc_contraction = "Mean-reverting outcome error", y_nc_covariate = "Covariate-dependent outcome shift")
EST_LAB <- c(did_ols = "Covariate-blind DiD", did_ipw = "IPW-DiD", drdid_imp = "Improved DR-DiD (parametric)",
             did_ps_rf = "IPW-DiD, random-forest propensity", did_ps_xgb = "IPW-DiD, boosting propensity",
             dml_lasso = "Lasso DR-DiD", dml_lasso_int = "Lasso DR-DiD with interactions", dml_rf = "Random-forest DR-DiD",
             dml_xgb = "Boosted DR-DiD", cf_did = "Causal forest", bcf_did = "Bayesian causal forest",
             match_full = "Full matching", simex_lasso = "SIMEX-lasso", simex_xgb = "SIMEX-boosting",
             misclass_lasso = "Known-probability correction, lasso", misclass_lasso_R = "Known-probability correction, lasso ($\\tau_R$)",
             misclass_xgb = "Known-probability correction, boosting", misclass_xgb_R = "Known-probability correction, boosting ($\\tau_R$)",
             bounds_lasso = "Reliability range, lasso", bounds_xgb = "Reliability range, boosting")
cell_parts <- function(id) {
  p <- strsplit(id, "__")
  data.frame(design = vapply(p, `[`, "", 1), mag = vapply(p, `[`, "", 2), fam = vapply(p, `[`, "", 3),
             N = as.numeric(sub("N", "", vapply(p, `[`, "", 4))))
}
# Reliability magnitudes ("r0.60") become "$r=0.6$"; the coercion only sees those, so it never warns.
mag_lab <- function(mag) {
  r <- as.numeric(sub("^r", "", if_else(grepl("^r[0-9]", mag), mag, NA_character_)))
  case_when(mag == "none" ~ "", grepl("^r[0-9]", mag) ~ paste0("$r=", r, "$"), .default = mag)
}
# Row label without the design (for tables grouped by design) or with it.
cell_short <- function(id) { z <- cell_parts(id); m <- mag_lab(z$mag)
  paste0(if_else(m == "", "", paste0(m, ", ")), z$fam, ", $n=", gsub(",", "{,}", big(z$N)), "$") }
cell_full <- function(id) paste0(DESIGN_LAB[cell_parts(id)$design], ", ", cell_short(id))
MAG_ORD <- c("none", "r0.90", "r0.75", "r0.60", "small", "medium", "large")
cell_key <- function(id) { z <- cell_parts(id)
  sprintf("%02d%02d%s%05d", match(z$design, names(DESIGN_LAB)), match(z$mag, MAG_ORD), z$fam, z$N) }
# "a - b" contrast codes to labels joined by a true minus.
con_lab <- function(x, lab) vapply(strsplit(x, " - ", fixed = TRUE),
  \(p) paste(if_else(p %in% names(lab), unname(lab[p]), p), collapse = " − "), "")
# Third amendment: implementation-check summary for Appendix H.
rck_text <- function() {
  a <- ex$amend3; if (is.null(a$checks)) return("pending")
  sprintf("%d of %d pass, SE/SD %s, coverage %s", sum(a$checks$pass %in% TRUE), sum(!is.na(a$checks$pass)),
          fmt(a$mc$se_ratio, 2), fmt(a$mc$cover))
}
