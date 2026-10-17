# tau_R truth for the R-form misclassification corrections (second amendment,
# Section 3). Re-scores existing production estimates; refits nothing.
#   cd simulation; Rscript runs/followup_tau_r.R [--pop_n=1000000]
# Writes amendment_2026-10-08/tau_r/:
#   tau_r_truth.csv       realised ATT and realised tau_R per (cell, rep)
#   tau_r_population.csv  population tau_R and ATT on one large draw per family and rate level
#   tau_r_scores.csv      bias and coverage of misclass_* rows against each target
#   tau_r_report.txt

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(digest)
})
source("R_scripts/1_sim_param.R"); source("R_scripts/2_data_sim.R"); source("R_scripts/3_registry.R")
source("R_scripts/35_population_targets.R"); source("runs/grid.R")
args <- commandArgs(trailingOnly = TRUE)
kv <- \(k, d) { h <- grep(paste0("^--", k, "="), args, value = TRUE); if (length(h)) sub(".*=", "", h[1]) else d }
POP_N <- as.integer(kv("pop_n", "1000000")); OUT <- "amendment_2026-10-08/tau_r"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
seed_for <- function(cell_id, rep) as.integer((digest::digest2int(cell_id) + rep) %% .Machine$integer.max)
z <- qnorm(0.975)

# f = kappa e (1 - e) / {q (1 - q)}, q = alpha0 + kappa e (33_misclass.R header).
f_of <- function(e, a0, a1) { k <- 1 - a0 - a1; q <- a0 + k * e; k * e * (1 - e) / (q * (1 - q)) }

grid <- build_grid() |> filter(me_design == "treat_misclass")
truth <- bind_rows(lapply(seq_len(nrow(grid)), \(i) {
  row <- grid[i, ]; a <- row$args[[1]]
  bind_rows(lapply(1:500, \(rep) {
    seed_r <- seed_for(row$cell_id, rep); set.seed(seed_r)
    gen <- data_gen_v2(N = row$N, params = default_params(row$panel_strength),
                       dgp_family = row$dgp_family, seed = seed_r, return_cf = TRUE)
    A <- wide_frame(gen$observed)$A                    # latent treatment, id order
    f <- f_of(gen$e_X, a$alpha0, a$alpha1)
    tibble(cell_id = row$cell_id, rep, att = mean(gen$theta_i[A == 1]), true_atet = gen$true_atet,
           tau_r = sum(A * gen$theta_i * f) / sum(A * f))
  }))
}))
if (!isTRUE(all.equal(truth$att, truth$true_atet))) {       # matches the stored original truth
  cli::cli_abort("Recomputed ATT differs from the stored truth: largest gap {signif(max(abs(truth$att - truth$true_atet)), 3)}.")
}
write.csv(truth, file.path(OUT, "tau_r_truth.csv"), row.names = FALSE)

# Population: tau_R = E[theta f e] / E[f e] and ATT = E[theta e] / E[e] on one large draw,
# with a delta-method SE for their difference (MC10 predicts zero).
lv <- distinct(grid, mag, args) |> mutate(a0 = vapply(args, `[[`, numeric(1), "alpha0"), a1 = vapply(args, `[[`, numeric(1), "alpha1"))
pop <- bind_rows(lapply(c("linear", "nonlinear"), \(fam) {
  set.seed(20261008); gen <- data_gen_v2(N = POP_N, params = default_params("low"), dgp_family = fam,
                                         seed = 20261008, return_cf = TRUE)
  th <- gen$theta_i; e <- gen$e_X
  bind_rows(lapply(seq_len(nrow(lv)), \(j) {
    fe <- f_of(e, lv$a0[j], lv$a1[j]) * e
    tr <- sum(th * fe) / sum(fe); ta <- sum(th * e) / sum(e)
    psi <- (th - tr) * fe / mean(fe) - (th - ta) * e / mean(e)
    tibble(dgp_family = fam, mag = lv$mag[j], n = POP_N, tau_r = tr, att = ta, diff = tr - ta,
           diff_se = sd(psi) / sqrt(POP_N), population_att = population_att(fam))
  }))
}))
write.csv(pop, file.path(OUT, "tau_r_population.csv"), row.names = FALSE)

# Re-score: R-form rows against tau_R and against the ATT; H-form rows against the ATT for reference.
wilson <- function(x) { n <- length(x); p <- mean(x); den <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / den) }
pr <- readRDS("sim_results_v2/_per_rep.rds") |>
  filter(cell_id %in% grid$cell_id, grepl("^misclass_", estimator), status == "ok") |>
  inner_join(select(truth, cell_id, rep, tau_r), by = c("cell_id", "rep"))
score <- function(g, tgt) { cov <- g$ci_lo <= tgt & tgt <= g$ci_hi; w <- wilson(cov)
  tibble(n = nrow(g), bias = mean(g$est - tgt), bias_mcse = sd(g$est - tgt) / sqrt(nrow(g)),
         rmse = sqrt(mean((g$est - tgt)^2)), coverage = mean(cov), cov_wilson_lo = w[1], cov_wilson_hi = w[2]) }
scores <- pr |> group_by(cell_id, estimator) |>
  group_modify(\(g, k) bind_rows(mutate(score(g, g$true_atet), target = "realised ATT"),
                                 mutate(score(g, g$tau_r), target = "realised tau_R"))) |>
  ungroup()
gap <- truth |> group_by(cell_id) |>
  summarise(gap_mean = mean(tau_r - att), gap_sd = sd(tau_r - att), gap_maxabs = max(abs(tau_r - att)), .groups = "drop")
scores <- left_join(scores, gap, by = "cell_id")
write.csv(scores, file.path(OUT, "tau_r_scores.csv"), row.names = FALSE)

txt <- capture.output({
  cat("tau_R truth and re-scoring,", format(Sys.time()), "\n\nPopulation (one draw of", POP_N, "per family):\n")
  print(as.data.frame(pop), digits = 4)
  cat("\nPer-replicate gap, realised tau_R minus realised ATT:\n"); print(as.data.frame(gap), digits = 3)
  cat("\nScores (R-form rows):\n")
  print(as.data.frame(filter(scores, grepl("_R$", estimator)) |>
    select(cell_id, estimator, target, bias, bias_mcse, coverage, cov_wilson_lo, cov_wilson_hi)), digits = 3)
})
writeLines(txt, file.path(OUT, "tau_r_report.txt")); cat(txt, sep = "\n")
