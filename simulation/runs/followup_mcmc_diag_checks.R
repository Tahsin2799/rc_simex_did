# Checks of R_scripts/37_mcmc_diag.R against chains whose answers are known
# (second amendment, Section 1). Writes amendment_2026-10-08/checks/mcmc_diag_checks.csv.
#   cd simulation; Rscript runs/followup_mcmc_diag_checks.R
source("R_scripts/37_mcmc_diag.R")
dir.create("amendment_2026-10-08/checks", recursive = TRUE, showWarnings = FALSE)
set.seed(20261008)
ar1 <- function(n, rho) { x <- numeric(n); x[1] <- rnorm(1) / sqrt(1 - rho^2)
  for (t in 2:n) x[t] <- rho * x[t - 1] + rnorm(1); x }
m <- 4; n <- 1000; S <- m * n; R <- 200
one <- function(rho, shift = 0) {
  x <- vapply(seq_len(m), \(j) ar1(n, rho) + shift * (j == 1) / sqrt(1 - rho^2), numeric(n))
  c(rhat = rhat_rank(x), bulk = ess_bulk(x), tail = ess_tail(x))
}
summ <- function(lab, rho, shift, expect_ess, rule) {
  v <- t(replicate(R, one(rho, shift)))
  data.frame(case = lab, rho, shift, reps = R, rhat_median = median(v[, "rhat"]), rhat_max = max(v[, "rhat"]),
             bulk_median = median(v[, "bulk"]), tail_median = median(v[, "tail"]),
             expected_ess = expect_ess, rule = rule)
}
res <- rbind(
  summ("iid", 0, 0, S, "R-hat < 1.01 in every rep; bulk-ESS within 10% of S"),
  summ("ar1_0.9", 0.9, 0, S * 0.1 / 1.9, "bulk-ESS within 15% of S(1-rho)/(1+rho)"),
  summ("ar1_0.9_shifted", 0.9, 2, NA, "R-hat > 1.1 in every rep"))
res$pass <- c(res$rhat_max[1] < 1.01 && abs(res$bulk_median[1] / S - 1) < 0.10,
              abs(res$bulk_median[2] / res$expected_ess[2] - 1) < 0.15,
              NA)
v3 <- t(replicate(R, one(0.9, 2))); res$pass[3] <- all(v3[, "rhat"] > 1.1)
write.csv(res, "amendment_2026-10-08/checks/mcmc_diag_checks.csv", row.names = FALSE)
print(res, digits = 4)
cat(if (all(res$pass)) "ALL PASS\n" else "FAILURE\n")
