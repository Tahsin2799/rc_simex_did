# Exact full-panel enumeration for the binary-support completion counterexample.
# No Monte Carlo or production estimator is used; dplyr supplies only if_else and case_when.
panel <- expand.grid(w = c(-1, 1), a = 0:1, y0 = 0:1, y1 = 0:1)
panel$e <- dplyr::if_else(panel$w == -1, 0.1, 0.9)
panel$p_y0 <- dplyr::case_when(
  panel$a == 0 & panel$w == -1 ~ 0,
  panel$a == 0                 ~ 0.1,
  panel$w == -1                ~ 0.1,
  .default                     = 0.9
)
panel$p_y1 <- dplyr::case_when(
  panel$a == 0 & panel$w == -1 ~ 0.8,
  panel$a == 0                 ~ 0,
  panel$w == -1                ~ 0.9,
  .default                     = 0.8
)
panel$prob <- 0.5 * dbinom(panel$a, 1, panel$e) *
  dbinom(panel$y0, 1, panel$p_y0) * dbinom(panel$y1, 1, panel$p_y1)
panel$d <- panel$y1 - panel$y0
weighted <- function(value, prob) sum(value * prob) / sum(prob)
conditional <- function(var, a, w) {
  use <- panel$a == a & panel$w == w
  weighted(panel[[var]][use], panel$prob[use])
}
mu0 <- vapply(c(-1, 1), \(w) conditional("d", 0, w), numeric(1))
mu1 <- vapply(c(-1, 1), \(w) conditional("d", 1, w), numeric(1))
b1 <- vapply(c(-1, 1), \(w) conditional("y0", 1, w), numeric(1))
r_latent <- mu0 + b1
smoothed <- function(z) {
  posterior_mass <- panel$prob * dnorm(z - panel$w)
  control <- panel$a == 0
  treated <- panel$a == 1
  mu <- weighted(panel$d[control], posterior_mass[control])
  baseline <- weighted(panel$y0[treated], posterior_mass[treated])
  c(mu_control = mu, baseline_treated = baseline,
    required_untreated_followup = mu + baseline)
}
at_zero <- smoothed(0)
treated_w_mass <- vapply(c(-1, 1), \(w)
  sum(panel$prob[panel$a == 1 & panel$w == w]), numeric(1))
latent_att <- weighted(mu1 - mu0, treated_w_mass)
checks <- data.frame(
  check = c("joint_probability_sum", "latent_cpt_difference_w_minus1",
            "latent_cpt_difference_w_plus1", "latent_binary_r_minus1",
            "latent_binary_r_plus1", "latent_att", "smoothed_control_trend_at_zero",
            "smoothed_treated_baseline_at_zero", "smoothed_binary_r_at_zero"),
  actual = c(sum(panel$prob), mu1 - mu0, r_latent, latent_att, at_zero),
  expected = c(1, 0, 0, 0.9, 0.8, 0, 0.71, 0.82, 1.53)
)
checks$absolute_error <- abs(checks$actual - checks$expected)
if (max(checks$absolute_error) >= 1e-12) cli::cli_abort("Binary-support example: largest error {signif(max(checks$absolute_error), 3)}.")
if (any(r_latent < 0 | r_latent > 1)) cli::cli_abort("Latent completion probabilities leave [0, 1]: range {signif(range(r_latent), 4)}.")
if (!all(vapply(c(-0.05, 0, 0.05), \(z) smoothed(z)[3] > 1, logical(1)))) {
  cli::cli_abort("Smoothing did not push the completion probability above 1 near W* = 0.")
}
write.csv(checks, "simulation/audit_2026-10-06/binary_support_checks.csv",
          row.names = FALSE)
print(checks, row.names = FALSE, digits = 12)
cat("PASS: binary CPT is feasible at the latent covariate and fails after smoothing.\n")
