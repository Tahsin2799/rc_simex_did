# Population quantities of the v2.0 DGP, by numerical integration.
#   source("R_scripts/1_sim_param.R"); source("R_scripts/35_population_targets.R")
#
# Production scores estimators against the replicate's REALISED treated-sample
# effect mean(theta_i[A == 1]). These give the population counterparts
# (post_prod_sim_steps.md Section 3.2 and A.7), including the effect of the
# truncation and rounding in data_gen_v2(), so both targets can be reported.

# E[min(1, max(0, X))], X ~ N(mu, sd).
.e_clamp01 <- function(mu, sd)
  integrate(\(x) x * dnorm(x, mu, sd), 0, 1, rel.tol = 1e-12)$value + pnorm(1, mu, sd, lower.tail = FALSE)
# Moments of max(lo, round(X)), X ~ N(mu, sd).
.round_trunc_moments <- function(mu, sd, lo) {
  k  <- lo:ceiling(mu + 40 * sd)
  pk <- pnorm(k + 0.5, mu, sd) - pnorm(k - 0.5, mu, sd)
  pk[1] <- pnorm(lo + 0.5, mu, sd)                    # everything that rounds to <= lo
  c(mean = sum(k * pk), var = sum(k^2 * pk) - sum(k * pk)^2)
}

# Population ATT. Effect modifiers c11..c15 are independent of A, so
# E[theta | A = 1] = theta (1 + sum_k g_k E z_k); z11, z13 have mean zero exactly,
# the truncated z12, z15 and rounded/truncated z14 do not.
population_att <- function(dgp_family, params = default_params("low")) {
  if (dgp_family == "linear") return(params$theta)
  m <- params$cov_moments; g <- params$gamma
  ez12 <- (.e_clamp01(m$c12_mean, m$c12_sd) - m$c12_mean) / m$c12_sd
  ez15 <- (.e_clamp01(m$c15_mean, m$c15_sd) - m$c15_mean) / m$c15_sd
  ez14 <- (.round_trunc_moments(m$c14_mean, m$c14_sd, 0L)[["mean"]] - m$c14_mean) / m$c14_sd
  params$theta * (1 + g$g12 * ez12 + g$g14 * ez14 + g$g15 * ez15)
}

# Population variance of each covariate that receives classical error, and the
# implied error variance at reliability r. me_x_classical() scales its error by
# the SAMPLE sd of the latent covariate; this is the fixed-population version.
population_me_variance <- function(reliability, covs = c("c1", "c3", "c4", "c6", "c8"),
                                   params = default_params("low")) {
  m <- params$cov_moments
  v <- vapply(covs, \(k) if (k == "c4") .round_trunc_moments(m$c4_mean, m$c4_sd, 2L)[["var"]]
              else m[[paste0(k, "_sd")]]^2, numeric(1))
  v * (1 / reliability - 1)
}
