# MCMC convergence diagnostics of Vehtari, Gelman, Simpson, Carpenter and
# Buerkner (2021, Bayesian Analysis 16:667-718): rank-normalised split R-hat
# (maximum of bulk and folded) and bulk- and tail-ESS. Follows the reference
# implementation in the `posterior` package, which is not installed here.
#   source("R_scripts/37_mcmc_diag.R")
# Input: an iterations x chains matrix of draws of one scalar.

# Autocovariance by FFT, "biased" (divided by N), as recommended by Geyer (1992).
.acov <- function(x) {
  N <- length(x); M <- stats::nextn(N); yc <- c(x - mean(x), rep(0, 2 * M - N))
  tr <- stats::fft(yc)
  Re(stats::fft(Conj(tr) * tr, inverse = TRUE))[seq_len(N)] / (N^2 * 2)
}

.split_chains <- function(x) {
  x <- as.matrix(x); n <- nrow(x); h <- n %/% 2
  if (n %% 2) x <- x[-(h + 1), , drop = FALSE]     # drop the middle draw when n is odd
  cbind(x[seq_len(h), , drop = FALSE], x[h + seq_len(h), , drop = FALSE])
}

# Rank-normalise jointly over all chains (Blom offset c = 3/8).
.z_scale <- function(x) {
  r <- rank(as.vector(x), ties.method = "average"); S <- length(r)
  matrix(stats::qnorm((r - 3 / 8) / (S - 2 * 3 / 8 + 1)), nrow = nrow(x))
}

.rhat_basic <- function(x) {
  n <- nrow(x); cm <- colMeans(x); cv <- apply(x, 2, stats::var)
  sqrt(((n * stats::var(cm)) / mean(cv) + n - 1) / n)
}

.ess_basic <- function(x) {
  x <- as.matrix(x); m <- ncol(x); n <- nrow(x)
  if (n < 3 || any(!is.finite(x)) || stats::var(as.vector(x)) == 0) return(NA_real_)
  acov <- vapply(seq_len(m), \(i) .acov(x[, i]), numeric(n))
  if (n == 1) acov <- matrix(acov, nrow = 1)
  chain_mean <- colMeans(x)
  mean_var <- mean(acov[1, ]) * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n + (if (m > 1) stats::var(chain_mean) else 0)
  rho <- numeric(n); t <- 0
  rho_even <- 1; rho[1] <- rho_even
  rho_odd <- 1 - (mean_var - mean(acov[2, ])) / var_plus; rho[2] <- rho_odd
  while (t < n - 5 && !is.nan(rho_even + rho_odd) && rho_even + rho_odd > 0) {   # Geyer: initial positive
    t <- t + 2
    rho_even <- 1 - (mean_var - mean(acov[t + 1, ])) / var_plus
    rho_odd  <- 1 - (mean_var - mean(acov[t + 2, ])) / var_plus
    if (rho_even + rho_odd >= 0) { rho[t + 1] <- rho_even; rho[t + 2] <- rho_odd }
  }
  max_t <- t
  if (rho_even > 0) rho[max_t + 1] <- rho_even
  t <- 0
  while (t <= max_t - 4) {                                                         # initial monotone
    t <- t + 2
    if (rho[t + 1] + rho[t + 2] > rho[t - 1] + rho[t]) {
      rho[t + 1] <- (rho[t - 1] + rho[t]) / 2; rho[t + 2] <- rho[t + 1]
    }
  }
  S <- m * n
  tau_hat <- max(-1 + 2 * sum(rho[seq_len(max_t)]) + rho[max_t + 1], 1 / log10(S))
  S / tau_hat
}

rhat_rank <- function(x) {
  s <- .split_chains(x)
  max(.rhat_basic(.z_scale(s)), .rhat_basic(.z_scale(abs(s - stats::median(s)))))
}
ess_bulk <- function(x) .ess_basic(.z_scale(.split_chains(x)))
ess_tail <- function(x) {
  s <- .split_chains(x)
  min(vapply(c(0.05, 0.95), \(p) .ess_basic((s <= stats::quantile(s, p)) + 0), numeric(1)))
}
