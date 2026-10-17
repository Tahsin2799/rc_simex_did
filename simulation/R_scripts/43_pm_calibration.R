# Posterior-mean calibration and SIMEX applied to it (sixth amendment pilot, amendment_2026-10-14b).
#   source("R_scripts/30_core.R"); source("R_scripts/34_simex_fold.R"); source("R_scripts/43_pm_calibration.R")
# Each error-prone covariate is replaced by its posterior mean E[W_c | W*_c, A] within the treatment group,
# under classical Gaussian error with variance s_c = (1 - r) Var(W*_c). The latent law of W_c given A is
# estimated nonparametrically by the grid NPMLE of a Gaussian location mixture (Kiefer-Wolfowitz; EM), so
# the posterior mean equals Tweedie's formula y + s (log p)'(y) under the fitted law. Coordinates are
# calibrated one at a time, which is exact when they are independent given A (as in the shape designs).
# Everything is cross-fitted: latent laws and s are estimated on the training folds.
# SIMEX on top: noise sqrt(lambda s) Z is added to W*, the noisy covariates are calibrated with error
# variance (1 + lambda) s under the same latent law (the latent law does not depend on lambda), the DR-DiD
# is refitted, and the curve is extrapolated quadratically to lambda = -1.

# Grid NPMLE of the mixing law for y = W + N(0, s). Support: a mix of quantile and equally spaced points.
npmle_fit <- function(y, s, K = 300L, iters = 400L) {
  u <- sort(unique(c(quantile(y, seq(0, 1, length.out = K %/% 2), names = FALSE),
                     seq(min(y), max(y), length.out = K %/% 2))))
  lg <- dnorm(outer(y, u, "-"), sd = sqrt(s), log = TRUE)
  L <- exp(lg - apply(lg, 1, max))                                   # row-scaled likelihoods
  w <- rep(1 / length(u), length(u))
  for (it in seq_len(iters)) w <- w * drop(crossprod(L, 1 / drop(L %*% w))) / length(y)
  list(u = u, w = w)
}
# Posterior mean of W given y = W + N(0, e) under a fitted mixing law (log-sum-exp for the tails).
pm_mean <- function(y, fit, e) {
  lg <- dnorm(outer(y, fit$u, "-"), sd = sqrt(e), log = TRUE) + rep(log(pmax(fit$w, 1e-300)), each = length(y))
  K <- exp(lg - apply(lg, 1, max))
  drop(K %*% fit$u) / rowSums(K)
}

# Cross-fitted DR-DiD with posterior-mean calibration, and SIMEX on top.
# Returns the curve Q (lambda = 0 is posterior-mean RC), the extrapolated estimate, and SEs from the
# influence values that treat the calibration as fixed.
pm_simex_curve <- function(X, A, dY, me_cols, reliability, learner = "glmnet", fold,
                           lambdas = c(0.5, 1, 2), B = 20L, trim = PS_TRIM, K_grid = 300L) {
  L <- FOLD_LEARNERS[[learner]]; X <- as.matrix(X); n <- nrow(X)
  me <- match(me_cols, colnames(X)); p <- length(me); J <- length(lambdas); lam <- c(0, lambdas)
  K <- max(fold); nk <- tabulate(fold, K); W <- X[, me, drop = FALSE]
  S <- t(vapply(seq_len(K), \(k) (1 - reliability) * apply(W[fold != k, , drop = FALSE], 2, var), numeric(p)))
  S <- matrix(S, K, p); s_row <- S[fold, , drop = FALSE]
  # Latent laws per (training folds k, group a, coordinate c), fitted on the unperturbed W*.
  fits <- lapply(seq_len(K), \(k) lapply(0:1, \(a) lapply(seq_len(p), \(c) {
    idx <- which(fold != k & A == a); npmle_fit(W[idx, c], S[k, c], K_grid) })))
  Z <- array(rnorm(n * p * J * B), c(n, p, J, B))
  zeta <- matrix(0, n, J + 1); Q <- numeric(J + 1)
  for (j in 0:J) {
    nb <- if (j == 0) 1L else B; Qb <- numeric(nb)
    for (b in seq_len(nb)) {
      Wn <- if (j == 0) W else W + Z[, , j, b] * sqrt(lam[j + 1] * s_row)
      tau_k <- numeric(K); phi <- numeric(n)
      for (k in seq_len(K)) {
        Xh <- X                                                      # calibrate every unit with fold k's laws
        for (c in seq_len(p)) for (a in 0:1) {
          ia <- which(A == a); Xh[ia, me[c]] <- pm_mean(Wn[ia, c], fits[[k]][[a + 1]][[c]], (1 + lam[j + 1]) * S[k, c])
        }
        tr <- which(fold != k); te <- which(fold == k); ctrl <- tr[A[tr] == 0]
        e_te <- L(Xh[tr, , drop = FALSE], A[tr], "ps")(Xh[te, , drop = FALSE])
        m0_te <- L(Xh[ctrl, , drop = FALSE], dY[ctrl], "or")(Xh[te, , drop = FALSE])
        fs <- .fold_score(A[te], dY[te], e_te, m0_te, trim)
        tau_k[k] <- fs$tau; phi[te] <- fs$psi * n / (K * nk[k])
      }
      Qb[b] <- mean(tau_k); zeta[, j + 1] <- zeta[, j + 1] + phi / nb
    }
    Q[j + 1] <- mean(Qb)
  }
  Zg <- cbind(1, lam, lam^2); cw <- drop(t(solve(crossprod(Zg), t(Zg))) %*% c(1, -1, 1))
  V <- .if_vcov(zeta)
  list(lam = lam, Q = Q, est = sum(cw * Q), se_fixed = sqrt(drop(t(cw) %*% V %*% cw)),
       pm_rc = Q[1], pm_rc_se_fixed = sqrt(V[1, 1]))
}
