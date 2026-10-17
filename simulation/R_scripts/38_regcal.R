# Regression calibration with a known reliability (second amendment, Section 2;
# Carroll, Ruppert, Stefanski and Crainiceanu 2006).
#   source("R_scripts/38_regcal.R")
#
# Each error-prone column W*_k of X is replaced by the best linear predictor of the
# latent W_k given all columns of X, computed within each treatment group
# (by_group = TRUE) or once for the whole sample. Classical error independent of
# (W, V, A) gives Cov(W_k, X_j | A) = Cov(W*_k, X_j | A) for j != k and
# Var(W_k | A) = Var(W*_k | A) - s_k, with error variance s_k = (1 - r) Var(W*_k) over
# the whole sample, the same calibration information SIMEX uses. A non-positive
# Var(W_k | A) is floored at 1% of Var(W*_k | A) and counted in attr(, "n_flag").
regcal <- function(X, A, me_cols, reliability, by_group = TRUE) {
  X <- as.matrix(X)
  if (!all(me_cols %in% colnames(X))) {
    cli::cli_abort("Error-prone columns {.val {setdiff(me_cols, colnames(X))}} are not in X ({.val {colnames(X)}}).")
  }
  s <- (1 - reliability) * apply(X[, me_cols, drop = FALSE], 2, stats::var)
  grp <- if (by_group) A else rep(0L, nrow(X))
  out <- X; n_flag <- 0L
  for (g in unique(grp)) {
    i <- which(grp == g); Z <- X[i, , drop = FALSE]
    keep <- union(me_cols, colnames(Z)[apply(Z, 2, stats::var) > 0])   # drop columns constant in the group
    Z <- Z[, keep, drop = FALSE]
    mu <- colMeans(Z); Szz <- stats::cov(Z)
    Swz <- Szz[me_cols, , drop = FALSE]
    for (k in seq_along(me_cols)) {
      v_lat <- Szz[me_cols[k], me_cols[k]] - s[k]
      if (v_lat <= 0) { v_lat <- 0.01 * Szz[me_cols[k], me_cols[k]]; n_flag <- n_flag + 1L }
      Swz[k, me_cols[k]] <- v_lat
    }
    B <- Swz %*% qr.solve(Szz, tol = 1e-10)                            # |me| x p
    out[i, me_cols] <- sweep(sweep(Z, 2, mu) %*% t(B), 2, mu[me_cols], "+")
  }
  attr(out, "n_flag") <- n_flag
  out
}
