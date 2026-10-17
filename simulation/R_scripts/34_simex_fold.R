# Fold-calibrated SIMEX for cross-fit panel DR-DiD, method version 3.0.
#   source("R_scripts/30_core.R"); source("R_scripts/34_simex_fold.R")
#
# Post-production amendment (2026-10-06). Implements the modified estimator of
# post_prod_sim_steps.md Section 11.C.S4.2 and the inference of S5-S6.1. It is
# NOT the production method: 32_simex.R stays frozen and its saved results are
# never relabelled as evaluations of this file.
#
# Differences from 32_simex.R, each tied to the derivation it implements:
#   * one fold partition shared by every (lambda, b) fit, so each row's
#     calibration and simulation draws are well defined (S4.2 step 2);
#   * calibration s_hat_{-k} = (1 - r) Var(W*) on the TRAINING rows of fold k
#     ("fold"), or a supplied fixed covariance ("fixed"), or the production
#     full-sample plug-in ("full", whose expansion is unproved, S4.1);
#   * fold-normalised Hajek weights, averaged over folds (S4.2 step 3; the
#     pooled normalisation needs max-over-fold rates, A.3);
#   * the standard normals Z are drawn once, up front, and can be returned;
#   * calibration gradient D from the Gaussian score (S10), calibration
#     influence h_s (S12), combined influence eta = zeta + h_s D' (S13);
#   * pointwise, approximation-aware (S29), simultaneous-band (S23-S24) and
#     outer (S31) intervals, with the Peano constant kappa(t) (S27) computed
#     exactly for any grid.
#
# Generic matrix interface: X is the numeric design matrix, `me_cols` names or
# indexes its error-prone columns. Nothing here assumes the c1..c20 layout.

SIMEX_FOLD_VERSION <- "simex_fold_3.0"

# Parametric learner (logit / OLS on the raw columns). Used by the
# theorem-compatible checks, where it is correctly specified for every
# pseudo-law; not a registered production learner.
.learner_glm <- function(X, y, type) {
  Xd  <- cbind(1, X)
  fit <- if (type == "ps") glm.fit(Xd, y, family = binomial()) else lm.fit(Xd, y)
  b   <- fit$coefficients; b[is.na(b)] <- 0
  if (type == "ps") \(newX) plogis(drop(cbind(1, newX) %*% b))
  else              \(newX) drop(cbind(1, newX) %*% b)
}
FOLD_LEARNERS <- c(LEARNERS, list(glm = .learner_glm))

# Fold assignment stratified by A, as in crossfit_drdid().
.strat_folds <- function(A, K) {
  fold <- integer(length(A))
  fold[A == 1] <- sample(rep(seq_len(K), length.out = sum(A == 1)))
  fold[A == 0] <- sample(rep(seq_len(K), length.out = sum(A == 0)))
  fold
}

# Whole clusters to folds (independent clusters, Appendix cluster extension).
.cluster_folds <- function(cluster, K) {
  ug <- unique(cluster)
  sample(rep(seq_len(K), length.out = length(ug)))[match(cluster, ug)]
}

# Variance of the estimators whose influence values are the columns of Psi
# (full-sample scale). Rows are independent when `cluster` is NULL; otherwise
# clusters are independent and the variance uses centred cluster sums with the
# G / (G - 1) small-sample factor. Returns Var(estimate), already divided by n.
.if_vcov <- function(Psi, cluster = NULL) {
  Psi <- as.matrix(Psi); n <- nrow(Psi)
  if (is.null(cluster)) return(cov(Psi) / n)
  Sg <- rowsum(Psi, cluster); G <- nrow(Sg)
  Sc <- sweep(Sg, 2, colMeans(Sg))
  G / (G - 1) * crossprod(Sc) / n^2
}

# Hajek DR-DiD score on one held-out fold, normalised within the fold.
.fold_score <- function(A, dY, e, m0, trim) {
  clipped <- e < trim[1] | e > trim[2]
  e  <- pmin(pmax(e, trim[1]), trim[2])
  r  <- dY - m0
  w1 <- A / mean(A)
  w0 <- (1 - A) * e / (1 - e); w0 <- w0 / mean(w0)
  tau1 <- mean(w1 * r); tau0 <- mean(w0 * r)
  list(tau = tau1 - tau0, psi = w1 * (r - tau1) - w0 * (r - tau0),
       n_clipped = sum(clipped))
}

# Cross-fit DR-DiD (no SIMEX) with the same fold-normalised score; returns the
# estimate, full-sample-scale influence values and held-out propensities.
fold_drdid <- function(X, A, dY, learner, fold, trim = PS_TRIM) {
  L <- FOLD_LEARNERS[[learner]]; K <- max(fold); nk <- tabulate(fold, K); n <- length(A)
  tau <- numeric(K); psi <- numeric(n); e_hat <- numeric(n)
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k); ctrl <- tr[A[tr] == 0]
    e_hat[te] <- L(X[tr, , drop = FALSE], A[tr], "ps")(X[te, , drop = FALSE])
    m0 <- L(X[ctrl, , drop = FALSE], dY[ctrl], "or")(X[te, , drop = FALSE])
    fs <- .fold_score(A[te], dY[te], e_hat[te], m0, trim)
    tau[k] <- fs$tau; psi[te] <- fs$psi * n / (K * nk[k])
  }
  list(est = mean(tau), psi = psi, e_hat = e_hat)
}

# Quadratic-extrapolation weights c(t) = L' z(t) for grid `lam` (J+1 x length(t)).
.c_weights <- function(lam, t) {
  Zg <- cbind(1, lam, lam^2)
  Lm <- solve(crossprod(Zg), t(Zg))                # 3 x (J+1)
  t(Lm) %*% rbind(1, t, t^2)
}

# Peano constant kappa(t) = int_t^{lam_J} |K_t(u)| du with
# K_t(u) = 1/2 sum_j c_j(t) (lam_j - u)_+^2 (S26-S27). K_t is quadratic between
# knots, so each piece is split at its real roots and integrated exactly.
peano_kappa <- function(lam, t) {
  vapply(t, \(tt) {
    cw <- drop(.c_weights(lam, tt))
    knots <- sort(unique(c(tt, lam[lam > tt])))
    total <- 0
    for (i in seq_len(length(knots) - 1L)) {
      a <- knots[i]; b <- knots[i + 1L]
      act <- lam >= b                               # terms with (lam_j - u)_+ > 0 on (a, b)
      p0 <- sum(cw[act] * lam[act]^2) / 2; p1 <- -sum(cw[act] * lam[act]); p2 <- sum(cw[act]) / 2
      antideriv <- \(u) p0 * u + p1 * u^2 / 2 + p2 * u^3 / 3
      rts <- if (abs(p2) > 1e-14) {
        disc <- p1^2 - 4 * p2 * p0
        if (disc > 0) (-p1 + c(-1, 1) * sqrt(disc)) / (2 * p2) else numeric(0)
      } else if (abs(p1) > 1e-14) -p0 / p1 else numeric(0)
      cuts <- sort(c(a, rts[rts > a & rts < b], b))
      total <- total + sum(abs(diff(antideriv(cuts))))
    }
    total
  }, numeric(1))
}

# ---- the curve -------------------------------------------------------------
# Returns every quantity the amendment's output contract asks for. All
# influence values are on the full-sample scale: est - target ~ mean(IF).
simex_fold_curve <- function(X, A, dY, me_cols, reliability,
                             lambdas = c(0.5, 1, 2), B = 20L, K = 5L,
                             learner = "glmnet",
                             calibration = c("fold", "fixed", "full"),
                             s_fixed = NULL, trim = PS_TRIM, fold = NULL,
                             cluster = NULL, save_draws = FALSE) {
  calibration <- rlang::arg_match(calibration)
  L <- FOLD_LEARNERS[[learner]]
  if (is.null(L)) cli::cli_abort("Unknown learner {.val {learner}}; registered: {.val {names(FOLD_LEARNERS)}}.")
  X <- as.matrix(X); n <- nrow(X)
  me <- if (is.character(me_cols)) match(me_cols, colnames(X)) else as.integer(me_cols)
  if (anyNA(me)) cli::cli_abort("Error-prone columns {.val {me_cols[is.na(me)]}} are not in X ({.val {colnames(X)}}).")
  p <- length(me); J <- length(lambdas); lam <- c(0, lambdas)
  if (is.null(fold)) fold <- if (is.null(cluster)) .strat_folds(A, K) else .cluster_folds(cluster, K)
  K <- max(fold); nk <- tabulate(fold, K)
  W <- X[, me, drop = FALSE]

  # Calibration: one covariance vector per fold, and the rowwise influence h_s.
  if (calibration == "fixed") {
    if (length(s_fixed) != p || any(s_fixed <= 0)) stop("s_fixed must be p positive variances.")
    S <- matrix(s_fixed, K, p, byrow = TRUE); h_s <- matrix(0, n, p)
  } else {
    v_full <- apply(W, 2, var)
    S <- if (calibration == "fold")
      t(vapply(seq_len(K), \(k) (1 - reliability) * apply(W[fold != k, , drop = FALSE], 2, var), numeric(p)))
    else matrix((1 - reliability) * v_full, K, p, byrow = TRUE)
    S <- matrix(S, K, p)
    Wc  <- sweep(W, 2, colMeans(W))
    h_s <- (1 - reliability) * sweep(Wc^2, 2, v_full)            # S12
  }
  s_row <- S[fold, , drop = FALSE]                                # n x p

  Z <- array(rnorm(n * p * J * B), c(n, p, J, B))                 # drawn once
  zeta <- matrix(0, n, J + 1); Q <- numeric(J + 1)
  Qjb  <- matrix(NA_real_, J + 1, B); Dhat <- matrix(0, J + 1, p)
  Drow <- array(0, c(n, J + 1, p))                                # rowwise phi*H, for D's MCSE
  n_clip <- matrix(0, J + 1, B)

  for (j in 0:J) {
    nb <- if (j == 0) 1L else B
    for (b in seq_len(nb)) {
      Xp <- X
      if (j > 0) Xp[, me] <- W + Z[, , j, b] * sqrt(lam[j + 1] * s_row)
      tau_k <- numeric(K); phi <- numeric(n)
      for (k in seq_len(K)) {
        tr <- which(fold != k); te <- which(fold == k)
        e_te  <- L(Xp[tr, , drop = FALSE], A[tr], "ps")(Xp[te, , drop = FALSE])
        ctrl  <- tr[A[tr] == 0]
        m0_te <- L(Xp[ctrl, , drop = FALSE], dY[ctrl], "or")(Xp[te, , drop = FALSE])
        fs <- .fold_score(A[te], dY[te], e_te, m0_te, trim)
        tau_k[k] <- fs$tau
        phi[te]  <- fs$psi * n / (K * nk[k])                      # full-sample IF scale
        n_clip[j + 1, b] <- n_clip[j + 1, b] + fs$n_clipped
      }
      Qjb[j + 1, b] <- mean(tau_k)
      zeta[, j + 1] <- zeta[, j + 1] + phi / nb
      if (j > 0) {
        H <- (Z[, , j, b]^2 - 1) / (2 * s_row)                    # S9
        Drow[, j + 1, ] <- Drow[, j + 1, ] + phi * H / nb
      }
    }
    Q[j + 1] <- mean(Qjb[j + 1, seq_len(nb)])
  }
  Dhat <- apply(Drow, c(2, 3), mean)                              # (J+1) x p, S10
  eta  <- zeta + h_s %*% t(Dhat)                                  # S13
  Lm   <- solve(crossprod(cbind(1, lam, lam^2)), t(cbind(1, lam, lam^2)))
  V_zeta <- .if_vcov(zeta, cluster); V_eta <- .if_vcov(eta, cluster)

  out <- list(
    method_version = SIMEX_FOLD_VERSION, calibration_design = calibration,
    learner = learner, reliability_ref = reliability, n = n, K = K, B = B,
    lam = lam, fold = fold, S_folds = S, Q = Q, Qjb = Qjb, coef = drop(Lm %*% Q),
    zeta = zeta, eta = eta, D = Dhat, D_mcse = apply(Drow, c(2, 3), sd) / sqrt(n),
    h_s = h_s, cluster = cluster, V_zeta = V_zeta, V_eta = V_eta,
    # Omega is the covariance of sqrt(n) (a_hat - a), so Var(a_hat) = Omega / n.
    Omega_fixed = n * Lm %*% V_zeta %*% t(Lm), Omega = n * Lm %*% V_eta %*% t(Lm),
    clip_frac = n_clip / n)
  if (save_draws) out$Z <- Z
  out
}

# ---- inference ---------------------------------------------------------------
# Pointwise interval for q2(t) (fixed-calibration and calibration-aware SEs),
# widened by M * kappa(t) for the true ATT (S29). M = 0 is exact-quadratic
# extrapolation; M is in outcome units per (lambda unit)^3 of the reference
# covariance and is a sensitivity assumption, not an estimate.
simex_fold_point <- function(cv, t = -1, M = 0, alpha = 0.05) {
  cw  <- drop(.c_weights(cv$lam, t))
  est <- sum(cw * cv$Q)
  se_fixed <- sqrt(drop(t(cw) %*% cv$V_zeta %*% cw))
  se_cal   <- sqrt(drop(t(cw) %*% cv$V_eta %*% cw))
  zq <- qnorm(1 - alpha / 2); bias_allow <- M * peano_kappa(cv$lam, t)
  c(est = est, se_fixed = se_fixed, se = se_cal, kappa = peano_kappa(cv$lam, t),
    M = M, ci_lo = est - zq * se_cal - bias_allow, ci_hi = est + zq * se_cal + bias_allow)
}

.kappa_cache <- new.env()

# Simultaneous band over t in [t_lo, t_hi] (S23) and the range / outer
# intervals (S24, S31). eps > 0 ridges a near-singular Omega (S5).
# Extrema are found on a dense grid and refined by optimize() on the bracketing
# cells; this is numerical, not a certified bound. `conservative` is the
# closed-form widening sup_t w_n(t) <= k sqrt(3 lambda_max(Omega) / n).
simex_fold_range <- function(cv, t_lo = -1, t_hi = 0, M = 0, alpha = 0.05,
                             eps = 0, calibration_aware = TRUE, grid = 2001L) {
  Om <- if (calibration_aware) cv$Omega else cv$Omega_fixed
  Om <- Om + eps * diag(3)
  k_a <- sqrt(qchisq(1 - alpha, 3)); a <- cv$coef; n <- cv$n
  q2 <- \(t) a[1] + a[2] * t + a[3] * t^2
  w  <- \(t) k_a * sqrt(pmax(colSums(rbind(1, t, t^2) * (Om %*% rbind(1, t, t^2))), 0) / n)
  bM <- \(t) if (M == 0) 0 * t else M * peano_kappa(cv$lam, t)
  lower <- \(t) q2(t) - w(t) - bM(t); upper <- \(t) q2(t) + w(t) + bM(t)
  tg <- seq(t_lo, t_hi, length.out = grid)
  # kappa on the grid depends only on (lam, grid), so it is computed once per
  # session; optimize() below still evaluates exact kappa at off-grid points.
  key <- paste(c(cv$lam, t_lo, t_hi, grid), collapse = ",")
  if (is.null(.kappa_cache[[key]])) .kappa_cache[[key]] <- peano_kappa(cv$lam, tg)
  kg <- .kappa_cache[[key]]
  ext <- function(f, maximum) {
    v <- q2(tg) + (if (maximum) 1 else -1) * (w(tg) + M * kg)   # f on the grid
    i <- if (maximum) which.max(v) else which.min(v)
    br <- c(tg[max(i - 1L, 1L)], tg[min(i + 1L, grid)])
    best <- v[i]
    if (br[2] > br[1]) {
      cand <- optimize(f, br, maximum = maximum)$objective
      best <- if (maximum) max(best, cand) else min(best, cand)
    }
    best
  }
  vtx <- if (abs(a[3]) > 0) -a[2] / (2 * a[3]) else NA_real_
  tq  <- c(t_lo, t_hi, if (is.finite(vtx) && vtx > t_lo && vtx < t_hi) vtx)
  wbar <- k_a * sqrt(3 * max(eigen(Om, symmetric = TRUE, only.values = TRUE)$values, 0) / n)
  c(range_lo = min(q2(tq)), range_hi = max(q2(tq)),               # exact extrema of q2-hat
    ci_lo = ext(lower, FALSE), ci_hi = ext(upper, TRUE),
    ci_lo_conservative = min(q2(tq)) - wbar - M * max(kg),
    ci_hi_conservative = max(q2(tq)) + wbar + M * max(kg),
    M = M, k_alpha = k_a)
}

# Simultaneous-band event for a known population curve (simulation checks):
# TRUE when |q2_hat(t) - q2(t)| <= w_n(t) on the whole grid.
simex_fold_band_covers <- function(cv, coef_true, t_lo = -1, t_hi = 0,
                                   alpha = 0.05, calibration_aware = TRUE, grid = 401L) {
  Om <- if (calibration_aware) cv$Omega else cv$Omega_fixed
  tg <- seq(t_lo, t_hi, length.out = grid); Zt <- rbind(1, tg, tg^2)
  w  <- sqrt(qchisq(1 - alpha, 3)) * sqrt(colSums(Zt * (Om %*% Zt)) / cv$n)
  all(abs(drop(crossprod(Zt, cv$coef - coef_true))) <= w)
}

# ---- data-frame wrappers (simulation design) -----------------------------------
# Known reliability r: extrapolate to t = -1. Returns .out() with the curve
# attached, so a harness can save all diagnostics.
simex_fold_drdid <- function(d, reliability, learner = "glmnet",
                             calibration = "fold", s_fixed = NULL, M = 0,
                             lambdas = c(0.5, 1, 2), B = 20L, K = 5L,
                             covs = ME_X_COVS) {
  X  <- design_matrix(d)
  cv <- simex_fold_curve(X, d$A, d$dY, covs, reliability, lambdas, B, K,
                         learner, calibration, s_fixed)
  pt <- simex_fold_point(cv, -1, M)
  out <- .out(pt[["est"]], pt[["se"]], pt[["ci_lo"]], pt[["ci_hi"]])
  attr(out, "se_fixed") <- pt[["se_fixed"]]
  attr(out, "curve") <- cv
  out
}
