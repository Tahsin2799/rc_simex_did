# Fold-calibrated regression calibration with a calibration-aware influence
# function (third amendment, amendment_2026-10-11/AMENDMENT_3.md).
#   source("R_scripts/30_core.R"); source("R_scripts/34_simex_fold.R"); source("R_scripts/39_regcal_fold.R")
#
# Group-wise regression calibration replaces the error-prone columns W* of
# Z = (W*, V) by their best linear predictor given Z within each treatment
# group a. With S_p = diag(s) on the error-prone coordinates (0 elsewhere),
# Sigma_a = Cov(Z | A = a) and mu_a = E[Z | A = a], the calibrated covariates are
#   Xhat = Z - S_p Sigma_a^{-1} (Z - mu_a),
# which equals 38_regcal.R's map. The calibration parameters
# eta = (mu_0, mu_1, Sigma_0, Sigma_1, s) are estimated on the training rows of
# each fold ("fold" calibration), so the held-out rows are independent of the
# calibration, as in the fold-calibrated SIMEX (Theorem 1 of the manuscript).
#
# The ATT functional at calibration eta is Q(eta) = E[A{D - m_eta(Xhat_1)}] / pi,
# where m_eta is the control outcome regression on the calibrated covariates.
# Because Xhat is an invertible affine map of Z within each group, m_eta(x)
# = m*_0(h_0^{-1}(x)), and the directional derivative of Q along delta is
#   dQ[delta] = -E[A G' {dh_1(Z)[delta] - dh_0(z')[delta]}] / pi,
# with G = grad_x m_eta(Xhat_1) and z' = h_0^{-1}(Xhat_1) the control-calibration
# preimage of a treated unit. dh_a(z)[delta] = -dR_a (z - mu_a) + R_a dmu_a with
# R_a = S_p Sigma_a^{-1} and dR_a = dS_p Sigma_a^{-1} - S_p Sigma_a^{-1} dSigma_a Sigma_a^{-1}.
# The calibration contribution of unit j is c_j = dQ[IF_eta(O_j)], and the
# combined influence value is chi_j = phi_j + c_j.
#
# G is the gradient of the FITTED control regression, taken by central
# differences of its predictions. That is exact for the parametric learner and
# for the lasso dictionary of 30_core.R (raw columns plus squares, a quadratic
# polynomial); it is not valid for tree learners, which are not used here.

RC_FOLD_VERSION <- "regcal_fold_1.0"

# Calibration parameters from the rows `rows`; s is pooled over those rows.
.rc_params <- function(Z, A, me, reliability, rows = seq_len(nrow(Z))) {
  Zr <- Z[rows, , drop = FALSE]; Ar <- A[rows]
  s <- (1 - reliability) * apply(Zr[, me, drop = FALSE], 2, var)
  grp <- lapply(0:1, \(a) {
    Za <- Zr[Ar == a, , drop = FALSE]
    Sig <- cov(Za)
    # Columns constant within the group (rare dummies) are dropped from the predictor, as
    # 38_regcal.R does; their rows and columns of the inverse are zero.
    keep <- union(me, which(diag(Sig) > 0))
    Si <- matrix(0, ncol(Za), ncol(Za)); Si[keep, keep] <- solve(Sig[keep, keep, drop = FALSE])
    v_obs <- diag(Sig)[me]
    s_eff <- pmin(s, 0.99 * v_obs)                     # latent variance floored at 1% (as 38_regcal.R)
    list(mu = colMeans(Za), Sig = Sig, Si = Si, s_eff = s_eff, n_flag = sum(s >= 0.99 * v_obs))
  })
  list(s = s, g = grp, me = me)
}

# Calibrated covariates for all rows of Z, given parameters `cp`.
.rc_apply <- function(Z, A, cp) {
  me <- cp$me; out <- Z
  for (a in 0:1) {
    i <- which(A == a); if (!length(i)) next
    ga <- cp$g[[a + 1]]
    U <- sweep(Z[i, , drop = FALSE], 2, ga$mu)
    out[i, me] <- Z[i, me, drop = FALSE] - sweep((U %*% ga$Si)[, me, drop = FALSE], 2, ga$s_eff, "*")
  }
  out
}

# Control-calibration preimage z' of treated calibrated points: V' = V and
# (I - D_s P_mm) w' = what - D_s P_mm mu0_m + D_s P_mv (v - mu0_v), P = Sigma_0^{-1}.
.rc_preimage <- function(Xhat1, cp) {
  me <- cp$me; g0 <- cp$g[[1]]; P <- g0$Si; ex <- setdiff(seq_len(ncol(Xhat1)), me)
  Ds <- diag(g0$s_eff, length(me))
  Mlhs <- diag(length(me)) - Ds %*% P[me, me, drop = FALSE]
  rhs <- t(Xhat1[, me, drop = FALSE]) - drop(Ds %*% P[me, me, drop = FALSE] %*% g0$mu[me])
  if (length(ex)) rhs <- rhs + Ds %*% P[me, ex, drop = FALSE] %*% (t(Xhat1[, ex, drop = FALSE]) - g0$mu[ex])
  zp <- Xhat1; zp[, me] <- t(solve(Mlhs, rhs)); zp
}

# Gradient of a fitted predictor in the error-prone coordinates (central differences).
.rc_grad <- function(pred, X, me, h = 0.01) {
  G <- matrix(0, nrow(X), ncol(X))
  for (k in me) {
    Xp <- X; Xm <- X; Xp[, k] <- Xp[, k] + h; Xm[, k] <- Xm[, k] - h
    G[, k] <- (pred(Xp) - pred(Xm)) / (2 * h)
  }
  G
}

# Generic directional derivative dQ[delta] from the accumulated pieces; delta is
# list(dmu = list(a0, a1), dSig = list(a0, a1), ds). Used by the checks.
rc_dQ <- function(parts, delta) {
  cp <- parts$cp; me <- cp$me; p <- length(cp$g[[1]]$mu)
  tot <- 0
  for (a in 0:1) {
    ga <- cp$g[[a + 1]]; Ua <- parts$U[[a + 1]]
    Sp <- matrix(0, p, p); Sp[cbind(me, me)] <- ga$s_eff
    dSp <- matrix(0, p, p); dSp[cbind(me, me)] <- delta$ds
    dR <- dSp %*% ga$Si - Sp %*% ga$Si %*% delta$dSig[[a + 1]] %*% ga$Si
    Ta <- -sum(diag(dR %*% Ua)) + sum(parts$gsum * (Sp %*% ga$Si %*% delta$dmu[[a + 1]]))
    tot <- tot + (if (a == 1) Ta else -Ta)
  }
  -tot / parts$n1
}

# Rowwise calibration contributions c_j = dQ[IF_eta(O_j)] (vectorised rc_dQ).
.rc_contrib <- function(Z, A, parts, reliability) {
  cp <- parts$cp; me <- cp$me; n <- nrow(Z); p <- ncol(Z)
  Wm <- Z[, me, drop = FALSE]; wbar <- colMeans(Wm); v <- apply(Wm, 2, var)
  ds <- (1 - reliability) * sweep(sweep(Wm, 2, wbar)^2, 2, v)          # n x |me|, IF of s
  cj <- numeric(n)
  for (a in 0:1) {
    ga <- cp$g[[a + 1]]; Ua <- parts$U[[a + 1]]; pia <- mean(A == a)
    Sp <- matrix(0, p, p); Sp[cbind(me, me)] <- ga$s_eff
    SiU <- ga$Si %*% Ua
    trS <- drop(ds %*% diag(SiU)[me])                                  # tr(dS_p Sigma^-1 U)
    K <- SiU %*% Sp %*% ga$Si; Ks <- (K + t(K)) / 2
    ua <- sweep(Z, 2, ga$mu)
    qf <- rowSums((ua %*% Ks) * ua) - sum(ga$Sig * Ks)                 # u'Ku - tr(Sigma K)
    ind <- (A == a) / pia
    lin <- drop(ua %*% t(Sp %*% ga$Si) %*% parts$gsum)                 # gsum' R_a u
    Ta <- -(trS - ind * qf) + ind * lin
    cj <- cj + (if (a == 1) Ta else -Ta)
  }
  -cj / parts$n1
}

# The estimator. calibration = "fold" (training-row calibration, the version
# covered by the expansion) or "full" (one full-sample calibration, as the
# comparator of the second amendment). Returns the estimate, both SEs, the
# influence values and diagnostics.
regcal_fold_drdid <- function(X, A, dY, me_cols, reliability, learner = "glmnet", fold,
                              calibration = c("fold", "full"), trim = PS_TRIM, h = 0.01) {
  calibration <- rlang::arg_match(calibration)
  Lr <- FOLD_LEARNERS[[learner]]; X <- as.matrix(X); n <- nrow(X); p <- ncol(X)
  me <- if (is.character(me_cols)) match(me_cols, colnames(X)) else as.integer(me_cols)
  K <- max(fold); nk <- tabulate(fold, K)
  cp_full <- .rc_params(X, A, me, reliability)
  tau <- numeric(K); phi <- numeric(n); e_hat <- numeric(n); n_flag <- 0L
  U <- list(matrix(0, p, p), matrix(0, p, p)); gsum <- numeric(p); n1 <- 0
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k)
    cp <- if (calibration == "fold") .rc_params(X, A, me, reliability, tr) else cp_full
    n_flag <- n_flag + cp$g[[1]]$n_flag + cp$g[[2]]$n_flag
    Xh <- .rc_apply(X, A, cp); colnames(Xh) <- colnames(X)
    ctrl <- tr[A[tr] == 0]
    e_hat[te] <- Lr(Xh[tr, , drop = FALSE], A[tr], "ps")(Xh[te, , drop = FALSE])
    mf <- Lr(Xh[ctrl, , drop = FALSE], dY[ctrl], "or")
    fs <- .fold_score(A[te], dY[te], e_hat[te], mf(Xh[te, , drop = FALSE]), trim)
    tau[k] <- fs$tau; phi[te] <- fs$psi * n / (K * nk[k])
    t1 <- te[A[te] == 1]
    if (length(t1)) {
      X1 <- Xh[t1, , drop = FALSE]
      G <- .rc_grad(mf, X1, me, h)
      zp <- .rc_preimage(X1, cp)
      U[[2]] <- U[[2]] + crossprod(sweep(X[t1, , drop = FALSE], 2, cp$g[[2]]$mu), G)
      U[[1]] <- U[[1]] + crossprod(sweep(zp, 2, cp$g[[1]]$mu), G)
      gsum <- gsum + colSums(G); n1 <- n1 + length(t1)
    }
  }
  parts <- list(cp = cp_full, U = U, gsum = gsum, n1 = n1)
  cj <- .rc_contrib(X, A, parts, reliability)
  chi <- phi + cj
  se_fixed <- sqrt(.if_vcov(phi)); se <- sqrt(.if_vcov(chi))
  list(method_version = RC_FOLD_VERSION, calibration = calibration, est = mean(tau),
       se = drop(se), se_fixed = drop(se_fixed), phi = phi, c = cj, chi = chi, parts = parts,
       calib_share = var(cj) / var(chi), e_hat = e_hat,
       clip_frac = mean(e_hat < trim[1] | e_hat > trim[2]), n_flag = n_flag)
}

# Estimator at a GIVEN calibration (full-sample map, no IF): used to check dQ
# against central differences of the refitted estimator.
regcal_at <- function(X, A, dY, cp, learner, fold, trim = PS_TRIM) {
  Xh <- .rc_apply(as.matrix(X), A, cp); colnames(Xh) <- colnames(X)
  fold_drdid(Xh, A, dY, learner, fold, trim)$est
}
