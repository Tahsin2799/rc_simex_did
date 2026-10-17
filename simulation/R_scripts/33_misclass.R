# Orthogonal known-rate correction for treatment misclassification under
# conditional parallel trends (paper Prop. 4, Theorem 3; PreReg v2.0 s4.3).
#   source("R_scripts/33_misclass.R")
#
# Observed A* with P(A*=1|A=1,X) = 1 - alpha1, P(A*=1|A=0,X) = alpha0,
# misclassification independent of (X, potential outcomes); kappa = 1 - a0 - a1.
# Observed-data nuisances:  q(X) = P(A*=1|X),  m_a(X) = E[dY | A*=a, X],
# c = m_1 - m_0.  Then  e = (q - a0)/kappa,  f = kappa e(1-e)/(q(1-q)),
# c = theta f, and
#   tau   = E[c K_H(q)] / E[L_H(q)],  K_H = e/f = q(1-q)/(1-a1-q),  L_H = e
#   tau_R = E[c K_R(q)] / E[L_R(q)],  K_R = e,                     L_R = f e
# (tau exact under heterogeneity; tau_R the fe-weighted ATET, exact under a
# homogeneous effect and less variable because it does not divide by f).
#
# Neyman-orthogonal scores (paper eq. 20-21):
#   phi_K = c K(q) + K(q)[A*(dY - m1)/q - (1-A*)(dY - m0)/(1-q)] + c K'(q)(A* - q)
#   phi_L = L(q) + L'(q)(A* - q)
# Estimator: Psi_hat = mean(phi_K) / mean(phi_L) with K-fold cross-fit
# nuisances; IF psi = (phi_K - Psi phi_L)/mean(phi_L); se = sd(psi)/sqrt(N).
# q_hat is clipped to [a0 + kappa*delta, 1 - a1 - kappa*delta], delta = PS_TRIM[1].
#
# Replaces the 2026-09-12 plug-in (linear OR by observed group, logit e*,
# bootstrap SE), which was misspecified because m_a = mu + theta p_a(X) is
# non-linear in X even when mu is linear.

.mc_funs <- function(alpha0, alpha1) {
  kappa <- 1 - alpha0 - alpha1
  e  <- \(q) (q - alpha0) / kappa
  K_H  <- \(q) q * (1 - q) / (1 - alpha1 - q)
  K_Hp <- \(q) ((1 - 2 * q) * (1 - alpha1 - q) + q * (1 - q)) / (1 - alpha1 - q)^2
  L_H  <- e;  L_Hp <- \(q) rep(1 / kappa, length(q))
  K_R  <- e;  K_Rp <- L_Hp
  G    <- \(q) (q - alpha0)^2 * (1 - alpha1 - q) / (kappa^2 * q * (1 - q))
  Gp   <- function(q) {
    n  <- (q - alpha0)^2 * (1 - alpha1 - q); dd <- kappa^2 * q * (1 - q)
    np <- 2 * (q - alpha0) * (1 - alpha1 - q) - (q - alpha0)^2
    dp <- kappa^2 * (1 - 2 * q)
    (np * dd - n * dp) / dd^2
  }
  list(kappa = kappa,
       H = list(K = K_H, Kp = K_Hp, L = L_H, Lp = L_Hp),
       R = list(K = K_R, Kp = K_Rp, L = G,   Lp = Gp))
}

# Cross-fit the three nuisances once; folds stratified on A*.
.mc_nuisances <- function(X, Astar, dY, learner = "xgb", K = 5L, q_lo, q_hi) {
  L <- LEARNERS[[learner]]; N <- length(Astar)
  fold <- integer(N)
  fold[Astar == 1] <- sample(rep(seq_len(K), length.out = sum(Astar == 1)))
  fold[Astar == 0] <- sample(rep(seq_len(K), length.out = sum(Astar == 0)))
  q_hat <- m1_hat <- m0_hat <- numeric(N)
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k)
    q_hat[te]  <- L(X[tr, , drop = FALSE], Astar[tr], "ps")(X[te, , drop = FALSE])
    t1 <- tr[Astar[tr] == 1]; t0 <- tr[Astar[tr] == 0]
    m1_hat[te] <- L(X[t1, , drop = FALSE], dY[t1], "or")(X[te, , drop = FALSE])
    m0_hat[te] <- L(X[t0, , drop = FALSE], dY[t0], "or")(X[te, , drop = FALSE])
  }
  list(q = pmin(pmax(q_hat, q_lo), q_hi), m1 = m1_hat, m0 = m0_hat)
}

.mc_estimate <- function(Astar, dY, nu, fk) {
  q <- nu$q; c_hat <- nu$m1 - nu$m0; N <- length(Astar)
  phi_K <- c_hat * fk$K(q) +
           fk$K(q) * (Astar * (dY - nu$m1) / q - (1 - Astar) * (dY - nu$m0) / (1 - q)) +
           c_hat * fk$Kp(q) * (Astar - q)
  phi_L <- fk$L(q) + fk$Lp(q) * (Astar - q)
  est <- mean(phi_K) / mean(phi_L)
  psi <- (phi_K - est * phi_L) / mean(phi_L)
  .out(est, sd(psi) / sqrt(N))
}

# target = "H" (tau, exact under heterogeneity) or "R" (tau_R, fe-weighted).
misclass_corr_known <- function(d, alpha0, alpha1, target = c("H", "R"),
                                learner = "xgb", K = 5L) {
  target <- rlang::arg_match(target)
  fn <- .mc_funs(alpha0, alpha1)
  q_lo <- alpha0 + fn$kappa * PS_TRIM[1]; q_hi <- 1 - alpha1 - fn$kappa * PS_TRIM[1]
  nu <- .mc_nuisances(design_matrix(d), d$A, d$dY, learner, K, q_lo, q_hi)
  .mc_estimate(d$A, d$dY, nu, fn[[target]])
}
