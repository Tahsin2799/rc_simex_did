# Population targets with nonparametric nuisances, by numerical integration (no estimator is fitted).
#   source("R_scripts/44_population_np.R")
# Used by runs/population_map.R (the paper's population map). Verified in October 2026 against the
# closed-form Gaussian results (regression calibration's bias b_q (Om_1 - Om_0)) and by grid convergence
# (four decimals when the grids are doubled and tripled); see references/contribution_candidates.md.
#
# Design, per error-prone coordinate c (independent given A): W_c | A = a ~ a dmu_c + sw_{a,c} L, with L a
# standardised log-normal of log-scale SD sl (sl = 0: standard Gaussian); W*_c = W_c + N(0, e); the untreated
# trend is f1(W1) + f2(W2), f1(w) = b1 w + bq w^2, f2(w) = b2 w (V and the noise do not affect the targets).
# With nonparametric nuisances the DR-DiD target minus the ATT is a sum over coordinates of
#   naive at error variance e:   E_1 f(W) - E_{y ~ W*_e | A=1} E_0[f(W) | W*_e = y]
#   calibration h_a:             E_1 f(W) - E_{y ~ W*_e | A=1} E_0[f(W) | W*_e = h_0^{-1}(h_1(y))]
# with h_a the within-group best linear predictor (regression calibration) or the exact posterior mean
# E_a[W | W*_e = y] (posterior-mean calibration). SIMEX extrapolates the naive target at e = s0 (1 + lambda),
# lambda in {0, 0.5, 1, 2}, quadratically to lambda = -1; SIMEX on a calibration does the same with the
# calibrated target, the calibration recomputed at each error variance.

NP <- local({
  gq <- local({ g <- seq(-8.5, 8.5, length.out = 2001); w <- dnorm(g); list(g = g, w = w / sum(w)) })
  gh <- local({ x <- seq(-7, 7, length.out = 161); w <- dnorm(x); list(x = x, w = w / sum(w)) })
  list(gq = gq, gh = gh)
})
np_stdL <- function(sl) if (sl == 0) NP$gq$g else (exp(sl * NP$gq$g) - exp(sl^2 / 2)) / sqrt((exp(sl^2) - 1) * exp(sl^2))
np_lat <- function(a, P) a * P$dmu + (if (a == 1) P$sw1 else P$sw0) * np_stdL(P$sl)
np_post <- function(y, f, a, P, e) {                     # E_a[f(W) | W*_e = y]
  w <- np_lat(a, P); K <- outer(y, w, \(yy, ww) dnorm(yy, ww, sqrt(e)))
  v <- drop(K %*% (f(w) * NP$gq$w)) / drop(K %*% NP$gq$w)
  ok <- is.finite(v)                                     # 0/0 where every kernel weight underflows (far tails)
  if (!all(ok)) v[!ok] <- approx(y[ok], v[ok], xout = y[!ok], rule = 2)$y
  v
}
# y grid: dense up to the 99.99th latent percentile, log-spaced through the far tail of skewed laws.
np_ygrid <- function(P, e) {
  w <- c(np_lat(0, P), np_lat(1, P)); ww <- c(NP$gq$w, NP$gq$w); o <- order(w)
  q_hi <- w[o][which(cumsum(ww[o]) / sum(ww) >= 0.9999)[1]]
  lo <- min(w) - 6 * sqrt(e); mid <- q_hi + 6 * sqrt(e); hi <- max(w) + 6 * sqrt(e)
  c(seq(lo, mid, length.out = 1200), if (hi > mid) mid - 1 + exp(seq(0, log(hi - mid + 1), length.out = 301))[-1])
}
np_treated_avg <- function(gfun, P, e) {                 # E over y ~ W*_e | A = 1 of gfun(y)
  y <- outer(np_lat(1, P), sqrt(e) * NP$gh$x, "+")
  sum(NP$gq$w * drop(matrix(gfun(as.vector(y)), nrow(y)) %*% NP$gh$w))
}
np_E1f <- function(f, P) sum(NP$gq$w * f(np_lat(1, P)))
np_blp <- function(a, P, e) { m <- a * P$dmu; v <- (if (a == 1) P$sw1 else P$sw0)^2
  k <- v / (v + e); list(f = \(y) m + k * (y - m), inv = \(x) m + (x - m) / k) }

# Naive, regression-calibration and posterior-mean-calibration targets (minus the ATT) for one coordinate.
np_terms <- function(f, P, e) {
  yg <- np_ygrid(P, e); pm0f <- np_post(yg, f, 0, P, e)
  m0 <- \(y) approx(yg, pm0f, xout = y, rule = 2)$y
  h1 <- np_blp(1, P, e); h0 <- np_blp(0, P, e)
  g0 <- cummax(np_post(yg, identity, 0, P, e)); g1 <- cummax(np_post(yg, identity, 1, P, e))
  c(naive = np_E1f(f, P) - np_treated_avg(m0, P, e),
    rc = np_E1f(f, P) - np_treated_avg(\(y) m0(h0$inv(h1$f(y))), P, e),
    pm = np_E1f(f, P) - np_treated_avg(\(y) m0(approx(g0, yg, xout = approx(yg, g1, xout = y, rule = 2)$y,
                                                         rule = 2, ties = mean)$y), P, e))
}
# Fisher non-Gaussianity of W*_e within group a: J = Var(W*) I(W*) - 1 (Brown's identity: the share of
# Var(W*) that the best linear predictor misses relative to the posterior mean is (e / Var(W*))^2 J).
np_J <- function(a, P, e) {
  yg <- np_ygrid(P, e); w <- np_lat(a, P)
  p <- drop(outer(yg, w, \(yy, ww) dnorm(yy, ww, sqrt(e))) %*% NP$gq$w)
  dy <- diff(yg); pm <- (p[-1] + p[-length(p)]) / 2; sc <- diff(log(p)) / dy
  v <- sum(NP$gq$w * (w - sum(NP$gq$w * w))^2) + e
  v * sum(sc^2 * pm * dy) / sum(pm * dy) - 1
}

# A design of the shape family (moments as in R_scripts/40_shape_design.R), group-variance ratio for W1.
np_design <- function(sl, r, bq, ratio = 1, b1 = 1, b2 = 0.8, dmu = c(0.6, 0.4), pi = 0.4) {
  varWs <- 1 + pi * (1 - pi) * dmu^2; s0 <- (1 - r) * varWs; sw <- sqrt(1 - s0)
  list(P1 = list(dmu = dmu[1], sw0 = sw[1], sw1 = sqrt(ratio) * sw[1], sl = sl),
       P2 = list(dmu = dmu[2], sw0 = sw[2], sw1 = sw[2], sl = sl),
       s0 = s0, f1 = \(w) b1 * w + bq * w^2, f2 = \(w) b2 * w)
}
NP_LAM <- c(0, 0.5, 1, 2)
NP_XW <- local({ Zg <- cbind(1, NP_LAM, NP_LAM^2); drop(t(solve(crossprod(Zg), t(Zg))) %*% c(1, -1, 1)) })
np_targets <- function(sl, r, bq, ratio = 1) {
  D <- np_design(sl, r, bq, ratio)
  C <- vapply(NP_LAM, \(l) np_terms(D$f1, D$P1, D$s0[1] * (1 + l)) + np_terms(D$f2, D$P2, D$s0[2] * (1 + l)), numeric(3))
  J0 <- np_J(0, D$P1, D$s0[1])
  data.frame(sl = sl, r = r, bq = bq, ratio = ratio,
             naive = C["naive", 1], rc = C["rc", 1], pm_rc = C["pm", 1],
             simex = sum(NP_XW * C["naive", ]), pm_simex = sum(NP_XW * C["pm", ]),
             J0_w1 = J0, J1_w1 = np_J(1, D$P1, D$s0[1]),
             gap_share_w1 = D$s0[1]^2 * J0)               # Var(W*_1 | A = 0) = 1 in this family
}
