# Shared estimator core, design v2.0.
#   source("R_scripts/30_core.R")
#
#   .out()            standard return vector
#   pin_threads()     one thread per process (called by the batch runner)
#   crossfit_drdid()  K-fold cross-fit Sant'Anna-Zhao (2020) panel DR-DiD
#                     score with pluggable nuisance learners; returns the
#                     estimate, the influence-function SE and the IF itself.

suppressPackageStartupMessages({
  library(glmnet)
  library(ranger)
  library(xgboost)
})

PS_TRIM <- c(0.02, 0.98)
Z975    <- qnorm(0.975)

# Every estimator returns this shape. lo/hi are identified-set endpoints for
# partial-identification estimators (NA otherwise); ci_lo/ci_hi is the
# interval whose coverage is scored (Wald by default).
.out <- function(est, se, ci_lo = est - Z975 * se, ci_hi = est + Z975 * se,
                 lo = NA_real_, hi = NA_real_) {
  c(est = est, se = se, ci_lo = ci_lo, ci_hi = ci_hi, lo = lo, hi = hi)
}

pin_threads <- function() {
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  options(ranger.num.threads = 1L)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1L); RhpcBLASctl::omp_set_num_threads(1L)
  }
  invisible(NULL)
}

# ---- nuisance learners: each returns a function(newX) -> predictions -------
# type = "ps" (binary, P(A=1|X)) or "or" (continuous mean).

# Lasso dictionary. The linear learner is fed a polynomial expansion of the
# design matrix, the standard "flexible" specification for linear nuisance
# learners in DML (Chernozhukov et al. 2018 Sec. 6 and the DoubleML 401(k)
# example: raw columns plus degree-2 polynomials of the continuous
# covariates; Zimmert 2020 and Knaus-Lechner-Strittmatter 2021 add all
# two-way interactions). Tree learners get the raw columns; they form
# interactions internally.
#   "raw"      27 raw columns. The v2.0 draft learner; with p = 27 << N this
#              is a cross-fit logit / OLS, not a lasso. Kept for the
#              2026-09-12 dictionary comparison (PreReg v2.0 s9).
#   "sq"       + squares of the 12 non-binary columns = 39. REGISTERED for
#              dml_lasso and for every lasso-based correction (PreReg s4.2, V13).
#   "sq_int"   + all 136 pairwise products of the 17 non-dummy columns = 175.
#              REGISTERED for dml_lasso_int, base panel only: at 4-6x the
#              cost of "sq" it is unaffordable inside the 61-refit SIMEX
#              wrapper (PreReg s10).
#   "quad"     + all pairwise products including dummies, ~370 after dropping
#              within-factor products (always zero). Not used: 12-16 s a fit.
# Non-binary columns are standardised on the training fold before expansion:
# products of raw-scale columns are badly conditioned and slow the
# coordinate descent by an order of magnitude. Columns constant on the
# training fold are dropped and the same map is applied at prediction. The
# penalty path is 50 values down to 0.01 lambda_max; the near-unpenalised
# tail of the glmnet default path costs most of the runtime and is never
# selected.
LASSO_DICT <- "sq"
LASSO_PATH <- list(nlambda = 50L, lambda.min.ratio = 0.01)

.dict_cols <- function(X, dict) {
  nm <- colnames(X); p <- ncol(X)
  is_bin <- apply(X, 2, \(v) all(v %in% c(0, 1)))
  is_dum <- grepl("_", nm, fixed = TRUE)          # dummy-coded factor levels
  out <- list(X)
  if (dict %in% c("sq", "sq_int", "quad")) {
    sq <- X[, !is_bin, drop = FALSE]^2
    colnames(sq) <- paste0(nm[!is_bin], "^2"); out <- c(out, list(sq))
  }
  if (dict %in% c("sq_int", "quad")) {
    idx <- if (dict == "quad") seq_len(p) else which(!is_dum)
    ij <- utils::combn(idx, 2)
    pr <- X[, ij[1, ], drop = FALSE] * X[, ij[2, ], drop = FALSE]
    colnames(pr) <- paste(nm[ij[1, ]], nm[ij[2, ]], sep = ":"); out <- c(out, list(pr))
  }
  do.call(cbind, out)
}

# Returns a transformer fitted on X: standardises non-binary columns, expands,
# drops constant columns; the returned function applies the same map to newX.
.dict_transform <- function(X, dict) {
  if (dict == "raw") return(\(newX) newX)
  is_bin <- apply(X, 2, \(v) all(v %in% c(0, 1)))
  mu <- colMeans(X); s <- apply(X, 2, sd); mu[is_bin] <- 0; s[is_bin] <- 1; s[s == 0] <- 1
  f <- \(newX) .dict_cols(scale(newX, center = mu, scale = s), dict)
  Z <- f(X)
  keep <- which(apply(Z, 2, \(v) any(v != v[1])))
  \(newX) f(newX)[, keep, drop = FALSE]
}

.learner_glmnet <- function(X, y, type, dict = LASSO_DICT) {
  fam <- if (type == "ps") "binomial" else "gaussian"
  tf  <- .dict_transform(X, dict)
  fit <- cv.glmnet(tf(X), y, family = fam, nfolds = 5, alpha = 1,
                   nlambda = LASSO_PATH$nlambda, lambda.min.ratio = LASSO_PATH$lambda.min.ratio)
  \(newX) as.numeric(predict(fit, newx = tf(newX), s = "lambda.min", type = "response"))
}

.learner_ranger <- function(X, y, type) {
  if (type == "ps") {
    fit <- ranger(x = X, y = factor(y, levels = c(0, 1)), num.trees = 500,
                  min.node.size = 5, probability = TRUE, num.threads = 1)
    \(newX) predict(fit, data = newX, num.threads = 1)$predictions[, "1"]
  } else {
    fit <- ranger(x = X, y = y, num.trees = 500, min.node.size = 5, num.threads = 1)
    \(newX) predict(fit, data = newX, num.threads = 1)$predictions
  }
}

# xgboost (hist). nrounds chosen once per call by early stopping on a 20 %
# holdout of the training fold, then refit on the full training fold.
.learner_xgb <- function(X, y, type, max_rounds = 500L) {
  obj <- if (type == "ps") "binary:logistic" else "reg:squarederror"
  params <- list(objective = obj, eta = 0.05, max_depth = 3, subsample = 0.7,
                 colsample_bytree = 0.8, min_child_weight = 5,
                 tree_method = "hist", nthread = 1)
  n <- nrow(X); hold <- sample.int(n, floor(0.2 * n))
  dtr <- xgb.DMatrix(X[-hold, , drop = FALSE], label = y[-hold])
  dho <- xgb.DMatrix(X[hold,  , drop = FALSE], label = y[hold])
  es  <- xgb.train(params = params, data = dtr, nrounds = max_rounds,
                   evals = list(ho = dho), early_stopping_rounds = 20L,
                   verbose = 0)
  # No best_iteration recorded: older xgboost raises, newer returns NULL or a
  # string that will not coerce. Each case falls back to 100 rounds.
  best <- tryCatch(xgb.attr(es, "best_iteration"), error = \(e) NULL)
  best <- withCallingHandlers(as.integer(best), warning = \(w) {
    if (grepl("coercion", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")
  })
  if (length(best) != 1L || is.na(best) || best < 1L) best <- 100L
  fit <- xgb.train(params = params, data = xgb.DMatrix(X, label = y),
                   nrounds = best, verbose = 0)
  \(newX) as.numeric(predict(fit, xgb.DMatrix(newX)))
}

LEARNERS <- list(glmnet = .learner_glmnet, ranger = .learner_ranger,
                 xgb = .learner_xgb,
                 glmnet_int = \(X, y, type) .learner_glmnet(X, y, type, dict = "sq_int"))

# ---- cross-fit DR-DiD -------------------------------------------------------
# Sant'Anna & Zhao (2020) panel score with Hajek-normalised weights:
#   tau = mean(w1 r) - mean(w0 r),  r = dY - m0(X),
#   w1 = A / mean(A),  w0 = [(1-A) e/(1-e)] / mean[(1-A) e/(1-e)].
# Influence function (nuisances treated as known, valid under cross-fitting
# and o(N^-1/4) nuisance rates):
#   psi = w1 (r - tau1) - w0 (r - tau0).
crossfit_drdid <- function(X, A, dY, learner = "xgb", K = 5L,
                           trim = PS_TRIM, return_if = FALSE) {
  L <- LEARNERS[[learner]]
  N <- length(A)
  fold <- integer(N)
  fold[A == 1] <- sample(rep(seq_len(K), length.out = sum(A == 1)))
  fold[A == 0] <- sample(rep(seq_len(K), length.out = sum(A == 0)))
  e_hat <- numeric(N); m0_hat <- numeric(N)
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k)
    e_hat[te]  <- L(X[tr, , drop = FALSE], A[tr], "ps")(X[te, , drop = FALSE])
    ctrl <- tr[A[tr] == 0]
    m0_hat[te] <- L(X[ctrl, , drop = FALSE], dY[ctrl], "or")(X[te, , drop = FALSE])
  }
  drdid_score(A, dY, e_hat, m0_hat, trim = trim, return_if = return_if)
}

drdid_score <- function(A, dY, e_hat, m0_hat, trim = PS_TRIM, return_if = FALSE) {
  e_hat <- pmin(pmax(e_hat, trim[1]), trim[2])
  N  <- length(A)
  r  <- dY - m0_hat
  w1 <- A / mean(A)
  w0 <- (1 - A) * e_hat / (1 - e_hat); w0 <- w0 / mean(w0)
  tau1 <- mean(w1 * r); tau0 <- mean(w0 * r)
  est <- tau1 - tau0
  psi <- w1 * (r - tau1) - w0 * (r - tau0)
  se  <- sd(psi) / sqrt(N)
  res <- .out(est, se)
  if (return_if) attr(res, "psi") <- psi
  res
}
