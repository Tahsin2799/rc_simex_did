# Natural-spline lasso with constant extrapolation (fifth amendment, amendment_2026-10-14/AMENDMENT_5.md).
#   source("R_scripts/30_core.R"); source("R_scripts/34_simex_fold.R"); source("R_scripts/41_spline_learner.R")
#   source("R_scripts/42_spline_clamped.R")
# The "glmnet_ns" learner of 41_spline_learner.R, except that each column of a new design matrix is
# clamped to the range of the training data before the spline basis is evaluated. Inside the training
# range the basis, the fit and every prediction are those of "glmnet_ns"; beyond it the fitted
# function is constant (its value at the boundary) instead of linear. Clamping the training matrix
# to its own range changes nothing, so with the same RNG state the cv.glmnet fit is identical.

.clamp_cols <- function(lo, hi) \(newX) sweep(sweep(newX, 2, lo, pmax), 2, hi, pmin)

.learner_glmnet_nsc <- function(X, y, type) {
  fam <- if (type == "ps") "binomial" else "gaussian"
  tf <- .ns_transform(X)
  clamp <- .clamp_cols(apply(X, 2, min), apply(X, 2, max))
  fit <- cv.glmnet(tf(X), y, family = fam, nfolds = 5, alpha = 1,
                   nlambda = LASSO_PATH$nlambda, lambda.min.ratio = LASSO_PATH$lambda.min.ratio)
  \(newX) as.numeric(predict(fit, newx = tf(clamp(newX)), s = "lambda.min", type = "response"))
}
FOLD_LEARNERS$glmnet_nsc <- .learner_glmnet_nsc
