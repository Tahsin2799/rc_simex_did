# Lasso on a natural cubic spline dictionary (fourth amendment, amendment_2026-10-11b/AMENDMENT_4.md).
#   source("R_scripts/30_core.R"); source("R_scripts/34_simex_fold.R"); source("R_scripts/41_spline_learner.R")
# Every non-binary column is expanded into a natural cubic spline basis with SPLINE_DF degrees of
# freedom (knots at training quantiles; linear beyond the boundary knots, as splines::ns does);
# binary columns enter as they are. The penalty and path are those of the "glmnet" learner. The
# fitted function is smooth, so central differences of its predictions give the gradient needed
# by regression calibration's calibration term to O(h^2), not exactly as for polynomial dictionaries.

SPLINE_DF <- 5L

.ns_transform <- function(X, df = SPLINE_DF) {
  is_bin <- apply(X, 2, \(v) all(v %in% c(0, 1)))
  bases <- lapply(seq_len(ncol(X)), \(j) if (is_bin[j]) NULL else splines::ns(X[, j], df = df))
  \(newX) do.call(cbind, lapply(seq_len(ncol(X)), \(j)
    if (is_bin[j]) newX[, j, drop = FALSE] else predict(bases[[j]], newX[, j])))
}

.learner_glmnet_ns <- function(X, y, type) {
  fam <- if (type == "ps") "binomial" else "gaussian"
  tf <- .ns_transform(X)
  fit <- cv.glmnet(tf(X), y, family = fam, nfolds = 5, alpha = 1,
                   nlambda = LASSO_PATH$nlambda, lambda.min.ratio = LASSO_PATH$lambda.min.ratio)
  \(newX) as.numeric(predict(fit, newx = tf(newX), s = "lambda.min", type = "response"))
}
FOLD_LEARNERS$glmnet_ns <- .learner_glmnet_ns
