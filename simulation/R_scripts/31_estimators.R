# Base estimator panel, design v2.0. Every function takes the WIDE frame
#   source("R_scripts/31_estimators.R")
# from wide_frame() (id, A, y0, y1, dY, c1..c20) and returns .out().
# Adjustment set is X = c1..c20 (design_matrix: 27 columns, less any column
# with no variation in the sample) for every covariate-adjusting estimator;
# no estimator receives y0 (PAPER_PLAN 3.2 F).

suppressPackageStartupMessages({
  library(DRDID)
  library(grf)
  library(bcf)
  library(MatchIt)
  library(optmatch)
  library(sandwich)
})

# 1. Covariate-blind DiD: dY on A. Correct SE (v1.3's long-form OLS ignored
#    the unit effect and over-covered; PAPER_PLAN 3.2 A).
did_ols <- function(d) {
  fit <- lm(dY ~ A, data = d)
  V   <- sandwich::vcovHC(fit, type = "HC1")
  .out(unname(coef(fit)["A"]), sqrt(V["A", "A"]))
}

# 2. Abadie (2005) IPW-DiD, logit PS, with the PS-estimation term in the IF
#    (DRDID::std_ipw_did_panel).
did_ipw <- function(d) {
  X <- cbind(1, design_matrix(d))
  r <- DRDID::std_ipw_did_panel(y1 = d$y1, y0 = d$y0, D = d$A, covariates = X,
                                boot = FALSE, inffunc = FALSE)
  .out(as.numeric(r$ATT), as.numeric(r$se))
}

# 3. Sant'Anna-Zhao improved DR-DiD (IPT + WLS), parametric benchmark.
#    ps.flag 2 means trust and nlm both failed and DRDID fell back to the
#    plain logit coefficients, so the number is standard DR-DiD and not the
#    improved one. Flag the rep rather than average it in under this label;
#    aggregate.R holds flagged reps out of the primary metrics (PreReg s5).
drdid_imp <- function(d) {
  X <- cbind(1, design_matrix(d))
  r <- DRDID::drdid_imp_panel(y1 = d$y1, y0 = d$y0, D = d$A, covariates = X,
                              boot = FALSE, inffunc = FALSE)
  out <- .out(as.numeric(r$ATT), as.numeric(r$se))
  if (isTRUE(r$ps.flag == 2)) attr(out, "status") <- "ipt_fallback"
  out
}

# 4-5. IPW-DiD with ML propensity only (no outcome model): RF (OOB) and xgb
#      (cross-fitted). Same Hajek score with m0 = 0.
did_ps_rf <- function(d) {
  X <- design_matrix(d)
  rf <- ranger(x = X, y = factor(d$A, levels = c(0, 1)), num.trees = 500,
               min.node.size = 5, probability = TRUE, num.threads = 1)
  drdid_score(d$A, d$dY, rf$predictions[, "1"], 0)
}

did_ps_xgb <- function(d, K = 5L) {
  X <- design_matrix(d); A <- d$A; N <- length(A)
  fold <- integer(N)
  fold[A == 1] <- sample(rep(seq_len(K), length.out = sum(A == 1)))
  fold[A == 0] <- sample(rep(seq_len(K), length.out = sum(A == 0)))
  e_hat <- numeric(N)
  for (k in seq_len(K)) {
    tr <- which(fold != k); te <- which(fold == k)
    e_hat[te] <- .learner_xgb(X[tr, ], A[tr], "ps")(X[te, ])
  }
  drdid_score(A, d$dY, e_hat, 0)
}

# 6-8. Cross-fit DR-DiD with ML nuisances (shared core, SZ panel score).
#      The lasso learner receives the degree-2 dictionary (LASSO_DICT).
#      `return_if = TRUE` attaches the per-unit influence values (attr "psi"),
#      which SIMEX and the bounds need (paper Theorem 1(iii), Theorem 2(i)).
dml_lasso <- function(d, return_if = FALSE)
  crossfit_drdid(design_matrix(d), d$A, d$dY, "glmnet", return_if = return_if)
# 6b. Same score, lasso on the squares-plus-interactions dictionary (base
#     panel only; see LASSO_DICT in 30_core.R).
dml_lasso_int <- function(d, return_if = FALSE)
  crossfit_drdid(design_matrix(d), d$A, d$dY, "glmnet_int", return_if = return_if)
dml_rf    <- function(d, return_if = FALSE)
  crossfit_drdid(design_matrix(d), d$A, d$dY, "ranger", return_if = return_if)
dml_xgb   <- function(d, return_if = FALSE)
  crossfit_drdid(design_matrix(d), d$A, d$dY, "xgb",    return_if = return_if)

# 9. Causal forest on dY, ATT via the grf doubly robust score.
cf_did <- function(d) {
  cf  <- causal_forest(X = design_matrix(d), Y = d$dY, W = d$A, num.threads = 1)
  ate <- average_treatment_effect(cf, target.sample = "treated")
  .out(as.numeric(ate[["estimate"]]), as.numeric(ate[["std.err"]]))
}

# 10. Bayesian causal forest on dY with logit pihat; posterior mean/SD of the
#     treated-average tau. MCMC halved relative to v1.3 (250 / 1000).
bcf_did <- function(d, nburn = 250L, nsim = 1000L) {
  X  <- design_matrix(d)
  ps <- glm(d$A ~ X, family = binomial)
  pihat <- pmin(pmax(as.numeric(fitted(ps)), PS_TRIM[1]), PS_TRIM[2])
  # bcf 2.0 prints MCMC banners regardless of verbose/no_output; swallow them.
  invisible(capture.output(
    fit <- bcf(y = d$dY, z = d$A, x_control = X, x_moderate = X, pihat = pihat,
               n_chains = 1, nburn = nburn, nsim = nsim,
               ntree_control = 100, ntree_moderate = 50,
               update_interval = nsim, no_output = TRUE, verbose = FALSE,
               n_threads = 1),
    type = "output"))
  draws <- rowMeans(fit$tau[, d$A == 1, drop = FALSE])
  .out(mean(draws), sd(draws),
       ci_lo = unname(quantile(draws, 0.025)), ci_hi = unname(quantile(draws, 0.975)))
}

# 11. Optimal full matching on logit PS; weighted dY on A, SE clustered by
#     matched stratum.
match_full <- function(d) {
  X  <- design_matrix(d)
  dd <- data.frame(A = d$A, dY = d$dY, X)
  ps <- as.numeric(fitted(glm(A ~ . - dY, data = dd, family = binomial)))
  options(optmatch_max_problem_size = Inf)
  m  <- matchit(A ~ . - dY, data = dd, distance = ps, method = "full",
                estimand = "ATT")
  md <- match.data(m)
  fit <- lm(dY ~ A, data = md, weights = weights)
  V   <- sandwich::vcovCL(fit, cluster = md$subclass)
  .out(unname(coef(fit)["A"]), sqrt(V["A", "A"]))
}
