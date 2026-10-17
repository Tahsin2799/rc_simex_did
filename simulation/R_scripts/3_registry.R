# Estimator registry, design v2.0. Base panel (12) runs in every cell; corrected
#   source("R_scripts/3_registry.R")
# estimators run only in the ME cells they target (see estimator_set()).

source("R_scripts/30_core.R")
source("R_scripts/31_estimators.R")
source("R_scripts/32_simex.R")
source("R_scripts/33_misclass.R")

base_estimators <- function() {
  list(
    did_ols    = did_ols,
    did_ipw    = did_ipw,
    drdid_imp  = drdid_imp,
    did_ps_rf  = did_ps_rf,
    did_ps_xgb = did_ps_xgb,
    dml_lasso  = dml_lasso,
    dml_lasso_int = dml_lasso_int,
    dml_rf     = dml_rf,
    dml_xgb    = dml_xgb,
    cf_did     = cf_did,
    bcf_did    = bcf_did,
    match_full = match_full
  )
}

# Corrected estimators need cell-level ME parameters (the assumed
# reliability, the misclassification rates). `me_args` is the wrapper's
# argument list from the grid; `r_L` is the lower reliability bound used by
# the partial-identification estimator (pre-registered at 0.6).
corrected_estimators <- function(me_design, me_args, r_L = 0.6) {
  out <- list()
  if (me_design %in% c("x_classical", "x_nc")) {
    # Known-reliability SIMEX uses the true reliability in x_classical cells
    # and the SAME nominal value in x_nc cells (deliberate misspecification:
    # the error there is treatment-correlated, not classical).
    r <- if (me_design == "x_classical") me_args$reliability else 0.75
    out$simex_lasso  <- \(d) simex_drdid(d, dml_lasso, reliability = r)
    out$simex_xgb    <- \(d) simex_drdid(d, dml_xgb,   reliability = r)
    out$bounds_lasso <- \(d) bounds_drdid(d, dml_lasso, r_L = r_L)
    out$bounds_xgb   <- \(d) bounds_drdid(d, dml_xgb,   r_L = r_L)
  }
  if (me_design == "treat_misclass") {
    # Orthogonal known-rate correction (paper Theorem 3): tau (exact under
    # heterogeneity, divides by f) and tau_R (fe-weighted, exact under a
    # homogeneous effect). Cross-fit xgboost nuisances m_1, m_0, q.
    # Two learners, as for SIMEX: the 2026-09-12 oracle check showed lasso
    # nuisances reproduce the oracle at N = 20000 while boosted ones carry
    # learner bias (paper Section 3.6 remark).
    a0 <- me_args$alpha0; a1 <- me_args$alpha1
    out$misclass_lasso   <- \(d) misclass_corr_known(d, a0, a1, "H", "glmnet")
    out$misclass_lasso_R <- \(d) misclass_corr_known(d, a0, a1, "R", "glmnet")
    out$misclass_xgb     <- \(d) misclass_corr_known(d, a0, a1, "H", "xgb")
    out$misclass_xgb_R   <- \(d) misclass_corr_known(d, a0, a1, "R", "xgb")
  }
  out
}

estimator_set <- function(me_design, me_args) {
  c(base_estimators(), corrected_estimators(me_design, me_args))
}
