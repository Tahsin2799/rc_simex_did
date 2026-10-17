# Second post-production amendment: follow-up studies

**Dated 8 October 2026. Status: FROZEN on 8 October 2026, before any follow-up replicate ran. Claude froze it under the author's instruction of 8 October 2026 ("10 - run it", "11 - do it", "12 - do it"). The author reviews it with the manuscript.**

This amendment was written after the v2.0 production results and the first amendment's confirmatory results ([AMENDMENT.md](../amendment_2026-10-06/AMENDMENT.md)) were seen. It is a local, dated record in the repository, not an external preregistration. It changes neither the frozen v2.0 design nor the first amendment. It closes three open items of [post_prod_sim_steps.md](../../post_prod_sim_steps.md): the BCF diagnosis (Section 6, P1), the comparator under a common information budget (Section 5) and the tau_R truth (Gate R).

## 1. BCF diagnosis pilot

**Question.** Is the under-coverage of `bcf_did` in the outcome-error cells caused by short MCMC (production used one chain with 250 burn-in iterations), or by the model?

| Item | Specification |
|---|---|
| Cells | `baseline__none__linear__N5000`; `y_diff__medium__linear__N5000` (production posterior SD / empirical SD 0.73); `y_diff__medium__nonlinear__N5000` (0.78). Production seeds, so the data equal production's |
| Replicates | 1-20 |
| Configurations | `prod`: nburn 250, nsim 1000. `long`: nburn 1000, nsim 1000. Four independent chains each, one `bcf()` call per chain with its own `random_seed`. Everything else as production `bcf_did`: logit `pihat` on the design matrix, 100/50 trees, `x_control = x_moderate` |
| Quantity | Treated-average tau per retained draw |
| Diagnostics | Rank-normalised split R-hat (maximum of bulk and folded) and bulk- and tail-ESS (5% and 95% quantiles) of Vehtari et al. (2021), implemented in [37_mcmc_diag.R](../R_scripts/37_mcmc_diag.R) and validated on AR(1) chains before use. Trace plots for replicate 1 of each cell |
| Comparisons | Posterior mean and SD per configuration, from the pooled chains and from chain 1 alone (the production setting); paired `long - prod` differences with MCSE; posterior SD over the empirical SD of the 500 production estimates of the same cell |

**Reading rule, fixed now.** MCMC length is judged not to explain the under-coverage if, under `prod`, rank-normalised R-hat < 1.01 and bulk-ESS > 400 in at least 90% of fits, and `long` changes the mean posterior SD by less than 5% or by less than two paired MCSEs. Otherwise MCMC length is judged to contribute.

**Compute.** 480 chain runs of about 110-180 s, roughly 19 CPU-hours on 6 workers.

## 2. Regression-calibration comparator under the same information

**Question.** How does a standard correction that uses exactly SIMEX's information (a known reliability and the sample variance of each error-prone covariate) compare with SIMEX?

**Method.** Regression calibration (Carroll et al., 2006). Each error-prone covariate is replaced by its best linear predictor given all observed covariates, computed within each treatment group. In that predictor, Cov(W, W\*) is the within-group Var(W\*) minus the error variance (1 - r) Var(W\*), with Var(W\*) taken over the whole sample as in the SIMEX calibration. The calibrated covariates then enter the same cross-fitted DR-DiD (`fold_drdid`: Hajek weights normalised within folds) with the learner SIMEX used. The SE comes from the DR-DiD influence values with the calibration treated as fixed, as production SIMEX did. A group's calibrated error-free variance that is not positive is floored at 1% of its observed variance and flagged.

Pooled calibration (one predictor for both groups) is an affine map of the observed covariates. For linear-index learners it spans the same space as the naive fit, so it is run only in the closed-form design, as a check. Each task also fits naive DR-DiD and DR-DiD on the latent covariates, with the same folds.

| Study | Cells | Learner | Pairs with |
|---|---|---|---|
| `rc_gauss` | closed form, reliability 0.6 known: n = 1,000 and 5,000. Regression calibration is exact here (Gaussian covariates, linear trend) | glm (both n); lasso (n = 1,000) | first amendment's `gauss_known` SIMEX rows (same data seeds) |
| `rc_dgp` | production design, `x_classical` r in {0.60, 0.75, 0.90} x {linear, nonlinear} x N in {1,000, 5,000}; production seeds | lasso (as `simex_fold`) | production `simex_lasso` and `dml_lasso`; revised SIMEX in the three first-amendment cells |
| `rc_xnc` | production design, `x_nc__large`, linear and nonlinear, N = 5,000; reliability 0.75, as production SIMEX used. The error model is wrong for both corrections | lasso | production `simex_lasso` |

There are 500 replicates per cell. The metrics are bias against the realised ATT (MCSE), empirical SD, SE/SD, coverage with Wilson intervals, width and RMSE. Paired differences (regression calibration minus SIMEX, minus naive and minus latent) are reported for the estimate, the absolute error and coverage, each with its MCSE, along with flag counts. There is no decision rule; every cell is reported.

**Compute.** Three or four fits per task, roughly 10 CPU-hours on 5 workers.

## 3. tau_R truth for the R-form misclassification corrections

This part only re-scores existing estimates; no estimator is refitted. `misclass_lasso_R` and `misclass_xgb_R` target tau_R = E[theta f e] / E[f e], with f = kappa e(1 - e) / {q(1 - q)} and q = alpha0 + kappa e. Production scored them against the realised ATT. For each of the 12 `treat_misclass` cells and each of the 500 replicates, the latent data are regenerated with the production seed, and the realised tau_R = sum(A theta f) / sum(A f) is computed with the latent A and the true e.

theta depends only on c11-c15, which are drawn independently of the covariates in e(X). The population tau_R therefore equals the population ATT (MC10), which is verified numerically on one large draw per family. The R-form rows are re-scored against the realised tau_R (bias with MCSE, coverage with Wilson intervals). These scores are reported beside the ATT-based scores, together with the per-replicate gap between realised tau_R and the realised ATT.

## 4. Reporting additions to the first amendment (no refits)

These are the runtime median, the 90th percentile and the total task-seconds of the confirmatory study, plus the median share of clipped propensities per cell. They are reported in Appendix E and Section 6.

## 5. Not included

The following are out of scope: BART with heteroskedastic errors and other BCF variants, regression calibration with an estimated reliability, and other comparators (for example Wen and Yan, 2025).

## 6. Deviations log

| Date | Deviation | Reason |
|---|---|---|
| 2026-10-08 | This file was last edited at 02:32. The BCF pilot started about 02:40 and the comparator at 02:44. The file was committed at 03:15 (commit fe9539b), after partial results had been seen while testing the aggregator: about 70 replicates per N = 5,000 comparator cell and 9 BCF fits. The design text above was not changed after the runs started | Timing record; the commit title "before any results" refers to the design, not the commit |
| 2026-10-08 03:09 | `followup_bcf_pilot.R` gained a `--reverse` option, so a second runner can use the cores freed when the comparator finishes. Each task re-checks its file before starting | Computational; tasks, seeds and outputs unchanged |
| 2026-10-08 | `followup_aggregate.R` tolerates configurations with no complete fit yet | Aggregation code only |
| 2026-10-08 07:00-11:41 | The machine shut down unexpectedly at about 07:00 (Windows events 41 and 6008; reboot at 11:11) while the BCF pilot had 396 of 480 chain files. All 396 files were checked readable. The runner was restarted at 11:41 with 11 workers, since the comparator had finished, and completed the remaining 84 chains by 12:02 | Computational; each task's seed depends only on cell, replicate, configuration and chain, and finished files are skipped, so no output differs from an uninterrupted run |
