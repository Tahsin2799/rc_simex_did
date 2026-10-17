# Third post-production amendment: covariate shape, and a calibration-aware variance for regression calibration

**Dated 11 October 2026. Status: FROZEN on 11 October 2026, before any study replicate ran. Claude froze it under the author's instruction of 11 October 2026 ("go on and do 1 and 2"). The author reviews it with the manuscript.**

This amendment was written after all earlier results were seen, including the second amendment's comparator ([AMENDMENT_2.md](../amendment_2026-10-08/AMENDMENT_2.md)). In every classical-error design run so far, group-wise regression calibration with SIMEX's information came within Monte Carlo error of the latent-data fit and beat SIMEX on mean absolute error. Its fixed-calibration intervals, however, under-covered in the closed-form design (0.89-0.93). This amendment answers two questions the author selected on 11 October 2026 (routes 1 and 2 of the 8 October readiness assessment in `PAPER_LOG.md`):

1. Does the ranking change in a design where linear calibration should fail?
2. Can regression calibration be given a calibration-aware variance by the same machinery as SIMEX, and does it fix the under-coverage?

It is a local, dated record, not an external preregistration. It changes no earlier design or result.

## 1. Method: fold-calibrated regression calibration (`regcal_fold_1.0`, [39_regcal_fold.R](../R_scripts/39_regcal_fold.R))

**Calibration.** Write Z = (W*, V), with error-prone columns W*. S_p is diag(s) on those coordinates and 0 elsewhere, with s = (1 - r) Var(W*) and r the known reliability. Sigma_a = Cov(Z | A = a) and mu_a = E[Z | A = a]. Group-wise regression calibration replaces Z by

  Xhat = Z - S_p Sigma_a^{-1} (Z - mu_a)  for units with A = a.

This is the second amendment's map (`38_regcal.R`), rewritten. Columns constant within a group are dropped from the predictor, and a latent variance that is not positive is floored at 1% of the observed variance and flagged, both as before. In the fold version, eta = (mu_0, mu_1, Sigma_0, Sigma_1, s) is estimated on the training rows of each fold. Both the training and held-out rows of that fold are then mapped with it, and the cross-fitted DR-DiD with fold-normalised weights follows (`fold_drdid` logic).

**Influence function.** The ATT functional at calibration eta is Q(eta) = E[A{D - m_eta(Xhat_1)}]/pi, where m_eta is the control outcome regression on the calibrated covariates. At the true nuisances the propensity drops out. Xhat is an invertible affine map of Z within each group, so m_eta(x) = m*_0(h_0^{-1}(x)), and the directional derivative is

  dQ[delta] = -E[A G' {dh_1(Z)[delta] - dh_0(z')[delta]}] / pi.

Here G is the gradient of m_eta at the treated unit's calibrated point, and z' = h_0^{-1}(Xhat_1). The perturbation is dh_a(z)[delta] = -dR_a (z - mu_a) + R_a dmu_a, with R_a = S_p Sigma_a^{-1}. The calibration contribution of unit j is c_j = dQ[IF_eta(O_j)], with the usual influence functions of the group means, group covariances and s. The combined influence value is chi_j = phi_j + c_j, and the SE is sd(chi)/sqrt(n).

G is the gradient of the fitted control regression, taken by central differences of its predictions. That is exact for the parametric learner and for the lasso dictionary (raw columns plus squares). Tree learners are not covered.

**Expansion.** The fold-specific calibration makes the argument of Theorem 1 (fold-calibrated SIMEX) carry over: condition on the training rows, expand at the realised calibration, then average the calibration error over folds. The conditions are the same: nuisances consistent with product rate o_p(n^{-1/2}) at nearby calibrations, differentiability of Q in eta, and mean-square continuity of the score in eta. The full statement and proof go into the manuscript's Appendix C. **Stated before the runs:** in the log-normal cells of Section 2, no polynomial lasso is consistent for the control regression, so neither this SE nor SIMEX's calibration-aware SE has a guarantee there. Their coverage in those cells is reported as an empirical check, not as a test of the theory.

**Checks passed before freezing** ([rc_fold_checks.R](../runs/rc_fold_checks.R) -> `checks/rc_fold_checks.csv`, `rc_fold_dQ_numeric.csv`, `rc_fold_mc.csv`):

- **The map:** it equals `38_regcal.R` to 1e-15, and the preimage inverts the control map to 1e-15.
- **Vectorised contributions:** they equal the generic directional derivative at influence-function directions, to 1e-15.
- **Central differences:** the analytic dQ matches central differences of the refitted estimator's outcome-regression part. The maximum relative error over 12 random directions in three designs, at n = 200,000 with the parametric learner, is 3.3e-4. The residual is the difference between the fold fits used for G and the single all-control fit used numerically. The tolerance of 1e-3 was set after the first run showed this source; the first tolerance, 1e-4, failed for that reason only.
- **Full DR-DiD derivative:** it differs from the analytic dQ by the propensity-weight path. That path has mean zero when the outcome model is correct, and in the Gaussian design it was at most 0.021 across four directions. In the log-normal designs the parametric learner is misspecified and the difference is large, which is the limitation stated above.
- **Large-sample bias:** the closed-form design at n = 200,000 gives bias -0.003 (SE 0.009).
- **Monte Carlo:** the closed-form design, n = 2,000, parametric learner, 400 replicates, fold calibration. Calibration-aware SE/SD 0.987, coverage 0.955; fixed-calibration SE/SD 0.889, coverage 0.908.

## 2. Study `shape`: covariate shape x trend

**Design** ([40_shape_design.R](../R_scripts/40_shape_design.R)). The closed-form design of the first amendment with its first two moments kept: A ~ Bern(0.4); W_c = A dmu_c + sqrt(v_c - s0_c) L_c, dmu = (0.6, 0.4), v = 1; W*_c = W_c + N(0, s0_c), s0 = (1 - r0) Var(W*), r0 = 0.6; V ~ N(0, 1); dY = beta'W + kappa_q (W_1^2 - E W_1^2) + 0.5 V + 0.5 A + N(0, 1), beta = (1, 0.8). The effect is constant, so the realised ATT is 0.5 in every replicate.

| Factor | Levels |
|---|---|
| shape | `gauss`: L_c ~ N(0, 1). `lnorm`: L_c a standardised log-normal, (exp(omega G) - exp(omega^2/2)) / sd, omega = 0.8 (skewness 3.7) |
| trend | `linear`: kappa_q = 0. `quad`: kappa_q = 0.5 (in the Gaussian case its variance is about one fifth of the linear part's) |
| n | 1,000 and 5,000 |

Means, variances and covariances of (W, W*, V, A) are the same for both shapes, so regression calibration's population map is the same in all four cells. Expectations stated before the runs, derived from the design and not from data:
- In `gauss_linear`, regression calibration is exact.
- In `gauss_quad` it is also exact, because the within-group posterior variance of W is constant and equal across groups, so the quadratic term's calibration error cancels in the contrast.
- In the `lnorm` cells the conditional mean of W given Z is nonlinear and heteroskedastic, so regression calibration is approximate.
- SIMEX's quadratic extrapolation error is unknown in every cell except `gauss_linear`, where q_2(-1) - tau = 0.16.

These are expectations, not hypotheses with a decision rule.

**Methods, on the same data and folds.** The learner is lasso ("glmnet", raw columns plus squares) in all eight cells, plus the parametric learner ("glm") in `gauss_linear` at both n. Learner seeds are set per method from (cell, rep, learner, method).
- `naive`: DR-DiD on Z.
- `latent`: DR-DiD on (W, V).
- `rc_group`: full-sample calibration, as the second amendment. Both SEs are saved; the fixed one is primary, as before.
- `rc_fold`: fold calibration, calibration-aware SE (primary) and fixed SE.
- `simex_fold`: `simex_fold_3.0`, fold calibration, grid (0, 0.5, 1, 2), B = 20, known reliability. Point estimate at t = -1, calibration-aware SE, intervals at M = 0 and M = 0.05 (the grid fixed in the first amendment).

## 3. Study `rc_dgp`: fold calibration on the original design

These are the first amendment's three production-design cells: classical error, nonlinear trends, reliability 0.6 at N = 1,000 and 5,000, and reliability 0.75 at N = 1,000. Lasso learner. The methods are `naive`, `latent`, `rc_group` and `rc_fold`. The production data seeds, the second amendment's folds and its learner seeds are reused, so `naive`, `latent` and `rc_group` reproduce the second amendment's rows exactly; the pilot confirmed this to 1e-15. The question is how large the calibration term is where linear calibration was already near exact.

## 4. Replicates, metrics, reporting

There are 500 replicates per cell: 10 `shape` cells and 3 `rc_dgp` cells, 6,500 tasks.

**Metrics:**
- bias against the realised ATT (MCSE);
- empirical SD;
- SE/SD for each saved SE;
- coverage with Wilson intervals, for each SE, and for SIMEX also at M = 0.05;
- width and RMSE;
- the calibration term's share of the variance, for `rc_fold` and `simex_fold`;
- median and maximum clipping share, and flag counts.

**Paired contrasts** for the estimate, the absolute error and coverage, each with MCSE: `rc_fold` minus `simex_fold`, minus `naive`, minus `latent` and minus `rc_group`; and `simex_fold` minus `latent`.

There is no decision rule. Every cell is reported, whichever correction it favours.

**Compute.** A two-replicate pilot per cell (26 tasks, written to `pilot/`, not part of the study) was used only for timing, errors and the reproducibility check; its estimates were not summarised. Lasso tasks took about 38 s at n = 1,000 and 118 s at n = 5,000 (SIMEX is 95% of that). The total is about 90 CPU-hours, about 9 hours on 10 workers.

## 5. Not included

Not included:
- boosting nuisances;
- estimated reliability;
- other error distributions or non-classical error;
- GSIMEX (Chen and Zhang, 2026) as a further comparator;
- the full-DR propensity-path term in the regression-calibration variance (zero at the truth when the outcome model is correct; omitted by design).

## 6. Deviations log

| Date | Deviation | Reason |
|---|---|---|
| 2026-10-11 02:20-02:40 | `runs/amend3_aggregate.R` was written while the study ran and tested only on the pilot output, with every metric suppressed. It computes the metrics and contrasts listed in Section 4 | Aggregation code only |
| 2026-10-11 11:07 | The study finished: 6,500 tasks, 0 errors, 8.88 h on 10 workers, no restarts | Record |
| 2026-10-11 afternoon | After the results were seen, `runs/amend3_rc_variance_diag.R` split regression calibration's calibration term into its group-mean, group-covariance and error-variance parts (closed-form design, n = 200,000) and computed tail statistics of the `rc_fold` estimates (`checks/rc_var_components.csv`, `rc_var_tails.csv`). These are explanatory diagnostics reported as such in Section 5.5 | Post hoc, labelled in the manuscript |
