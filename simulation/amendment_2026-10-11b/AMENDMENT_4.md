# Fourth post-production amendment: a design where regression calibration should fail

**Dated 11 October 2026. Status: FROZEN on 11 October 2026, before any study replicate ran. Claude froze it under the author's instruction of 11 October 2026 ("for route 1, see if we can find a design where regression calibration should fail"). The author reviews it with the manuscript.**

This amendment was written after all earlier results were seen, including the third amendment ([AMENDMENT_3.md](../amendment_2026-10-11/AMENDMENT_3.md)). There, regression calibration beat SIMEX wherever its calibration model held. With skewed covariates and a quadratic trend at reliability 0.6, both corrections failed. **The design below was chosen deliberately, from population calculations, as one where regression calibration's target is further from the ATT than SIMEX's.** It is not a representative design. Its purpose is to show that such designs exist, and to check whether finite-sample estimates and intervals follow the population calculation. Section 1 reports that calculation in full, including the regions where it favours regression calibration.

## 1. Population calculations made before the design was chosen (no estimator fitted)

**Gaussian covariates with group-specific dispersion** ([amend4_population_map.R](../runs/amend4_population_map.R), closed form, `population_map.csv`). In a Gaussian design where W given (Z, A) has a common variance in both groups, regression calibration is exact for any trend, because the posterior of W is Gaussian with the same variance in both groups (Appendix C.8). Its target departs from the ATT only if Var(W | A) differs between groups and the trend is nonlinear; the departure is bq (Om_1 - Om_0). The grid was: dispersion ratio 1-6, mean shift 0-0.6, bq in {0.5, 1}, reliability in {0.6, 0.75, 0.9} (108 settings). Regression calibration's bias was smaller than SIMEX's extrapolation error in every setting; the ratio was at most 0.32.

**Skewed (log-normal) covariates** ([amend4_population_lnorm.R](../runs/amend4_population_lnorm.R), numerical integration, `population_map_lnorm.csv`), with nonparametric nuisances. This reproduces the Gaussian closed form and the earlier Monte Carlo targets in the Gaussian rows.

| log-scale SD | reliability | trend | naive | regression calibration | SIMEX q_2(-1) |
|---|---|---|---|---|---|
| 0.8 | 0.6 | quadratic | 0.475 | -0.221 | 0.260 |
| 1.0 | 0.75 | linear | 0.407 | 0.170 | 0.224 |
| 1.0 | 0.75 | quadratic | 0.379 | 0.016 | 0.187 |
| 0.8 | 0.9 | linear | 0.188 | 0.089 | 0.050 |
| 1.0 | 0.9 | linear | 0.229 | **0.132** | **0.072** |
| 1.0 | 0.9 | quadratic | 0.192 | 0.053 | 0.046 |

(All entries are biases against the ATT, 0.5.)

SIMEX's population target is closer to the ATT than regression calibration's only at high reliability with skewed covariates. The mechanism: a flexible outcome regression learns the control group's nonlinear posterior mean of W, which does not transfer to the treated group. With a linear outcome model and a linear trend, regression calibration stays exact, because the best linear predictor preserves group means. Whether the lasso learners reproduce the nonparametric targets is therefore part of the question.

**SIMEX's true derivative bound** ([amend4_population_M.R](../runs/amend4_population_M.R), `population_M.csv`). It is 0.009 in the Gaussian reliability-0.9 cell. In the log-normal cells it is above 10, so the allowance M kappa(-1) would be uninformative. The numerical integration also degrades as the error variance tends to zero, so these values are indicative only. Intervals are reported at the pre-declared M in {0, 0.05}, as before.

## 2. Study `shape4`

**Design.** [40_shape_design.R](../R_scripts/40_shape_design.R) with omega = 1.0 (log-normal, skewness 6.2) and reliability 0.9 known, otherwise as in the third amendment: means, variances and covariances as in the closed-form design, trend linear or with 0.5 (W_1^2 - E W_1^2), constant effect, ATT 0.5.

| Cell | Shape | Trend | n | Learners |
|---|---|---|---|---|
| `lnorm_linear_r090_n5000` | log-normal | linear | 5,000 | lasso (raw + squares); lasso (natural splines, 5 df) |
| `lnorm_quad_r090_n5000` | log-normal | quadratic | 5,000 | both |
| `gauss_linear_r090_n5000` | Gaussian | linear | 5,000 | lasso (raw + squares), reference |
| `lnorm_linear_r090_n1000` | log-normal | linear | 1,000 | both |

The natural-spline lasso (`41_spline_learner.R`, `glmnet_ns`) approximates the nonparametric nuisances assumed in Section 1. Its gradient for regression calibration's calibration term is taken by central differences, which is accurate to about 3e-4 here but not exact.

**Methods** on the same data and folds, as in the third amendment: `naive`, `latent`, `rc_group`, `rc_fold` (calibration-aware SE) and `simex_fold` (intervals at M = 0 and 0.05). Seeds come from (cell, rep, learner, method).

**Expectations from Section 1** (nonparametric nuisances):
- `lnorm_linear_r090`: regression calibration bias near 0.13, SIMEX near 0.07.
- `lnorm_quad_r090`: both near 0.05.
- `gauss_linear_r090`: regression calibration exact, SIMEX bias near 0.005.

The lasso with squares may sit closer to regression calibration's exact linear case. There is no decision rule; every cell is reported.

**Metrics and contrasts:** as in the third amendment (`runs/amend3_aggregate.R`, pointed at this root). 500 replicates per cell, 3,500 tasks. A two-replicate pilot per cell (14 tasks, `pilot/`) was used only for timing and errors; its estimates were not summarised. The compute is about 109 CPU-hours, about 11 hours on 10 workers.

## 3. Deviations log

| Date | Deviation | Reason |
|---|---|---|
| 2026-10-12 22:10 to 2026-10-13 07:32 | The controller stopped after 2,738 tasks without writing its closing line. It was restarted on 13 Oct at 06:19 with the same command; the skip-finished logic left 762 tasks, which ran in 1.22 h. All 3,500 task files are readable, the progress log has 3,500 lines for 3,500 distinct tasks, and seeds depend only on (cell, rep, learner, method), so the restart does not change any output | Interruption; no protocol change |
| 2026-10-14 | The study used 335 task-hours, not the 109 planned. In the log-normal n = 5,000 cells with the lasso with squares, `simex_fold` had a heavy runtime tail (median 109-115 s, 95th percentile 3,700-3,800 s). The slow replicates are those with extreme draws: their naive and regression-calibration estimates are also more dispersed. glmnet convergence warnings were not recorded | Compute record |
| 2026-10-14 | **Numerical correction to Section 1.** The log-normal population values for $\omega=1.0$ were recomputed with a non-uniform evaluation grid in `runs/explore_candidates_population.R`; they converge to four decimals when the grids are doubled and tripled. The original uniform grid was too coarse in the log-normal tail, and its values were 0.01-0.025 too low. Corrected values (naive / regression calibration / SIMEX): linear, $r=0.75$: 0.417 / 0.180 / 0.238; quadratic, $r=0.75$: 0.397 / 0.036 / 0.211; linear, $r=0.9$: 0.241 / 0.144 / 0.083; quadratic, $r=0.9$: 0.216 / 0.078 / 0.072. The $\omega=0.8$, $r=0.9$ row agrees to 0.001; the $\omega=0.8$, $r=0.6$ row was not recomputed. The Section 2 expectations shift accordingly: linear, regression calibration 0.144 and SIMEX 0.083; quadratic, both about 0.075. The study's estimates (0.139 and 0.079 in the linear spline cell) lie closer to the corrected values | Numerical correction, found after the results were seen |
| 2026-10-14 | After the results were seen, `runs/amend4_tail_diag.R` computed robust summaries (median error, MAD), the error by how far a test-fold W_1 exceeds its training folds' maximum, a refit and decomposition of the worst replicate, and task runtimes (`checks/`). In the spline-lasso quadratic cell, 46 of 500 replicates have an error above 1 in at least one method, 43 of them for `latent`. These come from a few controls whose W_1 lies beyond the training range: the natural spline is linear there, so it underpredicts the quadratic trend, and its logit extrapolates to the upper trim. In replicate 389, one control (W_1 = 59.6 against a training maximum of 12.7) moved the estimate to -32.6. These diagnostics are explanatory and are labelled post hoc in the manuscript | Post hoc |
