# Fifth post-production amendment: the spline cells rerun with constant extrapolation

**Dated 14 October 2026. Status: FROZEN on 14 October 2026, before any study replicate ran. The author asked for the rerun on 14 October 2026 ("i still want to run it") after Claude had recommended against it (reasons in Section 1). Claude chose the learner change and froze this record. The author reviews it with the manuscript.**

This amendment was written after all of the fourth amendment's results had been seen ([AMENDMENT_4.md](../amendment_2026-10-11b/AMENDMENT_4.md), including its deviations log). It is a **sensitivity analysis**. The fourth amendment's spline results remain the primary results and are reported in full, whatever this study shows.

## 1. Why, and what was known before freezing

In the fourth amendment's spline-lasso cell with a quadratic trend (`lnorm_quad_r090_n5000`, learner `glmnet_ns`), 46 of 500 replicates have an error above 1 in at least one method, 43 of them for `latent`. The mean-based metrics are dominated by these replicates (RMSE: latent 2.37, `rc_fold` 1.98, `simex_fold` 6.14). The post hoc medians are latent -0.040, naive 0.216, `rc_fold` 0.072 and `simex_fold` 0.049.

The post hoc diagnosis ([amend4_tail_diag.R](../runs/amend4_tail_diag.R), `../amendment_2026-10-11b/checks/`):
- With log-normal covariates (omega = 1.0), some test-fold W_1 exceeds its training folds' maximum in every replicate.
- The natural spline is linear beyond its boundary knots, which `splines::ns` places at the range of the training data. Beyond that range, the outcome model underpredicts a quadratic trend, and the propensity logit can extrapolate past the upper trim of 0.98.
- In replicate 389, one control (W_1 = 59.6, against a training maximum of 12.7) had a residual of 1,318 and a normalised weight of 125, which moved the estimate to -32.6.
- Replicates with a test-fold excess above 10 account for 16 of the 43 latent failures; none occurred with an excess below 2.

Claude recommended against a rerun:
- the fourth amendment's headline comes from the linear spline cells, which have no replicate with an error above 1 at n = 5,000;
- the oracle fails in the same way, so the cell is not evidence about either correction;
- any learner change is chosen after seeing the failure.

The author chose to run it as a sensitivity analysis.

**What had been computed with the new learner before freezing:** the code check in Section 3 (predictions only, no estimator) and a six-task pilot used for timing and errors. The pilot's estimates were not summarised or read. The change was chosen from the mechanism above and was not tried on the failing replicates.

## 2. The learner change

`glmnet_nsc` ([42_spline_clamped.R](../R_scripts/42_spline_clamped.R)) is `glmnet_ns`, except that each column of a new design matrix is clamped to the training data's range before the spline basis is evaluated:
- inside the training range, the basis, the fit and every prediction are unchanged;
- beyond it, the fitted function is constant at its boundary value, as a tree learner would extrapolate;
- one change applies to both nuisances and to every method, including `latent`.

**Alternatives considered and not chosen:**
- *Boundary knots over the full sample (training and test covariates):* a cubic piece is still fitted over a range that holds no training data.
- *A tighter propensity trim:* `PS_TRIM` = (0.02, 0.98) is fixed for the whole paper.
- *Winsorising at training quantiles, such as the 1st and 99th percentiles:* this changes fits inside the range as well.
- *Adding W_1^2 to the spline dictionary:* this makes the learner closer to the lasso with squares, which the fourth amendment already runs, and changes what the spline learner is meant to approximate.

**What clamping can and cannot fix:**
- Clamping bounds the propensity of an out-of-range unit at its boundary value.
- It does not fix the outcome model. Under a quadratic trend, constant extrapolation underpredicts more than linear extrapolation does, so an out-of-range control's residual grows.
- Large errors can therefore remain whenever an out-of-range control has a large boundary propensity.
- In the linear cells, linear extrapolation of the outcome was correct, so clamping makes those few predictions worse. In the code check, the out-of-range outcome predictions moved by up to 8.4.

## 3. Study `shape4`, learner `glmnet_nsc`

**Cells:** `lnorm_linear_r090_n5000`, `lnorm_quad_r090_n5000` and `lnorm_linear_r090_n1000`, the fourth amendment's three spline cells. 500 replicates each, 1,500 tasks ([amend5_study.R](../runs/amend5_study.R), copied from `amend4_study.R`).

**Pairing:**
- Data and folds come from `seed_of(cell, rep, "data")` and `seed_of(cell, rep, "folds")`, as before.
- Method seeds use the learner label `glmnet_ns`, so the cross-validation folds of `cv.glmnet` and the SIMEX noise draws are also those of the fourth amendment.
- Each replicate is therefore paired with its fourth-amendment counterpart, and they differ only in predictions at covariate values outside a training fold's range.

**Methods:** `naive`, `latent`, `rc_group`, `rc_fold` and `simex_fold` (M = 0 and 0.05), unchanged.

**Code check** ([amend5_learner_check.R](../runs/amend5_learner_check.R), run 14 October 2026: PASS).
- Design: fold 1 of replicate 1 in each cell; observed and latent covariates; both nuisances; same RNG state.
- Inside the training range, `glmnet_nsc` and `glmnet_ns` predictions are identical (largest difference 0).
- Beyond it, `glmnet_nsc` predictions equal `glmnet_ns` predictions at the clamped point.

**Pilot** (`pilot/`, 2 replicates per cell, 6 workers): 0 errors; 40-42 s per task at n = 1,000, 129-133 s at n = 5,000. Expected compute is about 48 task-hours, about 5 hours on 10 workers.

**Expectations, written before the study** (no decision rule; every cell is reported):
- Quadratic cell:
  - fewer replicates with an error above 1 than the fourth amendment's 46, but probably not none (Section 2);
  - mean errors closer to the fourth amendment's medians than its means were.
- Linear cells: estimates within a few hundredths of the fourth amendment's, with SIMEX still closer to the ATT than regression calibration at n = 5,000 (fourth amendment: 0.079 against 0.139).

**Metrics and contrasts:**
- The fourth amendment's, from `runs/amend3_aggregate.R` pointed at this root, together with the robust summaries of `runs/amend4_tail_diag.R` (median error, MAD, count of errors above 1).
- Per method and cell, a paired comparison with the fourth amendment's `glmnet_ns` rows, replicate by replicate: mean and median of the change in error, and replicates whose error crosses 1 in either direction.
- The comparison script `runs/amend5_compare.R` will be written while the study runs and tested on the pilot output only.

## 4. Deviations log

| Date | Deviation | Reason |
|---|---|---|
| 2026-10-14 05:15-10:34 | The study ran in 5.32 h on 10 workers: 1,500/1,500 tasks, 0 errors, no restarts | Record |
| 2026-10-14 12:18 | `runs/amend5_compare.R` was written after the study finished, not while it ran as Section 3 planned, but before any study estimate was read. It was tested on the pilot output first, and it computes the metrics listed in Section 3 | Timing of the analysis code only |

## 5. Results (14 October 2026)

`study/amend3_*` (standard metrics) and `study/amend5_compare.{csv,txt}` (paired with the fourth amendment's `glmnet_ns` rows, by replicate).

- **Quadratic cell.** Errors above 1 fell for every method but did not disappear:

  | Method | Old | New |
  |---|---|---|
  | latent | 43 | 40 |
  | naive | 14 | 8 |
  | `rc_fold` | 17 | 9 |
  | `rc_group` | 21 | 10 |
  | `simex_fold` | 24 | 22 |

  RMSE roughly halved (`rc_fold` 1.98 → 0.78; `simex_fold` 6.14 → 2.16; latent 2.37 → 1.17). Medians are unchanged: latent −0.041, `rc_fold` 0.074, `simex_fold` 0.051. As Section 2 anticipated, clamping bounds the propensities but not the outcome residuals of out-of-range units, which still dominate.
- **Linear cells, $n=5{,}000$.** Changes of −0.002 to −0.006 (paired MCSE ≤ 0.003). Regression calibration 0.136, SIMEX 0.073: SIMEX is still closer to the ATT, as expected.
- **Linear cells, $n=1{,}000$.** Changes of −0.017 to −0.032. The latent RMSE rose from 0.147 to 0.195 and SIMEX had five errors above 1 instead of two, which is consistent with constant extrapolation being wrong for a linear trend.
- **Reading.** The sensitivity analysis leaves the fourth amendment's conclusions unchanged. The quadratic spline cell remains a learner failure that no change to the propensity extrapolation fixes.
