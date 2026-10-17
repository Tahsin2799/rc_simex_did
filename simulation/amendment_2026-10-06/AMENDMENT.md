# Post-production amendment to the v2.0 simulation study

**Dated 6 October 2026. Status: FROZEN on 6 October 2026, before any confirmatory replicate was run, under the author's written delegation. The author reviews it with the manuscript.**

This amendment was written **after** the v2.0 production results were seen. It is a local, dated record in the repository, not an external preregistration. It does not change the frozen design: tag `be5bff47e5ab11802d8146f044cdf2e6075139cc`, [ADEMP_prereg_v2.md](../ADEMP_prereg_v2.md), 64 cells × 500 replicates. Production results stay as recorded, including 1,120 `ipt_fallback` rows and one `match_full` error.

## 1. Why an amendment

The post-production audit ([post_prod_sim_steps.md](../../post_prod_sim_steps.md), Sections 3-4 and 11.C) found that the theory does not justify the production SIMEX standard error. The cause is full-sample calibration combined with pooled normalisation, so the production intervals are not validated by the theorem. The audit also found that exact quadratic extrapolation needs an explicit approximation allowance. Section 11.C.S4.2 derives a modified estimator that has a complete conditional expansion. This amendment adds that estimator and the studies needed to evaluate it. It does not re-score the production method.

## 2. Changes to reporting only (no refits)

[runs/aggregate.R](../runs/aggregate.R) was corrected on 6 October 2026, and [runs/report_tables.R](../runs/report_tables.R) writes [tables_v2/](../tables_v2/). The backup `sim_results_v2/_summary_pre_2026-10-06.rds` keeps the pre-fix summary.

| Repair | Effect |
|---|---|
| Bias MCSE from sd(est - truth) | Largest change 0.0000245 |
| Wilson intervals for coverage; actual interval denominators | 0 % and 100 % cells no longer report an MCSE of zero |
| CI-width MCSE; runtime median and 90th percentile | New columns |
| `bounds_*` midpoint removed from bias/RMSE | Range estimators are reported by endpoint bias, width, set inclusion and interval coverage |
| Population ATT beside the realised target | Nonlinear population ATT = 0.199787. Truncation and rounding of c12, c14 and c15 move it off 0.2 |
| Paired comparisons on clean and all-usable bases | `tables_v2/paired_comparisons.csv` |

Headline numbers are unchanged; see `tables_v2/headline_N5000.csv`.

## 3. New method and the work already done

- **Method:** `simex_fold_3.0`, specified in [SPEC_simex_fold.md](SPEC_simex_fold.md) and implemented in [34_simex_fold.R](../R_scripts/34_simex_fold.R). It is separate from the frozen `32_simex.R`.
- **Implementation checks (complete).** [checks/checks_report.txt](checks/checks_report.txt): 48 of 48 pass.
  - Exact algebra.
  - Large-n closed-form agreement for Q_j and the calibration gradient D.
  - A 400-replicate Monte Carlo in a theorem-compatible Gaussian design with a correctly specified parametric learner. SE/SD = 0.97-0.98, and coverage of q2(-1) is 0.94.
  - In the same design, exact quadratic extrapolation (M = 0) covers the true ATT only 38-40 % of the time. The true third-derivative bound restores coverage, but its allowance (2.22) is 14 times the actual extrapolation bias (0.16). A global M bound can be very conservative.
- **Diagnostic pilot (complete; not confirmatory).** [pilot/pilot_report.txt](pilot/pilot_report.txt) covers `x_classical__r0.60__nonlinear`, lasso, production seeds, 40 replicates at N = 1,000 and 2 timing replicates at N = 5,000.
  - Replaying production `simex_lasso` reproduced the saved estimate and SE exactly in 3 of 3 replicates.
  - No failures; at most 0.1 % of propensities clipped.
  - Runtime matches production: median 54 s at N = 1,000 and 144-159 s at N = 5,000 per fit.
  - The calibration term's share of the q2(-1) influence variance had a median of 0.16 % (N = 1,000) and 0.03 % (N = 5,000). In this design, the calibration correction to the SE is negligible. The open question is extrapolation bias: realised-target bias was about -0.035 to -0.040 in all three arms.
  - With 40 replicates, coverage is resolved only to about ±0.07, so no coverage claim follows.

## 4. Confirmatory design (frozen)

The pilot raised no change to the method. The only change before freezing is computational: kappa(t) on the range grid is cached, and the 48 checks were re-run with identical results. The runner is [runs/confirm_simex_fold.R](../runs/confirm_simex_fold.R) and the aggregator is [runs/confirm_aggregate.R](../runs/confirm_aggregate.R). There are 500 replicates per cell. Arms `fixed` (population reference covariance) and `fold` (estimated per training fold) share each replicate's data and folds. SIMEX uses B = 20, K = 5 and grid (0, 0.5, 1, 2). The reference reliability is r_L = 0.6 throughout.

| Study | Cells | Learners | Question | Targets |
|---|---|---|---|---|
| `gauss_known` | closed-form design, r0 = r_L = 0.6, n = 1,000 and 5,000 | glm (both n); lasso (n = 1,000) | Does the expansion deliver nominal coverage when its assumptions hold? | q2(-1) (exact), ATT, the whole q2 curve on [-1, 0] |
| `gauss_range` | closed-form design, r0 = 0.75, reference 0.6 (reliability treated as unknown in [0.6, 1]), n = 1,000 and 5,000 | as above | Range interval (S24) and outer ATT interval (S31) | range of q2 on [-1, 0], ATT at t0 = -0.625 |
| `dgp_simex` | production DGP, `x_classical` nonlinear: r0.60 N1000, r0.60 N5000, r0.75 N1000; production seeds | lasso (`sq` dictionary) | Behaviour on the production design, paired with production `simex_lasso` / `bounds_lasso` rows | realised ATT (registered), population ATT |
| `oracle` | production DGP, `x_classical` r0.60 and `treat_misclass` large, linear and nonlinear, N5000; production seeds | `dml_lasso`, `dml_xgb` on observed and on latent data | Same-sample naive-minus-oracle gap | realised ATT |

**Intervals recorded per replicate:**
- Pointwise Wald at t_point (-1, or -(1 - r0)/(1 - r_L) when r0 differs from r_L), with the fixed-calibration and calibration-aware SEs.
- S29 widening for M ∈ {0, 0.005, 0.01, 0.02, 0.05}, plus the true M in the closed-form design.
- S24 / S31 simultaneous-band range and outer intervals for the same M values, with the closed-form conservative version.

M is in dY units per (lambda unit)^3 of the reference covariance. It is a sensitivity grid; no M is estimated.

**Seeds and outputs:**
- Data seeds are production's `seed_for(cell, rep)` for production cells and `digest2int(cell | rep | data)` otherwise. Fold and method seeds are name-based.
- Rowwise influence values are saved for replicates 1-20 only, to limit disk use.
- Errors are recorded per task, never retried silently.

**Primary metrics** are bias, empirical SD, SE/SD ratio, coverage with Wilson intervals, width, the calibration share of variance, clipping, runtime, and paired differences (fold - fixed; new arm - production). Coverage of a target is reported only for the interval that targets it.

**Compute:** about 346,000 task-seconds, roughly 8.7 h wall on 11 workers, from measured per-task times.

**Not included (decided at freeze):**
- Estimated reliability with its own influence h_r.
- External calibration samples (S21).
- xgboost SIMEX.
- Error-model stress cells, BCF diagnosis, nonlinear assignment, and the correlated-heterogeneity tau_R study.

The manuscript reports these as limitations rather than results.

## 5. Deviations log

| Date | Deviation | Reason |
|---|---|---|
| 2026-10-06 | Reporting repairs in Section 2 | Audit findings; no estimates changed |
| 2026-10-06 | Added `simex_fold_3.0` and its checks and pilot | Post-production mathematical repair (Section 11.C) |
| 2026-10-06 | Before freezing: kappa cached in `simex_fold_range`; added `gauss_range` and `oracle` studies to the proposal | Runtime; plan Section 6 P0 items (range coverage, same-sample oracle) |
| 2026-10-06 (after launch, before any `bihs_semi` replicate) | Added study `bihs_semi`: BIHS-calibrated semi-synthetic design ([03_semisynth_population.R](../../application/R_v2/03_semisynth_population.R)), 500 reps; villages resampled, known ATT 0.25, reliability 0.7 on log assets and log land, cluster folds and village-clustered variance. Added `cluster` support to `34_simex_fold.R` (identical results without clusters; 50/50 checks) | Plan Section 6 P2 and the cluster extension needed for BIHS. Frozen studies unchanged |
| 2026-10-06 (after launch) | Runner dispatch changed to `chunk.size = 1`; run restarted, keeping completed replicate files | The default chunking gave single workers about 24 h of N = 5,000 tasks. No estimate is affected |
| 2026-10-07 (after the run, results seen) | `confirm_aggregate.R` adds paired contrasts for `bihs_semi` (SIMEX - naive, SIMEX - latent, naive - latent, per arm) to `confirm_paired.csv`. The 7,000 task files and every existing row are unchanged | Reporting only. The registered naive and latent biases have MCSE about 0.006, too coarse to resolve differences of about 0.02; the paired MCSE is about 0.0015 |
| 2026-10-08 (results seen) | `confirm_aggregate.R` adds runtime median and 90th percentile, total task-hours (`confirm_compute.csv`) and the median clipping share per cell | Reporting only, under the second amendment, Section 4 |
| 2026-10-08 | Follow-up studies (BCF pilot, regression-calibration comparator, tau_R truth) specified and frozen in a separate, second amendment: [AMENDMENT_2.md](../amendment_2026-10-08/AMENDMENT_2.md) | Open plan items; this amendment's frozen design is unchanged |

## 6. Sign-off

Confirmatory design frozen on 6 October 2026, before any confirmatory replicate ran. Claude froze it under the author's written delegation ("I am giving you the permission of taking the liberty of finishing it"). The author reviews it together with the manuscript. Any later change goes in the deviations log.
