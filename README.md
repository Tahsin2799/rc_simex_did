# Regression calibration and SIMEX for mismeasured covariates in doubly robust difference-in-differences

Code and saved results for the paper of the same title (Tahsin Shahriar). They reproduce every number, table and figure in the paper and its supplement.

## Contents

```
simulation/
  R_scripts/            data-generating processes, DR-DiD estimators, SIMEX and regression calibration
                        with cross-fitted calibration, population targets
  runs/                 study runners, aggregators and the numerical checks of Appendix H.2
  audit_2026-10-06/     two checks of the identification results (Appendix B)
  ADEMP_prereg_v2.md    design of the original 64-cell study, fixed before it ran
  amendment_*/          the dated protocol of each later study (AMENDMENT*.md), with its aggregated results
  tables_v2/            aggregated results of the original study
  population_map/       exact population targets and their checks (Section 4)
  bias_aware/           intervals with the bias-aware critical value, recomputed from the saved estimates (Section 5.6)
application/
  R_v2/                 BIHS application: panel construction, analysis, semi-synthetic population, estimated gap shares
  output/               aggregate tables of the application (no unit-level data)
manuscript/             sources of the paper and supplement, the exhibit builder and the build scripts
renv.lock, renv/        the package environment (R 4.5.2)
```

## Environment

Install R 4.5.2, open R in this folder and run `renv::restore()`. The `.Rprofile` files activate the environment, including from `simulation/`.

## Rebuild the paper from the saved results

From this folder:

```sh
Rscript manuscript/R/build_exhibits.R      # tables, figures and inline numbers -> manuscript/exhibits, manuscript/figures
Rscript manuscript/R/build_submission.R    # main text, supplement and the combined document (PDF)
```

The builder reads only the aggregated results in `simulation/` and `application/output/`, so this step takes about a minute and needs no data. The PDFs need a LaTeX installation with `pdflatex` and `bibtex`.

## Rerun the studies

Every simulation script runs from `simulation/`. The full pipeline, in order, with the output of each step in comments:

```sh
cd simulation
# Original study
Rscript runs/launch.R                      # -> sim_results_v2/
Rscript runs/aggregate.R                   # -> _summary.rds
Rscript runs/report_tables.R               # -> tables_v2/
# First amendment (amendment_2026-10-06/)
Rscript runs/simex_fold_checks.R           # -> checks/
Rscript runs/pilot_simex_fold.R            # -> pilot/
Rscript runs/confirm_simex_fold.R          # -> confirm/ (bihs_semi needs 03_semisynth_population.R)
Rscript runs/confirm_aggregate.R           # -> confirm_metrics_long.csv
# Second amendment (amendment_2026-10-08/)
Rscript runs/followup_mcmc_diag_checks.R   # -> checks/
Rscript runs/followup_rc_checks.R
Rscript runs/followup_bcf_pilot.R          # -> bcf/
Rscript runs/followup_rc.R                 # -> rc/
Rscript runs/followup_aggregate.R
Rscript runs/followup_tau_r.R              # -> tau_r/
# Third amendment (amendment_2026-10-11/)
Rscript runs/rc_fold_checks.R              # -> checks/
Rscript runs/amend3_study.R                # -> study/
Rscript runs/amend3_aggregate.R
Rscript runs/amend3_rc_variance_diag.R     # post hoc (Section 5.2)
# Fourth amendment (amendment_2026-10-11b/)
Rscript runs/amend4_population_map.R       # population calculations
Rscript runs/amend4_population_lnorm.R
Rscript runs/amend4_population_M.R
Rscript runs/amend4_study.R                # -> study/
Rscript runs/amend3_aggregate.R --root=amendment_2026-10-11b/study
Rscript runs/amend4_tail_diag.R            # post hoc -> checks/
# Fifth amendment (amendment_2026-10-14/)
Rscript runs/amend5_learner_check.R
Rscript runs/amend5_study.R                # -> study/
Rscript runs/amend3_aggregate.R --root=amendment_2026-10-14/study
Rscript runs/amend5_compare.R              # -> study/amend5_compare.csv
# Post hoc (population targets, gap share, bias-aware critical value)
Rscript runs/population_map.R              # -> population_map/
Rscript runs/population_map_convergence.R
Rscript runs/population_map_mc_check.R
Rscript runs/gap_bias_check.R
Rscript runs/gap_diagnostic_check.R
Rscript runs/gap_original_design.R
Rscript runs/gap_resolution.R
Rscript runs/bias_aware_intervals.R        # -> bias_aware/
# Numerical checks of the derivations (Appendix H.2)
Rscript runs/math_foundation_repair_checks.R
Rscript runs/math_simex_repair_checks.R
Rscript runs/math_misclass_repair_checks.R
Rscript audit_2026-10-06/check_identification_heatflow.R
Rscript audit_2026-10-06/check_binary_support.R
```

Seeds are derived from the cell and the replicate, so any replicate can be replayed on its own. Raw replicate files are not included; the runners regenerate them, and the aggregators rebuild the saved results. The compute is substantial: the original study used about 114 estimator-days, and the first, third, fourth and fifth amendments 113, 89, 335 and 53 task-hours.

`runs/launch.R` is the driver for the original study (`--jobs` sets the number of workers; finished batches are skipped, so rerunning the same command resumes). `runs/launch.sh` is the same driver for macOS and Linux.

## Data

The application uses the Bangladesh Integrated Household Survey (rounds 2011–12, 2015 and 2018–19), distributed by the International Food Policy Research Institute (IFPRI). The microdata are not included and may not be redistributed. To rerun the application, obtain the rounds from IFPRI and set `BIHS_DATA_ROOT` to the folder that holds them:

```sh
BIHS_DATA_ROOT=/path/to/bihs Rscript application/R_v2/01_build_panel.R
Rscript application/R_v2/02_analysis.R
Rscript application/R_v2/03_semisynth_population.R
Rscript application/R_v2/04_rc_gap.R
```

The scripts write unit-level files to `application/derived/`, which should not be shared, and aggregate tables to `application/output/`. The semi-synthetic study of the first amendment also needs `application/derived/semisynth_population.rds` from the third script.

## Licence

MIT (see `LICENSE`).
