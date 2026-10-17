# ADEMP Pre-Registration, version 2.0

**Title.** Measurement error in machine-learning difference-in-differences: what theory settles, what it does not, and a SIMEX correction with reliability bounds.

**Authors.** Tahsin Shahriar (corresponding).

**Version.** 2.0, drafted 2026-09-12. Supersedes v1.3 (`ADEMP_prereg.md`, filed 2026-05-17). This is a **re-scoping after a failed production run**: the v1.3 design was launched on 2026-05-20 and produced no usable output in sixteen days (post-mortem in `SIM_ISSUES_2026-09-11.md`); the corrected compute requirement for v1.3 is 4,593 CPU-days against a budget of roughly 333. Every change from v1.3 is listed in §11 with its reason. The scientific aim is narrower and the theoretical content is larger.

**Status.** Internal design specification, not externally registered. The author decided on 2026-09-12 not to file this study on OSF and to keep the record in this repository instead. The document is therefore binding by convention rather than by timestamp: it is fixed before the production run, every change is dated in §11, and the manuscript claims only that the design was fixed in advance, not that it was registered.

**Date fixed.** [date of the commit that precedes the first production batch]

**Identifier.** Git tag `sim-launch-v2.0` in the replication repository.

**Associated manuscript.** `manuscript/paper_jci.Rmd`.

**Code.** `simulation/R_scripts/` (v2 files: `1_sim_param.R`, `2_data_sim.R`, `30_core.R`, `31_estimators.R`, `32_simex.R`, `33_misclass.R`, `3_registry.R`), `simulation/runs/` (`grid.R`, `run_batch.R`, `launch.sh`, `aggregate.R`, `smoke.R`). v1 code is retained unmodified under `legacy_v1/` for provenance. Git SHA at filing recorded in §10.

Template: Siepe, Heck, Mikus, Pawel & Heycke (2024); framework: Morris, White & Crowther (2019).

---

## 1. Aims

**Primary aim.** For the two-period balanced-panel ATET, characterise which measurement-error (ME) structures the doubly-robust DiD score is *provably* insensitive to, which it is not, and how machine-learning nuisance estimation changes the picture. Four propositions (§4.1) settle the outcome-ME cases analytically; the simulation is reserved for the three cases where the answer depends on the nuisance learner: classical covariate ME, treatment-correlated covariate ME, and treatment misclassification.

**Secondary aim.** Evaluate two corrections that are new in this setting: (i) SIMEX applied to the cross-fit DR-DiD score with ML nuisances, under known reliability; (ii) partial-identification bounds on the ATET when the reliability is only known to lie in $[r_L, 1]$, with Imbens-Manski confidence intervals. And (iii) a covariate-conditional known-rate correction for treatment misclassification.

**Tertiary aim.** Provide a reliability-sensitivity analysis for the BIHS PESP application.

Phase, in the Heinze et al. (2024) taxonomy: Phase II for the comparison, Phase I/II for the corrections (new use of an existing method, first evaluation).

**Not aims.** Quantile estimands; staggered adoption; synthetic-control estimators (the design has one pre-period).

---

## 2. Data-Generating Mechanisms

### 2.1 Covariates

Unchanged from v1.3 §2.1: twenty covariates in four layers, moments from BIHS Round 1 (2011). Layer A ($c_1$-$c_5$) and Layer B ($c_6$-$c_{10}$) are confounders; Layer C ($c_{11}$-$c_{15}$) are effect modifiers under the non-linear family; Layer D ($c_{16}$-$c_{20}$) is noise. The ML design matrix has 27 columns ($k-1$ dummies for $c_9, c_{10}, c_{20}$).

### 2.2 Treatment

Unchanged: $A_i \sim \mathrm{Bernoulli}(\mathrm{logit}^{-1}(\kappa_0 + \kappa^\top X_i^{(A,B)}))$, $\kappa_0 = -0.55$, $P(A=1) \approx 0.29$ (realised 0.282-0.290 in the v2 smoke test).

### 2.3 Outcome (changed)

$$
Y_{i0} = \lambda_0 + \lambda_a A_i + \beta^\top X_i^{(A,B)} + n_\ell(X_i) + \alpha_i + \varepsilon_{i0},
\qquad
\Delta Y_i = \lambda_t + \delta^\top X_i^{(A,B)} + n_t(X_i) + \theta_i A_i + (\varepsilon_{i1} - \varepsilon_{i0}),
$$

with $\delta = \beta/2$ (**new**), $n_\ell = n_t = 0$ under the linear family, and under the non-linear family $n_\ell(X) = 0.19\,h(X)$, $n_t(X) = 0.10\,h(X)$ where $h$ is the standardised surface $0.10 z_1 z_6 + 0.15\sin z_3 + 0.10(z_4^2-1) + 0.10 c_2 z_3 + 0.08 c_7 z_6$ (scaled to unit SD), and $\theta_i = \theta(1 + \sum_k g_k z_{k})$ over Layer C as in v1.3.

**Why $\delta \ne 0$.** In v1.3 covariates entered only the level, so they cancelled in $\Delta Y$, unconditional parallel trends held, and a covariate-blind DiD was unbiased (verified: 400 reps, mean 0.202). Under that DGP covariate ME cannot bias any $\Delta Y$-based estimator, which would have made the paper's central cells uninformative. With $\delta = \beta/2$, conditional parallel trends holds and unconditional fails: the covariate-blind contrast has bias $E[\delta^\top X \mid A=1] - E[\delta^\top X \mid A=0] \approx -0.087$ (300 reps, $N=1000$), i.e. about $-0.44\,\theta$, while the parametric DR-DiD is unbiased (mean 0.199, empirical SD 0.072, mean model SE 0.074).

Parameters: $\lambda_0 = 5$, $\lambda_t = 0.2$, $\lambda_a = 0$, $\theta = 0.2$, $\sigma_\alpha = 0.7$, $\sigma_\varepsilon = 0.7$ (panel-strength "low" only; see §11).

### 2.4 Realised ATET

$\tau = E[\theta_i \mid A_i = 1]$, computed exactly per replicate as the treated mean of $\theta_i$ (this equals the counterfactual contrast because $\varepsilon$ cancels). Under the linear family $\tau = 0.2$.

### 2.5 Measurement-error designs (changed)

Let $W = (c_1, c_3, c_4, c_6, c_8)$ be the error-prone continuous confounders and $V$ the rest.

| Design | Wrapper | Parameterisation | Cells |
|---|---|---|---|
| Baseline | none | — | 1 |
| **Classical covariate ME** | `me_x_classical` | reliability $r = \mathrm{Var}(W_k)/\mathrm{Var}(W_k^*)$, $U_k \sim N(0, \mathrm{sd}(W_k)^2(1/r-1))$, same draw both periods | $r \in \{0.90, 0.75, 0.60\}$ |
| **Treatment-correlated covariate ME** | `me_x_nc` | $W_k^* = W_k + \mathrm{sd}(W_k)(\delta_a A + N(0, 0.1^2))$, per-variable scaling (v1.3 added an unscaled scalar) | $\delta_a \in \{0.1, 0.3, 0.5\}$ |
| **Treatment misclassification** | `me_treat_misclass` | $(\alpha_0, \alpha_1)$, $X$-independent | $\{(.02,.02), (.05,.10), (.10,.20)\}$ |
| Non-linear covariate-dependent outcome ME | `me_y_nc_covariate` | $Y_1^* = Y_1 + \rho\,\sigma_\varepsilon\,h(X) + N(0, 0.35^2)$ at $t=1$ only; $h$ as in §2.3 (v1.3 used raw child age, mean 28) | $\rho \in \{0.25, 0.5, 1\}$ |
| Classical outcome ME | `me_y_classical` | $\sigma = 0.7$ | 1 (proposition check) |
| Differential outcome ME | `me_y_diff` | $(\sigma_0, \sigma_1) = (0.7, 1.75)$ | 1 (proposition check) |
| Contraction outcome ME | `me_y_nc_contraction` | $(\gamma, \sigma) = (-0.3, 0.35)$ | 1 (proposition check) |

Sixteen ME cells. ME is applied to the observed data only.

### 2.6 Grid

$N \in \{1000, 5000\}$ × DGP family $\in \{\text{linear}, \text{nonlinear}\}$ × 16 ME cells = **64 cells**. $N = 5000$ is the ceiling because the BIHS child panel has about 2,300 children per wave; the v1.3 $N = 10000$ tier is dropped (§11).

---

## 3. Estimand

ATET in the post period, $\tau = E[Y_{i1}(1) - Y_{i1}(0) \mid A_i = 1]$. For the partial-identification estimator the target is the same $\tau$; its performance measure is coverage of $\tau$ by the Imbens-Manski interval and the width of the estimated identified set.

---

## 4. Methods

### 4.1 Propositions (pre-specified; proofs in the manuscript)

Notation: $\Delta Y = \mu(X) + \theta(X) A + \eta$ with $E[\eta \mid X, A] = 0$ (conditional PT); $e(X) = P(A=1\mid X)$; the DR-DiD functional is $\tau(e, m_0) = E[(w_1 - w_0)(\Delta Y - m_0(X))]$ with $w_1 = A/E[A]$, $w_0 \propto (1-A)e(X)/(1-e(X))$.

- **P1 (outcome ME).** (a) $\Delta Y^* = \Delta Y + U$, $U \perp (A, X)$, $E[U]=0$: $\tau$ unchanged; only variance rises. (b) $U$ heteroskedastic in $A$: same. (c) $\Delta Y^* = (1+\gamma)\Delta Y + U$: every estimator linear in $\Delta Y$ that identifies $\tau$ without ME identifies $(1+\gamma)\tau$; no covariate adjustment changes this. (d) $\Delta Y^* = \Delta Y + \rho g(X)$ for arbitrary $g$: the DR-DiD functional with the **correct propensity score** is unchanged for any $m_0$, because $E[(w_1 - w_0) g(X)] = 0$ when $w_0$ balances $X$; with a misspecified PS the functional is unchanged only if $m_0$ absorbs $g$. (This is what makes the `y_nc_covariate` cells discriminating between nuisance learners.)
- **P2 (covariate ME: no protection).** With $W^* = W + U$, $U \perp (W, V, A, \eta)$, and $X^* = (W^*, V)$, the DR-DiD functional evaluated with nuisances that are correct *for $X^*$* equals $\tau + b(r)$, $b(r) = E\big[E[\mu(X)\mid X^*, A=1] - E[\mu(X)\mid X^*, A=0] \mid A=1\big]$ (homogeneous effect), which is zero iff $r = 1$, or $\mu$ does not depend on $W$, or $A \perp W \mid W^*$. Double robustness protects against misspecification of $(e^*, m_0^*)$, not against $X^* \neq X$. In the linear-Gaussian case $b(r) = \delta_W^\top E[E[\epsilon_W \mid W^*, A=1] - E[\epsilon_W \mid W^*, A=0] \mid A=1]$ where $W = E[W\mid W^*] + \epsilon_W$ and $\mathrm{Var}(\epsilon_W) = (I - \Lambda(r))\Sigma_W$; $b$ is continuous in $r$ with $b(1) = 0$.
- **P3 (SIMEX and bounds; paper Theorems 1-2, proofs in paper Appendix A.7-A.8).** Let $\tau^*(\lambda)$ be the DR-DiD functional on $W^* + \lambda^{1/2} U'$, $U' \sim$ the ME law. If $\tau^*$ is quadratic on $[-1, \lambda_{\max}]$, SIMEX extrapolated to $\lambda = -1$ identifies $\tau$. For fixed $\lambda$ and $B$, the SIMEX estimator is a fixed linear map of averages of cross-fit DR-DiD estimators, hence asymptotically linear under the Chang (2020) nuisance-rate conditions; its variance is estimated by the influence-function estimator of the extrapolated linear combination (exact for the combination), of which the Stefanski-Cook estimator is an approximation. Under unknown common reliability $r \in [r_L, 1]$, the identified set is $\{\tau^*(\lambda(r)) : r \in [r_L, 1]\}$ with $\lambda(r) = -(1-r)/(1-r_L)$, an interval by continuity; its endpoint estimators are of the same form, and the Imbens-Manski (2004) / Stoye (2009) interval has asymptotic coverage at least $1-\alpha$ for $\tau$. Deviations of $\tau^*$ from quadratic are a bias of the same order as in every SIMEX application and are reported empirically.
- **P4 (misclassification under conditional PT).** With non-differential known rates $(\alpha_0, \alpha_1)$, $\kappa = 1 - \alpha_0 - \alpha_1$: $q(X) := P(A^* = 1 \mid X) = \alpha_0 + \kappa e(X)$ and $E[\Delta Y \mid A^*=1, X] - E[\Delta Y \mid A^*=0, X] = \theta(X) f(X)$ with $f(X) = \kappa\, e(1-e)/\{q(1-q)\}$. $f$ is a non-constant function of $e$; when $e \equiv \pi$ it equals the Aigner (1973) factor $\kappa\pi(1-\pi)/\{\pi^*(1-\pi^*)\}$, which equals $\kappa$ only if $\pi^*(1-\pi^*) = \pi(1-\pi)$ (the earlier draft's claim that it reduces to $\kappa$ was wrong). Every DR-DiD-family estimator with a nuisance correct for $(A^*, X)$ converges to $\tau^\circ = E[\theta f \mid A^* = 1]$, so the bias depends on the overlap structure, not only on the rates. Identification: $e = (q - \alpha_0)/\kappa$, $\theta = c/f$, $\tau = E[c\,H(q)]/E[e]$ with $H(q) = q(1-q)/(1-\alpha_1-q)$; $\tau_R = E[c\,e]/E[fe]$ is the $fe$-weighted ATET. Both are estimated by the Neyman-orthogonal cross-fit scores of paper Theorem 3 (see §4.3).

Positioning: Negi & Negi (2025, JAE) treat differential one-sided misclassification with a parametric two-step; Denteh & Kédagni (2022/2026) give unconditional bounds under known rates. P4 is the covariate-conditional case and the question the simulation asks is how ML nuisances on $A^*$ behave.

### 4.2 Base estimator panel (12)

All receive the 27-column design matrix; none receives $y_0$. The lasso learner does not see the 27 columns directly: it is fed the standard "flexible" dictionary for linear nuisance learners in double machine learning, the raw columns plus squares of the 12 non-binary columns (39 columns), with continuous columns standardised on the training fold, the penalty chosen by 5-fold cross-validation inside each training fold (`lambda.min`, 50-point path down to $0.01\,\lambda_{\max}$), and dummy-coded factors entering linearly (V13; precedent: Chernozhukov et al. 2018 Sec. 6 and the DoubleML 401(k) example use raw columns plus degree-2 polynomials of the continuous covariates for the lasso and raw columns for trees; Chetverikov, Liao and Chernozhukov 2021 give near-oracle rates for the cross-validated penalty). A second lasso entry adds all pairwise products of the 17 non-dummy columns (175 columns; the interaction dictionary of Zimmert 2020 and Knaus, Lechner and Strittmatter 2021) and runs in the base panel only, because at four to six times the cost of the squares dictionary it is not affordable inside the 61-refit SIMEX wrapper (§10). Tree learners receive the raw columns; they form interactions internally. All return $(\hat\tau, \widehat{se}, ci_{lo}, ci_{hi})$; coverage is scored on the interval.

| # | Label | Implementation | Interval |
|---|---|---|---|
| 1 | `did_ols` | $\Delta Y$ on $A$, HC1 | Wald |
| 2 | `did_ipw` | `DRDID::std_ipw_did_panel` (logit; IF includes PS-estimation term) | Wald |
| 3 | `drdid_imp` | `DRDID::drdid_imp_panel` (IPT + WLS) | Wald |
| 4 | `did_ps_rf` | Hajek IPW score, ranger OOB PS, 500 trees, min node 5 | Wald, IF |
| 5 | `did_ps_xgb` | Hajek IPW score, 5-fold cross-fit xgboost PS | Wald, IF |
| 6 | `dml_lasso` | cross-fit Sant'Anna-Zhao score, `cv.glmnet` nuisances on the squares dictionary (39 columns) | Wald, IF |
| 6b | `dml_lasso_int` | same, squares-plus-interactions dictionary (175 columns); base panel only | Wald, IF |
| 7 | `dml_rf` | same, ranger nuisances | Wald, IF |
| 8 | `dml_xgb` | same, xgboost nuisances (hist, depth 3, $\eta = 0.05$, rounds by 20 % holdout early stopping, refit) | Wald, IF |
| 9 | `cf_did` | `grf::causal_forest` on $\Delta Y$, ATT by DR score | Wald |
| 10 | `bcf_did` | `bcf` on $\Delta Y$, logit `pihat`, 250 burn / 1000 draws, 100 + 50 trees | posterior 2.5-97.5 % |
| 11 | `match_full` | optimal full matching on logit PS, weighted $\Delta Y$ on $A$, SE clustered by stratum | Wald |

Cross-fit $K = 5$, stratified on $A$. PS trimmed to $[0.02, 0.98]$. One thread per process.

### 4.3 Corrected estimators (cell-specific)

| Label | Cells | Definition |
|---|---|---|
| `simex_lasso`, `simex_xgb` | `x_classical`, `x_nc` | SIMEX around #6 / #8, $\lambda \in \{0.5, 1, 2\}$, $B = 20$, quadratic extrapolant; **SE = influence-function estimator of the extrapolated linear combination** (paper Theorem 1(iii)); the Stefanski-Cook SE is stored as an attribute for comparison only; assumed reliability = true $r$ in `x_classical`, fixed nominal $0.75$ in `x_nc` (deliberate misspecification) |
| `bounds_lasso`, `bounds_xgb` | `x_classical`, `x_nc` | identified set over $r \in [0.6, 1]$ from one SIMEX curve (reference $r_L = 0.6$, $B = 20$): min / max of the fitted quadratic over a 41-point grid on $[-1, 0]$ **plus its exact vertex when interior** (paper Theorem 2); endpoint SEs by the same influence-function construction at $c(\lambda_L), c(\lambda_U)$; Imbens-Manski 95 % interval; reports (lo, hi, ci) |
| `misclass_lasso`, `misclass_lasso_R`, `misclass_xgb`, `misclass_xgb_R` | `treat_misclass` | **Neyman-orthogonal known-rate correction** (paper Theorem 3) with true $(\alpha_0, \alpha_1)$: cross-fit nuisances $m_1, m_0$ (outcome change by observed group) and $q = P(A^*=1 \mid X)$, $K = 5$ folds stratified on $A^*$, lasso (`glmnet`, squares dictionary as in #6) or boosted (`xgboost`) learner, $\hat q$ clipped to $[\alpha_0 + \kappa\delta, 1 - \alpha_1 - \kappa\delta]$; scores $\phi_K, \phi_L$ of paper eq. (20)-(21); `*_lasso` / `*_xgb` target $\tau$ ($K = H$, exact under heterogeneity), `*_R` variants target $\tau_R$ ($K = e$, $fe$-weighted, exact under a homogeneous effect); IF SE. Replaces the plug-in (linear OR by group, logit $e^*$, bootstrap) drafted earlier the same day, which is misspecified because $m_a = \mu + \theta p_a(X)$ is non-linear in $X$ |

Hyperparameters are frozen at the tagged commit. No tuning after it; any change is a dated entry in §11.

### 4.4 Estimators dropped from v1.3 (with reason)

augsynth (structurally broken at $T_0 = 1$: NA jackknife SE, runtime exponent $> 2.7$); `causalweight::didDML` lasso and RF (replaced by the own cross-fit implementations so that all ML DR-DiD entries share one score; removes the repeated-cross-section label caveat); BCF with RF PS (duplicate); CiC and QDiD (different estimand and identification; $B = 50$ bootstrap gave unreliable coverage); Mahalanobis and NN-PS matching (full matching represents the family); gbm variants (replaced by xgboost).

---

## 5. Performance Measures

Per (cell, estimator) with Monte Carlo SEs (Morris et al. 2019): bias, empirical SE, model SE, RMSE, coverage of the reported interval, interval width, failure and flag rates, mean runtime. For `bounds_xgb` additionally the mean identified-set width and whether $\tau$ lies in the estimated set. Computed by `runs/aggregate.R`.

**Degenerate columns (added 2026-09-13).** `design_matrix()` drops any column with no variation in the sample. `c20` is a Layer D noise variable, absent from both the treatment and the outcome index, whose level 2 has prevalence 0.003; its dummy is therefore empty in $0.997^{1000} = 5.0$ % of $N = 1000$ draws and in $3 \times 10^{-7}$ of $N = 5000$ draws. `lm` and `glm` alias such a column, but `DRDID::std_ipw_did_panel` inverts $X'WX$ by Cholesky factorisation and errors, which cost `did_ipw` six of the 64 $N = 1000$ replicates of the 2026-09-13 dry run (§9). Dropping a constant column is an exact reparameterisation, so no estimate changes: rerunning those six returns finite estimates, and the five replicates flagged below return bit-identical values.

**Replicate outcomes (added 2026-09-13).** Every replicate is one of three, and `runs/aggregate.R` reports all three:

- *ok*: clean run. Bias, empirical SE, model SE, RMSE, coverage, width and runtime are computed on these replicates alone.
- *flagged*: a usable estimate, but the solver warned or fell back to something that is no longer the estimator on the label. Reported as `flag_rate`; the same metrics over ok plus flagged replicates are reported as `bias_all`, `rmse_all`, `coverage_all`.
- *failed*: no usable estimate. Reported as `fail_rate`.

Morris et al. (2019, §5.4) ask for the failure mechanism and for the metrics computed both ways. Only `drdid_imp` is flagged before launch. `DRDID::drdid_imp_panel` returns `ps.flag = 2` when `trust` and `nlm` both fail on the calibration loss and it reverts to the plain logit coefficients, so the estimate is standard DR-DiD rather than the improved estimator that carries the label. This happens in 6.5 % of $N = 1000$ replicates (200-draw check, 2026-09-13), almost all of them where `c20` level 2 has exactly one observation, and it is not removed by the column drop above. Moving from the clean basis to all usable replicates shifts the bias from $-0.0010$ (187 replicates) to $+0.0005$ (200), against an MCSE of 0.0051.

Headline tables: (i) `x_classical` × $N = 5000$ × both families, all estimators including corrections; (ii) `treat_misclass` likewise; (iii) proposition-check cells as one compact table showing the predicted bias against the observed.

No hypothesis tests across estimators; descriptive with MCSEs.

---

## 6. Number of replicates

$n_{sim} = 500$ per cell. MCSE for 95 % coverage $= \sqrt{0.95 \cdot 0.05/500} = 0.0097 < 0.01$. MCSE for bias with empirical SD 0.035 at $N = 5000$ is 0.0016. No downgrade rule; no upgrade.

---

## 7. Reproducibility

Seed $= \mathrm{digest2int}(\text{cell\_id}) + r$; estimator $i$ in replicate $r$ is seeded at seed $+ i$. Jobs are (cell, batch of 25 reps) in independent single-threaded R processes; any batch is replayable from `Rscript runs/run_batch.R --cell=<id> --batch=<k>`. `renv.lock` at the repo root; Zenodo release at submission.

---

## 8. Analysis plan

1. `runs/launch.sh` runs all pending batches (resumable). 2. `runs/aggregate.R` binds batches and writes `_summary.rds`. 3. Tables per §5; nested-zip plots per ME design faceted by $N$ and family. 4. Proposition checks: observed bias vs. predicted (0, 0, $\gamma\tau$, 0-with-correct-PS). 5. SIMEX: bias and coverage vs. naive; bounds: coverage of $\tau$ and set width vs. $r$. 6. Misclassification: naive attenuation across estimators vs. P4's $X$-dependent factor; corrected estimator bias and coverage. Anything else is labelled exploratory.

---

## 9. Pre-run checks (done 2026-09-12; dry run 2026-09-13)

- Smoke test (`runs/smoke.R`, $N = 1000$): all 12 base and 4 corrected estimators return finite estimates in the baseline, `x_classical` and `treat_misclass` cells, both families.
- DGP check: 300 reps, `drdid_imp` mean 0.1986 (truth 0.2), empirical SD 0.0723, mean SE 0.0735; `did_ols` mean 0.1127, PT gap $-0.086$.
- Misclassification correction, 20 draws at $(0.05, 0.10)$: naive adjusted 0.160, ratio correction 0.195 (SD 0.074), per-unit form 0.220 (SD 0.129) — ratio form adopted.
- Orthogonal misclassification correction, oracle diagnostic (one draw, $N = 20000$, linear family, $(\alpha_0, \alpha_1) = (0.05, 0.10)$, truth 0.2): DR-DiD on the true $A$ 0.185; naive DR-DiD on $A^*$ 0.153; corrected with oracle nuisances 0.183 ($\tau$) / 0.183 ($\tau_R$); with lasso nuisances 0.183 / 0.182; with boosted nuisances 0.172 / 0.171 (nuisance RMSE for $m_1$: lasso 0.085, boosted 0.143). The correction recovers the true-treatment estimator; the residual with boosting is learner bias. Small-sample check, 20 draws at $N = 1000$: naive 0.163; lasso-corrected 0.182 ($\tau$) / 0.181 ($\tau_R$), empirical SD 0.077, mean SE 0.086; boosted-corrected 0.160 / 0.154. Both learners are registered.
- Lasso dictionary comparison (V13), 35 replicates, $N = 1000$, same draws for every entry (MCSE of bias $\approx 0.015$, of coverage $\approx 0.05$). Bias / empirical SD / mean SE / coverage:

| cell | lasso, raw 27 | lasso, squares (registered `dml_lasso`) | lasso, squares + interactions (`dml_lasso_int`) | `dml_xgb` | `drdid_imp` |
|---|---|---|---|---|---|
| baseline, linear | $-0.026$ / 0.086 / 0.072 / 0.91 | $-0.029$ / 0.084 / 0.072 / 0.91 | $-0.044$ / 0.083 / 0.072 / 0.89 | $-0.046$ / 0.087 / 0.072 / 0.83 | $-0.012$ / 0.090 / 0.074 / 0.89 |
| baseline, nonlinear | $-0.026$ / 0.088 / 0.072 / 0.91 | $-0.029$ / 0.086 / 0.072 / 0.91 | $-0.042$ / 0.084 / 0.072 / 0.89 | $-0.046$ / 0.088 / 0.073 / 0.83 | $-0.011$ / 0.092 / 0.074 / 0.89 |
| `x_classical` $r = 0.75$, nonlinear | $-0.046$ / 0.087 / 0.073 / 0.89 | $-0.050$ / 0.087 / 0.073 / 0.86 | $-0.060$ / 0.086 / 0.072 / 0.83 | $-0.060$ / 0.086 / 0.074 / 0.86 | $-0.030$ / 0.092 / 0.074 / 0.89 |
| `y_nc_covariate` medium, nonlinear | $-0.028$ / 0.101 / 0.086 / 0.91 | $-0.034$ / 0.098 / 0.083 / 0.91 | $-0.040$ / 0.094 / 0.081 / 0.89 | $-0.041$ / 0.106 / 0.084 / 0.86 | $-0.010$ / 0.108 / 0.088 / 0.89 |

  Reading: at $N = 1000$ the three dictionaries are within Monte Carlo error of one another; the ordering (bias grows with the size of the dictionary, toward the boosted learner) is the regularisation component of the learner bias already documented for `dml_xgb`, and the non-linear surface of the DGP is mild enough at this $N$ that the extra terms buy no approximation. The dictionary is registered on specification grounds (it is what the DML literature means by a lasso nuisance learner), not because the pilot favoured it, and the pilot is recorded here so that the production comparison of `dml_lasso` with `dml_lasso_int` is read as pre-specified. Single draws at $N = 5000$, nonlinear baseline (truth 0.199): raw 0.197, squares 0.198, squares + interactions 0.184, boosted 0.183.
- Dry run (2026-09-13, 64 cells $\times$ 2 replicates, 3 jobs in parallel): 1,813 of 1,824 estimator runs clean, and the 11 exceptions are one defect with two faces, both at $N = 1000$ and both from `c20` level 2 (§5). Six `did_ipw` errors where the dummy is empty; five `drdid_imp` glm fallbacks where it holds one or two observations. Rerunning the 11 after the fix gives six finite `did_ipw` estimates and five `drdid_imp` values identical to those recorded. `runs/smoke.R` returns identical estimates, standard errors, intervals and bounds in all four of its cells.
- SIMEX / bounds re-validation after the IF-variance change, 30 replicates, $N = 2000$, linear family, $r = 0.6$, lasso nuisances, $B = 20$: naive bias $-0.012$ (empirical SD 0.053, mean SE 0.051, coverage 0.97); SIMEX bias $+0.003$ (empirical SD 0.061, mean IF SE 0.055, coverage 0.97; Stefanski-Cook SE 0.053); bounds: mean estimated set $[0.185, 0.207]$, Imbens-Manski coverage 0.93 (MCSE 0.05), mean interval width 0.21.
- SIMEX / bounds validation, 50 replicates, $N = 2000$, linear family, $r = 0.6$, $B = 10$ (MCSE of bias $\approx 0.008$, of coverage $\approx 0.04$):

| estimator | bias | empirical SD | mean SE | coverage | set |
|---|---|---|---|---|---|
| oracle lasso (true $X$) | $+0.006$ | 0.058 | 0.051 | 0.92 | |
| naive lasso on $X^*$ | $-0.018$ | 0.057 | 0.051 | 0.92 | |
| naive xgb on $X^*$ | $-0.039$ | 0.056 | 0.051 | 0.92 | |
| covariate-blind OLS | $-0.077$ | 0.055 | 0.051 | 0.68 | |
| SIMEX lasso | $-0.0005$ | 0.065 | 0.053 | 0.88 | |
| SIMEX xgb | $-0.027$ | 0.064 | 0.053 | 0.82 | |
| bounds xgb | (mid) $-0.034$ | 0.053 | | IM 0.96 | $[0.156, 0.175]$ |

  Reading: SIMEX removes the measurement-error component of the bias (lasso: $-0.018 \to -0.0005$) but not the learner's own finite-sample bias, which does not scale with $\lambda$ and is extrapolated as a constant (xgb: $-0.039 \to -0.027$; the naive xgb carries about $-0.02$ of learner bias relative to lasso at $N = 2000$). The Stefanski-Cook SE under-estimates by about 20 % at $B = 10$. Consequences pre-registered here: $B$ raised to 20; the SE replaced by the influence-function estimator (paper Theorem 1(iii)); `bounds_lasso` added alongside `bounds_xgb`; the learner-bias remark is stated in the manuscript after Theorem 1 and the simulation reports the naive-minus-oracle gap per learner as the quantity SIMEX is not designed to remove.

---

## 10. Computational environment and budget

Timing, 2026-09-12 (late, after V13), one replicate, single thread, Apple Silicon (aarch64-apple-darwin20), R 4.5.2, `runs/smoke.R` (seconds). Base-panel column at $N = 1000$ from the `x_classical` $r = 0.75$ linear cell, at $N = 5000$ from the baseline linear cell; correction rows from the `x_classical` and `treat_misclass` linear cells at each $N$.

Two readings were taken at $N = 5000$ for the four covariate-ME corrections, one inside the smoke run and one from a separate script that overlapped with it, so both had a second R process on the machine: 89.9 / 51.6 / 83.6 / 42.2 s and 97.7 / 49.2 / 97.1 / 49.0 s for `simex_lasso` / `simex_xgb` / `bounds_lasso` / `bounds_xgb`. The table and the budget use the slower reading and are therefore mildly conservative. The base-panel rows are from the two cells that ran with nothing else on the machine and agree to within a second (157.3 s and 158.5 s).

| estimator | N = 1000 | N = 5000 |
|---|---|---|
| did_ols, did_ipw, drdid_imp | 0.0 | 0.2 |
| did_ps_rf | 0.4 | 3.0 |
| did_ps_xgb | 0.1 | 0.6 |
| dml_lasso (squares dictionary) | 0.4 | 1.5 |
| dml_lasso_int (squares + interactions) | 2.3 | 5.3 |
| dml_rf | 2.8 | 21.3 |
| dml_xgb | 0.2 | 0.7 |
| cf_did | 4.2 | 37.8 |
| bcf_did | 40.1 | 64.2 |
| match_full | 0.8 | 22.7 |
| **base panel (12)** | **51.3** | **157.3** |
| simex_lasso ($B = 20$) | 25.3 | 97.7 |
| bounds_lasso ($B = 20$) | 26.4 | 97.1 |
| simex_xgb ($B = 20$) | 12.4 | 49.2 |
| bounds_xgb ($B = 20$) | 12.5 | 49.0 |
| **four covariate-ME corrections** | **76.6** | **293.0** |
| misclass_lasso, misclass_lasso_R (each) | 0.6 | 1.8 |
| misclass_xgb, misclass_xgb_R (each) | 0.3 | 1.1 |
| **four misclassification corrections** | **1.7** | **5.8** |

Budget at $n_{sim} = 500$: 32 cells per $N$; corrections in 12 covariate-ME cells and 6 misclassification cells per $N$.
$N = 1000$: $32 \cdot 500 \cdot 51.3 + 12 \cdot 500 \cdot 76.6 + 6 \cdot 500 \cdot 1.7 \approx 1.29 \times 10^6$ s $= 14.9$ CPU-days.
$N = 5000$: $32 \cdot 500 \cdot 157.3 + 12 \cdot 500 \cdot 293.0 + 6 \cdot 500 \cdot 5.8 \approx 4.29 \times 10^6$ s $= 49.7$ CPU-days.
**Total $\approx$ 65 CPU-days ($\approx$ 1,550 core-hours)**, against 45 before V13. Composition: base panel 39 CPU-days (of which BCF 19 and the new `dml_lasso_int` 1.4), the two lasso SIMEX/bounds entries 17, the two boosted ones 9, all misclassification entries 0.3. Seven jobs on the 8-core Mac: 10-14 days wall (single-process timings do not account for contention on the efficiency cores). Oracle Ampere A1, 16 OCPU: about 9 days, about \$35. Always-free 4 OCPU: about 35 days. BCF remains the single most expensive entry and is retained at 250/1000 draws as registered.

Cost lever, identified but not adopted. The lasso penalty is currently chosen by cross-validation separately in each of the five cross-fitting folds, so one estimate costs five cross-validated fits per nuisance. Selecting it once per estimator call and holding it fixed across folds (the DoubleML default when tuning is requested without `tune_on_folds`) would remove most of that work, the more so inside SIMEX, where the learner is refitted 61 times. The saving is not measured and is therefore not claimed as a number; it is recorded here so that, if the budget binds before launch, the change is a pre-stated deviation rather than a post-hoc tweak. It is not the default because per-fold tuning keeps the sample splitting exact, which is what the cross-fitting argument behind the paper's asymptotic-linearity lemma (Lemma 2) assumes.

Driver: `runs/launch.sh` on macOS/Linux, `runs/launch.R` on Windows (flat (cell, batch) jobs, one single-threaded R process each, resumable; tested end-to-end 2026-09-12 with two parallel jobs, and the two drivers verified on 2026-09-13 to write identical `batch_*.rds` apart from wall-clock runtimes). `launch.R` starts its workers with `--no-init-file` and an explicit `R_LIBS` so that renv's autoloader never runs inside a worker, which is the Windows failure of 2026-05-20 (`SIM_ISSUES_2026-09-11.md` I1, I4); it also kills any worker exceeding the per-job timeout and any survivor at exit (I2). Filled at launch: git SHA, `renv.lock` SHA-256, production hardware.

---

## 11. Deviations from v1.3 (all pre-launch; V1-V13 dated 2026-09-12, V14 2026-09-13)

| # | Change | Reason |
|---|---|---|
| V1 | Outcome DGP gains covariate-dependent trends $\delta = \beta/2$ and $n_t$ | v1.3 had trend-free covariates, so unconditional PT held and covariate ME had no bias channel (verified by MC) |
| V2 | `did_ols` on $\Delta Y$ with HC1 SE | v1.3 long-form OLS SE ignored $\alpha_i$; coverage 1.00 (verified) |
| V3 | `me_y_nc_covariate` standardised, non-linear, in $\sigma_\varepsilon$ units | v1.3 used raw age (mean 28), swamping the effect and testing confounding rather than ME |
| V4 | `me_x_nc` scaled per variable | v1.3 added one scalar to variables with SDs from 0.2 to 17 |
| V5 | `me_x_classical` parameterised by reliability | interpretable; needed by SIMEX and bounds |
| V6 | Estimator panel 18 → 12 base + 4 corrected | §4.4; corrections are the new contribution |
| V7 | $N \in \{1000, 5000\}$ | $N = 10000$ is outside the application's range and drove 60 % of cost |
| V8 | Panel strength fixed at "low" | the factor only rescaled $\sigma_\varepsilon$ for $\Delta Y$ estimators |
| V9 | Outcome-ME designs reduced to one proposition-check cell each | analytically settled (P1) |
| V10 | $n_{sim}$ 1000 → 500 | MCSE still $< 0.01$; budget |
| V11 | Runner: flat (cell, batch) jobs, one process each, no callr | post-mortem I1-I4, I12-I15 |
| V12 | No $y_0$ in any adjustment set | consistency across estimators |
| V13 | Lasso learner receives the squares dictionary (39 columns, standardised, per-fold 5-fold CV); `dml_lasso_int` (squares + all pairwise products of the non-dummy columns, 175 columns) added to the base panel; the raw-column lasso of the 2026-09-12 draft is dropped | 2026-09-12 lit review on lasso nuisance learners in DML (Chernozhukov et al. 2018 Sec. 6; Bach et al. 2024 DoubleML 401(k) example; Zimmert 2020; Knaus, Lechner and Strittmatter 2021; Chetverikov, Liao and Chernozhukov 2021; Bodory, Huber and Laffers 2022 use raw covariates only when the DGP is linear in them): a lasso on raw covariates with $p = 27 \ll N$ is a cross-fit logit/OLS, not the linear DML learner the literature means; the degree-2 polynomial dictionary is the standard specification. Full pairwise expansion (370 columns) and the interaction dictionary inside SIMEX are ruled out by cost (§10). Pilot comparison in §9 |

| V14 | `design_matrix()` drops columns with no variation in the sample; `drdid_imp` flags DRDID's `ps.flag = 2` glm fallback as its own replicate outcome; `runs/aggregate.R` computes the primary metrics on clean replicates and repeats them over all usable ones | 2026-09-13 dry run (§9): `c20` level 2 (prevalence 0.003) is empty in 5 % of $N = 1000$ draws, which errors `did_ipw` in `chol()`, and leaves one or two observations often enough that `drdid_imp` reports the standard estimator under the improved estimator's label. §5 |

D3 (callr timeouts) of v1.3 is retired with the estimators that needed it; the shell driver applies a per-job `timeout` instead.

---

## 12. References (additions to v1.3)

- Carroll, R. J., Ruppert, D., Stefanski, L. A. & Crainiceanu, C. M. (2006). *Measurement Error in Nonlinear Models*, 2nd ed., ch. 5 (SIMEX).
- Cook, J. R. & Stefanski, L. A. (1994). Simulation-extrapolation estimation in parametric measurement error models. *JASA* 89, 1314-1328.
- Stefanski, L. A. & Cook, J. R. (1995). Simulation-extrapolation: the measurement error jackknife. *JASA* 90, 1247-1256.
- Imbens, G. W. & Manski, C. F. (2004). Confidence intervals for partially identified parameters. *Econometrica* 72, 1845-1857.
- Stoye, J. (2009). More on confidence intervals for partially identified parameters. *Econometrica* 77, 1299-1315.
- Chang, N.-C. (2020). Double/debiased machine learning for difference-in-differences models. *Econometrics Journal* 23, 177-191.
- Negi, A. & Negi, D. S. (2025). Difference-in-differences with a misclassified treatment. *Journal of Applied Econometrics* 40(4), 411-423.
- Denteh, A. & Kédagni, D. (2022, rev. 2026). Misclassification in difference-in-differences models. arXiv:2207.11890.
- Sant'Anna, P. H. C. & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics* 219, 101-122.
