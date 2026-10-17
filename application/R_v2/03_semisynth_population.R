# Population for the BIHS-calibrated semi-synthetic study (plan Section 6, P2).
#   Rscript application/R_v2/03_semisynth_population.R      (from repo root; seconds)
# Writes application/derived/semisynth_population.rds (git-ignored: it holds
# covariate rows derived from restricted microdata) and a non-disclosive
# summary to application/output/semisynth_population.csv.
#
# Each replicate (simulation/runs/confirm_simex_fold.R, study "bihs_semi")
# resamples whole villages from the analytic cohort, keeps their real joint
# covariate distribution, and then simulates:
#   A  ~ Bern(logit^{-1}(X b_A)),  b_A from a logit of the real A on X;
#   dY = X b_Y + nl * z_assets * z_land + tau A + u_village + e,
#        b_Y from OLS of the real dY on X among controls, residual variance
#        split 20 % village / 80 % child, tau = 0.25 school years;
#   W* = W + U for log assets and log land, U Gaussian with fixed variance
#        (1/r0 - 1) Var(W), r0 = 0.7.
# The observed BIHS covariates play the role of the error-free truth here.

suppressPackageStartupMessages({ library(dplyr) })
an <- readRDS("application/derived/analysis_child_panel.rds")
COVS <- c("age11", "female", "hh_size", "n_child_5_15", "max_adult_yrs", "head_female",
          "head_age", "head_yrs", "log_assets", "log_land")
X <- cbind(as.matrix(an[, COVS]), model.matrix(~ factor(division), an)[, -1, drop = FALSE])
colnames(X) <- make.names(colnames(X))
X <- X[, apply(X, 2, \(v) any(v != v[1])), drop = FALSE]
ME <- c("log_assets", "log_land")
fitA <- glm(an$A ~ X, family = binomial())
ctrl <- an$A == 0
fitY <- lm(an$dY[ctrl] ~ X[ctrl, ])
sigma <- summary(fitY)$sigma
mu <- colMeans(X[, ME]); sdv <- apply(X[, ME], 2, sd)
r0 <- 0.7
pop <- list(X = X, cluster = an$cluster, bA = coef(fitA), bY = coef(fitY), sigma = sigma,
            icc = 0.2, nl = 0.25 * sigma, tau = 0.25, me = ME, r0 = r0,
            mu_me = mu, sd_me = sdv, s0 = (1 / r0 - 1) * apply(X[, ME], 2, var))
pop$bA[is.na(pop$bA)] <- 0; pop$bY[is.na(pop$bY)] <- 0
e <- plogis(drop(cbind(1, X) %*% pop$bA))
saveRDS(pop, "application/derived/semisynth_population.rds")
summ <- data.frame(item = c("children", "villages", "covariate columns", "tau (years)",
                            "outcome residual sd", "village share of residual variance",
                            "nonlinear coefficient", "reliability r0 (log assets, log land)",
                            "propensity p01", "propensity median", "propensity p99", "mean propensity"),
                   value = c(nrow(X), length(unique(an$cluster)), ncol(X), pop$tau, sigma, pop$icc,
                             pop$nl, r0, quantile(e, .01), median(e), quantile(e, .99), mean(e)))
write.csv(summ, "application/output/semisynth_population.csv", row.names = FALSE)
print(summ, digits = 4, row.names = FALSE)
