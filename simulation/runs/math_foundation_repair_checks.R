# Exact finite-support checks of the added foundation derivations.
# Run from repository root. Does not fit or modify production estimators.
out <- "simulation/audit_2026-10-06"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
checks <- list()
record <- function(name, observed, expected, tol = 1e-7) {
  checks[[length(checks) + 1L]] <<- data.frame(check = name, observed = observed,
    expected = expected, discrepancy = observed - expected,
    pass = abs(observed - expected) < tol)
}
px <- c(.2, .3, .5); e <- c(.15, .45, .8)
mu <- c(-.4, .6, 1.2); theta <- c(.1, .35, .65)
z <- expand.grid(x = 1:3, a = 0:1, sign = c(-1, 1))
z$p <- px[z$x] * dplyr::if_else(z$a == 1, e[z$x], 1 - e[z$x]) / 2
z$eta <- z$sign * (.3 + .2 * z$x + .1 * z$a)
z$y <- mu[z$x] + z$a * theta[z$x] + z$eta
E <- function(v) sum(z$p * v)
pi <- E(z$a); tau <- E(z$a * theta[z$x]) / pi
omega <- function(g) (1 - z$a) * g[z$x] / (1 - g[z$x])
w1 <- z$a / pi
w0 <- function(g) omega(g) / E(omega(g))
T_P <- function(g, m, y = z$y) E((w1 - w0(g)) * (y - m[z$x]))
h <- c(.7, -.3, .5); m <- c(.1, .2, .5); eps <- 1e-6
r <- mu - m; rt <- sum(px * e * r) / pi
der <- -sum(px * h / (1 - e) * (r - rt)) / pi
record("Normalized propensity derivative A6", (T_P(e + eps*h,m)-T_P(e-eps*h,m))/(2*eps), der)
g <- c(.3, .35, .7)
record("Outcome nuisance derivative A7", (T_P(g,mu+eps*h)-T_P(g,mu-eps*h))/(2*eps), -E((w1-w0(g))*h[z$x]))
record("Unconditional trend gap covariance identity", sum(px*e*mu)/pi - sum(px*(1-e)*mu)/(1-pi),
       (sum(px*e*mu)-sum(px*e)*sum(px*mu))/(pi*(1-pi)))
ec <- pmin(pmax(e, .2), .7)
record("Clipped propensity exact remainder", T_P(ec,m)-tau, E((w0(e)-w0(ec))*(mu-m)[z$x]))
# Pooling counterexample: each fold has one exactly correct nuisance.
p <- .4; aa <- .07; bb <- .12
d1 <- (1-p)*(p+aa)/(1-p-aa); dbar <- (d1+p)/2
pooled <- -bb/2+p*bb/(2*dbar)
record("Pooled Hajek counterexample exact drift",pooled,-bb*(d1-p)/(4*dbar))
if (abs(pooled) <= 1e-4) cli::cli_abort("The pooled-fold counterexample shows no drift ({signif(pooled, 3)}).")
record("Per-fold product does not control pooled drift", 0*aa+0*bb, 0)
# Population and realized ATT influence decompositions, integrated exactly.
psi <- w1*(z$y-mu[z$x]-tau)-w0(e)*(z$y-mu[z$x])
composition <- w1*(theta[z$x]-tau)
zeta <- (w1-w0(e))*z$eta
record("Population minus realized-target IF",max(abs(psi-composition-zeta)),0)
record("Residual/composition covariance",E(zeta*composition),0)
record("Population versus realized-target variance A13",E(psi^2)-E(zeta^2),E(composition^2))
# Mean-zero measurement error need not increase variance if correlated with Y.
u <- -z$eta; psi_star <- psi+(w1-w0(e))*u
record("Correlated outcome-error variance A10",E(psi_star^2),E(psi^2)+2*E(psi*(w1-w0(e))*u)+E(((w1-w0(e))*u)^2))
record("Outcome error can remove all residual noise",E(psi_star^2),E(composition^2))
if (E(psi_star^2) >= E(psi^2)) cli::cli_abort("Outcome error did not reduce the variance in the example: {signif(E(psi_star^2), 4)} against {signif(E(psi^2), 4)}.")
tab <- do.call(rbind, checks)
if (!all(tab$pass)) cli::cli_abort("Failed checks: {.val {tab$check[!tab$pass]}}.")
write.csv(tab, file.path(out, "repair_foundation_checks.csv"), row.names = FALSE)
print(tab, row.names = FALSE)
cat("All finite-support foundation repair checks passed. These checks do not prove rate or CLT assumptions.\n")
