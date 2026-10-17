# Independent Monte Carlo check of the population map (14 October 2026, after the map was computed).
# The map integrates numerically (R_scripts/44_population_np.R). Here the same targets are computed by
# brute force: draw a very large sample of one design, estimate the control regression E_0[f(W) | W*]
# by averaging f(W) within narrow bins of W*, and average over the treated units. Only the design
# (np_design) is shared with the map; the integration is not. Per coordinate, summed over both.
#   cd simulation; Rscript runs/population_map_mc_check.R
suppressPackageStartupMessages(library(dplyr))
source("R_scripts/44_population_np.R")
map <- read.csv("population_map/population_map.csv")
std_ln <- \(n, sl) if (sl == 0) rnorm(n) else (exp(sl * rnorm(n)) - exp(sl^2 / 2)) / sqrt((exp(sl^2) - 1) * exp(sl^2))

coord_terms <- function(f, P, e, n0 = 6e6, n1 = 4e6, bins = 6000L) {
  w0 <- P$sw0 * std_ln(n0, P$sl); y0 <- w0 + rnorm(n0, 0, sqrt(e))
  w1 <- P$dmu + P$sw1 * std_ln(n1, P$sl); y1 <- w1 + rnorm(n1, 0, sqrt(e))
  br <- unique(quantile(y0, seq(0, 1, length.out = bins + 1), names = FALSE))
  b <- findInterval(y0, br, all.inside = TRUE)
  yc <- tapply(y0, b, mean); mf <- tapply(f(w0), b, mean)                 # binned control regression
  m0 <- \(y) approx(yc, mf, xout = y, rule = 2)$y
  h1 <- np_blp(1, P, e); h0 <- np_blp(0, P, e)
  c(naive = mean(f(w1)) - mean(m0(y1)), rc = mean(f(w1)) - mean(m0(h0$inv(h1$f(y1)))))
}
pts <- data.frame(sl = c(0, 0.25, 0.25, 1), r = c(0.5, 0.5, 0.5, 0.9), bq = c(0, 0, 0.5, 0), ratio = 1)
res <- bind_rows(lapply(seq_len(nrow(pts)), \(i) {
  p <- pts[i, ]; set.seed(7000 + i); D <- np_design(p$sl, p$r, p$bq, p$ratio)
  mc <- coord_terms(D$f1, D$P1, D$s0[1]) + coord_terms(D$f2, D$P2, D$s0[2])
  q <- map[map$sl == p$sl & map$r == p$r & map$bq == p$bq & map$ratio == p$ratio, ]
  data.frame(p, naive_mc = mc[["naive"]], naive_map = q$naive, rc_mc = mc[["rc"]], rc_map = q$rc)
}))
write.csv(res, "population_map/population_map_mc_check.csv", row.names = FALSE)
print(res, digits = 4)
