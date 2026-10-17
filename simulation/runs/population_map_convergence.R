# Grid-convergence check of the population map (14 October 2026, after the map was computed).
# Recomputes selected designs with the latent quadrature, the noise quadrature and the y grid all
# doubled, and compares with population_map.csv. Chosen designs: mild skewness at low reliability,
# the group-variance ratio 2 with a quadratic trend, and the fourth amendment's design as a reference.
#   cd simulation; Rscript runs/population_map_convergence.R
suppressPackageStartupMessages(library(dplyr))
source("R_scripts/44_population_np.R")
map <- read.csv("population_map/population_map.csv")
pts <- data.frame(sl = c(0.25, 0.25, 0.5, 0.25, 0.75, 1), r = c(0.5, 0.8, 0.8, 0.5, 0.5, 0.9),
                  bq = c(0, 0, 0.5, 0.5, 0.5, 0), ratio = c(1, 2, 2, 2, 1, 1))
NP <- local({
  gq <- local({ g <- seq(-8.5, 8.5, length.out = 4001); w <- dnorm(g); list(g = g, w = w / sum(w)) })
  gh <- local({ x <- seq(-7, 7, length.out = 321); w <- dnorm(x); list(x = x, w = w / sum(w)) })
  list(gq = gq, gh = gh)
})
np_ygrid <- function(P, e) {
  w <- c(np_lat(0, P), np_lat(1, P)); ww <- c(NP$gq$w, NP$gq$w); o <- order(w)
  q_hi <- w[o][which(cumsum(ww[o]) / sum(ww) >= 0.9999)[1]]
  lo <- min(w) - 6 * sqrt(e); mid <- q_hi + 6 * sqrt(e); hi <- max(w) + 6 * sqrt(e)
  c(seq(lo, mid, length.out = 2400), if (hi > mid) mid - 1 + exp(seq(0, log(hi - mid + 1), length.out = 601))[-1])
}
cols <- c("naive", "rc", "pm_rc", "simex", "pm_simex")
res <- bind_rows(lapply(seq_len(nrow(pts)), \(i) {
  p <- pts[i, ]; fine <- np_targets(p$sl, p$r, p$bq, p$ratio)[, cols]
  base <- map[map$sl == p$sl & map$r == p$r & map$bq == p$bq & map$ratio == p$ratio, cols]
  cbind(p, setNames(fine, paste0(cols, "_fine")), max_abs_diff = max(abs(unlist(fine) - unlist(base))))
}))
write.csv(res, "population_map/population_map_convergence.csv", row.names = FALSE)
print(res, digits = 4)
