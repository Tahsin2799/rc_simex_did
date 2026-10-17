# Numerical check of Corollary C (bias of regression calibration with a linear trend, eq. C.23):
#   Q(eta_0) - tau = - sum_c beta_c E_1[ Delta_0c( h_0c^{-1}( h_1c(W*_c) ) ) ],
# against the directly computed target (R_scripts/44_population_np.R), in three log-normal designs.
#   cd simulation; Rscript runs/gap_bias_check.R    -> population_map/gap_bias_check.csv
suppressPackageStartupMessages(library(dplyr))
source("R_scripts/44_population_np.R")
dir.create("population_map", showWarnings = FALSE)
one <- function(sl, r) {
  D <- np_design(sl, r, 0); rep_sum <- 0; direct_sum <- 0
  for (cc in 1:2) {
    P <- if (cc == 1) D$P1 else D$P2; e <- D$s0[cc]; b <- if (cc == 1) 1 else 0.8
    yg <- np_ygrid(P, e); g0 <- np_post(yg, identity, 0, P, e)
    h1 <- np_blp(1, P, e); h0 <- np_blp(0, P, e)
    gap0 <- \(y) approx(yg, g0, xout = y, rule = 2)$y - h0$f(y)
    rep_sum <- rep_sum - b * np_treated_avg(\(y) gap0(h0$inv(h1$f(y))), P, e)
    direct_sum <- direct_sum + unname(np_terms(\(w) b * w, P, e)["rc"])
  }
  data.frame(sl = sl, r = r, representation = rep_sum, direct = direct_sum, abs_diff = abs(rep_sum - direct_sum))
}
res <- rbind(one(1, 0.9), one(0.5, 0.75), one(0.8, 0.6))
write.csv(res, "population_map/gap_bias_check.csv", row.names = FALSE)
print(res, digits = 6)
if (max(res$abs_diff) >= 1e-6) cli::cli_abort("(C.23) and the direct bias disagree: largest gap {signif(max(res$abs_diff), 3)}.")
