# The paper's population map (Section 5): targets of regression calibration, SIMEX and posterior-mean
# calibration (alone and with SIMEX) with nonparametric nuisances, over latent shape x trend x reliability x
# group-variance ratio. Deterministic; no estimator is fitted. Computed on 14 October 2026, after the
# simulation studies, as one systematic map in place of the sequence of designs used so far.
#   cd simulation; Rscript runs/population_map.R
# Writes population_map/population_map.csv.
suppressPackageStartupMessages(library(dplyr))
source("R_scripts/44_population_np.R")
OUT <- "population_map"; dir.create(OUT, showWarnings = FALSE)

# Check against values verified earlier (closed form and converged grids).
chk <- rbind(np_targets(1, 0.9, 0), np_targets(0, 0.6, 0.5, ratio = 2), np_targets(0.8, 0.6, 0.5))
ref <- data.frame(rc = c(0.14378, 0.03410, -0.21910), simex = c(0.08337, 0.34888, 0.26200))
ref_gap <- max(abs(c(chk$rc - ref$rc, chk$simex - ref$simex)))
if (!isTRUE(ref_gap < 5e-4)) cli::cli_abort("Population map does not reproduce the verified reference values: largest gap {signif(ref_gap, 3)}.")
cat("check passed\n")

grid <- expand.grid(sl = c(0, 0.25, 0.5, 0.75, 1), bq = c(0, 0.5), r = c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95), ratio = c(1, 2))
t0 <- Sys.time()
map <- bind_rows(lapply(seq_len(nrow(grid)), \(i) do.call(np_targets, as.list(grid[i, ]))))
map <- map |> mutate(skewness = if_else(sl == 0, 0, (exp(sl^2) + 2) * sqrt(exp(sl^2) - 1)),
                     better = if_else(abs(rc) <= abs(simex), "RC", "SIMEX"))
write.csv(map, file.path(OUT, "population_map.csv"), row.names = FALSE)
cat(sprintf("%d designs in %.1f min\n", nrow(map), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
options(width = 200); print(as.data.frame(map), digits = 3)
