# Simulation grid, design v2.0 (PreReg v2.0 s2.6-2.7).
#   source("runs/grid.R")
#   N in {1000, 5000} x DGP family {linear, nonlinear} x 16 ME cells = 64 cells
#   panel strength fixed at "low" (v1.3's second regime only rescaled
#   sigma_eps for dY-based estimators; dropped, see PreReg v2.0 s11).

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})

me_grid <- function() {
  tribble(
    ~me_design,          ~mag,     ~fn_name,               ~args,
    "baseline",          "none",   NA_character_,          list(),
    # covariate ME (the paper's focus); reliability r = Var(X)/Var(X*)
    "x_classical",       "r0.90",  "me_x_classical",       list(reliability = 0.90),
    "x_classical",       "r0.75",  "me_x_classical",       list(reliability = 0.75),
    "x_classical",       "r0.60",  "me_x_classical",       list(reliability = 0.60),
    "x_nc",              "small",  "me_x_nc",              list(delta_a = 0.10, sigma = 0.10),
    "x_nc",              "medium", "me_x_nc",              list(delta_a = 0.30, sigma = 0.10),
    "x_nc",              "large",  "me_x_nc",              list(delta_a = 0.50, sigma = 0.10),
    # treatment misclassification (Lewbel; Negi & Negi 2025; Denteh-Kedagni)
    "treat_misclass",    "small",  "me_treat_misclass",    list(alpha0 = 0.02, alpha1 = 0.02),
    "treat_misclass",    "medium", "me_treat_misclass",    list(alpha0 = 0.05, alpha1 = 0.10),
    "treat_misclass",    "large",  "me_treat_misclass",    list(alpha0 = 0.10, alpha1 = 0.20),
    # non-linear covariate-dependent outcome ME (rho in units of sigma_eps)
    "y_nc_covariate",    "small",  "me_y_nc_covariate",    list(rho = 0.25, sigma = 0.35),
    "y_nc_covariate",    "medium", "me_y_nc_covariate",    list(rho = 0.50, sigma = 0.35),
    "y_nc_covariate",    "large",  "me_y_nc_covariate",    list(rho = 1.00, sigma = 0.35),
    # proposition-check cells (analytically settled; one magnitude each)
    "y_classical",       "medium", "me_y_classical",       list(sigma_eps_y = 0.70),
    "y_diff",            "medium", "me_y_diff",            list(sigma0 = 0.70, sigma1 = 1.75),
    "y_nc_contraction",  "medium", "me_y_nc_contraction",  list(gamma = -0.30, sigma = 0.35),
  )
}

build_grid <- function(N = c(1000L, 5000L), dgp_family = c("linear", "nonlinear")) {
  g <- crossing(expand_grid(N = N, dgp_family = dgp_family), me_grid()) |>
    mutate(panel_strength = "low",
           cell_id = paste(me_design, mag, dgp_family, sprintf("N%d", N), sep = "__")) |>
    select(cell_id, N, dgp_family, panel_strength, me_design, mag, fn_name, args)
  dup <- unique(g$cell_id[duplicated(g$cell_id)])
  if (length(dup)) cli::cli_abort(c(
    "Duplicate cell id{?s} in the grid: {.val {dup}}.",
    i = "Two rows of `me_grid()` share a (me_design, mag) pair, so their batch files would overwrite each other."
  ))
  g
}

# Flat job table: one row per (cell, batch). n_reps split into batches of
# `batch_size` reps; job k covers reps ((k-1)*batch_size + 1) : (k*batch_size).
build_jobs <- function(grid = build_grid(), n_reps = 500L, batch_size = 25L) {
  nb <- ceiling(n_reps / batch_size)
  crossing(grid |> select(cell_id, N), batch = seq_len(nb)) |>
    mutate(rep_from = (batch - 1L) * batch_size + 1L,
           rep_to   = pmin(batch * batch_size, n_reps)) |>
    arrange(N, cell_id, batch)   # cheap N first so early failures surface fast
}
