# BIHS child panel for the application (paper1_rc_simex_did/notes/post_prod_sim_steps.md, Section 7).
#   BIHS_DATA_ROOT=<raw root> Rscript application/R_v2/01_build_panel.R      (from repo root)
#
# Reads the raw BIHS rounds read-only. Writes the child-level analysis file to
# application/derived/ (git-ignored: derived from restricted microdata) and
# aggregate, non-disclosive tables to application/output/.
#
# Cohort (fixed before any outcome model was fitted):
#   children aged 6-10 in round 1 (2011-12), currently attending school in 2011,
#   with an education-module record showing no primary stipend (code 5) now or
#   before 2011; present in the same household line in round 2 (2015) and
#   round 3 (2018-19), with consistent sex and plausible ageing.
# Exposure: current receipt of the primary stipend (code 5) in round 2.
# Comparison: no code 5 currently or before 2015 in round 2.
# Outcome: completed years of schooling, round 1 (pre) and round 3 (post).
# Missing education-module rows are "unknown", never "untreated".

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(tidyr)
})
root <- Sys.getenv("BIHS_DATA_ROOT")
if (!dir.exists(root)) cli::cli_abort("BIHS data root {.path {root}} not found; set {.envvar BIHS_DATA_ROOT}.")
out_d <- "application/derived"; out_t <- "application/output"
dir.create(out_d, showWarnings = FALSE, recursive = TRUE)
dir.create(out_t, showWarnings = FALSE, recursive = TRUE)
rd <- \(p) zap_labels(read_dta(file.path(root, p)))
files <- c(r1a = "bihs-r1-og/001_mod_a_male.dta", r1b1 = "bihs-r1-og/003_mod_b1_male.dta",
           r1b2 = "bihs-r1-og/004_mod_b2_male.dta", r1d1 = "bihs-r1-og/006_mod_d1_male.dta",
           r1g = "bihs-r1-og/010_mod_g_male.dta", r2b1 = "bihs-r2-og/003_r2_male_mod_b1.dta",
           r2b2 = "bihs-r2-og/004_r2_mod_b2_male.dta",
           r3b1 = "bihs-r3-og/BIHSRound3/Male/010_bihs_r3_male_mod_b1.dta",
           r3b2 = "bihs-r3-og/BIHSRound3/Male/011_bihs_r3_male_mod_b2.dta")
checksums <- data.frame(module = names(files), file = unname(files),
                        md5 = unname(tools::md5sum(file.path(root, files))))
write.csv(checksums, file.path(out_t, "raw_checksums.csv"), row.names = FALSE)
m <- lapply(files, rd)

# Highest class passed -> completed years. 22/33 are exam candidates (10/12
# completed only on passing, so the preceding level), 66/67 preschool and 99
# never attended are 0; 71-76 (professional/other) are undefined here.
years_of <- function(x) {
  y <- if_else(x %in% 0:16, x, NA_real_)
  y[x %in% c(22)] <- 9; y[x %in% c(33)] <- 11; y[x %in% c(66, 67, 99)] <- 0
  y
}
yes <- \(x) case_when(x == 1 ~ 1L, x == 2 ~ 0L, .default = NA_integer_)
# A household with no adult schooling record gets NA, not the -Inf of an empty max().
max_or_na <- \(y) if (all(is.na(y))) NA_real_ else max(y, na.rm = TRUE)
pid <- \(a01, mid) paste(a01, mid, sep = "_")

# ---- round 1 baseline -----------------------------------------------------------
hh <- m$r1a |> transmute(a01, division = div, district = dcode, cluster = paste(uncode, vcode_n, sep = "_"),
                         ftf_sample = as.integer(hh_type == 1), head_mid = a11)
ros <- m$r1b1 |> mutate(age = b1_02, female = as.integer(b1_01 == 2), yrs = years_of(b1_08))
hhv <- ros |> group_by(a01) |>
  summarise(hh_size = n(), n_child_5_15 = sum(age >= 5 & age <= 15, na.rm = TRUE),
            max_adult_yrs = max_or_na(yrs[age >= 18]), .groups = "drop")
head <- ros |> inner_join(hh |> select(a01, head_mid), by = "a01") |> filter(mid == head_mid) |>
  transmute(a01, head_female = female, head_age = age, head_yrs = yrs)
# Assets: current value of owned items; an all-missing household stays missing.
assets <- m$r1d1 |> filter(d1_03 == 1) |> group_by(a01) |>
  summarise(asset_value = if (all(is.na(d1_10))) NA_real_ else sum(d1_10, na.rm = TRUE), .groups = "drop")
# Land: area (decimals) of plots owned by a member or the household jointly.
land <- m$r1g |> filter(!is.na(g09a), g09a < 70 | g09a == 71) |> group_by(a01) |>
  summarise(land_dec = sum(g02, na.rm = TRUE), .groups = "drop")
b2_1 <- m$r1b2 |> transmute(a01, mid, prog11 = b2_09, past11a = b2_12a, past11b = b2_12b, b2_r1 = 1L)

base <- ros |> filter(age >= 6, age <= 10) |>
  transmute(a01, mid, person_id = pid(a01, mid), age11 = age, female, yrs11 = yrs,
            enrolled11 = yes(b1_09)) |>
  left_join(hh, by = "a01") |> left_join(hhv, by = "a01") |> left_join(head, by = "a01") |>
  left_join(assets, by = "a01") |> left_join(land, by = "a01") |> left_join(b2_1, by = c("a01", "mid")) |>
  mutate(land_dec = coalesce(land_dec, 0),                         # no owned plot recorded
         log_assets = log1p(asset_value), log_land = log1p(land_dec))
if (anyDuplicated(base$person_id)) {
  cli::cli_abort("{sum(duplicated(base$person_id))} duplicated person_id{?s} in the round 1 cohort.")
}

# ---- links to rounds 2 and 3 -------------------------------------------------------
link <- function(b1, b2, mid_b2, suffix) {
  r <- b1 |> filter(mem_stat == 0) |>
    transmute(a01_base = floor(a01), a01_now = a01, mid, female_now = as.integer(b1_01 == 2),
              age_now = b1_02, yrs_now = years_of(b1_08), enrolled_now = yes(b1_09))
  e <- b2 |> transmute(a01_now = a01, mid = .data[[mid_b2]], prog = b2_09, past_a = b2_12a,
                       past_b = b2_12b, b2_present = 1L)
  r <- r |> left_join(e, by = c("a01_now", "mid"))
  # A member must appear once; ambiguous split-household matches are dropped.
  r <- r |> group_by(a01_base, mid) |> filter(n() == 1) |> ungroup()
  names(r)[-(1:3)] <- paste0(names(r)[-(1:3)], suffix)
  r |> select(-a01_now)
}
l2 <- link(m$r2b1, m$r2b2, "mid", "15"); l3 <- link(m$r3b1, m$r3b2, "mid_b2", "18")
panel <- base |>
  left_join(l2, by = c("a01" = "a01_base", "mid")) |>
  left_join(l3, by = c("a01" = "a01_base", "mid"))

has5 <- \(...) Reduce(`|`, lapply(list(...), \(v) !is.na(v) & v == 5))
panel <- panel |> mutate(
  linked15 = !is.na(age_now15), linked18 = !is.na(age_now18),
  consistent = linked15 & linked18 & female_now15 == female & female_now18 == female &
               between(age_now15 - age11, 2, 6) & between(age_now18 - age11, 5, 10),
  hist11 = case_when(is.na(b2_r1) ~ "unknown (no R1 education record)",
                     has5(prog11, past11a, past11b) ~ "stipend by 2011",
                     TRUE ~ "no stipend by 2011"),
  hist15 = case_when(is.na(b2_present15) ~ "unknown (no R2 education record)",
                     !is.na(prog15) & prog15 == 5 ~ "current stipend 2015",
                     has5(past_a15, past_b15) ~ "past stipend only (2011-15)",
                     TRUE ~ "no stipend 2015"))

# ---- sample flow ------------------------------------------------------------------
flow <- list(); step <- function(label, d) { flow[[length(flow) + 1L]] <<- data.frame(step = label,
  children = nrow(d), households = n_distinct(d$a01), clusters = n_distinct(d$cluster)); d }
d <- step("R1 roster, age 6-10", panel)
d <- step("attending school in 2011", filter(d, enrolled11 %in% 1L))
d <- step("R1 education record, no stipend now or before 2011", filter(d, hist11 == "no stipend by 2011"))
d <- step("present in R2 and R3 household lines", filter(d, linked15, linked18))
d <- step("consistent sex and ageing across rounds", filter(d, consistent))
d <- step("R2 exposure known (current stipend or none ever)",
          filter(d, hist15 %in% c("current stipend 2015", "no stipend 2015")))
d <- step("outcome observed 2011 and 2018", filter(d, !is.na(yrs11), !is.na(yrs_now18)))
cov_names <- c("age11", "female", "hh_size", "n_child_5_15", "max_adult_yrs", "head_female",
               "head_age", "head_yrs", "log_assets", "log_land", "ftf_sample", "division")
d <- step("baseline covariates complete", filter(d, if_all(all_of(cov_names), \(v) !is.na(v))))
flow <- do.call(rbind, flow)
analysis <- d |> mutate(A = as.integer(hist15 == "current stipend 2015"),
                        dY = yrs_now18 - yrs11, enrolled18 = enrolled_now18)
write.csv(flow, file.path(out_t, "sample_flow.csv"), row.names = FALSE)

# Exposure-history table, before the exposure restriction.
hist_tab <- panel |> filter(enrolled11 %in% 1L) |> count(hist11, hist15, name = "children")
write.csv(hist_tab, file.path(out_t, "exposure_history.csv"), row.names = FALSE)
# Cross-check for DAG 3: receipt in 2015 against attendance in 2015 (all linked).
rev_tab <- panel |> filter(hist11 == "no stipend by 2011", linked15) |>
  count(hist15, enrolled_now15, name = "children")
write.csv(rev_tab, file.path(out_t, "receipt_by_attendance_2015.csv"), row.names = FALSE)

saveRDS(analysis, file.path(out_d, "analysis_child_panel.rds"))
cat("Sample flow:\n"); print(flow, row.names = FALSE)
cat("\nExposure histories (attending in 2011):\n"); print(hist_tab, row.names = FALSE)
cat("\nReceipt by 2015 attendance:\n"); print(rev_tab, row.names = FALSE)
cat("\nTreated:", sum(analysis$A), " controls:", sum(1 - analysis$A),
    " clusters:", n_distinct(analysis$cluster), "\n")
print(summary(analysis[, c("dY", cov_names[cov_names != "division"])]))
