# ── Stage 03: cohort assembly ─────────────────────────────────────────────────
# In  : 5_Data/1_Interim/*.rds  (stage 01-02)
# Out : 5_Data/2_Processed/*.rds
#       6_Results/2_Cohort/{Cohort_Filter_Flow, Cohort_Counts,
#                           Cohort_ParameterSpecific_Counts}.csv
#
# Assembly order:
#   1. maximal-effort CPETs (RER >= 1.0), split into HCM and non-HCM registries
#   2. drop visits with CAD / COPD / ILD
#   3. pair each CPET with the nearest echo within +/- 7 days
#   4. phenotypic HCM entry (wall >= 1.5 cm, apical, or gene+ with >= 1.3 cm)
#   5. baseline visit per patient, then analytic eligibility
#   6. non-HCM reference screening
#   7. longitudinal subcohort (repeat CPETs)
#   8. outcomes, medications, comorbidity flags
#   9. drop post-myectomy CPETs from the longitudinal frame
#  10. counts for the manuscript
#
# The cohort definitions themselves live in R/cohort.R and are called from here
# and nowhere else, so study entry and eligibility have exactly one definition.

source("R/00_setup.R")
stage_banner("Stage 03: cohort assembly")

analytic_cohort_mode <- cfg$cohort$analytic_cohort_mode
filter_labels <- analytic_filter_labels(analytic_cohort_mode)

cpx_derived <- load_data("cpx_derived")

# ── 1. Maximal-effort registry ────────────────────────────────────────────────
CPX_effort_all   <- prepare_cpx_registry(cpx_derived, rer_min = cfg$cohort$rer_min)
CPX_hcm_effort   <- CPX_effort_all %>% filter(MCl3 == "HCM")
CPX_nonhcm_effort <- CPX_effort_all %>% filter(!has_any_hcm_flag)

CPX_clean        <- format_cpx_registry(CPX_hcm_effort, keep_diagnosis = FALSE)
CPX_nonhcm_clean <- format_cpx_registry(CPX_nonhcm_effort, keep_diagnosis = TRUE)

filter_flow <- tibble(
  stage = c("Initial CPX import", "Assigned stable patient ID", "HCM diagnosis + RER >= 1.0"),
  n_tests = c(nrow(cpx_derived), nrow(cpx_derived), nrow(CPX_hcm_effort)),
  n_patients = c(n_distinct(cpx_derived$MRN), n_distinct(cpx_derived$ID),
                 n_distinct(CPX_hcm_effort$ID))
)

# ── 2. Comorbidity exclusions, at the visit level ─────────────────────────────
Comorbidities_table     <- load_data("comorbidities_raw")
Comorbidities_to_filter <- comorbidity_clean_visits(Comorbidities_table)

CPX_comorbidity_filtered        <- apply_comorbidity_exclusions(CPX_clean, Comorbidities_to_filter)
CPX_nonhcm_comorbidity_filtered <- apply_comorbidity_exclusions(CPX_nonhcm_clean, Comorbidities_to_filter)

filter_flow <- bind_rows(filter_flow, tibble(
  stage = "After excluding CAD/COPD/ILD",
  n_tests = nrow(CPX_comorbidity_filtered),
  n_patients = n_distinct(CPX_comorbidity_filtered$ID)
))

# ── 3. Echo alignment ─────────────────────────────────────────────────────────
echo_all <- prepare_echo_measures(load_data("echo_raw"))

cpx_echo        <- align_cpx_to_echo(CPX_comorbidity_filtered, echo_all, cfg$cohort$echo_window_days)
cpx_nonhcm_echo <- align_cpx_to_echo(CPX_nonhcm_comorbidity_filtered, echo_all, cfg$cohort$echo_window_days)

filter_flow <- bind_rows(filter_flow, tibble(
  stage = "Echo aligned (within +/- 7d)",
  n_tests = nrow(cpx_echo),
  n_patients = n_distinct(cpx_echo$ID)
))

# ── 4. Phenotypic HCM entry ───────────────────────────────────────────────────
genetics_flags      <- derive_genetics_flags(load_data("genetics_raw"))
cpx_echo_pre_entry  <- derive_hcm_pre_entry(cpx_echo, genetics_flags)

# Anyone who reached the HCM registry can never serve as a reference subject.
hcm_registry_mrns <- cpx_echo_pre_entry %>% distinct(MRN)
cpx_nonhcm_echo   <- cpx_nonhcm_echo %>% anti_join(hcm_registry_mrns, by = "MRN")

cpx_echo <- apply_hcm_entry_window(cpx_echo_pre_entry)

filter_flow <- bind_rows(filter_flow, tibble(
  stage = "HCM inclusion: phenotypic criteria (wall ≥1.5 cm, apical, or gene+ ≥1.3 cm)",
  n_tests = nrow(cpx_echo),
  n_patients = n_distinct(cpx_echo$ID)
))

# ── 5. Baseline visit and analytic eligibility ────────────────────────────────
# Fields coerced to numeric once here, so no downstream model has to repeat it.
num_vars <- c(
  "VO2_FRIEND2_PP", "VO2_FRIEND_PP", "VO2_WASSERMAN_PP", "VeVco2_slope",
  "pk.RER", "age", "BMI", "LBMI", "e_e_ave", "la_vol_index", "tr_max_vel",
  "med_peak_e_vel", "lat_peak_e_vel", "ef_modsp4", "lvot_max_gradient",
  "lv_septal_thickness", "mv_e_a", "mv_dec_time", "E_vel", "e_prime_ave"
)

# cross_sectional_eligible = both E/e' and LAVi measured (the 261-patient
# consistency cohort). analysis_base_eligible = the 441-patient base cohort that
# the primary parameter-specific analyses draw from.
baseline_df <- cpx_echo %>%
  derive_baseline_visits(num_vars) %>%
  add_cross_sectional_eligible(analytic_cohort_mode) %>%
  add_parameter_eligibility()

# The common cohort is exactly the intersection of the two parameters it
# requires, so the two definitions cannot drift apart.
if (analytic_cohort_mode == "ee_lavi") {
  stopifnot(identical(
    baseline_df$cross_sectional_eligible,
    baseline_df$eligible_e_e_ave & baseline_df$eligible_la_vol_index
  ))
}

filter_flow <- bind_rows(filter_flow, tibble(
  stage = filter_labels$stage,
  n_tests = NA_integer_,
  n_patients = sum(baseline_df$cross_sectional_eligible, na.rm = TRUE)
))

# ── 6. Non-HCM reference screening ────────────────────────────────────────────
baseline_nonhcm_df <- derive_nonhcm_baseline(cpx_nonhcm_echo, num_vars)

# ── 7. Longitudinal subcohort ─────────────────────────────────────────────────
# Spans every base-eligible patient; each longitudinal model then restricts to
# the patients with its own baseline parameter measured, exactly as the
# cross-sectional models do.
base_eligible_ids   <- baseline_df %>% filter(analysis_base_eligible) %>% distinct(ID)
cross_sectional_ids <- baseline_df %>% filter(cross_sectional_eligible) %>% distinct(ID)

longitudinal    <- build_longitudinal_subcohort(cpx_echo, base_eligible_ids, num_vars)
longitudinal_ids <- longitudinal$ids
longitudinal_cpx <- longitudinal$cpx
long_df          <- longitudinal$frame

# ── 8. Outcomes, medications, comorbidity flags ───────────────────────────────
outcomes        <- prepare_outcomes_table(load_data("outcomes_raw"))
outcomes_index  <- summarise_outcome_index(outcomes)
outcomes_dedup  <- summarise_outcome_dates(outcomes)

baseline_df <- add_outcome_endpoints(baseline_df, outcomes_dedup, outcomes_index)

medications <- derive_medication_flags(load_data("medications_raw"))
baseline_df        <- left_join_preserving_rows(baseline_df, medications, by = c("MRN", "cpx_test_date"))
long_df            <- left_join_preserving_rows(long_df, medications, by = c("MRN", "cpx_test_date"))
baseline_nonhcm_df <- left_join_preserving_rows(baseline_nonhcm_df, medications, by = c("MRN", "cpx_test_date"))

# ── 9. Post-myectomy CPETs leave the longitudinal frame ───────────────────────
n_long_pre_censor     <- nrow(long_df)
n_long_pts_pre_censor <- n_distinct(long_df$ID)
long_df <- censor_longitudinal_post_myectomy(long_df, outcomes_dedup)

message(sprintf("Post-myectomy CPETs dropped: %d obs (%d -> %d), %d -> %d patients",
                n_long_pre_censor - nrow(long_df), n_long_pre_censor, nrow(long_df),
                n_long_pts_pre_censor, n_distinct(long_df$ID)))

comorbidity_adjustment_flags <- derive_comorbidity_adjustment_flags(Comorbidities_table)
join_keys <- c("MRN", "cpx_test_date")
baseline_df        <- left_join_preserving_rows(baseline_df, comorbidity_adjustment_flags, join_keys)
long_df            <- left_join_preserving_rows(long_df, comorbidity_adjustment_flags, join_keys)
baseline_nonhcm_df <- left_join_preserving_rows(baseline_nonhcm_df, comorbidity_adjustment_flags, join_keys)

# ── 10. Counts reported in the manuscript ─────────────────────────────────────
long_summary <- long_df %>%
  group_by(ID) %>%
  summarise(
    n_tests = n(),
    max_fu_yrs = max(time_yrs, na.rm = TRUE),
    inter_test_days = if (n() > 1) median(diff(sort(days_from_baseline)), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

analytic_gradient   <- as_num(baseline_df$lvot_max_gradient[baseline_df$cross_sectional_eligible])
n_gradient_analytic <- sum(!is.na(analytic_gradient))
n_obst_analytic     <- sum(analytic_gradient >= 30, na.rm = TRUE)

cross_sectional_n <- sum(baseline_df$cross_sectional_eligible, na.rm = TRUE)
base_analytic_n   <- sum(baseline_df$analysis_base_eligible, na.rm = TRUE)
hf_cohort_n       <- sum(baseline_df$outcomes_analytic, na.rm = TRUE)
hf_events_n       <- sum(baseline_df$hf_composite == 1 & baseline_df$outcomes_analytic, na.rm = TRUE)

# Structural invariants of the nested cohorts (not expected counts, which go
# stale; these say only that the nesting still holds).
stopifnot(
  hf_cohort_n <= cross_sectional_n,
  cross_sectional_n <= nrow(baseline_df),
  hf_events_n <= hf_cohort_n,
  all(long_df$ID %in% base_eligible_ids$ID),
  !any(baseline_nonhcm_df$MRN %in% hcm_registry_mrns$MRN)
)

q <- function(x, p) unname(stats::quantile(x, p, na.rm = TRUE))

cohort_counts <- tibble::tribble(
  ~quantity, ~value,
  "Screened HCM cohort (patients)",                      nrow(baseline_df),
  "Base analytic cohort (patients)",                     base_analytic_n,
  "Base analytic cohort with follow-up (patients)",
    sum(baseline_df$analysis_base_eligible & baseline_df$outcomes_followup_ok, na.rm = TRUE),
  "Base analytic cohort composite events",
    sum(baseline_df$analysis_base_eligible & baseline_df$outcomes_followup_ok &
          baseline_df$hf_composite == 1, na.rm = TRUE),
  "Common E/e' + LAVi cohort (patients)",                cross_sectional_n,
  "Common-cohort outcomes subset (patients)",            hf_cohort_n,
  "Common-cohort outcomes subset composite events",      hf_events_n,
  "Analytic cohort with measured resting LVOT gradient", n_gradient_analytic,
  "Obstructive HCM (resting LVOT >= 30 mm Hg)",          n_obst_analytic,
  "Non-obstructive HCM (resting LVOT < 30 mm Hg)",       n_gradient_analytic - n_obst_analytic,
  "Resting LVOT gradient unavailable",                   cross_sectional_n - n_gradient_analytic,
  "Longitudinal patients",                               nrow(long_summary),
  "Longitudinal observations",                           nrow(long_df),
  "Longitudinal observations before myectomy censoring", n_long_pre_censor,
  "Longitudinal patients before myectomy censoring",     n_long_pts_pre_censor,
  "Longitudinal follow-up, median years",                median(long_summary$max_fu_yrs, na.rm = TRUE),
  "Longitudinal follow-up, Q1 years",                    q(long_summary$max_fu_yrs, 0.25),
  "Longitudinal follow-up, Q3 years",                    q(long_summary$max_fu_yrs, 0.75),
  "CPETs per patient, median",                           median(long_summary$n_tests, na.rm = TRUE),
  "CPETs per patient, Q1",                               q(long_summary$n_tests, 0.25),
  "CPETs per patient, Q3",                               q(long_summary$n_tests, 0.75),
  "Inter-test interval, median days",                    median(long_summary$inter_test_days, na.rm = TRUE),
  "Inter-test interval, Q1 days",                        q(long_summary$inter_test_days, 0.25),
  "Inter-test interval, Q3 days",                        q(long_summary$inter_test_days, 0.75)
)

# What each parameter's own available-case sample contains. These are the
# denominators the manuscript reports.
parameter_specific_counts <- summarise_parameter_specific_counts(
  baseline_df, cpx_echo, outcomes_dedup, num_vars
)

write_result(cohort_counts, "2_Cohort", "Cohort_Counts.csv")
write_result(filter_flow, "2_Cohort", "Cohort_Filter_Flow.csv")
write_result(parameter_specific_counts, "2_Cohort", "Cohort_ParameterSpecific_Counts.csv")

print(as.data.frame(cohort_counts), row.names = FALSE)
print(as.data.frame(parameter_specific_counts), row.names = FALSE)

# ── ASE threshold crossings ───────────────────────────────────────────────────
baseline_df <- add_ase_grading(
  baseline_df,
  eprime_bands = cfg$eprime_reference$bands,
  min_classifiable_age = cfg$eprime_reference$min_age
)
long_df <- left_join_preserving_rows(long_df, baseline_df %>% select(ID, fp_class), by = "ID")

# ── Save for the downstream stages ────────────────────────────────────────────
for (name in c(
  "echo_all", "cpx_echo", "genetics_flags", "Comorbidities_to_filter",
  "baseline_df", "baseline_nonhcm_df", "cross_sectional_ids",
  "longitudinal_ids", "longitudinal_cpx", "long_df", "long_summary",
  "outcomes", "outcomes_index", "outcomes_dedup", "medications",
  "comorbidity_adjustment_flags", "filter_flow"
)) {
  save_data(get(name), name, where = "processed")
}

save_data(
  list(
    num_vars = num_vars,
    analytic_cohort_mode = analytic_cohort_mode,
    filter_labels = filter_labels,
    n_long_pre_censor = n_long_pre_censor,
    n_long_pts_pre_censor = n_long_pts_pre_censor
  ),
  "cohort_parameters",
  where = "processed"
)

message("Stage 03 complete.")
