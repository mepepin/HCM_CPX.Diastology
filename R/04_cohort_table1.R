# ── Stage 04: reference comparison, cohort flow, and Table 1 ──────────────────
# In  : 5_Data/2_Processed/{baseline_df, baseline_nonhcm_df, filter_flow,
#                           cohort_parameters}.rds
# Out : 6_Results/2_Cohort/
#         Figure1_CohortFlow.pdf                  (main Figure 1)
#         Table1_Manuscript.docx                  (main Table 1)
#         FigureS_PropensityOverlap.pdf           (Figure S1)
#         Table_Manuscript_Descriptives.csv       (in-text numbers)
#         Table_Comparison_Balance_OverlapWeights.csv
#         Table_Comparison_Covariates_WeightedUnweighted.csv
#         Comparison_Estimand_Statement.md        (pasted into the Methods)
#
# The HCM cohort is compared with the non-HCM reference group using
# propensity-score overlap weights on age, sex, and BMI, which keeps every
# eligible patient rather than discarding unmatched comparators.

source("R/00_setup.R")
stage_banner("Stage 04: reference comparison, cohort flow, Table 1")

OUT_DIR <- results_dir("2_Cohort")

baseline_df <- load_data("baseline_df", "processed")
baseline_nonhcm_df <- load_data("baseline_nonhcm_df", "processed")
long_df <- load_data("long_df", "processed")
filter_flow <- load_data("filter_flow", "processed")
cohort_parameters <- load_data("cohort_parameters", "processed")
analytic_cohort_mode <- cohort_parameters$analytic_cohort_mode
analytic_filter_stage <- cohort_parameters$filter_labels$stage
analytic_filter_exclusion <- cohort_parameters$filter_labels$exclusion
analytic_filter_description <- cohort_parameters$filter_labels$description
pdftoppm <- Sys.which("pdftoppm")

# ══ Figure 1: Controls ════════════════════════════════════════════════════════
# (legacy chunk: figure1_matched_controls)
# ── 1:1 matched controls (age/sex/BMI) ────────────────────────────────────────
# Outcome-blind candidate grid with exact sex and age/BMI balance. Select the
# largest match for which every prespecified post-match SMD is <0.10.
# Table 1 and the reference comparison describe the whole base analytic
# cohort. A measured resting LVOT gradient is needed only to label obstruction,
# not to enter the comparison, so patients without one are retained and appear
# in Table 1 as a separate group.
hcm_intro_df <- baseline_df %>%
  filter(analysis_base_eligible) %>%
  mutate(
    hcm_primary_binary = factor(
      ifelse(
        hcm_combo_n_abnormal == 0,
        "HCM (Normal diastology)",
        "HCM (Abnormal diastology)"
      ),
      levels = c("HCM (Normal diastology)", "HCM (Abnormal diastology)")
    )
  )

nonhcm_control_pool <- baseline_nonhcm_df %>%
  filter(control_qc_eligible, !is.na(age), !is.na(Sex), !is.na(BMI))

stopifnot(
  nrow(nonhcm_control_pool) > 0,
  all(!nonhcm_control_pool$has_any_hcm_flag),
  all(nonhcm_control_pool$control_no_exclusionary_diagnosis),
  all(nonhcm_control_pool$ef_modsp4 >= 50),
  all(nonhcm_control_pool$lv_max_wall_thickness < 1.3),
  # An unmeasured resting gradient is permitted (see derive_nonhcm_baseline);
  # no measured gradient may reach the obstruction threshold.
  all(nonhcm_control_pool$lvot_max_gradient < 30, na.rm = TRUE),
  !any(is.na(nonhcm_control_pool$ef_modsp4)),
  !any(is.na(nonhcm_control_pool$lv_max_wall_thickness))
)

# Primary comparison: propensity-score overlap weighting, which keeps every
# eligible patient instead of discarding unmatched comparators. See
# 8_Docs/ANALYSIS_PLAN_2026-09.md.
comparison_weights <- overlap_weight_comparison(hcm_intro_df, nonhcm_control_pool)
comparison_hcm_df <- comparison_weights$cases
comparison_nonhcm_df <- comparison_weights$controls

write.csv(
  comparison_weights$balance %>%
    mutate(
      n_hcm = comparison_weights$n_cases,
      n_reference = comparison_weights$n_controls,
      effective_n_hcm = round(comparison_weights$ess_cases, 1),
      effective_n_reference = round(comparison_weights$ess_controls, 1)
    ),
  file.path(OUT_DIR, "Table_Comparison_Balance_OverlapWeights.csv"),
  row.names = FALSE
)

# Weighting balances the covariates in the propensity model by construction;
# assert it rather than assume it, and confirm nobody was dropped.
stopifnot(
  max(abs(comparison_weights$balance$`SMD (overlap-weighted)`)) < 0.01,
  comparison_weights$n_cases == sum(stats::complete.cases(
    hcm_intro_df[, c("age", "BMI")]
  ) & hcm_intro_df$Sex %in% c("Male", "Female")),
  comparison_weights$ess_cases <= comparison_weights$n_cases,
  comparison_weights$ess_controls <= comparison_weights$n_controls
)
message(sprintf(
  "Reference comparison: %d HCM and %d non-HCM references retained (effective n %.0f and %.0f); max |SMD| %.4f after weighting (%.3f before).",
  comparison_weights$n_cases, comparison_weights$n_controls,
  comparison_weights$ess_cases, comparison_weights$ess_controls,
  max(abs(comparison_weights$balance$`SMD (overlap-weighted)`)),
  max(abs(comparison_weights$balance$`SMD (unweighted)`))
))

comparison_hcm_df <- comparison_hcm_df %>%
  mutate(
    lvot_max_gradient_num = as_num(lvot_max_gradient),
    obstructive_hcm = lvot_max_gradient_num >= 30,
    figure1_hist_group = case_when(
      is.na(lvot_max_gradient_num) ~ "HCM, gradient not measured",
      obstructive_hcm ~ "Obstructive HCM",
      TRUE ~ "Non-obstructive HCM"
    )
  )

comparison_nonhcm_df <- comparison_nonhcm_df %>%
  mutate(figure1_hist_group = "CON")

matched_hcm_n <- nrow(comparison_hcm_df)
matched_nonhcm_n <- nrow(comparison_nonhcm_df)
hcm_normal_n <- sum(
  comparison_hcm_df$hcm_primary_binary == "HCM (Normal diastology)",
  na.rm = TRUE
)
hcm_abnormal_n <- sum(
  comparison_hcm_df$hcm_primary_binary == "HCM (Abnormal diastology)",
  na.rm = TRUE
)
n_hcm_obstructive <- sum(comparison_hcm_df$obstructive_hcm, na.rm = TRUE)
n_hcm_nonobstructive <- sum(!comparison_hcm_df$obstructive_hcm, na.rm = TRUE)

case_control_hist_df <- bind_rows(
  comparison_hcm_df %>%
    filter(figure1_hist_group != "HCM, gradient not measured") %>%
    mutate(Cohort = figure1_hist_group),
  comparison_nonhcm_df %>%
    mutate(
      Cohort = figure1_hist_group,
      HCM_Phenotype = NA_character_,
      hcm_primary_binary = NA
    )
) %>%
  mutate(Cohort = factor(
    Cohort,
    levels = c("CON", "Non-obstructive HCM", "Obstructive HCM")
  ))

eprime_reference_df <- comparison_nonhcm_df %>%
  filter(!is.na(age), !is.na(e_prime_ave)) %>%
  transmute(age = as.numeric(age), e_prime_ave = as.numeric(e_prime_ave))

if (nrow(eprime_reference_df) >= 10) {
  eprime_age_model <- lm(
    e_prime_ave ~ ns(age, df = 3),
    data = eprime_reference_df
  )
  case_control_hist_df <- case_control_hist_df %>%
    mutate(
      avg_eprime_age_expected = ifelse(
        !is.na(age),
        predict(eprime_age_model, newdata = tibble(age = as.numeric(age))),
        NA_real_
      ),
      avg_eprime_age_cal = e_prime_ave - avg_eprime_age_expected
    )
} else {
  case_control_hist_df <- case_control_hist_df %>%
    mutate(avg_eprime_age_expected = NA_real_, avg_eprime_age_cal = NA_real_)
}

# ══ Manuscript descriptives on the base analytic cohort ══════════════════════
# The manuscript reports cohort descriptives that no other table carried in a
# consumable form, so the Results text had to be written against numbers a
# reader could not trace. They are generated here instead, as quantity/value
# pairs, and the manuscript quotes this file.
fmt_med_iqr <- function(x, digits = 1) {
  x <- as_num(x)
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  sprintf(
    paste0("%.", digits, "f (IQR %.", digits, "f to %.", digits, "f)"),
    q[2], q[1], q[3]
  )
}
fmt_n_pct <- function(n, total) sprintf("%d (%.1f%%)", n, 100 * n / total)

# Exercise modality is recorded in the workbook as "Treadmill.(0).or.cycle.(1)"
# and is consumed by the FRIEND/Wasserman derivations, but it was dropped from
# the analysis frames, so the Methods had no traceable source for the test
# counts. Joined back here for reporting only.
cpx_modality <- load_data("cpx_derived", "interim") %>%
  transmute(
    MRN = as.character(MRN),
    cpx_test_date = as.Date(cpx_test_date),
    cpx_mode = as_num(.data[["Treadmill.(0).or.cycle.(1)"]])
  ) %>%
  distinct(MRN, cpx_test_date, .keep_all = TRUE)

descr_df <- hcm_intro_df %>%
  mutate(
    .mrn_chr = as.character(MRN),
    .test_date = as.Date(cpx_test_date)
  ) %>%
  left_join(
    cpx_modality,
    by = c(".mrn_chr" = "MRN", ".test_date" = "cpx_test_date")
  )
stopifnot(nrow(descr_df) == nrow(hcm_intro_df))

n_base <- nrow(descr_df)
gradient_recorded <- sum(!is.na(as_num(descr_df$lvot_max_gradient)))
obstructive_n <- sum(as_num(descr_df$lvot_max_gradient) >= 30, na.rm = TRUE)
ee_measured <- sum(!is.na(descr_df$e_e_ave))
lavi_measured <- sum(!is.na(descr_df$la_vol_index))
trv_measured <- sum(!is.na(descr_df$tr_max_vel))
both_measured <- descr_df %>% filter(!is.na(e_e_ave), !is.na(la_vol_index))
reference_pool_gradient <- as_num(nonhcm_control_pool$lvot_max_gradient)
reference_contributing_gradient <- as_num(comparison_nonhcm_df$lvot_max_gradient)

manuscript_descriptives <- tibble::tribble(
  ~quantity, ~value,
  "Base analytic cohort, n", as.character(n_base),
  "Age, mean +/- SD, years", sprintf("%.1f +/- %.1f", mean(as_num(descr_df$age), na.rm = TRUE), sd(as_num(descr_df$age), na.rm = TRUE)),
  "Male sex, n (%)", fmt_n_pct(sum(descr_df$Sex == "Male", na.rm = TRUE), n_base),
  "BMI, median (IQR), kg/m2", fmt_med_iqr(descr_df$BMI),
  "LVEF measured, n", as.character(sum(!is.na(as_num(descr_df$ef_modsp4)))),
  "LVEF, median (IQR), %", fmt_med_iqr(descr_df$ef_modsp4),
  "Resting LVOT gradient recorded, n", as.character(gradient_recorded),
  "Resting LVOT gradient >= 30 mm Hg among recorded, n (%)", fmt_n_pct(obstructive_n, gradient_recorded),
  "Apical HCM, n", as.character(sum(descr_df$apical_hcm == 1, na.rm = TRUE)),
  "Morphology unrecorded, n", as.character(sum(is.na(descr_df$apical_hcm))),
  "Peak VO2, median (IQR), % predicted", fmt_med_iqr(descr_df$VO2_FRIEND2_PP),
  "VE/VCO2 slope, median (IQR)", fmt_med_iqr(descr_df$VeVco2_slope),
  "Baseline CPET on treadmill, n", as.character(sum(descr_df$cpx_mode == 0, na.rm = TRUE)),
  "Baseline CPET on cycle ergometer, n", as.character(sum(descr_df$cpx_mode == 1, na.rm = TRUE)),
  "Baseline CPET modality not recorded or invalid, n", as.character(sum(is.na(descr_df$cpx_mode) | !descr_df$cpx_mode %in% c(0, 1))),
  "Average e-prime measured, n", as.character(sum(!is.na(as_num(descr_df$e_prime_ave)))),
  "Average e-prime measured but age < 20, not classifiable, n", as.character(sum(!is.na(as_num(descr_df$e_prime_ave)) & as_num(descr_df$age) < 20, na.rm = TRUE)),
  "Average e-prime classifiable by ASE age bands, n", as.character(sum(descr_df$eprime_classifiable, na.rm = TRUE)),
  "Average e-prime below the age-specific limit, n (%) of classifiable", fmt_n_pct(sum(descr_df$eprime_below_age_limit, na.rm = TRUE), sum(descr_df$eprime_classifiable, na.rm = TRUE)),
  "Age < 18 years, n", as.character(sum(as_num(descr_df$age) < 18, na.rm = TRUE)),
  "Age < 20 years, n", as.character(sum(as_num(descr_df$age) < 20, na.rm = TRUE)),
  "E/e-prime measured, n", as.character(ee_measured),
  "E/e-prime > 14, n (%) of measured", fmt_n_pct(sum(descr_df$abn_ee, na.rm = TRUE), ee_measured),
  "LAVi measured, n", as.character(lavi_measured),
  "LAVi > 34 mL/m2, n (%) of measured", fmt_n_pct(sum(descr_df$abn_lavi, na.rm = TRUE), lavi_measured),
  "TRVmax measured, n", as.character(trv_measured),
  "TRVmax > 2.8 m/s, n (%) of measured", fmt_n_pct(sum(descr_df$abn_trv, na.rm = TRUE), trv_measured),
  "Both E/e-prime and LAVi measured, n", as.character(nrow(both_measured)),
  "Both elevated, n", as.character(sum(both_measured$abn_ee & both_measured$abn_lavi)),
  "LAVi elevated only, n", as.character(sum(!both_measured$abn_ee & both_measured$abn_lavi)),
  "E/e-prime elevated only, n", as.character(sum(both_measured$abn_ee & !both_measured$abn_lavi)),
  "Neither elevated, n", as.character(sum(!both_measured$abn_ee & !both_measured$abn_lavi)),
  "Reference pool resting LVOT gradient, median (IQR), mm Hg", fmt_med_iqr(reference_pool_gradient, 2),
  "Reference pool resting LVOT gradient, maximum, mm Hg", sprintf("%.2f", max(reference_pool_gradient, na.rm = TRUE)),
  "Contributing references resting LVOT gradient, maximum, mm Hg", sprintf("%.2f", max(reference_contributing_gradient, na.rm = TRUE))
)

# The ASE bands must partition the patients who have a measured e-prime:
# classifiable plus under-age must equal measured, or a band boundary is wrong.
stopifnot(
  nrow(manuscript_descriptives) == 35,
  !any(is.na(manuscript_descriptives$value)),
  sum(descr_df$eprime_classifiable, na.rm = TRUE) +
    sum(!is.na(as_num(descr_df$e_prime_ave)) & as_num(descr_df$age) < 20, na.rm = TRUE) ==
    sum(!is.na(as_num(descr_df$e_prime_ave))),
  all(descr_df$eprime_age_limit[descr_df$eprime_classifiable] %in% c(9, 7, 6.5)),
  # Modality must account for every patient in the cohort.
  sum(descr_df$cpx_mode == 0, na.rm = TRUE) +
    sum(descr_df$cpx_mode == 1, na.rm = TRUE) +
    sum(is.na(descr_df$cpx_mode) | !descr_df$cpx_mode %in% c(0, 1)) == nrow(descr_df)
)
write.csv(
  manuscript_descriptives,
  file.path(OUT_DIR, "Table_Manuscript_Descriptives.csv"),
  row.names = FALSE
)
print(as.data.frame(manuscript_descriptives), row.names = FALSE)

# ══ Figure 1: CONSORT-style cohort assembly ═══════════════════════════════════
# (legacy chunk: figure_consort)
# Two-arm CONSORT. The vertical spine is HCM cohort assembly; the upstream
# measurement-completeness step defines the cross-sectional analytic cohort,
# and each diastolic parameter branches off it with its own denominators. The bottom band
# shows the secondary weighted reference comparison as the join of the HCM and non-HCM
# screening arms, so the Table 1 population is traceable.
#
# Every count is computed from the live pipeline objects. Nothing here is a
# hard-coded constant, so the figure cannot drift from the analysis.

consort_pal <- c(
  ink = "#2C3E50",
  excl_fill = "#F7F1E8",
  excl_line = "#A36A2B",
  excl_ink = "#5B3A12",
  node_fill = "#FFFFFF",
  key_fill = "#EEF3F7",
  comp_fill = "#F4F6F4",
  comp_line = "#5A6E5A",
  rule = "#B8C2CC"
)

# Integer formatter that does not pad to a common width the way format() does,
# so "n = 823" never renders as "n =  823".
consort_int <- function(x) formatC(x, format = "d", big.mark = ",")

# Publication wording for the pipeline stages. filter_flow keeps its internal
# stage names for the CONSORT tally; these are display labels only, written out
# so the figure reads as Methods prose rather than as pipeline shorthand.
consort_stage_label <- function(stage) {
  dplyr::case_when(
    grepl("Initial CPX import", stage) ~
      "Cardiopulmonary exercise tests in registry",
    grepl("HCM diagnosis", stage) ~
      "HCM diagnosis and peak RER ≥1.0",
    grepl("excluding CAD/COPD/ILD", stage) ~
      paste0(
        "No coronary artery disease, chronic obstructive\n",
        "pulmonary disease, or interstitial lung disease"
      ),
    grepl("Echo aligned", stage) ~
      "Resting echocardiogram within ±7 days of CPET",
    grepl("phenotypic criteria", stage) ~
      paste0(
        "HCM phenotypic entry criteria met (maximum wall\n",
        "≥1.5 cm, apical, or variant with wall ≥1.3 cm)"
      ),
    TRUE ~ stringr::str_wrap(
      stringr::str_replace_all(
        stage,
        c("E/e'" = "E/e′", "\\+/-" = "±", ">=" = "≥")
      ),
      width = 46
    )
  )
}

# ── main HCM spine ────────────────────────────────────────────────────────────
# The common-cohort row is a consistency-check subset rather than a step in
# the assembly, so it is not part of the spine; the parameter branches below
# carry the analytic denominators.
consort_steps <- filter_flow %>%
  filter(
    stage != "Assigned stable patient ID",
    !stringr::str_detect(stage, "consistency check")
  ) %>%
  mutate(
    excluded = lag(n_patients, default = NA_integer_) - n_patients,
    excluded = ifelse(is.na(excluded) | excluded <= 0, NA_integer_, excluded),
    stage_raw = stage,
    stage = consort_stage_label(stage),
    y = rev(seq(6.0, 12.0, length.out = n())),
    label = sprintf(
      "%s\nn = %s patients%s",
      stage,
      consort_int(n_patients),
      ifelse(
        is.na(n_tests),
        "",
        sprintf("\n(%s CPET records)", consort_int(n_tests))
      )
    )
  )

# Box height follows the number of text lines so a long stage label cannot
# overflow its rectangle into the node below.
consort_steps <- consort_steps %>%
  mutate(
    n_label_lines = stringr::str_count(label, "\n") + 1L,
    half_height = pmax(0.50, 0.17 + 0.115 * n_label_lines)
  )

consort_excl <- consort_steps %>%
  filter(!is.na(excluded)) %>%
  mutate(
    y_excl = y + 0.6,
    label = if_else(
      stringr::str_detect(stage_raw, "Baseline analytic filter"),
      sprintf(
        "Excluded: %s patients\n%s",
        consort_int(excluded),
        stringr::str_wrap(
          stringr::str_replace_all(analytic_filter_exclusion, "E/e'", "E/e′"),
          width = 26
        )
      ),
      sprintf("Excluded: %s patients", consort_int(excluded))
    ),
    n_label_lines = stringr::str_count(label, "\\n") + 1L,
    half_height = 0.20 + 0.11 * (n_label_lines - 1L)
  )

# ── base analytic cohort and parameter-specific branches ──────────────────────
cross_n <- sum(baseline_df$cross_sectional_eligible, na.rm = TRUE)
base_n <- sum(baseline_df$analysis_base_eligible, na.rm = TRUE)
cross_trv_n <- sum(
  baseline_df$cross_sectional_eligible & !is.na(baseline_df$tr_max_vel),
  na.rm = TRUE
)
cross_terminal <- tibble(
  x = 0,
  y = 4.58,
  label = sprintf(
    paste0(
      "BASE ANALYTIC COHORT\n",
      "n = %d patients with complete age, sex, and BMI\n",
      "and ≥1 cardiopulmonary exercise outcome\n",
      "(both E/e′ and LAVi measured in %d)"
    ),
    base_n,
    cross_n
  )
)

# ── parameter-specific available-case samples ─────────────────────────────────
# Each parameter is analysed in every patient in whom it was measured, so the
# branches carry their own cross-sectional, longitudinal, and outcome counts.
parameter_branches <- purrr::map_dfr(
  seq_along(primary_parameters),
  function(i) {
    param <- primary_parameters[i]
    keep <- baseline_df[[paste0("eligible_", param)]]
    ids <- baseline_df$ID[keep]
    lng <- long_df %>% filter(ID %in% ids)
    out_keep <- keep & baseline_df$outcomes_followup_ok
    tibble(
      x = c(-2.70, 0.30, 3.30)[i],
      y = 3.12,
      label = sprintf(
        paste0(
          "%s available-case sample\n",
          "n = %d cross-sectional\n",
          "%d longitudinal (%d CPETs)\n",
          "%d outcomes (%d first events)"
        ),
        c("E/e′", "LAVi", "TRVmax")[i],
        sum(keep, na.rm = TRUE),
        n_distinct(lng$ID),
        nrow(lng),
        sum(out_keep, na.rm = TRUE),
        sum(out_keep & baseline_df$hf_composite == 1, na.rm = TRUE)
      )
    )
  }
)

# ── secondary weighted reference comparison (bottom band) ────────────────────
# Left node: HCM cases eligible for matching. Right nodes: the non-HCM
# screening arm. Both feed the matched box, which is the Table 1 population.
n_nonhcm_screened <- nrow(baseline_nonhcm_df)
n_nonhcm_qc <- sum(baseline_nonhcm_df$control_qc_eligible, na.rm = TRUE)
n_nonhcm_pool <- nrow(nonhcm_control_pool)
n_hcm_match_pool <- nrow(hcm_intro_df)
n_matched_hcm <- nrow(comparison_hcm_df)
n_matched_con <- nrow(comparison_nonhcm_df)

comp_nodes <- tibble(
  x = c(-2.45, 2.45, 2.45),
  y = c(1.55, 1.55, 0.30),
  label = c(
    sprintf(
      paste0(
        "HCM cohort entering the comparison\n",
        "n = %d of %d in the base analytic cohort"
      ),
      n_hcm_match_pool,
      base_n
    ),
    sprintf(
      paste0(
        "Non-HCM comparator candidates\n",
        "n = %s patients\n",
        "same effort, comorbidity, and echo-alignment\n",
        "criteria; no HCM code or registry record"
      ),
      consort_int(n_nonhcm_screened)
    ),
    sprintf(
      paste0(
        "Passed diagnostic and structural QC\n",
        "n = %d (LVEF ≥50%%, maximum wall <1.3 cm, and\n",
        "resting LVOT gradient <30 mm Hg where measured)\n",
        "n = %d with complete age, sex, and BMI"
      ),
      n_nonhcm_qc,
      n_nonhcm_pool
    )
  ),
  half_width = c(1.72, 1.86, 1.86),
  half_height = c(0.44, 0.44, 0.56)
)

matched_node <- tibble(
  x = 0,
  y = -1.18,
  label = sprintf(
    paste0(
      "WEIGHTED BASELINE COMPARISON (Table 1)\n",
      "all %d HCM and %d non-HCM references retained (total n = %d)\n",
      "propensity-score overlap weights on age, sex, and BMI\n",
      "effective n = %.0f and %.0f; |SMD| < 0.001 after weighting"
    ),
    n_matched_hcm,
    n_matched_con,
    n_matched_hcm + n_matched_con,
    comparison_weights$ess_cases,
    comparison_weights$ess_controls
  )
)

box_half_height <- 0.50
spine_half_width <- 1.92
cross_half_width <- 1.98
cross_half_height <- 0.62
subcohort_half_width <- 1.40
subcohort_half_height <- 0.62
matched_half_width <- 2.42
matched_half_height <- 0.58

y_split <- 3.88 # horizontal bus feeding the parameter branches
y_band <- 2.38 # rule separating primary analyses from the secondary comparison
y_join <- -0.58 # horizontal bus feeding the matched-comparison box

p_consort <- ggplot() +
  # main spine
  geom_rect(
    data = consort_steps,
    aes(
      xmin = -spine_half_width,
      xmax = spine_half_width,
      ymin = y - half_height,
      ymax = y + half_height
    ),
    fill = consort_pal[["node_fill"]],
    color = consort_pal[["ink"]],
    linewidth = 0.5
  ) +
  geom_text(
    data = consort_steps,
    aes(x = 0, y = y, label = label),
    family = "Arial",
    size = 2.85,
    lineheight = 1.02
  ) +
  geom_segment(
    data = consort_steps %>% arrange(desc(y)) %>% slice_head(n = -1),
    aes(
      x = 0,
      xend = 0,
      y = y - half_height,
      yend = y - 1.5 + half_height
    ),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    color = consort_pal[["ink"]],
    linewidth = 0.4
  ) +
  # exclusion call-outs
  geom_segment(
    data = consort_excl,
    aes(x = 0, xend = 2.28, y = y_excl, yend = y_excl),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
    color = consort_pal[["excl_line"]],
    linewidth = 0.35
  ) +
  geom_rect(
    data = consort_excl,
    aes(
      xmin = 2.28,
      xmax = 4.98,
      ymin = y_excl - half_height,
      ymax = y_excl + half_height
    ),
    fill = consort_pal[["excl_fill"]],
    color = consort_pal[["excl_line"]],
    linewidth = 0.4
  ) +
  geom_text(
    data = consort_excl,
    aes(x = 3.63, y = y_excl, label = label),
    family = "Arial",
    size = 2.7,
    lineheight = 1.02,
    color = consort_pal[["excl_ink"]]
  ) +
  # analytic filter -> cross-sectional parent cohort
  geom_segment(
    aes(
      x = 0,
      xend = 0,
      y = min(consort_steps$y) - max(consort_steps$half_height),
      yend = cross_terminal$y + cross_half_height
    ),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    color = consort_pal[["ink"]],
    linewidth = 0.4
  ) +
  geom_rect(
    data = cross_terminal,
    aes(
      xmin = x - cross_half_width,
      xmax = x + cross_half_width,
      ymin = y - cross_half_height,
      ymax = y + cross_half_height
    ),
    fill = consort_pal[["key_fill"]],
    color = consort_pal[["ink"]],
    linewidth = 0.7
  ) +
  geom_text(
    data = cross_terminal,
    aes(x = x, y = y, label = label),
    family = "Arial",
    size = 2.85,
    lineheight = 1.05,
    fontface = "bold"
  ) +
  # parameter-specific branches
  geom_segment(
    aes(
      x = 0,
      xend = 0,
      y = cross_terminal$y - cross_half_height,
      yend = y_split
    ),
    color = consort_pal[["ink"]],
    linewidth = 0.4
  ) +
  geom_segment(
    aes(
      x = min(parameter_branches$x),
      xend = max(parameter_branches$x),
      y = y_split,
      yend = y_split
    ),
    color = consort_pal[["ink"]],
    linewidth = 0.4
  ) +
  geom_segment(
    data = parameter_branches,
    aes(x = x, xend = x, y = y_split, yend = y + subcohort_half_height),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    color = consort_pal[["ink"]],
    linewidth = 0.4
  ) +
  geom_rect(
    data = parameter_branches,
    aes(
      xmin = x - subcohort_half_width,
      xmax = x + subcohort_half_width,
      ymin = y - subcohort_half_height,
      ymax = y + subcohort_half_height
    ),
    fill = consort_pal[["node_fill"]],
    color = consort_pal[["ink"]],
    linewidth = 0.5
  ) +
  geom_text(
    data = parameter_branches,
    aes(x = x, y = y, label = label),
    family = "Arial",
    size = 2.6,
    lineheight = 1.02
  ) +
  # rule + caption separating the secondary weighted comparison
  annotate(
    "segment",
    x = -4.15,
    xend = 4.98,
    y = y_band,
    yend = y_band,
    color = consort_pal[["rule"]],
    linewidth = 0.4,
    linetype = "22"
  ) +
  annotate(
    "text",
    x = -4.15,
    y = y_band - 0.08,
    label = "Secondary reference comparison (not used for inference)",
    hjust = 0,
    vjust = 1,
    family = "Arial",
    size = 2.5,
    fontface = "italic",
    color = consort_pal[["comp_line"]]
  ) +
  # comparator arm
  geom_rect(
    data = comp_nodes,
    aes(
      xmin = x - half_width,
      xmax = x + half_width,
      ymin = y - half_height,
      ymax = y + half_height
    ),
    fill = consort_pal[["comp_fill"]],
    color = consort_pal[["comp_line"]],
    linewidth = 0.45
  ) +
  geom_text(
    data = comp_nodes,
    aes(x = x, y = y, label = label),
    family = "Arial",
    size = 2.5,
    lineheight = 1.02,
    color = "#243024"
  ) +
  # non-HCM candidates -> QC
  geom_segment(
    aes(x = 2.45, xend = 2.45, y = 1.55 - 0.44, yend = 0.30 + 0.56),
    arrow = arrow(length = unit(0.14, "cm"), type = "closed"),
    color = consort_pal[["comp_line"]],
    linewidth = 0.4
  ) +
  # both arms -> weighted comparison
  geom_segment(
    aes(x = -2.45, xend = -2.45, y = 1.55 - 0.44, yend = y_join),
    color = consort_pal[["comp_line"]],
    linewidth = 0.4
  ) +
  geom_segment(
    aes(x = 2.45, xend = 2.45, y = 0.30 - 0.56, yend = y_join),
    color = consort_pal[["comp_line"]],
    linewidth = 0.4
  ) +
  geom_segment(
    aes(x = -2.45, xend = 2.45, y = y_join, yend = y_join),
    color = consort_pal[["comp_line"]],
    linewidth = 0.4
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = 0,
      y = y_join,
      yend = matched_node$y + matched_half_height
    ),
    arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
    color = consort_pal[["comp_line"]],
    linewidth = 0.4
  ) +
  geom_rect(
    data = matched_node,
    aes(
      xmin = x - matched_half_width,
      xmax = x + matched_half_width,
      ymin = y - matched_half_height,
      ymax = y + matched_half_height
    ),
    fill = consort_pal[["comp_fill"]],
    color = consort_pal[["comp_line"]],
    linewidth = 0.7
  ) +
  geom_text(
    data = matched_node,
    aes(x = x, y = y, label = label),
    family = "Arial",
    size = 2.6,
    lineheight = 1.05,
    fontface = "bold",
    color = "#243024"
  ) +
  scale_x_continuous(limits = c(-4.25, 5.1), expand = c(0, 0)) +
  scale_y_continuous(
    limits = c(-1.95, max(consort_steps$y) + 0.75),
    expand = c(0, 0)
  ) +
  labs(title = "Figure 1. Cohort assembly and analysis populations") +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 11,
      margin = margin(b = 6)
    ),
    plot.margin = margin(8, 8, 8, 8)
  )

# Render-time guards: the parameter branches and the matched arms must be
# consistent with the cohorts they are drawn from, so a layout edit can never
# silently ship a figure whose arithmetic disagrees with the analysis.
stopifnot(
  cross_n <= base_n,
  n_distinct(long_df$ID) <= base_n,
  sum(baseline_df$analysis_base_eligible & baseline_df$outcomes_followup_ok, na.rm = TRUE) <= base_n,
  n_hcm_match_pool <= base_n,
  n_matched_hcm <= n_hcm_match_pool,
  n_matched_con <= n_nonhcm_pool,
  n_nonhcm_pool <= n_nonhcm_qc,
  n_nonhcm_qc <= n_nonhcm_screened
)

legend_consort_pdf <- paste0(
  "**Figure 1. Cohort assembly and analysis populations.** ",
  "The vertical spine reports sequential HCM cohort entry, with exclusions at each step shown at right. ",
  "Patients meeting phenotypic entry criteria with complete age, sex, and body mass index and at least one cardiopulmonary exercise outcome form the base analytic cohort. ",
  "Each resting LV diastolic parameter is then analyzed in every patient in whom that parameter was measured, so the three branches carry their own cross-sectional, longitudinal, and clinical-outcome denominators; ",
  "requiring two or three parameters simultaneously would remove patients from the analysis of the parameter that was measured least often. ",
  "The subset with both E/e' and LAVi measured is retained as a prespecified consistency check. ",
  "The secondary baseline comparison joins the HCM cohort to phenotypically screened non-HCM comparator candidates and supports Table 1 only. ",
  "Every eligible patient is retained and comparability is achieved by propensity-score overlap weighting on age, sex, and body mass index rather than by discarding unmatched patients; ",
  "absolute standardized differences are below 0.001 after weighting (0.92 for age before), and effective sample sizes are reported in the box. ",
  "One-to-one caliper matching of the same pools is reported as a sensitivity analysis. ",
  "All counts are computed by the analysis pipeline at render time. ",
  "CPET = cardiopulmonary exercise testing; LAVi = left atrial volume index; LVEF = left ventricular ejection fraction; ",
  "LVOT = left ventricular outflow tract; QC = quality control; RER = respiratory exchange ratio; ",
  "TRVmax = maximum tricuspid regurgitant velocity."
)

save_jacc_with_embedded_legend(
  p_consort,
  file.path(OUT_DIR, "Figure1_CohortFlow.pdf"),
  legend_text = legend_consort_pdf,
  w = 7.5,
  h = 10.0,
  legend_fontsize = 8.5
)
save_jacc_with_embedded_legend(
  p_consort,
  file.path(OUT_DIR, "Figure1_CohortFlow.pdf"),
  legend_text = legend_consort_pdf,
  w = 7.5,
  h = 10.0,
  legend_fontsize = 8.5
)

# (legacy chunk: figure_consort_display)
print(p_consort)

# ══ Figure 1: Export ══════════════════════════════════════════════════════════
# (legacy chunk: figure1_plot)
# ── Fig 1: CON vs nonobstructive vs obstructive HCM across 6 LVDD params ───────
# assemble plot df, then build/save the multi-panel comparison figure
validation_df <- case_control_hist_df %>%
  transmute(
    Cohort,
    lv_max_wall_thickness = as_num(lv_max_wall_thickness),
    lv_septal_thickness = as_num(lv_septal_thickness),
    lvot_max_gradient = as_num(lvot_max_gradient),
    la_vol_index = as_num(la_vol_index),
    ef_modsp4 = as_num(ef_modsp4),
    avg_eprime_age_cal = as_num(avg_eprime_age_cal)
  ) %>%
  mutate(Cohort = factor(
    Cohort,
    levels = c("CON", "Non-obstructive HCM", "Obstructive HCM")
  ))

# Draw order: obstructive HCM first (behind), then non-obstructive HCM, then CON
cohort_draw_order <- c("Obstructive HCM", "Non-obstructive HCM", "CON")

cohort_fill_palette <- c(
  "CON" = "darkgray",
  "Non-obstructive HCM" = "#C9A2A0",
  "Obstructive HCM" = "#8BAFC4"
)

# ══ Supplemental: propensity-score overlap and the weighted comparison ═══════
# Owner decision, card 2.2 of 8_Docs/DECISION_CARDS_2026-09.md. The unweighted
# age difference between the HCM cohort and the reference pool is large, so the
# common support is shown rather than asserted: panel A is the propensity
# distribution by group as observed, panel B the same distributions under the
# overlap weights that the comparison actually uses.
overlap_ps_df <- bind_rows(
  comparison_weights$cases %>%
    transmute(group = "HCM", ps = .ps, w = .weight),
  comparison_weights$controls %>%
    transmute(group = "Reference", ps = .ps, w = .weight)
) %>%
  mutate(group = factor(group, levels = c("Reference", "HCM")))

stopifnot(
  nrow(overlap_ps_df) ==
    comparison_weights$n_cases + comparison_weights$n_controls,
  all(overlap_ps_df$ps > 0 & overlap_ps_df$ps < 1),
  all(overlap_ps_df$w >= 0)
)

overlap_palette <- c(
  Reference = unname(cohort_fill_palette[["CON"]]),
  HCM = unname(cohort_fill_palette[["Non-obstructive HCM"]])
)

overlap_density_panel <- function(df, weighted, subtitle) {
  mapping <- if (weighted) {
    aes(x = ps, fill = group, colour = group, weight = w)
  } else {
    aes(x = ps, fill = group, colour = group)
  }
  ggplot(df, mapping) +
    geom_density(alpha = 0.45, linewidth = 0.4, adjust = 1.1) +
    scale_fill_manual(values = overlap_palette, name = NULL) +
    scale_colour_manual(values = overlap_palette, name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      x = "Propensity score (P[HCM] given age, sex, BMI)",
      y = "Density",
      subtitle = subtitle
    ) +
    theme_jacc() +
    theme(legend.position = "bottom")
}

p_overlap <- (
  overlap_density_panel(
    overlap_ps_df,
    weighted = FALSE,
    subtitle = sprintf(
      "Unweighted\n(n = %d HCM, %d ref)",
      comparison_weights$n_cases,
      comparison_weights$n_controls
    )
  ) |
    overlap_density_panel(
      overlap_ps_df,
      weighted = TRUE,
      subtitle = sprintf(
        "Overlap-weighted\n(effective n = %.0f HCM, %.0f ref)",
        comparison_weights$ess_cases,
        comparison_weights$ess_controls
      )
    )
) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 11), legend.position = "bottom")

show_and_save_jacc(
  p_overlap,
  file.path(OUT_DIR, "FigureS_PropensityOverlap.pdf"),
  w = 6.9,
  h = 3.6
)

# Weighted and unweighted covariate moments behind the SMDs.
overlap_covariate_labels <- c(
  age = "Age, years",
  sex_num = "Female sex, proportion",
  BMI = "BMI, kg/m2"
)
comparison_covariate_table <- comparison_weights$covariate_summary %>%
  mutate(
    Covariate = unname(overlap_covariate_labels[variable]),
    .before = 1
  ) %>%
  select(-variable)
stopifnot(!any(is.na(comparison_covariate_table$Covariate)))

write.csv(
  comparison_covariate_table,
  file.path(OUT_DIR, "Table_Comparison_Covariates_WeightedUnweighted.csv"),
  row.names = FALSE
)

# The estimand follows from the weights, so state it in the words the Methods
# should use rather than leaving a reader to infer it from the word "weighted".
comparison_estimand_statement <- sprintf(
  paste0(
    "The reference comparison is weighted by overlap weights (1 - e for HCM ",
    "patients, e for references, where e is the propensity score from a ",
    "logistic model on age, sex and body mass index). The estimand is ",
    "therefore the average treatment effect in the overlap population: the ",
    "contrast that would be observed among patients whose age, sex and body ",
    "mass index are compatible with membership of either group, rather than ",
    "in the HCM cohort as a whole or in the reference pool as a whole. All ",
    "%d eligible HCM patients and %d references contribute, with effective ",
    "sample sizes of %.0f and %.0f; the weights balance the three covariates ",
    "exactly, reducing the largest standardised mean difference from %.2f to ",
    "%.3f."
  ),
  comparison_weights$n_cases,
  comparison_weights$n_controls,
  comparison_weights$ess_cases,
  comparison_weights$ess_controls,
  max(abs(comparison_weights$balance$`SMD (unweighted)`)),
  max(abs(comparison_weights$balance$`SMD (overlap-weighted)`))
)
writeLines(
  c(
    "# Estimand statement for the reference comparison",
    "",
    "Generated by R/04_cohort_table1.R; paste into the Methods.",
    "",
    comparison_estimand_statement
  ),
  file.path(OUT_DIR, "Comparison_Estimand_Statement.md")
)
message(comparison_estimand_statement)

# ══ Table 1: Baseline Characteristics ═════════════════════════════════════════
# (legacy chunk: table1_inline)
# ── Table 1: baseline characteristics (gtsummary) ─────────────────────────────
make_yn_factor <- function(x) {
  factor(
    case_when(x == 1 ~ "Yes", x == 0 ~ "No", TRUE ~ NA_character_),
    levels = c("No", "Yes")
  )
}

table1_casecontrol_data <- bind_rows(
  comparison_nonhcm_df %>%
    mutate(
      table1_group = "Non-HCM reference",
      HCM_Phenotype = NA_character_
    ),
  comparison_hcm_df %>%
    mutate(
      # Owner decision (2026-09-18): Table 1 has no separate group for patients
      # without a recorded resting LVOT gradient. Obstruction requires a
      # measured gradient >= 30 mm Hg, so an unmeasured gradient is counted as
      # non-obstructive, consistent with the obstructive_hcm row below.
      table1_group = if_else(
        figure1_hist_group == "HCM, gradient not measured",
        "Non-obstructive HCM",
        figure1_hist_group
      ),
      HCM_Phenotype = as.character(HCM_Phenotype)
    )
) %>%
  mutate(
    table1_group = factor(
      table1_group,
      levels = c(
        "Non-HCM reference",
        "Non-obstructive HCM",
        "Obstructive HCM"
      )
    ),
    across(
      any_of(c(
        "VeVco2_slope",
        "ivsd",
        "lvpwd",
        "la_vol",
        "pulm_sys_vel",
        "pulm_dias_vel",
        "HRmax_PP"
      )),
      ~ as_num(.)
    ),
    # Demographics & clinical
    beta_blocker_use = make_yn_factor(bb_any),
    ndhp_ccb_use = make_yn_factor(ndhp_ccb_any),
    disopyramide_use = make_yn_factor(as_num(disopyramide)),
    acei_arb_use = make_yn_factor(as_num(`ACEI/ARB`)),
    diuretic_use = make_yn_factor(as_num(Diuretics)),
    statin_use = make_yn_factor(as_num(Statin)),
    diabetes_pre_test = make_yn_factor(dm_pre_test),
    hypertension = make_yn_factor(htn_pre_test),
    HCM_Phenotype = factor(
      case_when(
        is.na(HCM_Phenotype) | HCM_Phenotype == "Burned-out" ~ NA_character_,
        TRUE ~ HCM_Phenotype
      ),
      levels = c("Asymmetric Septal", "Symmetric", "Apical")
    ),
    mv_dec_time_ms = 1000 * as_num(mv_dec_time),
    obstructive_hcm = make_yn_factor(
      ifelse(!is.na(lvot_max_gradient) & lvot_max_gradient >= 30, 1, 0)
    )
  ) %>%
  select(
    table1_group,
    # Demographics
    age,
    Sex,
    BMI,
    # HCM phenotype
    HCM_Phenotype,
    # Comorbidities
    hypertension,
    diabetes_pre_test,
    # Medications
    beta_blocker_use,
    ndhp_ccb_use,
    disopyramide_use,
    acei_arb_use,
    diuretic_use,
    statin_use,
    # LV structure (ivsd removed — duplicate of lv_septal_thickness)
    lvot_max_gradient,
    lv_septal_thickness,
    lvpwd,
    ef_modsp4,
    # Diastolic indices
    e_prime_ave,
    med_peak_e_vel,
    lat_peak_e_vel,
    e_e_ave,
    la_vol_index,
    la_vol,
    tr_max_vel,
    mv_e_a,
    mv_dec_time_ms,
    # CPET
    pk.RER,
    HRmax_PP,
    VO2_FRIEND2_PP,
    VeVco2_slope,
    HRR
  )

# Every HCM patient must land in exactly one of the two HCM columns, and no
# patient may fall out of the factor levels.
stopifnot(
  !any(is.na(table1_casecontrol_data$table1_group)),
  sum(table1_casecontrol_data$table1_group %in%
    c("Non-obstructive HCM", "Obstructive HCM")) == nrow(comparison_hcm_df),
  sum(table1_casecontrol_data$table1_group == "Obstructive HCM") ==
    sum(as_num(comparison_hcm_df$lvot_max_gradient) >= 30, na.rm = TRUE)
)

table1_nonobstructive_n <- sum(
  table1_casecontrol_data$table1_group == "Non-obstructive HCM"
)
table1_nonobstructive_gradient_n <- sum(
  table1_casecontrol_data$table1_group == "Non-obstructive HCM" &
    !is.na(as_num(table1_casecontrol_data$lvot_max_gradient))
)

table1_overall <- table1_casecontrol_data %>%
  tbl_summary(
    by = table1_group,
    label = list(
      # Demographics
      age ~ "Age, years",
      Sex ~ "Female sex",
      BMI ~ "Body mass index, kg/m\u00b2",
      # HCM phenotype
      HCM_Phenotype ~ "HCM phenotype",
      # Comorbidities
      hypertension ~ "Hypertension",
      diabetes_pre_test ~ "Diabetes mellitus",
      # Medications
      beta_blocker_use ~ "\u03b2-blocker",
      ndhp_ccb_use ~ "Non-DHP CCB",
      disopyramide_use ~ "Disopyramide",
      acei_arb_use ~ "ACEi/ARB",
      diuretic_use ~ "Diuretic",
      statin_use ~ "Statin",
      # LV structure
      lvot_max_gradient ~ "Resting LVOT gradient, mm Hg",
      lv_septal_thickness ~ "IVSd, cm",
      lvpwd ~ "LVPWd, cm",
      ef_modsp4 ~ "LVEF, %",
      # Diastolic function
      e_prime_ave ~ "Average e\u2032, cm/s",
      med_peak_e_vel ~ "Septal e\u2032, cm/s",
      lat_peak_e_vel ~ "Lateral e\u2032, cm/s",
      e_e_ave ~ "E/e\u2032 (average)",
      la_vol_index ~ "LAVi, mL/m\u00b2",
      la_vol ~ "LA volume, mL",
      tr_max_vel ~ "TRVmax, cm/s",
      mv_e_a ~ "E/A ratio",
      mv_dec_time_ms ~ "MV deceleration time, ms",
      # CPET
      pk.RER ~ "Peak RER",
      HRmax_PP ~ "Max HR, % predicted",
      VO2_FRIEND2_PP ~ "Peak V\u0307O2, FRIEND 2.0 %pred",
      VeVco2_slope ~ "V\u0307E/V\u0307CO2 slope",
      HRR ~ "Heart rate recovery (1 min)"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)",
      all_dichotomous() ~ "{n} ({p}%)"
    ),
    type = list(
      Sex ~ "dichotomous",
      beta_blocker_use ~ "dichotomous",
      ndhp_ccb_use ~ "dichotomous",
      disopyramide_use ~ "dichotomous",
      acei_arb_use ~ "dichotomous",
      diuretic_use ~ "dichotomous",
      statin_use ~ "dichotomous",
      hypertension ~ "dichotomous",
      diabetes_pre_test ~ "dichotomous"
    ),
    value = list(
      Sex ~ "Female",
      beta_blocker_use ~ "Yes",
      ndhp_ccb_use ~ "Yes",
      disopyramide_use ~ "Yes",
      acei_arb_use ~ "Yes",
      diuretic_use ~ "Yes",
      statin_use ~ "Yes",
      hypertension ~ "Yes",
      diabetes_pre_test ~ "Yes"
    ),
    digits = list(
      tr_max_vel ~ 0,
      pk.RER ~ 2,
      all_continuous() ~ 1
    ),
    missing = "no"
  )

# --- Convert to tibble, insert section headers, build flextable ---
t1_tib <- table1_overall %>%
  as_tibble(col_labels = FALSE)

# Strip any markdown bold formatting from labels (e.g., **Age, years** -> Age, years)
names(t1_tib) <- paste0("V", seq_len(ncol(t1_tib)))
t1_tib <- t1_tib %>%
  mutate(V1 = gsub("\\*\\*(.+?)\\*\\*", "\\1", V1))

# Column headers with group Ns
group_ns <- table1_casecontrol_data %>%
  count(table1_group) %>%
  arrange(table1_group)
col_headers <- c(
  "Characteristic",
  paste0(group_ns$table1_group, "\n(N = ", group_ns$n, ")")
)

# Helper: create a section header row
make_section_row <- function(label) {
  row <- as.list(rep("", ncol(t1_tib)))
  names(row) <- names(t1_tib)
  row[[1]] <- label
  as_tibble(row)
}

# Remove the gtsummary-generated "HCM phenotype" label row (section header replaces it)
t1_tib <- t1_tib %>% filter(V1 != "HCM phenotype")

# Replace CON column (V2) with "—" for HCM-specific morphology rows.
nonhcm_dash_vars <- c(
  "Asymmetric Septal",
  "Symmetric",
  "Apical"
)
t1_tib <- t1_tib %>%
  mutate(
    V2 = ifelse(V1 %in% nonhcm_dash_vars, "\u2014", V2)
  )

# Section headers mapped to first variable in each group
section_map <- list(
  "Demographics" = "Age, years",
  "HCM Phenotype" = "Asymmetric Septal",
  "Comorbidities" = "Hypertension",
  "Medications" = "\u03b2-blocker",
  "LV Structure & Function" = "Resting LVOT gradient, mm Hg",
  "Diastolic Function" = "Average e\u2032, cm/s",
  "Cardiopulmonary Exercise Testing" = "Peak RER"
)

# Insert section headers in reverse order to preserve row indices
for (sec_name in rev(names(section_map))) {
  first_var <- section_map[[sec_name]]
  idx <- which(t1_tib$V1 == first_var)
  if (length(idx) == 0) {
    # Try partial match in case of encoding differences
    idx <- grep(substr(first_var, 1, 8), t1_tib$V1, fixed = TRUE)
  }
  if (length(idx) > 0) {
    idx <- idx[1]
    t1_tib <- bind_rows(
      t1_tib[seq_len(idx - 1), ],
      make_section_row(sec_name),
      t1_tib[idx:nrow(t1_tib), ]
    )
  }
}

# Identify section header row indices
section_indices <- which(t1_tib$V1 %in% names(section_map))

# Build flextable
thin_border <- officer::fp_border(color = "black", width = 0.8)

table1_ft <- flextable(t1_tib) %>%
  set_header_labels(values = setNames(col_headers, names(t1_tib))) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  bold(part = "header") %>%
  line_spacing(space = 1.0, part = "all") %>%
  # Compact padding for body
  padding(
    padding.top = 1,
    padding.bottom = 1,
    padding.left = 3,
    padding.right = 3,
    part = "body"
  ) %>%
  padding(
    padding.top = 2,
    padding.bottom = 2,
    padding.left = 3,
    padding.right = 3,
    part = "header"
  ) %>%
  # Three-line publication border style
  border_remove() %>%
  hline_top(border = thin_border, part = "header") %>%
  hline_bottom(border = thin_border, part = "header") %>%
  hline_bottom(border = thin_border, part = "body") %>%
  # Section header formatting: bold, italic, flush left
  bold(i = section_indices, j = 1, part = "body") %>%
  italic(i = section_indices, j = 1, part = "body") %>%
  padding(i = section_indices, j = 1, padding.left = 3, part = "body") %>%
  # Indent variable rows under section headers
  padding(
    i = setdiff(seq_len(nrow(t1_tib)), section_indices),
    j = 1,
    padding.left = 14,
    part = "body"
  ) %>%
  # Center-align data columns
  align(j = 2:ncol(t1_tib), align = "center", part = "all") %>%
  # Column widths
  width(j = 1, width = 2.2) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption(paste0(
    "Table 1. Baseline Characteristics by Cohort. Obstructive HCM denotes a ",
    "resting LVOT gradient \u2265 30 mm Hg; patients without a recorded resting ",
    "gradient are included with non-obstructive HCM, so the gradient summary ",
    "for that column describes the ",
    table1_nonobstructive_gradient_n, " of ", table1_nonobstructive_n,
    " patients with a recorded value."
  ))

table1_ft <- apply_subscript_format(table1_ft, t1_tib, "V1")

table1_ft %>%
  save_as_docx(path = file.path(OUT_DIR, "Table1_Manuscript.docx"))


message("Stage 04 complete.")
