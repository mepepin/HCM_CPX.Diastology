# Cohort assembly helpers: the only definitions of registry preparation,
# comorbidity exclusion, echo alignment, HCM entry, baseline selection, and
# analytic eligibility. Ported from _Scripts/R/utils.R (lines 360-511).

# ── Cohort assembly: single source of truth ───────────────────────────────────
# These helpers hold the ONLY definitions of registry preparation, comorbidity
# exclusion, HCM study entry, baseline-visit selection, and analytic
# eligibility. Every cohort is built by calling them — the primary cohort
# (RER >= 1.0) and the all-RER sensitivity cohort differ in exactly one
# argument, `rer_min`. Nothing is re-derived per analysis, so an edit to study
# entry or eligibility propagates everywhere instead of drifting between
# hand-rolled copies.

# CPX-level preparation. `rer_min = NULL` keeps submaximal tests (all-RER).
prepare_cpx_registry <- function(cpx_with_id, rer_min = 1.0) {
  out <- cpx_with_id
  if (!is.null(rer_min)) {
    out <- out %>% dplyr::filter(pk.RER >= rer_min)
  }
  out %>%
    dplyr::mutate(
      diag_primary = trimws(as.character(MCl3)),
      diag_secondary = trimws(as.character(MCl4)),
      has_any_hcm_flag = stringr::str_detect(
        stringr::str_to_upper(
          paste0(
            dplyr::coalesce(diag_primary, ""),
            " ",
            dplyr::coalesce(diag_secondary, "")
          )
        ),
        "HCM"
      )
    )
}

# Matching on MRN alone would let a comorbidity-positive visit survive via a
# different clean visit from the same patient, so the join is MRN + test date.
apply_comorbidity_exclusions <- function(cpx_clean_df, comorbidities_to_filter) {
  cpx_clean_df %>%
    dplyr::mutate(
      MRN = as.character(MRN),
      cpx_test_date = as.Date(cpx_test_date)
    ) %>%
    dplyr::semi_join(comorbidities_to_filter, by = c("MRN", "cpx_test_date"))
}

# Genetics merge + AHA/ACC phenotypic entry criteria, per aligned visit.
derive_hcm_pre_entry <- function(cpx_echo_df, genetics_flags) {
  n_pre_join <- nrow(cpx_echo_df)
  out <- cpx_echo_df %>% dplyr::left_join(genetics_flags, by = "MRN")
  stopifnot(nrow(out) == n_pre_join)

  out %>%
    dplyr::mutate(
      HCM_Phenotype = dplyr::coalesce(
        genetics_phenotype,
        as.character(HCM_Phenotype)
      ),
      apical_hcm = as.integer(HCM_Phenotype == "Apical"),
      hcm_selection_reason = dplyr::case_when(
        !is.na(lv_max_wall_thickness) & lv_max_wall_thickness >= 1.5 ~
          "Wall thickness >= 1.5 cm",
        HCM_Phenotype == "Apical" ~ "Apical HCM",
        # Gene+ with wall >= 1.3 cm meets AHA/ACC phenotypic HCM criteria
        pathogenic_variant == 1 &
          !is.na(lv_max_wall_thickness) &
          lv_max_wall_thickness >= 1.3 ~
          "Wall thickness >= 1.3 cm (gene+)",
        pathogenic_variant == 1 ~ "Pathogenic variant carrier",
        TRUE ~ "Excluded unclear HCM designation"
      ),
      hcm_selection_include = hcm_selection_reason %in%
        c(
          "Wall thickness >= 1.5 cm",
          "Apical HCM",
          "Wall thickness >= 1.3 cm (gene+)"
        )
    )
}

# Establish HCM entry at the first aligned examination meeting phenotypic
# criteria, then retain all subsequent eligible visits. This prevents ordinary
# longitudinal measurement variation in wall thickness from dropping later
# CPX/echo pairs after a patient has met the study definition.
apply_hcm_entry_window <- function(pre_entry_df) {
  hcm_index <- pre_entry_df %>%
    dplyr::filter(hcm_selection_include) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(
      hcm_index_date = min(as.Date(cpx_test_date)),
      .groups = "drop"
    )

  pre_entry_df %>%
    dplyr::inner_join(hcm_index, by = "ID") %>%
    dplyr::filter(as.Date(cpx_test_date) >= hcm_index_date) %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(cpx_test_date, .by_group = TRUE) %>%
    dplyr::mutate(
      days_from_baseline = as.numeric(
        as.Date(cpx_test_date) - dplyr::first(as.Date(cpx_test_date))
      )
    ) %>%
    dplyr::ungroup()
}

# Earliest aligned visit per patient, with the numeric coercions every
# downstream model relies on.
derive_baseline_visits <- function(cpx_echo_df, num_vars) {
  cpx_echo_df %>%
    dplyr::group_by(ID) %>%
    dplyr::arrange(days_from_baseline) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(
      dplyr::any_of(num_vars),
      ~ as.numeric(as.character(.))
    ))
}

# Analytic eligibility: complete age/sex/BMI + the measurements specified by
# `analytic_cohort_mode` + >=1 CPX outcome (peak VO2, VE/VCO2, or HRR).
add_cross_sectional_eligible <- function(df, analytic_cohort_mode) {
  out <- df %>%
    dplyr::mutate(
      cross_sectional_eligible = !is.na(as_num(age)) &
        Sex %in% c("Male", "Female") &
        !is.na(as_num(BMI)) &
        !is.na(as_num(e_e_ave)) &
        !is.na(as_num(la_vol_index)) &
        (analytic_cohort_mode != "complete_three" |
          !is.na(as_num(tr_max_vel))) &
        (!is.na(as_num(VO2_FRIEND2_PP)) |
          !is.na(as_num(VeVco2_slope)) |
          !is.na(as_num(HRR)))
    )

  required <- if (analytic_cohort_mode == "complete_three") {
    c("e_e_ave", "la_vol_index", "tr_max_vel")
  } else {
    c("e_e_ave", "la_vol_index")
  }
  stopifnot(
    out %>%
      dplyr::filter(cross_sectional_eligible) %>%
      dplyr::summarise(all_complete = all(!is.na(
        dplyr::pick(dplyr::all_of(required))
      ))) %>%
      dplyr::pull(all_complete)
  )

  out
}


# ── Analytic cohort labels ────────────────────────────────────────────────────
analytic_filter_labels <- function(analytic_cohort_mode) {
  stopifnot(analytic_cohort_mode %in% c("complete_three", "ee_lavi"))
  list(
    # Not a pipeline exclusion any more: each parameter is analysed in its own
    # available-case sample, and this cohort is the consistency check.
    stage = switch(
      analytic_cohort_mode,
      complete_three = "Common E/e' + LAVi + TRVmax cohort (consistency check)",
      ee_lavi = "Common E/e' + LAVi cohort (consistency check)"
    ),
    exclusion = switch(
      analytic_cohort_mode,
      complete_three = "E/e', LAVi, and/or TRVmax unavailable",
      ee_lavi = "E/e' and/or LAVi unavailable"
    ),
    description = switch(
      analytic_cohort_mode,
      complete_three = "complete E/e', LAVi, and TRVmax",
      ee_lavi = "complete E/e' and LAVi; TRVmax available-case"
    )
  )
}

# ── CPX registry formatting ───────────────────────────────────────────────────
# Shared by HCM and non-HCM registries: longitudinal time axis, analysis
# variable names, and removal of implausible %-predicted VO2 (outside 0-200).
# Derived fields (age, BMI, BSA, LBM, LBMI, VO2_*_PP, HRmax_PP) come from
# derive_cpx_fields() in R/derive.R.
format_cpx_registry <- function(df, keep_diagnosis = FALSE) {
  out <- df %>%
    # days since patient's first test = longitudinal time axis
    group_by(ID) %>%
    arrange(cpx_test_date, .by_group = TRUE) %>%
    mutate(
      days_from_baseline = as.numeric(cpx_test_date - first(cpx_test_date))
    ) %>%
    ungroup() %>%
    select(
      ID,
      MRN,
      cpx_test_date,
      age,
      Sex = `Sex:.M(0)/F(1)`,
      BMI,
      BSA,
      LBM,
      LBMI,
      HCM_Phenotype = MCl4,
      registry_diagnosis_comment = Comments_ddx,
      registry_hcm_field = HCM,
      CPX.sequential,
      days_from_baseline,
      pk.RER,
      VO2_FRIEND_PP,
      VO2_WASSERMAN_PP,
      VO2_FRIEND2_PP,
      HRmax_PP,
      HRR = `HR.recovery.(1min)`,
      VeVco2_slope = `ve/vco2.slope`,
      diag_primary,
      diag_secondary,
      has_any_hcm_flag
    ) %>%
    # %-pred VO2 outside (0,200] = data entry error
    filter(
      (is.na(VO2_FRIEND2_PP) | (VO2_FRIEND2_PP >= 0 & VO2_FRIEND2_PP <= 200)) &
        (is.na(VO2_FRIEND_PP) | (VO2_FRIEND_PP >= 0 & VO2_FRIEND_PP <= 200))
    )

  # dx text cols only needed to vet controls
  if (!keep_diagnosis) {
    out <- out %>% select(-diag_primary, -diag_secondary, -has_any_hcm_flag)
  }

  out
}

# CPX records at which CAD, COPD, and ILD were all absent (NA = absent).
# Keyed on MRN + test date so a clean visit cannot retain a different,
# comorbidity-positive visit from the same patient.
comorbidity_clean_visits <- function(comorbidities_table) {
  comorbidities_table %>%
    transmute(
      MRN = as.character(MRN),
      cpx_test_date = as.Date(cpx_test_date),
      cad_pre_test,
      copd_pre_test,
      interstitial_lung_dz_pre_test
    ) %>%
    filter(
      (cad_pre_test != 1 | is.na(cad_pre_test)) &
        (copd_pre_test != 1 | is.na(copd_pre_test)) &
        (interstitial_lung_dz_pre_test != 1 |
          is.na(interstitial_lung_dz_pre_test))
    ) %>%
    distinct(MRN, cpx_test_date)
}

# ── Echocardiography ──────────────────────────────────────────────────────────
echo_vars_core <- c(
  "e_e_ave",
  "e_e_lat",
  "e_e_med",
  "mv_a_dur",
  "mv_a_point",
  "mv_dec_time",
  "mv_e_a",
  "tr_max_vel",
  "la_vol_index",
  "med_peak_e_vel",
  "lat_peak_e_vel",
  "ef_modsp4",
  "la_vol",
  "pulm_dias_vel",
  "pulm_sys_vel",
  # `lv_max_pg` is the LV outflow/intracavitary Doppler gradient. It closely
  # reproduces 4 * (lv_v1_max / 100)^2 in the source extract. The generic
  # `max_pg` field is retained only for provenance checks and is not used to
  # classify obstruction.
  "lv_max_pg",
  "lv_v1_max",
  "max_pg",
  "ivsd",
  "lvpwd",
  "ivs_lvpw"
)

prepare_echo_measures <- function(echo_raw) {
  echo_raw %>%
    mutate(MRN = as.character(mrn), echo_date = as.Date(echo_date)) %>%
    select(MRN, echo_date, all_of(echo_vars_core)) %>%
    mutate(across(all_of(echo_vars_core), ~ as.numeric(as.character(.))))
}

# The Doppler pressure gradient should reproduce the simplified Bernoulli
# calculation from LV velocity; exported so the obstruction definition is
# independently auditable.
summarize_lvot_provenance <- function(echo_all) {
  echo_all %>%
    transmute(
      lv_max_pg,
      lv_v1_max,
      bernoulli_pg = 4 * (lv_v1_max / 100)^2,
      abs_difference = abs(lv_max_pg - bernoulli_pg)
    ) %>%
    filter(complete.cases(lv_max_pg, bernoulli_pg)) %>%
    summarise(
      n_pairs = n(),
      pearson_r = cor(lv_max_pg, bernoulli_pg),
      median_absolute_difference_mmhg = median(abs_difference),
      p95_absolute_difference_mmhg = quantile(abs_difference, 0.95)
    )
}

# Nearest-in-time echo per CPX within +/- window_days, plus derived indices:
# e' = mean(septal, lateral) if both else whichever exists; max wall =
# max(septal, posterior wall).
align_cpx_to_echo <- function(cpx_df, echo_all, window_days) {
  cpx_df %>%
    mutate(MRN = as.character(MRN), cpx_test_date = as.Date(cpx_test_date)) %>%
    left_join(echo_all, by = "MRN", relationship = "many-to-many") %>%
    mutate(
      delta_days = as.numeric(cpx_test_date - echo_date),
      abs_delta_days = abs(delta_days)
    ) %>%
    filter(!is.na(delta_days), abs_delta_days <= window_days) %>%
    # closest echo per test
    group_by(ID, CPX.sequential, cpx_test_date) %>%
    slice_min(order_by = abs_delta_days, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      E_vel = coalesce(med_peak_e_vel, lat_peak_e_vel),
      e_prime_ave = ifelse(
        !is.na(med_peak_e_vel) & !is.na(lat_peak_e_vel),
        (med_peak_e_vel + lat_peak_e_vel) / 2,
        coalesce(med_peak_e_vel, lat_peak_e_vel)
      ),
      lvot_max_gradient = lv_max_pg,
      lv_septal_thickness = ivsd,
      lv_max_wall_thickness = ifelse(
        is.na(ivsd) & is.na(lvpwd),
        NA_real_,
        pmax(ivsd, lvpwd, na.rm = TRUE)
      ),
      Sex = factor(Sex, levels = c(0, 1), labels = c("Male", "Female"))
    )
}

# ── Genetics ──────────────────────────────────────────────────────────────────
derive_genetics_flags <- function(genetics_raw) {
  names(genetics_raw) <- str_replace_all(names(genetics_raw), "\\.", " ")
  genetics_raw %>%
    transmute(
      MRN = as.character(MRN),
      genetics_result_raw = str_squish(str_replace_all(
        as.character(Column1),
        " ",
        " "
      )),
      genetics_result_raw = na_if(genetics_result_raw, ""),
      genetics_phenotype = na_if(str_squish(as.character(`HCM Phenotype`)), ""),
      genetics_result_available = as.integer(!is.na(genetics_result_raw)),
      pathogenic_variant = case_when(
        is.na(genetics_result_raw) ~ 0L,
        str_detect(str_to_upper(genetics_result_raw), "^NONE") ~ 0L,
        str_to_upper(genetics_result_raw) == "ASH" ~ 0L,
        TRUE ~ 1L
      )
    ) %>%
    distinct(MRN, .keep_all = TRUE)
}

# ── Non-HCM comparator screening ──────────────────────────────────────────────
# A no-HCM flag alone is insufficient for a reference group in a tertiary CPET
# registry. Exclude explicit cardiomyopathic, ischemic, valvular, congenital,
# pulmonary-vascular, transplant/device, and inherited arrhythmic diagnoses,
# then require preserved EF, non-hypertrophic wall thickness, and a validated
# resting LV Doppler gradient below the obstruction threshold.
control_exclusion_pattern <- paste(
  c(
    "HCM", "NICM", "ICM", "DCM", "ARVC", "LVNC", "CARDIOM", "MYOPATH",
    "(^|[^A-Z])CM([^A-Z]|$)", "AMYLOID", "SARCOID", "HEART FAILURE",
    "HFR?EF", "HFPEF", "CAD", "CORON", "ISCHEM", "LVH", "LV DILAT",
    "VALV", "AORT", "MITRAL", "TRICUSPID", "PULMONIC", "BAV",
    "(^|[^A-Z])(AVR|MVR|PVR|TVR)([^A-Z]|$)", "ROSS", "CONGEN", "TGA",
    "FONTAN", "TETRAL", "COARCT", "(^|[^A-Z])(VSD|ASD)([^A-Z]|$)",
    "SHUNT", "PULM(ONARY)? HYPERT", "PAH", "TRANSPLANT", "SCD", "ICD",
    "PACEMAKER", "CHANNELOPATH", "CPVT", "DANON", "LAMP2", "LMNA",
    "FLNC", "DESMOPLAKIN", "PLN"
  ),
  collapse = "|"
)

derive_nonhcm_baseline <- function(cpx_nonhcm_echo, num_vars) {
  cpx_nonhcm_echo %>%
    group_by(ID) %>%
    arrange(days_from_baseline) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      across(any_of(num_vars), ~ as.numeric(as.character(.))),
      control_diagnosis_text = str_to_upper(str_squish(paste(
        coalesce(diag_primary, ""),
        coalesce(diag_secondary, ""),
        coalesce(registry_diagnosis_comment, "")
      )))
    ) %>%
    mutate(
      control_no_exclusionary_diagnosis = !str_detect(
        control_diagnosis_text,
        control_exclusion_pattern
      ),
      # Preserved EF and non-hypertrophic walls are required outright. A
      # resting gradient is disqualifying only when it was measured and
      # obstructive: an unmeasured gradient is not evidence of obstruction in a
      # patient who already has LVEF >=50% and maximum wall <1.3 cm, and
      # requiring one removed 80 otherwise-eligible comparators.
      control_structurally_normal_echo = !is.na(ef_modsp4) & ef_modsp4 >= 50 &
        !is.na(lv_max_wall_thickness) & lv_max_wall_thickness < 1.3 &
        (is.na(lvot_max_gradient) | lvot_max_gradient < 30),
      control_qc_eligible = control_no_exclusionary_diagnosis &
        control_structurally_normal_echo
    )
}

# ── Clinical outcomes ─────────────────────────────────────────────────────────
intervention_binary_cols <- c(
  "pre_septal_reduction_surgery",
  "post_septal_reduction_surgery",
  "pre_ablation_surgery",
  "post_ablation_surgery",
  "pre_defibrillator",
  "post_defibrillator",
  "pre_pacemaker",
  "post_pacemaker",
  "pre_heart_transplant",
  "post_heart_transplant"
)

intervention_date_cols <- c(
  "pre_septal_reduction_surgery_date",
  "post_septal_reduction_surgery_date",
  "pre_ablation_surgery_date",
  "post_ablation_surgery_date",
  "pre_defibrillator_date",
  "post_defibrillator_date",
  "pre_pacemaker_date",
  "post_pacemaker_date",
  "pre_heart_transplant_date",
  "post_heart_transplant_date"
)

outcome_binary_cols <- c(
  "death",
  "post_acute_heart_failure",
  "post_chronic_heart_failure",
  "post_heart_transplant",
  "post_afib_flut",
  intervention_binary_cols
)

outcome_date_cols <- c(
  "death_date_clean",
  "post_acute_heart_failure_date",
  "post_chronic_heart_failure_date",
  "post_heart_transplant_date",
  "post_afib_flut_date",
  "last_enc_date",
  intervention_date_cols
)

outcome_yrs_cols <- c(
  "death_yrs",
  "post_acute_heart_failure_yrs",
  "post_chronic_heart_failure_yrs",
  "post_heart_transplant_yrs",
  "post_afib_flut_yrs"
)

# Source outcome fields, one row per source CPX row, with the source-sheet
# composite summaries retained for reference.
prepare_outcomes_table <- function(outcomes_raw) {
  outcomes_raw %>%
    mutate(
      MRN = as.character(MRN),
      cpx_test_date = coerce_excel_date(cpx_test_date)
    ) %>%
    select(any_of(c(
      "MRN",
      "cpx_test_date",
      outcome_binary_cols,
      outcome_date_cols,
      outcome_yrs_cols
    ))) %>%
    mutate(across(any_of(outcome_binary_cols), ~ as_num(.))) %>%
    mutate(across(any_of(outcome_date_cols), coerce_excel_date)) %>%
    mutate(across(any_of(outcome_yrs_cols), ~ as_num(.))) %>%
    # composite HF: event = any component (max); time = earliest (min). coalesce
    # uses Inf as the "no event" sentinel, flipped back to NA after pmin
    mutate(
      hf_composite = pmax(
        coalesce(post_acute_heart_failure, 0),
        coalesce(post_chronic_heart_failure, 0),
        coalesce(post_heart_transplant, 0),
        coalesce(death, 0),
        na.rm = TRUE
      ),
      hf_composite_yrs = pmin(
        coalesce(post_acute_heart_failure_yrs, Inf),
        coalesce(post_chronic_heart_failure_yrs, Inf),
        coalesce(post_heart_transplant_yrs, Inf),
        coalesce(death_yrs, Inf),
        na.rm = TRUE
      ),
      hf_composite_yrs = ifelse(
        is.infinite(hf_composite_yrs),
        NA_real_,
        hf_composite_yrs
      ),
      hf_composite_date_num = pmin(
        coalesce(as.numeric(post_acute_heart_failure_date), Inf),
        coalesce(as.numeric(post_chronic_heart_failure_date), Inf),
        coalesce(as.numeric(post_heart_transplant_date), Inf),
        coalesce(as.numeric(death_date_clean), Inf),
        na.rm = TRUE
      ),
      hf_composite_date_num = ifelse(
        is.infinite(hf_composite_date_num),
        NA_real_,
        hf_composite_date_num
      ),
      hf_or_death_no_transplant = pmax(
        coalesce(post_acute_heart_failure, 0),
        coalesce(post_chronic_heart_failure, 0),
        coalesce(death, 0),
        na.rm = TRUE
      ),
      hf_or_death_no_transplant_date_num = pmin(
        coalesce(as.numeric(post_acute_heart_failure_date), Inf),
        coalesce(as.numeric(post_chronic_heart_failure_date), Inf),
        coalesce(as.numeric(death_date_clean), Inf),
        na.rm = TRUE
      ),
      hf_or_death_no_transplant_date_num = ifelse(
        is.infinite(hf_or_death_no_transplant_date_num),
        NA_real_,
        hf_or_death_no_transplant_date_num
      )
    )
}

# Per-test intervention index (MRN + date key).
summarise_outcome_index <- function(outcomes) {
  outcomes %>%
    select(any_of(c(
      "MRN",
      "cpx_test_date",
      intervention_binary_cols,
      intervention_date_cols
    ))) %>%
    group_by(MRN, cpx_test_date) %>%
    summarise(
      across(any_of(intervention_binary_cols), max_numeric_or_na),
      across(any_of(intervention_date_cols), min_date_or_na),
      .groups = "drop"
    )
}

# One row per patient containing source dates only. The pre/post procedure
# flags are not used to establish temporal order because their labels are
# inconsistent across index CPX rows; event status is rebuilt from the actual
# date relative to each patient's selected baseline CPX (add_outcome_endpoints).
summarise_outcome_dates <- function(outcomes) {
  outcomes %>%
    group_by(MRN) %>%
    summarise(
      acute_hf_source_date = min_date_or_na(post_acute_heart_failure_date),
      chronic_hf_source_date = min_date_or_na(post_chronic_heart_failure_date),
      transplant_source_date = min_date_or_na(c(
        pre_heart_transplant_date,
        post_heart_transplant_date
      )),
      death_source_date = min_date_or_na(death_date_clean),
      myectomy_source_date = min_date_or_na(c(
        pre_septal_reduction_surgery_date,
        post_septal_reduction_surgery_date
      )),
      last_enc_date = max_date_or_na(last_enc_date),
      .groups = "drop"
    )
}

# Primary composite: acute HF event, heart transplantation, or death, each
# counted only when dated after the baseline CPX. Chronic-HF diagnosis is
# carried separately and is not an acute clinical endpoint.
add_outcome_endpoints <- function(baseline_df, outcomes_dedup, outcomes_index) {
  n_pre_join <- nrow(baseline_df)
  out <- baseline_df %>%
    left_join(outcomes_dedup, by = "MRN") %>%
    left_join(outcomes_index, by = c("MRN", "cpx_test_date")) %>%
    mutate(
      cpx_test_date = as.Date(cpx_test_date),
      acute_hf_date_num = ifelse(
        !is.na(acute_hf_source_date) & acute_hf_source_date > cpx_test_date,
        as.numeric(acute_hf_source_date),
        NA_real_
      ),
      chronic_hf_date_num = ifelse(
        !is.na(chronic_hf_source_date) & chronic_hf_source_date > cpx_test_date,
        as.numeric(chronic_hf_source_date),
        NA_real_
      ),
      transplant_date_num = ifelse(
        !is.na(transplant_source_date) & transplant_source_date > cpx_test_date,
        as.numeric(transplant_source_date),
        NA_real_
      ),
      death_date_num = ifelse(
        !is.na(death_source_date) & death_source_date > cpx_test_date,
        as.numeric(death_source_date),
        NA_real_
      ),
      myectomy_date_num = ifelse(
        !is.na(myectomy_source_date) & myectomy_source_date > cpx_test_date,
        as.numeric(myectomy_source_date),
        NA_real_
      ),
      hf_composite_date_num = pmin(
        coalesce(acute_hf_date_num, Inf),
        coalesce(transplant_date_num, Inf),
        coalesce(death_date_num, Inf),
        na.rm = TRUE
      ),
      hf_composite_date_num = ifelse(
        is.infinite(hf_composite_date_num),
        NA_real_,
        hf_composite_date_num
      ),
      hf_composite = as.integer(!is.na(hf_composite_date_num)),
      hf_composite_yrs = (hf_composite_date_num - as.numeric(cpx_test_date)) /
        365.25,
      hf_or_death_no_transplant_date_num = pmin(
        coalesce(acute_hf_date_num, Inf),
        coalesce(death_date_num, Inf),
        na.rm = TRUE
      ),
      hf_or_death_no_transplant_date_num = ifelse(
        is.infinite(hf_or_death_no_transplant_date_num),
        NA_real_,
        hf_or_death_no_transplant_date_num
      ),
      hf_or_death_no_transplant = as.integer(
        !is.na(hf_or_death_no_transplant_date_num)
      ),
      censor_time_yrs = as.numeric(last_enc_date - cpx_test_date) / 365.25,
      # Follow-up adequate for a time-to-event analysis, independent of which
      # diastolic parameter is the exposure.
      outcomes_followup_ok = (hf_composite == 1 & hf_composite_yrs > 0) |
        (hf_composite == 0 & !is.na(censor_time_yrs) & censor_time_yrs > 0),
      outcomes_analytic = cross_sectional_eligible & outcomes_followup_ok
    )
  stopifnot(nrow(out) == n_pre_join)
  out
}

# ── Medications and comorbidity adjustment flags ──────────────────────────────
# Reconcile the summary "Betablocker" flag with the individual drug columns:
# any drug = 1 -> on; all drugs explicitly 0 -> off; else NA. Mismatches are
# flagged in bb_discordant. Non-DHP calcium-channel blockers likewise.
bb_drug_cols <- c(
  "carvedilol",
  "metoprolol",
  "bisoprolol",
  "atenolol",
  "nebivolol",
  "propranolol",
  "nadolol",
  "sotalol",
  "other.lol"
)
non_dhp_ccb_cols <- c("diltiazem", "verapamil")
medication_flag_cols <- c(
  "Betablocker",
  bb_drug_cols,
  "Calcium.C.blocker",
  non_dhp_ccb_cols,
  "disopyramide",
  "ACEI/ARB",
  "Diuretics",
  "Statin"
)

derive_medication_flags <- function(medications_raw) {
  medications_raw %>%
    mutate(
      MRN = as.character(MRN),
      cpx_test_date = coerce_excel_date(cpx_test_date)
    ) %>%
    select(any_of(c("MRN", "cpx_test_date", medication_flag_cols))) %>%
    mutate(across(any_of(medication_flag_cols), ~ as_num(.))) %>%
    group_by(MRN, cpx_test_date) %>%
    summarise(
      across(any_of(medication_flag_cols), max_numeric_or_na),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      med_data_available = as.integer(any(
        !is.na(c_across(any_of(medication_flag_cols)))
      )),
      bb_drug_positive_n = sum(c_across(any_of(bb_drug_cols)) == 1, na.rm = TRUE),
      bb_drug_nonmissing_n = sum(!is.na(c_across(any_of(bb_drug_cols)))),
      bb_drug_all_zero = bb_drug_nonmissing_n == length(bb_drug_cols) &&
        sum(c_across(any_of(bb_drug_cols)) == 0, na.rm = TRUE) ==
          length(bb_drug_cols),
      bb_any = case_when(
        bb_drug_positive_n > 0 ~ 1,
        !is.na(Betablocker) & Betablocker == 1 ~ 1,
        !is.na(Betablocker) & Betablocker == 0 & bb_drug_positive_n == 0 ~ 0,
        is.na(Betablocker) & bb_drug_all_zero ~ 0,
        TRUE ~ NA_real_
      ),
      bb_discordant = as.integer(
        !is.na(Betablocker) &&
          ((Betablocker == 0 & bb_drug_positive_n > 0) ||
            (Betablocker == 1 & bb_drug_all_zero))
      ),
      non_dhp_ccb_any = case_when(
        sum(c_across(any_of(non_dhp_ccb_cols)) == 1, na.rm = TRUE) > 0 ~ 1,
        sum(!is.na(c_across(any_of(non_dhp_ccb_cols)))) ==
          length(non_dhp_ccb_cols) &&
          sum(c_across(any_of(non_dhp_ccb_cols)) == 0, na.rm = TRUE) ==
            length(non_dhp_ccb_cols) ~ 0,
        TRUE ~ NA_real_
      ),
      ndhp_ccb_any = non_dhp_ccb_any
    ) %>%
    ungroup() %>%
    select(-bb_drug_positive_n, -bb_drug_nonmissing_n, -bb_drug_all_zero)
}

derive_comorbidity_adjustment_flags <- function(comorbidities_table) {
  comorbidities_table %>%
    mutate(
      MRN = as.character(MRN),
      cpx_test_date = coerce_excel_date(cpx_test_date),
      dm_pre_test = as_num(dm_pre_test),
      htn_pre_test = as_num(htn_pre_test)
    ) %>%
    select(any_of(c("MRN", "cpx_test_date", "dm_pre_test", "htn_pre_test"))) %>%
    group_by(MRN, cpx_test_date) %>%
    summarise(
      dm_pre_test = max_numeric_or_na(dm_pre_test),
      htn_pre_test = max_numeric_or_na(htn_pre_test),
      .groups = "drop"
    )
}

left_join_preserving_rows <- function(x, y, by) {
  n_pre <- nrow(x)
  out <- left_join(x, y, by = by)
  stopifnot(nrow(out) == n_pre)
  out
}

# ── Longitudinal subcohort ────────────────────────────────────────────────────
# Exclude CPETs on or after surgical myectomy (actual procedure date), then
# require >= 2 distinct test dates spanning >= 0.5 years.
censor_longitudinal_post_myectomy <- function(long_df, outcomes_dedup) {
  long_df %>%
    left_join(
      outcomes_dedup %>% select(MRN, myectomy_source_date),
      by = "MRN"
    ) %>%
    filter(is.na(myectomy_source_date) | as.Date(cpx_test_date) < myectomy_source_date) %>%
    group_by(ID) %>%
    filter(
      n_distinct(cpx_test_date) >= 2,
      max(time_yrs, na.rm = TRUE) - min(time_yrs, na.rm = TRUE) >= 0.5
    ) %>%
    ungroup()
}

# ── ASE diastolic grading + filling-pressure class ────────────────────────────
# 3 primary indices vs ASE cutoffs (E/e'>14, LAVi>34, TRV>280) -> combo code/class.
# fp elevated if: >=2 primary abnl, OR 1 primary + >=1 supportive, OR (<=1 primary
# measured) + >=2 supportive (restrictive E/A, short DT, low PV S/D)
# eprime_bands: age-specific lower reference limits for average e-prime, in
# cm/s, as a list of list(max_age, limit) applied in order. Patients younger
# than min_classifiable_age are left unclassified rather than pushed into the
# youngest band, because the reference limits do not cover them.
add_ase_grading <- function(baseline_df,
                            eprime_bands = list(
                              list(max_age = 39.999, limit = 9.0),
                              list(max_age = 65.0, limit = 7.0),
                              list(max_age = 200, limit = 6.5)
                            ),
                            min_classifiable_age = 20) {
  eprime_limit_for <- function(age) {
    age <- as_num(age)
    out <- rep(NA_real_, length(age))
    for (band in rev(eprime_bands)) {
      out[!is.na(age) & age <= band$max_age] <- band$limit
    }
    out[is.na(age) | age < min_classifiable_age] <- NA_real_
    out
  }
  baseline_df %>%
    mutate(
      eprime_age_limit = eprime_limit_for(age),
      eprime_classifiable = !is.na(eprime_age_limit) &
        !is.na(as_num(e_prime_ave)),
      eprime_below_age_limit = eprime_classifiable &
        as_num(e_prime_ave) < eprime_age_limit,
      pv_sd_ratio = ifelse(
        !is.na(pulm_sys_vel) & !is.na(pulm_dias_vel) & pulm_dias_vel > 0,
        pulm_sys_vel / pulm_dias_vel,
        NA_real_
      ),
      abn_ee = !is.na(e_e_ave) & e_e_ave > 14,
      abn_lavi = !is.na(la_vol_index) & la_vol_index > 34,
      abn_trv = !is.na(tr_max_vel) & tr_max_vel > 280,
      n_primary_available = (!is.na(e_e_ave)) +
        (!is.na(la_vol_index)) +
        (!is.na(tr_max_vel)),
      n_primary_abnormal = abn_ee + abn_lavi + abn_trv,
      hcm_combo_code = case_when(
        n_primary_available >= 1 ~ paste0(
          as.integer(abn_ee),
          as.integer(abn_lavi),
          as.integer(abn_trv)
        ),
        TRUE ~ NA_character_
      ),
      hcm_combo_class = case_when(
        hcm_combo_code == "000" ~ "None abnormal",
        hcm_combo_code == "100" ~ "E/e' only",
        hcm_combo_code == "010" ~ "LAVi only",
        hcm_combo_code == "001" ~ "TRVmax only",
        hcm_combo_code == "110" ~ "E/e' + LAVi",
        hcm_combo_code == "101" ~ "E/e' + TRVmax",
        hcm_combo_code == "011" ~ "LAVi + TRVmax",
        hcm_combo_code == "111" ~ "All 3 abnormal",
        TRUE ~ NA_character_
      ),
      hcm_combo_class = factor(
        hcm_combo_class,
        levels = c(
          "None abnormal",
          "E/e' only",
          "LAVi only",
          "TRVmax only",
          "E/e' + LAVi",
          "E/e' + TRVmax",
          "LAVi + TRVmax",
          "All 3 abnormal"
        )
      ),
      hcm_combo_n_abnormal = ifelse(
        n_primary_available >= 1,
        n_primary_abnormal,
        NA_integer_
      ),
      hcm_combo_abnormal_group = factor(
        hcm_combo_n_abnormal,
        levels = 0:3,
        labels = c("0", "1", "2", "3")
      ),
      restrictive_pattern = !is.na(mv_e_a) & mv_e_a >= 2,
      # Source values are stored in seconds; apply the 150-ms threshold in
      # consistent units. This legacy classifier is supplemental only.
      short_dt = !is.na(mv_dec_time) & (mv_dec_time * 1000) < 150,
      low_pv_sd = !is.na(pv_sd_ratio) & pv_sd_ratio < 1,
      n_supportive_elevated = restrictive_pattern + short_dt + low_pv_sd,
      elevated_by_primary = n_primary_abnormal >= 2,
      elevated_by_mixed = n_primary_abnormal == 1 & n_supportive_elevated >= 1,
      elevated_by_supportive_only = n_primary_available <= 1 &
        n_supportive_elevated >= 2,
      fp_class = case_when(
        elevated_by_primary ~ "Elevated",
        elevated_by_mixed ~ "Elevated",
        elevated_by_supportive_only ~ "Elevated",
        TRUE ~ "Not Elevated"
      ),
      fp_class = factor(fp_class, levels = c("Not Elevated", "Elevated"))
    )
}

# ── Parameter-specific eligibility ────────────────────────────────────────────
# Base eligibility is what every analysis needs regardless of which diastolic
# parameter is the exposure: complete age, sex, BMI and at least one CPET
# outcome. Each parameter then adds only its own measurement, so requiring one
# parameter never removes patients from another parameter's analysis.
# See 8_Docs/ANALYSIS_PLAN_2026-09.md.
primary_parameters <- c("e_e_ave", "la_vol_index", "tr_max_vel")
primary_parameter_labels <- c(
  e_e_ave = "E/e' (average)",
  la_vol_index = "LAVi",
  tr_max_vel = "TRVmax"
)

add_parameter_eligibility <- function(df) {
  out <- df %>%
    dplyr::mutate(
      analysis_base_eligible = !is.na(as_num(age)) &
        Sex %in% c("Male", "Female") &
        !is.na(as_num(BMI)) &
        (!is.na(as_num(VO2_FRIEND2_PP)) |
          !is.na(as_num(VeVco2_slope)) |
          !is.na(as_num(HRR)))
    )
  for (p in primary_parameters) {
    out[[paste0("eligible_", p)]] <- out$analysis_base_eligible &
      !is.na(as_num(out[[p]]))
  }
  out
}

# Repeated-test subcohort for a given set of patient IDs: >1 aligned test, an
# observed peak VO2 and time, and the longitudinal time axis. Used both for the
# primary longitudinal cohort and for the parameter-specific counts, so the two
# cannot drift apart.
build_longitudinal_subcohort <- function(cpx_echo, ids, num_vars) {
  repeat_ids <- cpx_echo %>%
    dplyr::semi_join(ids, by = "ID") %>%
    dplyr::count(ID, name = "num_aligned_tests") %>%
    dplyr::filter(num_aligned_tests > 1)

  cpx <- cpx_echo %>% dplyr::semi_join(repeat_ids, by = "ID")

  frame <- cpx %>%
    dplyr::filter(!is.na(VO2_FRIEND2_PP), !is.na(days_from_baseline)) %>%
    dplyr::mutate(
      dplyr::across(dplyr::any_of(num_vars), ~ as.numeric(as.character(.))),
      time_yrs = days_from_baseline / 365.25,
      cohort = "HCM"
    )

  list(ids = repeat_ids, cpx = cpx, frame = frame)
}

# Cross-sectional, longitudinal, and outcome counts each parameter would have
# in its own available-case sample, alongside the common-cohort counts.
# Aggregate only; writes nothing.
summarise_parameter_specific_counts <- function(baseline_df, cpx_echo, outcomes_dedup, num_vars) {
  count_row <- function(label, keep) {
    ids <- baseline_df %>% dplyr::filter(keep) %>% dplyr::distinct(ID)
    lng <- build_longitudinal_subcohort(cpx_echo, ids, num_vars)$frame %>%
      censor_longitudinal_post_myectomy(outcomes_dedup)
    out_keep <- keep & baseline_df$outcomes_followup_ok
    tibble::tibble(
      sample = label,
      n_cross_sectional = sum(keep, na.rm = TRUE),
      n_peak_vo2 = sum(keep & !is.na(as_num(baseline_df$VO2_FRIEND2_PP)), na.rm = TRUE),
      n_vevco2 = sum(keep & !is.na(as_num(baseline_df$VeVco2_slope)), na.rm = TRUE),
      n_longitudinal_patients = dplyr::n_distinct(lng$ID),
      n_longitudinal_observations = nrow(lng),
      n_outcomes_patients = sum(out_keep, na.rm = TRUE),
      n_outcomes_events = sum(out_keep & baseline_df$hf_composite == 1, na.rm = TRUE)
    )
  }

  rows <- lapply(primary_parameters, function(p) {
    count_row(primary_parameter_labels[[p]], baseline_df[[paste0("eligible_", p)]])
  })
  dplyr::bind_rows(
    count_row("Common cohort (E/e' + LAVi)", baseline_df$cross_sectional_eligible),
    dplyr::bind_rows(rows),
    count_row("Base (age/sex/BMI + >=1 CPET outcome)", baseline_df$analysis_base_eligible)
  )
}

# ── Missing-indicator handling for adjustment covariates ──────────────────────
# Owner decision, card 1.5 of 8_Docs/DECISION_CARDS_2026-09.md. The expanded
# adjustment set is only ~40% complete, so requiring complete cases on it
# discards more than half the analytic sample and makes covariate availability,
# not the exposure, define the population.
#
# Two devices, both of which keep every patient:
#
#   * A binary covariate becomes a three-level factor (No / Yes / Unknown), so
#     "not recorded" is its own stratum instead of a dropped patient. This is
#     what Dr. Haddad described as the honest handling for a yes/no field: we
#     do not pretend to know which way an unrecorded value falls.
#
#   * A continuous covariate becomes a presence indicator plus a gated value,
#     present * value. Because the value is multiplied by presence, whatever
#     constant fills the missing entries contributes only through the indicator
#     and cancels out of the exposure's coefficient: fitting b*(d*x) + g*d with
#     x filled by c where d = 0 gives the same exposure estimate for any c.
#     Mean imputation of a continuous covariate would instead assert a value
#     that was never measured.
#
# Returns the augmented data and the covariate names to put in the model.
# binary_style: "factor" emits one three-level factor per variable, which is
# the natural form for lm(), lme() and coxph(). "dummy" emits two 0/1 numeric
# columns instead, <v>_yes and <v>_unknown, with "No" as the reference. The two
# are mathematically identical parameterisations; "dummy" exists because
# fit_gamm_binary_panel() coerces every covariate except Sex to numeric, which
# would turn a factor into all-NA and silently drop every row.
add_missing_indicator_covariates <- function(df,
                                             binary_vars = character(0),
                                             gated_numeric_vars = character(0),
                                             binary_style = c("factor", "dummy")) {
  binary_style <- match.arg(binary_style)
  new_covariates <- character(0)

  for (v in binary_vars) {
    stopifnot(v %in% names(df))
    observed <- setdiff(unique(df[[v]]), NA)
    if (!all(observed %in% c(0, 1))) {
      stop(
        sprintf(
          "add_missing_indicator_covariates(): '%s' is not 0/1 coded (values: %s)",
          v, paste(sort(observed), collapse = ", ")
        ),
        call. = FALSE
      )
    }
    if (binary_style == "factor") {
      nm <- paste0(v, "_3lvl")
      df[[nm]] <- factor(
        dplyr::case_when(
          is.na(df[[v]]) ~ "Unknown",
          df[[v]] == 1 ~ "Yes",
          TRUE ~ "No"
        ),
        levels = c("No", "Yes", "Unknown")
      )
      new_covariates <- c(new_covariates, nm)
    } else {
      yes_nm <- paste0(v, "_yes")
      unknown_nm <- paste0(v, "_unknown")
      df[[yes_nm]] <- as.integer(!is.na(df[[v]]) & df[[v]] == 1)
      df[[unknown_nm]] <- as.integer(is.na(df[[v]]))
      new_covariates <- c(new_covariates, yes_nm, unknown_nm)
    }
  }

  for (v in gated_numeric_vars) {
    stopifnot(v %in% names(df))
    present_nm <- paste0(v, "_present")
    gated_nm <- paste0(v, "_gated")
    df[[present_nm]] <- as.integer(!is.na(df[[v]]))
    # The fill value is arbitrary precisely because it is gated away.
    df[[gated_nm]] <- dplyr::if_else(is.na(df[[v]]), 0, as.numeric(df[[v]])) *
      df[[present_nm]]
    new_covariates <- c(new_covariates, gated_nm, present_nm)
  }

  list(data = df, covariates = new_covariates)
}

# A covariate that takes one value in a given model frame carries no
# information and breaks lm()/coxph() contrasts, which is a real possibility
# for an "Unknown" level inside a small subset. Drop those rather than let the
# fit fail or return an NA coefficient, and report what was dropped.
drop_degenerate_covariates <- function(df, covariates, verbose = FALSE) {
  keep <- vapply(covariates, function(v) {
    if (!v %in% names(df)) return(FALSE)
    length(setdiff(unique(df[[v]]), NA)) >= 2
  }, logical(1))
  if (verbose && any(!keep)) {
    message(sprintf(
      "  dropped single-valued covariate(s) in this frame: %s",
      paste(covariates[!keep], collapse = ", ")
    ))
  }
  covariates[keep]
}
