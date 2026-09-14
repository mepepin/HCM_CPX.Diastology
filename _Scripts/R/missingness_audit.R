# Missingness and cohort-selection audit helpers.
#
# These functions are intentionally independent of manuscript rendering. They
# write de-identified aggregate tables only and never export MRNs or patient IDs.

audit_as_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

audit_is_observed <- function(x) {
  if (inherits(x, "Date")) return(!is.na(x))
  if (is.character(x) || is.factor(x)) {
    return(!is.na(x) & trimws(as.character(x)) != "")
  }
  !is.na(x)
}

build_missingness_table <- function(data, variable_spec, cohort_label) {
  purrr::pmap_dfr(variable_spec, function(variable, label, role, type) {
    if (!variable %in% names(data)) return(tibble::tibble())
    observed <- audit_is_observed(data[[variable]])
    tibble::tibble(
      Cohort = cohort_label,
      Variable = variable,
      Label = label,
      Role = role,
      Type = type,
      N_total = nrow(data),
      N_observed = sum(observed),
      N_missing = sum(!observed),
      Percent_missing = 100 * mean(!observed)
    )
  })
}

build_missingness_by_group <- function(
  data,
  variable_spec,
  group_var,
  group_labels,
  cohort_label
) {
  if (!group_var %in% names(data)) return(tibble::tibble())
  purrr::map_dfr(names(group_labels), function(group_value) {
    group_data <- data[
      !is.na(data[[group_var]]) &
        as.character(data[[group_var]]) == group_value,
      ,
      drop = FALSE
    ]
    build_missingness_table(
      group_data,
      variable_spec,
      paste0(cohort_label, ": ", unname(group_labels[[group_value]]))
    )
  })
}

audit_binary_smd <- function(p1, p0) {
  denominator <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  (p1 - p0) / denominator
}

audit_continuous_smd <- function(x1, x0) {
  pooled_sd <- sqrt((stats::var(x1) + stats::var(x0)) / 2)
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  (mean(x1) - mean(x0)) / pooled_sd
}

build_selection_smd_table <- function(
  data,
  group_var,
  variable_spec,
  comparison_label,
  group1_label,
  group0_label
) {
  stopifnot(group_var %in% names(data))
  group_value <- audit_as_numeric(data[[group_var]])

  purrr::pmap_dfr(variable_spec, function(variable, label, type) {
    if (!variable %in% names(data)) return(tibble::tibble())
    raw_value <- data[[variable]]
    value <- if (type == "binary") {
      audit_as_numeric(raw_value)
    } else {
      audit_as_numeric(raw_value)
    }
    observed <- !is.na(value) & !is.na(group_value) & group_value %in% c(0, 1)
    x1 <- value[observed & group_value == 1]
    x0 <- value[observed & group_value == 0]

    if (type == "binary") {
      summary1 <- if (length(x1)) mean(x1) else NA_real_
      summary0 <- if (length(x0)) mean(x0) else NA_real_
      smd <- audit_binary_smd(summary1, summary0)
      statistic <- "Proportion"
    } else {
      summary1 <- if (length(x1)) mean(x1) else NA_real_
      summary0 <- if (length(x0)) mean(x0) else NA_real_
      smd <- if (length(x1) > 1 && length(x0) > 1) {
        audit_continuous_smd(x1, x0)
      } else {
        NA_real_
      }
      statistic <- "Mean"
    }

    tibble::tibble(
      Comparison = comparison_label,
      Variable = variable,
      Label = label,
      Type = type,
      Statistic = statistic,
      Group_1 = group1_label,
      Group_1_N_observed = length(x1),
      Group_1_value = summary1,
      Group_0 = group0_label,
      Group_0_N_observed = length(x0),
      Group_0_value = summary0,
      SMD = smd,
      Absolute_SMD = abs(smd)
    )
  })
}

run_missingness_selection_audit <- function(
  baseline_df,
  long_df,
  output_dir
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  variable_spec <- tibble::tribble(
    ~variable, ~label, ~role, ~type,
    "age", "Age", "Primary covariate", "continuous",
    "Sex", "Sex", "Primary covariate", "categorical",
    "BMI", "BMI", "Primary covariate", "continuous",
    "e_e_ave", "Average E/e'", "Cohort-defining exposure", "continuous",
    "la_vol_index", "LAVi", "Cohort-defining exposure", "continuous",
    "tr_max_vel", "TRVmax", "Available-case exposure", "continuous",
    "VO2_FRIEND2_PP", "Peak VO2 (FRIEND 2.0 % predicted)", "Outcome", "continuous",
    "VeVco2_slope", "VE/VCO2 slope", "Outcome", "continuous",
    "HRR", "Heart-rate recovery", "Secondary outcome", "continuous",
    "lv_max_wall_thickness", "Maximum LV wall thickness", "Structural sensitivity covariate", "continuous",
    "lv_septal_thickness", "LV septal thickness", "Structural sensitivity covariate", "continuous",
    "lvot_max_gradient", "Resting maximum LVOT gradient", "Structural sensitivity covariate", "continuous",
    "apical_hcm", "Apical HCM morphology", "Morphology sensitivity covariate", "binary",
    "HCM_Phenotype", "HCM phenotype", "Descriptive phenotype", "categorical",
    "bb_any", "Beta-blocker use", "Medication sensitivity covariate", "binary",
    "ndhp_ccb_any", "Non-DHP calcium-channel blocker use", "Medication sensitivity covariate", "binary",
    "htn_pre_test", "Hypertension", "Clinical sensitivity covariate", "binary",
    "dm_pre_test", "Diabetes", "Clinical sensitivity covariate", "binary",
    "ef_modsp4", "LVEF", "Descriptive echocardiographic variable", "continuous",
    "e_prime_ave", "Average e'", "Descriptive diastolic variable", "continuous",
    "hf_composite", "Primary clinical composite status", "Outcome", "binary",
    "last_enc_date", "Last encounter date", "Follow-up", "date"
  )

  analytic_df <- baseline_df |>
    dplyr::filter(cross_sectional_eligible)
  outcomes_df <- baseline_df |>
    dplyr::filter(outcomes_analytic)
  longitudinal_ids <- long_df |>
    dplyr::distinct(ID)
  longitudinal_baseline_df <- baseline_df |>
    dplyr::semi_join(longitudinal_ids, by = "ID")
  trv_measured_df <- analytic_df |>
    dplyr::filter(!is.na(tr_max_vel))

  missingness_overall <- dplyr::bind_rows(
    build_missingness_table(baseline_df, variable_spec, "Phenotypic HCM entry"),
    build_missingness_table(analytic_df, variable_spec, "E/e' + LAVi analytic cohort"),
    build_missingness_table(trv_measured_df, variable_spec, "TRVmax-measured subset"),
    build_missingness_table(longitudinal_baseline_df, variable_spec, "Longitudinal patient subset"),
    build_missingness_table(outcomes_df, variable_spec, "Clinical-outcomes cohort")
  )

  utils::write.csv(
    missingness_overall,
    file.path(output_dir, "Table_Missingness_ByVariableAndCohort.csv"),
    row.names = FALSE
  )

  missingness_by_event <- build_missingness_by_group(
    analytic_df,
    variable_spec,
    "hf_composite",
    c(`0` = "no primary event", `1` = "primary event"),
    "E/e' + LAVi analytic cohort"
  )
  utils::write.csv(
    missingness_by_event,
    file.path(output_dir, "Table_Missingness_ByEventStatus.csv"),
    row.names = FALSE
  )

  requirement_flags <- baseline_df |>
    dplyr::transmute(
      ID,
      demographics_complete = !is.na(audit_as_numeric(age)) &
        Sex %in% c("Male", "Female") &
        !is.na(audit_as_numeric(BMI)),
      ee_measured = !is.na(audit_as_numeric(e_e_ave)),
      lavi_measured = !is.na(audit_as_numeric(la_vol_index)),
      ee_lavi_measured = ee_measured & lavi_measured,
      cpet_outcome_available = !is.na(audit_as_numeric(VO2_FRIEND2_PP)) |
        !is.na(audit_as_numeric(VeVco2_slope)) |
        !is.na(audit_as_numeric(HRR)),
      cross_sectional_eligible,
      trv_measured = cross_sectional_eligible &
        !is.na(audit_as_numeric(tr_max_vel)),
      outcomes_analytic
    )

  cohort_flow <- tibble::tibble(
    Stage_ID = c(
      "hcm_entry",
      "demographics_complete",
      "ee_measured",
      "lavi_measured",
      "cross_sectional",
      "trv_measured",
      "outcomes",
      "longitudinal"
    ),
    Parent_stage_ID = c(
      NA_character_,
      "hcm_entry",
      "demographics_complete",
      "ee_measured",
      "lavi_measured",
      "cross_sectional",
      "cross_sectional",
      "cross_sectional"
    ),
    Branch = c(
      "Upstream analytic flow",
      "Upstream analytic flow",
      "Upstream analytic flow",
      "Upstream analytic flow",
      "Upstream analytic flow",
      "Nested cross-sectional subset",
      "Nested cross-sectional subset",
      "Nested cross-sectional subset"
    ),
    Stage = c(
      "Phenotypic HCM entry",
      "Complete age, sex, and BMI",
      "+ E/e' measured",
      "+ LAVi measured",
      "+ at least one CPET outcome: cross-sectional cohort",
      "TRVmax measured within cross-sectional cohort",
      "Clinical-outcomes subcohort",
      "Longitudinal subcohort after intervention censoring"
    ),
    N_patients = c(
      nrow(baseline_df),
      sum(requirement_flags$demographics_complete),
      sum(requirement_flags$demographics_complete & requirement_flags$ee_measured),
      sum(requirement_flags$demographics_complete &
        requirement_flags$ee_measured & requirement_flags$lavi_measured),
      sum(requirement_flags$cross_sectional_eligible),
      sum(requirement_flags$trv_measured),
      sum(requirement_flags$outcomes_analytic),
      dplyr::n_distinct(long_df$ID)
    ),
    Parent_N = c(
      NA_integer_,
      nrow(baseline_df),
      sum(requirement_flags$demographics_complete),
      sum(requirement_flags$demographics_complete & requirement_flags$ee_measured),
      sum(requirement_flags$demographics_complete &
        requirement_flags$ee_measured & requirement_flags$lavi_measured),
      sum(requirement_flags$cross_sectional_eligible),
      sum(requirement_flags$cross_sectional_eligible),
      sum(requirement_flags$cross_sectional_eligible)
    )
  ) |>
    dplyr::mutate(Excluded_from_parent = Parent_N - N_patients)
  utils::write.csv(
    cohort_flow,
    file.path(output_dir, "Table_AnalyticCohort_RequirementFlow.csv"),
    row.names = FALSE
  )

  comparison_spec <- tibble::tribble(
    ~variable, ~label, ~type,
    "age", "Age", "continuous",
    "female", "Female sex", "binary",
    "BMI", "BMI", "continuous",
    "e_e_ave", "Average E/e'", "continuous",
    "la_vol_index", "LAVi", "continuous",
    "apical_hcm", "Apical HCM", "binary",
    "lv_max_wall_thickness", "Maximum LV wall thickness", "continuous",
    "lv_septal_thickness", "LV septal thickness", "continuous",
    "lvot_max_gradient", "Resting maximum LVOT gradient", "continuous",
    "obstructive_hcm", "Resting LVOT gradient >=30 mm Hg", "binary",
    "VO2_FRIEND2_PP", "Peak VO2 (FRIEND 2.0 % predicted)", "continuous",
    "VeVco2_slope", "VE/VCO2 slope", "continuous",
    "HRR", "Heart-rate recovery", "continuous",
    "hf_composite", "Primary clinical composite event", "binary"
  )

  comparison_df <- baseline_df |>
    dplyr::mutate(
      female = dplyr::case_when(
        Sex == "Female" ~ 1,
        Sex == "Male" ~ 0,
        TRUE ~ NA_real_
      ),
      obstructive_hcm = dplyr::case_when(
        is.na(audit_as_numeric(lvot_max_gradient)) ~ NA_real_,
        audit_as_numeric(lvot_max_gradient) >= 30 ~ 1,
        TRUE ~ 0
      ),
      analytic_included = as.integer(cross_sectional_eligible),
      ee_lavi_measured = as.integer(!is.na(audit_as_numeric(e_e_ave)) &
        !is.na(audit_as_numeric(la_vol_index))),
      trv_measured = as.integer(!is.na(audit_as_numeric(tr_max_vel)))
    )

  selection_smd <- dplyr::bind_rows(
    build_selection_smd_table(
      comparison_df,
      "ee_lavi_measured",
      comparison_spec,
      "E/e' + LAVi measurement selection",
      "Both measured",
      "Either missing"
    ),
    build_selection_smd_table(
      comparison_df,
      "analytic_included",
      comparison_spec,
      "Final cross-sectional cohort selection",
      "Included",
      "Excluded"
    ),
    build_selection_smd_table(
      comparison_df |> dplyr::filter(cross_sectional_eligible),
      "trv_measured",
      comparison_spec,
      "TRVmax measurement within E/e' + LAVi cohort",
      "TRVmax measured",
      "TRVmax missing"
    )
  )
  utils::write.csv(
    selection_smd,
    file.path(output_dir, "Table_CohortSelection_StandardizedDifferences.csv"),
    row.names = FALSE
  )

  # Parsimonious diagnostic model for whether TRVmax was measured. This is a
  # selection audit only; it does not establish a missing-at-random mechanism.
  trv_model_df <- analytic_df |>
    dplyr::transmute(
      trv_measured = as.integer(!is.na(tr_max_vel)),
      age = audit_as_numeric(age),
      Sex = factor(Sex, levels = c("Male", "Female")),
      BMI = audit_as_numeric(BMI),
      e_e_z = as.numeric(scale(audit_as_numeric(e_e_ave))),
      lavi_z = as.numeric(scale(audit_as_numeric(la_vol_index))),
      wall_thickness_z = as.numeric(scale(audit_as_numeric(lv_max_wall_thickness)))
    ) |>
    tidyr::drop_na()

  trv_model <- stats::glm(
    trv_measured ~ age + Sex + BMI + e_e_z + lavi_z + wall_thickness_z,
    data = trv_model_df,
    family = stats::binomial()
  )
  trv_model_table <- broom::tidy(
    trv_model,
    exponentiate = TRUE,
    conf.int = TRUE
  ) |>
    dplyr::transmute(
      Term = term,
      Odds_ratio = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      P_value = p.value,
      N = stats::nobs(trv_model),
      TRVmax_measured = sum(trv_model_df$trv_measured == 1),
      TRVmax_missing = sum(trv_model_df$trv_measured == 0)
    )
  utils::write.csv(
    trv_model_table,
    file.path(output_dir, "Table_TRVmax_Measurement_LogisticModel.csv"),
    row.names = FALSE
  )

  trv_model_fit <- tibble::tibble(
    N = stats::nobs(trv_model),
    TRVmax_measured = sum(trv_model_df$trv_measured == 1),
    TRVmax_missing = sum(trv_model_df$trv_measured == 0),
    Null_deviance = trv_model$null.deviance,
    Residual_deviance = trv_model$deviance,
    Likelihood_ratio_chisq = trv_model$null.deviance - trv_model$deviance,
    Degrees_of_freedom = trv_model$df.null - trv_model$df.residual,
    Likelihood_ratio_P = stats::pchisq(
      trv_model$null.deviance - trv_model$deviance,
      df = trv_model$df.null - trv_model$df.residual,
      lower.tail = FALSE
    ),
    McFadden_R2 = 1 - as.numeric(stats::logLik(trv_model)) /
      as.numeric(stats::logLik(
        stats::glm(
          trv_measured ~ 1,
          data = trv_model_df,
          family = stats::binomial()
        )
      ))
  )
  utils::write.csv(
    trv_model_fit,
    file.path(output_dir, "Table_TRVmax_Measurement_ModelFit.csv"),
    row.names = FALSE
  )

  trv_selection_smd <- selection_smd |>
    dplyr::filter(
      Comparison == "TRVmax measurement within E/e' + LAVi cohort"
    )
  ee_lavi_selection_smd <- selection_smd |>
    dplyr::filter(Comparison == "E/e' + LAVi measurement selection")
  stage1_summary <- tibble::tibble(
    Metric = c(
      "Phenotypic HCM entry N",
      "E/e' + LAVi cross-sectional N",
      "Excluded before cross-sectional cohort N",
      "TRVmax measured within cross-sectional cohort N",
      "TRVmax missing within cross-sectional cohort N",
      "TRVmax missing within cross-sectional cohort percent",
      "Clinical-outcomes subcohort N",
      "Longitudinal subcohort N",
      "TRVmax measurement model likelihood-ratio P",
      "TRVmax measurement model McFadden R2",
      "Largest absolute SMD: TRVmax measured vs missing",
      "Largest absolute SMD: E/e' + LAVi measured vs either missing"
    ),
    Value = c(
      nrow(baseline_df),
      nrow(analytic_df),
      nrow(baseline_df) - nrow(analytic_df),
      nrow(trv_measured_df),
      nrow(analytic_df) - nrow(trv_measured_df),
      100 * (nrow(analytic_df) - nrow(trv_measured_df)) / nrow(analytic_df),
      nrow(outcomes_df),
      dplyr::n_distinct(long_df$ID),
      trv_model_fit$Likelihood_ratio_P,
      trv_model_fit$McFadden_R2,
      max(trv_selection_smd$Absolute_SMD, na.rm = TRUE),
      max(ee_lavi_selection_smd$Absolute_SMD, na.rm = TRUE)
    )
  )
  utils::write.csv(
    stage1_summary,
    file.path(output_dir, "Stage1_Missingness_Selection_Summary.csv"),
    row.names = FALSE
  )

  list(
    variable_spec = variable_spec,
    missingness_overall = missingness_overall,
    missingness_by_event = missingness_by_event,
    cohort_flow = cohort_flow,
    selection_smd = selection_smd,
    trv_model = trv_model,
    trv_model_table = trv_model_table,
    trv_model_fit = trv_model_fit,
    stage1_summary = stage1_summary
  )
}
