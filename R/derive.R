# ── CPX derived fields ────────────────────────────────────────────────────────
# derive_cpx_fields() is the ONLY place the pipeline defines the CPX-derived
# variables used downstream:
#
#   age, BMI, BSA, LBM, LBMI,
#   VO2_FRIEND_PP, VO2_FRIEND2_PP, VO2_WASSERMAN_PP, HRmax_PP
#
# Every later stage (cohort filters, adjustment sets, matching, reference
# thresholds, outcomes) uses these columns and nothing else, so a change made
# here propagates to the entire analysis on the next `Rscript run_all.R`.
#
# DERIVED-FIELD CORRECTION (applied 2026-09-16 on the project owner's written
# authorization; owner remains responsible for the decision)
#   The `CPX` sheet of CPX_input.xlsx stores these quantities in calculated
#   columns whose formulas reference cells in OTHER rows, so a stored value
#   generally belongs to a different patient's test (verified 2026-09-15: the
#   stored age matched the referenced row's age in 5,160/5,160 rows and never
#   the row's own patient; stored age agreed with the independent age fields in
#   the Outcomes_2025 and medications sheets in 0.3% and 0.4% of rows).
#   Every quantity below is therefore recomputed from that row's own raw
#   inputs, using the workbook's own formulas:
#
#     age                = (cpx_test_date - DOB) / 365.25
#     weight (kg)        = Weight.(pounds) x 0.453592
#     height (cm)        = height.(in.) x 2.54            (0 in. -> missing)
#     BMI                = kg / (cm/100)^2
#     BSA                = sqrt(cm x kg / 3600)           (Mosteller)
#     LBM (NHANES)       male:   14.729 - 0.071*age + 0.21*cm  + 0.468*kg
#                        female: -14.292 - 0.046*age + 0.201*cm + 0.347*kg
#     LBMI               = LBM / (cm/100)^2
#     FRIEND1 predicted  = 79.9 - 0.39*age - 13.7*sex - 0.127*lb  (treadmill only)
#     FRIEND2 predicted  = 45.2 - 0.35*age - 10.9*(sex+1) - 0.15*lb
#                          + 0.68*in - 0.46*(mode+1)
#     predicted peak HR  = 209.3 - 0.72*age
#     Wasserman          five-step chain (cycle factor, predicted weight,
#                          over/under weight, ml/min, treadmill correction)
#     %-predicted        = observed / predicted x 100
#
#   Sex is coded 0 = male, 1 = female and mode 0 = treadmill, 1 = cycle, as in
#   the workbook; note FRIEND2 uses (sex+1) and (mode+1) while FRIEND1 and the
#   LBM/Wasserman branches use the raw coding.
#
#   ONE DELIBERATE DEVIATION: age uses 365.25 days per year (exact elapsed
#   years), where the workbook formula divides by 365. The difference is under
#   0.1 years. Owner decision, 2026-09-16.
#
#   Missing raw inputs stay missing; nothing is imputed or zero-filled. The
#   cached workbook columns are overwritten with the recomputed values so no
#   downstream code can read a misaligned value by accident.
#
#   Cross-file discrepancies (63 date-of-birth and ~55 sex mismatches between
#   the CPET and echo extracts, of ~4,900 comparable rows) are left as they
#   are: the CPET record is primary and the counts are reported, per owner
#   decision.
#
#   R/02_derive.R writes an aggregate agreement audit to
#   6_Results/1_Audit/CPX_Derived_Age_Agreement.csv after every run.

derive_cpx_fields <- function(cpx) {
  lb <- as_num(cpx[["Weight.(pounds)"]])
  inch <- as_num(cpx[["height.(in.)"]])
  sex <- as_num(cpx[["Sex:.M(0)/F(1)"]])
  mode <- as_num(cpx[["Treadmill.(0).or.cycle.(1)"]])
  vo2 <- as_num(cpx[["VO2.(ml/min/kg)"]])
  peak_hr <- as_num(cpx[["peak.HR.(bpm)"]])

  age_years <- as.numeric(as.Date(cpx$cpx_test_date) - as.Date(cpx$DOB)) / 365.25
  weight_kg <- lb * 0.453592
  height_cm <- ifelse(!is.na(inch) & inch == 0, NA_real_, inch * 2.54)

  lbm <- dplyr::case_when(
    is.na(sex) | is.na(age_years) | is.na(height_cm) | is.na(weight_kg) ~ NA_real_,
    sex == 0 ~ 14.729 - 0.071 * age_years + 0.21 * height_cm + 0.468 * weight_kg,
    sex == 1 ~ -14.292 - 0.046 * age_years + 0.201 * height_cm + 0.347 * weight_kg,
    TRUE ~ NA_real_
  )

  # FRIEND 2015 (treadmill only) and FRIEND 2018 reference equations
  friend1_pred <- ifelse(
    !is.na(mode) & mode == 0,
    79.9 - 0.39 * age_years - 13.7 * sex - 0.127 * lb,
    NA_real_
  )
  friend2_pred <- 45.2 - 0.35 * age_years - 10.9 * (sex + 1) -
    0.15 * lb + 0.68 * inch - 0.46 * (mode + 1)
  hr_pred <- 209.3 - 0.72 * age_years

  # Wasserman/Hansen predicted peak VO2
  cycle_factor <- dplyr::case_when(
    is.na(sex) | is.na(age_years) ~ NA_real_,
    sex == 0 ~ 50.72 - 0.372 * age_years,
    sex == 1 ~ 22.78 - 0.17 * age_years,
    TRUE ~ NA_real_
  )
  wasserman_pred_weight <- dplyr::case_when(
    is.na(sex) | is.na(height_cm) ~ NA_real_,
    sex == 0 ~ 0.79 * height_cm - 60.7,
    sex == 1 ~ 0.65 * height_cm - 42.8,
    TRUE ~ NA_real_
  )
  weight_class <- dplyr::case_when(
    is.na(weight_kg) | is.na(wasserman_pred_weight) ~ NA_real_,
    weight_kg > wasserman_pred_weight ~ 1,
    weight_kg == wasserman_pred_weight ~ 0,
    weight_kg < wasserman_pred_weight ~ -1,
    TRUE ~ NA_real_
  )
  wasserman_ml_min <- dplyr::case_when(
    is.na(sex) | is.na(weight_class) | is.na(cycle_factor) ~ NA_real_,
    sex == 0 & weight_class == -1 ~
      ((wasserman_pred_weight + weight_kg) / 2) * cycle_factor,
    sex == 0 & weight_class == 0 ~ weight_kg * cycle_factor,
    sex == 0 & weight_class == 1 ~
      wasserman_pred_weight * cycle_factor + 6 * (weight_kg - wasserman_pred_weight),
    sex == 1 & weight_class == -1 ~
      ((wasserman_pred_weight + weight_kg + 86) / 2) * cycle_factor,
    sex == 1 & weight_class == 0 ~ (weight_kg + 43) * cycle_factor,
    sex == 1 & weight_class == 1 ~
      (wasserman_pred_weight + 43) * cycle_factor + 6 * (weight_kg - wasserman_pred_weight),
    TRUE ~ NA_real_
  )
  wasserman_corrected <- dplyr::case_when(
    is.na(mode) | is.na(wasserman_ml_min) ~ NA_real_,
    mode == 0 ~ wasserman_ml_min * 1.11,
    mode == 1 ~ wasserman_ml_min,
    TRUE ~ NA_real_
  )
  wasserman_pred_vo2 <- wasserman_corrected / weight_kg

  pct_predicted <- function(observed, predicted) {
    ifelse(
      is.na(observed) | is.na(predicted) | !is.finite(predicted) | predicted <= 0,
      NA_real_,
      observed / predicted * 100
    )
  }

  cpx %>%
    dplyr::mutate(
      age = age_years,
      BSA = sqrt((height_cm * weight_kg) / 3600), # Mosteller
      BMI = weight_kg / ((height_cm / 100)^2),
      LBM = lbm,
      LBMI = LBM / (height_cm / 100)^2,
      VO2_FRIEND_PP = pct_predicted(vo2, friend1_pred),
      VO2_FRIEND2_PP = pct_predicted(vo2, friend2_pred),
      VO2_WASSERMAN_PP = pct_predicted(vo2, wasserman_pred_vo2),
      HRmax_PP = pct_predicted(peak_hr, hr_pred),
      # Overwrite the misaligned workbook columns so nothing downstream can
      # read a stale value.
      `Weight.(kg)` = weight_kg,
      `height.(cm)` = height_cm,
      LBM.NHANES.no.race = lbm,
      `FRIEND1.%.predicted.VO2` = VO2_FRIEND_PP,
      `FRIEND2.%predicted.VO2` = VO2_FRIEND2_PP,
      `Wasserman.%predicted.VO2.(2005)` = VO2_WASSERMAN_PP,
      `%predicted.HR.FRIEND` = HRmax_PP
    )
}

# Aggregate agreement between the derived age and the independent age fields
# in other sheets (Outcomes_2025 `age_y`, medications `age`). Returns counts
# only; no identifiers or row-level values.
audit_cpx_derived_age <- function(cpx_derived, outcomes_raw, medications_raw,
                                  tolerance_years = 0.1) {
  reference_age <- function(sheet, age_col, label) {
    sheet %>%
      dplyr::transmute(
        MRN = as.character(MRN),
        cpx_test_date = coerce_excel_date(cpx_test_date),
        reference_age = as_num(.data[[age_col]])
      ) %>%
      dplyr::filter(!is.na(cpx_test_date), !is.na(reference_age)) %>%
      dplyr::distinct(MRN, cpx_test_date, .keep_all = TRUE) %>%
      dplyr::mutate(reference = label)
  }

  refs <- dplyr::bind_rows(
    reference_age(outcomes_raw, "age_y", "Outcomes_2025 age_y"),
    reference_age(medications_raw, "age", "medications age")
  )

  cpx_derived %>%
    dplyr::transmute(
      MRN = as.character(MRN),
      cpx_test_date = as.Date(cpx_test_date),
      derived_age = as_num(age)
    ) %>%
    dplyr::inner_join(refs, by = c("MRN", "cpx_test_date"), relationship = "many-to-many") %>%
    dplyr::filter(!is.na(derived_age)) %>%
    dplyr::group_by(reference) %>%
    dplyr::summarise(
      n_compared = dplyr::n(),
      n_within_tolerance = sum(abs(derived_age - reference_age) <= tolerance_years),
      pct_within_tolerance = round(100 * n_within_tolerance / n_compared, 1),
      median_abs_difference_years = round(stats::median(abs(derived_age - reference_age)), 2),
      tolerance_years = tolerance_years,
      .groups = "drop"
    )
}
