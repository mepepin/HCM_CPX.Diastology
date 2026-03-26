library(openxlsx)
library(dplyr)
library(tidyr)

project_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
if (basename(project_root) == "_Scripts") {
  project_root <- dirname(project_root)
}

input_file <- file.path(project_root, "1_Input", "CPX_input.xlsx")
output_dir <- file.path(project_root, "2_Output", "Cycle")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "Cycle_Linked_Tests.xlsx")

CPX_raw <- read.xlsx(input_file, sheet = "CPX", detectDates = TRUE)
Outcomes_raw <- read.xlsx(input_file, sheet = "Outcomes_2025", detectDates = TRUE)
Comorbidities_raw <- read.xlsx(input_file, sheet = "co-morbidities", detectDates = TRUE)

cycling_pattern <- "cycl|bike|erg|supine|upright|ice"
noncycling_modes <- c("walking", "running", "vo2 study")

CPX_cycle <- CPX_raw %>%
  mutate(
    mode_clean = tolower(trimws(dplyr::coalesce(Mode, ""))),
    protocol_clean = tolower(trimws(dplyr::coalesce(Protocol, ""))),
    cycle_flag = suppressWarnings(as.numeric(`Treadmill.(0).or.cycle.(1)`)),
    is_cycling_test =
      ((cycle_flag == 1) & !(mode_clean %in% noncycling_modes)) |
      grepl(cycling_pattern, mode_clean) |
      grepl(cycling_pattern, protocol_clean),
    Exercise_mode = if_else(
      nzchar(mode_clean),
      Mode,
      if_else(cycle_flag == 1, "Cycling", NA_character_)
    ),
    Cycling_form = case_when(
      grepl("supine", protocol_clean) ~ "Supine bike",
      grepl("upright", protocol_clean) ~ "Upright bike",
      grepl("(^|[^a-z])ice([^a-z]|$)", protocol_clean) ~ "ICE/research bike",
      grepl("manual", protocol_clean) ~ "Manual bike",
      is_cycling_test ~ "Other/unspecified cycling",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(is_cycling_test) %>%
  mutate(MRN = as.character(MRN)) %>%
  group_by(MRN) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup() %>%
  mutate(
    cpx_test_date = as.Date(cpx_test_date),
    BSA = sqrt((`height.(cm)` * `Weight.(kg)`) / 3600),
    BMI = `Weight.(kg)` / ((`height.(cm)` / 100) ^ 2),
    LBM = LBM.NHANES.no.race,
    LBMI = LBM / ((`height.(cm)` / 100) ^ 2)
  ) %>%
  select(
    ID, MRN, cpx_test_date, age,
    Sex = `Sex:.M(0)/F(1)`,
    Exercise_mode, Cycling_form, Protocol, CPX.sequential,
    BMI, BSA, LBM, LBMI,
    pk.RER,
    VO2_kg = `VO2.(ml/min/kg)`,
    VO2_FRIEND_PP = `FRIEND1.%.predicted.VO2`,
    VO2_FRIEND2_PP = `FRIEND2.%predicted.VO2`,
    VO2_WASSERMAN_PP = `Wasserman.%predicted.VO2.(2005)`,
    HRmax_PP = `%predicted.HR.FRIEND`,
    HRR = `HR.recovery.(1min)`,
    VeVco2_slope = `ve/vco2.slope`,
    OUES, Bike_Watts = `Bike.Watts`, Arrhythmia, Notes
  ) %>%
  filter(
    is.na(VO2_FRIEND2_PP) | dplyr::between(VO2_FRIEND2_PP, 0, 200),
    is.na(VO2_FRIEND_PP) | dplyr::between(VO2_FRIEND_PP, 0, 200),
    is.na(VO2_WASSERMAN_PP) | dplyr::between(VO2_WASSERMAN_PP, 0, 200)
  )

Outcomes_clean <- Outcomes_raw %>%
  transmute(
    MRN = as.character(MRN),
    cpx_test_date = as.Date(cpx_test_date, origin = "1899-12-30"),
    last_enc_date = as.Date(last_enc_date, origin = "1899-12-30"),
    death = as.numeric(death),
    death_yrs = as.numeric(death_yrs),
    post_acute_heart_failure = as.numeric(post_acute_heart_failure),
    post_acute_heart_failure_yrs = as.numeric(post_acute_heart_failure_yrs),
    post_cardiac_arrest = as.numeric(post_cardiac_arrest),
    post_cardiac_arrest_yrs = as.numeric(post_cardiac_arrest_yrs),
    post_heart_transplant = as.numeric(post_heart_transplant),
    post_heart_transplant_yrs = as.numeric(post_heart_transplant_yrs),
    post_ventricular_tachycardia = as.numeric(post_ventricular_tachycardia),
    post_ventricular_tachycardia_yrs = as.numeric(post_ventricular_tachycardia_yrs),
    post_ventricular_fib_flut_tachy = as.numeric(post_ventricular_fib_flut_tachy),
    post_ventricular_fib_flut_tachy_yrs = as.numeric(post_ventricular_fib_flut_tachy_yrs),
    post_defibrillator = as.numeric(post_defibrillator),
    post_defibrillator_yrs = as.numeric(post_defibrillator_yrs),
    post_afib_flut = as.numeric(post_afib_flut),
    post_afib_flut_yrs = as.numeric(post_afib_flut_yrs),
    post_septal_reduction_surgery = as.numeric(post_septal_reduction_surgery),
    post_septal_reduction_surgery_yrs = as.numeric(post_septal_reduction_surgery_yrs)
  ) %>%
  group_by(MRN, cpx_test_date) %>%
  summarise(
    last_enc_date = if (all(is.na(last_enc_date))) as.Date(NA) else max(last_enc_date, na.rm = TRUE),
    across(
      c(
        death, post_acute_heart_failure, post_cardiac_arrest, post_heart_transplant,
        post_ventricular_tachycardia, post_ventricular_fib_flut_tachy,
        post_defibrillator, post_afib_flut, post_septal_reduction_surgery
      ),
      ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)
    ),
    across(
      c(
        death_yrs, post_acute_heart_failure_yrs, post_cardiac_arrest_yrs,
        post_heart_transplant_yrs, post_ventricular_tachycardia_yrs,
        post_ventricular_fib_flut_tachy_yrs, post_defibrillator_yrs,
        post_afib_flut_yrs, post_septal_reduction_surgery_yrs
      ),
      ~ if (all(is.na(.x))) NA_real_ else min(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

Comorbidities_clean <- Comorbidities_raw %>%
  transmute(
    MRN = as.character(MRN),
    cpx_test_date = as.Date(cpx_test_date),
    hf_dx = as.numeric(hf_pre_test),
    phtn_dx = as.numeric(
      pmax(
        dplyr::coalesce(as.numeric(pah_pre_test), 0),
        dplyr::coalesce(as.numeric(pahtn_pre_test), 0),
        dplyr::coalesce(as.numeric(pulm_htn_pre_test), 0),
        na.rm = TRUE
      )
    )
  ) %>%
  group_by(MRN, cpx_test_date) %>%
  summarise(
    hf_dx = if (all(is.na(hf_dx))) NA_real_ else max(hf_dx, na.rm = TRUE),
    phtn_dx = if (all(is.na(phtn_dx))) NA_real_ else max(phtn_dx, na.rm = TRUE),
    .groups = "drop"
  )

linked_cycle_tests <- CPX_cycle %>%
  left_join(Comorbidities_clean, by = c("MRN", "cpx_test_date")) %>%
  left_join(Outcomes_clean, by = c("MRN", "cpx_test_date")) %>%
  mutate(
    Disease_group = case_when(
      phtn_dx == 1 ~ "PHTN",
      hf_dx == 1 ~ "HF",
      TRUE ~ "No HF/PHTN"
    ),
    Sex = case_when(
      Sex == 0 ~ "Male",
      Sex == 1 ~ "Female",
      TRUE ~ NA_character_
    ),
    ventricular_arrhythmia_composite = pmax(
      dplyr::coalesce(post_ventricular_tachycardia, 0),
      dplyr::coalesce(post_ventricular_fib_flut_tachy, 0),
      na.rm = TRUE
    ),
    ventricular_arrhythmia_yrs = pmin(
      if_else(post_ventricular_tachycardia == 1, post_ventricular_tachycardia_yrs, Inf),
      if_else(post_ventricular_fib_flut_tachy == 1, post_ventricular_fib_flut_tachy_yrs, Inf),
      na.rm = TRUE
    ),
    primary_composite = pmax(
      dplyr::coalesce(death, 0),
      dplyr::coalesce(post_cardiac_arrest, 0),
      dplyr::coalesce(post_heart_transplant, 0),
      na.rm = TRUE
    ),
    primary_composite_yrs = pmin(
      if_else(death == 1, death_yrs, Inf),
      if_else(post_cardiac_arrest == 1, post_cardiac_arrest_yrs, Inf),
      if_else(post_heart_transplant == 1, post_heart_transplant_yrs, Inf),
      na.rm = TRUE
    ),
    ventricular_arrhythmia_yrs = ifelse(is.infinite(ventricular_arrhythmia_yrs), NA_real_, ventricular_arrhythmia_yrs),
    primary_composite_yrs = ifelse(is.infinite(primary_composite_yrs), NA_real_, primary_composite_yrs),
    follow_up_yrs = as.numeric(last_enc_date - cpx_test_date) / 365.25
  ) %>%
  filter(!is.na(death)) %>%
  arrange(MRN, cpx_test_date, CPX.sequential) %>%
  mutate(
    across(
      c(
        death, post_acute_heart_failure, post_cardiac_arrest, post_heart_transplant,
        post_ventricular_tachycardia, post_ventricular_fib_flut_tachy, post_defibrillator,
        post_afib_flut, post_septal_reduction_surgery, primary_composite,
        ventricular_arrhythmia_composite, hf_dx, phtn_dx
      ),
      ~ case_when(.x == 1 ~ "Yes", .x == 0 ~ "No", TRUE ~ NA_character_)
    )
  ) %>%
  rename(
    Patient_ID = ID,
    Test_Date = cpx_test_date,
    Test_Sequence = CPX.sequential,
    Peak_RER = pk.RER,
    Peak_VO2_mL_kg_min = VO2_kg,
    Peak_VO2_FRIEND1_pct = VO2_FRIEND_PP,
    Peak_VO2_FRIEND2_pct = VO2_FRIEND2_PP,
    Peak_VO2_Wasserman2005_pct = VO2_WASSERMAN_PP,
    Peak_HR_pct_predicted = HRmax_PP,
    HR_Recovery_1min = HRR,
    VE_VCO2_Slope = VeVco2_slope,
    Last_Encounter_Date = last_enc_date,
    Death = death,
    Death_Years = death_yrs,
    Acute_Heart_Failure = post_acute_heart_failure,
    Acute_Heart_Failure_Years = post_acute_heart_failure_yrs,
    Cardiac_Arrest = post_cardiac_arrest,
    Cardiac_Arrest_Years = post_cardiac_arrest_yrs,
    Heart_Transplant = post_heart_transplant,
    Heart_Transplant_Years = post_heart_transplant_yrs,
    Ventricular_Tachycardia = post_ventricular_tachycardia,
    Ventricular_Tachycardia_Years = post_ventricular_tachycardia_yrs,
    Ventricular_Fib_Flutter_Tachy = post_ventricular_fib_flut_tachy,
    Ventricular_Fib_Flutter_Tachy_Years = post_ventricular_fib_flut_tachy_yrs,
    Defibrillator = post_defibrillator,
    Defibrillator_Years = post_defibrillator_yrs,
    AFib_Flutter = post_afib_flut,
    AFib_Flutter_Years = post_afib_flut_yrs,
    Septal_Reduction_Surgery = post_septal_reduction_surgery,
    Septal_Reduction_Surgery_Years = post_septal_reduction_surgery_yrs,
    Heart_Failure_PreTest = hf_dx,
    Pulmonary_Hypertension_PreTest = phtn_dx,
    Ventricular_Arrhythmia_Composite = ventricular_arrhythmia_composite,
    Ventricular_Arrhythmia_Composite_Years = ventricular_arrhythmia_yrs,
    Primary_Composite = primary_composite,
    Primary_Composite_Years = primary_composite_yrs,
    FollowUp_Years = follow_up_yrs
  ) %>%
  select(
    Patient_ID, MRN, Test_Date, Test_Sequence, Disease_group,
    age, Sex, Exercise_mode, Cycling_form, Protocol,
    BMI, BSA, LBM, LBMI, Peak_RER, Peak_VO2_mL_kg_min,
    Peak_VO2_FRIEND1_pct, Peak_VO2_FRIEND2_pct, Peak_VO2_Wasserman2005_pct,
    Peak_HR_pct_predicted, HR_Recovery_1min, VE_VCO2_Slope, OUES, Bike_Watts,
    Arrhythmia, Notes, Heart_Failure_PreTest, Pulmonary_Hypertension_PreTest,
    Last_Encounter_Date, FollowUp_Years,
    Primary_Composite, Primary_Composite_Years,
    Ventricular_Arrhythmia_Composite, Ventricular_Arrhythmia_Composite_Years,
    Death, Death_Years,
    Acute_Heart_Failure, Acute_Heart_Failure_Years,
    Cardiac_Arrest, Cardiac_Arrest_Years,
    Heart_Transplant, Heart_Transplant_Years,
    Ventricular_Tachycardia, Ventricular_Tachycardia_Years,
    Ventricular_Fib_Flutter_Tachy, Ventricular_Fib_Flutter_Tachy_Years,
    Defibrillator, Defibrillator_Years,
    AFib_Flutter, AFib_Flutter_Years,
    Septal_Reduction_Surgery, Septal_Reduction_Surgery_Years
  )

summary_overall <- tibble(
  Metric = c(
    "Export date",
    "Linked cycling CPX tests",
    "Unique patients",
    "Median follow-up, years",
    "Primary composite events",
    "Deaths"
  ),
  Value = c(
    format(Sys.Date()),
    nrow(linked_cycle_tests),
    n_distinct(linked_cycle_tests$MRN),
    format(round(median(linked_cycle_tests$FollowUp_Years, na.rm = TRUE), 2), nsmall = 2),
    sum(linked_cycle_tests$Primary_Composite == "Yes", na.rm = TRUE),
    sum(linked_cycle_tests$Death == "Yes", na.rm = TRUE)
  )
)

summary_by_group <- linked_cycle_tests %>%
  count(Disease_group, name = "Tests") %>%
  left_join(
    linked_cycle_tests %>%
      group_by(Disease_group) %>%
      summarise(
        Patients = n_distinct(MRN),
        Primary_Composite_Events = sum(Primary_Composite == "Yes", na.rm = TRUE),
        Deaths = sum(Death == "Yes", na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Disease_group"
  )

data_dictionary <- tibble(
  Column = names(linked_cycle_tests),
  Description = c(
    "Stable patient identifier generated from MRN within the cycle cohort",
    "Medical record number",
    "CPX test date",
    "Sequential CPX number from source sheet",
    "Phenotype assignment used in the cycle report",
    "Age at CPX, years",
    "Biologic sex from source sheet",
    "Exercise mode field from source sheet",
    "Derived cycling form",
    "Protocol label from source sheet",
    "Body mass index, kg/m^2",
    "Body surface area, m^2",
    "Lean body mass, kg",
    "Lean body mass index, kg/m^2",
    "Peak respiratory exchange ratio",
    "Peak VO2, mL/kg/min",
    "Peak VO2 percent predicted, FRIEND 1",
    "Peak VO2 percent predicted, FRIEND 2",
    "Peak VO2 percent predicted, Wasserman 2005",
    "Peak heart rate percent predicted",
    "Heart rate recovery at 1 minute",
    "VE/VCO2 slope",
    "Oxygen uptake efficiency slope",
    "Peak bike workload, watts",
    "Arrhythmia field from source sheet",
    "Notes from source sheet",
    "Pre-test heart failure flag from comorbidity sheet",
    "Pre-test pulmonary hypertension flag from comorbidity sheet",
    "Last encounter date from outcomes sheet",
    "Observed follow-up from CPX to last encounter, years",
    "Composite of death, cardiac arrest, or transplant",
    "Years from CPX to primary composite event",
    "Composite of ventricular tachycardia or ventricular fibrillation/flutter/tachycardia",
    "Years from CPX to ventricular arrhythmia composite event",
    "Death after CPX",
    "Years from CPX to death",
    "Acute heart failure hospitalization after CPX",
    "Years from CPX to acute heart failure",
    "Cardiac arrest after CPX",
    "Years from CPX to cardiac arrest",
    "Heart transplant after CPX",
    "Years from CPX to heart transplant",
    "Ventricular tachycardia after CPX",
    "Years from CPX to ventricular tachycardia",
    "Ventricular fibrillation/flutter/tachycardia after CPX",
    "Years from CPX to ventricular fibrillation/flutter/tachycardia",
    "Defibrillator after CPX",
    "Years from CPX to defibrillator",
    "Atrial fibrillation/flutter after CPX",
    "Years from CPX to atrial fibrillation/flutter",
    "Septal reduction surgery after CPX",
    "Years from CPX to septal reduction surgery"
  )
)

wb <- createWorkbook(creator = "OpenAI Codex")

addWorksheet(wb, "Summary", gridLines = FALSE)
addWorksheet(wb, "Linked Tests", gridLines = FALSE)
addWorksheet(wb, "Data Dictionary", gridLines = FALSE)

title_style <- createStyle(
  fontSize = 14, textDecoration = "bold", fgFill = "#DCE6F1",
  halign = "left", border = "Bottom", borderColour = "#4F81BD"
)
section_style <- createStyle(
  fontSize = 11, textDecoration = "bold", fgFill = "#EAF2F8",
  border = "Bottom", borderColour = "#A9CCE3"
)
header_style <- createStyle(
  textDecoration = "bold", fgFill = "#1F4E78", fontColour = "#FFFFFF",
  halign = "center", valign = "center", border = "TopBottomLeftRight",
  borderColour = "#FFFFFF"
)
body_style <- createStyle(valign = "top")
date_style <- createStyle(numFmt = "yyyy-mm-dd")
decimal_style <- createStyle(numFmt = "0.00")
integer_style <- createStyle(numFmt = "0")
wrap_style <- createStyle(valign = "top", wrapText = TRUE)

writeData(wb, "Summary", x = "Ergocycle CPX linked outcomes export", startCol = 1, startRow = 1)
addStyle(wb, "Summary", title_style, rows = 1, cols = 1, gridExpand = TRUE)

writeData(wb, "Summary", x = "Overall", startCol = 1, startRow = 3)
addStyle(wb, "Summary", section_style, rows = 3, cols = 1, gridExpand = TRUE)
writeData(wb, "Summary", x = summary_overall, startCol = 1, startRow = 4, headerStyle = header_style)

writeData(wb, "Summary", x = "By phenotype", startCol = 1, startRow = 12)
addStyle(wb, "Summary", section_style, rows = 12, cols = 1, gridExpand = TRUE)
writeData(wb, "Summary", x = summary_by_group, startCol = 1, startRow = 13, headerStyle = header_style)
setColWidths(wb, "Summary", cols = 1:5, widths = c(28, 16, 14, 24, 12))
freezePane(wb, "Summary", firstActiveRow = 4)

writeDataTable(
  wb, "Linked Tests", x = linked_cycle_tests,
  tableStyle = "TableStyleMedium2", withFilter = TRUE
)
addStyle(
  wb, "Linked Tests", body_style,
  rows = 2:(nrow(linked_cycle_tests) + 1),
  cols = 1:ncol(linked_cycle_tests), gridExpand = TRUE, stack = TRUE
)

date_cols <- which(names(linked_cycle_tests) %in% c("Test_Date", "Last_Encounter_Date"))
decimal_cols <- which(names(linked_cycle_tests) %in% c(
  "BMI", "BSA", "LBM", "LBMI", "Peak_RER", "Peak_VO2_mL_kg_min",
  "Peak_VO2_FRIEND1_pct", "Peak_VO2_FRIEND2_pct", "Peak_VO2_Wasserman2005_pct",
  "Peak_HR_pct_predicted", "HR_Recovery_1min", "VE_VCO2_Slope", "OUES",
  "FollowUp_Years", "Primary_Composite_Years", "Ventricular_Arrhythmia_Composite_Years",
  "Death_Years", "Acute_Heart_Failure_Years", "Cardiac_Arrest_Years",
  "Heart_Transplant_Years", "Ventricular_Tachycardia_Years",
  "Ventricular_Fib_Flutter_Tachy_Years", "Defibrillator_Years",
  "AFib_Flutter_Years", "Septal_Reduction_Surgery_Years"
))
integer_cols <- which(names(linked_cycle_tests) %in% c("Patient_ID", "Test_Sequence", "age", "Bike_Watts"))
wrap_cols <- which(names(linked_cycle_tests) %in% c("Protocol", "Arrhythmia", "Notes"))

if (length(date_cols) > 0) {
  addStyle(
    wb, "Linked Tests", date_style,
    rows = 2:(nrow(linked_cycle_tests) + 1), cols = date_cols,
    gridExpand = TRUE, stack = TRUE
  )
}
if (length(decimal_cols) > 0) {
  addStyle(
    wb, "Linked Tests", decimal_style,
    rows = 2:(nrow(linked_cycle_tests) + 1), cols = decimal_cols,
    gridExpand = TRUE, stack = TRUE
  )
}
if (length(integer_cols) > 0) {
  addStyle(
    wb, "Linked Tests", integer_style,
    rows = 2:(nrow(linked_cycle_tests) + 1), cols = integer_cols,
    gridExpand = TRUE, stack = TRUE
  )
}
if (length(wrap_cols) > 0) {
  addStyle(
    wb, "Linked Tests", wrap_style,
    rows = 2:(nrow(linked_cycle_tests) + 1), cols = wrap_cols,
    gridExpand = TRUE, stack = TRUE
  )
}

setColWidths(wb, "Linked Tests", cols = 1:ncol(linked_cycle_tests), widths = "auto")
setColWidths(wb, "Linked Tests", cols = wrap_cols, widths = 30)
freezePane(wb, "Linked Tests", firstActiveRow = 2, firstActiveCol = 6)

writeDataTable(
  wb, "Data Dictionary", x = data_dictionary,
  tableStyle = "TableStyleLight9", withFilter = TRUE
)
setColWidths(wb, "Data Dictionary", cols = 1:2, widths = c(34, 88))
addStyle(
  wb, "Data Dictionary", wrap_style,
  rows = 2:(nrow(data_dictionary) + 1), cols = 1:2,
  gridExpand = TRUE, stack = TRUE
)
freezePane(wb, "Data Dictionary", firstActiveRow = 2)

saveWorkbook(wb, output_file, overwrite = TRUE)

cat(sprintf("Wrote %s with %d linked cycling tests.\n", output_file, nrow(linked_cycle_tests)))
