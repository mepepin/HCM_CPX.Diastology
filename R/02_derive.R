# ── Stage 02: CPX derived fields ──────────────────────────────────────────────
# In  : 5_Data/1_Interim/cpx_raw.rds
# Out : 5_Data/1_Interim/cpx_derived.rds
#
# Age, BMI, BSA, LBM, LBMI and the %-predicted VO2/HR fields are recomputed from
# same-row raw inputs by derive_cpx_fields() in R/derive.R. They are recomputed
# rather than read from the workbook because the cached values in those columns
# belong to other rows (see the note in 2_Config/config.yml).

source("R/00_setup.R")
stage_banner("Stage 02: CPX derived fields")

cpx_raw <- load_data("cpx_raw")

cpx_derived <- derive_cpx_fields(cpx_raw)
stopifnot(nrow(cpx_derived) == nrow(cpx_raw))

save_data(cpx_derived, "cpx_derived")

message("Stage 02 complete.")
