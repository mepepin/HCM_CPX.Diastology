# ── Stage 01: import the raw workbooks ────────────────────────────────────────
# In  : 1_Input/{CPX_input.xlsx, mark_extract_*.xlsx, HCM_Genetics.xlsx}
# Out : 5_Data/1_Interim/{cpx_raw, comorbidities_raw, echo_raw, genetics_raw,
#                         outcomes_raw, medications_raw}.rds
#
# Sheets are read exactly as supplied; no value is modified here. The stable
# patient ID is assigned from MRN so that nothing downstream needs the MRN
# except the joins between sheets.

source("R/00_setup.R")
stage_banner("Stage 01: import raw source tables")

read_sheet <- function(path, sheet) {
  openxlsx::read.xlsx(path, sheet = sheet, detectDates = TRUE)
}

input_files <- c(
  cpx_workbook      = cfg$paths$cpx_workbook,
  echo_extract      = cfg$paths$echo_extract,
  genetics_workbook = cfg$paths$genetics_workbook
)
missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs)) {
  stop("Input file(s) not found: ", paste(missing_inputs, collapse = ", "), call. = FALSE)
}

# Record which files this run read, so any result can be tied back to the exact
# input workbooks that produced it.
write_result(
  tibble(
    input = names(input_files),
    path = unname(input_files),
    size_bytes = file.size(input_files),
    modified = format(file.mtime(input_files), "%Y-%m-%d %H:%M:%S"),
    sha256 = vapply(input_files, function(p) digest::digest(file = p, algo = "sha256"), character(1)),
    read_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  "1_Audit", "Input_File_Manifest.csv"
)

# The CPX sheet was once sorted after its formulas were written, which left
# formulas pointing at other rows. Excel caches formula results, so that damage
# is invisible once a sheet is read. Check the stored formulas and stop if a
# cross-row reference turns up anywhere not already recorded in the config.
assert_input_formula_integrity(
  unname(input_files),
  out_dir  = results_dir("1_Audit"),
  accepted = cfg$input_integrity$accepted_cross_row
)

# One row per test. cur_group_id() orders groups by MRN, so the ID mapping is
# deterministic across runs.
cpx_raw <- read_sheet(cfg$paths$cpx_workbook, "CPX") %>%
  group_by(MRN) %>%
  mutate(ID = cur_group_id()) %>%
  ungroup()

save_data(cpx_raw, "cpx_raw")
save_data(read_sheet(cfg$paths$cpx_workbook, "co-morbidities"), "comorbidities_raw")
save_data(read_sheet(cfg$paths$cpx_workbook, "Outcomes_2025"), "outcomes_raw")
save_data(read_sheet(cfg$paths$cpx_workbook, "medications"), "medications_raw")
save_data(read_sheet(cfg$paths$echo_extract, 1), "echo_raw")
save_data(read_sheet(cfg$paths$genetics_workbook, 1), "genetics_raw")

message(sprintf("CPX sheet: %d tests from %d patients",
                nrow(cpx_raw), n_distinct(cpx_raw$ID)))
message("Stage 01 complete.")
