# ── Stage 08: assemble manuscript deliverables ────────────────────────────────
# Inputs : 2_Config/manuscript_manifest.csv  (destination, output_file, source_file)
#          6_Results/<section>/...           (outputs of stages 03-07)
# Outputs: 6_Results/6_Manuscript/main/        main figures and tables
#          6_Results/6_Manuscript/supplement/  supplemental figures and tables
#          6_Results/6_Manuscript/MANIFEST.csv source, destination, and SHA-256
#
# Also runs the de-identification check over the whole results tree.
#
# The manifest is the single list of submission files. To rename or renumber a
# supplemental item, edit 2_Config/manuscript_manifest.csv and rerun this stage.
# The stage fails if any listed source file is missing.
source("R/00_setup.R")
stage_banner("Stage 08: assemble manuscript deliverables")

# ── Central Illustration ─────────────────────────────────────────────────────
# Built here, before the bundle is assembled, so it is regenerated on every run
# from the result files the other stages wrote. See R/central_illustration.R.
ci_dir <- results_dir("7_Central_Illustration")
ci_cohort <- load_data("baseline_df", "processed") %>% filter(analysis_base_eligible)
ci_years <- paste0(
  format(min(ci_cohort$cpx_test_date), "%Y"), "\u2013",
  format(max(ci_cohort$cpx_test_date), "%Y")
)
ci_values <- build_central_illustration(cfg$paths$results, ci_dir, ci_years)
write_central_illustration_legend(ci_values, file.path(ci_dir, "Central_Illustration_Legend.md"))
message("Central Illustration written to ", ci_dir)

manifest <- read.csv("2_Config/manuscript_manifest.csv", stringsAsFactors = FALSE)
stopifnot(
  all(c("destination", "output_file", "source_file") %in% names(manifest)),
  all(manifest$destination %in% c("main", "supplement")),
  !anyDuplicated(paste(manifest$destination, manifest$output_file))
)

manifest$source_path <- file.path(cfg$paths$results, manifest$source_file)
missing_sources <- manifest$source_path[!file.exists(manifest$source_path)]
if (length(missing_sources)) {
  stop("Manifest source file(s) not found:\n  ", paste(missing_sources, collapse = "\n  "), call. = FALSE)
}

manuscript_dir <- results_dir("6_Manuscript")
for (dest in unique(manifest$destination)) {
  dest_dir <- file.path(manuscript_dir, dest)
  unlink(dest_dir, recursive = TRUE) # rebuilt from the manifest on every run
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
}

manifest$output_path <- file.path(manuscript_dir, manifest$destination, manifest$output_file)
copied <- file.copy(manifest$source_path, manifest$output_path, overwrite = TRUE, copy.date = TRUE)
stopifnot(all(copied))

manifest$sha256 <- vapply(manifest$output_path, function(p) digest::digest(file = p, algo = "sha256"), character(1))
write.csv(
  manifest[, c("destination", "output_file", "source_file", "sha256")],
  file.path(manuscript_dir, "MANIFEST.csv"),
  row.names = FALSE
)

# ── Numbers ledger ───────────────────────────────────────────────────────────
# Every headline number the manuscript reports, read back out of the result
# files rather than transcribed, so the ledger cannot drift from the analysis.
# 8_Docs/NUMBERS_LEDGER_2026-09.md pairs these with the values currently in the
# draft; this file is the authoritative right-hand column.
read_result <- function(relative_path) {
  path <- file.path(cfg$paths$results, relative_path)
  if (!file.exists(path)) {
    stop("Numbers ledger source missing: ", relative_path, call. = FALSE)
  }
  utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

cohort_counts_tbl <- read_result("2_Cohort/Cohort_Counts.csv")
pick_count <- function(quantity) {
  hit <- cohort_counts_tbl$value[cohort_counts_tbl$quantity == quantity]
  stopifnot(length(hit) == 1)
  hit
}

parameter_counts_tbl <- read_result("2_Cohort/Cohort_ParameterSpecific_Counts.csv")
balance_tbl <- read_result("2_Cohort/Table_Comparison_Balance_OverlapWeights.csv")
rcs_tbl <- read_result("3_CrossSectional/Table2_Figure3_RCS_Summary.csv")
lmm_tbl <- read_result("4_Longitudinal/Table3_Longitudinal_LMM_Primary.csv")
cox_tbl <- read_result("5_Outcomes/TableS_Figure5_Cox_Summary.csv")

numbers_ledger <- dplyr::bind_rows(
  tibble::tibble(
    Section = "Cohort",
    Quantity = c(
      "Screened HCM patients",
      "Base analytic cohort",
      "Base analytic cohort with follow-up",
      "Composite events in the base cohort",
      "Longitudinal patients",
      "Longitudinal CPETs",
      "Longitudinal follow-up, median years"
    ),
    Value = as.character(c(
      pick_count("Screened HCM cohort (patients)"),
      pick_count("Base analytic cohort (patients)"),
      pick_count("Base analytic cohort with follow-up (patients)"),
      pick_count("Base analytic cohort composite events"),
      pick_count("Longitudinal patients"),
      pick_count("Longitudinal observations"),
      sprintf("%.2f", as.numeric(pick_count("Longitudinal follow-up, median years")))
    )),
    Source = "2_Cohort/Cohort_Counts.csv"
  ),
  parameter_counts_tbl %>%
    dplyr::transmute(
      Section = "Denominators by parameter",
      Quantity = sprintf("%s: cross-sectional / longitudinal patients / outcomes / events", sample),
      Value = sprintf(
        "%d / %d / %d / %d",
        n_cross_sectional, n_longitudinal_patients, n_outcomes_patients, n_outcomes_events
      ),
      Source = "2_Cohort/Cohort_ParameterSpecific_Counts.csv"
    ),
  tibble::tibble(
    Section = "Reference comparison",
    Quantity = c(
      "HCM patients / references contributing",
      "Effective sample sizes (overlap-weighted)",
      "Largest |SMD| unweighted then weighted"
    ),
    Value = c(
      sprintf("%d / %d", balance_tbl$n_hcm[1], balance_tbl$n_reference[1]),
      sprintf("%.0f / %.0f", balance_tbl$effective_n_hcm[1], balance_tbl$effective_n_reference[1]),
      sprintf(
        "%.2f then %.4f",
        max(abs(balance_tbl[["SMD (unweighted)"]])),
        max(abs(balance_tbl[["SMD (overlap-weighted)"]]))
      )
    ),
    Source = "2_Cohort/Table_Comparison_Balance_OverlapWeights.csv"
  ),
  rcs_tbl %>%
    dplyr::transmute(
      Section = "Cross-sectional splines",
      Quantity = sprintf("%s ~ %s: N, overall P, nonlinearity q", Outcome, Parameter),
      Value = sprintf("%d, %s, %s", N, `Overall p`, `Nonlinearity q`),
      Source = "3_CrossSectional/Table2_Figure3_RCS_Summary.csv"
    ),
  lmm_tbl %>%
    dplyr::transmute(
      Section = "Trajectories",
      Quantity = sprintf("%s ~ %s: patients, interaction beta, P", Outcome, Parameter),
      Value = sprintf(
        "%d, %.3f, %s",
        Patients, `Interaction beta per year per 1 SD`, `P value`
      ),
      Source = "4_Longitudinal/Table3_Longitudinal_LMM_Primary.csv"
    ),
  cox_tbl %>%
    dplyr::transmute(
      Section = "Outcomes",
      Quantity = sprintf("%s: N / events, HR (95%% CI), P, BH q", Parameter),
      Value = sprintf("%s, %s, %s, %s", Sample, `HR (95% CI)`, `P value`, `BH q value`),
      Source = "5_Outcomes/TableS_Figure5_Cox_Summary.csv"
    )
)

stopifnot(
  nrow(numbers_ledger) > 25,
  !any(is.na(numbers_ledger$Value)),
  all(nzchar(numbers_ledger$Value)),
  !any(grepl("NA", numbers_ledger$Value, fixed = TRUE))
)

write.csv(
  numbers_ledger,
  file.path(manuscript_dir, "Numbers_Ledger.csv"),
  row.names = FALSE
)
message(sprintf("Numbers ledger: %d headline values.", nrow(numbers_ledger)))

# Nothing patient-level may reach the exported results (see R/validate.R).
assert_results_deidentified(cfg$paths$results)


message("Stage 08 complete.")
