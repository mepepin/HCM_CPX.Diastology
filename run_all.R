#!/usr/bin/env Rscript
# ── Run the HCM diastology analysis pipeline ──────────────────────────────────
# Usage (from the project root):
#   Rscript run_all.R            run every stage in order
#   Rscript run_all.R 03 05      run only the listed stages
#
# Stages live in R/ and run in order:
#   01_import      read the three source workbooks
#   02_derive      CPX derived fields (age, BMI, %-predicted VO2)
#   03_cohort      cohort assembly and the counts the manuscript reports
#   04_cohort_table1   reference comparison, Figure 1, Table 1, Figure S1
#   05_cross_sectional Figures 2, 3, S2 and the cross-sectional tables
#   06_longitudinal    Figure 4, Figure S3 and the mixed models
#   07_outcomes        Figure 5, Figure S4 and the Cox models
#   08_manuscript      Central Illustration, submission bundle, numbers ledger
#
# Each stage runs in a fresh R process, so a stage can only use what earlier
# stages wrote to disk. Output goes to 7_Logs/<stage>.log and the run stops at
# the first failing stage.

cfg <- yaml::read_yaml("2_Config/config.yml")
stage_dir <- "R"
log_dir <- cfg$paths$logs

stage_files <- sort(list.files(stage_dir, pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
stage_files <- stage_files[basename(stage_files) != "00_setup.R"]
if (!length(stage_files)) stop("No stage scripts found. Run from the project root.")

selected <- commandArgs(trailingOnly = TRUE)
if (length(selected)) {
  stage_files <- stage_files[substr(basename(stage_files), 1, 2) %in% selected]
  if (!length(stage_files)) stop("None of the requested stages exist: ", paste(selected, collapse = ", "))
}

dir.create(log_dir, showWarnings = FALSE)
rscript <- file.path(R.home("bin"), "Rscript")
run_start <- Sys.time()

for (stage in stage_files) {
  log_file <- file.path(log_dir, sub("\\.R$", ".log", basename(stage)))
  message(format(Sys.time(), "%H:%M:%S"), "  running ", basename(stage))
  t0 <- Sys.time()
  status <- system2(rscript, c("--no-save", "--no-restore", shQuote(stage)),
                    stdout = log_file, stderr = log_file)
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (status != 0) {
    message("\n---- last lines of ", log_file, " ----")
    message(paste(tail(readLines(log_file, warn = FALSE), 25), collapse = "\n"))
    stop(sprintf("Stage %s failed (exit status %d).", basename(stage), status), call. = FALSE)
  }
  message(sprintf("          finished in %.1f min", minutes))
}

writeLines(capture.output(sessioninfo::session_info()), file.path(log_dir, "session_info.txt"))
message(sprintf("Pipeline complete in %.1f min.",
                as.numeric(difftime(Sys.time(), run_start, units = "mins"))))
