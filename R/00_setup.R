# ── Shared setup ──────────────────────────────────────────────────────────────
# Sourced at the top of every stage script. Loads packages, reads
# 2_Config/config.yml, sources the helper files in this folder, and defines the
# handful of I/O helpers the stages share.
#
# Each stage runs in its own R process, so a stage can only use what an earlier
# stage wrote to disk.

if (!file.exists("2_Config/config.yml")) {
  stop("Run from the project root (the folder containing run_all.R).", call. = FALSE)
}

cfg <- yaml::read_yaml("2_Config/config.yml")

options(warn = 1)

if (nzchar(cfg$paths$r_library) && dir.exists(cfg$paths$r_library)) {
  .libPaths(c(normalizePath(cfg$paths$r_library), .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(ggtext)
  library(splines)
  library(lme4)
  library(nlme)
  library(mgcv)
  library(broom)
  library(survival)
  library(scales)
  library(viridis)
  library(ggrepel)
  library(openxlsx)
  library(gtsummary)
  library(flextable)
  library(digest)
})

options(
  hcm.rcs_boot_cache_dir  = file.path(cfg$paths$cache, "rcs_piecewise_boot"),
  hcm.rcs_bootstrap_reps  = cfg$rcs$bootstrap_reps,
  hcm.rcs_bootstrap_seed  = cfg$rcs$bootstrap_seed,
  hcm.gamm_stats_dir      = file.path(cfg$paths$results, "4_Longitudinal")
)

# Helper files: definitions, small utilities, and the model/figure helpers that
# more than one stage calls.
for (f in sort(list.files("R", pattern = "^[a-z_]+\\.R$", full.names = TRUE))) {
  source(f)
}

# Figures are written explicitly with ggsave()/cairo_pdf(). Keep a Cairo PNG
# device open as the current device so text measured while laying out legends
# uses the same font metrics, and printed plots go to a temp file rather than
# Rplots.pdf.
grDevices::png(tempfile(fileext = ".png"), width = 6.9, height = 4.5,
               units = "in", res = 150, type = "cairo")

# Jittered point layers use this seed (see position_jitter(seed = ...)).
JITTER_SEED <- cfg$figures$jitter_seed

# ── Stage I/O ─────────────────────────────────────────────────────────────────
# Patient-level objects are .rds under 5_Data/; aggregate outputs go to
# 6_Results/<section>/. Both helpers create directories as needed.
save_data <- function(x, name, where = c("interim", "processed")) {
  where <- match.arg(where)
  dir <- cfg$paths[[where]]
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file.path(dir, paste0(name, ".rds")))
  invisible(x)
}

load_data <- function(name, where = c("interim", "processed")) {
  where <- match.arg(where)
  path <- file.path(cfg$paths[[where]], paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing ", path, ". Run the earlier stages first.", call. = FALSE)
  }
  readRDS(path)
}

results_dir <- function(section) {
  dir <- file.path(cfg$paths$results, section)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

write_result <- function(x, section, file) {
  path <- file.path(results_dir(section), file)
  write.csv(x, path, row.names = FALSE)
  invisible(path)
}

stage_banner <- function(title) {
  message(strrep("=", 78))
  message(title, "  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "]")
  message(strrep("=", 78))
}
