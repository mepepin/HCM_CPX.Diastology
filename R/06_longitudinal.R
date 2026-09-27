# ── Stage 06: longitudinal analyses ───────────────────────────────────────────
# In  : 5_Data/2_Processed/{baseline_df, long_df, cohort_parameters}.rds
# Out : 6_Results/4_Longitudinal/
#         Figure4_Longitudinal_GAMM.pdf              (main Figure 4)
#         FigureS_LAVi_Binary_Trajectories.pdf       (Figure S3)
#         Table3_Longitudinal_LMM_Primary.csv/.docx  (Table S2, primary inference)
#         Table3_Figure4_GAMM_Summary.csv/.docx      (GAMM refit)
#         Table_Longitudinal_FollowupHorizon.csv     (per-model follow-up)
#         Table_Consistency_CommonCohort_LMM.csv, Table_Sensitivity_AdultsOnly_LMM.csv,
#         Table_Sensitivity_MissingIndicator_LMM.csv,
#         Table_Sensitivity_MultipleImputation_LMM.csv  (Methods sensitivities)
#       5_Data/2_Processed/z_source.rds              (LAVi Z-scores used by stage 07)
#
# Inference on trajectory modification comes from the linear mixed models; the
# generalized additive mixed models are a prespecified exploratory check on
# whether those trajectories are linear, and are not corrected for multiplicity.

source("R/00_setup.R")
stage_banner("Stage 06: longitudinal analyses")

OUT_DIR <- results_dir("4_Longitudinal")
AUDIT_DIR <- results_dir("1_Audit")

# ── Inputs from earlier stages ────────────────────────────────────────────────
baseline_df <- load_data("baseline_df", "processed")
long_df <- load_data("long_df", "processed")
cohort_parameters <- load_data("cohort_parameters", "processed")
analytic_cohort_mode <- cohort_parameters$analytic_cohort_mode

# ══ GAMLSS Z-Score Derivation ═════════════════════════════════════════════════
# (legacy chunk: figure4_gamlss_zscore)
# ── age/sex-referenced LAVi Z-scores (GAMLSS) ─────────────────────────────────
# model both mean (mu) and SD (sigma) as smooth fns of age + sex via pb() splines,
# then z = (obs - mu_hat)/sigma_hat. falls back to plain scale() if the fit fails.
# these z-scores feed the Cox models downstream
derive_lavi_gamlss_results <- function(df) {
  z_dat <- df %>%
    transmute(
      ID,
      age = suppressWarnings(as.numeric(age)),
      Sex = factor(Sex),
      la_vol_index = suppressWarnings(as.numeric(la_vol_index))
    ) %>%
    filter(!is.na(ID), !is.na(age), !is.na(Sex), !is.na(la_vol_index))

  if (nrow(z_dat) == 0) {
    return(list(
      la_vol_index = list(
        model = NULL,
        data = tibble(ID = numeric(), seamlss_z = numeric())
      )
    ))
  }

  fit_obj <- tryCatch(
    gamlss(
      la_vol_index ~ pb(age) + Sex,
      sigma.formula = ~ pb(age) + Sex,
      family = NO,
      data = z_dat,
      trace = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(fit_obj)) {
    mu_hat <- fitted(fit_obj, what = "mu")
    sigma_hat <- pmax(fitted(fit_obj, what = "sigma"), 1e-6)
    z_dat$seamlss_z <- (z_dat$la_vol_index - mu_hat) / sigma_hat
  } else {
    z_dat$seamlss_z <- as.numeric(scale(z_dat$la_vol_index))
  }

  list(
    la_vol_index = list(model = fit_obj, data = z_dat %>% select(ID, seamlss_z))
  )
}

# Derive z-scores from HCM baseline data
gamlss_combined_baseline <- baseline_df
gamlss_results <- derive_lavi_gamlss_results(gamlss_combined_baseline)

# ══ Continuous GAMM Trajectory Panels ═════════════════════════════════════════
# (legacy chunk: figure4_gamm)
# ── Fig 4: GAMM trajectories of exercise capacity by baseline diastolic burden ─
# mixed GAM (smooth time x stratum + random patient intercept) on the repeated-
# test cohort. stratify by baseline LVDD burden (#abnormal indices) and by each
# binary index. cap follow-up at 10y to avoid sparse-tail extrapolation
# z_source = per-patient LAVi z for the Cox models
z_source <- gamlss_results[["la_vol_index"]]$data %>%
  select(ID, dd_zscore = seamlss_z) %>%
  filter(!is.na(dd_zscore)) %>%
  group_by(ID) %>%
  summarise(dd_zscore = first(dd_zscore), .groups = "drop")

lvdd_burden_lookup <- baseline_df %>%
  transmute(
    ID,
    lvdd_burden = case_when(
      hcm_combo_n_abnormal %in% 1:3 ~ as.character(hcm_combo_n_abnormal),
      TRUE ~ NA_character_
    )
  ) %>%
  distinct(ID, .keep_all = TRUE)

longitudinal_lvdd_df <- long_df %>%
  left_join(lvdd_burden_lookup, by = "ID") %>%
  filter(
    !is.na(lvdd_burden),
    !is.na(time_yrs),
    !is.na(age),
    !is.na(BMI),
    is.finite(time_yrs),
    time_yrs <= 10
  ) %>%
  mutate(
    ID_fac = factor(ID),
    Sex = factor(Sex, levels = c("Male", "Female")),
    lvdd_burden = factor(lvdd_burden, levels = c("1", "2", "3"))
  )

binary_baseline_lookup <- baseline_df %>%
  mutate(
    abn_ee = case_when(!is.na(e_e_ave) ~ e_e_ave > 14, TRUE ~ NA),
    abn_lavi = case_when(
      !is.na(la_vol_index) ~ la_vol_index > 34,
      TRUE ~ NA
    ),
    abn_trv = case_when(!is.na(tr_max_vel) ~ tr_max_vel > 280, TRUE ~ NA)
  ) %>%
  select(ID, abn_ee, abn_lavi, abn_trv) %>%
  distinct(ID, .keep_all = TRUE)

longitudinal_binary_df <- long_df %>%
  left_join(binary_baseline_lookup, by = "ID") %>%
  filter(
    !is.na(time_yrs),
    !is.na(age),
    !is.na(BMI),
    is.finite(time_yrs),
    time_yrs <= 10
  ) %>%
  mutate(
    ID_fac = factor(ID),
    Sex = factor(Sex, levels = c("Male", "Female"))
  )

continuous_baseline_lookup <- baseline_df %>%
  transmute(
    ID,
    baseline_e_e_ave = e_e_ave,
    baseline_la_vol_index = la_vol_index,
    baseline_tr_max_vel = tr_max_vel
  ) %>%
  distinct(ID, .keep_all = TRUE)

longitudinal_continuous_df <- long_df %>%
  left_join(continuous_baseline_lookup, by = "ID") %>%
  filter(
    !is.na(time_yrs),
    !is.na(age),
    !is.na(BMI),
    is.finite(time_yrs),
    time_yrs <= 10
  ) %>%
  mutate(
    ID_fac = factor(ID),
    Sex = factor(Sex, levels = c("Male", "Female"))
  )

figure4_core_covariates <- c("age", "Sex", "BMI")
figure4_medication_covariates <- c("bb_any", "ndhp_ccb_any")
figure4_additional_sensitivity_covariates <- c("dm_pre_test")
figure4_medication_adjusted_covariates <- c(
  figure4_core_covariates,
  figure4_medication_covariates
)
figure4_full_covariates <- c(
  figure4_medication_adjusted_covariates,
  figure4_additional_sensitivity_covariates
)


figure4_index_specs <- list(
  list(
    var = "e_e_ave",
    baseline_var = "baseline_e_e_ave",
    status_var = "abn_ee",
    label_short = "E/e'",
    label_long = "E/e' (average)",
    threshold_subtitle = "E/e' ≤ 14 versus above 14",
    status_reference_label = "≤ 14",
    status_high_label = "Above 14",
    suffix = "Ee"
  ),
  list(
    var = "la_vol_index",
    baseline_var = "baseline_la_vol_index",
    status_var = "abn_lavi",
    label_short = "LAVi",
    label_long = "LAVi (mL/m\u00B2)",
    threshold_subtitle = "LAVi ≤ 34 mL/m\u00B2 versus above 34 mL/m\u00B2",
    status_reference_label = "≤ 34",
    status_high_label = "Above 34",
    suffix = "LAVi"
  ),
  list(
    var = "tr_max_vel",
    baseline_var = "baseline_tr_max_vel",
    status_var = "abn_trv",
    label_short = "TRVmax",
    label_long = "TRVmax (m/s)",
    threshold_subtitle = "TRVmax ≤ 2.8 m/s versus above 2.8 m/s",
    status_reference_label = "≤ 2.8",
    status_high_label = "Above 2.8",
    suffix = "TRV"
  )
)

# The medication and diabetes tiers model covariate missingness with
# three-level (No / Yes / Unknown) factors instead of requiring recorded
# status, which previously dropped roughly 60% of observations and made the
# unadjusted-versus-adjusted comparison a mixture of confounding adjustment
# and cohort restriction (owner decision, card 1.5 of
# 8_Docs/DECISION_CARDS_2026-09.md). The covariate set is unchanged; only the
# handling of unrecorded values differs. Table_Sensitivity_MissingIndicator_LMM
# carries the explicit complete-case contrast.
figure4_binary_indicator_vars <- c(
  figure4_medication_covariates,
  figure4_additional_sensitivity_covariates
)
# Dummy form here: fit_gamm_binary_panel() coerces its covariates to numeric.
figure4_binary_indicator_build <- add_missing_indicator_covariates(
  longitudinal_binary_df,
  binary_vars = figure4_binary_indicator_vars,
  binary_style = "dummy"
)
longitudinal_binary_indicator_df <- figure4_binary_indicator_build$data
figure4_indicator_name <- function(vars) {
  as.vector(rbind(paste0(vars, "_yes"), paste0(vars, "_unknown")))
}

figure4_model_sets <- list(
  list(key = "unadjusted", label = "Unadjusted", covariates = character(0)),
  list(
    key = "demographic",
    label = "Age, sex, BMI adjusted",
    covariates = figure4_core_covariates
  ),
  list(
    key = "medication",
    label = "Age, sex, BMI + medication adjusted",
    covariates = c(
      figure4_core_covariates,
      figure4_indicator_name(figure4_medication_covariates)
    )
  ),
  list(
    key = "full",
    label = "Age, sex, BMI + medication + diabetes adjusted",
    covariates = c(
      figure4_core_covariates,
      figure4_indicator_name(figure4_binary_indicator_vars)
    )
  )
)

figure4_model_results <- setNames(
  lapply(figure4_model_sets, function(x) {
    list(VO2_FRIEND2_PP = list(), VeVco2_slope = list())
  }),
  vapply(figure4_model_sets, function(x) x$key, character(1))
)

for (model_set in figure4_model_sets) {
  for (spec in figure4_index_specs) {
    figure4_model_results[[model_set$key]][["VO2_FRIEND2_PP"]][[spec$var]] <-
      fit_gamm_binary_panel(
        longitudinal_binary_indicator_df,
        "VO2_FRIEND2_PP",
        spec,
        covariates = model_set$covariates,
        stats_tag = model_set$key
      )
    figure4_model_results[[model_set$key]][["VeVco2_slope"]][[spec$var]] <-
      fit_gamm_binary_panel(
        longitudinal_binary_indicator_df,
        "VeVco2_slope",
        spec,
        covariates = model_set$covariates,
        stats_tag = model_set$key
      )
  }
}

figure4_results <- figure4_model_results[["demographic"]]
figure4_unadjusted_results <- figure4_model_results[["unadjusted"]]
figure4_model_label_order <- vapply(
  figure4_model_sets,
  function(x) x$label,
  character(1)
)
figure4_main_spec <- figure4_index_specs[[which(
  vapply(figure4_index_specs, function(spec) spec$var, character(1)) ==
    "la_vol_index"
)]]

figure4_main_results <- list(
  VO2_FRIEND2_PP = setNames(
    lapply(figure4_index_specs, function(spec) {
      fit_gamm_continuous_panel(
        longitudinal_continuous_df,
        "VO2_FRIEND2_PP",
        spec,
        covariates = figure4_core_covariates
      )
    }),
    vapply(figure4_index_specs, function(spec) spec$var, character(1))
  ),
  VeVco2_slope = setNames(
    lapply(figure4_index_specs, function(spec) {
      fit_gamm_continuous_panel(
        longitudinal_continuous_df,
        "VeVco2_slope",
        spec,
        covariates = figure4_core_covariates
      )
    }),
    vapply(figure4_index_specs, function(spec) spec$var, character(1))
  )
)

figure4_main_vo2_limits <- collect_prediction_limits(figure4_main_results[[
  "VO2_FRIEND2_PP"
]])
figure4_main_vevco2_limits <- collect_prediction_limits(figure4_main_results[[
  "VeVco2_slope"
]])

if (
  all(vapply(figure4_main_results[["VO2_FRIEND2_PP"]], is.null, logical(1))) &&
    all(vapply(figure4_main_results[["VeVco2_slope"]], is.null, logical(1)))
) {
  stop("Figure 4 continuous trajectory panels could not be generated.")
}

figure4_panel_titles <- c("A", "B", "C", "D", "E", "F")

figure4_vo2_panels <- lapply(seq_along(figure4_index_specs), function(i) {
  spec <- figure4_index_specs[[i]]
  make_continuous_trajectory_panel(
    figure4_main_results[["VO2_FRIEND2_PP"]][[spec$var]],
    y_label = if (i == 1) label_peak_vo2_friend_pred_unicode else NULL,
    y_limits = figure4_main_vo2_limits,
    title = paste0("(", figure4_panel_titles[[i]], ") ", spec$label_short),
    show_y = i == 1,
    show_x = FALSE,
    show_legend = i == 1,
    stats_label = format_continuous_panel_stats_label(
      figure4_main_results[["VO2_FRIEND2_PP"]][[spec$var]]
    ),
    legend_position = if (i == 1) c(0.04, 0.06) else "right",
    legend_direction = "vertical",
    legend_justification = if (i == 1) c(0, 0) else NULL,
    stats_fill_alpha = 0.62
  )
})

figure4_vevco2_panels <- lapply(seq_along(figure4_index_specs), function(i) {
  spec <- figure4_index_specs[[i]]
  make_continuous_trajectory_panel(
    figure4_main_results[["VeVco2_slope"]][[spec$var]],
    y_label = if (i == 1) label_vevco2_unicode else NULL,
    y_limits = figure4_main_vevco2_limits,
    title = paste0("(", figure4_panel_titles[[i + 3]], ") ", spec$label_short),
    show_y = i == 1,
    show_x = TRUE,
    show_legend = FALSE,
    stats_label = format_continuous_panel_stats_label(
      figure4_main_results[["VeVco2_slope"]][[spec$var]]
    ),
    legend_position = "right",
    legend_direction = "vertical",
    stats_fill_alpha = 0.62
  )
})

fig4 <- wrap_plots(
  c(figure4_vo2_panels, figure4_vevco2_panels),
  ncol = 3,
  byrow = TRUE
)

figure4_main_width <- 7.45
figure4_main_height <- 4.15
figure4_preprint_width <- 7.45
figure4_preprint_height <- 4.30
figure4_legend_fontsize <- 9

show_and_save_jacc(
  fig4,
  file.path(OUT_DIR, "Figure4_Longitudinal_GAMM.pdf"),
  w = figure4_main_width,
  h = figure4_main_height
)
save_jacc(
  fig4,
  file.path(OUT_DIR, "Figure4_Longitudinal_GAMM_Preprint.pdf"),
  w = figure4_preprint_width,
  h = figure4_preprint_height
)
save_jacc(
  fig4,
  file.path(OUT_DIR, "Figure4_Longitudinal_GAMM_Preprint.pdf"),
  w = figure4_preprint_width,
  h = figure4_preprint_height
)

legend_fig4_pdf <- paste0(
  "**Figure 4. Longitudinal Cardiopulmonary Trajectories Across LV Diastolic Parameters.** ",
  "Age-, sex-, and BMI-adjusted generalized additive mixed models related each continuous baseline parameter to repeated measures over time; ",
  "lines show predictions at the medians of the lower and upper halves. ",
  "**(A)** ",
  label_peak_vo2_friend_pred_md,
  " trajectory across baseline E/e’. ",
  "**(B)** ",
  label_peak_vo2_friend_pred_md,
  " trajectory across baseline LAVi. ",
  "**(C)** ",
  label_peak_vo2_friend_pred_md,
  " trajectory across baseline TRV$_max$. ",
  "**(D)** ",
  label_vevco2_md,
  " trajectory across baseline E/e’. ",
  "**(E)** ",
  label_vevco2_md,
  " trajectory across baseline LAVi. ",
  "**(F)** ",
  label_vevco2_md,
  " trajectory across baseline TRV$_max$."
)

save_jacc_with_embedded_legend(
  fig4,
  file.path(OUT_DIR, "Figure4_Longitudinal_GAMM.pdf"),
  legend_text = legend_fig4_pdf,
  w = figure4_main_width,
  h = figure4_main_height,
  legend_fontsize = figure4_legend_fontsize
)
save_jacc_with_embedded_legend(
  fig4,
  file.path(OUT_DIR, "Figure4_Longitudinal_GAMM.pdf"),
  legend_text = legend_fig4_pdf,
  w = figure4_main_width,
  h = figure4_main_height,
  legend_fontsize = figure4_legend_fontsize
)

# ══ Model Summaries ═══════════════════════════════════════════════════════════
# (legacy chunk: figure4_summary_table)
# ── table: GAMM smooth-term stats per index (primary, core-adjusted models) ────
figure4_summary_table_raw <- bind_rows(
  bind_rows(lapply(figure4_index_specs, function(spec) {
    extract_continuous_model_rows(
      figure4_main_results[["VO2_FRIEND2_PP"]][[spec$var]],
      label_peak_vo2_friend_pred
    )
  })),
  bind_rows(lapply(figure4_index_specs, function(spec) {
    extract_continuous_model_rows(
      figure4_main_results[["VeVco2_slope"]][[spec$var]],
      label_vevco2
    )
  }))
)

figure4_summary_table_display <- figure4_summary_table_raw %>%
  transmute(
    Outcome = case_when(
      Outcome == label_peak_vo2_friend_pred ~ label_peak_vo2,
      TRUE ~ Outcome
    ),
    `Baseline Parameter` = Parameter,
    `GAMM Term` = `Model component`,
    `Effective df` = ifelse(is.na(edf), "\u2014", sprintf("%.2f", edf)),
    Statistic = `Statistic`,
    `P value` = `P value`,
    p_num
  )

outcome_break_rows <- figure4_summary_table_display %>%
  count(Outcome, name = "n_rows") %>%
  pull(n_rows) %>%
  cumsum()
parameter_break_rows <- figure4_summary_table_display %>%
  count(Outcome, `Baseline Parameter`, name = "n_rows") %>%
  pull(n_rows) %>%
  cumsum()

figure4_summary_table <- figure4_summary_table_display %>%
  select(-p_num)

figure4_ft <- figure4_summary_table_display %>%
  flextable(
    col_keys = c(
      "Outcome",
      "Baseline Parameter",
      "GAMM Term",
      "Effective df",
      "Statistic",
      "P value"
    )
  ) %>%
  merge_v(j = c("Outcome", "Baseline Parameter")) %>%
  valign(
    j = c("Outcome", "Baseline Parameter"),
    valign = "top",
    part = "body"
  ) %>%
  align(
    j = c("Outcome", "Baseline Parameter", "GAMM Term"),
    align = "left",
    part = "all"
  ) %>%
  align(
    j = c("Effective df", "Statistic", "P value"),
    align = "center",
    part = "all"
  ) %>%
  bold(i = ~ p_num < 0.05, j = "P value", part = "body") %>%
  format_pub_table(
    caption = "Table 3. Key longitudinal GAMM terms for Figure 4 continuous parameter-specific analyses"
  ) %>%
  hline(
    i = parameter_break_rows,
    border = officer::fp_border(color = "grey55", width = 0.4),
    part = "body"
  ) %>%
  hline(
    i = outcome_break_rows,
    border = officer::fp_border(color = "black", width = 0.8),
    part = "body"
  )

figure4_ft <- apply_subscript_format(
  figure4_ft,
  figure4_summary_table_display,
  "Outcome"
)
figure4_ft <- apply_scientific_p_format(
  figure4_ft,
  figure4_summary_table_display,
  value_col = "P value",
  p_col = "p_num",
  digits = 2,
  threshold = 0.05
)


write.csv(
  figure4_summary_table,
  file.path(OUT_DIR, "Table3_Figure4_GAMM_Summary.csv"),
  row.names = FALSE
)
figure4_ft %>%
  save_as_docx(path = file.path(OUT_DIR, "Table3_Figure4_GAMM_Summary.docx"))

# ══ Complete Model Statistics ═════════════════════════════════════════════════
# (legacy chunk: figure4_complete_table)
# ── table: full GAMM coefficients (parametric + smooth) across all indices ─────
figure4_complete_table_raw <- bind_rows(
  bind_rows(lapply(figure4_index_specs, function(spec) {
    extract_complete_gamm_rows(
      figure4_results[["VO2_FRIEND2_PP"]][[spec$var]],
      label_peak_vo2_friend_pred
    )
  })),
  bind_rows(lapply(figure4_index_specs, function(spec) {
    extract_complete_gamm_rows(
      figure4_results[["VeVco2_slope"]][[spec$var]],
      label_vevco2
    )
  }))
) %>%
  mutate(
    outcome_order = match(Outcome, c(label_peak_vo2_friend_pred, label_vevco2)),
    parameter_order = match(
      Parameter,
      vapply(figure4_index_specs, function(spec) spec$label_long, character(1))
    ),
    term_class_order = match(`Term class`, c("Smooth", "Parametric"))
  ) %>%
  arrange(outcome_order, parameter_order, term_class_order, Term) %>%
  select(-outcome_order, -parameter_order, -term_class_order)

figure4_complete_table_display <- figure4_complete_table_raw %>%
  transmute(
    Outcome = case_when(
      Outcome == label_peak_vo2_friend_pred ~ label_peak_vo2,
      TRUE ~ Outcome
    ),
    `Baseline Parameter` = Parameter,
    Sample = sprintf("%d / %d", Observations, Patients),
    Fit = `Fit mode`,
    `Term Class` = `Term class`,
    Term,
    `Effect summary` = case_when(
      `Term class` == "Smooth" ~ sprintf(
        "edf = %.2f; ref df = %.2f",
        edf,
        `Ref df`
      ),
      TRUE ~ sprintf("\u03B2 = %.2f; SE = %.2f", Estimate, `Std. Error`)
    ),
    Statistic = sprintf("%s = %.2f", `Statistic type`, `Statistic value`),
    `P value` = ifelse(
      `P value raw` < 0.001,
      "<0.001",
      sprintf("%.3f", `P value raw`)
    ),
    p_num = `P value raw`
  )

figure4_complete_term_break_rows <- figure4_complete_table_display %>%
  count(Outcome, `Baseline Parameter`, `Term Class`, name = "n_rows") %>%
  pull(n_rows) %>%
  cumsum()
figure4_complete_parameter_break_rows <- figure4_complete_table_display %>%
  count(Outcome, `Baseline Parameter`, name = "n_rows") %>%
  pull(n_rows) %>%
  cumsum()
figure4_complete_outcome_break_rows <- figure4_complete_table_display %>%
  count(Outcome, name = "n_rows") %>%
  pull(n_rows) %>%
  cumsum()

figure4_complete_table <- figure4_complete_table_display %>%
  select(-p_num)

figure4_complete_ft <- figure4_complete_table_display %>%
  flextable(
    col_keys = c(
      "Outcome",
      "Baseline Parameter",
      "Sample",
      "Fit",
      "Term Class",
      "Term",
      "Effect summary",
      "Statistic",
      "P value"
    )
  ) %>%
  merge_v(
    j = c("Outcome", "Baseline Parameter", "Sample", "Fit", "Term Class")
  ) %>%
  valign(
    j = c("Outcome", "Baseline Parameter", "Sample", "Fit", "Term Class"),
    valign = "top",
    part = "body"
  ) %>%
  align(
    j = c(
      "Outcome",
      "Baseline Parameter",
      "Fit",
      "Term Class",
      "Term",
      "Effect summary"
    ),
    align = "left",
    part = "all"
  ) %>%
  align(
    j = c("Sample", "Statistic", "P value"),
    align = "center",
    part = "all"
  ) %>%
  bold(i = ~ p_num < 0.05, j = "P value", part = "body") %>%
  format_pub_table(
    caption = "Supplemental Table. Complete longitudinal GAMM statistics underlying the supplemental parameter-specific Figure 4 binary-status analyses"
  ) %>%
  hline(
    i = figure4_complete_term_break_rows,
    border = officer::fp_border(color = "grey75", width = 0.3),
    part = "body"
  ) %>%
  hline(
    i = figure4_complete_parameter_break_rows,
    border = officer::fp_border(color = "grey55", width = 0.4),
    part = "body"
  ) %>%
  hline(
    i = figure4_complete_outcome_break_rows,
    border = officer::fp_border(color = "black", width = 0.8),
    part = "body"
  )

figure4_complete_ft <- apply_subscript_format(
  figure4_complete_ft,
  figure4_complete_table_display,
  "Outcome"
)


write.csv(
  figure4_complete_table,
  file.path(OUT_DIR, "TableS_Figure4_GAMM_Complete.csv"),
  row.names = FALSE
)
figure4_complete_ft %>%
  save_as_docx(path = file.path(OUT_DIR, "TableS_Figure4_GAMM_Complete.docx"))

# (legacy chunk: figure4_validation)
# ── Fig 4 sanity check: GAMM N, EDF, time-smooth P per panel ──────────────────
figure4_validation_row <- function(result, outcome_label, parameter_label) {
  if (is.null(result)) {
    return(NULL)
  }
  tibble(
    Outcome = outcome_label,
    Parameter = parameter_label,
    Observations = result$n_observations,
    Patients = result$n_patients,
    Quartiles = paste0(
      result$quartile_summary$quartile,
      ": ",
      result$quartile_summary$n_patients,
      collapse = ", "
    )
  )
}

figure4_main_validation_df <- bind_rows(
  bind_rows(lapply(figure4_index_specs, function(spec) {
    figure4_validation_row(
      figure4_main_results[["VO2_FRIEND2_PP"]][[spec$var]],
      label_peak_vo2,
      spec$label_short
    )
  })),
  bind_rows(lapply(figure4_index_specs, function(spec) {
    figure4_validation_row(
      figure4_main_results[["VeVco2_slope"]][[spec$var]],
      label_vevco2,
      spec$label_short
    )
  }))
)

cat("**Figure 4 Validation (Main Continuous Quartile Analysis)**\n\n")
for (i in seq_len(nrow(figure4_main_validation_df))) {
  row_i <- figure4_main_validation_df[i, ]
  cat(sprintf(
    "- %s | %s: **%d observations across %d patients**; quartile patients = **%s**\n",
    row_i$Outcome,
    row_i$Parameter,
    row_i$Observations,
    row_i$Patients,
    row_i$Quartiles
  ))
}
cat(sprintf(
  "- Baseline LAVi Z-score derivation available for downstream Cox models: **%s**\n",
  ifelse(nrow(z_source) > 0, "yes", "no")
))

# ══ Linear Mixed-Model Confirmation ═══════════════════════════════════════════
# (legacy chunk: figure4_linear_mixed_confirmation)
# A parsimonious random-intercept model confirms the direction and significance
# of each time-by-baseline-parameter relationship without estimating nonlinear
# smooths. Baseline parameters are standardized across patients (not rows).
extract_continuous_gamm_interaction_p <- function(result) {
  if (is.null(result)) return(NA_real_)
  baseline_var <- result$spec$baseline_var %||% result$spec$var
  smooth_terms <- c(
    paste0("ti(time_yrs,", baseline_var, ")"),
    paste0("ti(", baseline_var, ",time_yrs)")
  )
  p_value <- extract_gamm_smooth_p(result, smooth_terms)
  if (is.na(p_value)) {
    p_value <- extract_gamm_parametric_p(
      result,
      c(
        paste0(baseline_var, ":time_yrs"),
        paste0("time_yrs:", baseline_var)
      )
    )
  }
  p_value
}

fit_linear_mixed_confirmation <- function(
  data,
  outcome,
  spec,
  random_structure = c("intercept", "intercept_slope"),
  extra_covariates = character(0)
) {
  random_structure <- match.arg(random_structure)
  baseline_var <- spec$baseline_var
  required <- c(
    "ID",
    "ID_fac",
    "time_yrs",
    outcome,
    baseline_var,
    figure4_core_covariates,
    extra_covariates
  )
  model_df <- data %>%
    select(all_of(required)) %>%
    drop_na() %>%
    filter(is.finite(time_yrs), is.finite(.data[[outcome]])) %>%
    group_by(ID) %>%
    filter(n_distinct(time_yrs) >= 2) %>%
    ungroup()

  patient_baseline_df <- model_df %>% distinct(ID, .keep_all = TRUE)
  patient_baseline <- patient_baseline_df[[baseline_var]]
  baseline_mean <- mean(patient_baseline, na.rm = TRUE)
  baseline_sd <- sd(patient_baseline, na.rm = TRUE)
  if (
    nrow(model_df) < 50 ||
      n_distinct(model_df$ID) < 20 ||
      !is.finite(baseline_sd) ||
      baseline_sd <= 0
  ) {
    return(NULL)
  }

  model_df <- model_df %>%
    mutate(
      baseline_parameter_z = (.data[[baseline_var]] - baseline_mean) /
        baseline_sd,
      ID_fac = droplevels(factor(ID_fac)),
      Sex = droplevels(factor(Sex, levels = c("Male", "Female")))
    )

  model_formula <- reformulate(
    c(
      "time_yrs * baseline_parameter_z",
      figure4_core_covariates,
      drop_degenerate_covariates(model_df, extra_covariates)
    ),
    response = outcome
  )
  random_formula <- if (random_structure == "intercept_slope") {
    list(ID_fac = nlme::pdDiag(~ 1 + time_yrs))
  } else {
    ~ 1 | ID_fac
  }
  fit <- tryCatch(
    nlme::lme(
      fixed = model_formula,
      random = random_formula,
      data = model_df,
      method = "REML",
      na.action = na.omit,
      control = nlme::lmeControl(
        opt = "optim",
        maxIter = 200,
        msMaxIter = 200
      )
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  coefficient_table <- as.data.frame(summary(fit)$tTable)
  coefficient_table$term <- rownames(coefficient_table)
  interaction_row <- coefficient_table %>%
    filter(term %in% c(
      "time_yrs:baseline_parameter_z",
      "baseline_parameter_z:time_yrs"
    )) %>%
    slice(1)
  if (nrow(interaction_row) != 1) return(NULL)

  critical_value <- qt(0.975, df = interaction_row$DF)
  list(
    fit = fit,
    formula = model_formula,
    data = model_df,
    estimate = interaction_row$Value,
    std_error = interaction_row$Std.Error,
    df = interaction_row$DF,
    ci_low = interaction_row$Value - critical_value * interaction_row$Std.Error,
    ci_high = interaction_row$Value + critical_value * interaction_row$Std.Error,
    p_value = interaction_row$`p-value`,
    baseline_mean = baseline_mean,
    baseline_sd = baseline_sd,
    random_structure = random_structure
  )
}

figure4_lmm_confirmation <- map_dfr(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome) {
    map_dfr(figure4_index_specs, function(spec) {
      fit_obj <- fit_linear_mixed_confirmation(
        longitudinal_continuous_df,
        outcome,
        spec
      )
      if (is.null(fit_obj)) return(tibble())
      gamm_p <- extract_continuous_gamm_interaction_p(
        figure4_main_results[[outcome]][[spec$var]]
      )
      tibble(
        Outcome = ifelse(
          outcome == "VO2_FRIEND2_PP",
          "Peak VO2 (% predicted)",
          "VE/VCO2 slope"
        ),
        Parameter = spec$label_short,
        Patients = n_distinct(fit_obj$data$ID),
        Observations = nrow(fit_obj$data),
        `Baseline mean` = fit_obj$baseline_mean,
        `Baseline SD` = fit_obj$baseline_sd,
        `Time x parameter beta per year per 1 SD` = fit_obj$estimate,
        `CI low` = fit_obj$ci_low,
        `CI high` = fit_obj$ci_high,
        `LMM P value` = fit_obj$p_value,
        `GAMM interaction P value` = gamm_p,
        `Significance classification concordant` = ifelse(
          is.na(gamm_p),
          NA,
          (fit_obj$p_value < 0.05) == (gamm_p < 0.05)
        )
      )
    })
  }
)

# Each parameter now has its own available-case sample, so denominators differ
# between parameters by design (the previous audit required E/e' and LAVi to
# share one denominator, which was a property of the common cohort). Verify
# instead that every reported denominator is exactly the complete-case count
# implied by its outcome, baseline parameter, and core covariates, and that no
# sample exceeds the base longitudinal cohort.
figure4_lmm_expected_n <- map_dfr(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome) {
    map_dfr(figure4_index_specs, function(spec) {
      expected <- longitudinal_continuous_df %>%
        select(all_of(c(
          "ID", "ID_fac", "time_yrs", outcome, spec$baseline_var,
          figure4_core_covariates
        ))) %>%
        drop_na() %>%
        filter(is.finite(time_yrs), is.finite(.data[[outcome]])) %>%
        group_by(ID) %>%
        filter(n_distinct(time_yrs) >= 2) %>%
        ungroup()
      tibble(
        Outcome = ifelse(
          outcome == "VO2_FRIEND2_PP",
          "Peak VO2 (% predicted)",
          "VE/VCO2 slope"
        ),
        Parameter = spec$label_short,
        expected_patients = n_distinct(expected$ID),
        expected_observations = nrow(expected)
      )
    })
  }
)

figure4_lmm_denominator_audit <- figure4_lmm_confirmation %>%
  select(Outcome, Parameter, Patients, Observations) %>%
  left_join(figure4_lmm_expected_n, by = c("Outcome", "Parameter")) %>%
  mutate(
    denominator_matches = Patients == expected_patients &
      Observations == expected_observations
  )
stopifnot(
  all(figure4_lmm_denominator_audit$denominator_matches),
  all(figure4_lmm_confirmation$Patients <=
    n_distinct(longitudinal_continuous_df$ID))
)
write.csv(
  figure4_lmm_denominator_audit,
  file.path(OUT_DIR, "Table_Longitudinal_Denominator_Audit.csv"),
  row.names = FALSE
)

write.csv(
  figure4_lmm_confirmation,
  paste0(
    file.path(OUT_DIR, ""),
    "Table4_LinearMixedModel_Confirmation.csv"
  ),
  row.names = FALSE
)

# Random-structure audit. The random-slope candidate is diagonal (pdDiag),
# avoiding an unsupported intercept-slope correlation. Models are compared by
# ML with identical fixed effects. The LRT is descriptive because the null
# variance lies on the boundary; AIC and convergence are reported alongside it.
fit_lmm_random_structure_audit <- function(fit_obj) {
  warning_messages <- character(0)
  fit_with_warnings <- function(...) {
    withCallingHandlers(
      tryCatch(nlme::lme(...), error = function(e) NULL),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }

  control <- nlme::lmeControl(
    opt = "optim",
    maxIter = 300,
    msMaxIter = 300,
    niterEM = 100,
    tolerance = 1e-7,
    msTol = 1e-7
  )
  ri_ml <- fit_with_warnings(
    fixed = fit_obj$formula,
    random = ~ 1 | ID_fac,
    data = fit_obj$data,
    method = "ML",
    na.action = na.omit,
    control = control
  )
  rs_ml <- fit_with_warnings(
    fixed = fit_obj$formula,
    random = list(ID_fac = nlme::pdDiag(~ 1 + time_yrs)),
    data = fit_obj$data,
    method = "ML",
    na.action = na.omit,
    control = control
  )

  if (is.null(ri_ml)) {
    return(list(
      ri = NULL,
      rs = rs_ml,
      warning = paste(unique(warning_messages), collapse = " | ")
    ))
  }

  extract_interaction <- function(model) {
    if (is.null(model)) {
      return(tibble(beta = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
    }
    tbl <- as.data.frame(summary(model)$tTable) %>%
      tibble::rownames_to_column("term") %>%
      filter(term %in% c(
        "time_yrs:baseline_parameter_z",
        "baseline_parameter_z:time_yrs"
      )) %>%
      slice(1)
    if (nrow(tbl) != 1) {
      return(tibble(beta = NA_real_, ci_low = NA_real_, ci_high = NA_real_, p = NA_real_))
    }
    critical <- qt(0.975, df = tbl$DF)
    tibble(
      beta = tbl$Value,
      ci_low = tbl$Value - critical * tbl$Std.Error,
      ci_high = tbl$Value + critical * tbl$Std.Error,
      p = tbl$`p-value`
    )
  }

  rs_vc <- if (!is.null(rs_ml)) nlme::VarCorr(rs_ml) else NULL
  slope_sd <- if (!is.null(rs_vc) && "time_yrs" %in% rownames(rs_vc)) {
    suppressWarnings(as.numeric(rs_vc["time_yrs", "StdDev"]))
  } else {
    NA_real_
  }
  residual_sd <- if (!is.null(rs_vc) && "Residual" %in% rownames(rs_vc)) {
    suppressWarnings(as.numeric(rs_vc["Residual", "StdDev"]))
  } else {
    NA_real_
  }

  comparison <- if (!is.null(rs_ml)) {
    tryCatch(anova(ri_ml, rs_ml), error = function(e) NULL)
  } else {
    NULL
  }
  lrt_p <- if (!is.null(comparison) && nrow(comparison) >= 2) {
    as.numeric(comparison$`p-value`[2])
  } else {
    NA_real_
  }

  list(
    ri = ri_ml,
    rs = rs_ml,
    ri_effect = extract_interaction(ri_ml),
    rs_effect = extract_interaction(rs_ml),
    slope_sd = slope_sd,
    residual_sd = residual_sd,
    lrt_p = lrt_p,
    warning = paste(unique(warning_messages), collapse = " | ")
  )
}

figure4_random_structure_audit <- map_dfr(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome) {
    map_dfr(figure4_index_specs, function(spec) {
      fit_obj <- fit_linear_mixed_confirmation(
        longitudinal_continuous_df,
        outcome,
        spec
      )
      if (is.null(fit_obj)) return(tibble())
      audit <- fit_lmm_random_structure_audit(fit_obj)
      rs_converged <- !is.null(audit$rs)
      ri_aic <- if (!is.null(audit$ri)) AIC(audit$ri) else NA_real_
      rs_aic <- if (rs_converged) AIC(audit$rs) else NA_real_
      delta_aic <- ri_aic - rs_aic
      stable_slope <- rs_converged &&
        is.finite(audit$slope_sd) &&
        audit$slope_sd > 1e-6 &&
        is.finite(audit$residual_sd) &&
        !nzchar(audit$warning)
      supports_slope <- stable_slope &&
        is.finite(delta_aic) &&
        delta_aic >= 2 &&
        is.finite(audit$lrt_p) &&
        audit$lrt_p < 0.05

      tibble(
        Outcome = ifelse(
          outcome == "VO2_FRIEND2_PP",
          "Peak VO2 (% predicted)",
          "VE/VCO2 slope"
        ),
        Parameter = spec$label_short,
        Patients = n_distinct(fit_obj$data$ID),
        Observations = nrow(fit_obj$data),
        `Random-intercept AIC` = ri_aic,
        `Random-slope AIC` = rs_aic,
        `Delta AIC (RI - RS)` = delta_aic,
        `Descriptive LRT P` = audit$lrt_p,
        `Random-slope SD` = audit$slope_sd,
        `Residual SD` = audit$residual_sd,
        `RI interaction beta` = audit$ri_effect$beta,
        `RI interaction CI low` = audit$ri_effect$ci_low,
        `RI interaction CI high` = audit$ri_effect$ci_high,
        `RI interaction P` = audit$ri_effect$p,
        `RS interaction beta` = audit$rs_effect$beta,
        `RS interaction CI low` = audit$rs_effect$ci_low,
        `RS interaction CI high` = audit$rs_effect$ci_high,
        `RS interaction P` = audit$rs_effect$p,
        `Random slope converged` = rs_converged,
        `Random slope stable` = stable_slope,
        `Random slope supported` = supports_slope,
        Warnings = audit$warning
      )
    })
  }
)

primary_random_structure_rows <- figure4_random_structure_audit %>%
  filter(Parameter %in% c("E/e'", "LAVi"))
figure4_random_structure_decision <- primary_random_structure_rows %>%
  group_by(Outcome) %>%
  summarise(
    `Primary models meeting rule` = sum(`Random slope supported`),
    `Primary models evaluated` = n(),
    Decision = if_else(
      `Primary models evaluated` == 2 &
        `Primary models meeting rule` == 2,
      "Random intercept plus independent random slope",
      "Random intercept"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Rule = paste0(
      "Within each outcome, use random slopes only if both primary E/e' and ",
      "LAVi models converge without warnings, estimate a nonzero slope SD, ",
      "improve AIC by at least 2, and have descriptive LRT P<0.05."
    )
  ) %>%
  select(
    Outcome,
    Decision,
    Rule,
    `Primary models meeting rule`,
    `Primary models evaluated`
  )


# The primary trajectory estimates use the selected structure, so every
# comparison in this stage must use the same one. Assembling it numerically
# here -- the published table formats its P values as text -- gives the
# sensitivity analyses a baseline that is the published primary rather than a
# differently specified model.
figure4_primary_numeric <- figure4_random_structure_audit %>%
  left_join(
    figure4_random_structure_decision %>% select(Outcome, Decision),
    by = "Outcome"
  ) %>%
  mutate(
    primary_structure = if_else(
      Decision == "Random intercept plus independent random slope",
      "intercept_slope",
      "intercept"
    ),
    primary_beta = if_else(
      primary_structure == "intercept_slope",
      `RS interaction beta`,
      `RI interaction beta`
    ),
    primary_p = if_else(
      primary_structure == "intercept_slope",
      `RS interaction P`,
      `RI interaction P`
    )
  ) %>%
  select(
    Outcome, Parameter, Patients, Observations,
    primary_structure, primary_beta, primary_p
  )
stopifnot(
  nrow(figure4_primary_numeric) == nrow(figure4_lmm_confirmation),
  !any(is.na(figure4_primary_numeric$primary_p)),
  n_distinct(figure4_primary_numeric$primary_structure) >= 1
)

figure4_primary_structure_for <- function(outcome_label) {
  st <- unique(
    figure4_primary_numeric$primary_structure[
      figure4_primary_numeric$Outcome == outcome_label
    ]
  )
  stopifnot(length(st) == 1)
  st
}

# ── Medication adjustment: complete case versus missing-indicator handling ────
# The sequential-adjustment ladder elsewhere in this stage requires recorded
# medication status, which drops roughly 60% of observations, so an
# unadjusted-versus-medication-adjusted comparison there mixes confounding
# adjustment with cohort restriction. Repeating the medication-adjusted
# trajectory model with three-level (No / Yes / Unknown) covariates keeps every
# patient, so the comparison isolates the adjustment (owner decision, card 1.5
# of 8_Docs/DECISION_CARDS_2026-09.md).
figure4_med_indicator_build <- add_missing_indicator_covariates(
  longitudinal_continuous_df,
  binary_vars = figure4_medication_covariates
)

figure4_medication_handling <- map_dfr(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome) {
    outcome_label <- ifelse(
      outcome == "VO2_FRIEND2_PP",
      "Peak VO2 (% predicted)",
      "VE/VCO2 slope"
    )
    map_dfr(figure4_index_specs, function(spec) {
      handling_specs <- list(
        list(
          label = "Complete case (medication recorded)",
          data = longitudinal_continuous_df,
          covariates = figure4_medication_covariates
        ),
        list(
          label = "Three-level No / Yes / Unknown",
          data = figure4_med_indicator_build$data,
          covariates = figure4_med_indicator_build$covariates
        )
      )
      map_dfr(handling_specs, function(handling) {
        fit_obj <- fit_linear_mixed_confirmation(
          handling$data,
          outcome,
          spec,
          random_structure = figure4_primary_structure_for(outcome_label),
          extra_covariates = handling$covariates
        )
        if (is.null(fit_obj)) return(tibble())
        tibble(
          Outcome = outcome_label,
          Parameter = spec$label_short,
          `Covariate handling` = handling$label,
          Patients = n_distinct(fit_obj$data$ID),
          Observations = nrow(fit_obj$data),
          `Interaction beta per year per 1 SD` = fit_obj$estimate,
          `CI low` = fit_obj$ci_low,
          `CI high` = fit_obj$ci_high,
          `P value` = fit_obj$p_value
        )
      })
    })
  }
)

# Modelling the missingness can only add observations, so each indicator-handled
# fit must be at least as large as its complete-case counterpart.
figure4_medication_handling_nesting <- figure4_medication_handling %>%
  select(Outcome, Parameter, `Covariate handling`, Observations) %>%
  tidyr::pivot_wider(
    names_from = `Covariate handling`,
    values_from = Observations
  )
stopifnot(
  nrow(figure4_medication_handling_nesting) > 0,
  all(
    figure4_medication_handling_nesting$`Three-level No / Yes / Unknown` >=
      figure4_medication_handling_nesting$`Complete case (medication recorded)`
  )
)

write.csv(
  figure4_medication_handling,
  file.path(OUT_DIR, "Table_Sensitivity_MissingIndicator_LMM.csv"),
  row.names = FALSE
)
print(as.data.frame(figure4_medication_handling), row.names = FALSE)

# ── Observed follow-up horizon per trajectory model ──────────────────────────
# The interaction is a per-year slope, so it is only supported over the window
# the patients were actually observed. The Results text has to state that
# window, and it differs by parameter and outcome because each model has its
# own available-case frame. Reported as the distribution of each patient's last
# included test.
figure4_followup_horizon <- map_dfr(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome) {
    outcome_label <- ifelse(
      outcome == "VO2_FRIEND2_PP",
      "Peak VO2 (% predicted)",
      "VE/VCO2 slope"
    )
    map_dfr(figure4_index_specs, function(spec) {
      fit_obj <- fit_linear_mixed_confirmation(
        longitudinal_continuous_df,
        outcome,
        spec,
        random_structure = figure4_primary_structure_for(outcome_label)
      )
      if (is.null(fit_obj)) return(tibble())
      per_patient <- fit_obj$data %>%
        group_by(ID) %>%
        summarise(last_yrs = max(time_yrs), .groups = "drop")
      tibble(
        Outcome = outcome_label,
        Parameter = spec$label_short,
        Patients = nrow(per_patient),
        Observations = nrow(fit_obj$data),
        `Follow-up median, years` = median(per_patient$last_yrs),
        `Follow-up Q1, years` = quantile(per_patient$last_yrs, 0.25),
        `Follow-up Q3, years` = quantile(per_patient$last_yrs, 0.75),
        `Follow-up maximum, years` = max(per_patient$last_yrs),
        `Patients followed >= 5 years` = sum(per_patient$last_yrs >= 5),
        `Patients followed >= 8 years` = sum(per_patient$last_yrs >= 8)
      )
    })
  }
)
stopifnot(
  nrow(figure4_followup_horizon) == nrow(figure4_primary_numeric),
  all(figure4_followup_horizon$`Follow-up maximum, years` <= 10 + 1e-8)
)
write.csv(
  figure4_followup_horizon,
  file.path(OUT_DIR, "Table_Longitudinal_FollowupHorizon.csv"),
  row.names = FALSE
)
print(as.data.frame(figure4_followup_horizon), row.names = FALSE)

# ── Consistency check: the same models in the common E/e' + LAVi cohort ───────
# Each parameter is analysed in its own available-case sample in the primary
# models, so denominators differ between parameters. This refits all three in
# the single cohort of patients who have both cohort-defining parameters
# measured, to show that analysing each parameter where it was measured does
# not change the trajectory conclusions (owner decision, card 2.3 of
# 8_Docs/DECISION_CARDS_2026-09.md).
figure4_common_cohort_ids <- baseline_df %>%
  filter(!is.na(e_e_ave), !is.na(la_vol_index)) %>%
  pull(ID)

figure4_common_continuous_df <- longitudinal_continuous_df %>%
  filter(ID %in% figure4_common_cohort_ids)

figure4_cohort_consistency <- map_dfr(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome) {
    outcome_label <- ifelse(
      outcome == "VO2_FRIEND2_PP",
      "Peak VO2 (% predicted)",
      "VE/VCO2 slope"
    )
    map_dfr(figure4_index_specs, function(spec) {
      primary_row <- figure4_primary_numeric %>%
        filter(Outcome == outcome_label, Parameter == spec$label_short)
      if (nrow(primary_row) != 1) return(tibble())
      common_fit <- fit_linear_mixed_confirmation(
        figure4_common_continuous_df,
        outcome,
        spec,
        random_structure = primary_row$primary_structure
      )
      if (is.null(common_fit)) return(tibble())
      tibble(
        Outcome = outcome_label,
        Parameter = spec$label_short,
        `Patients available-case` = primary_row$Patients,
        `Patients common cohort` = n_distinct(common_fit$data$ID),
        `Observations available-case` = primary_row$Observations,
        `Observations common cohort` = nrow(common_fit$data),
        `Interaction beta available-case` = primary_row$primary_beta,
        `Interaction beta common cohort` = common_fit$estimate,
        `P value available-case` = primary_row$primary_p,
        `P value common cohort` = common_fit$p_value,
        `Random-effects structure` = primary_row$primary_structure,
        `Conclusion unchanged` =
          (primary_row$primary_p < 0.05) == (common_fit$p_value < 0.05)
      )
    })
  }
)

# The common cohort is nested inside every parameter's available-case sample,
# so it can never be larger; a violation would mean the wrong frame was used.
stopifnot(
  nrow(figure4_cohort_consistency) == nrow(figure4_primary_numeric),
  all(figure4_cohort_consistency$`Patients common cohort` <=
    figure4_cohort_consistency$`Patients available-case`),
  !any(is.na(figure4_cohort_consistency$`P value common cohort`))
)

write.csv(
  figure4_cohort_consistency,
  file.path(OUT_DIR, "Table_Consistency_CommonCohort_LMM.csv"),
  row.names = FALSE
)
print(as.data.frame(figure4_cohort_consistency), row.names = FALSE)

# ── Sensitivity: multiple imputation of the missing exposures ─────────────────
# Labelled sensitivity analysis (owner decision, card 1.2). Read the header of
# mi_lmm_exposure_sensitivity() in R/imputation.R before quoting
# these numbers: there is no substantive-model-compatible imputation available
# for a mixed model, so this imputation model is not derived from the LMM and
# instead conditions on each patient's own outcome slope. TRVmax is never
# imputed as an exposure and appears only as a gated auxiliary.
if (isTRUE(cfg$imputation$enabled)) {
  figure4_mi_exposures <- list(
    list(var = "e_e_ave", baseline_var = "baseline_e_e_ave", label = "E/e'"),
    list(var = "la_vol_index", baseline_var = "baseline_la_vol_index", label = "LAVi")
  )

  figure4_mi_sensitivity <- map_dfr(
    c("VO2_FRIEND2_PP", "VeVco2_slope"),
    function(outcome) {
      outcome_label <- ifelse(
        outcome == "VO2_FRIEND2_PP",
        "Peak VO2 (% predicted)",
        "VE/VCO2 slope"
      )
      map_dfr(figure4_mi_exposures, function(mi_spec) {
        spec_template <- figure4_index_specs[[which(
          vapply(figure4_index_specs, function(x) x$var, character(1)) == mi_spec$var
        )]]
        reference <- figure4_primary_numeric %>%
          filter(Outcome == outcome_label, Parameter == spec_template$label_short)
        stopifnot(nrow(reference) == 1)
        mi_lmm_exposure_sensitivity(
          longitudinal_continuous_df,
          outcome_var = outcome,
          outcome_label = outcome_label,
          exposure = mi_spec$var,
          exposure_label = spec_template$label_short,
          baseline_var = mi_spec$baseline_var,
          auxiliary_baseline_vars = setdiff(
            c("baseline_e_e_ave", "baseline_la_vol_index", "baseline_tr_max_vel"),
            mi_spec$baseline_var
          ),
          fitter = fit_linear_mixed_confirmation,
          spec_template = spec_template,
          random_structure = reference$primary_structure,
          m = cfg$imputation$m,
          seed = cfg$imputation$seed,
          available_case = list(
            patients = reference$Patients,
            estimate = reference$primary_beta,
            p_value = reference$primary_p
          )
        )
      })
    }
  )

  stopifnot(
    nrow(figure4_mi_sensitivity) == 4,
    all(figure4_mi_sensitivity$`Exposure measured` +
      figure4_mi_sensitivity$`Exposure imputed` ==
      figure4_mi_sensitivity$`Patients in imputation frame`),
    all(figure4_mi_sensitivity$Imputations == cfg$imputation$m)
  )

  write.csv(
    figure4_mi_sensitivity,
    file.path(OUT_DIR, "Table_Sensitivity_MultipleImputation_LMM.csv"),
    row.names = FALSE
  )
  print(as.data.frame(figure4_mi_sensitivity), row.names = FALSE)
}

# ── Sensitivity: adults only ──────────────────────────────────────────────────
# Same rationale as the cross-sectional adults-only refit: the percent-
# predicted outcome rests on adult reference equations (owner decision, card
# 1.3). Patients below the age floor at their baseline CPET are excluded
# entirely, so a patient contributes all of their visits or none and the
# trajectory cohort stays internally consistent.
min_age_sensitivity <- cfg$cohort$min_age_sensitivity

if (!is.null(min_age_sensitivity)) {
  adult_baseline_ids <- baseline_df %>%
    filter(!is.na(age), age >= min_age_sensitivity) %>%
    pull(ID)

  figure4_adult_continuous_df <- longitudinal_continuous_df %>%
    filter(ID %in% adult_baseline_ids)

  figure4_age_sensitivity <- map_dfr(
    c("VO2_FRIEND2_PP", "VeVco2_slope"),
    function(outcome) {
      map_dfr(figure4_index_specs, function(spec) {
        outcome_label <- ifelse(
          outcome == "VO2_FRIEND2_PP",
          "Peak VO2 (% predicted)",
          "VE/VCO2 slope"
        )
        primary_row <- figure4_primary_numeric %>%
          filter(Outcome == outcome_label, Parameter == spec$label_short)
        if (nrow(primary_row) != 1) return(tibble())
        adult_fit <- fit_linear_mixed_confirmation(
          figure4_adult_continuous_df,
          outcome,
          spec,
          random_structure = primary_row$primary_structure
        )
        if (is.null(adult_fit)) return(tibble())
        tibble(
          Outcome = primary_row$Outcome,
          Parameter = spec$label_short,
          `Age floor (years)` = min_age_sensitivity,
          `Patients all ages` = primary_row$Patients,
          `Patients adults only` = n_distinct(adult_fit$data$ID),
          `Observations all ages` = primary_row$Observations,
          `Observations adults only` = nrow(adult_fit$data),
          `Interaction beta all ages` = primary_row$primary_beta,
          `Interaction beta adults only` = adult_fit$estimate,
          `P value all ages` = primary_row$primary_p,
          `P value adults only` = adult_fit$p_value,
          `Conclusion unchanged` =
            (primary_row$primary_p < 0.05) == (adult_fit$p_value < 0.05)
        )
      })
    }
  )

  stopifnot(
    nrow(figure4_age_sensitivity) == nrow(figure4_primary_numeric),
    all(figure4_age_sensitivity$`Patients adults only` <=
      figure4_age_sensitivity$`Patients all ages`),
    !any(is.na(figure4_age_sensitivity$`P value adults only`))
  )

  write.csv(
    figure4_age_sensitivity,
    file.path(OUT_DIR, "Table_Sensitivity_AdultsOnly_LMM.csv"),
    row.names = FALSE
  )
  print(as.data.frame(figure4_age_sensitivity), row.names = FALSE)
}

# Primary longitudinal presentation: one compact row per outcome/parameter.
# E/e' and LAVi share the fixed parent cohort; TRVmax remains available-case.
figure4_lmm_primary_table <- figure4_lmm_confirmation %>%
  transmute(
    Outcome,
    Parameter,
    `Analysis sample` = "Parameter available-case sample",
    Patients,
    Observations,
    `Interaction beta per year per 1 SD` =
      `Time x parameter beta per year per 1 SD`,
    `95% CI` = sprintf("%.2f to %.2f", `CI low`, `CI high`),
    `P value` = case_when(
      `LMM P value` < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", `LMM P value`)
    ),
    `Primary inference` = if_else(
      `LMM P value` < 0.05,
      "Evidence of trajectory modification",
      "No evidence of trajectory modification"
    )
  )

write.csv(
  figure4_lmm_primary_table,
  paste0(
    file.path(OUT_DIR, ""),
    "Table3_Longitudinal_LMM_Primary.csv"
  ),
  row.names = FALSE
)

figure4_lmm_primary_ft <- flextable(figure4_lmm_primary_table) %>%
  align(j = c("Patients", "Observations", "P value"), align = "center") %>%
  align(
    j = c(
      "Interaction beta per year per 1 SD",
      "95% CI"
    ),
    align = "center"
  ) %>%
  colformat_double(
    j = "Interaction beta per year per 1 SD",
    digits = 2
  ) %>%
  format_pub_table(
    caption = paste0(
      "Table 3. Primary longitudinal linear mixed-effects models. ",
      "The interaction estimates the annual difference in outcome trajectory ",
      "per 1-SD higher baseline parameter; models adjust for age, sex, and BMI ",
      "and include a patient random intercept."
    )
  ) %>%
  autofit()

save_as_docx(
  figure4_lmm_primary_ft,
  path = paste0(
    file.path(OUT_DIR, ""),
    "Table3_Longitudinal_LMM_Primary.docx"
  )
)

# Aggregate visit-density audit. No identifiers are exported.
longitudinal_patient_visit_counts <- longitudinal_continuous_df %>%
  filter(!is.na(VO2_FRIEND2_PP), !is.na(time_yrs)) %>%
  group_by(ID) %>%
  summarise(
    visits = n_distinct(cpx_test_date),
    follow_up_span_yrs = max(time_yrs) - min(time_yrs),
    .groups = "drop"
  )

longitudinal_visit_audit <- tibble(
  Metric = c(
    "Patients",
    "Peak VO2 observations",
    "VE/VCO2 observations",
    "Patients with exactly 2 visits",
    "Patients with exactly 3 visits",
    "Patients with >=4 visits",
    "Median visits per patient",
    "25th percentile visits per patient",
    "75th percentile visits per patient",
    "Median follow-up span",
    "25th percentile follow-up span",
    "75th percentile follow-up span",
    "Maximum follow-up span"
  ),
  Value = c(
    n_distinct(longitudinal_continuous_df$ID),
    sum(!is.na(longitudinal_continuous_df$VO2_FRIEND2_PP)),
    sum(!is.na(longitudinal_continuous_df$VeVco2_slope)),
    sum(longitudinal_patient_visit_counts$visits == 2),
    sum(longitudinal_patient_visit_counts$visits == 3),
    sum(longitudinal_patient_visit_counts$visits >= 4),
    median(longitudinal_patient_visit_counts$visits),
    quantile(longitudinal_patient_visit_counts$visits, 0.25),
    quantile(longitudinal_patient_visit_counts$visits, 0.75),
    median(longitudinal_patient_visit_counts$follow_up_span_yrs),
    quantile(longitudinal_patient_visit_counts$follow_up_span_yrs, 0.25),
    quantile(longitudinal_patient_visit_counts$follow_up_span_yrs, 0.75),
    max(longitudinal_patient_visit_counts$follow_up_span_yrs)
  ),
  Unit = c(
    "patients", "observations", "observations", "patients", "patients",
    "patients", "visits", "visits", "visits", "years", "years", "years",
    "years"
  )
)
write.csv(
  longitudinal_visit_audit,
  paste0(
    file.path(OUT_DIR, ""),
    "TableS_Longitudinal_VisitDistribution.csv"
  ),
  row.names = FALSE
)


write.csv(
  figure4_random_structure_audit,
  paste0(
    file.path(OUT_DIR, ""),
    "TableS_Longitudinal_RandomStructure_Audit.csv"
  ),
  row.names = FALSE
)
write.csv(
  figure4_random_structure_decision,
  paste0(
    file.path(OUT_DIR, ""),
    "TableS_Longitudinal_RandomStructure_Decision.csv"
  ),
  row.names = FALSE
)

# Final primary estimates use one prespecified random structure per outcome.
# This avoids choosing a structure separately for individual parameters while
# honoring consistent evidence of slope heterogeneity within an outcome.
figure4_lmm_primary_table <- figure4_random_structure_audit %>%
  left_join(
    figure4_random_structure_decision %>% select(Outcome, Decision),
    by = "Outcome"
  ) %>%
  mutate(
    `Analysis sample` = "Parameter available-case sample",
    `Interaction beta per year per 1 SD` = if_else(
      Decision == "Random intercept plus independent random slope",
      `RS interaction beta`,
      `RI interaction beta`
    ),
    `CI low` = if_else(
      Decision == "Random intercept plus independent random slope",
      `RS interaction CI low`,
      `RI interaction CI low`
    ),
    `CI high` = if_else(
      Decision == "Random intercept plus independent random slope",
      `RS interaction CI high`,
      `RI interaction CI high`
    ),
    `P numeric` = if_else(
      Decision == "Random intercept plus independent random slope",
      `RS interaction P`,
      `RI interaction P`
    ),
    `95% CI` = sprintf("%.2f to %.2f", `CI low`, `CI high`),
    `P value` = case_when(
      `P numeric` < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", `P numeric`)
    )
  ) %>%
  # One prespecified family: 3 parameters x 2 primary CPET outcomes, matching
  # the cross-sectional family, with inference on the adjusted value (owner
  # decision, card 1.9 of 8_Docs/DECISION_CARDS_2026-09.md).
  mutate(
    `Q numeric` = p.adjust(`P numeric`, method = "BH"),
    `Q value` = case_when(
      `Q numeric` < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", `Q numeric`)
    ),
    `Primary inference` = if_else(
      `Q numeric` < 0.05,
      "Evidence of trajectory modification",
      "No evidence of trajectory modification"
    )
  ) %>%
  transmute(
    Outcome,
    Parameter,
    `Analysis sample`,
    Patients,
    Observations,
    `Random-effects structure` = Decision,
    `Interaction beta per year per 1 SD`,
    `95% CI`,
    `P value`,
    `Q value`,
    `Primary inference`
  )

write.csv(
  figure4_lmm_primary_table,
  file.path(OUT_DIR, "Table3_Longitudinal_LMM_Primary.csv"),
  row.names = FALSE
)

figure4_lmm_primary_ft <- flextable(figure4_lmm_primary_table) %>%
  align(
    j = c(
      "Patients",
      "Observations",
      "Interaction beta per year per 1 SD",
      "95% CI",
      "P value",
      "Q value"
    ),
    align = "center"
  ) %>%
  colformat_double(
    j = "Interaction beta per year per 1 SD",
    digits = 2
  ) %>%
  format_pub_table(
    caption = paste0(
      "Table 3. Primary longitudinal linear mixed-effects models. ",
      "The interaction estimates the annual difference in outcome trajectory ",
      "per 1-SD higher baseline parameter; all models adjust for age, sex, and BMI. ",
      "Q values are Benjamini-Hochberg adjusted within the prespecified family ",
      "of six parameter-outcome comparisons, and inference is on the adjusted value."
    )
  ) %>%
  autofit()

save_as_docx(
  figure4_lmm_primary_ft,
  path = file.path(OUT_DIR, "Table3_Longitudinal_LMM_Primary.docx")
)

# Direct selection sensitivity for the two cohort-defining parameters: compare
# each parameter's own available-case repeated-test sample with its nested
# TRVmax-complete subset using the identical linear mixed-model specification.
# The first sample is NOT the common E/e' + LAVi cohort -- the fitter drops
# missing values per parameter, so E/e' and LAVi have different denominators
# here. The common cohort is analysed separately in the consistency table
# above.
figure4_lmm_cohort_selection <- map_dfr(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome) {
    outcome_label <- ifelse(
      outcome == "VO2_FRIEND2_PP",
      "Peak VO2 (% predicted)",
      "VE/VCO2 slope"
    )
    selected_decision <- figure4_random_structure_decision %>%
      filter(Outcome == outcome_label) %>%
      pull(Decision)
    selected_structure <- ifelse(
      selected_decision == "Random intercept plus independent random slope",
      "intercept_slope",
      "intercept"
    )
    map_dfr(
      figure4_index_specs[vapply(
        figure4_index_specs,
        function(spec) spec$var %in% c("e_e_ave", "la_vol_index"),
        logical(1)
      )],
      function(spec) {
        sample_specs <- list(
          list(
            label = "Parameter available-case sample",
            data = longitudinal_continuous_df
          ),
          list(
            label = "TRVmax-complete longitudinal subset",
            data = longitudinal_continuous_df %>%
              filter(!is.na(baseline_tr_max_vel))
          )
        )
        map_dfr(sample_specs, function(sample_spec) {
          fit_obj <- fit_linear_mixed_confirmation(
            sample_spec$data,
            outcome,
            spec,
            random_structure = selected_structure
          )
          if (is.null(fit_obj)) return(tibble())
          tibble(
            Outcome = outcome_label,
            Parameter = spec$label_short,
            `Analysis sample` = sample_spec$label,
            `Random-effects structure` = selected_decision,
            Patients = n_distinct(fit_obj$data$ID),
            Observations = nrow(fit_obj$data),
            `Time x parameter beta per year per 1 SD` = fit_obj$estimate,
            `CI low` = fit_obj$ci_low,
            `CI high` = fit_obj$ci_high,
            `P value` = fit_obj$p_value
          )
        })
      }
    )
  }
)

write.csv(
  figure4_lmm_cohort_selection,
  paste0(
    file.path(AUDIT_DIR, ""),
    "Table_CohortSelection_Longitudinal_EeLAVi.csv"
  ),
  row.names = FALSE
)


# ══ Supplemental Figures: Parameter-Specific Binary Trajectories ══════════════
# (legacy chunk: figureS_diastolic_continuous)
# ── Fig S: GAMM trajectories split by each binary index (supplement to Fig 4) ──
for (spec in figure4_index_specs) {
  result <- build_binary_subfigure(spec)
  if (is.null(result)) {
    next
  }

  pdf_name <- paste0("FigureS_", spec$suffix, "_Binary_Trajectories")
  output_pdf <- paste0(file.path(OUT_DIR, ""), pdf_name, ".pdf")
  manuscript_pdf <- paste0(file.path(OUT_DIR, ""), pdf_name, ".pdf")

  show_and_save_jacc(result$fig, output_pdf, w = 7.2, h = 3.1)

  legend_text <- paste0(
    "**Supplemental Figure.** Parameter-specific longitudinal cardiopulmonary trajectories stratified by binary baseline ",
    spec$label_long,
    " status. ",
    "**(A)** GAMM-derived ",
    label_peak_vo2_md,
    " trajectories. ",
    "**(B)** GAMM-derived ",
    label_vevco2_md,
    " trajectories. ",
    "**(C)** Summary of the key status and interaction terms for both models. ",
    "Baseline threshold groups were defined as ",
    spec$threshold_subtitle,
    ". ",
    "Longitudinal models were adjusted for age, sex, and BMI."
  )

  save_jacc_with_embedded_legend(
    result$fig,
    manuscript_pdf,
    legend_text = legend_text,
    w = 7.2,
    h = 3.1
  )

  summary_name <- paste0("TableS_", spec$suffix, "_Binary_GAMM_Summary.csv")
  write.csv(
    result$summary,
    paste0(file.path(OUT_DIR, ""), summary_name),
    row.names = FALSE
  )

  # (Legacy copies of these outputs under *_Tertile_* filenames were aliases of
  #  the binary-threshold figures, not tertile analyses, and are not produced.)
}

message("Stage 06 complete.")
