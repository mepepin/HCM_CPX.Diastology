# ── Stage 07: clinical outcomes ───────────────────────────────────────────────
# In  : 5_Data/2_Processed/{baseline_df, long_df, outcomes*, medications,
#                           z_source, cohort_parameters}.rds
# Out : 6_Results/5_Outcomes/
#         Figure5_HeartFailure_Outcomes.pdf              (main Figure 5)
#         Figure5_HeartFailure_Outcomes_Sensitivity.pdf  (Figure S4)
#         Table_Figure5_Event_Components.csv             (main Table 2)
#         TableS_Figure5_Cox_Summary.csv/.docx           (Table S3)
#         TableS_Cox_PH_Diagnostics.csv                  (Table S4)
#         Table_Consistency_MeasurementCohorts_Cox.csv   (Table S5)
#         Table_Manuscript_Outcome_Values.csv            (in-text values)
#         Table_Figure5_Myectomy_Cox.csv                 (secondary end point)
#         Table_Sensitivity_*_Cox.csv                    (Methods sensitivities)
#
# Each exposure is fitted in its own Cox model per SD, adjusted for age, sex,
# and BMI; peak VO2, VE/VCO2 slope, septal thickness and the resting LVOT
# gradient are comparators fitted in parallel models, never jointly.

source("R/00_setup.R")
stage_banner("Stage 07: clinical outcomes")

OUT_DIR <- results_dir("5_Outcomes")
AUDIT_DIR <- results_dir("1_Audit")


# ── Inputs from earlier stages ────────────────────────────────────────────────
baseline_df <- load_data("baseline_df", "processed")
long_df <- load_data("long_df", "processed")
outcomes <- load_data("outcomes", "processed")
outcomes_dedup <- load_data("outcomes_dedup", "processed")
outcomes_index <- load_data("outcomes_index", "processed")
medications <- load_data("medications", "processed")
z_source <- load_data("z_source", "processed")
outcomes_raw <- load_data("outcomes_raw", "interim")
cohort_parameters <- load_data("cohort_parameters", "processed")
analytic_cohort_mode <- cohort_parameters$analytic_cohort_mode

# ══ Kaplan-Meier Curves and Cox Models ════════════════════════════════════════
# (legacy chunk: figure5_survival)
# ── Fig 5: Cox PH for composite HF endpoint ───────────────────────────────────
# time = test -> event/censor (yrs); censor at last encounter. predictors scaled
# to interpretable units (per 10 mmHg, per 0.5cm, per SD z) so HRs are comparable.
# also derives a transplant-free endpoint + intervention timing for sensitivity
surv_df <- baseline_df %>%
  # Every base-eligible patient with adequate follow-up: an endpoint recorded
  # and (for events) a usable event time. Each Kaplan-Meier stratification and
  # each Cox model then reduces to its own exposure's available cases, so no
  # exposure loses patients because a different parameter was unmeasured
  # (8_Docs/ANALYSIS_PLAN_2026-09.md). No imputation anywhere.
  filter(analysis_base_eligible, outcomes_followup_ok) %>%
  mutate(
    hf_composite_date = as.Date(hf_composite_date_num, origin = "1970-01-01"),
    hf_or_death_no_transplant_date = as.Date(
      hf_or_death_no_transplant_date_num,
      origin = "1970-01-01"
    ),
    major_intervention_date_num = pmin(
      coalesce(as.numeric(post_septal_reduction_surgery_date), Inf),
      coalesce(as.numeric(post_ablation_surgery_date), Inf),
      coalesce(as.numeric(post_heart_transplant_date), Inf),
      na.rm = TRUE
    ),
    major_intervention_date_num = ifelse(
      is.infinite(major_intervention_date_num),
      NA_real_,
      major_intervention_date_num
    ),
    major_intervention_date = as.Date(
      major_intervention_date_num,
      origin = "1970-01-01"
    ),
    event_time_yrs = as.numeric(hf_composite_date - cpx_test_date) / 365.25,
    event_time_no_transplant_yrs = as.numeric(
      hf_or_death_no_transplant_date - cpx_test_date
    ) /
      365.25,
    myectomy_event_time_yrs = (myectomy_date_num - as.numeric(cpx_test_date)) /
      365.25,
    myectomy_event = as.integer(
      !is.na(myectomy_event_time_yrs) & myectomy_event_time_yrs > 0
    ),
    major_intervention_time_yrs = as.numeric(
      major_intervention_date - cpx_test_date
    ) /
      365.25,
    censor_time_yrs = as.numeric(as.Date(last_enc_date) - cpx_test_date) /
      365.25,
    lvot_max_gradient_10 = lvot_max_gradient / 10,
    lv_septal_thickness_0_5 = lv_septal_thickness / 0.5,
    vo2_friend2_10lower = -VO2_FRIEND2_PP / 10,
    vevco2_slope_5 = VeVco2_slope / 5,
    e_e_ave_z = as.numeric(scale(e_e_ave)),
    la_vol_index_z = as.numeric(scale(la_vol_index)),
    tr_max_vel_z = as.numeric(scale(tr_max_vel)),
    lvot_gradient_z = as.numeric(scale(lvot_max_gradient)),
    septal_thickness_z = as.numeric(scale(lv_septal_thickness)),
    vo2_friend2_z = as.numeric(scale(-VO2_FRIEND2_PP)),
    vevco2_slope_z = as.numeric(scale(VeVco2_slope)),
    follow_up_yrs = case_when(
      hf_composite == 1 &
        !is.na(event_time_yrs) &
        event_time_yrs > 0 ~ event_time_yrs,
      !is.na(censor_time_yrs) & censor_time_yrs > 0 ~ censor_time_yrs,
      TRUE ~ NA_real_
    ),
    follow_up_yrs = ifelse(
      is.na(follow_up_yrs) | follow_up_yrs <= 0,
      0.01,
      follow_up_yrs
    ),
    myectomy_follow_up_yrs = case_when(
      myectomy_event == 1 ~ myectomy_event_time_yrs,
      !is.na(censor_time_yrs) & censor_time_yrs > 0 ~ censor_time_yrs,
      TRUE ~ NA_real_
    )
  ) %>%
  filter(!is.na(follow_up_yrs), follow_up_yrs > 0)

# Outcomes cohort follow-up summary (for manuscript comment #19)
# surv_df already has follow_up_yrs computed and filtered; use it directly
outcomes_fu_summary <- surv_df %>%
  summarise(
    n_total = n(),
    n_events = sum(hf_composite == 1, na.rm = TRUE),
    n_censored = sum(hf_composite == 0, na.rm = TRUE),
    med_fu = median(follow_up_yrs, na.rm = TRUE),
    q1_fu = quantile(follow_up_yrs, 0.25, na.rm = TRUE),
    q3_fu = quantile(follow_up_yrs, 0.75, na.rm = TRUE)
  )
cat(sprintf(
  "Outcomes cohort: N=%d, events=%d (%.1f%%), censored=%d\n",
  outcomes_fu_summary$n_total,
  outcomes_fu_summary$n_events,
  100 * outcomes_fu_summary$n_events / outcomes_fu_summary$n_total,
  outcomes_fu_summary$n_censored
))
cat(sprintf(
  "Median follow-up: %.1f years (IQR %.1f–%.1f)\n",
  outcomes_fu_summary$med_fu,
  outcomes_fu_summary$q1_fu,
  outcomes_fu_summary$q3_fu
))

# Detailed component accounting. "Any occurrence" is non-mutually exclusive;
# "first composite event" assigns each patient to the event that determined the
# primary time-to-event endpoint and explicitly preserves same-day ties.
surv_df <- surv_df %>%
  rowwise() %>%
  mutate(
    first_event_date_num = min(
      c(acute_hf_date_num, transplant_date_num, death_date_num),
      na.rm = TRUE
    ),
    first_event_date_num = ifelse(
      is.infinite(first_event_date_num),
      NA_real_,
      first_event_date_num
    ),
    first_event_ties = sum(
      c(acute_hf_date_num, transplant_date_num, death_date_num) ==
        first_event_date_num,
      na.rm = TRUE
    ),
    first_event_component = case_when(
      is.na(first_event_date_num) ~ "No composite event",
      first_event_ties > 1 ~ "Multiple components on same day",
      acute_hf_date_num == first_event_date_num ~ "Acute heart failure",
      transplant_date_num == first_event_date_num ~ "Heart transplantation",
      death_date_num == first_event_date_num ~ "All-cause death",
      TRUE ~ "Unclassified"
    )
  ) %>%
  ungroup()

event_component_table <- tibble(
  Event = c(
    "Primary composite",
    "Acute heart failure",
    "Heart transplantation",
    "All-cause death",
    "Multiple components on same day",
    "Chronic heart-failure diagnosis (not in composite)",
    "Surgical myectomy (separate secondary outcome)"
  ),
  any_n = c(
    sum(surv_df$hf_composite == 1, na.rm = TRUE),
    sum(!is.na(surv_df$acute_hf_date_num)),
    sum(!is.na(surv_df$transplant_date_num)),
    sum(!is.na(surv_df$death_date_num)),
    sum(surv_df$first_event_component == "Multiple components on same day"),
    sum(!is.na(surv_df$chronic_hf_date_num)),
    sum(surv_df$myectomy_event == 1, na.rm = TRUE)
  ),
  first_n = c(
    sum(surv_df$hf_composite == 1, na.rm = TRUE),
    sum(surv_df$first_event_component == "Acute heart failure"),
    sum(surv_df$first_event_component == "Heart transplantation"),
    sum(surv_df$first_event_component == "All-cause death"),
    sum(surv_df$first_event_component == "Multiple components on same day"),
    NA_integer_,
    NA_integer_
  )
) %>%
  mutate(
    `Any occurrence, n (%)` = sprintf(
      "%d (%.1f%%)",
      any_n,
      100 * any_n / nrow(surv_df)
    ),
    `First composite event, n (%)` = ifelse(
      is.na(first_n),
      "—",
      sprintf("%d (%.1f%%)", first_n, 100 * first_n / nrow(surv_df))
    )
  ) %>%
  select(Event, `Any occurrence, n (%)`, `First composite event, n (%)`)

write.csv(
  event_component_table,
  file.path(OUT_DIR, "Table_Figure5_Event_Components.csv"),
  row.names = FALSE
)
event_component_table %>%
  flextable() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  autofit() %>%
  save_as_docx(
    path = file.path(OUT_DIR, "Table_Figure5_Event_Components.docx")
  )

# Source-semantics audit. This is intentionally aggregate and records what the
# supplied workbook can and cannot establish. A dated acute-HF variable is
# present, but no hospitalization/admission field or data dictionary is present;
# therefore the defensible label is "acute heart failure event."
acute_hf_source_qc <- outcomes %>%
  group_by(MRN) %>%
  summarise(
    acute_hf_flag = max_numeric_or_na(post_acute_heart_failure),
    acute_hf_date = min_date_or_na(post_acute_heart_failure_date),
    .groups = "drop"
  )

endpoint_definition_audit <- tibble(
  Audit_item = c(
    "Dated acute heart failure field",
    "Hospitalization/admission field in Outcomes_2025",
    "Full source extract: unique patients with acute-HF flag = 1",
    "Full source extract: unique patients with an acute-HF date",
    "Full source extract: flag = 1 but date unavailable",
    "Full source extract: date present but flag not equal to 1",
    "Defensible manuscript label",
    "Hospitalization wording"
  ),
  Result = c(
    "post_acute_heart_failure_date",
    ifelse(
      any(grepl("hospital|admission", names(outcomes_raw), ignore.case = TRUE)),
      "Present",
      "Absent"
    ),
    as.character(sum(acute_hf_source_qc$acute_hf_flag == 1, na.rm = TRUE)),
    as.character(sum(!is.na(acute_hf_source_qc$acute_hf_date))),
    as.character(sum(
      acute_hf_source_qc$acute_hf_flag == 1 &
        is.na(acute_hf_source_qc$acute_hf_date),
      na.rm = TRUE
    )),
    as.character(sum(
      !is.na(acute_hf_source_qc$acute_hf_date) &
        (is.na(acute_hf_source_qc$acute_hf_flag) |
          acute_hf_source_qc$acute_hf_flag != 1),
      na.rm = TRUE
    )),
    "Hospitalization for acute heart failure",
    # The extract carries flags and dates, not codes or encounter types. The
    # owner attested on 2026-09-18 that events were ascertained from ICD-coded
    # EHR diagnoses and the dates of hospitalization and death, which supplies
    # the data dictionary this row previously said was missing. The code lists
    # themselves still have to come from the extraction query.
    "Established by owner attestation (2026-09-18): ICD-coded EHR hospitalization and death dates; code lists pending from the extraction query"
  )
)

write.csv(
  endpoint_definition_audit,
  file.path(OUT_DIR, "Table_Endpoint_Definition_Audit.csv"),
  row.names = FALSE
)

fig5_palette <- c(
  reference = "#C4CBD3",
  highlight = "#C88B8B",
  forest = "#2F4858",
  text = "#1F2328",
  refline = "#B0B8C1"
)

make_km_panel <- function(
  data,
  strata_var,
  title,
  labels = c("Normal", "Abnormal"),
  inplot_legend = FALSE
) {
  surv_sub <- data %>%
    filter(!is.na(.data[[strata_var]])) %>%
    mutate(
      km_group = factor(
        ifelse(.data[[strata_var]], labels[2], labels[1]),
        levels = labels
      )
    )

  if (nrow(surv_sub) < 20 || n_distinct(surv_sub$km_group) < 2) {
    return(NULL)
  }

  km_fit <- survfit(
    Surv(follow_up_yrs, hf_composite) ~ km_group,
    data = surv_sub
  )
  km_data <- broom::tidy(km_fit) %>%
    mutate(
      strata = str_remove(strata, "^km_group="),
      strata = factor(strata, levels = labels)
    ) %>%
    filter(!is.na(strata))

  lr_test <- survdiff(
    Surv(follow_up_yrs, hf_composite) ~ km_group,
    data = surv_sub
  )
  lr_p <- 1 - pchisq(lr_test$chisq, df = length(lr_test$n) - 1)
  km_colors <- setNames(
    c(fig5_palette["reference"], fig5_palette["highlight"]),
    labels
  )

  # KM curve
  p_km <- ggplot(
    km_data,
    aes(x = time, y = estimate, color = strata, fill = strata, group = strata)
  ) +
    geom_step(linewidth = 0.9) +
    geom_ribbon(
      aes(ymin = conf.low, ymax = conf.high),
      alpha = 0.15,
      stat = "identity",
      show.legend = FALSE
    ) +
    scale_color_manual(values = km_colors, name = NULL) +
    scale_fill_manual(values = km_colors, guide = "none") +
    scale_x_continuous(
      limits = c(0, 10),
      breaks = seq(0, 10, by = 2),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_continuous(limits = c(0, 1), labels = percent) +
    labs(
      title = title,
      subtitle = paste0("Log-rank p = ", format.pval(lr_p, digits = 2)),
      x = "Years from Baseline CPET",
      y = "Event-Free Survival"
    ) +
    theme_jacc(base_size = 8.5) +
    theme(
      legend.position = if (inplot_legend) c(0.045, 0.06) else "none",
      legend.justification = c(0, 0),
      legend.direction = "vertical",
      legend.background = element_rect(
        fill = alpha("white", 0.9),
        color = "grey75",
        linewidth = 0.25
      ),
      legend.key = element_rect(fill = alpha("white", 0), color = NA),
      legend.margin = margin(3, 5, 3, 5),
      legend.key.width = unit(0.55, "cm"),
      legend.key.height = unit(0.28, "cm"),
      legend.text = element_text(size = 6.8),
      plot.title = element_text(size = 8.8),
      plot.subtitle = element_text(size = 7.4),
      axis.line = element_line(color = fig5_palette["text"], linewidth = 0.35),
      axis.text = element_text(color = fig5_palette["text"]),
      axis.title.y = element_text(
        color = fig5_palette["text"],
        margin = margin(r = 3)
      ),
      axis.title.x = element_text(size = 7.0, color = fig5_palette["text"]),
      axis.text.x = element_text(size = 6.5, color = fig5_palette["text"]),
      axis.ticks.x = element_line(
        color = fig5_palette["text"],
        linewidth = 0.25
      ),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    guides(
      color = guide_legend(
        override.aes = list(linewidth = 0.95, alpha = 1),
        keywidth = unit(0.55, "cm"),
        keyheight = unit(0.28, "cm")
      )
    )

  p_km
}

p_km_ee <- make_km_panel(
  surv_df,
  "abn_ee",
  "E/e' > 14",
  labels = c("E/e' ≤ 14", "E/e' > 14"),
  inplot_legend = TRUE
)
p_km_lavi <- make_km_panel(
  surv_df,
  "abn_lavi",
  "LAVi > 34 mL/m²",
  labels = c("LAVi ≤ 34", "LAVi > 34"),
  inplot_legend = TRUE
)
p_km_vevco2 <- make_km_panel(
  surv_df %>%
    mutate(
      high_vevco2_30 = case_when(
        is.na(VeVco2_slope) ~ NA,
        VeVco2_slope > 30 ~ TRUE,
        TRUE ~ FALSE
      )
    ),
  "high_vevco2_30",
  "V̇E/V̇CO₂ Slope > 30",
  labels = c("V̇E/V̇CO₂ ≤ 30", "V̇E/V̇CO₂ > 30"),
  inplot_legend = TRUE
)
p_km_vo2 <- make_km_panel(
  surv_df %>%
    mutate(
      low_peak_vo2_80 = case_when(
        is.na(VO2_FRIEND2_PP) ~ NA,
        VO2_FRIEND2_PP < 80 ~ TRUE,
        TRUE ~ FALSE
      )
    ),
  "low_peak_vo2_80",
  "Peak V̇O₂ < 80% Predicted",
  labels = c("Peak V̇O₂ ≥ 80%", "Peak V̇O₂ < 80%"),
  inplot_legend = TRUE
)

surv_model_df <- surv_df %>%
  mutate(
    pathogenic_variant_model = if_else(pathogenic_variant == 1, 1L, 0L)
  ) %>%
  filter(
    !is.na(age),
    !is.na(Sex),
    !is.na(BMI)
  )

cox_z_specs <- list(
  ee = list(var = "e_e_ave_z", label = "EEprime"),
  lavi = list(var = "la_vol_index_z", label = "LAVi"),
  trv = list(var = "tr_max_vel_z", label = "TRVmax"),
  lvot = list(var = "lvot_gradient_z", label = "LVOTgradient"),
  septal = list(var = "septal_thickness_z", label = "SeptalThickness"),
  vo2 = list(var = "vo2_friend2_z", label = "PeakVO2"),
  vevco2 = list(var = "vevco2_slope_z", label = "VEVCO2")
)

cox_z_results <- imap(cox_z_specs, function(spec, nm) {
  df <- surv_model_df %>% filter(!is.na(.data[[spec$var]]))
  model <- if (nrow(df) >= 20 && sum(df$hf_composite == 1, na.rm = TRUE) >= 5) {
    coxph(
      as.formula(paste0(
        "Surv(follow_up_yrs, hf_composite) ~ ",
        spec$var,
        " + age + Sex + BMI"
      )),
      data = df
    )
  } else {
    NULL
  }
  if (!is.null(model)) {
    write.csv(
      tidy(model, exponentiate = TRUE, conf.int = TRUE),
      paste0(file.path(OUT_DIR, "Table3_Cox_"), spec$label, ".csv"),
      row.names = FALSE
    )
  }
  list(model = model, data = df)
})

surv_ee <- cox_z_results$ee$data
cox_ee <- cox_z_results$ee$model
surv_lavi <- cox_z_results$lavi$data
cox_lavi <- cox_z_results$lavi$model
surv_trv <- cox_z_results$trv$data
cox_trv <- cox_z_results$trv$model
surv_lvot <- cox_z_results$lvot$data
cox_lvot <- cox_z_results$lvot$model
surv_septal <- cox_z_results$septal$data
cox_septal <- cox_z_results$septal$model
surv_vo2 <- cox_z_results$vo2$data
cox_vo2 <- cox_z_results$vo2$model
surv_vevco2 <- cox_z_results$vevco2$data
cox_vevco2 <- cox_z_results$vevco2$model

# Outcome-model denominator ladder; this is reporting only and does not alter
# the prespecified primary age/sex/BMI models above.
cox_denominator_ladder <- tidyr::crossing(
  Parameter = c("E/e'", "LAVi", "TRVmax"),
  Tier = c(
    "Primary: parameter + age/sex/BMI",
    "Sensitivity: + apical morphology",
    "Sensitivity: + validated LVOT gradient",
    "Sensitivity: expanded structure/medications/hypertension"
  )
) %>%
  mutate(
    parameter_var = recode(Parameter,
      "E/e'" = "e_e_ave_z",
      "LAVi" = "la_vol_index_z",
      "TRVmax" = "tr_max_vel_z"
    ),
    required_vars = map2(parameter_var, Tier, function(parameter, tier) {
      vars <- c(parameter, "follow_up_yrs", "hf_composite", "age", "Sex", "BMI")
      if (grepl("apical", tier)) vars <- c(vars, "apical_hcm")
      if (grepl("validated LVOT", tier)) vars <- c(vars, "lvot_max_gradient")
      if (grepl("expanded", tier)) {
        vars <- c(
          vars,
          "lv_septal_thickness",
          "lvot_max_gradient",
          "bb_any",
          "ndhp_ccb_any",
          "htn_pre_test"
        )
      }
      unique(vars)
    }),
    N = map_int(required_vars, function(vars) {
      sum(complete.cases(surv_df[, vars, drop = FALSE]))
    }),
    Events = map2_int(required_vars, N, function(vars, n) {
      cc <- complete.cases(surv_df[, vars, drop = FALSE])
      sum(surv_df$hf_composite[cc] == 1, na.rm = TRUE)
    })
  ) %>%
  select(Parameter, Tier, N, Events)

write.csv(
  cox_denominator_ladder,
  file.path(OUT_DIR, "Table_Cox_Denominator_Ladder.csv"),
  row.names = FALSE
)

# Deprecated composite classifier is intentionally excluded from the primary
# outcome models; retain null placeholders for legacy supplemental code paths.
cox_fp <- NULL
surv_fp_model <- tibble()

# Crude (unadjusted) Cox models for the three non-significant predictors in Figure 5E,
# used in the supplemental crude-vs-adjusted comparison figure.
cox_lvot_crude <- if (!is.null(cox_lvot)) {
  tryCatch(
    coxph(
      Surv(follow_up_yrs, hf_composite) ~ lvot_gradient_z,
      data = surv_lvot
    ),
    error = function(e) NULL
  )
} else {
  NULL
}

cox_vo2_crude <- if (!is.null(cox_vo2)) {
  tryCatch(
    coxph(Surv(follow_up_yrs, hf_composite) ~ vo2_friend2_z, data = surv_vo2),
    error = function(e) NULL
  )
} else {
  NULL
}

cox_vevco2_crude <- if (!is.null(cox_vevco2)) {
  tryCatch(
    coxph(
      Surv(follow_up_yrs, hf_composite) ~ vevco2_slope_z,
      data = surv_vevco2
    ),
    error = function(e) NULL
  )
} else {
  NULL
}

# Secondary septal-reduction-therapy outcome. Only surgical myectomy has a
# validated procedure date in this extract; the generic ablation field is not
# used because it cannot distinguish alcohol septal ablation from arrhythmia
# ablation. Models are parsimoniously adjusted for validated baseline LVOT
# gradient, the principal clinical determinant of myectomy referral.
myectomy_specs <- list(
  ee = list(var = "e_e_ave_z", label = "E/e'"),
  lavi = list(var = "la_vol_index_z", label = "LAVi"),
  trv = list(var = "tr_max_vel_z", label = "TRVmax")
)

myectomy_cox_results <- imap(myectomy_specs, function(spec, nm) {
  df <- surv_df %>%
    filter(
      !is.na(myectomy_follow_up_yrs),
      myectomy_follow_up_yrs > 0,
      !is.na(.data[[spec$var]]),
      !is.na(lvot_gradient_z)
    )
  model <- if (
    nrow(df) >= 20 && sum(df$myectomy_event == 1, na.rm = TRUE) >= 10
  ) {
    coxph(
      as.formula(paste0(
        "Surv(myectomy_follow_up_yrs, myectomy_event) ~ ",
        spec$var,
        " + lvot_gradient_z"
      )),
      data = df
    )
  } else {
    NULL
  }
  list(model = model, data = df, spec = spec)
})

myectomy_forest_rows <- imap_dfr(myectomy_cox_results, function(x, nm) {
  if (is.null(x$model)) {
    return(tibble())
  }
  tidy(x$model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == x$spec$var) %>%
    transmute(
      exposure = x$spec$label,
      N = nrow(x$data),
      Events = sum(x$data$myectomy_event == 1, na.rm = TRUE),
      estimate,
      conf.low,
      conf.high,
      p.value
    )
})

write.csv(
  myectomy_forest_rows,
  file.path(OUT_DIR, "Table_Figure5_Myectomy_Cox.csv"),
  row.names = FALSE
)

# HCM morphology sensitivity for the primary clinical composite. Preserve the
# primary age/sex/BMI adjustment and add apical-vs-known-nonapical morphology;
# unknown morphology is not recoded as nonapical.
morphology_cox_rows <- imap_dfr(
  myectomy_specs,
  function(spec, nm) {
    df <- surv_df %>%
      filter(
        !is.na(.data[[spec$var]]),
        !is.na(age),
        !is.na(Sex),
        !is.na(BMI),
        !is.na(apical_hcm)
      )
    model <- if (
      nrow(df) >= 20 && sum(df$hf_composite == 1, na.rm = TRUE) >= 10
    ) {
      coxph(
        as.formula(paste0(
          "Surv(follow_up_yrs, hf_composite) ~ ",
          spec$var,
          " + age + Sex + BMI + apical_hcm"
        )),
        data = df
      )
    } else {
      NULL
    }
    if (is.null(model)) {
      return(tibble())
    }
    tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term %in% c(spec$var, "apical_hcm")) %>%
      transmute(
        Parameter = spec$label,
        Term = ifelse(term == "apical_hcm", "Apical vs nonapical HCM", spec$label),
        N = nrow(df),
        Events = sum(df$hf_composite == 1, na.rm = TRUE),
        HR = estimate,
        CI_low = conf.low,
        CI_high = conf.high,
        p_value = p.value
      )
  }
)

write.csv(
  morphology_cox_rows,
  file.path(OUT_DIR, "Table_Sensitivity_Cox_Morphology.csv"),
  row.names = FALSE
)

# Same-sample robustness for E/e' and LAVi, plus available-case TRVmax. Every
# unadjusted/adjusted comparison for a given parameter uses identical patients
# and events, so a change in the estimate cannot be attributed to denominator
# loss. E/e' and LAVi share the fixed primary cohort by design.
primary_cox_exposures <- c(
  e_e_ave_z = "E/e'",
  la_vol_index_z = "LAVi",
  tr_max_vel_z = "TRVmax"
)

fit_same_sample_cox <- function(data, exposure, covariates, model_label) {
  formula <- reformulate(
    c(exposure, covariates),
    response = "Surv(follow_up_yrs, hf_composite)"
  )
  fit <- tryCatch(
    coxph(formula, data = data, ties = "efron", x = TRUE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  exposure_row <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == exposure)
  if (nrow(exposure_row) != 1) return(NULL)

  ph_global_p <- tryCatch(
    unname(cox.zph(fit)$table["GLOBAL", "p"]),
    error = function(e) NA_real_
  )
  morphology_row <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "apical_hcm")

  exposure_row %>%
    transmute(
      Model = model_label,
      N = nrow(data),
      Events = sum(data$hf_composite == 1, na.rm = TRUE),
      `Model coefficients` = length(coef(fit)),
      `Events per coefficient` = Events / `Model coefficients`,
      HR = estimate,
      `CI low` = conf.low,
      `CI high` = conf.high,
      `P value` = p.value,
      `Global PH P value` = ph_global_p,
      `Apical morphology HR` = ifelse(
        nrow(morphology_row) == 1,
        morphology_row$estimate,
        NA_real_
      ),
      `Apical morphology P value` = ifelse(
        nrow(morphology_row) == 1,
        morphology_row$p.value,
        NA_real_
      )
    )
}

primary_cox_common_required <- c(
  "follow_up_yrs",
  "hf_composite",
  "e_e_ave_z",
  "la_vol_index_z",
  "age",
  "Sex",
  "BMI"
)

primary_cox_same_sample_rows <- map_dfr(
  names(primary_cox_exposures),
  function(exposure) {
    parameter_df <- surv_df %>%
      filter(complete.cases(across(all_of(c(
        primary_cox_common_required,
        exposure
      )))))
    bind_rows(
      fit_same_sample_cox(
        parameter_df,
        exposure,
        character(0),
        "Unadjusted"
      ),
      fit_same_sample_cox(
        parameter_df,
        exposure,
        c("age", "Sex", "BMI"),
        "Age/sex/BMI adjusted"
      )
    ) %>%
      mutate(
        Parameter = unname(primary_cox_exposures[[exposure]]),
        `Analysis sample` = ifelse(
          exposure == "tr_max_vel_z" && analytic_cohort_mode == "ee_lavi",
          "TRVmax available-case outcomes subset",
          "Primary E/e' + LAVi outcomes cohort"
        ),
        .before = 1
      )
  }
)

known_morphology_same_sample_rows <- map_dfr(
  names(primary_cox_exposures),
  function(exposure) {
    known_morphology_cox_df <- surv_df %>%
      filter(complete.cases(across(all_of(c(
        primary_cox_common_required,
        exposure,
        "apical_hcm"
      )))))
    bind_rows(
      fit_same_sample_cox(
        known_morphology_cox_df,
        exposure,
        c("age", "Sex", "BMI"),
        "Age/sex/BMI adjusted"
      ),
      fit_same_sample_cox(
        known_morphology_cox_df,
        exposure,
        c("age", "Sex", "BMI", "apical_hcm"),
        "Age/sex/BMI + apical morphology"
      )
    ) %>%
      mutate(
        Parameter = unname(primary_cox_exposures[[exposure]]),
        `Analysis sample` = "Known-morphology outcomes subset",
        .before = 1
      )
  }
)

cox_same_sample_sensitivity <- bind_rows(
  primary_cox_same_sample_rows,
  known_morphology_same_sample_rows
)

cox_same_sample_n_audit <- cox_same_sample_sensitivity %>%
  group_by(`Analysis sample`, Parameter) %>%
  summarise(
    patient_counts = n_distinct(N),
    event_counts = n_distinct(Events),
    .groups = "drop"
  )
stopifnot(
  all(cox_same_sample_n_audit$patient_counts == 1),
  all(cox_same_sample_n_audit$event_counts == 1)
)

write.csv(
  cox_same_sample_sensitivity,
  paste0(
    file.path(OUT_DIR, ""),
    "Table_Sensitivity_Cox_SameSample.csv"
  ),
  row.names = FALSE
)

# ── Manuscript values: observed follow-up on the survival frame ──────────────
# The Results report the observed time to first event or censoring. It was not
# exported anywhere, so the sentence had no traceable source.
manuscript_outcome_values <- tibble(
  quantity = c(
    "Patients in the survival frame, n",
    "Composite events, n (%)",
    "Observed time to first event or censoring, median (IQR), years"
  ),
  value = c(
    as.character(nrow(surv_df)),
    sprintf(
      "%d (%.1f%%)",
      sum(surv_df$hf_composite == 1, na.rm = TRUE),
      100 * mean(surv_df$hf_composite == 1, na.rm = TRUE)
    ),
    sprintf(
      "%.1f (IQR %.1f to %.1f)",
      median(surv_df$follow_up_yrs, na.rm = TRUE),
      quantile(surv_df$follow_up_yrs, 0.25, na.rm = TRUE),
      quantile(surv_df$follow_up_yrs, 0.75, na.rm = TRUE)
    )
  )
)
stopifnot(!any(is.na(manuscript_outcome_values$value)))
write.csv(
  manuscript_outcome_values,
  file.path(OUT_DIR, "Table_Manuscript_Outcome_Values.csv"),
  row.names = FALSE
)
print(as.data.frame(manuscript_outcome_values), row.names = FALSE)

# ── Disclosure: measurement completeness by outcome status ────────────────────
# Available-case analysis assumes that whether a parameter was measured is
# unrelated to the outcome given covariates. It is not, here, so the pattern is
# published rather than assumed away (owner decision, card 1.2 of
# 8_Docs/DECISION_CARDS_2026-09.md). "Measured" uses the same definition as the
# Cox models -- the parameter and the adjustment covariates all present -- so
# the totals in this table are exactly the Cox denominators, which is asserted
# below.
measurement_disclosure_params <- c(
  e_e_ave = "E/e'",
  la_vol_index = "LAVi",
  tr_max_vel = "TRVmax"
)

measurement_by_outcome <- map_dfr(
  names(measurement_disclosure_params),
  function(param) {
    d <- surv_df %>%
      filter(!is.na(hf_composite)) %>%
      mutate(
        measured = as.integer(
          !is.na(.data[[param]]) & !is.na(age) & !is.na(Sex) & !is.na(BMI)
        ),
        event = as.integer(hf_composite == 1)
      )
    crude <- glm(measured ~ event, data = d, family = binomial())
    adjusted <- glm(
      measured ~ event + age + Sex + BMI,
      data = d,
      family = binomial()
    )
    crude_row <- broom::tidy(crude, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "event")
    adj_row <- broom::tidy(adjusted, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "event")
    n_event <- sum(d$event == 1)
    n_no_event <- sum(d$event == 0)
    measured_event <- sum(d$measured == 1 & d$event == 1)
    measured_no_event <- sum(d$measured == 1 & d$event == 0)
    tibble(
      Parameter = unname(measurement_disclosure_params[[param]]),
      `Patients with follow-up` = nrow(d),
      Events = n_event,
      `Measured among events` = measured_event,
      `% measured among events` = 100 * measured_event / n_event,
      `Non-events` = n_no_event,
      `Measured among non-events` = measured_no_event,
      `% measured among non-events` = 100 * measured_no_event / n_no_event,
      `Difference (percentage points)` =
        100 * measured_event / n_event - 100 * measured_no_event / n_no_event,
      `Crude OR measured given event` = crude_row$estimate,
      `Crude CI low` = crude_row$conf.low,
      `Crude CI high` = crude_row$conf.high,
      `Crude P value` = crude_row$p.value,
      `Adjusted OR measured given event` = adj_row$estimate,
      `Adjusted CI low` = adj_row$conf.low,
      `Adjusted CI high` = adj_row$conf.high,
      `Adjusted P value` = adj_row$p.value,
      `Total measured` = measured_event + measured_no_event
    )
  }
)

# Tie the disclosure table to the models it describes: each parameter's total
# measured count must equal the N its own Cox model used. If these ever drift
# apart, one of the two definitions has changed.
measurement_expected_n <- c(
  `E/e'` = nrow(surv_ee),
  LAVi = nrow(surv_lavi),
  TRVmax = nrow(surv_trv)
)
stopifnot(
  nrow(measurement_by_outcome) == length(measurement_disclosure_params),
  all(measurement_by_outcome$`Total measured` ==
    measurement_expected_n[measurement_by_outcome$Parameter]),
  all(measurement_by_outcome$`% measured among events` <= 100),
  !any(is.na(measurement_by_outcome$`Adjusted P value`))
)

write.csv(
  measurement_by_outcome,
  file.path(OUT_DIR, "Table_Measurement_By_Outcome.csv"),
  row.names = FALSE
)
print(as.data.frame(measurement_by_outcome), row.names = FALSE)

# ── Consistency check: the nested measurement cohorts ─────────────────────────
# The primary Cox models analyse each parameter in its own available-case
# sample. Measurement is not independent of the outcome here (card 1.2), so the
# same three exposures are refitted in the two nested cohorts that require the
# parameters to have been measured together: the common E/e' + LAVi cohort and
# the all-three-measured cohort (owner decision, card 2.3 of
# 8_Docs/DECISION_CARDS_2026-09.md). One fitter and one adjustment set
# throughout, so the samples are the only thing that differs.
#
# For TRVmax the common cohort and the all-three-measured cohort are the same
# set of patients by construction, so those two rows agree by definition; the
# sample definition column states the requirement explicitly rather than
# leaving identical denominators looking like a duplication error.
cox_consistency_samples <- list(
  list(
    label = "Parameter available-case sample",
    definition = "Exposure, follow-up and covariates measured",
    require = character(0)
  ),
  list(
    label = "Common E/e' + LAVi cohort",
    definition = "Additionally requires both cohort-defining parameters",
    require = c("e_e_ave_z", "la_vol_index_z")
  ),
  list(
    label = "All-three-measured cohort",
    definition = "Additionally requires E/e', LAVi and TRVmax",
    require = c("e_e_ave_z", "la_vol_index_z", "tr_max_vel_z")
  )
)

cox_cohort_consistency <- map_dfr(
  names(primary_cox_exposures),
  function(exposure) {
    map_dfr(cox_consistency_samples, function(sample_spec) {
      sample_df <- surv_df %>%
        filter(complete.cases(across(all_of(unique(c(
          "follow_up_yrs",
          "hf_composite",
          exposure,
          "age",
          "Sex",
          "BMI",
          sample_spec$require
        ))))))
      fit <- fit_same_sample_cox(
        sample_df,
        exposure,
        c("age", "Sex", "BMI"),
        "Age/sex/BMI adjusted"
      )
      if (is.null(fit)) return(tibble())
      fit %>%
        transmute(
          Parameter = unname(primary_cox_exposures[[exposure]]),
          `Analysis sample` = sample_spec$label,
          `Sample definition` = sample_spec$definition,
          N,
          Events,
          HR,
          `CI low`,
          `CI high`,
          `P value`,
          `Significant at 0.05` = `P value` < 0.05
        )
    })
  }
)

# Each cohort is nested in the one before it, so N must not increase down the
# three samples for any parameter. A rise would mean a filter was dropped.
cox_consistency_nesting <- cox_cohort_consistency %>%
  mutate(
    sample_order = match(
      `Analysis sample`,
      vapply(cox_consistency_samples, function(x) x$label, character(1))
    )
  ) %>%
  arrange(Parameter, sample_order) %>%
  group_by(Parameter) %>%
  summarise(
    nested = all(diff(N) <= 0),
    rows = n(),
    .groups = "drop"
  )
stopifnot(
  nrow(cox_cohort_consistency) ==
    length(primary_cox_exposures) * length(cox_consistency_samples),
  all(cox_consistency_nesting$nested),
  all(cox_consistency_nesting$rows == length(cox_consistency_samples)),
  !any(is.na(cox_cohort_consistency$`P value`))
)

write.csv(
  cox_cohort_consistency,
  file.path(OUT_DIR, "Table_Consistency_MeasurementCohorts_Cox.csv"),
  row.names = FALSE
)
print(as.data.frame(cox_cohort_consistency), row.names = FALSE)

# ── Sensitivity: multiple imputation of the missing exposures ─────────────────
# Labelled sensitivity analysis, not a primary estimate (owner decision, card
# 1.2 of 8_Docs/DECISION_CARDS_2026-09.md). It answers what the hazard ratios
# would be if E/e' and LAVi were missing at random given age, sex, BMI, TRVmax
# availability and the outcome. TRVmax is never imputed; see R/
# imputation.R for why and for how it still contributes.
if (isTRUE(cfg$imputation$enabled)) {
  mi_cox_rows <- mi_cox_exposure_sensitivity(
    surv_df,
    exposures = c(e_e_ave = "E/e'", la_vol_index = "LAVi"),
    m = cfg$imputation$m,
    seed = cfg$imputation$seed
  )

  # Put the available-case estimate beside it, so the comparison is explicit
  # rather than left to the reader to assemble from two tables.
  mi_available_case <- tibble(
    Parameter = c("E/e'", "LAVi"),
    `Available-case N` = c(nrow(surv_ee), nrow(surv_lavi)),
    `Available-case HR` = c(
      unname(exp(coef(cox_ee)["e_e_ave_z"])),
      unname(exp(coef(cox_lavi)["la_vol_index_z"]))
    ),
    `Available-case P value` = c(
      summary(cox_ee)$coefficients["e_e_ave_z", "Pr(>|z|)"],
      summary(cox_lavi)$coefficients["la_vol_index_z", "Pr(>|z|)"]
    )
  )

  mi_cox_sensitivity <- mi_cox_rows %>%
    left_join(mi_available_case, by = "Parameter") %>%
    mutate(
      `Conclusion unchanged` =
        (`P value` < 0.05) == (`Available-case P value` < 0.05)
    )

  # The imputation frame is every patient with follow-up and complete
  # covariates, so it must be at least as large as either available-case
  # sample, and measured + imputed must account for all of it.
  stopifnot(
    nrow(mi_cox_sensitivity) == 2,
    all(mi_cox_sensitivity$`Exposure measured` +
      mi_cox_sensitivity$`Exposure imputed` ==
      mi_cox_sensitivity$`Patients in imputation frame`),
    all(mi_cox_sensitivity$`Patients in imputation frame` >=
      mi_cox_sensitivity$`Available-case N`),
    all(mi_cox_sensitivity$`Exposure measured` ==
      mi_cox_sensitivity$`Available-case N`),
    all(mi_cox_sensitivity$Imputations == cfg$imputation$m),
    all(is.finite(mi_cox_sensitivity$HR))
  )

  write.csv(
    mi_cox_sensitivity,
    file.path(OUT_DIR, "Table_Sensitivity_MultipleImputation_Cox.csv"),
    row.names = FALSE
  )
  print(as.data.frame(mi_cox_sensitivity), row.names = FALSE)
}

# ── Sensitivity: adults only ──────────────────────────────────────────────────
# Percent-predicted peak VO2 rests on adult reference equations, so every Cox
# exposure is refitted with patients below the age floor excluded (owner
# decision, card 1.3 of 8_Docs/DECISION_CARDS_2026-09.md). Both columns come
# from the same fitter and the same age/sex/BMI adjustment, so the only
# difference between them is the age restriction.
min_age_sensitivity <- cfg$cohort$min_age_sensitivity

cox_age_sensitivity_exposures <- c(
  la_vol_index_z = "LAVi Z-score",
  e_e_ave_z = "E/e'",
  tr_max_vel_z = "TRVmax",
  lvot_gradient_z = "Resting LVOT gradient",
  septal_thickness_z = "LV septal thickness",
  vo2_friend2_z = "Peak V̇O2 (per 1 SD lower)",
  vevco2_slope_z = "V̇E/V̇CO2 slope"
)

if (!is.null(min_age_sensitivity)) {
  cox_age_sensitivity <- map_dfr(
    names(cox_age_sensitivity_exposures),
    function(exposure) {
      exposure_df <- surv_df %>%
        filter(complete.cases(across(all_of(c(
          "follow_up_yrs",
          "hf_composite",
          exposure,
          "age",
          "Sex",
          "BMI"
        )))))
      adult_df <- exposure_df %>% filter(age >= min_age_sensitivity)
      all_ages_fit <- fit_same_sample_cox(
        exposure_df, exposure, c("age", "Sex", "BMI"), "All ages"
      )
      adults_fit <- fit_same_sample_cox(
        adult_df, exposure, c("age", "Sex", "BMI"), "Adults only"
      )
      if (is.null(all_ages_fit) || is.null(adults_fit)) return(tibble())
      tibble(
        Parameter = unname(cox_age_sensitivity_exposures[[exposure]]),
        `Age floor (years)` = min_age_sensitivity,
        `N all ages` = all_ages_fit$N,
        `Events all ages` = all_ages_fit$Events,
        `N adults only` = adults_fit$N,
        `Events adults only` = adults_fit$Events,
        `Patients excluded` = all_ages_fit$N - adults_fit$N,
        `HR all ages` = all_ages_fit$HR,
        `HR adults only` = adults_fit$HR,
        `CI low adults only` = adults_fit$`CI low`,
        `CI high adults only` = adults_fit$`CI high`,
        `P value all ages` = all_ages_fit$`P value`,
        `P value adults only` = adults_fit$`P value`,
        `Conclusion unchanged` =
          (all_ages_fit$`P value` < 0.05) == (adults_fit$`P value` < 0.05)
      )
    }
  )

  stopifnot(
    nrow(cox_age_sensitivity) == length(cox_age_sensitivity_exposures),
    all(cox_age_sensitivity$`N adults only` <= cox_age_sensitivity$`N all ages`),
    all(cox_age_sensitivity$`Events adults only` <=
      cox_age_sensitivity$`Events all ages`),
    !any(is.na(cox_age_sensitivity$`P value adults only`))
  )

  write.csv(
    cox_age_sensitivity,
    file.path(OUT_DIR, "Table_Sensitivity_AdultsOnly_Cox.csv"),
    row.names = FALSE
  )
  print(as.data.frame(cox_age_sensitivity), row.names = FALSE)
}

# Direct selection sensitivity for the cohort-defining E/e' and LAVi Cox
# models: same exposure and adjustment set, once in the primary outcomes cohort
# and once after additionally requiring TRVmax.
cox_cohort_selection_rows <- map_dfr(
  c("e_e_ave_z", "la_vol_index_z"),
  function(exposure) {
    base_df <- surv_df %>%
      filter(complete.cases(across(all_of(primary_cox_common_required))))
    sample_specs <- list(
      list(
        label = "Primary E/e' + LAVi outcomes cohort",
        data = base_df
      ),
      list(
        label = "TRVmax-complete outcomes subset",
        data = base_df %>% filter(!is.na(tr_max_vel_z))
      )
    )
    map_dfr(sample_specs, function(sample_spec) {
      fit_row <- fit_same_sample_cox(
        sample_spec$data,
        exposure,
        c("age", "Sex", "BMI"),
        "Age/sex/BMI adjusted"
      )
      if (is.null(fit_row)) return(tibble())
      fit_row %>%
        mutate(
          Parameter = unname(primary_cox_exposures[[exposure]]),
          `Analysis sample` = sample_spec$label,
          `Effect scale` = "Per 1 SD in the primary outcomes cohort",
          .before = 1
        )
    })
  }
)

write.csv(
  cox_cohort_selection_rows,
  paste0(
    file.path(AUDIT_DIR, ""),
    "Table_CohortSelection_Outcomes_EeLAVi.csv"
  ),
  row.names = FALSE
)

forest_rows <- bind_rows(
  if (!is.null(cox_lavi)) {
    tidy(cox_lavi, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "la_vol_index_z") %>%
      transmute(exposure = "LAVi", estimate, conf.low, conf.high, p.value)
  },
  if (!is.null(cox_ee)) {
    tidy(cox_ee, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "e_e_ave_z") %>%
      transmute(exposure = "E/e'", estimate, conf.low, conf.high, p.value)
  },
  if (!is.null(cox_trv)) {
    tidy(cox_trv, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "tr_max_vel_z") %>%
      transmute(exposure = "TRVmax", estimate, conf.low, conf.high, p.value)
  },
  if (!is.null(cox_lvot)) {
    tidy(cox_lvot, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "lvot_gradient_z") %>%
      transmute(
        exposure = "Resting LVOT gradient",
        estimate,
        conf.low,
        conf.high,
        p.value
      )
  },
  if (!is.null(cox_septal)) {
    tidy(cox_septal, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "septal_thickness_z") %>%
      transmute(
        exposure = "LV septal thickness",
        estimate,
        conf.low,
        conf.high,
        p.value
      )
  },
  if (!is.null(cox_vo2)) {
    tidy(cox_vo2, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "vo2_friend2_z") %>%
      transmute(
        exposure = "Peak V̇O2 (per SD lower)",
        estimate,
        conf.low,
        conf.high,
        p.value
      )
  },
  if (!is.null(cox_vevco2)) {
    tidy(cox_vevco2, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "vevco2_slope_z") %>%
      transmute(
        exposure = "V̇E/V̇CO2 slope",
        estimate,
        conf.low,
        conf.high,
        p.value
      )
  }
) %>%
  arrange(desc(estimate)) %>%
  mutate(
    exposure_wrapped = case_when(
      str_detect(as.character(exposure), regex("^Peak", ignore_case = TRUE)) ~
        "Lower peak V̇O₂",
      str_detect(as.character(exposure), regex("LVOT gradient", ignore_case = TRUE)) ~
        "Resting LVOT\ngradient",
      TRUE ~ str_wrap(as.character(exposure), width = 16)
    ),
    hr_label = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    forest_color = ifelse(
      !is.na(p.value) & p.value < 0.05,
      fig5_palette["highlight"],
      fig5_palette["reference"]
    ),
    exposure_wrapped = factor(
      exposure_wrapped,
      levels = rev(unique(exposure_wrapped))
    ),
    exposure = factor(exposure, levels = rev(unique(exposure)))
  )

if (nrow(forest_rows) > 0) {
  xmax_data <- max(forest_rows$conf.high, na.rm = TRUE)
  xmax <- xmax_data * 1.25
  xmin <- max(0.5, min(forest_rows$conf.low, na.rm = TRUE) * 0.85)

  p_hf_forest <- ggplot(forest_rows, aes(x = estimate, y = exposure_wrapped)) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = fig5_palette["refline"],
      linewidth = 0.5
    ) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high, color = forest_color),
      height = 0.14,
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    geom_point(
      size = 2.8,
      shape = 21,
      stroke = 0.4,
      aes(fill = forest_color, color = forest_color),
      show.legend = FALSE
    ) +
    geom_text(
      aes(x = pmin(conf.high * 1.035, xmax_data * 1.13), label = hr_label),
      hjust = 0,
      size = 2.2,
      color = fig5_palette["text"]
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_x_continuous(
      limits = c(xmin, xmax),
      breaks = c(0.5, 1, 1.5, 2, 3, 4),
      expand = expansion(mult = c(0.01, 0.005))
    ) +
    labs(
      title = "Adjusted Hazard Ratios",
      x = "Hazard Ratio per SD",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_jacc() +
    theme(
      legend.position = "none",
      plot.margin = margin(4, 0, 4, 10),
      axis.text.y = element_text(
        angle = 0,
        hjust = 1,
        size = 6.8,
        lineheight = 0.62,
        color = fig5_palette["text"]
      ),
      axis.text.x = element_text(color = fig5_palette["text"]),
      axis.title.x = element_text(color = fig5_palette["text"]),
      axis.line = element_line(color = fig5_palette["text"], linewidth = 0.35)
    )

  # Event-count panel for the lower-left position in the original Figure 5
  # geometry. Counts identify the event that first determined the primary
  # composite; surgical myectomy is separate and is not part of that composite.
  figure5_event_rows <- tibble(
    event = c(
      "Primary composite",
      "Acute heart failure",
      "Heart transplantation",
      "All-cause death",
      "Same-day multiple components",
      "Surgical myectomy (separate)"
    ),
    n = c(
      sum(surv_df$hf_composite == 1, na.rm = TRUE),
      sum(surv_df$first_event_component == "Acute heart failure"),
      sum(surv_df$first_event_component == "Heart transplantation"),
      sum(surv_df$first_event_component == "All-cause death"),
      sum(surv_df$first_event_component == "Multiple components on same day"),
      sum(surv_df$myectomy_event == 1, na.rm = TRUE)
    )
  ) %>%
    mutate(y = rev(seq_len(n())))

  p_event_counts <- ggplot(figure5_event_rows) +
    geom_hline(
      yintercept = figure5_event_rows$y - 0.5,
      color = "#E3E7EA",
      linewidth = 0.3
    ) +
    geom_text(
      aes(x = 0, y = y, label = event),
      hjust = 0,
      size = 2.35,
      color = fig5_palette["text"]
    ) +
    geom_text(
      aes(x = 1, y = y, label = n),
      hjust = 1,
      size = 2.45,
      fontface = "bold",
      color = ifelse(
        figure5_event_rows$event == "Primary composite",
        fig5_palette["highlight"],
        fig5_palette["forest"]
      )
    ) +
    annotate(
      "text",
      x = 1,
      y = max(figure5_event_rows$y) + 0.72,
      label = "n",
      hjust = 1,
      size = 2.35,
      fontface = "bold",
      color = fig5_palette["text"]
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(
      limits = c(0.35, max(figure5_event_rows$y) + 1.05),
      expand = c(0, 0)
    ) +
    labs(
      title = "Clinical Event Counts",
      subtitle = sprintf("First composite event; cohort N = %d", nrow(surv_df))
    ) +
    theme_void(base_family = "Arial", base_size = 8.5) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 8.8,
        hjust = 0.5,
        margin = margin(b = 2)
      ),
      plot.subtitle = element_text(
        size = 7.0,
        hjust = 0.5,
        color = "grey40",
        margin = margin(b = 5)
      ),
      plot.margin = margin(4, 7, 4, 4)
    )

  # Panel E shows the primary composite only. The myectomy (septal reduction)
  # models were previously appended here as three extra rows; they are reported
  # in Table_Figure5_Myectomy_Cox.csv instead, at the owner's request, so the
  # panel carries one endpoint rather than two.
  combined_forest_rows <- forest_rows %>%
    transmute(
      exposure = as.character(exposure_wrapped),
      estimate,
      conf.low,
      conf.high,
      p.value,
      hr_label,
      forest_color
    )

  combined_order <- as.character(forest_rows$exposure_wrapped)
  combined_forest_rows <- combined_forest_rows %>%
    mutate(exposure = factor(exposure, levels = rev(unique(combined_order))))

  combined_xmax_data <- max(combined_forest_rows$conf.high, na.rm = TRUE)
  combined_xmin <- max(
    0.4,
    min(combined_forest_rows$conf.low, na.rm = TRUE) * 0.85
  )
  combined_xmax <- combined_xmax_data * 1.28

  p_outcomes_forest <- ggplot(
    combined_forest_rows,
    aes(x = estimate, y = exposure)
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = fig5_palette["refline"],
      linewidth = 0.5
    ) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high, color = forest_color),
      height = 0.12,
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    geom_point(
      aes(fill = forest_color, color = forest_color),
      size = 2.45,
      shape = 21,
      stroke = 0.35,
      show.legend = FALSE
    ) +
    geom_text(
      aes(
        x = pmin(conf.high * 1.03, combined_xmax_data * 1.10),
        label = hr_label
      ),
      hjust = 0,
      size = 1.82,
      color = fig5_palette["text"]
    ) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_x_continuous(
      limits = c(combined_xmin, combined_xmax),
      breaks = c(0.5, 1, 1.5, 2),
      expand = expansion(mult = c(0.01, 0.005))
    ) +
    labs(
      title = "Adjusted Hazard Ratios",
      subtitle = "Primary composite; adjusted for age, sex and BMI",
      x = "Hazard Ratio per SD",
      y = NULL
    ) +
    coord_cartesian(clip = "off") +
    theme_jacc(base_size = 8.5) +
    theme(
      legend.position = "none",
      plot.margin = margin(4, 0, 4, 10),
      plot.title = element_text(size = 8.8),
      plot.subtitle = element_text(size = 6.8),
      axis.text.y = element_text(
        angle = 0,
        hjust = 1,
        size = 6.4,
        lineheight = 0.78,
        color = fig5_palette["text"]
      ),
      axis.text.x = element_text(color = fig5_palette["text"]),
      axis.title.x = element_text(color = fig5_palette["text"]),
      axis.line = element_line(color = fig5_palette["text"], linewidth = 0.35)
    )

  km_placeholder <- function(label) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = label, size = 3.1) +
      theme_void()
  }

  fig5 <- wrap_plots(
    A = if (!is.null(p_km_ee)) {
      p_km_ee
    } else {
      km_placeholder("E/e' KM unavailable")
    },
    B = if (!is.null(p_km_lavi)) {
      p_km_lavi
    } else {
      km_placeholder("LAVi KM unavailable")
    },
    C = if (!is.null(p_km_vevco2)) {
      p_km_vevco2
    } else {
      km_placeholder("V̇E/V̇CO₂ slope KM unavailable")
    },
    D = if (!is.null(p_km_vo2)) {
      p_km_vo2
    } else {
      km_placeholder("Peak V̇O₂ KM unavailable")
    },
    E = patchwork::free(p_outcomes_forest, side = "l", type = "space"),
    design = "
ABC
DEE
"
  ) +
    plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(face = "bold", size = 11))

  show_and_save_jacc(
    fig5,
    file.path(OUT_DIR, "Figure5_HeartFailure_Outcomes.pdf"),
    w = 7.0,
    h = 6.8
  )

  figure5_preprint_width <- 7.2
  figure5_preprint_height <- 5.0
  save_jacc(
    fig5,
    file.path(OUT_DIR, "Figure5_HeartFailure_Outcomes_Preprint.pdf"),
    w = figure5_preprint_width,
    h = figure5_preprint_height
  )
  save_jacc(
    fig5,
    file.path(OUT_DIR, "Figure5_HeartFailure_Outcomes_Preprint.pdf"),
    w = figure5_preprint_width,
    h = figure5_preprint_height
  )

  legend_fig5_pdf <- paste0(
    "**Figure 5. Clinical Outcomes Across LV Diastolic Parameters.** ",
    "The primary endpoint was a dated acute heart failure event, heart transplantation, or death; ",
    "primary Cox models adjusted for age, sex, and BMI and report hazard ratios per standard deviation; peak V̇O$_2$ was inverse-coded so its hazard ratio is per SD lower fitness. ",
    sprintf(
      paste0(
        "Among %d patients, %d first composite events comprised %d acute heart-failure events, ",
        "%d heart transplantation, %d all-cause deaths, and %d same-day multiple-component event; ",
        "%d surgical myectomies occurred and were treated as a separate secondary outcome. "
      ),
      nrow(surv_df),
      sum(surv_df$hf_composite == 1, na.rm = TRUE),
      sum(surv_df$first_event_component == "Acute heart failure"),
      sum(surv_df$first_event_component == "Heart transplantation"),
      sum(surv_df$first_event_component == "All-cause death"),
      sum(surv_df$first_event_component == "Multiple components on same day"),
      sum(surv_df$myectomy_event == 1, na.rm = TRUE)
    ),
    "**(A)** Kaplan-Meier event-free survival stratified by E/e\u2019 \u226414 versus >14. ",
    "**(B)** Kaplan-Meier event-free survival stratified by LAVi \u226434 mL/m\u00B2 versus >34 mL/m\u00B2. ",
    "**(C)** Kaplan-Meier event-free survival stratified by V̇E/V̇CO$_2$ slope \u226430 versus >30. ",
    "**(D)** Kaplan-Meier event-free survival stratified by peak V̇O$_2$ <80% versus ≥80% FRIEND 2.0 predicted. ",
    "**(E)** Adjusted hazard ratios (95% CI) for the listed parameters; filled circles indicate P < 0.05. ",
    sprintf(
      paste0(
        "Forest-model N/events were E/e’ %d/%d, LAVi %d/%d, TRVmax %d/%d, ",
        "resting LVOT gradient %d/%d, LV septal thickness %d/%d, ",
        "peak V̇O$_2$ %d/%d, and V̇E/V̇CO$_2$ slope %d/%d. "
      ),
      nrow(surv_ee), sum(surv_ee$hf_composite == 1, na.rm = TRUE),
      nrow(surv_lavi), sum(surv_lavi$hf_composite == 1, na.rm = TRUE),
      nrow(surv_trv), sum(surv_trv$hf_composite == 1, na.rm = TRUE),
      nrow(surv_lvot), sum(surv_lvot$hf_composite == 1, na.rm = TRUE),
      nrow(surv_septal), sum(surv_septal$hf_composite == 1, na.rm = TRUE),
      nrow(surv_vo2), sum(surv_vo2$hf_composite == 1, na.rm = TRUE),
      nrow(surv_vevco2), sum(surv_vevco2$hf_composite == 1, na.rm = TRUE)
    ),
    "The time-constant V̇E/V̇CO$_2$ hazard ratio violated proportional hazards and is interpreted through a supplemental log-time-varying coefficient analysis. ",
    "Septal reduction therapy is analysed separately in the supplement with LVOT-gradient-adjusted models; generic ablation was omitted because procedure type could not be validated."
  )

  save_jacc_with_embedded_legend(
    fig5,
    file.path(OUT_DIR, "Figure5_HeartFailure_Outcomes.pdf"),
    legend_text = legend_fig5_pdf,
    w = 7.0,
    h = 6.8
  )
  save_jacc_with_embedded_legend(
    fig5,
    file.path(OUT_DIR, "Figure5_HeartFailure_Outcomes.pdf"),
    legend_text = legend_fig5_pdf,
    w = 7.0,
    h = 6.8
  )

}

# ══ Supplemental Table: Proportional-Hazards Diagnostics ══════════════════════
# (legacy chunk: figure5_ph_diagnostics)
# Schoenfeld-residual test for proportional-hazards assumption on each Cox
# model fit in Figure 5. Reports global and per-covariate chi-square p-values.
# If the global test rejects PH at alpha = 0.05 for a key diastolic
# predictor, the discussion section flags this and recommends a time-
# interaction sensitivity analysis.

ph_specs <- list(
  list(label = "LAVi Z-score", model = cox_lavi, surv_data = surv_lavi),
  list(label = "E/e' Z-score", model = cox_ee, surv_data = surv_ee),
  list(label = "TRVmax Z-score", model = cox_trv, surv_data = surv_trv),
  list(label = "LVOT Z-score", model = cox_lvot, surv_data = surv_lvot),
  list(label = "Septal Z-score", model = cox_septal, surv_data = surv_septal),
  list(label = "Peak VO2 Z-score", model = cox_vo2, surv_data = surv_vo2),
  list(label = "VE/VCO2 Z-score", model = cox_vevco2, surv_data = surv_vevco2)
)

ph_rows <- purrr::map_dfr(ph_specs, function(spec) {
  if (is.null(spec$model)) {
    return(NULL)
  }
  zph <- tryCatch(
    survival::cox.zph(spec$model, transform = "km"),
    error = function(e) NULL
  )
  if (is.null(zph)) {
    return(NULL)
  }
  tab <- as.data.frame(zph$table)
  tab$term <- rownames(tab)
  tibble(
    Model = spec$label,
    Term = tab$term,
    `Chi-square` = sprintf("%.2f", tab$chisq),
    df = tab$df,
    `P value` = ifelse(tab$p < 0.001, "<0.001", sprintf("%.3f", tab$p)),
    p_num = tab$p
  )
})

write.csv(
  ph_rows %>% select(-p_num),
  file.path(OUT_DIR, "TableS_Cox_PH_Diagnostics.csv"),
  row.names = FALSE
)

ph_violations <- ph_rows %>%
  filter(Term == "GLOBAL", !is.na(p_num), p_num < 0.05) %>%
  pull(Model)


cat("\n\n")
if (length(ph_violations) > 0) {
  cat(sprintf(
    "**Proportional-hazards violations detected** in the following model(s): %s. ",
    paste(ph_violations, collapse = ", ")
  ))
  cat(
    "These models should be interpreted with caution; results are reported with ",
    "a time-stratified sensitivity in the discussion.\n"
  )
} else {
  cat(
    "**No global PH violations detected** (all global tests P >= 0.05). ",
    "Hazard ratios are interpretable as time-constant effects.\n"
  )
}


# ══ Supplemental Table: Figure 5D Cox Model Summary ═══════════════════════════
# (legacy chunk: figure5_summary_table)
# ── table: Cox HRs (95% CI, P) for each predictor/endpoint ────────────────────
figure5_summary_table_raw <- bind_rows(
  if (!is.null(cox_fp)) {
    tidy(cox_fp, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(
        term %in%
          c(
            "fp_classElevated",
            "age",
            "SexFemale",
            "BMI",
            "bb_any",
            "ndhp_ccb_any"
          )
      ) %>%
      transmute(
        model_order = 1L,
        term_order = match(
          term,
          c(
            "fp_classElevated",
            "age",
            "SexFemale",
            "BMI",
            "bb_any",
            "ndhp_ccb_any"
          )
        ),
        Model = "Composite LVDD model",
        Parameter = case_when(
          term == "fp_classElevated" ~ "LVDD",
          term == "age" ~ "Age",
          term == "SexFemale" ~ "Female sex",
          term == "BMI" ~ "BMI",
          term == "bb_any" ~ "β-blocker use",
          term == "ndhp_ccb_any" ~ "Non-DHP CCB use",
          TRUE ~ term
        ),
        Scale = case_when(
          term == "fp_classElevated" ~ "Elevated vs not elevated",
          term == "age" ~ "Per 1-year increase",
          term == "SexFemale" ~ "Female vs male",
          term == "BMI" ~ "Per 1 kg/m² increase",
          term %in% c("bb_any", "ndhp_ccb_any") ~ "Use vs no use",
          TRUE ~ "As modeled"
        ),
        Sample = sprintf(
          "%d / %d",
          nrow(surv_fp_model),
          sum(surv_fp_model$hf_composite == 1, na.rm = TRUE)
        ),
        estimate,
        conf.low,
        conf.high,
        p_num = p.value
      )
  },
  if (!is.null(cox_lavi)) {
    tidy(cox_lavi, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "la_vol_index_z") %>%
      transmute(
        model_order = 2L,
        term_order = 1L,
        Model = "LAVi model",
        Parameter = "LAVi Z-score",
        Scale = "Per 1 SD higher",
        Sample = sprintf(
          "%d / %d",
          nrow(surv_lavi),
          sum(surv_lavi$hf_composite == 1, na.rm = TRUE)
        ),
        estimate,
        conf.low,
        conf.high,
        p_num = p.value
      )
  },
  if (!is.null(cox_ee)) {
    tidy(cox_ee, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "e_e_ave_z") %>%
      transmute(
        model_order = 3L,
        term_order = 1L,
        Model = "E/e' model",
        Parameter = "E/e'",
        Scale = "Per 1 SD higher",
        Sample = sprintf(
          "%d / %d",
          nrow(surv_ee),
          sum(surv_ee$hf_composite == 1, na.rm = TRUE)
        ),
        estimate,
        conf.low,
        conf.high,
        p_num = p.value
      )
  },
  if (!is.null(cox_trv)) {
    tidy(cox_trv, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "tr_max_vel_z") %>%
      transmute(
        model_order = 4L,
        term_order = 1L,
        Model = "TRVmax model",
        Parameter = "TRVmax",
        Scale = "Per 1 SD higher",
        Sample = sprintf(
          "%d / %d",
          nrow(surv_trv),
          sum(surv_trv$hf_composite == 1, na.rm = TRUE)
        ),
        estimate,
        conf.low,
        conf.high,
        p_num = p.value
      )
  },
  if (!is.null(cox_lvot)) {
    tidy(cox_lvot, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "lvot_gradient_z") %>%
      transmute(
        model_order = 5L,
        term_order = 1L,
        Model = "LVOT gradient model",
        Parameter = "Resting LVOT gradient",
        Scale = "Per 1 SD higher",
        Sample = sprintf(
          "%d / %d",
          nrow(surv_lvot),
          sum(surv_lvot$hf_composite == 1, na.rm = TRUE)
        ),
        estimate,
        conf.low,
        conf.high,
        p_num = p.value
      )
  },
  if (!is.null(cox_septal)) {
    tidy(cox_septal, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "septal_thickness_z") %>%
      transmute(
        model_order = 6L,
        term_order = 1L,
        Model = "LV septal thickness model",
        Parameter = "LV septal thickness",
        Scale = "Per 1 SD higher",
        Sample = sprintf(
          "%d / %d",
          nrow(surv_septal),
          sum(surv_septal$hf_composite == 1, na.rm = TRUE)
        ),
        estimate,
        conf.low,
        conf.high,
        p_num = p.value
      )
  },
  if (!is.null(cox_vo2)) {
    tidy(cox_vo2, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "vo2_friend2_z") %>%
      transmute(
        model_order = 7L,
        term_order = 1L,
        Model = "Peak V̇O2 model",
        Parameter = "Peak V̇O2",
        Scale = "Per 1 SD lower",
        Sample = sprintf(
          "%d / %d",
          nrow(surv_vo2),
          sum(surv_vo2$hf_composite == 1, na.rm = TRUE)
        ),
        estimate,
        conf.low,
        conf.high,
        p_num = p.value
      )
  },
  if (!is.null(cox_vevco2)) {
    tidy(cox_vevco2, exponentiate = TRUE, conf.int = TRUE) %>%
      filter(term == "vevco2_slope_z") %>%
      transmute(
        model_order = 8L,
        term_order = 1L,
        Model = "V̇E/V̇CO2 slope model",
        Parameter = "V̇E/V̇CO2 slope",
        Scale = "Per 1 SD higher",
        Sample = sprintf(
          "%d / %d",
          nrow(surv_vevco2),
          sum(surv_vevco2$hf_composite == 1, na.rm = TRUE)
        ),
        estimate,
        conf.low,
        conf.high,
        p_num = p.value
      )
  }
) %>%
  arrange(model_order, term_order)

# Multiplicity family for the outcome analyses: one exposure per Cox model
# (owner decision, card 1.1 of 8_Docs/DECISION_CARDS_2026-09.md). The composite
# LVDD model also reports its adjustment covariates, which are not hypotheses in
# this family, so they are named here and excluded from the adjustment rather
# than assumed absent.
figure5_covariate_parameters <- c(
  "Age",
  "Female sex",
  "BMI",
  "β-blocker use",
  "Non-DHP CCB use"
)

# NA-safe: p.adjust() would count NA entries in its denominator, so the
# adjustment is computed on the exposure rows alone.
adjust_exposure_q <- function(p, is_exposure) {
  q <- rep(NA_real_, length(p))
  q[is_exposure] <- p.adjust(p[is_exposure], method = "BH")
  q
}

figure5_summary_table_display <- figure5_summary_table_raw %>%
  mutate(
    is_exposure = !Parameter %in% figure5_covariate_parameters,
    q_num = adjust_exposure_q(p_num, is_exposure)
  ) %>%
  transmute(
    Model,
    Parameter,
    Scale,
    Sample,
    `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    `P value` = format_p_scientific(p_num, digits = 2),
    `BH q value` = ifelse(
      is_exposure,
      format_p_scientific(q_num, digits = 2),
      "not in family"
    ),
    p_num,
    q_num,
    is_exposure
  )
# BH across exposures is only valid if each contributes exactly one test, so
# assert that rather than a hard-coded row count.
stopifnot(
  sum(figure5_summary_table_display$is_exposure) ==
    n_distinct(figure5_summary_table_display$Model[
      figure5_summary_table_display$is_exposure
    ]),
  !any(is.na(figure5_summary_table_display$q_num[
    figure5_summary_table_display$is_exposure
  ])),
  all(
    figure5_summary_table_display$q_num[
      figure5_summary_table_display$is_exposure
    ] >=
      figure5_summary_table_display$p_num[
        figure5_summary_table_display$is_exposure
      ]
  )
)

figure5_summary_break_rows <- figure5_summary_table_display %>%
  count(Model, name = "n_rows") %>%
  pull(n_rows) %>%
  cumsum()

figure5_summary_export <- figure5_summary_table_display %>%
  select(-p_num, -q_num, -is_exposure)

figure5_summary_ft <- figure5_summary_table_display %>%
  flextable(
    col_keys = c(
      "Model",
      "Parameter",
      "Scale",
      "Sample",
      "HR (95% CI)",
      "P value",
      "BH q value"
    )
  ) %>%
  set_header_labels(
    Model = "Cox model",
    Parameter = "Figure 5D parameter",
    Scale = "Contrast / scale",
    Sample = "N / Events"
  ) %>%
  merge_v(j = c("Model", "Sample")) %>%
  valign(j = c("Model", "Sample"), valign = "top", part = "body") %>%
  align(j = c("Model", "Parameter", "Scale"), align = "left", part = "all") %>%
  align(
    j = c("Sample", "HR (95% CI)", "P value", "BH q value"),
    align = "center",
    part = "all"
  ) %>%
  bold(i = ~ p_num < 0.05, j = "P value", part = "body") %>%
  bold(i = ~ !is.na(q_num) & q_num < 0.05, j = "BH q value", part = "body") %>%
  format_pub_table(
    caption = paste0(
      "Supplemental Table. Multivariable Cox statistics for parameters displayed in Figure 5D. ",
      "Each parameter was entered in a separate Cox model adjusted for age, sex, and BMI. ",
      "Peak V̇O₂ Z-score is sign-reversed (`scale(-V̇O₂)`) so that HR > 1 indicates lower fitness. ",
      "Q values are Benjamini-Hochberg adjusted across the seven exposures tested against the composite endpoint."
    )
  ) %>%
  hline(
    i = figure5_summary_break_rows,
    border = officer::fp_border(color = "grey55", width = 0.4),
    part = "body"
  )


write.csv(
  figure5_summary_export,
  file.path(OUT_DIR, "TableS_Figure5_Cox_Summary.csv"),
  row.names = FALSE
)
figure5_summary_ft %>%
  save_as_docx(path = file.path(OUT_DIR, "TableS_Figure5_Cox_Summary.docx"))


# ══ Medication-Adjusted and Intervention-Censored Sensitivities ═══════════════
# (legacy chunk: figure5_sensitivity)
# ── Fig 5 sensitivity: re-fit Cox censoring at major interventions ────────────
# myectomy/ablation/transplant alter physiology; censor follow-up there to check
# the HRs aren't artifacts of post-intervention events
summarize_intervention_row <- function(
  df,
  label,
  pre_var,
  post_var,
  post_date_var
) {
  post_time_yrs <- as.numeric(df[[post_date_var]] - df$cpx_test_date) / 365.25
  post_time_yrs <- post_time_yrs[is.finite(post_time_yrs) & post_time_yrs > 0]
  denom_n <- max(nrow(df), 1)

  tibble(
    Intervention = label,
    `Pre-index n (%)` = sprintf(
      "%d (%.1f%%)",
      sum(df[[pre_var]] == 1, na.rm = TRUE),
      100 * sum(df[[pre_var]] == 1, na.rm = TRUE) / denom_n
    ),
    `Post-index n (%)` = sprintf(
      "%d (%.1f%%)",
      sum(df[[post_var]] == 1, na.rm = TRUE),
      100 * sum(df[[post_var]] == 1, na.rm = TRUE) / denom_n
    ),
    `Median years to post event [IQR]` = ifelse(
      length(post_time_yrs) == 0,
      "NA",
      sprintf(
        "%.1f [%.1f, %.1f]",
        median(post_time_yrs, na.rm = TRUE),
        quantile(post_time_yrs, 0.25, na.rm = TRUE),
        quantile(post_time_yrs, 0.75, na.rm = TRUE)
      )
    )
  )
}

intervention_summary_table <- bind_rows(
  summarize_intervention_row(
    baseline_df,
    "Septal reduction surgery",
    "pre_septal_reduction_surgery",
    "post_septal_reduction_surgery",
    "post_septal_reduction_surgery_date"
  ),
  summarize_intervention_row(
    baseline_df,
    "Ablation surgery",
    "pre_ablation_surgery",
    "post_ablation_surgery",
    "post_ablation_surgery_date"
  ),
  summarize_intervention_row(
    baseline_df,
    "Defibrillator",
    "pre_defibrillator",
    "post_defibrillator",
    "post_defibrillator_date"
  ),
  summarize_intervention_row(
    baseline_df,
    "Pacemaker",
    "pre_pacemaker",
    "post_pacemaker",
    "post_pacemaker_date"
  ),
  summarize_intervention_row(
    baseline_df,
    "Heart transplant",
    "pre_heart_transplant",
    "post_heart_transplant",
    "post_heart_transplant_date"
  )
)

write.csv(
  intervention_summary_table,
  file.path(OUT_DIR, "Table_Sensitivity_Interventions_Summary.csv"),
  row.names = FALSE
)
intervention_summary_table %>%
  flextable() %>%
  bold(part = "header") %>%
  align(
    j = c("Intervention", "Median years to post event [IQR]"),
    align = "left",
    part = "all"
  ) %>%
  theme_booktabs() %>%
  autofit() %>%
  save_as_docx(
    path = file.path(OUT_DIR, "Table_Sensitivity_Interventions_Summary.docx")
  )

surv_sensitivity_df <- surv_df %>%
  mutate(
    hf_or_death_no_transplant_event = as.integer(
      hf_or_death_no_transplant == 1 &
        !is.na(event_time_no_transplant_yrs) &
        event_time_no_transplant_yrs > 0
    ),
    follow_up_no_transplant_yrs = case_when(
      hf_or_death_no_transplant_event == 1 ~ event_time_no_transplant_yrs,
      !is.na(censor_time_yrs) & censor_time_yrs > 0 ~ censor_time_yrs,
      TRUE ~ NA_real_
    ),
    follow_up_no_transplant_yrs = ifelse(
      is.na(follow_up_no_transplant_yrs) | follow_up_no_transplant_yrs <= 0,
      0.01,
      follow_up_no_transplant_yrs
    ),
    event_before_intervention = hf_or_death_no_transplant_event == 1 &
      (is.na(major_intervention_time_yrs) |
        major_intervention_time_yrs <= 0 |
        event_time_no_transplant_yrs <= major_intervention_time_yrs),
    follow_up_intervention_censored_yrs = pmin(
      ifelse(
        hf_or_death_no_transplant_event == 1 &
          !is.na(event_time_no_transplant_yrs) &
          event_time_no_transplant_yrs > 0,
        event_time_no_transplant_yrs,
        Inf
      ),
      ifelse(
        !is.na(censor_time_yrs) & censor_time_yrs > 0,
        censor_time_yrs,
        Inf
      ),
      ifelse(
        !is.na(major_intervention_time_yrs) & major_intervention_time_yrs > 0,
        major_intervention_time_yrs,
        Inf
      ),
      na.rm = TRUE
    ),
    follow_up_intervention_censored_yrs = ifelse(
      is.infinite(follow_up_intervention_censored_yrs) |
        follow_up_intervention_censored_yrs <= 0,
      0.01,
      follow_up_intervention_censored_yrs
    ),
    hf_or_death_intervention_censored_event = as.integer(
      event_before_intervention
    )
  )

tidy_sensitivity_cox <- function(
  model,
  term,
  exposure_label,
  scenario_label,
  n_obs,
  n_events
) {
  if (is.null(model)) {
    return(tibble())
  }
  tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == !!term) %>%
    transmute(
      Scenario = scenario_label,
      Exposure = exposure_label,
      N = n_obs,
      Events = n_events,
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p_value = p.value
    )
}

run_figure5_sensitivity_suite <- function(
  data,
  scenario_label,
  time_var,
  event_var,
  med_terms = c("bb_any", "ndhp_ccb_any")
) {
  out_rows <- list()
  # The medication adjustment is a parameter so this suite can be run either
  # complete-case or with the three-level missing-indicator covariates (owner
  # decision, card 1.5 of 8_Docs/DECISION_CARDS_2026-09.md).
  med_rhs <- paste(med_terms, collapse = " + ")

  surv_fp_sens <- data %>% filter(!is.na(fp_class))
  if (
    nrow(surv_fp_sens) >= 20 &&
      sum(surv_fp_sens[[event_var]] == 1, na.rm = TRUE) >= 5
  ) {
    fit <- tryCatch(
      coxph(
        as.formula(paste0(
          "Surv(",
          time_var,
          ", ",
          event_var,
          ") ~ fp_class + age + Sex + BMI + ", med_rhs
        )),
        data = surv_fp_sens
      ),
      error = function(e) NULL
    )
    out_rows[[length(out_rows) + 1]] <- tidy_sensitivity_cox(
      fit,
      "fp_classElevated",
      "LVDD",
      scenario_label,
      nrow(surv_fp_sens),
      sum(surv_fp_sens[[event_var]] == 1, na.rm = TRUE)
    )
  }

  surv_z_sens <- data %>%
    left_join(z_source, by = "ID") %>%
    filter(!is.na(dd_zscore))
  if (
    nrow(surv_z_sens) >= 20 &&
      sum(surv_z_sens[[event_var]] == 1, na.rm = TRUE) >= 5
  ) {
    fit <- tryCatch(
      coxph(
        as.formula(paste0(
          "Surv(",
          time_var,
          ", ",
          event_var,
          ") ~ dd_zscore + age + Sex + BMI + ", med_rhs
        )),
        data = surv_z_sens
      ),
      error = function(e) NULL
    )
    out_rows[[length(out_rows) + 1]] <- tidy_sensitivity_cox(
      fit,
      "dd_zscore",
      "LAVi",
      scenario_label,
      nrow(surv_z_sens),
      sum(surv_z_sens[[event_var]] == 1, na.rm = TRUE)
    )
  }

  surv_ee_sens <- data %>% filter(!is.na(e_e_ave_z))
  if (
    nrow(surv_ee_sens) >= 20 &&
      sum(surv_ee_sens[[event_var]] == 1, na.rm = TRUE) >= 5
  ) {
    fit <- tryCatch(
      coxph(
        as.formula(paste0(
          "Surv(",
          time_var,
          ", ",
          event_var,
          ") ~ e_e_ave_z + age + Sex + BMI + ", med_rhs
        )),
        data = surv_ee_sens
      ),
      error = function(e) NULL
    )
    out_rows[[length(out_rows) + 1]] <- tidy_sensitivity_cox(
      fit,
      "e_e_ave_z",
      "E/e'",
      scenario_label,
      nrow(surv_ee_sens),
      sum(surv_ee_sens[[event_var]] == 1, na.rm = TRUE)
    )
  }

  surv_trv_sens <- data %>% filter(!is.na(tr_max_vel_z))
  if (
    nrow(surv_trv_sens) >= 20 &&
      sum(surv_trv_sens[[event_var]] == 1, na.rm = TRUE) >= 5
  ) {
    fit <- tryCatch(
      coxph(
        as.formula(paste0(
          "Surv(",
          time_var,
          ", ",
          event_var,
          ") ~ tr_max_vel_z + age + Sex + BMI + ", med_rhs
        )),
        data = surv_trv_sens
      ),
      error = function(e) NULL
    )
    out_rows[[length(out_rows) + 1]] <- tidy_sensitivity_cox(
      fit,
      "tr_max_vel_z",
      "TRVmax",
      scenario_label,
      nrow(surv_trv_sens),
      sum(surv_trv_sens[[event_var]] == 1, na.rm = TRUE)
    )
  }

  sensitivity_specs <- tribble(
    ~filter_var          , ~term                , ~rhs                                                    , ~label                ,
    "lvot_gradient_z"    , "lvot_gradient_z"    , paste0("lvot_gradient_z + age + Sex + BMI + ", med_rhs)    , "Resting LVOT gradient" ,
    "septal_thickness_z" , "septal_thickness_z" , paste0("septal_thickness_z + age + Sex + BMI + ", med_rhs) , "LV septal thickness" ,
    "vo2_friend2_z"      , "vo2_friend2_z"      , paste0("vo2_friend2_z + age + Sex + BMI + ", med_rhs)      , "Peak V̇O2 (per SD lower)" ,
    "vevco2_slope_z"     , "vevco2_slope_z"     , paste0("vevco2_slope_z + age + Sex + BMI + ", med_rhs)     , "V̇E/V̇CO2 slope"
  )

  for (i in seq_len(nrow(sensitivity_specs))) {
    spec <- sensitivity_specs[i, ]
    dat_sub <- data %>% filter(!is.na(.data[[spec$filter_var]]))
    if (
      nrow(dat_sub) < 20 || sum(dat_sub[[event_var]] == 1, na.rm = TRUE) < 5
    ) {
      next
    }
    fit <- tryCatch(
      coxph(
        as.formula(paste0(
          "Surv(",
          time_var,
          ", ",
          event_var,
          ") ~ ",
          spec$rhs
        )),
        data = dat_sub
      ),
      error = function(e) NULL
    )
    out_rows[[length(out_rows) + 1]] <- tidy_sensitivity_cox(
      fit,
      spec$term,
      spec$label,
      scenario_label,
      nrow(dat_sub),
      sum(dat_sub[[event_var]] == 1, na.rm = TRUE)
    )
  }

  bind_rows(out_rows)
}

# Medication-adjusted sensitivity: the reported version models medication
# missingness with three-level (No / Yes / Unknown) covariates rather than
# requiring recorded status (owner decision, card 1.5 of
# 8_Docs/DECISION_CARDS_2026-09.md). Requiring it halved every denominator
# without moving the hazard ratios, so the complete-case comparison is kept
# below as an explicit supplemental contrast rather than as the headline.
surv_bb_df <- surv_sensitivity_df %>%
  filter(med_data_available == 1, !is.na(bb_any), !is.na(ndhp_ccb_any))

surv_med_indicator_build <- add_missing_indicator_covariates(
  surv_sensitivity_df,
  binary_vars = c("bb_any", "ndhp_ccb_any")
)

figure5_sensitivity_scenarios <- list(
  list(
    label = "Main composite + expanded covariates",
    time_var = "follow_up_yrs",
    event_var = "hf_composite"
  ),
  list(
    label = "HF/death without transplant + expanded covariates",
    time_var = "follow_up_no_transplant_yrs",
    event_var = "hf_or_death_no_transplant_event"
  ),
  list(
    label = "HF/death censored at intervention + expanded covariates",
    time_var = "follow_up_intervention_censored_yrs",
    event_var = "hf_or_death_intervention_censored_event"
  )
)

run_figure5_sensitivity_scenarios <- function(data, med_terms) {
  map_dfr(figure5_sensitivity_scenarios, function(scenario) {
    run_figure5_sensitivity_suite(
      data,
      scenario$label,
      scenario$time_var,
      scenario$event_var,
      med_terms = med_terms
    )
  })
}

figure5_sensitivity_hr_table <- run_figure5_sensitivity_scenarios(
  surv_med_indicator_build$data,
  surv_med_indicator_build$covariates
)

figure5_sensitivity_complete_case_rows <- run_figure5_sensitivity_scenarios(
  surv_bb_df,
  c("bb_any", "ndhp_ccb_any")
)

if (nrow(figure5_sensitivity_hr_table) > 0) {
  figure5_sensitivity_hr_table <- figure5_sensitivity_hr_table %>%
    mutate(
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
      `P value` = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
    )

  write.csv(
    figure5_sensitivity_hr_table,
    file.path(OUT_DIR, "Table_Sensitivity_Figure5_BetaBlockerAndIntervention.csv"),
    row.names = FALSE
  )

  figure5_sensitivity_hr_table %>%
    select(Scenario, Exposure, N, Events, `HR (95% CI)`, `P value`) %>%
    flextable() %>%
    bold(part = "header") %>%
    align(j = c("Scenario", "Exposure"), align = "left", part = "all") %>%
    theme_booktabs() %>%
    autofit() %>%
    save_as_docx(
      path = file.path(OUT_DIR, "Table_Sensitivity_Figure5_BetaBlockerAndIntervention.docx")
    )

  # ── Reported versus complete-case medication adjustment ────────────────────
  # Side-by-side contrast of the two handlings. The reported suite models
  # medication missingness; the complete-case suite is what requiring recorded
  # status produces on the same models (card 1.5).
  figure5_missing_indicator_comparison <- bind_rows(
    figure5_sensitivity_complete_case_rows %>%
      select(Scenario, Exposure, N, Events, HR, CI_low, CI_high, p_value) %>%
      mutate(`Covariate handling` = "Complete case (medication recorded)"),
    figure5_sensitivity_hr_table %>%
      select(Scenario, Exposure, N, Events, HR, CI_low, CI_high, p_value) %>%
      mutate(
        `Covariate handling` = "Three-level No / Yes / Unknown (reported)"
      )
  ) %>%
    arrange(Scenario, Exposure, `Covariate handling`) %>%
    mutate(
      `HR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", HR, CI_low, CI_high),
      `P value` = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))
    )

  # Modelling the missingness can only add patients, never remove them, so
  # every indicator-handled fit must be at least as large as its complete-case
  # counterpart. A violation would mean the two arms differ in more than the
  # covariate handling.
  figure5_missing_indicator_nesting <- figure5_missing_indicator_comparison %>%
    select(Scenario, Exposure, `Covariate handling`, N) %>%
    tidyr::pivot_wider(names_from = `Covariate handling`, values_from = N) %>%
    filter(
      !is.na(`Complete case (medication recorded)`),
      !is.na(`Three-level No / Yes / Unknown (reported)`)
    )
  stopifnot(
    nrow(figure5_missing_indicator_nesting) > 0,
    all(
      figure5_missing_indicator_nesting$`Three-level No / Yes / Unknown` >=
        figure5_missing_indicator_nesting$`Complete case (medication recorded)`
    )
  )

  write.csv(
    figure5_missing_indicator_comparison %>%
      select(
        Scenario, Exposure, `Covariate handling`, N, Events,
        `HR (95% CI)`, `P value`
      ),
    file.path(OUT_DIR, "Table_Sensitivity_MissingIndicator_Cox.csv"),
    row.names = FALSE
  )
  print(
    as.data.frame(
      figure5_missing_indicator_comparison %>%
        filter(Scenario == "Main composite + expanded covariates") %>%
        select(Exposure, `Covariate handling`, N, Events, `HR (95% CI)`, `P value`)
    ),
    row.names = FALSE
  )

  figure5_sensitivity_plot_df <- figure5_sensitivity_hr_table %>%
    mutate(
      Scenario = factor(
        Scenario,
        levels = c(
          "Main composite + expanded covariates",
          "HF/death without transplant + expanded covariates",
          "HF/death censored at intervention + expanded covariates"
        )
      ),
      Exposure = factor(Exposure, levels = rev(unique(Exposure)))
    )

  p_fig5_sensitivity <- ggplot(
    figure5_sensitivity_plot_df,
    aes(x = HR, y = Exposure)
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = fig5_palette["refline"],
      linewidth = 0.45
    ) +
    geom_errorbarh(
      aes(xmin = CI_low, xmax = CI_high),
      height = 0.12,
      linewidth = 0.65,
      color = fig5_palette["forest"]
    ) +
    geom_point(
      size = 2.2,
      shape = 21,
      stroke = 0.35,
      fill = fig5_palette["highlight"],
      color = fig5_palette["forest"]
    ) +
    facet_wrap(~Scenario, ncol = 1, scales = "free_y") +
    scale_x_continuous(breaks = c(0.5, 0.75, 1, 1.5, 2, 3, 4)) +
    labs(x = "Hazard Ratio", y = NULL) +
    theme_jacc(base_size = 8.3) +
    theme(
      strip.text = element_text(face = "bold", size = 7.6),
      axis.text.y = element_text(size = 7.0, color = fig5_palette["text"]),
      axis.text.x = element_text(color = fig5_palette["text"]),
      axis.line = element_line(color = fig5_palette["text"], linewidth = 0.3)
    )

  show_and_save_jacc(
    p_fig5_sensitivity,
    file.path(OUT_DIR, "Figure5_HeartFailure_Outcomes_Sensitivity.pdf"),
    w = 7.0,
    h = 6.0
  )
}

figure5_primary_compare <- if (exists("forest_rows") && nrow(forest_rows) > 0) {
  forest_rows %>%
    transmute(
      Scenario = "Primary composite (expanded covariates)",
      Exposure = as.character(exposure),
      `HR (95% CI)` = sprintf(
        "%.2f (%.2f-%.2f)",
        estimate,
        conf.low,
        conf.high
      ),
      `P value` = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
} else {
  tibble()
}

figure5_sensitivity_comparison <- bind_rows(
  figure5_primary_compare,
  figure5_sensitivity_hr_table %>%
    select(Scenario, Exposure, `HR (95% CI)`, `P value`)
)

if (nrow(figure5_sensitivity_comparison) > 0) {
  write.csv(
    figure5_sensitivity_comparison,
    file.path(OUT_DIR, "Table_Sensitivity_Figure5_Comparison.csv"),
    row.names = FALSE
  )
}

# (legacy chunk: figure5_validation)
# ── Cohort sizes ──────────────────────────────────────────────────────────────
hf_follow_n_live <- nrow(surv_df)
hf_event_n_live <- sum(surv_df$hf_composite == 1, na.rm = TRUE)
hf_model_n <- nrow(surv_model_df)
hf_model_events <- sum(surv_model_df$hf_composite == 1, na.rm = TRUE)

# Helper: extract HR, 95% CI, and P from a fitted coxph object for a named term
extract_cox <- function(fit, term_name) {
  if (is.null(fit)) {
    return(list(hr = NA, lo = NA, hi = NA, p = NA, n = NA, events = NA))
  }
  td <- tryCatch(
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE),
    error = function(e) NULL
  )
  if (is.null(td)) {
    return(list(hr = NA, lo = NA, hi = NA, p = NA, n = NA, events = NA))
  }
  row <- td[td$term == term_name, ]
  if (nrow(row) == 0) {
    return(list(hr = NA, lo = NA, hi = NA, p = NA, n = NA, events = NA))
  }
  nd <- tryCatch(nrow(fit$model), error = function(e) NA)
  ne <- tryCatch(sum(fit$y[, 2]), error = function(e) NA)
  list(
    hr = row$estimate,
    lo = row$conf.low,
    hi = row$conf.high,
    p = row$p.value,
    n = nd,
    events = ne
  )
}

fmt_hr <- function(x) sprintf("HR %.2f (%.2f\u2013%.2f)", x$hr, x$lo, x$hi)
fmt_p <- function(x) {
  if (is.na(x$p)) {
    return("P = NA")
  }
  if (x$p < 0.001) {
    sprintf("P = %s", format_p_scientific(x$p))
  } else {
    sprintf("P = %.3f", x$p)
  }
}
fmt_row <- function(label, x) {
  flag <- if (!is.na(x$p) && x$p < 0.05) " \u2713" else ""
  sprintf(
    "- **%s** | N = %d, events = %d | %s, %s%s\n",
    label,
    x$n,
    x$events,
    fmt_hr(x),
    fmt_p(x),
    flag
  )
}

r_lavi <- extract_cox(cox_lavi, "la_vol_index_z")
r_ee <- extract_cox(cox_ee, "e_e_ave_z")
r_trv <- extract_cox(cox_trv, "tr_max_vel_z")
r_lvot <- extract_cox(cox_lvot, "lvot_gradient_z")
r_septal <- extract_cox(cox_septal, "septal_thickness_z")
r_vo2 <- extract_cox(cox_vo2, "vo2_friend2_z")
r_vevco2 <- extract_cox(cox_vevco2, "vevco2_slope_z")

cat("::: {.manuscript-check}\n")
cat("**Figure 5 Validation — Cohort Sizes**\n\n")
cat(sprintf(
  "- Full outcomes cohort (KM): **%d patients**, **%d events** (%.1f%%)\n",
  hf_follow_n_live,
  hf_event_n_live,
  100 * hf_event_n_live / hf_follow_n_live
))
cat(sprintf(
  "- Cox model cohort (complete covariates): **%d patients**, **%d events** (%.1f%%)\n",
  hf_model_n,
  hf_model_events,
  100 * hf_model_events / hf_model_n
))
cat("\n")

cat("**Figure 5 Validation — Cox Model HRs (per 1 SD; \u2713 = P < 0.05)**\n\n")
cat(fmt_row("LAVi Z-score", r_lavi))
cat(fmt_row("E/e\u2019 Z-score", r_ee))
cat(fmt_row("TRVmax Z-score", r_trv))
cat(fmt_row("Resting LVOT gradient Z-score", r_lvot))
cat(fmt_row("LV septal thickness Z-score", r_septal))
cat(fmt_row("Peak V\u0307O\u2082 Z-score (per SD decrease)", r_vo2))
cat(fmt_row("V\u0307E/V\u0307CO\u2082 slope Z-score", r_vevco2))
cat(":::\n")


message("Stage 07 complete.")

message("Stage 07 complete.")
