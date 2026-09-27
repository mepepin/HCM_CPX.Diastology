# ── Stage 05: cross-sectional models ──────────────────────────────────────────
# In  : 5_Data/2_Processed/{baseline_df, cohort_parameters}.rds
# Out : 6_Results/3_CrossSectional/
#         Figure2_ASE2025_FillingPressure.pdf          (main Figure 2)
#         Figure3_Nonlinear_DiastolicIndices_RCS.pdf   (main Figure 3)
#         Figure_S1_AllDiastolic_RCS.pdf               (Figure S2)
#         Table_Sensitivity_SameSample_Adjustment.csv  (Table S1)
#         Table2_Figure3_RCS_Summary.csv/.docx         (spline P and q values)
#         Stats_CrossSectional_LRT_*.csv               (partial F, bootstrap q)
#         Table_Manuscript_CrossSectional_Values.csv   (in-text values)
#         Table_Consistency_CommonCohort_RCS.csv, Table_Sensitivity_AdultsOnly_RCS.csv,
#         Table_Sensitivity_MultipleImputation_RCS.csv (Methods sensitivities)
#
# Each diastolic index is modelled against each fitness outcome with 4-knot
# restricted cubic splines adjusted for age, sex, and BMI; the per-SD linear
# estimates quoted in the Results come from the same-sample adjustment table.

source("R/00_setup.R")
stage_banner("Stage 05: cross-sectional models")

OUT_DIR <- results_dir("3_CrossSectional")
AUDIT_DIR <- results_dir("1_Audit")

# ── Inputs from earlier stages ────────────────────────────────────────────────
baseline_df <- load_data("baseline_df", "processed")
medications <- load_data("medications", "processed")
analytic_cohort_mode <- load_data("cohort_parameters", "processed")$analytic_cohort_mode

# ══ Figure 2: HCM Primary-Parameter Phenotypes ════════════════════════════════
# (legacy chunk: figure2_phenotypes)
# ── Panel A: UpSet plot of primary parameter combinations ────────────────────
df_combo <- baseline_df %>%
  filter(cross_sectional_eligible) %>%
  mutate(
    eprime_age_lln = case_when(
      age >= 20 & age < 40 ~ 9.0,
      age >= 40 & age <= 65 ~ 7.0,
      age > 65 ~ 6.5,
      TRUE ~ NA_real_
    ),
    eprime_age_margin = e_prime_ave - eprime_age_lln,
    hcm_combo_code = paste0(
      as.integer(e_e_ave > 14),
      as.integer(la_vol_index > 34)
    ),
    hcm_combo_class = factor(
      case_when(
        hcm_combo_code == "00" ~ "Neither above threshold",
        hcm_combo_code == "10" ~ "E/e' only",
        hcm_combo_code == "01" ~ "LAVi only",
        hcm_combo_code == "11" ~ "E/e' + LAVi",
        TRUE ~ NA_character_
      ),
      levels = c(
        "Neither above threshold",
        "E/e' only",
        "LAVi only",
        "E/e' + LAVi"
      )
    ),
    hcm_combo_abnormal_group = factor(
      stringr::str_count(hcm_combo_code, "1"),
      levels = 0:3,
      labels = c("0", "1", "2", "3")
    )
  )
combo_palette <- c(
  "0" = "#E5E8EC",
  "1" = alpha("#C9A2A0", 0.72),
  "2" = "#B97E79",
  "3" = "#8E5C59"
)

combo_palette_fig3 <- c(
  "0" = "#CBD3DA",
  "1" = alpha("#C9A2A0", 0.65),
  "2" = "#B97E79",
  "3" = "#8E5C59"
)

combo_levels <- levels(df_combo$hcm_combo_class)
combo_parameter_names <- c("E/e'", "LAVi")
combo_signature <- tibble(
  hcm_combo_code = c("00", "10", "01", "11"),
  hcm_combo_class = factor(
    c(
      "Neither above threshold",
      "E/e' only",
      "LAVi only",
      "E/e' + LAVi"
    ),
    levels = combo_levels
  )
) %>%
  mutate(
    hcm_combo_abnormal_group = factor(
      stringr::str_count(hcm_combo_code, "1"),
      levels = 0:3,
      labels = names(combo_palette)
    ),
    `E/e'` = ifelse(substr(hcm_combo_code, 1, 1) == "1", "+", "-"),
    LAVi = ifelse(substr(hcm_combo_code, 2, 2) == "1", "+", "-"),
    x_id = row_number()
  )

param_levels <- c("LAVi", "E/e'")
param_positions <- c("LAVi" = 1.20, "E/e'" = 2.28)
x_limits_upset <- c(0.5, nrow(combo_signature) + 0.5)

intersection_counts <- combo_signature %>%
  select(x_id, hcm_combo_class, hcm_combo_abnormal_group) %>%
  left_join(df_combo %>% count(hcm_combo_class), by = "hcm_combo_class") %>%
  mutate(n = coalesce(n, 0L))

upset_matrix <- combo_signature %>%
  pivot_longer(
    cols = all_of(combo_parameter_names),
    names_to = "Parameter",
    values_to = "Status"
  ) %>%
  mutate(
    Parameter = factor(Parameter, levels = param_levels),
    param_id = unname(param_positions[as.character(Parameter)]),
    dot_fill = ifelse(Status == "+", "active", "inactive")
  )

upset_connectors <- upset_matrix %>%
  filter(Status == "+") %>%
  group_by(x_id, hcm_combo_class, hcm_combo_abnormal_group) %>%
  summarize(
    ymin = min(param_id),
    ymax = max(param_id),
    n_active = n(),
    .groups = "drop"
  ) %>%
  filter(n_active >= 2)

p_upset_top <- intersection_counts %>%
  ggplot(aes(x = x_id, y = n, fill = hcm_combo_abnormal_group)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.2) +
  geom_text(aes(label = n), vjust = -0.25, size = 3.1) +
  scale_fill_manual(values = combo_palette, drop = FALSE) +
  scale_x_continuous(
    limits = x_limits_upset,
    breaks = combo_signature$x_id,
    labels = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = "Patients") +
  theme_jacc() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 6, 4, 4)
  )

p_upset_matrix <- ggplot(upset_matrix, aes(x = x_id, y = param_id)) +
  geom_segment(
    data = upset_connectors,
    aes(x = x_id, xend = x_id, y = ymin, yend = ymax),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_point(aes(color = dot_fill), shape = 16, size = 2.5) +
  scale_color_manual(
    values = c(active = "black", inactive = "grey75"),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = x_limits_upset,
    breaks = combo_signature$x_id,
    labels = NULL
  ) +
  scale_y_continuous(
    breaks = unname(param_positions[param_levels]),
    labels = param_levels,
    limits = c(0.70, 2.78),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL) +
  theme_jacc() +
  theme(
    legend.position = "none",
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6.9, margin = margin(r = 4)),
    panel.grid = element_blank(),
    plot.margin = margin(0, 6, 0, 4)
  )

p_bar <- wrap_elements(
  full = patchwork::patchworkGrob(
    p_upset_top / p_upset_matrix + plot_layout(heights = c(4.25, 1.05))
  )
)

p_upset_header <- patchwork::wrap_elements(
  full = grid::grobTree(
    grid::textGrob(
      "(B) E/e' and LAVi Threshold Combinations",
      x = grid::unit(0, "npc"),
      y = grid::unit(1, "npc") - grid::unit(2, "pt"),
      just = c("left", "top"),
      gp = grid::gpar(fontface = "bold", fontsize = 11, col = "#23313F")
    ),
    grid::textGrob(
      sprintf(
        "%s cohort: n = %d",
        if_else(
          analytic_cohort_mode == "complete_three",
          "Complete E/e', LAVi, and TRVmax",
          "Complete E/e' and LAVi"
        ),
        nrow(df_combo)
      ),
      x = grid::unit(0, "npc"),
      y = grid::unit(1, "npc") - grid::unit(14, "pt"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 8.4, col = "grey40")
    )
  )
)

# ── Panel B: HCM diastolic abnormality prevalence ────────────────────────────
fig2_prevalence_thresholds <- tribble(
  ~Group                  , ~variable        , ~Parameter       , ~direction , ~cutoff , ~ThresholdLabel        ,
  "Primary parameters"    , "e_e_ave"        , "E/e'"           , ">"        ,  14     , "> 14"                  ,
  "Primary parameters"    , "la_vol_index"   , "LAVi"           , ">"        ,  34     , "> 34 mL/m\u00b2"       ,
  "Primary parameters"    , "tr_max_vel"     , "TRVmax"         , ">"        , 280     , "> 2.8 m/s"             ,
  "Supportive parameters" , "eprime_age_margin", "Age-calibrated average e'", "<", 0, "below age-specific reference limit",
  "Supportive parameters" , "mv_e_a"         , "MV E/A"         , ">="       ,   2     , "\u2265 2.0"           ,
  "Supportive parameters" , "mv_dec_time"    , "MV Decel. time" , "<"        , 150     , "< 150 ms"             ,
  "Supportive parameters" , "pv_sd_ratio"    , "PV S/D Ratio"   , "<"        ,   1     , "< 1.0"
)

fig2_eval_abnormal <- function(x, direction, cutoff) {
  case_when(
    direction == ">" ~ x > cutoff,
    direction == ">=" ~ x >= cutoff,
    direction == "<" ~ x < cutoff,
    direction == "<=" ~ x <= cutoff,
    TRUE ~ NA
  )
}

fig2_prepare_values <- function(x, variable) {
  x_num <- as_num(x)
  if (variable == "mv_dec_time") {
    med_val <- stats::median(x_num, na.rm = TRUE)
    if (is.finite(med_val) && med_val < 10) {
      x_num <- x_num * 1000
    }
  }
  x_num
}

fig2_prevalence_summary <- bind_rows(lapply(
  seq_len(nrow(fig2_prevalence_thresholds)),
  function(i) {
    spec <- fig2_prevalence_thresholds[i, ]
    values <- fig2_prepare_values(df_combo[[spec$variable]], spec$variable)
    available_n <- sum(!is.na(values))
    abnormal <- fig2_eval_abnormal(values, spec$direction, spec$cutoff)
    abnormal_n <- sum(abnormal, na.rm = TRUE)
    abnormal_pct <- if (available_n > 0) {
      100 * abnormal_n / available_n
    } else {
      NA_real_
    }

    tibble(
      Group = spec$Group,
      Parameter = spec$Parameter,
      ThresholdLabel = spec$ThresholdLabel,
      available_n = available_n,
      abnormal_n = abnormal_n,
      abnormal_pct = abnormal_pct,
      label = sprintf("%d/%d (%.1f%%)", abnormal_n, available_n, abnormal_pct)
    )
  }
)) %>%
  mutate(
    Group = factor(
      Group,
      levels = c("Primary parameters", "Supportive parameters")
    ),
    Parameter = factor(
      Parameter,
      levels = c(
        "E/e'",
        "LAVi",
        "TRVmax",
        "Age-calibrated average e'",
        "MV E/A",
        "MV Decel. time",
        "PV S/D Ratio"
      )
    )
  )

fig2_prev_palette <- c(
  "Primary parameters" = "#8E5C59",
  "Supportive parameters" = "#B8C1C9"
)
fig2_prev_label_palette <- c(
  "Primary parameters" = "#8E5C59",
  "Supportive parameters" = "#66717A"
)

fig2_prev_axis_max <- max(
  70,
  ceiling((max(fig2_prevalence_summary$abnormal_pct, na.rm = TRUE) + 10) / 10) *
    10
)
fig2_prev_breaks <- seq(0, fig2_prev_axis_max, by = 20)
fig2_prev_plot_max <- fig2_prev_axis_max + 18

fig2_prev_order <- c(
  "Age-calibrated average e'",
  "LAVi",
  "E/e'",
  "PV S/D Ratio",
  "TRVmax",
  "MV Decel. time",
  "MV E/A"
)

fig2_prev_axis_labels <- c(
  "E/e'" = "E/e'\n> 14",
  "LAVi" = "LAVi\n> 34 mL/m\u00b2",
  "TRVmax" = "TRVmax\n> 2.8 m/s",
  "Age-calibrated average e'" = "Age-calibrated average e'\nbelow age-specific limit",
  "MV E/A" = "MV E/A\n\u2265 2.0",
  "MV Decel. time" = "MV Decel. time\n< 150 ms",
  "PV S/D Ratio" = "PV S/D Ratio\n< 1.0"
)

fig2_prevalence_plot_df <- fig2_prevalence_summary %>%
  mutate(
    Parameter = factor(as.character(Parameter), levels = rev(fig2_prev_order)),
    ParameterLabel = factor(
      fig2_prev_axis_labels[as.character(Parameter)],
      levels = rev(unname(fig2_prev_axis_labels[fig2_prev_order]))
    ),
    pct_label_x = pmin(abnormal_pct + 2.5, fig2_prev_axis_max - 4.5),
    pct_label = sprintf("%.1f%%", abnormal_pct),
    count_label_x = fig2_prev_axis_max + 9,
    count_label = sprintf("%d/%d", abnormal_n, available_n)
  )

p_prev_header <- patchwork::wrap_elements(
  full = grid::grobTree(
    grid::textGrob(
      "(A) LV Diastolic Parameter Prevalence",
      x = grid::unit(0, "npc"),
      y = grid::unit(1, "npc") - grid::unit(2, "pt"),
      just = c("left", "top"),
      gp = grid::gpar(
        fontface = "bold",
        fontsize = 11,
        col = "#23313F",
        lineheight = 0.95
      )
    )
  )
)

p_prev <- ggplot(
  fig2_prevalence_plot_df,
  aes(x = ParameterLabel, y = abnormal_pct, fill = Group)
) +
  geom_col(width = 0.68, color = NA) +
  geom_text(
    aes(y = pct_label_x, label = pct_label, color = Group),
    hjust = 0,
    size = 2.7,
    fontface = "bold",
    show.legend = FALSE
  ) +
  geom_text(
    aes(y = count_label_x, label = count_label),
    hjust = 0,
    size = 2.35,
    color = "#56616B",
    show.legend = FALSE
  ) +
  geom_hline(
    yintercept = 0,
    inherit.aes = FALSE,
    linewidth = 0.3,
    color = "#3E4650"
  ) +
  scale_fill_manual(values = fig2_prev_palette, guide = "none") +
  scale_color_manual(values = fig2_prev_label_palette, guide = "none") +
  scale_y_continuous(
    limits = c(0, fig2_prev_plot_max),
    breaks = fig2_prev_breaks,
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = "Threshold prevalence") +
  coord_flip(clip = "off") +
  theme_jacc(base_size = 8.1) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 6.9, margin = margin(r = 4)),
    axis.text.x = element_text(
      size = 6.6,
      color = "#2D3440",
      margin = margin(r = 7)
    ),
    axis.text.y = element_text(
      size = 6.3,
      lineheight = 0.92,
      color = "#2D3440"
    ),
    axis.line.y = element_line(color = "#3E4650", linewidth = 0.3),
    axis.ticks.y = element_line(color = "#3E4650", linewidth = 0.25),
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E3E7EB", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1, 10, 3, 4)
  )

p_prev_panel <- wrap_elements(full = p_prev)

# ── Combination-phenotype Kruskal-Wallis tests ───────────────────────────────
# Computed once here and reused by the panel headers below and by the exported
# manuscript values, so the figure labels and the Results text cannot diverge.
# These were previously computed inline as figure labels only, leaving the
# manuscript without a traceable source.
combo_kw_vo2 <- kruskal.test(VO2_FRIEND2_PP ~ hcm_combo_class, data = df_combo)
combo_kw_vevco2 <- kruskal.test(VeVco2_slope ~ hcm_combo_class, data = df_combo)

# ── Panel C: Peak VO2 by combination phenotype ───────────────────────────────
p_box_kw_p <- format.pval(combo_kw_vo2$p.value, digits = 2)

p_box_header <- patchwork::wrap_elements(
  full = grid::grobTree(
    grid::textGrob(
      expression(bold("(C) Peak " * dot(V) * O[2] * " by Combination Phenotype")),
      x = grid::unit(0, "npc"),
      y = grid::unit(1, "npc") - grid::unit(2, "pt"),
      just = c("left", "top"),
      gp = grid::gpar(
        fontface = "bold",
        fontfamily = "Helvetica",
        fontsize = 11,
        col = "#23313F"
      )
    )
  )
)

p_box_top <- df_combo %>%
  filter(!is.na(VO2_FRIEND2_PP)) %>%
  ggplot(aes(
    x = hcm_combo_class,
    y = VO2_FRIEND2_PP,
    fill = hcm_combo_abnormal_group
  )) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    linewidth = 0.3,
    color = alpha("#56606A", 0.8)
  ) +
  geom_point(position = position_jitter(width = 0.15, seed = JITTER_SEED), alpha = 0.18, size = 0.8, color = "#5E6872") +
  scale_fill_manual(values = combo_palette_fig3, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = expression("Peak " * dot(V) * O[2] * " (FRIEND 2.0 %)" )
  ) +
  theme_jacc() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(
      family = "Helvetica",
      margin = margin(r = 7)
    ),
    plot.margin = margin(4, 6, 0, 4)
  )

p_box_sig <- ggplot(upset_matrix, aes(x = hcm_combo_class, y = param_id)) +
  geom_segment(
    data = upset_connectors,
    aes(
      x = hcm_combo_class,
      xend = hcm_combo_class,
      y = ymin,
      yend = ymax
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_point(aes(color = dot_fill), shape = 16, size = 2.5) +
  scale_color_manual(
    values = c(active = "black", inactive = "grey75"),
    drop = FALSE
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = unname(param_positions[param_levels]),
    labels = param_levels,
    limits = c(0.70, 2.78),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL) +
  theme_jacc() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 6.9, margin = margin(r = 4)),
    panel.grid = element_blank(),
    plot.margin = margin(0, 6, 4, 4)
  )

p_box <- p_box_top / p_box_sig + plot_layout(heights = c(4.2, 1.35))
p_box_panel <- wrap_elements(full = patchwork::patchworkGrob(p_box))

# ── Panel D: VE/VCO2 by combination phenotype ────────────────────────────────
p_vevco2_kw_p <- format.pval(combo_kw_vevco2$p.value, digits = 2)

p_vevco2_header <- patchwork::wrap_elements(
  full = grid::grobTree(
    grid::textGrob(
      expression(bold("(D) " * dot(V) * E * "/" * dot(V) * C * O[2] * " by Combination Phenotype")),
      x = grid::unit(0, "npc"),
      y = grid::unit(1, "npc") - grid::unit(2, "pt"),
      just = c("left", "top"),
      gp = grid::gpar(
        fontface = "bold",
        fontfamily = "Helvetica",
        fontsize = 11,
        col = "#23313F"
      )
    )
  )
)

p_vevco2_top <- df_combo %>%
  filter(!is.na(VeVco2_slope)) %>%
  ggplot(aes(
    x = hcm_combo_class,
    y = VeVco2_slope,
    fill = hcm_combo_abnormal_group
  )) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    linewidth = 0.3,
    color = alpha("#56606A", 0.8)
  ) +
  geom_point(position = position_jitter(width = 0.15, seed = JITTER_SEED), alpha = 0.18, size = 0.8, color = "#5E6872") +
  scale_fill_manual(values = combo_palette_fig3, drop = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(
    x = NULL,
    y = expression(dot(V) * E * "/" * dot(V) * C * O[2] * " Slope")
  ) +
  theme_jacc() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(
      family = "Helvetica",
      margin = margin(r = 7)
    ),
    plot.margin = margin(4, 6, 0, 4)
  )

p_vevco2_sig <- ggplot(upset_matrix, aes(x = hcm_combo_class, y = param_id)) +
  geom_segment(
    data = upset_connectors,
    aes(
      x = hcm_combo_class,
      xend = hcm_combo_class,
      y = ymin,
      yend = ymax
    ),
    inherit.aes = FALSE,
    linewidth = 0.55,
    color = "black"
  ) +
  geom_point(aes(color = dot_fill), shape = 16, size = 2.5) +
  scale_color_manual(
    values = c(active = "black", inactive = "grey75"),
    drop = FALSE
  ) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    breaks = unname(param_positions[param_levels]),
    labels = param_levels,
    limits = c(0.70, 2.78),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = NULL, y = NULL) +
  theme_jacc() +
  theme(
    legend.position = "none",
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(
      size = 6.9,
      margin = margin(r = 4),
      color = "transparent"
    ),
    panel.grid = element_blank(),
    plot.margin = margin(0, 6, 4, 4)
  )

p_vevco2 <- p_vevco2_top / p_vevco2_sig + plot_layout(heights = c(4.2, 1.35))
p_vevco2_panel <- wrap_elements(full = patchwork::patchworkGrob(p_vevco2))

# ── Compose Figure 2 ─────────────────────────────────────────────────────────
fig_ase <- ((p_prev_header / p_prev_panel + plot_layout(heights = c(0.10, 1))) |
  (p_upset_header / p_bar + plot_layout(heights = c(0.14, 1)))) /
  ((p_box_header / p_box_panel + plot_layout(heights = c(0.14, 1))) |
    (p_vevco2_header / p_vevco2_panel + plot_layout(heights = c(0.14, 1)))) +
  plot_layout(heights = c(0.98, 0.92), widths = c(0.92, 1.08)) &
  theme(plot.tag = element_text(face = "bold", size = 11))

figure2_main_width <- 7.2
figure2_main_height <- 8.1
figure2_preprint_width <- 7.2
figure2_preprint_height <- 6.20
figure2_preprint_legend_fontsize <- 9

show_and_save_jacc(
  fig_ase,
  file.path(OUT_DIR, "Figure2_ASE2025_FillingPressure.pdf"),
  w = figure2_main_width,
  h = figure2_main_height
)

legend_fig2_pdf <- paste0(
  "**Figure 2. Baseline HCM LV Diastolic Parameters and Cardiopulmonary Performance.** ",
  if (analytic_cohort_mode == "complete_three") {
    sprintf(
      "Complete E/e\u2019, LAVi, and TRVmax were required (n=%d). ",
      nrow(df_combo)
    )
  } else {
    sprintf(
      "Complete E/e\u2019 and LAVi were required (n=%d); TRVmax was available in %d. ",
      nrow(df_combo),
      sum(!is.na(df_combo$tr_max_vel))
    )
  },
  paste0(
    "This figure displays combinations of E/e\u2019 and LAVi threshold ",
    "crossings and therefore requires both parameters; the cross-sectional ",
    "models analyse each parameter in its own available-case sample. "
  ),
  "**(A)** Parameter prevalence using observed denominators. ",
  "**(B)** Combinations of E/e\u2019 and LAVi threshold crossings. ",
  "**(C)** ",
  label_peak_vo2_md,
  " (FRIEND 2.0 % predicted) stratified by combinatorial ",
  "filling-pressure phenotype, with Kruskal-Wallis P value. ",
  "**(D)** ",
  label_vevco2_md,
  " stratified by the same combinatorial phenotype, ",
  "with Kruskal-Wallis P value."
)

save_jacc_with_embedded_legend(
  fig_ase,
  file.path(OUT_DIR, "Figure2_ASE2025_FillingPressure.pdf"),
  legend_text = legend_fig2_pdf,
  w = figure2_main_width,
  h = figure2_main_height
)
save_jacc_with_embedded_legend(
  fig_ase,
  file.path(OUT_DIR, "Figure2_ASE2025_FillingPressure.pdf"),
  legend_text = legend_fig2_pdf,
  w = figure2_main_width,
  h = figure2_main_height
)
save_jacc(
  fig_ase,
  file.path(OUT_DIR, "Figure2_ASE2025_FillingPressure_Preprint.pdf"),
  w = figure2_preprint_width,
  h = figure2_preprint_height
)
save_jacc(
  fig_ase,
  file.path(OUT_DIR, "Figure2_ASE2025_FillingPressure_Preprint.pdf"),
  w = figure2_preprint_width,
  h = figure2_preprint_height
)

# ── Manuscript validation ────────────────────────────────────────────────────
cat("**Figure 2 Validation**\n\n")
fig2_combo_counts <- table(df_combo$hcm_combo_class)
cat(sprintf(
  "- Neither E/e' nor LAVi above threshold: **%d**\n",
  fig2_combo_counts["Neither above threshold"]
))
cat(sprintf("- E/e' only: **%d**\n", fig2_combo_counts["E/e' only"]))
cat(sprintf("- LAVi only: **%d**\n", fig2_combo_counts["LAVi only"]))
cat(sprintf("- E/e' + LAVi: **%d**\n", fig2_combo_counts["E/e' + LAVi"]))
cat(sprintf(
  "\n- %s KW p = **%s**\n",
  label_peak_vo2_md,
  p_box_kw_p
))
cat(sprintf(
  "- %s KW p = **%s**\n",
  label_vevco2_md,
  p_vevco2_kw_p
))

# (legacy chunk: figure2_stratified_table)
# ── table: baseline characteristics by diastolic phenotype combo ──────────────
tbl_ase <- df_combo %>%
  select(
    hcm_combo_class,
    age,
    Sex,
    BMI,
    LBMI,
    HCM_Phenotype,
    VO2_FRIEND2_PP,
    VO2_WASSERMAN_PP,
    VeVco2_slope,
    HRR,
    lvot_max_gradient,
    lv_septal_thickness,
    e_prime_ave,
    med_peak_e_vel,
    lat_peak_e_vel,
    e_e_ave,
    la_vol_index,
    tr_max_vel,
    mv_e_a,
    mv_dec_time,
    ef_modsp4
  ) %>%
  tbl_summary(
    by = hcm_combo_class,
    label = list(
      age ~ "Age, years",
      Sex ~ "Sex",
      BMI ~ "BMI, kg/m\u00B2",
      LBMI ~ "LBMI, kg/m\u00B2",
      HCM_Phenotype ~ "HCM Phenotype",
      VO2_FRIEND2_PP ~ "Peak V\u0307O\u2082 (FRIEND 2.0 %pred)",
      VO2_WASSERMAN_PP ~ "Peak V\u0307O\u2082 (Wasserman %pred)",
      VeVco2_slope ~ "V\u0307E/V\u0307CO\u2082 Slope",
      HRR ~ "HR Recovery (1 min)",
      lvot_max_gradient ~ "Resting LVOT Gradient, mm Hg",
      lv_septal_thickness ~ "LV Septal Thickness, cm",
      e_prime_ave ~ "Average e', cm/s",
      med_peak_e_vel ~ "Septal e', cm/s",
      lat_peak_e_vel ~ "Lateral e', cm/s",
      e_e_ave ~ "E/e' (average)",
      la_vol_index ~ "LAVi, mL/m\u00B2",
      tr_max_vel ~ "TRVmax, cm/s",
      mv_e_a ~ "E/A Ratio",
      mv_dec_time ~ "MV Deceleration Time, ms",
      ef_modsp4 ~ "LVEF, %"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    missing = "no"
  ) %>%
  add_p(
    test = list(
      all_continuous() ~ "kruskal.test",
      all_categorical() ~ "fisher.test"
    ),
    test.args = all_tests("fisher.test") ~ list(simulate.p.value = TRUE)
  ) %>%
  bold_labels() %>%
  modify_caption(
    "**Table: Patient Characteristics by E/e'--LAVi Parameter Pattern**"
  )


# Table_ASE2025 export removed — not part of the manuscript

# ══ RCS Helpers ═══════════════════════════════════════════════════════════════
# (legacy chunk: figure3_helpers)
# ── Fig 3 setup: adjustment sets + analytic subsets ───────────────────────────
# Each index is analyzed in its own available-case sample with core adjustment.
# Additional complete-case requirements apply only to the expanded sensitivity.
rcs_adjustment_vars <- c("age", "Sex", "BMI")
rcs_adjustment_numeric <- c("age", "BMI")
figure3_adjustment_covariates <- c("bb_any", "ndhp_ccb_any", "htn_pre_test")
figure3_structural_covariates <- c(
  "lv_septal_thickness",
  "lvot_max_gradient"
)
figure3_adjustment_vars <- c(
  rcs_adjustment_vars,
  figure3_structural_covariates,
  figure3_adjustment_covariates
)
figure3_adjustment_numeric <- c(
  rcs_adjustment_numeric,
  figure3_structural_covariates,
  figure3_adjustment_covariates
)
figure3_adjustment_label <- "age, sex, and BMI"

# Missing-indicator handling for the expanded set (owner decision, card 1.5).
# The binaries become three-level factors and the resting gradient becomes a
# gated value plus a presence indicator, so the expanded sensitivity no longer
# requires covariate completeness. Septal thickness is kept as an ordinary
# covariate because it is measured in essentially everyone (1 missing of 441).
figure3_indicator_binary_covariates <- figure3_adjustment_covariates
figure3_indicator_gated_covariates <- c("lvot_max_gradient")
figure3_indicator_required_covariates <- c(
  rcs_adjustment_vars,
  "lv_septal_thickness"
)

# Each parameter is analysed in its own available-case sample: the frame holds
# every base-eligible patient (complete age/sex/BMI + >=1 CPET outcome) and
# fit_rcs_model() restricts to complete cases on the outcome, that parameter,
# and the adjustment set. Requiring one parameter therefore never drops
# patients from another parameter's model. The common E/e' + LAVi cohort is
# retained below as a consistency check. See 8_Docs/ANALYSIS_PLAN_2026-09.md.
cross_sectional_analysis_df <- baseline_df %>%
  filter(analysis_base_eligible) %>%
  mutate(
    across(
      any_of(c(
        "age",
        "BMI",
        "VO2_FRIEND2_PP",
        "VeVco2_slope",
        "HRR",
        hcm_primary_indices
      )),
      ~ as_num(.)
    ),
    Sex = factor(Sex, levels = c("Male", "Female"))
  )

write.csv(
  cross_sectional_analysis_df,
  file.path(cfg$paths$processed, "Data_CrossSectional_Analysis.csv"),
  row.names = FALSE
)

figure3_model_df <- cross_sectional_analysis_df

# ── Manuscript values: combination tests and per-parameter SD units ──────────
# The per-SD effect estimates are reported in the manuscript, so the SD each
# one is expressed in has to be quotable. It is computed on the same frame the
# peak-VO2 model uses for that parameter.
manuscript_xs_values <- bind_rows(
  tibble::tibble(
    quantity = c(
      "Kruskal-Wallis P, peak VO2 across E/e-prime-LAVi combinations",
      "Kruskal-Wallis P, VE/VCO2 slope across E/e-prime-LAVi combinations"
    ),
    value = c(
      format.pval(combo_kw_vo2$p.value, digits = 3),
      format.pval(combo_kw_vevco2$p.value, digits = 3)
    )
  ),
  purrr::map_dfr(hcm_primary_indices, function(idx) {
    frame <- cross_sectional_analysis_df %>%
      filter(complete.cases(across(all_of(c("VO2_FRIEND2_PP", idx, rcs_adjustment_vars)))))
    x <- as_num(frame[[idx]])
    tibble::tibble(
      quantity = sprintf("SD of %s, peak-VO2 model frame (n = %d)", idx, nrow(frame)),
      value = sprintf("%.2f", stats::sd(x, na.rm = TRUE))
    )
  })
)
stopifnot(nrow(manuscript_xs_values) == 5, !any(is.na(manuscript_xs_values$value)))
write.csv(
  manuscript_xs_values,
  file.path(OUT_DIR, "Table_Manuscript_CrossSectional_Values.csv"),
  row.names = FALSE
)
print(as.data.frame(manuscript_xs_values), row.names = FALSE)

figure3_indicator_build <- add_missing_indicator_covariates(
  cross_sectional_analysis_df,
  binary_vars = figure3_indicator_binary_covariates,
  gated_numeric_vars = figure3_indicator_gated_covariates
)
cross_sectional_indicator_df <- figure3_indicator_build$data
figure3_indicator_adjustment_vars <- c(
  figure3_indicator_required_covariates,
  figure3_indicator_build$covariates
)

cat(sprintf(
  "Figure 3 HCM analytic subset: %d patients\n",
  nrow(figure3_model_df)
))

# Explicit denominator ladder so readers can distinguish covariate adjustment
# from complete-case selection. Primary inference uses the first tier.
figure3_denominator_ladder <- tidyr::crossing(
  Outcome = c("Peak VO2", "VE/VCO2 slope"),
  Parameter = c("E/e'", "LAVi", "TRVmax"),
  Tier = c(
    "Primary: parameter + age/sex/BMI",
    "Sensitivity: + septal thickness/LVOT gradient",
    "Sensitivity: + structure/medications/hypertension"
  )
) %>%
  mutate(
    outcome_var = recode(Outcome,
      "Peak VO2" = "VO2_FRIEND2_PP",
      "VE/VCO2 slope" = "VeVco2_slope"
    ),
    parameter_var = recode(Parameter,
      "E/e'" = "e_e_ave",
      "LAVi" = "la_vol_index",
      "TRVmax" = "tr_max_vel"
    ),
    required_vars = pmap(
      list(outcome_var, parameter_var, Tier),
      function(outcome, parameter, tier) {
        base <- c(outcome, parameter, "age", "Sex", "BMI")
        if (grepl("septal thickness", tier)) {
          base <- c(base, figure3_structural_covariates)
        }
        if (grepl("structure/medications", tier)) {
          base <- c(
            base,
            figure3_structural_covariates,
            figure3_adjustment_covariates
          )
        }
        unique(base)
      }
    ),
    N = map_int(required_vars, function(vars) {
      sum(complete.cases(cross_sectional_analysis_df[, vars, drop = FALSE]))
    })
  ) %>%
  select(Outcome, Parameter, Tier, N)

write.csv(
  figure3_denominator_ladder,
  file.path(OUT_DIR, "Table_Model_Denominator_Ladder.csv"),
  row.names = FALSE
)

# ══ Restricted Cubic Spline Models ════════════════════════════════════════════
# (legacy chunk: figure3_models)
# ── cross-sectional RCS: index -> exercise outcome ────────────────────────────
# restricted cubic spline (4 knots; 3 exposure df, including 2 nonlinear df)
# per primary index, adjusted; LRT tests overall association and nonlinearity.
fit_figure3_outcome <- function(outcome_var, model_df = figure3_model_df,
                                run_piecewise = TRUE) {
  model_list <- list()
  lrt_list <- list()

  for (idx in hcm_primary_indices) {
    fit_obj <- fit_rcs_model(
      model_df,
      outcome_var,
      idx,
      nk = 4, # 4 knots = 2 df, standard for n this size
      adjust_vars = rcs_adjustment_vars,
      adjust_numeric = rcs_adjustment_numeric,
      run_piecewise = run_piecewise
    )
    if (is.null(fit_obj)) {
      next
    }
    model_list[[idx]] <- fit_obj
    lrt_list[[idx]] <- extract_rcs_stats(fit_obj, idx, outcome_var)
  }

  # Multiplicity adjustment is added only after both primary outcomes have
  # been fitted, so the declared family contains all six primary tests.
  lrt_table <- bind_rows(lrt_list)
  panel_stats <- tibble(index = character(), stats_label = character())

  if (length(model_list) == 0) {
    return(list(
      models = model_list,
      lrt_table = lrt_table,
      panel_stats = panel_stats
    ))
  }

  pred_curves <- bind_rows(lapply(names(model_list), function(idx) {
    m <- model_list[[idx]]
    lo <- quantile(m$df[[idx]], 0.02, na.rm = TRUE)
    hi <- quantile(m$df[[idx]], 0.98, na.rm = TRUE)
    pred_grid_rcs(
      idx,
      m$df,
      m$spline,
      m$knots,
      lo,
      hi,
      adjust_vars = m$adjust_vars
    ) %>%
      mutate(
        label = hcm_primary_labels[var],
        ase_cutoff = ase_cutoffs_hcm$ase_cutoff[match(
          var,
          ase_cutoffs_hcm$variable
        )]
      )
  }))

  point_cloud <- bind_rows(lapply(names(model_list), function(idx) {
    m <- model_list[[idx]]
    lo <- quantile(m$df[[idx]], 0.02, na.rm = TRUE)
    hi <- quantile(m$df[[idx]], 0.98, na.rm = TRUE)
    m$df %>%
      transmute(x = .data[[idx]], y = .data[[outcome_var]], var = idx) %>%
      filter(is.finite(x), is.finite(y), x >= lo, x <= hi)
  }))

  piecewise_lines <- bind_rows(lapply(names(model_list), function(idx) {
    m <- model_list[[idx]]
    if (is.null(m$piecewise)) {
      return(NULL)
    }
    lo <- quantile(m$df[[idx]], 0.02, na.rm = TRUE)
    hi <- quantile(m$df[[idx]], 0.98, na.rm = TRUE)
    pred_grid_piecewise(
      idx,
      m$df,
      m$piecewise$fit,
      m$piecewise$best_knot,
      lo,
      hi,
      adjust_vars = m$adjust_vars
    )
  }))
  if (is.null(piecewise_lines) || nrow(piecewise_lines) == 0) {
    piecewise_lines <- tibble(var = character(), x = numeric(), fit = numeric())
  }

  list(
    models = model_list,
    lrt_table = lrt_table,
    panel_stats = panel_stats,
    pred_curves = pred_curves,
    point_cloud = point_cloud,
    piecewise_lines = piecewise_lines
  )
}

# ── Consistency check: the same models in the common E/e' + LAVi cohort ───────
# Prespecified check that analysing each parameter in its own available-case
# sample does not change the conclusions drawn from the common cohort. The
# breakpoint bootstrap is not needed here, so it is skipped.
figure3_common_cohort_df <- figure3_model_df %>% filter(cross_sectional_eligible)

figure3_cohort_consistency <- bind_rows(lapply(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome_var) {
    common <- fit_figure3_outcome(
      outcome_var,
      model_df = figure3_common_cohort_df,
      run_piecewise = FALSE
    )$lrt_table %>%
      transmute(
        outcome, index, label,
        n_common_cohort = n,
        overall_p_common_cohort = overall_p,
        nonlinear_p_common_cohort = nonlinear_p
      )
    available <- fit_figure3_outcome(outcome_var, run_piecewise = FALSE)$lrt_table %>%
      transmute(
        outcome, index, label,
        n_available_case = n,
        overall_p_available_case = overall_p,
        nonlinear_p_available_case = nonlinear_p
      )
    available %>%
      left_join(common, by = c("outcome", "index", "label")) %>%
      mutate(n_added = n_available_case - n_common_cohort)
  }
))
write.csv(
  figure3_cohort_consistency,
  file.path(OUT_DIR, "Table_Consistency_CommonCohort_RCS.csv"),
  row.names = FALSE
)
print(as.data.frame(figure3_cohort_consistency), row.names = FALSE)

# ── Sensitivity: adults only ──────────────────────────────────────────────────
# The FRIEND and Wasserman reference equations were developed in adult
# populations, so percent-predicted peak VO2 is an extrapolation for the few
# patients tested below cohort$min_age_sensitivity. Every primary spline model
# is refitted with those patients excluded (owner decision, card 1.3 of
# 8_Docs/DECISION_CARDS_2026-09.md). The breakpoint bootstrap is not needed to
# compare P values, so it is skipped.
min_age_sensitivity <- cfg$cohort$min_age_sensitivity

if (!is.null(min_age_sensitivity)) {
  figure3_adult_df <- figure3_model_df %>%
    filter(!is.na(age), age >= min_age_sensitivity)

  figure3_age_sensitivity <- bind_rows(lapply(
    c("VO2_FRIEND2_PP", "VeVco2_slope"),
    function(outcome_var) {
      adults <- fit_figure3_outcome(
        outcome_var,
        model_df = figure3_adult_df,
        run_piecewise = FALSE
      )$lrt_table %>%
        transmute(
          outcome, index, label,
          n_adults_only = n,
          overall_p_adults_only = overall_p,
          nonlinear_p_adults_only = nonlinear_p
        )
      all_ages <- fit_figure3_outcome(outcome_var, run_piecewise = FALSE)$lrt_table %>%
        transmute(
          outcome, index, label,
          n_all_ages = n,
          overall_p_all_ages = overall_p,
          nonlinear_p_all_ages = nonlinear_p
        )
      all_ages %>%
        left_join(adults, by = c("outcome", "index", "label")) %>%
        mutate(
          age_floor_years = min_age_sensitivity,
          n_excluded = n_all_ages - n_adults_only,
          overall_conclusion_unchanged =
            (overall_p_all_ages < 0.05) == (overall_p_adults_only < 0.05),
          nonlinear_conclusion_unchanged =
            (nonlinear_p_all_ages < 0.05) == (nonlinear_p_adults_only < 0.05)
        )
    }
  ))

  # Benjamini-Hochberg within the same six-test family as the primary models
  # (card 1.9), so the adults-only column is judged the same way.
  figure3_age_sensitivity <- figure3_age_sensitivity %>%
    mutate(
      overall_q_all_ages = p.adjust(overall_p_all_ages, method = "BH"),
      overall_q_adults_only = p.adjust(overall_p_adults_only, method = "BH"),
      nonlinear_q_all_ages = p.adjust(nonlinear_p_all_ages, method = "BH"),
      nonlinear_q_adults_only = p.adjust(nonlinear_p_adults_only, method = "BH")
    )

  # Every primary model must appear in both columns: a dropped fit would look
  # like agreement if it silently vanished.
  stopifnot(
    nrow(figure3_age_sensitivity) ==
      2 * length(intersect(hcm_primary_indices, unique(figure3_age_sensitivity$index))),
    !any(is.na(figure3_age_sensitivity$overall_p_adults_only)),
    !any(is.na(figure3_age_sensitivity$n_adults_only)),
    all(figure3_age_sensitivity$n_adults_only <= figure3_age_sensitivity$n_all_ages)
  )

  write.csv(
    figure3_age_sensitivity,
    file.path(OUT_DIR, "Table_Sensitivity_AdultsOnly_RCS.csv"),
    row.names = FALSE
  )
  print(as.data.frame(figure3_age_sensitivity), row.names = FALSE)
}

figure3_vo2 <- fit_figure3_outcome("VO2_FRIEND2_PP")
figure3_vevco2 <- fit_figure3_outcome("VeVco2_slope")
figure3_hrr <- fit_figure3_outcome("HRR")

# ── Sensitivity: multiple imputation of the missing exposures ─────────────────
# Labelled sensitivity analysis (owner decision, card 1.2 of
# 8_Docs/DECISION_CARDS_2026-09.md). Substantive-model-compatible imputation,
# so the spline the analysis tests is the spline the imputation assumed; see
# R/imputation.R. TRVmax is never imputed as an exposure and appears
# only as a gated auxiliary, so it has no row here.
if (isTRUE(cfg$imputation$enabled)) {
  figure3_mi_primary_lrt <- bind_rows(
    figure3_vo2$lrt_table %>% mutate(mi_outcome = "VO2_FRIEND2_PP"),
    figure3_vevco2$lrt_table %>% mutate(mi_outcome = "VeVco2_slope")
  )

  figure3_mi_specs <- tidyr::crossing(
    mi_outcome = c("VO2_FRIEND2_PP", "VeVco2_slope"),
    mi_exposure = c("e_e_ave", "la_vol_index")
  )

  figure3_mi_sensitivity <- purrr::pmap_dfr(
    figure3_mi_specs,
    function(mi_outcome, mi_exposure) {
      reference <- figure3_mi_primary_lrt %>%
        filter(.data$mi_outcome == !!mi_outcome, .data$index == !!mi_exposure)
      stopifnot(nrow(reference) == 1)
      mi_rcs_exposure_sensitivity(
        cross_sectional_analysis_df,
        outcome_var = mi_outcome,
        outcome_label = ifelse(
          mi_outcome == "VO2_FRIEND2_PP",
          "Peak VO2 (% predicted)",
          "VE/VCO2 slope"
        ),
        exposure = mi_exposure,
        exposure_label = reference$label,
        auxiliary_params = setdiff(primary_parameters, mi_exposure),
        m = cfg$imputation$m,
        seed = cfg$imputation$seed,
        available_case = list(
          n = reference$n,
          overall_p = reference$overall_p,
          nonlinear_p = reference$nonlinear_p
        )
      )
    }
  )

  stopifnot(
    nrow(figure3_mi_sensitivity) == nrow(figure3_mi_specs),
    all(figure3_mi_sensitivity$`Exposure measured` +
      figure3_mi_sensitivity$`Exposure imputed` ==
      figure3_mi_sensitivity$`Patients in imputation frame`),
    all(figure3_mi_sensitivity$Imputations == cfg$imputation$m)
  )

  write.csv(
    figure3_mi_sensitivity,
    file.path(OUT_DIR, "Table_Sensitivity_MultipleImputation_RCS.csv"),
    row.names = FALSE
  )
  print(as.data.frame(figure3_mi_sensitivity), row.names = FALSE)
}

# One prespecified multiplicity family: 3 parameters x 2 primary CPET outcomes.
figure3_primary_lrt_all <- bind_rows(
  figure3_vo2$lrt_table,
  figure3_vevco2$lrt_table
) %>%
  annotate_rcs_lrt_table() %>%
  mutate(
    nonlinearity_adjustment_family =
      "Six primary parameter-outcome comparisons"
  )

figure3_vo2$lrt_table <- figure3_primary_lrt_all %>%
  filter(outcome == "VO2_FRIEND2_PP")
figure3_vevco2$lrt_table <- figure3_primary_lrt_all %>%
  filter(outcome == "VeVco2_slope")
figure3_hrr$lrt_table <- annotate_rcs_lrt_table(figure3_hrr$lrt_table) %>%
  mutate(
    nonlinearity_adjustment_family =
      "Three exploratory heart-rate-recovery comparisons"
  )

add_figure3_panel_stats <- function(result_obj) {
  result_obj$panel_stats <- if (nrow(result_obj$lrt_table) > 0) {
    result_obj$lrt_table %>%
      rowwise() %>%
      mutate(stats_label = format_rcs_stats_label(cur_data(), compact = TRUE)) %>%
      ungroup() %>%
      select(index, stats_label)
  } else {
    tibble(index = character(), stats_label = character())
  }
  result_obj
}

figure3_vo2 <- add_figure3_panel_stats(figure3_vo2)
figure3_vevco2 <- add_figure3_panel_stats(figure3_vevco2)
figure3_hrr <- add_figure3_panel_stats(figure3_hrr)

# Verify that every plotted N is exactly the complete-case denominator implied
# by its outcome, parameter, and age/sex/BMI adjustment set.
figure3_denominator_audit <- figure3_primary_lrt_all %>%
  rowwise() %>%
  mutate(
    expected_n = sum(complete.cases(
      figure3_model_df[, c(outcome, index, rcs_adjustment_vars), drop = FALSE]
    )),
    population = paste0(hcm_primary_labels[[index]], " available-case sample"),
    denominator_matches = n == expected_n
  ) %>%
  ungroup()
stopifnot(all(figure3_denominator_audit$denominator_matches))
write.csv(
  figure3_denominator_audit %>%
    select(outcome, index, population, n, expected_n, denominator_matches),
  paste0(
    file.path(OUT_DIR, ""),
    "Table_Figure3_Denominator_Audit.csv"
  ),
  row.names = FALSE
)

# Morphology sensitivity: repeat the six main index-outcome models after
# excluding apical HCM. An interaction analysis is not used because the apical
# subgroup is too small for stable nonlinear effect-modification estimates.
figure3_nonapical_df <- figure3_model_df %>%
  filter(!is.na(apical_hcm), apical_hcm == 0)

figure3_nonapical_sensitivity <- map_dfr(
  c("VO2_FRIEND2_PP", "VeVco2_slope"),
  function(outcome_var) {
    map_dfr(hcm_primary_indices, function(idx) {
      fit <- fit_rcs_model(
        figure3_nonapical_df,
        outcome_var,
        idx,
        nk = 4,
        adjust_vars = rcs_adjustment_vars,
        adjust_numeric = rcs_adjustment_numeric,
        run_piecewise = FALSE
      )
      if (is.null(fit)) {
        return(tibble())
      }
      extract_rcs_stats(fit, idx, outcome_var) %>%
        transmute(
          outcome,
          index,
          n_nonapical = n,
          overall_p_nonapical = overall_p,
          nonlinear_p_nonapical = nonlinear_p
        )
    })
  }
)

write.csv(
  figure3_nonapical_sensitivity,
  file.path(OUT_DIR, "Table_Sensitivity_NonApical_RCS.csv"),
  row.names = FALSE
)

write.csv(
  figure3_vo2$lrt_table,
  file.path(OUT_DIR, "Stats_CrossSectional_LRT_FRIEND2.csv"),
  row.names = FALSE
)
write.csv(
  figure3_vevco2$lrt_table,
  file.path(OUT_DIR, "Stats_CrossSectional_LRT_VeVco2.csv"),
  row.names = FALSE
)
write.csv(
  figure3_hrr$lrt_table,
  file.path(OUT_DIR, "Stats_CrossSectional_LRT_HRR.csv"),
  row.names = FALSE
)

build_figure3_panels <- function(
  result_obj,
  y_label,
  line_color,
  fill_color,
  stats_position = "bottom_right",
  stats_fill_alpha = 0.62
) {
  if (is.null(result_obj$pred_curves) || nrow(result_obj$pred_curves) == 0) {
    return(list())
  }
  y_limits <- range(
    c(result_obj$pred_curves$lo, result_obj$pred_curves$hi),
    na.rm = TRUE
  )
  lapply(seq_along(hcm_primary_indices), function(i) {
    idx <- hcm_primary_indices[i]
    stats_label <- result_obj$panel_stats %>%
      filter(index == idx) %>%
      pull(stats_label)
    panel_stats_position <- if (length(stats_position) >= i) {
      stats_position[[i]]
    } else {
      stats_position[[1]]
    }
    make_rcs_panel(
      pred_df = filter(result_obj$pred_curves, var == idx),
      idx = idx,
      y_label = if (i == 1) y_label else NULL,
      line_color = line_color,
      fill_color = fill_color,
      y_limits = y_limits,
      show_y = i == 1,
      stats_label = if (length(stats_label) > 0) stats_label[[1]] else NULL,
      point_df = NULL,
      stats_position = panel_stats_position,
      # Retain piecewise calculations in supplemental diagnostics, but do not
      # draw them in the primary figure because they can visually imply an
      # unsupported threshold when nonlinearity is not detected.
      piecewise_df = NULL,
      stats_fill_alpha = stats_fill_alpha
    )
  })
}

friend2_panels <- build_figure3_panels(
  figure3_vo2,
  label_peak_vo2_friend_plot,
  jacc_cols["navy"],
  jacc_cols["blue"],
  stats_position = rep("bottom_left", 3),
  stats_fill_alpha = 0.92
)
vevco2_panels <- build_figure3_panels(
  figure3_vevco2,
  label_vevco2_title_plot,
  jacc_cols["red"],
  jacc_cols["orange"],
  stats_position = rep("upper_left", 3),
  stats_fill_alpha = 0.92
)
hrr_panels <- build_figure3_panels(
  figure3_hrr,
  "Heart Rate Recovery (1 min)",
  jacc_cols["teal"],
  jacc_cols["blue"],
  stats_position = rep("bottom_left", 3),
  stats_fill_alpha = 0.62
)

if (length(friend2_panels) == 3 && length(vevco2_panels) == 3) {
  fig_rcs_combined <- wrap_plots(friend2_panels, ncol = 3) /
    wrap_plots(vevco2_panels, ncol = 3) +
    plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")") &
    theme(plot.tag = element_text(face = "bold", size = 11))

  show_and_save_jacc(
    fig_rcs_combined,
    file.path(OUT_DIR, "Figure3_Nonlinear_DiastolicIndices_RCS.pdf"),
    w = 7.0,
    h = 5.8
  )

  figure3_preprint_width <- 7.2
  figure3_preprint_height <- 4.45
  save_jacc(
    fig_rcs_combined,
    file.path(OUT_DIR, "Figure3_Nonlinear_DiastolicIndices_RCS_Preprint.pdf"),
    w = figure3_preprint_width,
    h = figure3_preprint_height
  )
  save_jacc(
    fig_rcs_combined,
    file.path(OUT_DIR, "Figure3_Nonlinear_DiastolicIndices_RCS_Preprint.pdf"),
    w = figure3_preprint_width,
    h = figure3_preprint_height
  )

  legend_fig3_pdf <- paste0(
    "**Figure 3. Left Ventricular Diastolic Parameters and Cardiopulmonary Fitness.** ",
    "Primary restricted cubic spline models adjusted for age, sex, and BMI, with each diastolic index analyzed in all patients in whom it was measured; panel annotations give the analyzed N. ",
    "Solid lines represent spline estimates with shaded 95% confidence intervals; ",
    "dashed vertical lines denote ASE abnormality thresholds; ",
    "Panel annotations report the analyzed N, overall association P value, and nonlinearity q value. ",
    "The overall P value tests whether the full spline adds information beyond age, sex, and BMI; the nonlinearity q value tests whether curvature improves fit beyond a straight-line relationship after adjustment across the six primary comparisons. ",
    "**(A)** ",
    label_peak_vo2_friend_pred_md,
    " as a function of E/e\u2019. ",
    "**(B)** ",
    label_peak_vo2_friend_pred_md,
    " as a function of LAVi. ",
    "**(C)** ",
    label_peak_vo2_friend_pred_md,
    " as a function of TRV$_max$. ",
    "**(D)** ",
    label_vevco2_md,
    " as a function of E/e\u2019. ",
    "**(E)** ",
    label_vevco2_md,
    " as a function of LAVi. ",
    "**(F)** ",
    label_vevco2_md,
    " as a function of TRV$_max$."
  )

  save_jacc_with_embedded_legend(
    fig_rcs_combined,
    file.path(OUT_DIR, "Figure3_Nonlinear_DiastolicIndices_RCS.pdf"),
    legend_text = legend_fig3_pdf,
    w = 7.0,
    h = 5.8
  )
  save_jacc_with_embedded_legend(
    fig_rcs_combined,
    file.path(OUT_DIR, "Figure3_Nonlinear_DiastolicIndices_RCS.pdf"),
    legend_text = legend_fig3_pdf,
    w = 7.0,
    h = 5.8
  )
}

# ══ LRTs and Curvature Statistics ═════════════════════════════════════════════
# (legacy chunk: figure3_summary_table)
# ── table: Fig 3 RCS LRT stats (overall + nonlinearity P per index/outcome) ───
figure3_summary_table_raw <- bind_rows(
  figure3_vo2$lrt_table %>% mutate(Outcome = label_peak_vo2_friend_pred),
  figure3_vevco2$lrt_table %>% mutate(Outcome = label_vevco2)
) %>%
  mutate(
    outcome_order = match(Outcome, c(label_peak_vo2_friend_pred, label_vevco2)),
    parameter_order = match(index, hcm_primary_indices)
  ) %>%
  arrange(outcome_order, parameter_order)

figure3_summary_table_display <- figure3_summary_table_raw %>%
  transmute(
    Outcome = case_when(
      Outcome == label_peak_vo2_friend_pred ~ label_peak_vo2,
      TRUE ~ Outcome
    ),
    Parameter = label,
    Population = ifelse(
      index == "tr_max_vel",
      "TRVmax available-case subset",
      "E/e' + LAVi cohort"
    ),
    N = n,
    `Overall p` = ifelse(
      overall_p < 0.001,
      "<0.001",
      sprintf("%.3f", overall_p)
    ),
    `Overall q` = ifelse(
      overall_q < 0.001,
      "<0.001",
      sprintf("%.3f", overall_q)
    ),
    `Nonlinearity q` = ifelse(
      nonlinear_q < 0.001,
      "<0.001",
      sprintf("%.3f", nonlinear_q)
    ),
    # Inference is on the adjusted value now that the family is adjusted
    # throughout (card 1.9).
    Interpretation = case_when(
      overall_q < 0.05 & nonlinear_q < 0.05 ~
        "Association with evidence of nonlinearity",
      overall_q < 0.05 ~
        "Association; no evidence of nonlinearity",
      TRUE ~ "No clear association"
    ),
    overall_p_num = overall_p,
    overall_q_num = overall_q,
    nonlinear_q_num = nonlinear_q
  )

# Model-selection and breakpoint diagnostics are retained for supplemental
# review without burdening the primary table or figure.
figure3_diagnostics_table <- figure3_summary_table_raw %>%
  transmute(
    Outcome = case_when(
      Outcome == label_peak_vo2_friend_pred ~ label_peak_vo2,
      TRUE ~ Outcome
    ),
    Parameter = label,
    Population = ifelse(
      index == "tr_max_vel",
      "TRVmax available-case subset",
      "E/e' + LAVi cohort"
    ),
    N = n,
    `Linear-component p` = linear_p,
    `Overall spline p` = overall_p,
    `Overall spline q (six-test family)` = overall_q,
    `Nonlinearity p` = nonlinear_p,
    `Nonlinearity q (six-test family)` = nonlinear_q,
    `Spline vs linear Delta AIC` = delta_AIC_nonlinear,
    `Piecewise knot` = piecewise_knot,
    `Piecewise bootstrap p` = piecewise_p_used,
    `Piecewise q` = piecewise_q,
    `Piecewise knot support range` = knot_support_range,
    `Curvature zone` = curvature_zone
  )

write.csv(
  figure3_diagnostics_table,
  paste0(
    file.path(OUT_DIR, ""),
    "TableS_Figure3_RCS_Diagnostics.csv"
  ),
  row.names = FALSE
)
write.csv(
  figure3_diagnostics_table,
  paste0(
    file.path(OUT_DIR, ""),
    "TableS_Figure3_RCS_Diagnostics.csv"
  ),
  row.names = FALSE
)

figure3_outcome_break_rows <- figure3_summary_table_display %>%
  count(Outcome, name = "n_rows") %>%
  pull(n_rows) %>%
  cumsum()

figure3_summary_table <- figure3_summary_table_display %>%
  select(-overall_p_num, -overall_q_num, -nonlinear_q_num)

figure3_ft <- figure3_summary_table_display %>%
  flextable(
    col_keys = c(
      "Outcome",
      "Parameter",
      "Population",
      "N",
      "Overall p",
      "Overall q",
      "Nonlinearity q",
      "Interpretation"
    )
  ) %>%
  merge_v(j = "Outcome") %>%
  valign(j = "Outcome", valign = "top", part = "body") %>%
  align(
    j = c("Outcome", "Parameter", "Population", "Interpretation"),
    align = "left",
    part = "all"
  ) %>%
  align(
    j = c(
      "N",
      "Overall p",
      "Overall q",
      "Nonlinearity q"
    ),
    align = "center",
    part = "all"
  ) %>%
  bold(i = ~ overall_q_num < 0.05, j = "Overall q", part = "body") %>%
  bold(i = ~ nonlinear_q_num < 0.05, j = "Nonlinearity q", part = "body") %>%
  format_pub_table(
    caption = paste0(
      "Table 2. Primary age-, sex-, and BMI-adjusted restricted cubic spline ",
      "results. Q values are Benjamini-Hochberg adjusted within the ",
      "prespecified family of six parameter-outcome comparisons, applied to ",
      "the overall association and the nonlinearity test alike."
    )
  ) %>%
  hline(
    i = figure3_outcome_break_rows,
    border = officer::fp_border(color = "black", width = 0.8),
    part = "body"
  )

figure3_ft <- apply_subscript_format(
  figure3_ft,
  figure3_summary_table_display,
  "Outcome"
)


write.csv(
  figure3_summary_table,
  file.path(OUT_DIR, "Table2_Figure3_RCS_Summary.csv"),
  row.names = FALSE
)
figure3_ft %>%
  save_as_docx(path = file.path(OUT_DIR, "Table2_Figure3_RCS_Summary.docx"))

# (legacy chunk: figure3_validation)
# ── Fig 3 sanity check: print analytic N + key effect sizes ───────────────────
cat("**Figure 3 Validation**\n\n")
cat(sprintf(
  "- Exact-date medication-adjusted analytic subset: **%d** patients\n",
  nrow(figure3_model_df)
))
cat(sprintf("- Adjustment set: **%s**\n", figure3_adjustment_label))

# ── Cross-correlation of diastolic indices ─────────────────────────────────
cat("\n**Diastolic index pairwise correlations (Spearman):**\n\n")
corr_df <- figure3_model_df %>%
  select(any_of(hcm_primary_indices)) %>%
  mutate(across(everything(), ~ as_num(.)))
pairs <- combn(hcm_primary_indices, 2, simplify = FALSE)
for (pr in pairs) {
  complete <- corr_df %>% filter(!is.na(.data[[pr[1]]]), !is.na(.data[[pr[2]]]))
  if (nrow(complete) >= 10) {
    rho <- cor(complete[[pr[1]]], complete[[pr[2]]], method = "spearman")
    cat(sprintf(
      "- %s vs %s: r = %.3f (n = %d)\n",
      hcm_primary_labels[pr[1]],
      hcm_primary_labels[pr[2]],
      rho,
      nrow(complete)
    ))
  }
}

# ── Per-model sample sizes ─────────────────────────────────────────────────
cat("\n**Per-model sample sizes:**\n\n")
for (outcome_name in c("VO2_FRIEND2_PP", "VeVco2_slope")) {
  result <- switch(
    outcome_name,
    VO2_FRIEND2_PP = figure3_vo2,
    VeVco2_slope = figure3_vevco2
  )
  for (idx in hcm_primary_indices) {
    m <- result$models[[idx]]
    if (!is.null(m)) {
      cat(sprintf(
        "- %s ~ %s: n = %d\n",
        outcome_name,
        hcm_primary_labels[idx],
        nrow(m$df)
      ))
    }
  }
}

# ── Verify hypertension covariate is present ───────────────────────────────
cat(sprintf(
  "\n- Hypertension (htn_pre_test) available: **%s** (non-missing: %d / %d)\n",
  ifelse("htn_pre_test" %in% names(figure3_model_df), "yes", "NO"),
  sum(!is.na(figure3_model_df$htn_pre_test)),
  nrow(figure3_model_df)
))

# ══ Same-Sample Adjustment and Morphology Sensitivities ═══════════════════════
# (legacy chunk: figure3_same_sample_sensitivity)
# These linear sensitivity models deliberately hold the analyzed rows fixed
# while covariates are added. This separates adjustment effects from losses due
# to complete-case selection; the nonlinear RCS models above remain primary.
figure3_outcome_specs <- c(
  VO2_FRIEND2_PP = "Peak VO2 (% predicted)",
  VeVco2_slope = "VE/VCO2 slope"
)
figure3_parameter_specs <- c(
  e_e_ave = "E/e'",
  la_vol_index = "LAVi",
  tr_max_vel = "TRVmax"
)
figure3_primary_sample_label <- if (
  analytic_cohort_mode == "complete_three"
) {
  "Primary complete-three cohort"
} else {
  "Parameter available-case sample"
}

fit_same_sample_lm <- function(data, outcome, parameter, covariates) {
  required <- unique(c(outcome, parameter, covariates))
  model_df <- data %>%
    select(all_of(required)) %>%
    drop_na() %>%
    mutate(parameter_z = as.numeric(scale(.data[[parameter]])))

  if (nrow(model_df) < 20 || sd(model_df$parameter_z) == 0) {
    return(NULL)
  }

  fit <- lm(
    reformulate(c("parameter_z", covariates), response = outcome),
    data = model_df
  )
  result <- broom::tidy(fit, conf.int = TRUE) %>%
    filter(term == "parameter_z")
  if (nrow(result) != 1) return(NULL)

  list(fit = fit, data = model_df, result = result)
}

figure3_same_sample_rows <- map_dfr(
  names(figure3_outcome_specs),
  function(outcome) {
    map_dfr(names(figure3_parameter_specs), function(parameter) {
      # Full primary cohort: same rows for crude and minimally adjusted fits.
      primary_vars <- c(outcome, parameter, rcs_adjustment_vars)
      primary_df <- cross_sectional_analysis_df %>%
        filter(complete.cases(across(all_of(primary_vars))))

      # Expanded tier: the same smaller set is used for all three fits.
      expanded_vars <- c(outcome, parameter, figure3_adjustment_vars)
      expanded_df <- cross_sectional_analysis_df %>%
        filter(complete.cases(across(all_of(expanded_vars))))

      # Indicator-handled expanded tier: the same expanded covariates, but
      # missingness is modelled instead of excluded, so the sample is the
      # primary one (card 1.5).
      indicator_vars <- c(
        outcome,
        parameter,
        figure3_indicator_required_covariates
      )
      indicator_df <- cross_sectional_indicator_df %>%
        filter(complete.cases(across(all_of(indicator_vars))))
      indicator_covariates <- drop_degenerate_covariates(
        indicator_df,
        figure3_indicator_adjustment_vars
      )

      model_specs <- list(
        list(
          sample = figure3_primary_sample_label,
          model = "Unadjusted",
          data = primary_df,
          covariates = character(0)
        ),
        list(
          sample = figure3_primary_sample_label,
          model = "Age/sex/BMI adjusted",
          data = primary_df,
          covariates = rcs_adjustment_vars
        ),
        list(
          sample = "Expanded-covariate complete subset",
          model = "Unadjusted",
          data = expanded_df,
          covariates = character(0)
        ),
        list(
          sample = "Expanded-covariate complete subset",
          model = "Age/sex/BMI adjusted",
          data = expanded_df,
          covariates = rcs_adjustment_vars
        ),
        list(
          sample = "Expanded-covariate complete subset",
          model = paste0(
            "Age/sex/BMI + septal thickness/LVOT gradient/",
            "medications/hypertension"
          ),
          data = expanded_df,
          covariates = figure3_adjustment_vars
        ),
        list(
          sample = "Expanded-covariate indicator-handled sample",
          model = "Unadjusted",
          data = indicator_df,
          covariates = character(0)
        ),
        list(
          sample = "Expanded-covariate indicator-handled sample",
          model = "Age/sex/BMI adjusted",
          data = indicator_df,
          covariates = rcs_adjustment_vars
        ),
        list(
          sample = "Expanded-covariate indicator-handled sample",
          model = paste0(
            "Age/sex/BMI + septal thickness/LVOT gradient/",
            "medications/hypertension (missing-indicator handled)"
          ),
          data = indicator_df,
          covariates = indicator_covariates
        )
      )

      map_dfr(model_specs, function(spec) {
        fit_obj <- fit_same_sample_lm(
          spec$data,
          outcome,
          parameter,
          spec$covariates
        )
        if (is.null(fit_obj)) return(tibble())
        fit_obj$result %>%
          transmute(
            Outcome = unname(figure3_outcome_specs[[outcome]]),
            Parameter = unname(figure3_parameter_specs[[parameter]]),
            `Analysis sample` = spec$sample,
            Model = spec$model,
            N = nrow(fit_obj$data),
            `Beta per 1 SD higher` = estimate,
            `CI low` = conf.low,
            `CI high` = conf.high,
            `P value` = p.value,
            `Adjusted R2` = summary(fit_obj$fit)$adj.r.squared
          )
      })
    })
  }
)

figure3_same_sample_rows <- figure3_same_sample_rows %>%
  mutate(
    `Inference tier` = case_when(
      `Analysis sample` == "Expanded-covariate complete subset" ~
        "Supplemental expanded-covariate sensitivity",
      `Analysis sample` == "Expanded-covariate indicator-handled sample" ~
        "Supplemental expanded-covariate sensitivity (missing-indicator handled)",
      TRUE ~ "Primary same-sample robustness"
    ),
    .after = `Analysis sample`
  )

# Every crude/adjusted comparison must retain an identical denominator within
# its outcome, parameter, and analysis-sample stratum.
figure3_same_sample_n_audit <- figure3_same_sample_rows %>%
  group_by(Outcome, Parameter, `Analysis sample`) %>%
  summarise(n_distinct_n = n_distinct(N), .groups = "drop")
stopifnot(all(figure3_same_sample_n_audit$n_distinct_n == 1))

if (analytic_cohort_mode == "complete_three") {
  figure3_primary_n_audit <- figure3_same_sample_rows %>%
    filter(`Analysis sample` == figure3_primary_sample_label) %>%
    summarise(n_distinct_n = n_distinct(N)) %>%
    pull(n_distinct_n)
  figure3_expanded_n_audit <- figure3_same_sample_rows %>%
    filter(`Analysis sample` == "Expanded-covariate complete subset") %>%
    summarise(n_distinct_n = n_distinct(N)) %>%
    pull(n_distinct_n)
  stopifnot(
    figure3_primary_n_audit == 1,
    figure3_expanded_n_audit == 1
  )
}

write.csv(
  figure3_same_sample_rows,
  paste0(
    file.path(OUT_DIR, ""),
    "Table_Sensitivity_SameSample_Adjustment.csv"
  ),
  row.names = FALSE
)

# Direct selection sensitivity: refit the E/e' and LAVi linear associations
# after additionally requiring a measured TRVmax. This quantifies what changes
# because of that requirement, without changing the exposure, outcome, or
# adjustment set.
figure3_cohort_selection_rows <- map_dfr(
  names(figure3_outcome_specs),
  function(outcome) {
    map_dfr(c("e_e_ave", "la_vol_index"), function(parameter) {
      sample_specs <- list(
        list(
          label = "Parameter available-case sample",
          data = cross_sectional_analysis_df
        ),
        list(
          label = "TRVmax-measured subset",
          data = cross_sectional_analysis_df %>% filter(!is.na(tr_max_vel))
        )
      )
      map_dfr(sample_specs, function(sample_spec) {
        fit_obj <- fit_same_sample_lm(
          sample_spec$data,
          outcome,
          parameter,
          rcs_adjustment_vars
        )
        if (is.null(fit_obj)) return(tibble())
        fit_obj$result %>%
          transmute(
            Outcome = unname(figure3_outcome_specs[[outcome]]),
            Parameter = unname(figure3_parameter_specs[[parameter]]),
            `Analysis sample` = sample_spec$label,
            N = nrow(fit_obj$data),
            `Beta per 1 SD higher` = estimate,
            `CI low` = conf.low,
            `CI high` = conf.high,
            `P value` = p.value
          )
      })
    })
  }
)

write.csv(
  figure3_cohort_selection_rows,
  paste0(
    file.path(AUDIT_DIR, ""),
    "Table_CohortSelection_CrossSectional_EeLAVi.csv"
  ),
  row.names = FALSE
)

# Morphology sensitivity on one fixed known-morphology subset. The exposure
# coefficient before and after adding apical phenotype can therefore be compared
# without a denominator change. Unknown morphology is never treated as nonapical.
figure3_morphology_rows <- map_dfr(
  names(figure3_outcome_specs),
  function(outcome) {
    map_dfr(names(figure3_parameter_specs), function(parameter) {
      required <- c(
        outcome,
        parameter,
        rcs_adjustment_vars,
        "apical_hcm"
      )
      morphology_df <- cross_sectional_analysis_df %>%
        filter(complete.cases(across(all_of(required))))

      model_specs <- list(
        list(
          model = "Age/sex/BMI adjusted",
          covariates = rcs_adjustment_vars
        ),
        list(
          model = "Age/sex/BMI + apical morphology",
          covariates = c(rcs_adjustment_vars, "apical_hcm")
        )
      )

      map_dfr(model_specs, function(spec) {
        fit_obj <- fit_same_sample_lm(
          morphology_df,
          outcome,
          parameter,
          spec$covariates
        )
        if (is.null(fit_obj)) return(tibble())
        morphology_term <- broom::tidy(fit_obj$fit, conf.int = TRUE) %>%
          filter(term == "apical_hcm")
        fit_obj$result %>%
          transmute(
            Outcome = unname(figure3_outcome_specs[[outcome]]),
            Parameter = unname(figure3_parameter_specs[[parameter]]),
            Model = spec$model,
            N = nrow(fit_obj$data),
            `Apical HCM, n` = sum(morphology_df$apical_hcm == 1),
            `Nonapical HCM, n` = sum(morphology_df$apical_hcm == 0),
            `Parameter beta per 1 SD higher` = estimate,
            `Parameter CI low` = conf.low,
            `Parameter CI high` = conf.high,
            `Parameter P value` = p.value,
            `Apical morphology beta` = ifelse(
              nrow(morphology_term) == 1,
              morphology_term$estimate,
              NA_real_
            ),
            `Apical morphology P value` = ifelse(
              nrow(morphology_term) == 1,
              morphology_term$p.value,
              NA_real_
            )
          )
      })
    })
  }
)

write.csv(
  figure3_morphology_rows,
  paste0(
    file.path(OUT_DIR, ""),
    "Table_Sensitivity_Morphology_SameSample.csv"
  ),
  row.names = FALSE
)


# ══ Figure S1: All-Parameter RCS ══════════════════════════════════════════════
# (legacy chunk: supplemental_s1)
# ── Fig S1: RCS curves for all diastolic params x outcomes (full grid) ─────────
supp_outcome_specs <- tribble(
  ~outcome         , ~outcome_label                       , ~curve_color        , ~fill_color         ,
  "VO2_FRIEND2_PP" , "Peak V̇O<sub>2</sub> (% Predicted)" , jacc_cols["navy"]   , jacc_cols["blue"]   ,
  "VeVco2_slope"   , "V̇E/V̇CO<sub>2</sub> Slope"         , jacc_cols["orange"] , jacc_cols["orange"]
)

supp_rcs_all_curves <- list()
supp_rcs_all_stats <- list()

for (i in seq_len(nrow(supp_outcome_specs))) {
  spec <- supp_outcome_specs[i, ]
  outcome_curves <- list()
  outcome_stats <- list()

  for (idx in diastolic_indices) {
    fit_obj <- fit_rcs_model(
      figure3_model_df,
      spec$outcome,
      idx,
      nk = 4,
      adjust_vars = rcs_adjustment_vars,
      adjust_numeric = rcs_adjustment_numeric,
      run_piecewise = FALSE
    )
    if (is.null(fit_obj)) {
      next
    }

    lo <- quantile(fit_obj$df[[idx]], 0.02, na.rm = TRUE)
    hi <- quantile(fit_obj$df[[idx]], 0.98, na.rm = TRUE)

    outcome_curves[[idx]] <- pred_grid_rcs(
      idx,
      fit_obj$df,
      fit_obj$spline,
      fit_obj$knots,
      lo,
      hi,
      adjust_vars = fit_obj$adjust_vars
    ) %>%
      mutate(
        outcome = spec$outcome,
        outcome_label = spec$outcome_label,
        label = label_map[var],
        panel_label = paste0(spec$outcome_label, " | ", label_map[var]),
        curve_color = spec$curve_color,
        fill_color = spec$fill_color
      ) %>%
      left_join(ase_cutoffs_all, by = c("var" = "variable"))

    stats_row <- tryCatch(
      extract_rcs_stats(fit_obj, idx, spec$outcome) %>%
        mutate(label = label_map[idx]),
      error = function(e) NULL
    )
    if (!is.null(stats_row)) outcome_stats[[idx]] <- stats_row
  }

  supp_rcs_all_curves[[i]] <- bind_rows(outcome_curves)

  if (length(outcome_stats) > 0) {
    supp_rcs_all_stats[[i]] <- annotate_rcs_lrt_table(bind_rows(
      outcome_stats
    )) %>%
      mutate(
        outcome_label = spec$outcome_label,
        panel_label = paste0(spec$outcome_label, " | ", label_map[index])
      )
  }
}

supp_rcs_all_curves <- bind_rows(supp_rcs_all_curves)

supp_rcs_stats_df <- if (length(supp_rcs_all_stats) > 0) {
  bind_rows(supp_rcs_all_stats) %>%
    rowwise() %>%
    mutate(stats_label = format_rcs_stats_label(cur_data())) %>%
    ungroup() %>%
    select(panel_label, stats_label) %>%
    mutate(x = Inf, y = -Inf)
} else {
  tibble(
    panel_label = character(),
    stats_label = character(),
    x = numeric(),
    y = numeric()
  )
}

if (nrow(supp_rcs_all_curves) > 0) {
  fig_s1_rcs_all <- ggplot(supp_rcs_all_curves, aes(x = x, y = fit)) +
    geom_ribbon(
      aes(ymin = lo, ymax = hi, fill = outcome_label),
      alpha = 0.14,
      color = NA
    ) +
    geom_line(aes(color = outcome_label), linewidth = 0.7) +
    geom_vline(
      data = supp_rcs_all_curves %>%
        distinct(panel_label, ase_cutoff) %>%
        filter(!is.na(ase_cutoff)),
      aes(xintercept = ase_cutoff),
      inherit.aes = FALSE,
      color = jacc_cols["blue"],
      linetype = "dashed",
      linewidth = 0.4
    ) +
    {
      if (nrow(supp_rcs_stats_df) > 0) {
        ggtext::geom_richtext(
          data = supp_rcs_stats_df,
          aes(x = x, y = y, label = stats_label),
          inherit.aes = FALSE,
          hjust = 1.02,
          vjust = -0.15,
          size = 1.9,
          lineheight = 1.1,
          label.size = 0.2,
          fill = alpha("white", 0.88),
          color = "black",
          label.colour = "black"
        )
      }
    } +
    facet_wrap(~panel_label, scales = "free", ncol = 4) +
    scale_color_manual(
      values = setNames(
        supp_outcome_specs$curve_color,
        supp_outcome_specs$outcome_label
      )
    ) +
    scale_fill_manual(
      values = setNames(
        supp_outcome_specs$fill_color,
        supp_outcome_specs$outcome_label
      )
    ) +
    labs(
      title = "Supplemental Restricted Cubic Splines Across All Available Diastolic Parameters",
      subtitle = paste0(
        "Age-, sex-, and BMI-adjusted spline fits for ",
        label_peak_vo2_plot,
        " and ",
        label_vevco2_plot
      ),
      x = NULL,
      y = "Predicted outcome"
    ) +
    theme_jacc() +
    theme(
      strip.text = ggtext::element_markdown(size = 6),
      axis.text = element_text(size = 6),
      axis.title = ggtext::element_markdown(size = 8),
      # subtitle + colour-legend labels carry <sub> markup; render as markdown
      # so V̇O2/V̇E·V̇CO2 subscripts typeset instead of printing literal tags
      plot.subtitle = ggtext::element_markdown(size = 9, hjust = 0.5, color = "grey40"),
      legend.text = ggtext::element_markdown(size = 8),
      legend.position = "bottom"
    )

  show_and_save_jacc(
    fig_s1_rcs_all,
    file.path(OUT_DIR, "Figure_S1_AllDiastolic_RCS.pdf"),
    w = 11,
    h = 10.5
  )
  save_jacc_with_embedded_legend(
    fig_s1_rcs_all,
    file.path(OUT_DIR, "Figure_S1_AllDiastolic_RCS.pdf"),
    legend_text = paste0(
      "**Figure S1.** Supplemental restricted cubic spline models relating all available diastolic parameters to ",
      label_peak_vo2_md,
      " and ",
      label_vevco2_md,
      ". ",
      "Each panel shows age-, sex-, and BMI-adjusted spline-estimated associations with shaded 95% confidence intervals. ",
      "Dashed blue vertical lines indicate ASE reference cutoffs where available. ",
      "Annotation boxes report patient count, overall association P value, nonlinearity FDR q-value ",
      "(BH-adjusted within each outcome across all 12 diastolic indices), and ΔAIC (spline vs. linear model)."
    ),
    w = 11,
    h = 10.5
  )
}


message("Stage 05 complete.")

message("Stage 05 complete.")
