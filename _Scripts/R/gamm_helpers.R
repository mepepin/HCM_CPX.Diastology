# ── Figure 4 rendering constants ─────────────────────────────────────────────
figure4_lvdd_palette <- c(
  "1 abnormal parameter"  = "#C9A2A0",
  "2 abnormal parameters" = "#B97E79",
  "All 3 abnormal"        = "#8E5C59"
)
figure4_lvdd_labels <- c(
  "1" = "1 abnormal parameter",
  "2" = "2 abnormal parameters",
  "3" = "All 3 abnormal"
)
figure4_status_palette  <- c("Normal"  = "#B8C1C9", "Abnormal" = "#B5423E")
figure4_quartile_palette <- c(
  # Preserve the original sequential-blue visual language. The requested
  # analysis displays two halves rather than four quartiles, but that analytical
  # change should not introduce a new gray/rose design.
  "Lower half" = "#B7CAD7",
  "Upper half" = "#2F5873"
)

# ── Shared prediction-grid covariate filler ───────────────────────────────────
fill_pred_covariates <- function(pred_df, model_df, covariates) {
  if ("age" %in% covariates)
    pred_df$age <- median(model_df$age, na.rm = TRUE)
  if ("Sex" %in% covariates)
    pred_df$Sex <- factor(names(which.max(table(model_df$Sex))),
                          levels = levels(model_df$Sex))
  if ("BMI" %in% covariates)
    pred_df$BMI <- median(model_df$BMI, na.rm = TRUE)

  for (covar in setdiff(covariates, c("age", "Sex", "BMI"))) {
    vals <- model_df[[covar]]
    pred_df[[covar]] <- if (is.factor(vals)) {
      factor(mode_nonmissing(vals), levels = levels(vals))
    } else if (is.numeric(vals) || is.integer(vals)) {
      uv <- unique(na.omit(vals))
      if (length(uv) && all(uv %in% c(0, 1))) as.numeric(mode_nonmissing(vals))
      else median(vals, na.rm = TRUE)
    } else {
      mode_nonmissing(vals)
    }
  }
  pred_df
}

# ── Shared model-fitting utilities ────────────────────────────────────────────
.out_stub <- function(outcome_var) {
  case_when(outcome_var == "VO2_FRIEND2_PP" ~ "VO2",
            outcome_var == "VeVco2_slope"   ~ "VEVCO2",
            TRUE                            ~ outcome_var)
}

.write_gamm_summary <- function(fit_summary, output_prefix) {
  write.csv(as.data.frame(fit_summary$s.table), paste0(output_prefix, "_SmoothTerms.csv"))
  write.csv(as.data.frame(fit_summary$p.table), paste0(output_prefix, "_ParametricTerms.csv"))
}

.gamm_prep_df <- function(data, outcome_var, extra_vars, model_covariates) {
  numeric_covariates <- setdiff(model_covariates, "Sex")
  data %>%
    mutate(across(any_of(c(outcome_var, "time_yrs", extra_vars, numeric_covariates)), as_num),
           ID_fac = factor(ID),
           Sex    = factor(Sex, levels = c("Male", "Female")))
}

.gamm_fit_with_fallback <- function(smooth_formula, fallback_formula, model_df) {
  # Primary fit: mgcv::gamm() with random intercept + random slope on time_yrs
  # per patient. Methods text describes this structure; previous bam()
  # implementation had random intercept only.
  fit_mode <- "smooth"
  fit <- tryCatch(
    mgcv::gamm(smooth_formula, data = model_df,
               random = list(ID_fac = nlme::pdDiag(~ 1 + time_yrs)),
               method = "REML",
               control = nlme::lmeControl(opt = "optim", maxIter = 200,
                                          msMaxIter = 200,
                                          niterEM = 100, tolerance = 1e-6,
                                          msTol = 1e-6)),
    error = function(e) {
      fit_mode <<- "fallback"
      tryCatch(
        mgcv::gamm(fallback_formula, data = model_df,
                   random = list(ID_fac = nlme::pdDiag(~ 1 + time_yrs)),
                   method = "REML",
                   control = nlme::lmeControl(opt = "optim", maxIter = 200,
                                              msMaxIter = 200,
                                              niterEM = 100, tolerance = 1e-6,
                                              msTol = 1e-6)),
        error = function(e2) {
          # Final fallback: drop random slope when lme refuses to converge.
          fit_mode <<- "fallback_intercept_only"
          mgcv::gamm(fallback_formula, data = model_df,
                     random = list(ID_fac = ~ 1),
                     method = "REML",
                     control = nlme::lmeControl(opt = "optim", maxIter = 200,
                                                msMaxIter = 200,
                                                niterEM = 100, tolerance = 1e-6,
                                                msTol = 1e-6))
        }
      )
    }
  )
  list(fit = fit, fit_mode = fit_mode)
}

.gamm_predict_ci <- function(fit, pred_df) {
  # For gamm objects, prediction is via the $gam slot at population-level
  # (random effects are integrated out automatically — no exclude= needed).
  gam_obj <- if (inherits(fit, "list") && !is.null(fit$gam)) fit$gam else fit
  pt <- predict(gam_obj, newdata = pred_df, type = "response", se.fit = TRUE)
  pred_df %>% mutate(fit = pt$fit, se = pt$se.fit,
                     lo  = fit - 1.96 * se, hi = fit + 1.96 * se)
}

# Helper: extract the gam component from gamm() output for summary/coef calls.
.gamm_gam <- function(fit_obj) {
  if (inherits(fit_obj, "list") && !is.null(fit_obj$gam)) fit_obj$gam
  else fit_obj
}

# Residual df: for gamm() output the df.residual lives on $gam; for bam/gam
# objects it lives at the top level. Used by downstream p-value extractors.
.gamm_dfresid <- function(fit_obj) {
  gam_obj <- .gamm_gam(fit_obj)
  tryCatch(gam_obj$df.residual, error = \(e) NA_real_)
}

# Relabel mgcv's by-ordered-factor smooth term names to the legacy convention
# used by downstream extract_*_model_rows code. mgcv outputs terms like
# "s(time_yrs):status_binary_ordAbnormal"; we map them back to
# "ti(time_yrs,status_abnormal)" so all extraction logic continues to work.
.relabel_by_smooth <- function(smooth_tbl, ordered_var, legacy_name) {
  if (!nrow(smooth_tbl)) return(smooth_tbl)
  pattern <- paste0("^s\\(time_yrs\\):", ordered_var, ".+$")
  smooth_tbl %>%
    mutate(term = ifelse(grepl(pattern, term, perl = TRUE), legacy_name, term))
}

# ── GAMM fitting functions ────────────────────────────────────────────────────
fit_gamm_binary_panel <- function(data, outcome_var, spec,
                                  covariates = figure4_core_covariates) {
  status_col         <- spec$status_var
  model_covariates   <- unique(covariates)
  numeric_covariates <- setdiff(model_covariates, "Sex")

  model_df <- .gamm_prep_df(data, outcome_var, character(0), model_covariates) %>%
    mutate(
      status_abnormal = case_when(
        is.na(.data[[status_col]])  ~ NA_real_,
        .data[[status_col]]         ~ 1,
        TRUE                        ~ 0
      ),
      status_binary = factor(
        case_when(is.na(status_abnormal) ~ NA_character_,
                  status_abnormal == 1  ~ "Abnormal",
                  TRUE                  ~ "Normal"),
        levels = c("Normal", "Abnormal")
      ),
      # Ordered factor enables the difference-smooth varying-coefficient
      # parameterization in mgcv: s(time_yrs, by=status_binary_ord) returns
      # a single smooth representing the trajectory deviation of Abnormal vs
      # Normal (the reference level).
      status_binary_ord = ordered(status_binary, levels = c("Normal", "Abnormal"))
    ) %>%
    filter(!is.na(.data[[outcome_var]]), !is.na(status_binary),
           !is.na(time_yrs), is.finite(time_yrs))

  if (length(model_covariates))
    model_df <- model_df %>% filter(if_all(all_of(model_covariates), ~ !is.na(.x)))

  if (nrow(model_df) < 50 || n_distinct(model_df$ID) < 20 ||
      n_distinct(model_df$status_binary) < 2) return(NULL)

  # Varying-coefficient form: reference smooth + difference smooth.
  # Random intercept + slope per patient supplied via gamm(random=...).
  rhs_smooth   <- paste(c("s(time_yrs, bs='cr', k=5)",
                           "s(time_yrs, by=status_binary_ord, bs='cr', k=5)",
                           "status_binary",
                           model_covariates), collapse = " + ")
  rhs_fallback <- paste(c("status_binary * time_yrs",
                           model_covariates), collapse = " + ")
  smooth_formula   <- as.formula(paste(outcome_var, "~", rhs_smooth))
  fallback_formula <- as.formula(paste(outcome_var, "~", rhs_fallback))

  res         <- .gamm_fit_with_fallback(smooth_formula, fallback_formula, model_df)
  fit         <- res$fit
  fit_gam     <- .gamm_gam(fit)
  fit_summary <- summary(fit_gam)
  output_prefix <- paste0("../2_Output/Stats_GAMM_", .out_stub(outcome_var),
                          "_", spec$suffix, "Binary")
  .write_gamm_summary(fit_summary, output_prefix)

  pred_df <- expand_grid(
    time_yrs      = seq(0, min(10, max(model_df$time_yrs, na.rm = TRUE)), by = 0.25),
    status_binary = factor(c("Normal", "Abnormal"), levels = c("Normal", "Abnormal"))
  ) %>%
    mutate(status_abnormal = ifelse(status_binary == "Abnormal", 1, 0),
           status_binary_ord = ordered(status_binary,
                                       levels = c("Normal", "Abnormal")),
           ID_fac = model_df$ID_fac[1]) %>%
    fill_pred_covariates(model_df, model_covariates) %>%
    .gamm_predict_ci(fit, .)

  smooth_tbl <- as_tibble(fit_summary$s.table, rownames = "term") %>%
    .relabel_by_smooth(ordered_var = "status_binary_ord",
                       legacy_name = "ti(time_yrs,status_abnormal)")

  list(
    fit = fit, fit_mode = res$fit_mode,
    smooth_terms    = smooth_tbl,
    parametric_terms = as_tibble(fit_summary$p.table, rownames = "term"),
    pred_df          = pred_df,
    n_observations   = nrow(model_df),
    n_patients       = n_distinct(model_df$ID),
    n_normal_patients   = model_df %>% distinct(ID, status_binary) %>% filter(status_binary == "Normal") %>% nrow(),
    n_abnormal_patients = model_df %>% distinct(ID, status_binary) %>% filter(status_binary == "Abnormal") %>% nrow(),
    outcome          = outcome_var,
    spec             = spec,
    status_caption   = spec$threshold_subtitle,
    covariates       = model_covariates
  )
}

fit_gamm_continuous_panel <- function(data, outcome_var, spec,
                                      covariates = figure4_core_covariates) {
  baseline_var     <- spec$baseline_var %||% spec$var
  model_covariates <- unique(covariates)

  model_df <- .gamm_prep_df(data, outcome_var, baseline_var, model_covariates) %>%
    filter(!is.na(.data[[outcome_var]]), !is.na(.data[[baseline_var]]),
           !is.na(time_yrs), is.finite(time_yrs))

  if (length(model_covariates))
    model_df <- model_df %>% filter(if_all(all_of(model_covariates), ~ !is.na(.x)))

  if (nrow(model_df) < 50 || n_distinct(model_df$ID) < 20 ||
      n_distinct(model_df[[baseline_var]]) < 8) return(NULL)

  rhs_smooth   <- paste(c("s(time_yrs, bs='cr', k=5)",
                           paste0("s(", baseline_var, ", bs='cr', k=5)"),
                           paste0("ti(time_yrs,", baseline_var, ", bs=c('cr','cr'), k=c(5,4))"),
                           model_covariates), collapse = " + ")
  rhs_fallback <- paste(c(paste0(baseline_var, " * time_yrs"),
                           model_covariates), collapse = " + ")
  smooth_formula   <- as.formula(paste(outcome_var, "~", rhs_smooth))
  fallback_formula <- as.formula(paste(outcome_var, "~", rhs_fallback))

  res         <- .gamm_fit_with_fallback(smooth_formula, fallback_formula, model_df)
  fit         <- res$fit
  fit_gam     <- .gamm_gam(fit)
  fit_summary <- summary(fit_gam)
  output_prefix <- paste0("../2_Output/Stats_GAMM_", .out_stub(outcome_var),
                          "_", spec$suffix, "Continuous")
  .write_gamm_summary(fit_summary, output_prefix)

  quartile_lookup <- model_df %>%
    distinct(ID, baseline_value = .data[[baseline_var]]) %>%
    mutate(
      quartile_num = ntile(baseline_value, 2),
      quartile = factor(
        quartile_num,
        levels = 1:2,
        labels = c("Lower half", "Upper half")
      )
    )
  quartile_summary <- quartile_lookup %>%
    group_by(quartile) %>%
    summarise(representative_value = median(baseline_value, na.rm = TRUE),
              q_min = min(baseline_value, na.rm = TRUE),
              q_max = max(baseline_value, na.rm = TRUE),
              n_patients = dplyr::n(), .groups = "drop")

  pred_df <- expand_grid(
    time_yrs = seq(0, min(10, max(model_df$time_yrs, na.rm = TRUE)), by = 0.25),
    quartile = quartile_summary$quartile
  ) %>%
    left_join(quartile_summary %>% select(quartile, representative_value), by = "quartile") %>%
    mutate(ID_fac = model_df$ID_fac[1])
  pred_df[[baseline_var]] <- pred_df$representative_value
  pred_df <- pred_df %>%
    fill_pred_covariates(model_df, model_covariates) %>%
    .gamm_predict_ci(fit, .)

  list(
    fit = fit, fit_mode = res$fit_mode,
    smooth_terms     = as_tibble(fit_summary$s.table, rownames = "term"),
    parametric_terms = as_tibble(fit_summary$p.table, rownames = "term"),
    pred_df          = pred_df,
    quartile_summary = quartile_summary,
    n_observations   = nrow(model_df),
    n_patients       = n_distinct(model_df$ID),
    outcome          = outcome_var,
    spec             = spec,
    covariates       = model_covariates
  )
}

fit_gamm_lvdd_panel <- function(data, outcome_var, covariates = figure4_core_covariates) {
  model_covariates <- unique(covariates)

  model_df <- .gamm_prep_df(data, outcome_var, character(0), model_covariates) %>%
    mutate(lvdd_burden = factor(as.character(lvdd_burden), levels = c("1", "2", "3")),
           # Ordered factor for difference-smooth varying-coefficient
           lvdd_burden_ord = ordered(lvdd_burden, levels = c("1", "2", "3"))) %>%
    filter(!is.na(.data[[outcome_var]]), !is.na(lvdd_burden),
           !is.na(time_yrs), is.finite(time_yrs))

  if (length(model_covariates))
    model_df <- model_df %>% filter(if_all(all_of(model_covariates), ~ !is.na(.x)))

  if (nrow(model_df) < 50 || n_distinct(model_df$ID) < 20 ||
      n_distinct(model_df$lvdd_burden) < 2) return(NULL)

  rhs_smooth   <- paste(c("s(time_yrs, bs='cr', k=5)",
                           "s(time_yrs, by=lvdd_burden_ord, bs='cr', k=5)",
                           "lvdd_burden",
                           model_covariates), collapse = " + ")
  rhs_fallback <- paste(c("lvdd_burden * time_yrs",
                           model_covariates), collapse = " + ")
  smooth_formula   <- as.formula(paste(outcome_var, "~", rhs_smooth))
  fallback_formula <- as.formula(paste(outcome_var, "~", rhs_fallback))

  res         <- .gamm_fit_with_fallback(smooth_formula, fallback_formula, model_df)
  fit         <- res$fit
  fit_gam     <- .gamm_gam(fit)
  fit_summary <- summary(fit_gam)
  output_prefix <- paste0("../2_Output/Stats_GAMM_", .out_stub(outcome_var), "_LVDDBurden")
  .write_gamm_summary(fit_summary, output_prefix)

  pred_df <- expand_grid(
    time_yrs    = seq(0, min(10, max(model_df$time_yrs, na.rm = TRUE)), by = 0.25),
    lvdd_burden = factor(c("1", "2", "3"), levels = c("1", "2", "3"))
  ) %>%
    mutate(lvdd_burden_ord = ordered(lvdd_burden, levels = c("1", "2", "3")),
           lvdd_group = factor(unname(figure4_lvdd_labels[as.character(lvdd_burden)]),
                               levels = unname(figure4_lvdd_labels)),
           ID_fac = model_df$ID_fac[1]) %>%
    fill_pred_covariates(model_df, model_covariates) %>%
    .gamm_predict_ci(fit, .)

  smooth_tbl_lvdd <- as_tibble(fit_summary$s.table, rownames = "term") %>%
    mutate(term = case_when(
      grepl("^s\\(time_yrs\\):lvdd_burden_ord2$", term) ~ "lvdd_burden2:time_yrs",
      grepl("^s\\(time_yrs\\):lvdd_burden_ord3$", term) ~ "lvdd_burden3:time_yrs",
      TRUE ~ term
    ))

  list(
    fit = fit, fit_mode = res$fit_mode,
    smooth_terms      = smooth_tbl_lvdd,
    parametric_terms  = as_tibble(fit_summary$p.table, rownames = "term"),
    pred_df           = pred_df,
    n_observations    = nrow(model_df),
    n_patients        = n_distinct(model_df$ID),
    n_group_patients  = model_df %>% distinct(ID, lvdd_burden) %>% count(lvdd_burden, name = "n_patients"),
    covariates        = model_covariates
  )
}

# ── Model term extraction ─────────────────────────────────────────────────────
extract_gamm_smooth_p <- function(result, terms) {
  if (is.null(result) || !nrow(result$smooth_terms)) return(NA_real_)
  df_resid  <- .gamm_dfresid(result$fit)
  term_row  <- result$smooth_terms %>% filter(term %in% terms) %>% slice(1)
  if (!nrow(term_row)) return(NA_real_)
  raw_p <- term_row$`p-value`[[1]]
  if (!is.na(raw_p) && (raw_p > 0 || is.na(df_resid))) return(raw_p)
  stats::pf(term_row$F[[1]], df1 = pmax(term_row$Ref.df[[1]], 1e-8),
            df2 = pmax(df_resid, 1e-8), lower.tail = FALSE)
}

extract_gamm_parametric_p <- function(result, terms) {
  if (is.null(result) || !nrow(result$parametric_terms)) return(NA_real_)
  df_resid <- .gamm_dfresid(result$fit)
  term_row <- result$parametric_terms %>% filter(term %in% terms) %>% slice(1)
  if (!nrow(term_row)) return(NA_real_)
  raw_p <- term_row$`Pr(>|t|)`[[1]]
  if (!is.na(raw_p) && (raw_p > 0 || is.na(df_resid))) return(raw_p)
  2 * stats::pt(abs(term_row$`t value`[[1]]), df = df_resid, lower.tail = FALSE)
}

extract_binary_model_rows <- function(result, outcome_label) {
  if (is.null(result)) return(NULL)
  df_resid <- .gamm_dfresid(result$fit)

  smooth_rows <- result$smooth_terms %>%
    filter(term %in% c("s(time_yrs)", "ti(time_yrs,status_abnormal)")) %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long,
      Observations = result$n_observations, Patients = result$n_patients,
      `Model component` = case_when(
        term == "s(time_yrs)"                     ~ "Time smooth",
        term == "ti(time_yrs,status_abnormal)"    ~
          paste0("Time × abnormal ", result$spec$label_short, " status"),
        TRUE ~ term
      ),
      edf = round(edf, 2), `F statistic` = round(F, 2),
      p_num = ifelse(is.na(`p-value`) | `p-value` > 0 | is.na(df_resid), `p-value`,
                     stats::pf(F, df1 = pmax(Ref.df, 1e-8), df2 = pmax(df_resid, 1e-8),
                               lower.tail = FALSE)),
      `P value` = as.character(ifelse(`p-value` < 0.001, "<0.001", sprintf("%.3f", `p-value`))),
      Interpretation = as.character(case_when(
        term == "ti(time_yrs,status_abnormal)" & p_num < 0.05 ~ "Trajectory differs by abnormal status",
        term == "ti(time_yrs,status_abnormal)"                 ~ "No clear trajectory modification by abnormal status",
        term == "s(time_yrs)"          & p_num < 0.05 ~ "Nonlinear temporal trend present",
        TRUE                                           ~ "Not statistically significant"
      ))
    )

  status_rows <- result$parametric_terms %>%
    filter(term == "status_binaryAbnormal") %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long,
      Observations = result$n_observations, Patients = result$n_patients,
      `Model component` = paste0("Abnormal ", result$spec$label_short, " status"),
      edf = NA_real_, `F statistic` = round(`t value`^2, 2),
      p_num = ifelse(is.na(`Pr(>|t|)`) | `Pr(>|t|)` > 0 | is.na(df_resid), `Pr(>|t|)`,
                     2 * stats::pt(abs(`t value`), df = df_resid, lower.tail = FALSE)),
      `P value` = as.character(ifelse(`Pr(>|t|)` < 0.001, "<0.001", sprintf("%.3f", `Pr(>|t|)`))),
      Interpretation = as.character(ifelse(p_num < 0.05, "Baseline abnormal-status association present",
                              "Not statistically significant"))
    )

  fallback_rows <- result$parametric_terms %>%
    filter(term %in% c("status_binaryAbnormal:time_yrs", "time_yrs:status_binaryAbnormal")) %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long,
      Observations = result$n_observations, Patients = result$n_patients,
      `Model component` = paste0("Time × abnormal ", result$spec$label_short,
                                 " status (linear fallback)"),
      edf = NA_real_, `F statistic` = round(`t value`^2, 2),
      p_num = ifelse(is.na(`Pr(>|t|)`) | `Pr(>|t|)` > 0 | is.na(df_resid), `Pr(>|t|)`,
                     2 * stats::pt(abs(`t value`), df = df_resid, lower.tail = FALSE)),
      `P value` = as.character(ifelse(`Pr(>|t|)` < 0.001, "<0.001", sprintf("%.3f", `Pr(>|t|)`))),
      Interpretation = as.character(ifelse(p_num < 0.05, "Trajectory differs by abnormal status",
                              "No clear trajectory modification by abnormal status"))
    )

  bind_rows(smooth_rows, status_rows, fallback_rows) %>%
    mutate(term_order = case_when(`Model component` == "Time smooth" ~ 1L,
                                  str_detect(`Model component`, "^Abnormal ") ~ 2L,
                                  TRUE ~ 3L)) %>%
    arrange(term_order) %>% select(-term_order)
}

extract_continuous_model_rows <- function(result, outcome_label) {
  if (is.null(result)) return(NULL)
  df_resid    <- .gamm_dfresid(result$fit)
  baseline_var <- result$spec$baseline_var %||% result$spec$var
  s_base      <- paste0("s(", baseline_var, ")")
  ti_terms    <- c(paste0("ti(time_yrs,", baseline_var, ")"),
                   paste0("ti(", baseline_var, ",time_yrs)"))

  smooth_rows <- result$smooth_terms %>%
    filter(term %in% c("s(time_yrs)", s_base, ti_terms)) %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long,
      `Model component` = case_when(
        term == "s(time_yrs)" ~ "Time smooth",
        term == s_base        ~ paste0("Baseline ", result$spec$label_short, " smooth"),
        term %in% ti_terms    ~ paste0("Time × baseline ", result$spec$label_short, " interaction"),
        TRUE ~ term
      ),
      edf = round(edf, 2), Statistic = sprintf("F = %.2f", F),
      p_num = ifelse(is.na(`p-value`) | `p-value` > 0 | is.na(df_resid), `p-value`,
                     stats::pf(F, df1 = pmax(Ref.df, 1e-8), df2 = pmax(df_resid, 1e-8),
                               lower.tail = FALSE))
    )

  fallback_rows <- result$parametric_terms %>%
    filter(term %in% c(baseline_var, paste0(baseline_var, ":time_yrs"),
                        paste0("time_yrs:", baseline_var))) %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long,
      `Model component` = case_when(
        term == baseline_var ~ paste0("Baseline ", result$spec$label_short, " (linear fallback)"),
        TRUE ~ paste0("Time × baseline ", result$spec$label_short, " interaction (linear fallback)")
      ),
      edf = NA_real_, Statistic = sprintf("t = %.2f", `t value`),
      p_num = ifelse(is.na(`Pr(>|t|)`) | `Pr(>|t|)` > 0 | is.na(df_resid), `Pr(>|t|)`,
                     2 * stats::pt(abs(`t value`), df = df_resid, lower.tail = FALSE))
    )

  bind_rows(smooth_rows, fallback_rows) %>%
    mutate(`P value` = format_p_mixed(p_num, sci_digits = 2, threshold = 0.05, decimal_digits = 3),
           term_order = case_when(`Model component` == "Time smooth" ~ 1L,
                                  str_detect(`Model component`, "^Baseline ") ~ 2L,
                                  TRUE ~ 3L)) %>%
    arrange(term_order) %>% select(-term_order)
}

label_figure4_term <- function(term, spec) {
  case_when(
    term == "(Intercept)"    ~ "Intercept",
    term == "age"            ~ "Age",
    term == "SexFemale"      ~ "Female sex",
    term == "BMI"            ~ "BMI",
    term == "bb_any"         ~ "β-blocker use",
    term == "ndhp_ccb_any"   ~ "Non-DHP CCB use",
    term == "dm_pre_test"    ~ "Diabetes",
    term == "time_yrs"       ~ "Time (linear)",
    term == "status_binaryAbnormal" ~
      paste0("Abnormal ", spec$label_short, " status"),
    term %in% c("status_binaryAbnormal:time_yrs", "time_yrs:status_binaryAbnormal") ~
      paste0("Time × abnormal ", spec$label_short, " status (linear)"),
    term == "s(time_yrs)"   ~ "Time smooth",
    term == "ti(time_yrs,status_abnormal)" ~
      paste0("Time × abnormal ", spec$label_short, " status"),
    term == "s(ID_fac)"     ~ "Subject random intercept",
    TRUE ~ term
  )
}

extract_complete_gamm_rows <- function(result, outcome_label) {
  if (is.null(result)) return(NULL)
  fit_mode_label <- ifelse(result$fit_mode == "smooth", "Smooth GAMM", "Linear fallback GAMM")
  smooth_rows <- result$smooth_terms %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long,
      `Fit mode` = fit_mode_label,
      Observations = result$n_observations, Patients = result$n_patients,
      `Term class` = "Smooth", Term = label_figure4_term(term, result$spec),
      Estimate = NA_real_, `Std. Error` = NA_real_,
      edf = edf, `Ref df` = Ref.df, `Statistic type` = "F", `Statistic value` = F,
      `P value raw` = `p-value`
    )
  parametric_rows <- result$parametric_terms %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long,
      `Fit mode` = fit_mode_label,
      Observations = result$n_observations, Patients = result$n_patients,
      `Term class` = "Parametric", Term = label_figure4_term(term, result$spec),
      Estimate = Estimate, `Std. Error` = `Std. Error`,
      edf = NA_real_, `Ref df` = NA_real_, `Statistic type` = "t",
      `Statistic value` = `t value`, `P value raw` = `Pr(>|t|)`
    )
  bind_rows(smooth_rows, parametric_rows)
}

extract_gamm_sensitivity_rows <- function(result, outcome_label, model_label) {
  if (is.null(result)) return(NULL)
  fit_mode_label <- ifelse(result$fit_mode == "smooth", "Smooth GAMM", "Linear fallback GAMM")
  status_rows <- result$parametric_terms %>%
    filter(term == "status_binaryAbnormal") %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long, Model = model_label,
      Sample = sprintf("%d / %d", result$n_observations, result$n_patients),
      Fit = fit_mode_label,
      `Relationship term` = paste0("Abnormal ", result$spec$label_short, " status"),
      `Effect summary` = sprintf("β = %.2f; SE = %.2f", Estimate, `Std. Error`),
      Statistic = sprintf("t = %.2f", `t value`),
      p_num = `Pr(>|t|)`,
      `P value` = as.character(ifelse(`Pr(>|t|)` < 0.001, "<0.001", sprintf("%.3f", `Pr(>|t|)`)))
    )
  smooth_interaction_rows <- result$smooth_terms %>%
    filter(term == "ti(time_yrs,status_abnormal)") %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long, Model = model_label,
      Sample = sprintf("%d / %d", result$n_observations, result$n_patients),
      Fit = fit_mode_label,
      `Relationship term` = paste0("Time × abnormal ", result$spec$label_short, " status"),
      `Effect summary` = sprintf("edf = %.2f; ref df = %.2f", edf, Ref.df),
      Statistic = sprintf("F = %.2f", F),
      p_num = `p-value`,
      `P value` = as.character(ifelse(`p-value` < 0.001, "<0.001", sprintf("%.3f", `p-value`)))
    )
  fallback_interaction_rows <- result$parametric_terms %>%
    filter(term %in% c("status_binaryAbnormal:time_yrs", "time_yrs:status_binaryAbnormal")) %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long, Model = model_label,
      Sample = sprintf("%d / %d", result$n_observations, result$n_patients),
      Fit = fit_mode_label,
      `Relationship term` = paste0("Time × abnormal ", result$spec$label_short, " status (linear)"),
      `Effect summary` = sprintf("β = %.2f; SE = %.2f", Estimate, `Std. Error`),
      Statistic = sprintf("t = %.2f", `t value`),
      p_num = `Pr(>|t|)`,
      `P value` = as.character(ifelse(`Pr(>|t|)` < 0.001, "<0.001", sprintf("%.3f", `Pr(>|t|)`)))
    )
  bind_rows(status_rows, smooth_interaction_rows, fallback_interaction_rows)
}

extract_gamm_covariate_rows <- function(result, outcome_label,
                                        model_label = "Age, sex, BMI adjusted") {
  if (is.null(result)) return(NULL)
  result$parametric_terms %>%
    filter(term %in% c("age", "SexFemale", "BMI", "bb_any", "ndhp_ccb_any", "dm_pre_test")) %>%
    transmute(
      Outcome = outcome_label, Parameter = result$spec$label_long, Model = model_label,
      Sample = sprintf("%d / %d", result$n_observations, result$n_patients),
      Covariate = label_figure4_term(term, result$spec),
      Estimate = Estimate, `Std. Error` = `Std. Error`,
      Statistic = sprintf("t = %.2f", `t value`),
      p_num = `Pr(>|t|)`,
      `P value` = as.character(ifelse(`Pr(>|t|)` < 0.001, "<0.001", sprintf("%.3f", `Pr(>|t|)`)))
    )
}

extract_lvdd_model_rows <- function(result, outcome_label) {
  if (is.null(result)) return(NULL)
  df_resid <- .gamm_dfresid(result$fit)

  smooth_rows <- result$smooth_terms %>%
    filter(term %in% c("s(time_yrs)", "ti(time_yrs,lvdd_burden)", "ti(lvdd_burden,time_yrs)")) %>%
    transmute(
      Outcome = outcome_label,
      `Model component` = case_when(
        term == "s(time_yrs)" ~ "Time smooth",
        term %in% c("ti(time_yrs,lvdd_burden)", "ti(lvdd_burden,time_yrs)") ~
          "Time × LVDD burden group",
        TRUE ~ term
      ),
      edf = round(edf, 2), Statistic = sprintf("F = %.2f", F),
      p_num = ifelse(is.na(`p-value`) | `p-value` > 0 | is.na(df_resid), `p-value`,
                     stats::pf(F, df1 = pmax(Ref.df, 1e-8), df2 = pmax(df_resid, 1e-8),
                               lower.tail = FALSE))
    )
  burden_rows <- result$parametric_terms %>%
    filter(term %in% c("lvdd_burden2", "lvdd_burden3")) %>%
    transmute(
      Outcome = outcome_label,
      `Model component` = case_when(
        term == "lvdd_burden2" ~ "2 abnormal parameters vs 1",
        term == "lvdd_burden3" ~ "All 3 abnormal vs 1",
        TRUE ~ term
      ),
      edf = NA_real_, Statistic = sprintf("t = %.2f", `t value`),
      p_num = ifelse(is.na(`Pr(>|t|)`) | `Pr(>|t|)` > 0 | is.na(df_resid), `Pr(>|t|)`,
                     2 * stats::pt(abs(`t value`), df = df_resid, lower.tail = FALSE))
    )
  fallback_rows <- result$parametric_terms %>%
    filter(term %in% c("lvdd_burden2:time_yrs", "time_yrs:lvdd_burden2",
                        "lvdd_burden3:time_yrs", "time_yrs:lvdd_burden3")) %>%
    transmute(
      Outcome = outcome_label,
      `Model component` = case_when(
        term %in% c("lvdd_burden2:time_yrs", "time_yrs:lvdd_burden2") ~
          "Time × 2 abnormal parameters vs 1 (linear)",
        TRUE ~ "Time × all 3 abnormal vs 1 (linear)"
      ),
      edf = NA_real_, Statistic = sprintf("t = %.2f", `t value`),
      p_num = ifelse(is.na(`Pr(>|t|)`) | `Pr(>|t|)` > 0 | is.na(df_resid), `Pr(>|t|)`,
                     2 * stats::pt(abs(`t value`), df = df_resid, lower.tail = FALSE))
    )
  bind_rows(smooth_rows, burden_rows, fallback_rows) %>%
    mutate(`P value` = format_p_mixed(p_num, sci_digits = 2, threshold = 0.05, decimal_digits = 3),
           term_order = case_when(
             `Model component` == "Time smooth"              ~ 1L,
             `Model component` == "2 abnormal parameters vs 1" ~ 2L,
             `Model component` == "All 3 abnormal vs 1"     ~ 3L,
             TRUE ~ 4L
           )) %>%
    arrange(term_order) %>% select(-term_order)
}

collect_prediction_limits <- function(result_list) {
  pred_tbl <- map_dfr(result_list, \(x) if (is.null(x)) NULL else x$pred_df)
  if (!nrow(pred_tbl)) return(c(0, 1))
  range(c(pred_tbl$lo, pred_tbl$hi), na.rm = TRUE)
}

# ── Stats-label formatters ────────────────────────────────────────────────────
format_plot_p_html <- function(p, digits = 2, decimal_digits = 3, sci_cutoff = 0.001) {
  if (is.na(p) || p == 0) return(ifelse(is.na(p), "NA", "0"))
  if (p >= sci_cutoff) return(sprintf(paste0("%.", decimal_digits, "f"), p))
  sci <- scientific_p_parts(p, digits = digits)
  paste0(sci$mantissa, " × 10<sup>", sci$exponent, "</sup>")
}

format_continuous_panel_stats_label <- function(result) {
  if (is.null(result)) return(NULL)
  baseline_var <- result$spec$baseline_var %||% result$spec$var
  s_base  <- paste0("s(", baseline_var, ")")
  ti_fwd  <- paste0("ti(time_yrs,", baseline_var, ")")
  ti_rev  <- paste0("ti(", baseline_var, ",time_yrs)")
  fb_fwd  <- paste0(baseline_var, ":time_yrs")
  fb_rev  <- paste0("time_yrs:", baseline_var)
  time_p        <- extract_gamm_smooth_p(result, "s(time_yrs)")
  baseline_p    <- extract_gamm_smooth_p(result, s_base)
  interaction_p <- extract_gamm_smooth_p(result, c(ti_fwd, ti_rev))
  if (is.na(baseline_p))    baseline_p    <- extract_gamm_parametric_p(result, baseline_var)
  if (is.na(interaction_p)) interaction_p <- extract_gamm_parametric_p(result, c(fb_fwd, fb_rev))
  paste0("Obs/patients: ", result$n_observations, "/", result$n_patients,
         "<br>Time P = ", format_plot_p_html(time_p),
         "<br>Baseline P = ", format_plot_p_html(baseline_p),
         "<br>Interaction P = ", format_plot_p_html(interaction_p))
}

format_lvdd_panel_stats_label <- function(result) {
  if (is.null(result)) return(NULL)
  time_p   <- extract_gamm_smooth_p(result,     "s(time_yrs)")
  burden2_p <- extract_gamm_parametric_p(result, "lvdd_burden2")
  burden3_p <- extract_gamm_parametric_p(result, "lvdd_burden3")
  smooth_p  <- extract_gamm_smooth_p(result,    c("ti(time_yrs,lvdd_burden)", "ti(lvdd_burden,time_yrs)"))
  int2_p    <- extract_gamm_parametric_p(result, c("lvdd_burden2:time_yrs", "time_yrs:lvdd_burden2"))
  int3_p    <- extract_gamm_parametric_p(result, c("lvdd_burden3:time_yrs", "time_yrs:lvdd_burden3"))
  lines <- c(
    paste0("Obs/patients: ", result$n_observations, "/", result$n_patients),
    paste0("Time P = ",   format_plot_p_html(time_p)),
    paste0("2 vs 1 P = ", format_plot_p_html(burden2_p)),
    paste0("3 vs 1 P = ", format_plot_p_html(burden3_p))
  )
  if (!is.na(smooth_p)) lines <- c(lines, paste0("Interaction P = ", format_plot_p_html(smooth_p)))
  else lines <- c(lines, paste0("2×Time P = ", format_plot_p_html(int2_p)),
                          paste0("3×Time P = ", format_plot_p_html(int3_p)))
  paste(lines, collapse = "<br>")
}

# ── Panel-building functions ──────────────────────────────────────────────────
empty_gamm_panel <- function(label) {
  ggplot() + annotate("text", x = 0.5, y = 0.5, label = label, size = 3.3) + theme_void()
}

extract_legend_grob <- function(plot_obj) {
  gt      <- ggplotGrob(plot_obj)
  grob_idx <- which(vapply(gt$grobs, \(x) if (is.null(x$name)) "" else x$name,
                           character(1)) == "guide-box")
  if (!length(grob_idx)) return(grid::nullGrob())
  gt$grobs[[grob_idx[1]]]
}

make_binary_trajectory_panel <- function(result, y_label, y_limits,
                                         show_y = FALSE, show_x = FALSE,
                                         show_legend = FALSE) {
  if (is.null(result)) return(empty_gamm_panel("Model unavailable"))
  status_display_labels <- c(
    "Normal" = result$spec$status_reference_label %||% "At/below threshold",
    "Abnormal" = result$spec$status_high_label %||% "Above threshold"
  )
  p <- ggplot(result$pred_df, aes(x = time_yrs, y = fit,
                                   color = status_binary, fill = status_binary)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(
      values = figure4_status_palette,
      labels = status_display_labels,
      drop = FALSE
    ) +
    scale_fill_manual(
      values = figure4_status_palette,
      labels = status_display_labels,
      drop = FALSE
    ) +
    scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.03, 0.08)),
                       labels = label_number(accuracy = 1)) +
    labs(title    = result$spec$label_short,
         subtitle = result$status_caption,
         x = if (show_x) "Years from Baseline CPET" else NULL,
         y = y_label, color = "Baseline threshold") +
    theme_jacc(base_size = 8.1) +
    theme(legend.position  = if (show_legend) "bottom" else "none",
          legend.justification = c(0.5, 0), legend.direction = "horizontal",
          legend.background = element_rect(fill = alpha("white", 0.90),
                                           color = "grey70", linewidth = 0.25),
          legend.key.height = unit(0.24, "cm"),
          legend.text  = ggtext::element_markdown(size = 6.0),
          legend.title = ggtext::element_markdown(size = 6.2),
          legend.margin = margin(2, 3, 2, 3),
          plot.title    = ggtext::element_markdown(size = 8.7, face = "bold"),
          plot.subtitle = ggtext::element_markdown(size = 6.2, color = "grey30"),
          axis.title    = ggtext::element_markdown(size = 7.5),
          axis.text     = element_text(size = 6.9))
  if (!show_x) p <- p + theme(axis.title.x = element_blank(),
                               axis.text.x  = element_blank(),
                               axis.ticks.x = element_blank())
  if (!show_y) p <- p + theme(axis.title.y = element_blank(),
                               axis.text.y  = element_blank(),
                               axis.ticks.y = element_blank())
  p
}

make_continuous_trajectory_panel <- function(result, y_label, y_limits, title,
                                             show_y = FALSE, show_x = FALSE,
                                             show_legend = FALSE, stats_label = NULL,
                                             legend_position = "bottom",
                                             legend_direction = "horizontal",
                                             legend_justification = NULL,
                                             stats_fill_alpha = 0.9) {
  if (is.null(result)) return(empty_gamm_panel("Model unavailable"))
  computed_lj <- legend_justification %||%
    (if (is.numeric(legend_position)) c(0, 0)
     else if (identical(legend_position, "right")) c(0, 0.5)
     else c(0.5, 0))
  p <- ggplot(result$pred_df, aes(x = time_yrs, y = fit, color = quartile, fill = quartile)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.14, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = figure4_quartile_palette, drop = FALSE) +
    scale_fill_manual(values  = figure4_quartile_palette, drop = FALSE) +
    scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.03, 0.08)),
                       labels = label_number(accuracy = 1)) +
    labs(title = title,
         x = if (show_x) "Years from Baseline CPET" else NULL,
         y = y_label, color = "Baseline half") +
    theme_jacc(base_size = 8.1) +
    theme(legend.position      = if (show_legend) legend_position else "none",
          legend.justification = computed_lj, legend.direction = legend_direction,
          legend.background = element_rect(fill = alpha("white", 0.90),
                                           color = "grey70", linewidth = 0.25),
          legend.key.width  = unit(0.48, "cm"), legend.key.height = unit(0.24, "cm"),
          legend.text  = ggtext::element_markdown(size = 6.0),
          legend.title = ggtext::element_markdown(size = 6.2),
          legend.margin = margin(2, 3, 2, 3),
          plot.title    = element_text(size = 8.7, face = "bold", hjust = 0),
          plot.subtitle = ggtext::element_markdown(size = 6.2, color = "grey30"),
          axis.title    = element_text(size = 7.5),
          axis.text     = element_text(size = 6.9))
  p <- add_plot_stats_box(p, stats_label = stats_label,
                          x_range = range(result$pred_df$time_yrs, na.rm = TRUE),
                          y_range = range(y_limits, na.rm = TRUE),
                          position = "upper_left", size = 2.25,
                          fill_alpha = stats_fill_alpha)
  if (!show_x) p <- p + theme(axis.title.x = element_blank(),
                               axis.text.x  = element_blank(),
                               axis.ticks.x = element_blank())
  if (!show_y) p <- p + theme(axis.title.y = element_blank(),
                               axis.text.y  = element_blank(),
                               axis.ticks.y = element_blank())
  p
}

make_lvdd_trajectory_panel <- function(result, y_label, y_limits, title,
                                       show_y = FALSE, show_x = FALSE,
                                       show_legend = FALSE, stats_label = NULL) {
  if (is.null(result)) return(empty_gamm_panel("Model unavailable"))
  p <- ggplot(result$pred_df, aes(x = time_yrs, y = fit,
                                   color = lvdd_group, fill = lvdd_group)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, linewidth = 0, show.legend = FALSE) +
    geom_line(linewidth = 1.0) +
    scale_color_manual(values = figure4_lvdd_palette, drop = FALSE) +
    scale_fill_manual(values  = figure4_lvdd_palette, drop = FALSE) +
    scale_x_continuous(limits = c(0, 10), breaks = seq(0, 10, by = 2),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.03, 0.08)),
                       labels = label_number(accuracy = 1)) +
    labs(title    = title,
         subtitle = "Baseline LVDD burden: 1, 2, or all 3 abnormal primary parameters",
         x = if (show_x) "Years from Baseline CPET" else NULL,
         y = y_label, color = "Baseline LVDD\nburden") +
    theme_jacc(base_size = 8.1) +
    theme(legend.position      = if (show_legend) "bottom" else "none",
          legend.justification = c(0.5, 0), legend.direction = "horizontal",
          legend.background = element_rect(fill = alpha("white", 0.90),
                                           color = "grey70", linewidth = 0.25),
          legend.key.height = unit(0.24, "cm"),
          legend.text  = ggtext::element_markdown(size = 6.0),
          legend.title = ggtext::element_markdown(size = 6.2),
          legend.margin = margin(2, 3, 2, 3),
          plot.title    = ggtext::element_markdown(size = 8.7, face = "bold"),
          plot.subtitle = ggtext::element_markdown(size = 6.2, color = "grey30"),
          axis.title    = ggtext::element_markdown(size = 7.5),
          axis.text     = element_text(size = 6.9))
  p <- add_plot_stats_box(p, stats_label = stats_label,
                          x_range = range(result$pred_df$time_yrs, na.rm = TRUE),
                          y_range = range(y_limits, na.rm = TRUE),
                          position = "upper_left", size = 2.35)
  if (!show_x) p <- p + theme(axis.title.x = element_blank(),
                               axis.text.x  = element_blank(),
                               axis.ticks.x = element_blank())
  if (!show_y) p <- p + theme(axis.title.y = element_blank(),
                               axis.text.y  = element_blank(),
                               axis.ticks.y = element_blank())
  p
}

# ── Binary subfigure assembly ─────────────────────────────────────────────────
build_binary_subfigure_summary <- function(result, outcome_label) {
  if (is.null(result)) return(NULL)
  extract_binary_model_rows(result, outcome_label) %>%
    transmute(
      Outcome = outcome_label,
      `Model component` = case_when(
        str_detect(`Model component`, regex("^Time.*abnormal", ignore_case = TRUE)) ~
          "Time x threshold group",
        str_detect(`Model component`, regex("abnormal", ignore_case = TRUE)) ~
          "Threshold group",
        TRUE ~ `Model component`
      ),
              edf = ifelse(is.na(edf), "—", sprintf("%.2f", edf)),
              `F statistic` = sprintf("%.2f", `F statistic`), `P value`)
}

make_parameter_table_grob <- function(df) {
  headers        <- c("Model component", "edf", "F statistic", "P value")
  outcome_levels <- unique(df$Outcome)

  layout_rows <- map_dfr(outcome_levels, function(ol) {
    bind_rows(
      tibble(row_type = "group", Outcome = ol, `Model component` = ol,
             edf = "", `F statistic` = "", `P value` = ""),
      df %>% filter(Outcome == ol) %>% mutate(row_type = "body")
    )
  })
  table_df   <- layout_rows
  row_heights <- ifelse(table_df$row_type == "group", 0.92, 1.08)

  panel_layout <- grid::grid.layout(
    nrow = nrow(table_df) + 1, ncol = 4,
    widths  = grid::unit(c(0.50, 0.12, 0.22, 0.16), "npc"),
    heights = grid::unit(c(1.18, row_heights), "null")
  )

  grobs <- list(
    grid::rectGrob(gp = grid::gpar(fill = "white", col = "#1F2328", lwd = 0.8),
                   vp = grid::viewport(layout = panel_layout,
                                       layout.pos.row = 1:(nrow(table_df) + 1),
                                       layout.pos.col = 1:4)),
    grid::segmentsGrob(x0 = grid::unit(0.02, "npc"), x1 = grid::unit(0.98, "npc"),
                       y0 = grid::unit(0.06, "npc"), y1 = grid::unit(0.06, "npc"),
                       gp = grid::gpar(col = "#1F2328", lwd = 0.9),
                       vp = grid::viewport(layout = panel_layout,
                                           layout.pos.col = 1:4, layout.pos.row = 1))
  )
  for (col_idx in seq_along(headers)) {
    grobs[[length(grobs) + 1]] <- grid::textGrob(
      headers[col_idx],
      x = if (col_idx == 1) grid::unit(0.03, "npc") else grid::unit(0.5, "npc"),
      y = grid::unit(0.5, "npc"),
      just = if (col_idx == 1) c("left","center") else c("center","center"),
      gp = grid::gpar(fontface = "bold", fontsize = 6.8, fontfamily = "Times New Roman",
                      col = "#1F2328", lineheight = 0.95),
      vp = grid::viewport(layout = panel_layout, layout.pos.row = 1, layout.pos.col = col_idx)
    )
  }
  for (row_idx in seq_len(nrow(table_df))) {
    table_row    <- row_idx + 1
    is_group_row <- identical(table_df$row_type[[row_idx]], "group")
    if (is_group_row) {
      grobs[[length(grobs) + 1]] <- grid::rectGrob(
        gp = grid::gpar(fill = "#F4F1EC", col = NA),
        vp = grid::viewport(layout = panel_layout, layout.pos.row = table_row, layout.pos.col = 1:4))
      grobs[[length(grobs) + 1]] <- grid::textGrob(
        table_df$Outcome[[row_idx]], x = grid::unit(0.03, "npc"), y = grid::unit(0.5, "npc"),
        just = c("left","center"),
        gp = grid::gpar(fontface = "bold", fontsize = 7.0, fontfamily = "Times New Roman",
                        col = "#23313F"),
        vp = grid::viewport(layout = panel_layout, layout.pos.row = table_row, layout.pos.col = 1:4))
      grobs[[length(grobs) + 1]] <- grid::segmentsGrob(
        x0 = grid::unit(0.02, "npc"), x1 = grid::unit(0.98, "npc"),
        y0 = grid::unit(0.02, "npc"), y1 = grid::unit(0.02, "npc"),
        gp = grid::gpar(col = "#9AA4AE", lwd = 0.6),
        vp = grid::viewport(layout = panel_layout, layout.pos.col = 1:4, layout.pos.row = table_row))
      next
    }
    row_values <- table_df[row_idx, c("Model component","edf","F statistic","P value")]
    for (col_idx in seq_along(headers)) {
      grobs[[length(grobs) + 1]] <- grid::textGrob(
        row_values[[col_idx]],
        x = if (col_idx == 1) grid::unit(0.03, "npc") else grid::unit(0.5, "npc"),
        y = grid::unit(0.5, "npc"),
        just = if (col_idx == 1) c("left","center") else c("center","center"),
        gp = grid::gpar(fontsize = if (col_idx == 1) 6.65 else 6.9,
                        fontfamily = "Times New Roman", col = "#1F2328", lineheight = 0.98),
        vp = grid::viewport(layout = panel_layout, layout.pos.row = table_row, layout.pos.col = col_idx))
    }
    grobs[[length(grobs) + 1]] <- grid::segmentsGrob(
      x0 = grid::unit(0.02, "npc"), x1 = grid::unit(0.98, "npc"),
      y0 = grid::unit(0.02, "npc"), y1 = grid::unit(0.02, "npc"),
      gp = grid::gpar(col = "#D0D7DE", lwd = 0.35),
      vp = grid::viewport(layout = panel_layout, layout.pos.col = 1:4, layout.pos.row = table_row))
  }
  grid::grobTree(children = do.call(grid::gList, grobs),
                 vp = grid::viewport(layout = panel_layout))
}

build_binary_subfigure <- function(spec, result_source = figure4_results,
                                   vo2_limits = NULL, vevco2_limits = NULL) {
  vo2_result    <- result_source[["VO2_FRIEND2_PP"]][[spec$var]]
  vevco2_result <- result_source[["VeVco2_slope"]][[spec$var]]
  if (is.null(vo2_result) && is.null(vevco2_result)) return(NULL)

  vo2_limits    <- vo2_limits    %||% collect_prediction_limits(result_source[["VO2_FRIEND2_PP"]])
  vevco2_limits <- vevco2_limits %||% collect_prediction_limits(result_source[["VeVco2_slope"]])

  p_vo2 <- if (!is.null(vo2_result))
    make_binary_trajectory_panel(vo2_result, y_label = label_peak_vo2_friend_unicode,
                                 y_limits = vo2_limits, show_y = TRUE, show_x = TRUE,
                                 show_legend = TRUE)
  else empty_gamm_panel("Peak V̇O2 unavailable")

  p_vevco2 <- if (!is.null(vevco2_result))
    make_binary_trajectory_panel(vevco2_result, y_label = label_vevco2_title_unicode,
                                 y_limits = vevco2_limits, show_y = TRUE, show_x = TRUE,
                                 show_legend = FALSE)
  else empty_gamm_panel("V̇E/V̇CO2 unavailable")

  panel_summary  <- bind_rows(build_binary_subfigure_summary(vo2_result,    label_peak_vo2),
                               build_binary_subfigure_summary(vevco2_result, label_vevco2))
  table_grob      <- make_parameter_table_grob(panel_summary)
  subtitle_source <- vo2_result %||% vevco2_result
  status_subtitle <- paste0("Baseline ", subtitle_source$status_caption)

  p_vo2 <- p_vo2 +
    labs(title = paste0("(A) ", label_peak_vo2_unicode), subtitle = status_subtitle) +
    theme(plot.title    = ggtext::element_markdown(size = 8.5, face = "bold", hjust = 0),
          plot.subtitle = ggtext::element_markdown(size = 6.0, color = "grey30", hjust = 0),
          legend.position      = c(0.04, 0.98), legend.justification = c(0, 1),
          legend.background    = element_rect(fill = alpha("white", 0.92),
                                              color = "grey65", linewidth = 0.25),
          legend.key.height    = unit(0.18, "cm"),
          legend.text          = ggtext::element_markdown(size = 5.5),
          legend.title         = ggtext::element_markdown(size = 5.7),
          legend.margin        = margin(1, 2, 1, 2),
          plot.margin          = margin(7, 4, 4, 4))
  p_vevco2 <- p_vevco2 +
    labs(title = paste0("(B) ", label_vevco2_unicode), subtitle = status_subtitle) +
    theme(plot.title    = ggtext::element_markdown(size = 8.5, face = "bold", hjust = 0),
          plot.subtitle = ggtext::element_markdown(size = 6.0, color = "grey30", hjust = 0),
          plot.margin   = margin(7, 4, 4, 4))

  table_panel <- grid::grobTree(
    grid::textGrob("(C) Longitudinal Model Summary",
                   x = unit(0.03, "npc"), y = unit(0.93, "npc"), just = c("left","top"),
                   gp = grid::gpar(fontface = "bold", fontsize = 8.4,
                                   fontfamily = "Times New Roman", col = "#1F2328")),
    grid::grobTree(table_grob,
                   vp = grid::viewport(x = 0.52, y = 0.54, width = 0.92, height = 0.62,
                                       just = c("center","center")))
  )

  fig <- (p_vo2 | p_vevco2 | wrap_elements(full = table_panel)) +
    plot_layout(widths = c(1.02, 1.02, 1.20))

  list(fig = fig, summary = panel_summary)
}
