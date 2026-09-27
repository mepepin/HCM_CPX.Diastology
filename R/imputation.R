# ── Multiple imputation: labelled sensitivity analysis only ───────────────────
# Owner decision, card 1.2 of 8_Docs/DECISION_CARDS_2026-09.md. The primary
# analyses stay available-case. Multiple imputation answers the separate
# question of what the estimates would be if the exposures were missing at
# random given the covariates and the outcome.
#
# Three constraints come from that card and are enforced here:
#
#   * Only E/e' and LAVi are imputed. TRVmax is never imputed, because whether
#     a tricuspid jet is measurable is itself tied to the value being measured;
#     treating it as missing at random would manufacture information. It still
#     contributes as an auxiliary through its fully observed gated form
#     (a presence indicator and present * value), so no TRVmax value is ever
#     invented while its information is still used.
#
#   * The outcome belongs in the imputation model. For a survival outcome that
#     means the event indicator plus the Nelson-Aalen estimate of the
#     cumulative hazard rather than raw follow-up time (White & Royston 2009);
#     omitting the outcome biases the imputed exposure toward the null.
#
#   * m and the seed are fixed in configuration so the sensitivity analysis is
#     reproducible run to run.

# Nelson-Aalen cumulative hazard evaluated at each patient's follow-up time.
nelson_aalen_cumhaz <- function(time, event) {
  fit <- survival::survfit(survival::Surv(time, event) ~ 1)
  stats::approx(
    x = c(0, fit$time),
    y = c(0, fit$cumhaz),
    xout = time,
    method = "constant",
    rule = 2
  )$y
}

# Rubin's rules for one scalar estimate. Done explicitly rather than through a
# pooling helper so the degrees of freedom and the fraction of missing
# information are visible and auditable.
rubin_pool <- function(estimates, variances) {
  m <- length(estimates)
  stopifnot(m > 1, length(variances) == m, all(is.finite(estimates)), all(variances > 0))
  q_bar <- mean(estimates)
  u_bar <- mean(variances)
  b <- stats::var(estimates)
  total_var <- u_bar + (1 + 1 / m) * b
  # Barnard-Rubin adjusted df would need the complete-data df; the classical
  # Rubin df is used here and reported alongside the FMI.
  r <- (1 + 1 / m) * b / u_bar
  df <- (m - 1) * (1 + 1 / r)^2
  fmi <- (r + 2 / (df + 3)) / (r + 1)
  se <- sqrt(total_var)
  tibble::tibble(
    estimate = q_bar,
    std_error = se,
    df = df,
    conf_low = q_bar - stats::qt(0.975, df) * se,
    conf_high = q_bar + stats::qt(0.975, df) * se,
    p_value = 2 * stats::pt(-abs(q_bar) / se, df),
    relative_increase_variance = r,
    fraction_missing_information = fmi
  )
}

# Multiply imputes the named exposures and refits each one's primary Cox model,
# pooling across imputations. Returns one row per exposure.
mi_cox_exposure_sensitivity <- function(surv_df,
                                        exposures,
                                        covariates = c("age", "Sex", "BMI"),
                                        auxiliary_gated = "tr_max_vel",
                                        time_var = "follow_up_yrs",
                                        event_var = "hf_composite",
                                        m = 20,
                                        seed = 20260916) {
  impute_vars <- names(exposures)

  frame <- surv_df %>%
    dplyr::filter(
      !is.na(.data[[time_var]]),
      !is.na(.data[[event_var]]),
      dplyr::if_all(dplyr::all_of(covariates), ~ !is.na(.x))
    )

  # Auxiliary in fully observed gated form: never imputed, still informative.
  aux_build <- add_missing_indicator_covariates(
    frame,
    gated_numeric_vars = auxiliary_gated
  )
  frame <- aux_build$data
  aux_vars <- aux_build$covariates

  frame <- frame %>%
    dplyr::mutate(
      .na_cumhaz = nelson_aalen_cumhaz(.data[[time_var]], .data[[event_var]]),
      .event = as.integer(.data[[event_var]] == 1)
    )

  mice_vars <- c(impute_vars, covariates, aux_vars, ".na_cumhaz", ".event")
  mice_df <- frame %>%
    dplyr::select(dplyr::all_of(mice_vars)) %>%
    as.data.frame()

  # Everything except the two exposures must be complete, or mice would impute
  # it too and the analysis would quietly stop being the one that was agreed.
  non_exposure <- setdiff(mice_vars, impute_vars)
  incomplete_non_exposure <- non_exposure[
    vapply(non_exposure, function(v) any(is.na(mice_df[[v]])), logical(1))
  ]
  if (length(incomplete_non_exposure)) {
    stop(
      "mi_cox_exposure_sensitivity(): only the exposures may be imputed, but ",
      paste(incomplete_non_exposure, collapse = ", "),
      " also contain missing values.",
      call. = FALSE
    )
  }

  method <- rep("", length(mice_vars))
  names(method) <- mice_vars
  method[impute_vars] <- "pmm"

  imputed <- mice::mice(
    mice_df,
    m = m,
    method = method,
    seed = seed,
    printFlag = FALSE
  )

  purrr::map_dfr(impute_vars, function(exposure) {
    fits <- lapply(seq_len(m), function(i) {
      completed <- mice::complete(imputed, i)
      completed[[time_var]] <- frame[[time_var]]
      completed[[event_var]] <- frame[[event_var]]
      # Standardised within each completed dataset, matching the primary
      # models' per-sample z-scores.
      completed$.exposure_z <- as.numeric(scale(completed[[exposure]]))
      fit <- survival::coxph(
        stats::reformulate(
          c(".exposure_z", covariates),
          response = sprintf("survival::Surv(%s, %s)", time_var, event_var)
        ),
        data = completed,
        ties = "efron"
      )
      coefs <- summary(fit)$coefficients
      list(
        estimate = unname(coefs[".exposure_z", "coef"]),
        variance = unname(coefs[".exposure_z", "se(coef)"])^2
      )
    })
    pooled <- rubin_pool(
      vapply(fits, function(x) x$estimate, numeric(1)),
      vapply(fits, function(x) x$variance, numeric(1))
    )
    pooled %>%
      dplyr::transmute(
        Parameter = unname(exposures[[exposure]]),
        `Patients in imputation frame` = nrow(frame),
        Events = sum(frame[[event_var]] == 1, na.rm = TRUE),
        `Exposure measured` = sum(!is.na(frame[[exposure]])),
        `Exposure imputed` = sum(is.na(frame[[exposure]])),
        Imputations = m,
        HR = exp(estimate),
        `CI low` = exp(conf_low),
        `CI high` = exp(conf_high),
        `P value` = p_value,
        `Rubin df` = df,
        `Fraction missing information` = fraction_missing_information
      )
  })
}

# ── Substantive-model-compatible imputation for the spline models ─────────────
# The cross-sectional substantive model contains a restricted cubic spline of
# the exposure. Imputing that exposure from a model without its own nonlinear
# terms would pull the imputations toward a straight line and attenuate the
# nonlinearity test, so the imputation here is substantive-model compatible
# (smcfcs: Bartlett et al. 2015) rather than ordinary FCS.
#
# Three details matter for the result to mean anything:
#
#   * Knots are computed once on the observed data and then held fixed, so
#     every imputation is tested against the same basis as the primary model
#     rather than against a basis that moves with the imputed values.
#
#   * The exposure is imputed on the log scale and back-transformed as a
#     passive variable. Imputing it directly produced negative E/e' values,
#     which would have made the spline extrapolate below its boundary knot.
#
#   * The basis is Hmisc::rcspline.eval(inclx = TRUE), whose first column is
#     the exposure itself. The linear model is therefore a true parameter
#     subset of the spline model, which is what mice::D1 needs to test
#     nonlinearity; a natural-spline basis spans the same space but its
#     coefficients are not nested, and pooling silently returned the overall
#     test twice.
#
# The breakpoint bootstrap is deliberately not run here: it searches for a
# knot, and card 1.2 keeps that search on the available-case analysis only.
mi_rcs_exposure_sensitivity <- function(analysis_df,
                                        outcome_var,
                                        outcome_label,
                                        exposure,
                                        exposure_label,
                                        auxiliary_params,
                                        covariates = c("age", "Sex", "BMI"),
                                        nk = 4,
                                        m = 20,
                                        seed = 20260916,
                                        numit = 10,
                                        available_case = NULL) {
  log_var <- paste0("log_", exposure)

  frame <- analysis_df %>%
    dplyr::filter(
      !is.na(.data[[outcome_var]]),
      dplyr::if_all(dplyr::all_of(covariates), ~ !is.na(.x))
    )
  observed <- frame %>% dplyr::filter(!is.na(.data[[exposure]]))
  if (nrow(observed) < 20) return(tibble::tibble())

  if (any(observed[[exposure]] <= 0)) {
    stop(
      sprintf(
        "mi_rcs_exposure_sensitivity(): '%s' has non-positive values, so the log-scale imputation is invalid.",
        exposure
      ),
      call. = FALSE
    )
  }

  knots <- build_rcs_basis_df(observed[[exposure]], idx = exposure, nk = nk)$knots

  # Auxiliaries in gated form so they are never themselves imputed.
  aux_build <- add_missing_indicator_covariates(
    frame,
    gated_numeric_vars = auxiliary_params
  )
  aux_vars <- aux_build$covariates

  mi_df <- aux_build$data %>%
    dplyr::mutate(
      !!log_var := ifelse(
        is.na(.data[[exposure]]),
        NA_real_,
        log(.data[[exposure]])
      )
    ) %>%
    dplyr::select(dplyr::all_of(c(
      outcome_var, log_var, exposure, covariates, aux_vars
    ))) %>%
    as.data.frame()

  non_imputed <- setdiff(names(mi_df), c(log_var, exposure))
  incomplete <- non_imputed[
    vapply(non_imputed, function(v) any(is.na(mi_df[[v]])), logical(1))
  ]
  if (length(incomplete)) {
    stop(
      "mi_rcs_exposure_sensitivity(): only the exposure may be imputed, but ",
      paste(incomplete, collapse = ", "), " also contain missing values.",
      call. = FALSE
    )
  }

  method <- rep("", ncol(mi_df))
  names(method) <- names(mi_df)
  method[log_var] <- "norm"
  method[exposure] <- sprintf("exp(%s)", log_var)

  predictor_matrix <- matrix(
    0,
    nrow = ncol(mi_df),
    ncol = ncol(mi_df),
    dimnames = list(names(mi_df), names(mi_df))
  )
  predictor_matrix[log_var, c(covariates, aux_vars)] <- 1

  smformula <- sprintf(
    "%s ~ Hmisc::rcspline.eval(%s, knots = c(%s), inclx = TRUE) + %s",
    outcome_var,
    exposure,
    paste(sprintf("%.10f", knots), collapse = ", "),
    paste(covariates, collapse = " + ")
  )

  set.seed(seed)
  imputed <- smcfcs::smcfcs(
    mi_df,
    smtype = "lm",
    smformula = smformula,
    method = method,
    predictorMatrix = predictor_matrix,
    m = m,
    numit = numit,
    rjlimit = 5000
  )

  fits <- lapply(imputed$impDatasets, function(completed) {
    basis <- build_rcs_basis_df(completed[[exposure]], idx = exposure, knots = knots)
    completed <- dplyr::bind_cols(completed, basis$basis_df)
    basis_names <- basis$basis_names
    list(
      base = stats::lm(
        stats::reformulate(covariates, response = outcome_var),
        data = completed
      ),
      linear = stats::lm(
        stats::reformulate(c(basis_names[1], covariates), response = outcome_var),
        data = completed
      ),
      spline = stats::lm(
        stats::reformulate(c(basis_names, covariates), response = outcome_var),
        data = completed
      )
    )
  })

  imputed_values <- unlist(lapply(imputed$impDatasets, function(z) z[[exposure]]))
  stopifnot(all(is.finite(imputed_values)), all(imputed_values > 0))

  mira_of <- function(which) mice::as.mira(lapply(fits, `[[`, which))
  overall <- mice::D1(mira_of("spline"), mira_of("base"))
  nonlinear <- mice::D1(mira_of("spline"), mira_of("linear"))

  # The reduced model must be a strict parameter subset, or D1 silently tests
  # the wrong contrast; nk knots give nk - 1 exposure terms.
  stopifnot(
    overall$result[1, "df1"] == nk - 1,
    nonlinear$result[1, "df1"] == nk - 2
  )

  out <- tibble::tibble(
    Outcome = outcome_label,
    Parameter = exposure_label,
    `Patients in imputation frame` = nrow(frame),
    `Exposure measured` = nrow(observed),
    `Exposure imputed` = nrow(frame) - nrow(observed),
    Imputations = m,
    `Overall P (imputed)` = overall$result[1, "P(>F)"],
    `Nonlinear P (imputed)` = nonlinear$result[1, "P(>F)"]
  )
  if (!is.null(available_case)) {
    # The imputation's observed subset must be exactly the primary model's
    # sample, or the two columns are not comparable.
    stopifnot(nrow(observed) == available_case$n)
    out <- out %>%
      dplyr::mutate(
        `N (available case)` = available_case$n,
        `Overall P (available case)` = available_case$overall_p,
        `Nonlinear P (available case)` = available_case$nonlinear_p,
        `Overall conclusion unchanged` =
          (`Overall P (imputed)` < 0.05) == (available_case$overall_p < 0.05),
        `Nonlinear conclusion unchanged` =
          (`Nonlinear P (imputed)` < 0.05) == (available_case$nonlinear_p < 0.05)
      )
  }
  out
}

# ── Imputation for the trajectory models ─────────────────────────────────────
# smcfcs supports lm, logistic, brlogistic, poisson, weibull, coxph and compet
# substantive models -- there is no mixed-model type -- so the trajectory
# family cannot have a substantive-model-compatible imputation in the sense
# used for the spline models. This is the documented compromise, and the
# limitation is real: the imputation model is not derived from the linear
# mixed model it feeds.
#
# What is done instead. The missing quantity is a patient-level baseline
# measurement, so imputation is at patient level (one row per patient), and the
# substantive model's outcome enters through the two subject-specific summaries
# that the interaction is actually about: each patient's first observed outcome
# and the slope of their own outcome over time. Omitting the slope would impute
# the exposure independently of the trajectory and drive the interaction toward
# zero by construction, which is the failure mode this analysis is meant to
# test rather than create. Imputation is on the log scale, as for the spline
# models, so no negative measurement is ever generated.
#
# `fitter` is passed in rather than looked up, so this function has no hidden
# dependency on a model-fitting helper defined in a stage script.
mi_lmm_exposure_sensitivity <- function(long_analysis_df,
                                        outcome_var,
                                        outcome_label,
                                        exposure,
                                        exposure_label,
                                        baseline_var,
                                        auxiliary_baseline_vars,
                                        fitter,
                                        spec_template,
                                        random_structure = "intercept",
                                        covariates = c("age", "Sex", "BMI"),
                                        m = 20,
                                        seed = 20260916,
                                        available_case = NULL) {
  log_var <- paste0("log_", baseline_var)

  # Same frame rule as the primary confirmation models: usable outcome, time
  # and covariates, and at least two distinct visit times per patient.
  frame <- long_analysis_df %>%
    dplyr::filter(
      !is.na(.data[[outcome_var]]),
      !is.na(time_yrs),
      is.finite(time_yrs),
      dplyr::if_all(dplyr::all_of(covariates), ~ !is.na(.x))
    ) %>%
    dplyr::group_by(ID) %>%
    dplyr::filter(dplyr::n_distinct(time_yrs) >= 2) %>%
    dplyr::ungroup()
  if (!nrow(frame)) return(tibble::tibble())

  patient_level <- frame %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(c(baseline_var, auxiliary_baseline_vars, covariates)), dplyr::first),
      .first_outcome = .data[[outcome_var]][which.min(time_yrs)],
      .own_slope = {
        y <- .data[[outcome_var]]
        t <- time_yrs
        if (stats::var(t) > 0) unname(stats::coef(stats::lm(y ~ t))[2]) else NA_real_
      },
      .n_visits = dplyr::n_distinct(time_yrs),
      .follow_up_span = max(time_yrs) - min(time_yrs),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.na(.own_slope))

  observed <- patient_level %>% dplyr::filter(!is.na(.data[[baseline_var]]))
  if (nrow(observed) < 20) return(tibble::tibble())
  if (any(observed[[baseline_var]] <= 0)) {
    stop(
      sprintf("mi_lmm_exposure_sensitivity(): '%s' has non-positive values.", baseline_var),
      call. = FALSE
    )
  }

  aux_build <- add_missing_indicator_covariates(
    patient_level,
    gated_numeric_vars = auxiliary_baseline_vars
  )
  aux_vars <- aux_build$covariates

  mice_input <- aux_build$data %>%
    dplyr::mutate(
      !!log_var := ifelse(
        is.na(.data[[baseline_var]]),
        NA_real_,
        log(.data[[baseline_var]])
      )
    ) %>%
    dplyr::select(dplyr::all_of(c(
      log_var, covariates, aux_vars,
      ".first_outcome", ".own_slope", ".n_visits", ".follow_up_span"
    ))) %>%
    as.data.frame()

  non_imputed <- setdiff(names(mice_input), log_var)
  incomplete <- non_imputed[
    vapply(non_imputed, function(v) any(is.na(mice_input[[v]])), logical(1))
  ]
  if (length(incomplete)) {
    stop(
      "mi_lmm_exposure_sensitivity(): only the baseline exposure may be imputed, but ",
      paste(incomplete, collapse = ", "), " also contain missing values.",
      call. = FALSE
    )
  }

  method <- rep("", ncol(mice_input))
  names(method) <- names(mice_input)
  method[log_var] <- "norm"

  imputed <- mice::mice(
    mice_input,
    m = m,
    method = method,
    seed = seed,
    printFlag = FALSE
  )

  results <- lapply(seq_len(m), function(i) {
    completed <- mice::complete(imputed, i)
    filled <- tibble::tibble(
      ID = patient_level$ID,
      !!baseline_var := exp(completed[[log_var]])
    )
    stopifnot(all(is.finite(filled[[baseline_var]])), all(filled[[baseline_var]] > 0))
    long_filled <- frame %>%
      dplyr::select(-dplyr::all_of(baseline_var)) %>%
      dplyr::inner_join(filled, by = "ID")
    fit <- fitter(
      long_filled,
      outcome_var,
      spec_template,
      random_structure = random_structure
    )
    if (is.null(fit)) return(NULL)
    list(estimate = fit$estimate, variance = fit$std_error^2, n = dplyr::n_distinct(fit$data$ID))
  })
  results <- results[!vapply(results, is.null, logical(1))]
  if (length(results) < 2) return(tibble::tibble())

  pooled <- rubin_pool(
    vapply(results, function(x) x$estimate, numeric(1)),
    vapply(results, function(x) x$variance, numeric(1))
  )

  out <- pooled %>%
    dplyr::transmute(
      Outcome = outcome_label,
      Parameter = exposure_label,
      `Patients in imputation frame` = nrow(patient_level),
      `Exposure measured` = nrow(observed),
      `Exposure imputed` = nrow(patient_level) - nrow(observed),
      Imputations = length(results),
      `Interaction beta (imputed)` = estimate,
      `CI low` = conf_low,
      `CI high` = conf_high,
      `P value (imputed)` = p_value,
      `Fraction missing information` = fraction_missing_information
    )
  if (!is.null(available_case)) {
    stopifnot(nrow(observed) == available_case$patients)
    out <- out %>%
      dplyr::mutate(
        `Patients (available case)` = available_case$patients,
        `Interaction beta (available case)` = available_case$estimate,
        `P value (available case)` = available_case$p_value,
        `Conclusion unchanged` =
          (`P value (imputed)` < 0.05) == (available_case$p_value < 0.05)
      )
  }
  out
}
