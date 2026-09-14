# ── RCS model infrastructure ──────────────────────────────────────────────────

make_rhs_formula <- function(terms) paste(terms, collapse = " + ")

build_rcs_newdata <- function(var, df, x_vals, adjust_vars) {
  newd <- tibble(!!var := x_vals)
  for (adj_var in adjust_vars) {
    if (!adj_var %in% names(df) || adj_var == var) next
    x <- df[[adj_var]]
    if (is.factor(x)) {
      newd[[adj_var]] <- factor(mode_nonmissing(x), levels = levels(x))
    } else if (is.character(x)) {
      newd[[adj_var]] <- mode_nonmissing(x)
    } else if (is.logical(x)) {
      mv <- mode_nonmissing(x)
      newd[[adj_var]] <- ifelse(is.na(mv), NA, as.logical(mv))
    } else if (is.numeric(x) || is.integer(x)) {
      xn <- x[!is.na(x)]
      newd[[adj_var]] <- if (!length(xn)) NA_real_
        else if (all(unique(xn) %in% c(0, 1))) {
          mv <- mode_nonmissing(xn)
          ifelse(is.na(mv), NA_real_, as.numeric(mv))
        } else median(xn, na.rm = TRUE)
    } else {
      newd[[adj_var]] <- mode_nonmissing(x)
    }
  }
  newd
}

build_rcs_basis_df <- function(x_vals, idx, nk = NULL, knots = NULL) {
  basis <- if (!is.null(knots)) Hmisc::rcspline.eval(x_vals, knots = knots, inclx = TRUE)
           else                  Hmisc::rcspline.eval(x_vals, nk = nk,    inclx = TRUE)
  basis_df <- as_tibble(basis)
  names(basis_df) <- paste0(idx, "_rcs", seq_len(ncol(basis_df)))
  list(basis_df = basis_df, knots = attr(basis, "knots"), basis_names = names(basis_df))
}

build_rcs_prediction_frame <- function(var, df, knots, x_vals, adjust_vars) {
  newd <- build_rcs_newdata(var, df, x_vals, adjust_vars = adjust_vars)
  basis_info <- build_rcs_basis_df(newd[[var]], idx = var, knots = knots)
  bind_cols(newd, basis_info$basis_df)
}

assert_prediction_frame_vars <- function(pred_frame, required_vars, context) {
  missing_vars <- setdiff(required_vars, names(pred_frame))
  if (length(missing_vars))
    stop(sprintf("%s prediction frame missing required covariates: %s",
                 context, paste(missing_vars, collapse = ", ")), call. = FALSE)
  invisible(pred_frame)
}

# Bootstrap-based inference for the data-dredged piecewise breakpoint test.
# The original implementation reported anova(linear, best_pw)$Pr(>F) without
# any correction for the search over `grid_n` candidate knots. The reported
# p-value is therefore anti-conservative (the null distribution of the test
# statistic when the best knot is selected by AIC is wider than F_{1, n-p}).
# This function constructs an empirical null distribution of the best-AIC
# piecewise F-statistic under H0: outcome is linear in `idx` (no breakpoint),
# and returns a bootstrap p-value.
#
# Caching: results are saved by data fingerprint to RCS_BOOT_CACHE_DIR so that
# rerenders do not repeat the simulation. Pass force = TRUE to invalidate.
RCS_BOOT_CACHE_DIR  <- "../2_Output/cache/rcs_piecewise_boot"
RCS_BOOTSTRAP_REPS  <- 2000
RCS_BOOTSTRAP_SEED  <- 20260514

bootstrap_piecewise_p <- function(df_idx, outcome_var, idx,
                                  linear_model, adjust_vars,
                                  observed_F, observed_knot,
                                  B = RCS_BOOTSTRAP_REPS,
                                  min_side_n = 20, grid_n = 41,
                                  cache_dir = RCS_BOOT_CACHE_DIR,
                                  force = FALSE) {
  if (!is.finite(observed_F)) return(list(p_boot = NA_real_, B = 0L, F_null = numeric(0)))

  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  cache_key <- digest::digest(list(
    outcome_var = outcome_var, idx = idx,
    n = nrow(df_idx),
    adjust_vars = sort(adjust_vars),
    B = B, grid_n = grid_n, min_side_n = min_side_n,
    seed = RCS_BOOTSTRAP_SEED,
    fitted = sum(fitted(linear_model)),
    rss    = sum(residuals(linear_model)^2)
  ))
  cache_path <- file.path(cache_dir, paste0("piecewise_", cache_key, ".rds"))

  if (!force && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    return(list(
      p_boot = mean(cached$F_null >= observed_F, na.rm = TRUE),
      B      = sum(!is.na(cached$F_null)),
      F_null = cached$F_null
    ))
  }

  set.seed(RCS_BOOTSTRAP_SEED)
  mu_hat   <- fitted(linear_model)
  sigma_hat <- sqrt(sum(residuals(linear_model)^2) /
                    linear_model$df.residual)
  x        <- df_idx[[idx]]
  candidate_knots <- unique(as.numeric(stats::quantile(
    x, probs = seq(0.15, 0.85, length.out = grid_n), na.rm = TRUE, type = 8)))
  candidate_knots <- candidate_knots[is.finite(candidate_knots)]
  candidate_knots <- candidate_knots[vapply(candidate_knots, function(k) {
    sum(x <= k, na.rm = TRUE) >= min_side_n && sum(x > k, na.rm = TRUE) >= min_side_n
  }, logical(1))]
  if (!length(candidate_knots)) return(list(p_boot = NA_real_, B = 0L, F_null = numeric(0)))

  null_formula_linear <- as.formula(paste0(outcome_var, " ~ ",
                                            make_rhs_formula(c(adjust_vars, idx))))

  F_null <- vapply(seq_len(B), function(b) {
    y_star <- mu_hat + rnorm(length(mu_hat), 0, sigma_hat)
    df_b   <- df_idx
    df_b[[outcome_var]] <- y_star
    lin_b  <- lm(null_formula_linear, data = df_b)
    pw_aic <- vapply(candidate_knots, function(k) {
      df_b$.hinge <- pmax(df_b[[idx]] - k, 0)
      fit_b <- lm(as.formula(paste0(outcome_var, " ~ ",
                                     make_rhs_formula(c(adjust_vars, idx, ".hinge")))),
                  data = df_b)
      AIC(fit_b)
    }, numeric(1))
    best_k <- candidate_knots[which.min(pw_aic)]
    df_b$.hinge <- pmax(df_b[[idx]] - best_k, 0)
    best_b <- lm(as.formula(paste0(outcome_var, " ~ ",
                                    make_rhs_formula(c(adjust_vars, idx, ".hinge")))),
                 data = df_b)
    aov_b  <- anova(lin_b, best_b)
    f_b    <- aov_b$F[2]
    if (is.null(f_b) || !is.finite(f_b)) NA_real_ else f_b
  }, numeric(1))

  saveRDS(list(F_null = F_null, observed_F = observed_F, observed_knot = observed_knot,
               cache_key = cache_key, timestamp = Sys.time()),
          cache_path)

  list(p_boot = mean(F_null >= observed_F, na.rm = TRUE),
       B      = sum(!is.na(F_null)),
       F_null = F_null)
}

fit_piecewise_model <- function(df_idx, outcome_var, idx, linear_model,
                                adjust_vars, min_side_n = 20, grid_n = 41) {
  x <- df_idx[[idx]]
  candidate_knots <- unique(as.numeric(stats::quantile(
    x, probs = seq(0.15, 0.85, length.out = grid_n), na.rm = TRUE, type = 8)))
  candidate_knots <- candidate_knots[is.finite(candidate_knots)]
  candidate_knots <- candidate_knots[vapply(candidate_knots, function(k) {
    sum(x <= k, na.rm = TRUE) >= min_side_n && sum(x > k, na.rm = TRUE) >= min_side_n
  }, logical(1))]
  if (!length(candidate_knots)) return(NULL)

  piecewise_candidates <- map_dfr(candidate_knots, function(knot) {
    df_pw  <- df_idx %>% mutate(.hinge = pmax(.data[[idx]] - knot, 0))
    fit_obj <- lm(as.formula(paste0(outcome_var, " ~ ",
                                    make_rhs_formula(c(adjust_vars, idx, ".hinge")))),
                  data = df_pw)
    tibble(knot = knot, fit = list(fit_obj), AIC = AIC(fit_obj), BIC = BIC(fit_obj))
  }) %>% arrange(AIC)

  best_row  <- piecewise_candidates %>% slice(1)
  best_fit  <- best_row$fit[[1]]
  pw_test   <- anova(linear_model, best_fit)
  observed_F <- pw_test$F[2]
  vc        <- vcov(best_fit)
  beta_linear <- unname(coef(best_fit)[idx])
  beta_hinge  <- unname(coef(best_fit)[".hinge"])
  slope_pre_se  <- sqrt(vc[idx, idx])
  slope_post_se <- sqrt(vc[idx, idx] + vc[".hinge", ".hinge"] + 2 * vc[idx, ".hinge"])
  support_tbl   <- piecewise_candidates %>%
    mutate(delta_AIC = AIC - min(AIC, na.rm = TRUE))

  # Bootstrap p-value to correct for AIC-driven knot selection
  boot_res <- bootstrap_piecewise_p(df_idx, outcome_var, idx, linear_model,
                                     adjust_vars, observed_F = observed_F,
                                     observed_knot = best_row$knot[[1]],
                                     min_side_n = min_side_n, grid_n = grid_n)

  list(
    fit             = best_fit,
    best_knot       = best_row$knot[[1]],
    delta_AIC       = AIC(best_fit) - AIC(linear_model),
    test            = pw_test,
    bootstrap       = boot_res,
    candidate_table = support_tbl %>% select(knot, AIC, BIC, delta_AIC),
    knot_support_low  = min(support_tbl$knot[support_tbl$delta_AIC <= 2], na.rm = TRUE),
    knot_support_high = max(support_tbl$knot[support_tbl$delta_AIC <= 2], na.rm = TRUE),
    slope_pre         = beta_linear,
    slope_post        = beta_linear + beta_hinge,
    slope_pre_lo      = beta_linear - 1.96 * slope_pre_se,
    slope_pre_hi      = beta_linear + 1.96 * slope_pre_se,
    slope_post_lo     = beta_linear + beta_hinge - 1.96 * slope_post_se,
    slope_post_hi     = beta_linear + beta_hinge + 1.96 * slope_post_se
  )
}

estimate_rcs_curvature <- function(df, model, idx, knots, lo, hi,
                                   adjust_vars, n = 121) {
  x_grid    <- seq(lo, hi, length.out = n)
  base_step <- diff(range(df[[idx]], na.rm = TRUE)) / 400
  base_step <- ifelse(is.finite(base_step) && base_step > 0, base_step, 0.01)
  beta <- coef(model)
  vc   <- vcov(model)

  curvature_tbl <- map_dfr(x_grid, function(x0) {
    step <- min(base_step, x0 - lo, hi - x0)
    if (!is.finite(step) || step <= 0) return(NULL)
    x_triplet  <- c(x0 - step, x0, x0 + step)
    pred_frame <- build_rcs_prediction_frame(idx, df, knots, x_triplet, adjust_vars)
    assert_prediction_frame_vars(pred_frame, adjust_vars, "RCS curvature")
    mm <- model.matrix(delete.response(terms(model)), data = pred_frame)
    d1 <- (mm[3, ] - mm[1, ])                        / (2 * step)
    d2 <- (mm[3, ] - 2 * mm[2, ] + mm[1, ])          / step^2
    slope      <- sum(d1 * beta)
    slope_se   <- sqrt(as.numeric(d1 %*% vc %*% d1))
    curvature  <- sum(d2 * beta)
    curv_se    <- sqrt(as.numeric(d2 %*% vc %*% d2))
    curv_z     <- ifelse(curv_se > 0, curvature / curv_se, NA_real_)
    tibble(
      x           = x0,
      slope       = slope,
      slope_lo    = slope - 1.96 * slope_se,
      slope_hi    = slope + 1.96 * slope_se,
      curvature   = curvature,
      curvature_lo = curvature - 1.96 * curv_se,
      curvature_hi = curvature + 1.96 * curv_se,
      curvature_p = ifelse(is.na(curv_z), NA_real_, 2 * stats::pnorm(-abs(curv_z)))
    )
  })

  if (!nrow(curvature_tbl)) return(tibble())
  curvature_tbl %>%
    mutate(curvature_q = p.adjust(curvature_p, method = "BH"),
           significant_curvature = !is.na(curvature_q) & curvature_q < 0.05)
}

summarize_curvature_regions <- function(curvature_tbl, digits = 1) {
  sig_tbl <- curvature_tbl %>% filter(significant_curvature) %>% arrange(x)
  if (!nrow(sig_tbl)) return("None")
  median_gap <- stats::median(diff(sort(unique(curvature_tbl$x))), na.rm = TRUE)
  if (!is.finite(median_gap) || median_gap <= 0) median_gap <- 0
  sig_tbl <- sig_tbl %>%
    mutate(
      curvature_direction = ifelse(curvature > 0, "Upward", "Downward"),
      region_id = cumsum(c(TRUE,
        (diff(x) > 1.5 * median_gap) |
        (curvature_direction[-1] != curvature_direction[-n()])))
    )
  region_tbl <- sig_tbl %>%
    group_by(region_id, curvature_direction) %>%
    summarise(start = min(x, na.rm = TRUE), end = max(x, na.rm = TRUE),
              n_points = dplyr::n(), .groups = "drop") %>%
    filter(n_points >= 2 | (end - start) > median_gap)
  if (!nrow(region_tbl)) return("None")
  paste(sprintf(paste0("%s %.", digits, "f–%.", digits, "f"),
                region_tbl$curvature_direction, region_tbl$start, region_tbl$end),
        collapse = "; ")
}

fit_rcs_model <- function(df, outcome_var, idx, nk = 4,
                          adjust_vars   = figure3_adjustment_vars,
                          adjust_numeric = figure3_adjustment_numeric,
                          run_piecewise = TRUE) {
  adjust_vars    <- intersect(adjust_vars,    names(df))
  adjust_numeric <- intersect(adjust_numeric, names(df))
  df_idx <- df %>%
    mutate(across(any_of(c(outcome_var, idx, adjust_numeric)), as_num),
           Sex = factor(Sex, levels = c("Male", "Female"))) %>%
    filter(if_all(all_of(c(outcome_var, idx, adjust_vars)), ~ !is.na(.x)))
  if (nrow(df_idx) < 20) return(NULL)

  base   <- lm(as.formula(paste0(outcome_var, " ~ ", make_rhs_formula(adjust_vars))), data = df_idx)
  linear <- lm(as.formula(paste0(outcome_var, " ~ ", make_rhs_formula(c(adjust_vars, idx)))), data = df_idx)
  basis_info  <- build_rcs_basis_df(df_idx[[idx]], idx = idx, nk = nk)
  df_model    <- bind_cols(df_idx, basis_info$basis_df)
  spline      <- lm(as.formula(paste(outcome_var, "~", make_rhs_formula(c(adjust_vars, basis_info$basis_names)))),
                    data = df_model)
  # Overall association for an RCS exposure must compare the covariate-only
  # model with the full spline model. The base-vs-linear comparison is retained
  # separately as a descriptive linear-component test.
  linear_test    <- anova(base, linear)
  overall_test   <- anova(base, spline)
  nonlinear_test <- anova(linear, spline)
  piecewise <- if (isTRUE(run_piecewise)) {
    fit_piecewise_model(
      df_idx,
      outcome_var,
      idx,
      linear_model = linear,
      adjust_vars = adjust_vars
    )
  } else {
    NULL
  }
  lo <- quantile(df_idx[[idx]], 0.02, na.rm = TRUE)
  hi <- quantile(df_idx[[idx]], 0.98, na.rm = TRUE)
  curvature_tbl  <- estimate_rcs_curvature(df_idx, spline, idx, basis_info$knots,
                                           lo, hi, adjust_vars = adjust_vars, n = 121)
  list(
    df = df_idx, base = base, linear = linear, spline = spline,
    linear_test = linear_test, overall_test = overall_test,
    nonlinear_test = nonlinear_test,
    knots = basis_info$knots, adjust_vars = adjust_vars,
    piecewise = piecewise, curvature = curvature_tbl,
    curvature_zone = summarize_curvature_regions(curvature_tbl)
  )
}

pred_grid_rcs <- function(var, df, model, knots, lo, hi, adjust_vars, n = 200) {
  pred_df <- build_rcs_prediction_frame(var, df, knots, seq(lo, hi, length.out = n),
                                        adjust_vars = adjust_vars)
  assert_prediction_frame_vars(pred_df, adjust_vars, "RCS spline")
  p <- predict(model, newdata = pred_df, se.fit = TRUE)
  tibble(x = pred_df[[var]], fit = p$fit, se = p$se.fit,
         lo = p$fit - 1.96 * p$se.fit, hi = p$fit + 1.96 * p$se.fit, var = var)
}

pred_grid_piecewise <- function(var, df, model, knot, lo, hi, adjust_vars, n = 200) {
  pred_df <- build_rcs_newdata(var, df, seq(lo, hi, length.out = n), adjust_vars) %>%
    mutate(.hinge = pmax(.data[[var]] - knot, 0))
  assert_prediction_frame_vars(pred_df, adjust_vars, "RCS piecewise")
  p <- predict(model, newdata = pred_df, se.fit = TRUE)
  tibble(x = pred_df[[var]], fit = p$fit, var = var)
}

extract_rcs_stats <- function(fit_obj, idx, outcome_key) {
  tibble(
    outcome = outcome_key, index = idx, label = hcm_primary_labels[idx],
    n = nrow(fit_obj$df),
    delta_AIC_linear    = AIC(fit_obj$linear) - AIC(fit_obj$base),
    delta_AIC_nonlinear = AIC(fit_obj$spline) - AIC(fit_obj$linear),
    linear_F   = fit_obj$linear_test$F[2],
    linear_p   = fit_obj$linear_test$`Pr(>F)`[2],
    overall_F  = fit_obj$overall_test$F[2],
    overall_p  = fit_obj$overall_test$`Pr(>F)`[2],
    nonlinear_F = fit_obj$nonlinear_test$F[2],
    nonlinear_p = fit_obj$nonlinear_test$`Pr(>F)`[2],
    piecewise_knot      = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$best_knot),
    piecewise_delta_AIC = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$delta_AIC),
    piecewise_F = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$test$F[2]),
    piecewise_p = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$test$`Pr(>F)`[2]),
    piecewise_p_boot = ifelse(is.null(fit_obj$piecewise) ||
                               is.null(fit_obj$piecewise$bootstrap),
                              NA_real_, fit_obj$piecewise$bootstrap$p_boot),
    piecewise_boot_B = ifelse(is.null(fit_obj$piecewise) ||
                                is.null(fit_obj$piecewise$bootstrap),
                              NA_integer_, fit_obj$piecewise$bootstrap$B),
    slope_pre_knot    = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$slope_pre),
    slope_pre_knot_lo = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$slope_pre_lo),
    slope_pre_knot_hi = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$slope_pre_hi),
    slope_post_knot    = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$slope_post),
    slope_post_knot_lo = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$slope_post_lo),
    slope_post_knot_hi = ifelse(is.null(fit_obj$piecewise), NA_real_, fit_obj$piecewise$slope_post_hi),
    knot_support_range = ifelse(
      is.null(fit_obj$piecewise), NA_character_,
      sprintf("%.1f-%.1f", fit_obj$piecewise$knot_support_low, fit_obj$piecewise$knot_support_high)
    ),
    curvature_zone     = fit_obj$curvature_zone,
    any_local_curvature = ifelse(nrow(fit_obj$curvature) > 0,
                                 any(fit_obj$curvature$significant_curvature, na.rm = TRUE), FALSE),
    R2_base   = summary(fit_obj$base)$r.squared,
    R2_linear = summary(fit_obj$linear)$r.squared,
    R2_spline = summary(fit_obj$spline)$r.squared
  )
}

annotate_rcs_lrt_table <- function(tbl) {
  if (is.null(tbl) || !nrow(tbl) ||
      !all(c("nonlinear_p", "piecewise_p") %in% names(tbl))) return(tibble())
  # Use the bootstrap-corrected piecewise p when available; fall back to the
  # uncorrected anova p only if bootstrap was not computed (e.g., very small
  # samples where no candidate knot survives min_side_n filter).
  tbl_out <- tbl %>% arrange(nonlinear_p)
  if ("piecewise_p_boot" %in% names(tbl_out)) {
    tbl_out <- tbl_out %>%
      mutate(piecewise_p_used = ifelse(is.na(piecewise_p_boot),
                                        piecewise_p, piecewise_p_boot))
  } else {
    tbl_out <- tbl_out %>% mutate(piecewise_p_used = piecewise_p)
  }
  tbl_out %>%
    mutate(
      nonlinear_q = p.adjust(nonlinear_p,       method = "BH"),
      piecewise_q = p.adjust(piecewise_p_used,  method = "BH"),
      meaningful_nonlinearity = nonlinear_q < 0.10 & delta_AIC_nonlinear <= -2,
      meaningful_piecewise    = piecewise_q < 0.10 & piecewise_delta_AIC <= -2
    )
}

format_rcs_stats_label <- function(stats_row, compact = FALSE) {
  fmt_sci <- function(p) {
    if (is.na(p)) return("NA")
    if (p >= 0.001) return(sprintf("%.3f", p))
    sci <- scientific_p_parts(p, digits = 1)
    paste0(sci$mantissa, " × 10<sup>", sci$exponent, "</sup>")
  }
  if (isTRUE(compact)) {
    n_label <- if (identical(as.character(stats_row$index), "tr_max_vel")) {
      paste0("N = ", stats_row$n, " (available case)")
    } else {
      paste0("N = ", stats_row$n)
    }
    return(paste0(
      n_label,
      "<br>Overall P = ", fmt_sci(stats_row$overall_p),
      "<br>Nonlinearity q = ", fmt_sci(stats_row$nonlinear_q)
    ))
  }
  paste0(
    "Patients: ", stats_row$n,
    "<br>Overall P = ", fmt_sci(stats_row$overall_p),
    "<br>NL q = ", fmt_sci(stats_row$nonlinear_q),
    "<br>Delta AIC = ",
    number(stats_row$delta_AIC_nonlinear, accuracy = 0.01)
  )
}

build_single_outcome_rcs_summary_table <- function(tbl, outcome_label) {
  tbl %>%
    mutate(
      Interpretation = case_when(
        meaningful_nonlinearity & meaningful_piecewise ~ "Meaningful nonlinearity with piecewise threshold",
        meaningful_nonlinearity                        ~ "Meaningful nonlinearity",
        meaningful_piecewise                           ~ "Piecewise threshold without strong spline nonlinearity",
        overall_p < 0.05                               ~ "Approximately linear association",
        TRUE                                           ~ "No meaningful association"
      ),
      overall_p_fmt       = ifelse(overall_p   < 0.001, "<0.001", sprintf("%.3f", overall_p)),
      nonlinear_q_fmt     = ifelse(nonlinear_q < 0.001, "<0.001", sprintf("%.3f", nonlinear_q)),
      piecewise_q_fmt     = ifelse(piecewise_q < 0.001, "<0.001", sprintf("%.3f", piecewise_q)),
      delta_AIC_nonlinear = round(delta_AIC_nonlinear, 2),
      piecewise_knot_fmt  = ifelse(is.na(piecewise_knot), "NA", sprintf("%.1f", piecewise_knot)),
      Outcome             = outcome_label,
      Parameter           = label,
      N                   = n,
      `Overall p`         = overall_p_fmt,
      `Nonlinearity q`    = nonlinear_q_fmt,
      `Spline vs linear Delta AIC` = delta_AIC_nonlinear,
      `Piecewise knot`    = piecewise_knot_fmt,
      `Piecewise q`       = piecewise_q_fmt,
      `Curvature zone`    = curvature_zone
    ) %>%
    select(Outcome, Parameter, N, `Overall p`, `Nonlinearity q`,
           `Spline vs linear Delta AIC`, `Piecewise knot`, `Piecewise q`,
           `Curvature zone`, Interpretation)
}

# ── RCS panel theming and plotting ────────────────────────────────────────────
theme_jacc_rcs <- function(base_size = 9.5) {
  theme_jacc(base_size = base_size) +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      legend.position  = "none",
      strip.background = element_rect(fill = "grey97", color = "grey75", linewidth = 0.3),
      strip.text       = element_text(face = "bold", size = base_size,
                                      margin = margin(4, 4, 4, 4)),
      panel.spacing    = grid::unit(1.0, "lines"),
      axis.title.y     = ggtext::element_markdown(margin = margin(r = 6)),
      axis.text        = element_text(size = base_size - 0.3),
      plot.margin      = margin(4, 10, 4, 4)
    )
}

make_rcs_panel <- function(pred_df, idx, y_label, line_color, fill_color, y_limits,
                           show_y = FALSE, stats_label = NULL, point_df = NULL,
                           stats_position = "bottom_right", piecewise_df = NULL,
                           stats_fill_alpha = 0.9) {
  p <- ggplot(pred_df, aes(x = x, y = fit)) +
    {if (!is.null(point_df) && nrow(point_df) > 0)
      geom_point(data = point_df, aes(x = x, y = y), inherit.aes = FALSE,
                 color = alpha("grey28", 0.18), size = 0.6, shape = 16)} +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.16, fill = alpha(fill_color, 0.82)) +
    {if (!is.null(piecewise_df) && nrow(piecewise_df) > 0)
      geom_line(data = piecewise_df, aes(x = x, y = fit), inherit.aes = FALSE,
                color = "black", linewidth = 0.7, linetype = "22")} +
    geom_line(linewidth = 1.05, color = line_color) +
    geom_vline(xintercept = unique(pred_df$ase_cutoff),
               color = jacc_cols["gray"], linetype = "dashed", linewidth = 0.55) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = y_limits, expand = expansion(mult = c(0.03, 0.08)),
                       labels = label_number(accuracy = 1)) +
    labs(title = hcm_primary_labels[idx], x = NULL, y = y_label) +
    theme_jacc_rcs() +
    theme(plot.title = element_text(face = "bold", size = 9.5, hjust = 0.5),
          plot.margin = margin(4, 4, 4, 4))

  p <- add_plot_stats_box(p, stats_label = stats_label,
                          x_range = range(pred_df$x, na.rm = TRUE),
                          y_range = range(y_limits, na.rm = TRUE),
                          position = stats_position, size = 2.5,
                          fill_alpha = stats_fill_alpha)

  if (!show_y) p <- p + theme(axis.title.y = element_blank(),
                               axis.text.y  = element_blank(),
                               axis.ticks.y = element_blank())
  p
}
