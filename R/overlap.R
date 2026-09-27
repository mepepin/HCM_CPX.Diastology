overlap_weight_comparison <- function(case_df, control_df,
                                      covariates = c("age", "sex_num", "BMI")) {
  prep <- function(df, is_case) {
    df %>%
      dplyr::mutate(
        .is_case = is_case,
        age = as_num(age),
        sex_num = ifelse(as.character(Sex) == "Female", 1L, 0L),
        BMI = as_num(BMI)
      )
  }
  cases <- prep(case_df, 1L)
  controls <- prep(control_df, 0L)
  complete_rows <- function(df) df[stats::complete.cases(df[, covariates]), , drop = FALSE]
  cases <- complete_rows(cases)
  controls <- complete_rows(controls)

  pooled <- dplyr::bind_rows(
    cases %>% dplyr::select(dplyr::all_of(c(".is_case", covariates))),
    controls %>% dplyr::select(dplyr::all_of(c(".is_case", covariates)))
  )
  ps_model <- stats::glm(
    stats::reformulate(covariates, response = ".is_case"),
    data = pooled, family = stats::binomial()
  )
  e <- stats::fitted(ps_model)
  e_cases <- e[seq_len(nrow(cases))]
  e_controls <- e[nrow(cases) + seq_len(nrow(controls))]

  cases$.weight <- 1 - e_cases      # overlap weight for the case arm
  controls$.weight <- e_controls    # overlap weight for the comparator arm
  # Keep the score itself: the supplemental overlap figure plots its
  # distribution by group before and after weighting (card 2.2).
  cases$.ps <- e_cases
  controls$.ps <- e_controls

  smd <- function(x, w_case, w_control) {
    xc <- x[seq_len(nrow(cases))]
    xk <- x[nrow(cases) + seq_len(nrow(controls))]
    m1 <- stats::weighted.mean(xc, w_case)
    m0 <- stats::weighted.mean(xk, w_control)
    v1 <- sum(w_case * (xc - m1)^2) / sum(w_case)
    v0 <- sum(w_control * (xk - m0)^2) / sum(w_control)
    (m1 - m0) / sqrt((v1 + v0) / 2)
  }
  ones_case <- rep(1, nrow(cases))
  ones_control <- rep(1, nrow(controls))
  balance <- tibble::tibble(
    variable = covariates,
    `SMD (unweighted)` = vapply(covariates, function(v)
      smd(pooled[[v]], ones_case, ones_control), numeric(1)),
    `SMD (overlap-weighted)` = vapply(covariates, function(v)
      smd(pooled[[v]], cases$.weight, controls$.weight), numeric(1))
  )
  ess <- function(w) sum(w)^2 / sum(w^2)

  # Weighted and unweighted covariate moments per arm, which is what a reviewer
  # checks the SMDs against (card 2.2).
  weighted_sd <- function(x, w) {
    m <- stats::weighted.mean(x, w)
    sqrt(sum(w * (x - m)^2) / sum(w))
  }
  covariate_summary <- purrr::map_dfr(covariates, function(v) {
    xc <- cases[[v]]
    xk <- controls[[v]]
    tibble::tibble(
      variable = v,
      `HCM mean (unweighted)` = mean(xc),
      `HCM SD (unweighted)` = stats::sd(xc),
      `Reference mean (unweighted)` = mean(xk),
      `Reference SD (unweighted)` = stats::sd(xk),
      `HCM mean (overlap-weighted)` = stats::weighted.mean(xc, cases$.weight),
      `HCM SD (overlap-weighted)` = weighted_sd(xc, cases$.weight),
      `Reference mean (overlap-weighted)` =
        stats::weighted.mean(xk, controls$.weight),
      `Reference SD (overlap-weighted)` =
        weighted_sd(xk, controls$.weight)
    )
  }) %>%
    dplyr::left_join(balance, by = "variable")

  list(
    cases = cases,
    controls = controls,
    balance = balance,
    covariate_summary = covariate_summary,
    n_cases = nrow(cases),
    n_controls = nrow(controls),
    ess_cases = ess(cases$.weight),
    ess_controls = ess(controls$.weight),
    ps_model = ps_model
  )
}
