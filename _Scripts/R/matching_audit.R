# Comparator matching and quality-control helpers.
#
# Matching is based only on prespecified baseline variables (age, sex, BMI).
# Clinical outcomes and diastolic/exercise measurements are never used to
# select a matching specification. Only aggregate diagnostics are exported.

extract_matchit_balance <- function(match_object) {
  cobalt::bal.tab(
    match_object,
    un = TRUE,
    m.threshold = 0.10,
    binary = "std"
  )$Balance |>
    tibble::as_tibble(rownames = "variable") |>
    dplyr::filter(variable %in% c("age", "sex_num", "BMI")) |>
    dplyr::transmute(
      variable,
      SMD_unmatched = .data[["Diff.Un"]],
      SMD_matched = .data[["Diff.Adj"]],
      Absolute_SMD_matched = abs(.data[["Diff.Adj"]])
    )
}

extract_matched_arms <- function(match_object, case_df, control_df) {
  matched_pool <- MatchIt::match.data(match_object)
  matched_cases <- matched_pool |>
    dplyr::filter(.is_case == 1L) |>
    dplyr::arrange(as.character(subclass))
  matched_controls <- matched_pool |>
    dplyr::filter(.is_case == 0L) |>
    dplyr::arrange(as.character(subclass))

  stopifnot(nrow(matched_cases) == nrow(matched_controls))

  case_rows <- matched_cases$.orig_row
  control_rows <- matched_controls$.orig_row
  list(
    cases = case_df[case_rows, , drop = FALSE] |>
      dplyr::mutate(match_id = dplyr::row_number()),
    controls = control_df[control_rows, , drop = FALSE] |>
      dplyr::mutate(match_id = dplyr::row_number()),
    pairs = tibble::tibble(
      match_id = seq_len(nrow(matched_cases)),
      case_row = case_rows,
      control_row = control_rows
    )
  )
}

build_matching_candidate_grid <- function() {
  calipers <- c(0.025, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25)
  orders <- c("closest", "farthest", "largest", "smallest")

  dplyr::bind_rows(
    tidyr::crossing(
      Method = "Propensity-score nearest neighbor",
      Caliper = calipers,
      Order = orders
    ),
    tidyr::crossing(
      Method = "Mahalanobis within propensity-score caliper",
      Caliper = calipers,
      Order = orders
    ),
    tibble::tibble(
      Method = "Unrestricted Mahalanobis benchmark",
      Caliper = NA_real_,
      Order = "default"
    )
  ) |>
    dplyr::mutate(
      Candidate_ID = sprintf("M%02d", dplyr::row_number()),
      .before = 1
    )
}

run_matching_candidate <- function(pool, candidate) {
  match_args <- list(
    formula = .is_case ~ age + sex_num + BMI,
    data = pool,
    method = "nearest",
    exact = ~ sex_num,
    replace = FALSE,
    ratio = 1,
    estimand = "ATT"
  )

  if (candidate$Method == "Unrestricted Mahalanobis benchmark") {
    match_args$distance <- "mahalanobis"
  } else {
    match_args$distance <- "glm"
    match_args$caliper <- candidate$Caliper
    match_args$std.caliper <- TRUE
    match_args$m.order <- candidate$Order
    if (
      candidate$Method ==
        "Mahalanobis within propensity-score caliper"
    ) {
      match_args$mahvars <- ~ age + BMI
    }
  }

  tryCatch(
    suppressWarnings(do.call(MatchIt::matchit, match_args)),
    error = function(e) NULL
  )
}

match_controls_balance_grid <- function(
  case_df,
  control_df,
  target_smd = 0.10
) {
  case_df <- case_df |>
    dplyr::mutate(.case_row = dplyr::row_number())
  control_df <- control_df |>
    dplyr::mutate(.control_row = dplyr::row_number())

  prep <- function(data) {
    data |>
      dplyr::mutate(
        age = as.numeric(age),
        sex_num = dplyr::case_when(
          as.character(Sex) == "Female" ~ 1L,
          as.character(Sex) == "Male" ~ 0L,
          TRUE ~ NA_integer_
        ),
        BMI = as.numeric(BMI)
      )
  }
  case_prepped <- prep(case_df) |>
    dplyr::filter(stats::complete.cases(age, sex_num, BMI))
  control_prepped <- prep(control_df) |>
    dplyr::filter(stats::complete.cases(age, sex_num, BMI))

  pool <- dplyr::bind_rows(
    case_prepped |>
      dplyr::transmute(
        .is_case = 1L,
        age,
        sex_num,
        BMI,
        .orig_row = .case_row
      ),
    control_prepped |>
      dplyr::transmute(
        .is_case = 0L,
        age,
        sex_num,
        BMI,
        .orig_row = .control_row
      )
  )

  candidate_grid <- build_matching_candidate_grid()
  candidate_objects <- vector("list", nrow(candidate_grid))
  candidate_summaries <- vector("list", nrow(candidate_grid))

  for (i in seq_len(nrow(candidate_grid))) {
    candidate <- candidate_grid[i, , drop = FALSE]
    match_object <- run_matching_candidate(pool, candidate)
    candidate_objects[[i]] <- match_object

    if (is.null(match_object)) {
      candidate_summaries[[i]] <- candidate |>
        dplyr::mutate(
          Matched_pairs = 0L,
          SMD_age = NA_real_,
          SMD_sex = NA_real_,
          SMD_BMI = NA_real_,
          Maximum_absolute_SMD = NA_real_,
          Balance_pass = FALSE
        )
      next
    }

    balance <- extract_matchit_balance(match_object)
    matched_n <- sum(MatchIt::match.data(match_object)$.is_case == 1L)
    get_smd <- function(variable_name) {
      value <- balance |>
        dplyr::filter(variable == variable_name) |>
        dplyr::pull(SMD_matched)
      if (length(value) == 1) value else NA_real_
    }
    max_smd <- max(balance$Absolute_SMD_matched, na.rm = TRUE)
    candidate_summaries[[i]] <- candidate |>
      dplyr::mutate(
        Matched_pairs = matched_n,
        SMD_age = get_smd("age"),
        SMD_sex = get_smd("sex_num"),
        SMD_BMI = get_smd("BMI"),
        Maximum_absolute_SMD = max_smd,
        Balance_pass = is.finite(max_smd) & max_smd < target_smd
      )
  }

  candidate_summary <- dplyr::bind_rows(candidate_summaries) |>
    dplyr::mutate(
      Method_preference = dplyr::case_when(
        Method == "Mahalanobis within propensity-score caliper" ~ 1L,
        Method == "Propensity-score nearest neighbor" ~ 2L,
        TRUE ~ 3L
      )
    )

  passing <- candidate_summary |>
    dplyr::filter(Balance_pass) |>
    dplyr::arrange(
      dplyr::desc(Matched_pairs),
      Maximum_absolute_SMD,
      Method_preference,
      Caliper,
      Candidate_ID
    )
  if (!nrow(passing)) {
    stop(
      "No prespecified matching candidate achieved all absolute SMDs < ",
      target_smd,
      call. = FALSE
    )
  }

  selected_id <- passing$Candidate_ID[[1]]
  selected_index <- match(selected_id, candidate_grid$Candidate_ID)
  selected_object <- candidate_objects[[selected_index]]
  matched_arms <- extract_matched_arms(
    selected_object,
    case_df,
    control_df
  )
  selected_balance <- extract_matchit_balance(selected_object)

  stopifnot(
    nrow(matched_arms$cases) == nrow(matched_arms$controls),
    all(selected_balance$Absolute_SMD_matched < target_smd)
  )

  list(
    cases = matched_arms$cases,
    controls = matched_arms$controls,
    pairs = matched_arms$pairs,
    match_n = nrow(matched_arms$cases),
    smd_pre = selected_balance |>
      dplyr::select(variable, SMD = SMD_unmatched),
    smd_post = selected_balance |>
      dplyr::select(variable, SMD = SMD_matched),
    match_object = selected_object,
    selected_specification = passing[1, , drop = FALSE],
    candidate_summary = candidate_summary |>
      dplyr::arrange(
        dplyr::desc(Balance_pass),
        dplyr::desc(Matched_pairs),
        Maximum_absolute_SMD
      )
  )
}
