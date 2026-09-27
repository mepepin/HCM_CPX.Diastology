# Small shared utilities. Ported from _Scripts/R/utils.R (lines 512-535).

# ── Misc helpers ──────────────────────────────────────────────────────────────
fmt_n   <- function(x) format(x, big.mark = ",")
fmt_pct <- function(x) sprintf("%.1f", x)

mode_nonmissing <- function(x) {
  x_nonmiss <- x[!is.na(x)]
  if (!length(x_nonmiss)) return(NA)
  names(sort(table(x_nonmiss), decreasing = TRUE))[1]
}

residualize_on_age <- function(df, idx) {
  df <- df %>% mutate(.row_id = row_number())
  dsub <- df %>% filter(!is.na(.data[[idx]]), !is.na(age))
  if (nrow(dsub) < 10) {
    df[[paste0(idx, "_resid")]] <- NA_real_
    return(df %>% select(-.row_id))
  }
  fit <- lm(as.formula(paste0(idx, " ~ ns(age, df = 3)")), data = dsub)
  dsub$.resid <- resid(fit)
  df <- df %>% left_join(dsub %>% select(.row_id, .resid), by = ".row_id")
  df[[paste0(idx, "_resid")]] <- df$.resid
  df %>% select(-.row_id, -.resid)
}


# Excel columns often arrive as character/factor; coerce without warning spam.
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# NA-safe reducers; Excel dates arrive either as serial numbers or as Date.
coerce_excel_date <- function(x) {
  if (inherits(x, "Date")) {
    return(as.Date(x))
  }
  x_num <- suppressWarnings(as.numeric(x))
  as.Date(openxlsx::convertToDate(x_num))
}

min_date_or_na <- function(x) {
  x <- as.Date(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) as.Date(NA) else min(x)
}

max_date_or_na <- function(x) {
  x <- as.Date(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) as.Date(NA) else max(x)
}

max_numeric_or_na <- function(x) {
  x <- as_num(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else max(x)
}

min_numeric_or_na <- function(x) {
  x <- as_num(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else min(x)
}
