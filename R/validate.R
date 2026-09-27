# ── De-identification guard for exported results ──────────────────────────────
# 6_Results/ is tracked in git and may be published, so nothing patient-level
# may reach it. Patient-level tables belong in 5_Data/ (git-ignored).
#
# assert_results_deidentified() scans every CSV under the results tree and
# stops the pipeline if it finds either
#   * a column named like an identifier (MRN, ID, DOB, birth_date, ...), or
#   * a cell that is entirely a 7-10 digit integer (the shape of a medical
#     record number), excluding yyyymmdd dates,
# outside the columns named in `allow` (file sizes, checksums, and similar
# legitimately long values).
#
# Called at the end of the last stage, so a run cannot finish while an
# identifier is sitting in an exported file.

identifier_column_pattern <- "^(mrn|id|patient_?id|subject_?id|studyid|study_id|dob|birth_?date|record_?id|ssn|accession)$"

assert_results_deidentified <- function(results_dir,
                                        allow = c("size_bytes", "sha256", "n_pairs")) {
  csv_files <- list.files(results_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  problems <- character(0)

  for (path in csv_files) {
    tab <- tryCatch(
      utils::read.csv(path, colClasses = "character", check.names = FALSE, nrows = -1),
      error = function(e) NULL
    )
    if (is.null(tab) || !ncol(tab)) next

    nms <- tolower(trimws(names(tab)))
    bad_cols <- names(tab)[grepl(identifier_column_pattern, nms)]
    if (length(bad_cols)) {
      problems <- c(problems, sprintf("%s: identifier column(s) %s",
                                      path, paste(bad_cols, collapse = ", ")))
    }

    for (j in seq_along(tab)) {
      if (tolower(names(tab)[j]) %in% tolower(allow)) next
      # Only whole-cell numbers of record-number length count. Digits embedded
      # in text (sheet names, file names) and yyyymmdd dates are not
      # identifiers, and flagging them made the check cry wolf.
      whole_number <- grepl("^\\s*[0-9]{7,10}\\s*$", tab[[j]])
      date_like <- grepl("^\\s*(19|20)[0-9]{6}\\s*$", tab[[j]])
      hits <- whole_number & !date_like
      if (any(hits, na.rm = TRUE)) {
        problems <- c(problems, sprintf("%s: column '%s' has %d value(s) shaped like a record number",
                                        path, names(tab)[j], sum(hits, na.rm = TRUE)))
      }
    }
  }

  if (length(problems)) {
    stop("De-identification check failed for exported results:\n  ",
         paste(problems, collapse = "\n  "),
         "\nMove patient-level tables to the processed-data directory (5_Data/2_Processed).",
         call. = FALSE)
  }

  message(sprintf("De-identification check passed: %d result CSVs, no identifiers.", length(csv_files)))
  invisible(length(csv_files))
}

# ── Input-workbook formula integrity ──────────────────────────────────────────
# openxlsx returns only the values Excel cached, so a formula that points at
# the wrong row is invisible to the pipeline. These helpers read the formulas
# themselves out of the .xlsx (a zip of XML) and report any formula whose
# references point at a row other than its own within the same sheet.
#
# This is the defect that produced the September 2026 derived-field
# misalignment: the CPX sheet was sorted after its formulas were written, so
# every formula stayed behind pointing at the row it used to occupy. The check
# exists so that a re-sorted or newly exported workbook cannot reach the
# analysis unnoticed.

excel_sheet_formula_summary <- function(path) {
  tmp <- tempfile()
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(path, exdir = tmp)

  book <- xml2::read_xml(file.path(tmp, "xl", "workbook.xml"))
  ns <- xml2::xml_ns(book)
  sheet_nodes <- xml2::xml_find_all(book, "//d1:sheets/d1:sheet", ns)
  sheet_names <- xml2::xml_attr(sheet_nodes, "name")
  sheet_rids <- xml2::xml_attr(sheet_nodes, "id")

  rels <- xml2::read_xml(file.path(tmp, "xl", "_rels", "workbook.xml.rels"))
  rel_nodes <- xml2::xml_find_all(rels, "//*[local-name()='Relationship']")
  rel_map <- stats::setNames(
    xml2::xml_attr(rel_nodes, "Target"),
    xml2::xml_attr(rel_nodes, "Id")
  )

  purrr::map_dfr(seq_along(sheet_names), function(i) {
    target <- rel_map[[sheet_rids[i]]]
    sheet_path <- file.path(tmp, "xl", sub("^/?xl/", "", target))
    if (!file.exists(sheet_path)) {
      return(tibble::tibble(
        workbook = basename(path), sheet = sheet_names[i], n_formulas = NA_integer_,
        n_cross_row = NA_integer_, n_ref_errors = NA_integer_,
        n_shared_followers = NA_integer_, cross_row_columns = "unreadable"
      ))
    }
    sheet <- xml2::read_xml(sheet_path)
    f_nodes <- xml2::xml_find_all(sheet, "//*[local-name()='f']")
    if (!length(f_nodes)) {
      return(tibble::tibble(
        workbook = basename(path), sheet = sheet_names[i], n_formulas = 0L,
        n_cross_row = 0L, n_ref_errors = 0L, n_shared_followers = 0L,
        cross_row_columns = ""
      ))
    }
    formula_text <- xml2::xml_text(f_nodes)
    cell_ref <- xml2::xml_attr(xml2::xml_parent(f_nodes), "r")
    own_row <- as.integer(sub("^[A-Z]+", "", cell_ref))
    own_col <- sub("[0-9]+$", "", cell_ref)
    followers <- formula_text == "" |
      (xml2::xml_attr(f_nodes, "t") %in% "shared" & formula_text == "")

    # Same-sheet references only: a reference qualified with a different sheet
    # name is a legitimate lookup, not a misaligned row.
    stripped <- gsub("'[^']*'!|\\b[A-Za-z_][A-Za-z0-9_.]*!", "", formula_text)
    referenced_rows <- stringr::str_extract_all(stripped, "\\$?[A-Z]{1,3}\\$?[0-9]+")
    cross_row <- vapply(seq_along(referenced_rows), function(k) {
      if (followers[k] || !length(referenced_rows[[k]])) return(FALSE)
      rows <- as.integer(sub("^\\$?[A-Z]{1,3}\\$?", "", referenced_rows[[k]]))
      any(rows != own_row[k], na.rm = TRUE)
    }, logical(1))

    tibble::tibble(
      workbook = basename(path),
      sheet = sheet_names[i],
      n_formulas = length(f_nodes),
      n_cross_row = sum(cross_row),
      n_ref_errors = sum(grepl("#REF!", formula_text, fixed = TRUE)),
      n_shared_followers = sum(followers),
      cross_row_columns = paste(sort(unique(own_col[cross_row])), collapse = " ")
    )
  })
}

# Runs the formula check over the configured input workbooks, writes the
# summary, and stops the pipeline if any same-sheet cross-row formula is found.
assert_input_formula_integrity <- function(paths, out_dir, accepted = list()) {
  summary_tbl <- purrr::map_dfr(paths, excel_sheet_formula_summary)

  accepted_cols <- function(workbook, sheet) {
    for (a in accepted) {
      if (identical(a$workbook, workbook) && identical(a$sheet, sheet)) {
        return(strsplit(trimws(a$columns), "\\s+")[[1]])
      }
    }
    character(0)
  }

  summary_tbl <- summary_tbl %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      unexpected_columns = paste(
        setdiff(
          strsplit(trimws(cross_row_columns), "\\s+")[[1]],
          c("", accepted_cols(workbook, sheet))
        ),
        collapse = " "
      ),
      status = dplyr::case_when(
        is.na(n_cross_row) ~ "UNREADABLE",
        n_cross_row == 0 ~ "PASS",
        nzchar(unexpected_columns) ~ "FAIL",
        TRUE ~ "KNOWN (recomputed from raw inputs)"
      )
    ) %>%
    dplyr::ungroup()

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_tbl, file.path(out_dir, "Input_Formula_Integrity.csv"), row.names = FALSE)

  offenders <- summary_tbl %>% dplyr::filter(status == "FAIL")
  if (nrow(offenders)) {
    stop(
      "Cross-row formulas found in column(s) not recorded in ",
      "2_Config/config.yml (input_integrity$accepted_cross_row). Cached values ",
      "in these columns may belong to other records:\n  ",
      paste(sprintf(
        "%s [%s]: unexpected columns %s (%d of %d formulas cross-row)",
        offenders$workbook, offenders$sheet, offenders$unexpected_columns,
        offenders$n_cross_row, offenders$n_formulas
      ), collapse = "\n  "),
      "\nEstablish what changed before extending the accepted list; affected ",
      "values must be recomputed from same-row raw inputs (R/derive.R).",
      call. = FALSE
    )
  }

  known <- summary_tbl %>% dplyr::filter(startsWith(status, "KNOWN"))
  if (nrow(known)) {
    message(sprintf(
      "Input formula check: known cross-row damage confirmed unchanged in %s (recomputed downstream).",
      paste(sprintf("%s[%s]", known$workbook, known$sheet), collapse = ", ")
    ))
  }

  ref_errors <- summary_tbl %>% dplyr::filter(!is.na(n_ref_errors), n_ref_errors > 0)
  if (nrow(ref_errors)) {
    message(sprintf(
      "Note: %d formula cell(s) contain #REF! errors (%s). Not used by the pipeline unless a listed column is read.",
      sum(ref_errors$n_ref_errors),
      paste(sprintf("%s[%s]", ref_errors$workbook, ref_errors$sheet), collapse = ", ")
    ))
  }

  message(sprintf(
    "Input formula check: %d sheet(s) scanned, %d formulas, no unexpected cross-row references.",
    nrow(summary_tbl), sum(summary_tbl$n_formulas, na.rm = TRUE)
  ))
  invisible(summary_tbl)
}
