# Publication styling: palette, notation labels, ggplot theme, table and
# P-value formatting, and figure export (PDF/PNG, XeLaTeX-typeset legends).
# Ported verbatim from _Scripts/R/utils.R (lines 1-359).

# ── JACC color palette ────────────────────────────────────────────────────────
jacc_cols <- c(
  navy   = "#2C3E50", red    = "#E74C3C", blue   = "#3498DB",
  green  = "#27AE60", teal   = "#16A085", orange = "#F39C12",
  gray   = "#7F8C8D", purple = "#8E44AD"
)

# ── Shared cardiopulmonary notation labels ────────────────────────────────────
label_peak_vo2                  <- "Peak V̇O2"
label_peak_vo2_md               <- "Peak V̇O$_2$"
label_peak_vo2_plot             <- "Peak V̇O<sub>2</sub>"
label_peak_vo2_unicode          <- "Peak V̇O₂"
label_peak_vo2_friend           <- paste0(label_peak_vo2,         " (FRIEND 2.0 %)")
label_peak_vo2_friend_md        <- paste0(label_peak_vo2_md,      " (FRIEND 2.0 %)")
label_peak_vo2_friend_plot      <- paste0(label_peak_vo2_plot,    " (FRIEND 2.0 %)")
label_peak_vo2_friend_unicode   <- paste0(label_peak_vo2_unicode, " (FRIEND 2.0 %)")
label_peak_vo2_friend_pred_unicode <- paste0(label_peak_vo2_unicode, " (FRIEND 2.0 % predicted)")
label_peak_vo2_friend_pred      <- paste0(label_peak_vo2,         " (FRIEND 2.0 % predicted)")
label_peak_vo2_friend_pred_md   <- paste0(label_peak_vo2_md,      " (FRIEND 2.0 % predicted)")
label_peak_vo2_friend_pred_plot <- paste0(label_peak_vo2_plot,    " (FRIEND 2.0 % predicted)")
label_vevco2                    <- "V̇E/V̇CO2 slope"
label_vevco2_md                 <- "V̇E/V̇CO$_2$ slope"
label_vevco2_plot               <- "V̇E/V̇CO<sub>2</sub> slope"
label_vevco2_unicode            <- "V̇E/V̇CO₂ slope"
label_vevco2_title              <- "V̇E/V̇CO2 Slope"
label_vevco2_title_md           <- "V̇E/V̇CO$_2$ Slope"
label_vevco2_title_plot         <- "V̇E/V̇CO<sub>2</sub> Slope"
label_vevco2_title_unicode      <- "V̇E/V̇CO₂ Slope"

# ── JACC publication theme ────────────────────────────────────────────────────
theme_jacc <- function(base_size = 9) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      plot.subtitle    = element_text(size = base_size, hjust = 0.5, color = "grey40"),
      axis.title       = element_text(face = "bold"),
      axis.text        = element_text(color = "black"),
      axis.line        = element_line(color = "black", linewidth = 0.4),
      axis.ticks       = element_line(color = "black", linewidth = 0.3),
      panel.grid       = element_blank(),
      strip.background = element_rect(fill = "grey94", color = NA),
      strip.text       = element_text(face = "bold", size = base_size),
      legend.position  = "top",
      legend.title     = element_blank(),
      legend.text      = element_text(size = base_size - 1),
      plot.margin      = margin(4, 8, 4, 4)
    )
}

# ── Publication-quality flextable formatter ───────────────────────────────────
format_pub_table <- function(ft, caption = NULL) {
  thin_border <- officer::fp_border(color = "black", width = 0.8)
  ft <- ft %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    line_spacing(space = 1.0, part = "all") %>%
    padding(padding.top = 1, padding.bottom = 1, padding.left = 3, padding.right = 3, part = "body") %>%
    padding(padding.top = 2, padding.bottom = 2, padding.left = 3, padding.right = 3, part = "header") %>%
    border_remove() %>%
    hline_top(border = thin_border, part = "header") %>%
    hline_bottom(border = thin_border, part = "header") %>%
    hline_bottom(border = thin_border, part = "body") %>%
    set_table_properties(layout = "autofit") %>%
    autofit()
  if (!is.null(caption)) ft <- ft %>% set_caption(caption)
  ft
}

# ── Flextable subscript helpers ───────────────────────────────────────────────
subscript_paragraph <- function(x) {
  if (is.na(x)) return(as_paragraph(as_chunk("")))
  marker <- "<<SUB2>>"
  text <- gsub("CO2", paste0("CO", marker), x,   fixed = TRUE)
  text <- gsub("O2",  paste0("O",  marker), text, fixed = TRUE)
  parts <- strsplit(text, marker, fixed = TRUE)[[1]]
  chunk_list <- list()
  for (i in seq_along(parts)) {
    if (nzchar(parts[[i]])) chunk_list[[length(chunk_list) + 1]] <- as_chunk(parts[[i]])
    if (i < length(parts))  chunk_list[[length(chunk_list) + 1]] <- as_sub("2")
  }
  if (!length(chunk_list)) chunk_list[[1]] <- as_chunk("")
  do.call(as_paragraph, chunk_list)
}

apply_subscript_format <- function(ft, data, columns) {
  for (col in intersect(columns, names(data))) {
    idx <- which(!is.na(data[[col]]) & grepl("O2|CO2", data[[col]]))
    if (!length(idx)) next
    for (i in idx) ft <- compose(ft, i = i, j = col, value = subscript_paragraph(data[[col]][i]))
  }
  ft
}

# ── P-value formatting ────────────────────────────────────────────────────────
superscript_int <- function(x) {
  sup_map <- c("0"="⁰","1"="¹","2"="²","3"="³","4"="⁴","5"="⁵","6"="⁶","7"="⁷","8"="⁸","9"="⁹","-"="⁻","+"="⁺")
  paste0(unname(sup_map[strsplit(as.character(x), "", fixed = TRUE)[[1]]]), collapse = "")
}

scientific_p_parts <- function(p, digits = 2) {
  if (is.na(p)) return(list(mantissa = NA_character_, exponent = NA_integer_))
  if (p == 0)   return(list(mantissa = "0",           exponent = NA_integer_))
  sci      <- sprintf(paste0("%.", digits, "e"), p)
  parts    <- strsplit(sci, "e", fixed = TRUE)[[1]]
  mantissa <- sub("0+$", "", parts[1])
  mantissa <- sub("\\.$", "", mantissa)
  list(mantissa = mantissa, exponent = as.integer(parts[2]))
}

format_p_scientific <- function(p, digits = 2) {
  vapply(p, function(val) {
    if (is.na(val)) return(NA_character_)
    if (val == 0)   return("0")
    sci <- sprintf(paste0("%.", digits, "e"), val)
    parts    <- strsplit(sci, "e", fixed = TRUE)[[1]]
    mantissa <- sub("0+$", "", sub("\\.$", "", parts[1]))
    paste0(mantissa, " × 10", superscript_int(as.integer(parts[2])))
  }, character(1))
}

format_p_mixed <- function(p, sci_digits = 2, threshold = 0.05, decimal_digits = 3) {
  vapply(p, function(val) {
    if (is.na(val))       return(NA_character_)
    if (val < threshold)  return(format_p_scientific(val, digits = sci_digits))
    sprintf(paste0("%.", decimal_digits, "f"), val)
  }, character(1))
}

format_p_exact <- function(p) {
  vapply(p, function(val) {
    if (is.na(val)) return(NA_character_)
    sprintf("%.15g", val)
  }, character(1))
}

scientific_p_paragraph <- function(p, digits = 2) {
  sci <- scientific_p_parts(p, digits = digits)
  if (is.na(sci$mantissa)) return(as_paragraph(as_chunk("")))
  if (is.na(sci$exponent)) return(as_paragraph(as_chunk(sci$mantissa)))
  as_paragraph(as_chunk(sci$mantissa), as_chunk(" x 10"), as_sup(as.character(sci$exponent)))
}

apply_scientific_p_format <- function(ft, data, value_col = "P value", p_col = "p_num",
                                      digits = 2, threshold = Inf) {
  if (!all(c(value_col, p_col) %in% names(data))) return(ft)
  idx <- which(!is.na(data[[p_col]]) & data[[p_col]] < threshold)
  if (!length(idx)) return(ft)
  for (i in idx) {
    ft <- compose(ft, i = i, j = value_col,
                  value = scientific_p_paragraph(data[[p_col]][i], digits = digits))
  }
  ft
}

# ── Export helpers ────────────────────────────────────────────────────────────
save_jacc <- function(plot, path, w = 6.9, h = 4.5) {
  ggsave(filename = path, plot = plot, device = cairo_pdf,
         width = w, height = h, units = "in", dpi = 600)
  if (grepl("\\.pdf$", path, ignore.case = TRUE)) {
    ggsave(filename = sub("\\.pdf$", ".png", path, ignore.case = TRUE),
           plot = plot, device = "png", width = w, height = h,
           units = "in", dpi = 600, bg = "white")
  }
}

draw_export_object <- function(x) {
  if (inherits(x, c("grob", "gTree", "gtable"))) grid::grid.draw(x)
  else print(x, newpage = FALSE)
}

make_embedded_legend_grob <- function(text, width_in, fontfamily = "Times New Roman",
                                      fontsize = 10, lineheight = 1.25, margin_in = 0.14) {
  legend_html <- text
  legend_html <- gsub("&",  "&amp;",   legend_html, fixed = TRUE)
  legend_html <- gsub("<",  "&lt;",    legend_html, fixed = TRUE)
  legend_html <- gsub(">",  "&gt;",    legend_html, fixed = TRUE)
  legend_html <- gsub("\\*\\*(.*?)\\*\\*", "<strong>\\1</strong>", legend_html, perl = TRUE)
  legend_html <- gsub("\\$_([A-Za-z0-9]+)\\$", "<sub>\\1</sub>",  legend_html, perl = TRUE)

  usable_width_in <- width_in - 2 * margin_in
  label_grob <- gridtext::textbox_grob(
    legend_html,
    x = grid::unit(margin_in / width_in, "npc"),
    y = grid::unit(1, "npc"),
    width = grid::unit(usable_width_in, "in"),
    hjust = 0, vjust = 1, halign = 0, valign = 1,
    gp = grid::gpar(fontsize = fontsize, fontfamily = fontfamily,
                    lineheight = lineheight, col = "#111111"),
    padding = grid::unit(c(0, 0, 0, 0), "pt"),
    margin  = grid::unit(c(0, 0, 0, 0), "pt"),
    r = grid::unit(0, "pt"),
    box_gp = grid::gpar(col = NA, fill = NA)
  )
  legend_height_in <- 0.06 +
    grid::convertHeight(grid::grobHeight(label_grob), "in", valueOnly = TRUE) +
    0.10
  wrapped_grob <- grid::grobTree(
    label_grob,
    vp = grid::viewport(x = 0, y = 0, width = 1, height = 1, just = c("left", "bottom"))
  )
  list(
    grob = grid::grobTree(
      wrapped_grob,
      vp = grid::viewport(x = 0, y = 1 - 0.06 / legend_height_in,
                          width = 1, height = 1, just = c("left", "top"))
    ),
    height_in = legend_height_in
  )
}

save_jacc_with_embedded_legend <- function(plot, path, legend_text, w = 6.9, h = 4.5,
                                           legend_fontsize = 10,
                                           legend_family = "Times New Roman") {
  # gridtext can wrap rich text but cannot fully justify it. Render the legend
  # through XeLaTeX so every embedded figure legend is flush left and fully
  # justified while the final line remains left aligned. The plot itself stays
  # vector-native and is included without scaling or restyling.
  xelatex <- Sys.which("xelatex")
  if (!nzchar(xelatex)) {
    stop("XeLaTeX is required to export left-aligned, fully justified legends.")
  }

  latex_escape_legend <- function(x) {
    x <- gsub(
      "\\*\\*(.*?)\\*\\*",
      "QQBOLDOPENQQ\\1QQBOLDCLOSEQQ",
      x,
      perl = TRUE
    )
    x <- gsub(
      "\\$_([A-Za-z0-9]+)\\$",
      "QQSUBOPENQQ\\1QQSUBCLOSEQQ",
      x,
      perl = TRUE
    )
    replacements <- c(
      "\\" = "\\textbackslash{}",
      "&" = "\\&",
      "%" = "\\%",
      "$" = "\\$",
      "#" = "\\#",
      "_" = "\\_",
      "{" = "\\{",
      "}" = "\\}",
      "~" = "\\textasciitilde{}",
      "^" = "\\textasciicircum{}"
    )
    for (token in names(replacements)) {
      x <- gsub(token, replacements[[token]], x, fixed = TRUE)
    }
    x <- gsub("QQBOLDOPENQQ", "\\textbf{", x, fixed = TRUE)
    x <- gsub("QQBOLDCLOSEQQ", "}", x, fixed = TRUE)
    x <- gsub("QQSUBOPENQQ", "\\textsubscript{", x, fixed = TRUE)
    x <- gsub("QQSUBCLOSEQQ", "}", x, fixed = TRUE)
    x
  }

  out_dir <- dirname(path)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  build_dir <- tempfile("jacc_legend_")
  dir.create(build_dir)
  on.exit(unlink(build_dir, recursive = TRUE, force = TRUE), add = TRUE)

  plot_pdf <- file.path(build_dir, "figure_body.pdf")
  grDevices::cairo_pdf(plot_pdf, width = w, height = h)
  grid::grid.newpage()
  draw_export_object(plot)
  grDevices::dev.off()

  margin_in <- 0.14
  usable_width_in <- w - 2 * margin_in
  lineheight_pt <- legend_fontsize * 1.25
  tex_path <- file.path(build_dir, "figure_with_legend.tex")
  tex_source <- paste0(
    "\\documentclass[varwidth=", sprintf("%.4fin", w), ",border=0pt]{standalone}\n",
    "\\usepackage{graphicx}\n",
    "\\usepackage{fontspec}\n",
    "\\usepackage{ragged2e}\n",
    "\\setmainfont{", legend_family, "}\n",
    "\\setlength{\\parindent}{0pt}\n",
    "\\begin{document}\n",
    "\\begin{minipage}{", sprintf("%.4fin", w), "}\n",
    "\\includegraphics[width=", sprintf("%.4fin", w),
    ",height=", sprintf("%.4fin", h), "]{\\detokenize{",
    normalizePath(plot_pdf, winslash = "/", mustWork = TRUE), "}}\\par\n",
    "\\vspace{0.06in}\n",
    "\\hspace*{", sprintf("%.4fin", margin_in), "}",
    "\\begin{minipage}{", sprintf("%.4fin", usable_width_in), "}\n",
    "\\setlength{\\parindent}{0pt}\n",
    "{\\fontsize{", sprintf("%.2f", legend_fontsize), "pt}{",
    sprintf("%.2f", lineheight_pt),
    paste0(
      "pt}\\selectfont\\sloppy",
      "\\setlength{\\emergencystretch}{3em}",
      "\\setlength{\\spaceskip}{0.333em plus 0.167em minus 0.05em}",
      "\\justifying\\setlength{\\parindent}{0pt} "
    ),
    latex_escape_legend(legend_text), "\\par}\n",
    "\\end{minipage}\\par\n",
    "\\vspace{0.10in}\n",
    "\\end{minipage}\n",
    "\\end{document}\n"
  )
  writeLines(tex_source, tex_path, useBytes = TRUE)

  latex_log <- file.path(build_dir, "xelatex.log")
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(build_dir)
  status <- system2(
    xelatex,
    c("-interaction=nonstopmode", "-halt-on-error", basename(tex_path)),
    stdout = latex_log,
    stderr = latex_log
  )
  if (!identical(status, 0L)) {
    stop(
      "XeLaTeX legend export failed:\n",
      paste(readLines(latex_log, warn = FALSE), collapse = "\n")
    )
  }
  setwd(old_wd)
  file.copy(
    file.path(build_dir, "figure_with_legend.pdf"),
    path,
    overwrite = TRUE
  )
  invisible(path)
}

show_and_save_jacc <- function(plot, path, w = 6.9, h = 4.5) {
  print(plot)
  save_jacc(plot, path, w = w, h = h)
  invisible(plot)
}

# ── Annotation helper (used by both RCS and GAMM panels) ─────────────────────
add_plot_stats_box <- function(p, stats_label, x_range, y_range,
                               position = "upper_left", size = 2.5, fill_alpha = 0.9) {
  if (is.null(stats_label) || !nzchar(stats_label)) return(p)
  x_pad <- 0.04 * diff(x_range)
  y_pad <- 0.05 * diff(y_range)
  pos <- switch(
    position,
    "upper_left"  = list(x = x_range[1] + x_pad, y = y_range[2] - y_pad, hjust = 0, vjust = 1),
    "bottom_left" = list(x = x_range[1] + x_pad, y = y_range[1] + y_pad, hjust = 0, vjust = 0),
    list(x = Inf, y = -Inf, hjust = 1.02, vjust = -0.15)
  )
  p + ggtext::geom_richtext(
    data = data.frame(x = pos$x, y = pos$y, lab = stats_label),
    aes(x = x, y = y, label = lab),
    inherit.aes = FALSE,
    hjust = pos$hjust, vjust = pos$vjust,
    size = size, lineheight = 1.1, label.size = 0.2,
    fill = alpha("white", fill_alpha), color = "black", label.colour = "black"
  )
}
