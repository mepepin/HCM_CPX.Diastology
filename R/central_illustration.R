# ── Central Illustration ──────────────────────────────────────────────────────
# JACC: Heart Failure requires a Central Illustration for original research: a
# single conceptual snapshot of the paper that must not duplicate content from
# the numbered figures, submitted as a TIF of at least 300 dpi, at least
# 7 inches wide, with lettering no smaller than 10 points at final size.
# The journal's illustrators redraw the final version, so this is the author
# submission: accurate, clean, and explicit enough to redraw from.
#
# Every number is read from a result file written by an earlier stage, so the
# illustration cannot drift from the analysis. Nothing is typed in.
#
# Colour: the three diastolic indices are identities, so they take the first
# three slots of the validated categorical palette (blue, orange, aqua), which
# pass the all-pairs colour-vision checks. Aqua sits below 3:1 contrast on
# white, so colour never carries meaning alone: every chip is paired with a
# dark text label. Peak VO2, the comparator, is neutral grey.

ci_colors <- c(
  ee = "#2a78d6", lavi = "#eb6834", trv = "#1baf7a",
  comparator = "#8a8983", ink = "#0b0b0b", ink2 = "#52514e",
  card = "#f3f3f1", rule = "#d9d8d4", surface = "#ffffff"
)

ci_read_values <- function(results_dir) {
  rd <- function(p) utils::read.csv(file.path(results_dir, p), check.names = FALSE,
                                    stringsAsFactors = FALSE)
  one <- function(x) { stopifnot(length(x) == 1, !is.na(x)); x }

  cc <- rd("2_Cohort/Cohort_Counts.csv")
  cnt <- function(q) one(as.numeric(cc$value[cc$quantity == q]))
  ps <- rd("2_Cohort/Cohort_ParameterSpecific_Counts.csv")
  pcount <- function(s, col) one(as.integer(ps[[col]][ps$sample == s]))

  ss <- rd("3_CrossSectional/Table_Sensitivity_SameSample_Adjustment.csv")
  beta <- function(outcome, param) {
    r <- ss[ss$Outcome == outcome & ss$Parameter == param &
              ss$`Analysis sample` == "Parameter available-case sample" &
              ss$Model == "Age/sex/BMI adjusted", ]
    one(r$`Beta per 1 SD higher`)
  }
  t2 <- rd("3_CrossSectional/TableS_Figure3_RCS_Diagnostics.csv")
  lrt <- rbind(rd("3_CrossSectional/Stats_CrossSectional_LRT_FRIEND2.csv"),
               rd("3_CrossSectional/Stats_CrossSectional_LRT_VeVco2.csv"))
  t3 <- rd("4_Longitudinal/Table3_Longitudinal_LMM_Primary.csv")
  lmm <- function(outcome, param, col) one(t3[[col]][t3$Outcome == outcome & t3$Parameter == param])
  cox <- rd("5_Outcomes/TableS_Figure5_Cox_Summary.csv")
  hr <- function(param) one(cox$`HR (95% CI)`[cox$Parameter == param])
  ov <- rd("5_Outcomes/Table_Manuscript_Outcome_Values.csv")
  fu <- one(ov$value[ov$quantity == "Observed time to first event or censoring, median (IQR), years"])

  list(
    n_base = cnt("Base analytic cohort (patients)"),
    n_long = cnt("Longitudinal patients"),
    n_out = cnt("Base analytic cohort with follow-up (patients)"),
    n_events = cnt("Base analytic cohort composite events"),
    fu_median = sub(" .*", "", fu),
    n_ee = pcount("E/e' (average)", "n_cross_sectional"),
    n_lavi = pcount("LAVi", "n_cross_sectional"),
    n_trv = pcount("TRVmax", "n_cross_sectional"),
    vo2 = c(ee = beta("Peak VO2 (% predicted)", "E/e'"),
            lavi = beta("Peak VO2 (% predicted)", "LAVi"),
            trv = beta("Peak VO2 (% predicted)", "TRVmax")),
    vev = c(ee = beta("VE/VCO2 slope", "E/e'"),
            lavi = beta("VE/VCO2 slope", "LAVi"),
            trv = beta("VE/VCO2 slope", "TRVmax")),
    overall_q_max = max(p.adjust(lrt$overall_p, "BH")),
    nonlin_q_min = min(t2$`Nonlinearity q (six-test family)`),
    lavi_slope = as.numeric(lmm("VE/VCO2 slope", "LAVi", "Interaction beta per year per 1 SD")),
    lavi_slope_q = lmm("VE/VCO2 slope", "LAVi", "Q value"),
    lavi_slope_n = as.integer(lmm("VE/VCO2 slope", "LAVi", "Patients")),
    vo2_traj_p_min = min(as.numeric(t3$`P value`[t3$Outcome == "Peak VO2 (% predicted)"])),
    hr_lavi = hr("LAVi Z-score"), hr_ee = hr("E/e'"), hr_trv = hr("TRVmax"),
    hr_vo2 = hr("Peak V̇O2")
  )
}

# Unicode minus for negative values; one decimal.
ci_num <- function(x, digits = 1, plus = TRUE) {
  s <- formatC(abs(x), format = "f", digits = digits)
  paste0(if (x < 0) "−" else if (plus) "+" else "", s)
}
ci_hr <- function(s) gsub("-", "–", s)  # en dash in CI ranges

draw_central_illustration <- function(v, years) {
  library(grid)
  C <- ci_colors
  VO2 <- "V̇O₂"; VEV <- "V̇E/V̇CO₂"
  LAB <- c(ee = "E/e′", lavi = "LAVi", trv = "TRVmax")
  txt <- function(label, x, y, size = 10, face = "plain", col = C[["ink"]],
                  just = "left", ...) {
    grid.text(label, x = unit(x, "in"), y = unit(y, "in"), just = just,
              gp = gpar(fontfamily = "Arial", fontsize = size, fontface = face,
                        col = col, lineheight = 1.1), ...)
  }
  card <- function(x, y, w, h) {
    grid.roundrect(unit(x, "in"), unit(y, "in"), unit(w, "in"), unit(h, "in"),
                   just = c("left", "bottom"), r = unit(0.06, "in"),
                   gp = gpar(fill = C[["card"]], col = NA))
  }
  chip <- function(x, y, key, d = 0.10) {
    grid.roundrect(unit(x, "in"), unit(y, "in"), unit(d, "in"), unit(d, "in"),
                   just = c("left", "center"), r = unit(0.02, "in"),
                   gp = gpar(fill = C[[key]], col = NA))
  }
  arrow_h <- function(x0, x1, y) {
    grid.lines(unit(c(x0, x1), "in"), unit(c(y, y), "in"),
               arrow = arrow(length = unit(0.07, "in"), type = "closed"),
               gp = gpar(col = C[["ink2"]], fill = C[["ink2"]], lwd = 1.2))
  }

  grid.newpage()
  grid.rect(gp = gpar(fill = C[["surface"]], col = NA))

  # Geometry (inches). Page 7.5 x 5.6; 0.2 in side margins.
  top <- 4.95; bot <- 1.10; h <- top - bot
  cx <- 0.20; cw <- 1.25                    # cohort strip
  w <- 1.70; xs <- c(1.67, 3.59, 5.51)      # three finding cards
  value_x <- function(x) x + 1.56           # right edge for numbers in a card

  # Title
  txt("Central Illustration", cx, 5.43, size = 11, face = "bold", col = C[["ink2"]])
  txt("Resting Diastolic Indices, Cardiopulmonary Fitness, and Heart Failure Outcomes in HCM",
      cx, 5.19, size = 12, face = "bold")

  # ── Cohort strip ──
  card(cx, bot, cw, h)
  x <- cx + 0.12
  txt(format(v$n_base, big.mark = ","), x, top - 0.40, size = 20, face = "bold")
  txt("patients\nwith HCM", x, top - 0.80, size = 10)
  txt(sprintf("Echo-aligned\nCPET,\n%s", years), x, top - 1.42, size = 10, col = C[["ink2"]])
  txt("Resting indices", x, top - 2.02, size = 10, face = "bold")
  for (i in 1:3) {
    k <- c("ee", "lavi", "trv")[i]
    n <- c(ee = v$n_ee, lavi = v$n_lavi, trv = v$n_trv)[[k]]
    yy <- top - 2.36 - (i - 1) * 0.42
    chip(x, yy + 0.07, k)
    txt(LAB[[k]], x + 0.17, yy + 0.07, size = 10, face = "bold", just = c("left", "center"))
    txt(sprintf("n = %d", n), x + 0.17, yy - 0.12, size = 10, col = C[["ink2"]], just = c("left", "center"))
  }

  # ── Three finding cards ──
  for (x in xs) card(x, bot, w, h)
  for (x in c(cx + cw, xs[1] + w, xs[2] + w)) arrow_h(x + 0.03, x + 0.19, top - 0.26)
  heads <- list(
    c("1  Exercise capacity", "Per SD higher index"),
    c("2  Change over time", sprintf("%d patients, serial CPET", v$n_long)),
    c("3  Heart failure events", sprintf("%d patients, %d events", v$n_out, v$n_events))
  )
  for (i in 1:3) {
    txt(heads[[i]][1], xs[i] + 0.12, top - 0.24, size = 11, face = "bold")
    txt(heads[[i]][2], xs[i] + 0.12, top - 0.47, size = 10, col = C[["ink2"]])
  }

  # Card 1: per-SD associations as numbers (Figure 3 holds the curves)
  x <- xs[1] + 0.12
  block <- function(y, title, unit_lab, vals) {
    txt(title, x, y, size = 10, face = "bold")
    txt(unit_lab, x, y - 0.20, size = 10, col = C[["ink2"]])
    for (i in 1:3) {
      k <- c("ee", "lavi", "trv")[i]
      yy <- y - 0.46 - (i - 1) * 0.23
      chip(x, yy, k, d = 0.09)
      txt(LAB[[k]], x + 0.15, yy, size = 10, just = c("left", "center"))
      txt(ci_num(vals[[k]]), value_x(xs[1]), yy, size = 10, face = "bold", just = c("right", "center"))
    }
  }
  block(top - 0.85, paste0("Lower peak ", VO2), "% predicted", v$vo2)
  block(top - 2.06, paste0("Higher ", VEV, " slope"), "slope units", v$vev)
  txt(sprintf("All 6 associations\nq ≤ %.3f; linear,\nno threshold", v$overall_q_max),
      x, bot + 0.30, size = 10, col = C[["ink2"]])

  # Card 2: schematic of trajectory modification (not data; Figure 4 holds data)
  x <- xs[2] + 0.12
  pushViewport(viewport(x = unit(x + 0.20, "in"), y = unit(top - 1.78, "in"),
                        width = unit(1.00, "in"), height = unit(0.92, "in"),
                        just = c("left", "bottom")))
  grid.lines(c(0, 0, 1), c(1, 0, 0), gp = gpar(col = C[["rule"]], lwd = 1))
  grid.lines(c(0.05, 0.95), c(0.30, 0.40), gp = gpar(col = C[["comparator"]], lwd = 2))
  grid.lines(c(0.05, 0.95), c(0.36, 0.90), gp = gpar(col = C[["lavi"]], lwd = 2))
  popViewport()
  txt("Higher\nLAVi", x + 1.25, top - 0.97, size = 10, just = c("left", "center"))
  txt("Lower\nLAVi", x + 1.25, top - 1.43, size = 10, just = c("left", "center"))
  txt(VEV, x + 0.08, top - 1.32, size = 10, col = C[["ink2"]], rot = 90, just = "center")
  txt("Time (schematic)", x + 0.70, top - 1.92, size = 10, col = C[["ink2"]], just = "center")
  txt(sprintf("Higher LAVi: steeper\nrise in %s slope", VEV), x, top - 2.30, size = 10, face = "bold")
  txt(sprintf("%s per year per SD\nq = %s; n = %d", ci_num(v$lavi_slope, 2), v$lavi_slope_q, v$lavi_slope_n),
      x, top - 2.78, size = 10)
  # "linear" matters: the mixed models found no linear modification, while the
  # generalized additive models suggested nonlinear modification by E/e' and
  # TRVmax, which the manuscript reports.
  txt(sprintf("Peak %s trajectory:\nno linear modification\nby any index", VO2), x, bot + 0.30, size = 10, col = C[["ink2"]])

  # Card 3: hazard ratios as numbers (Figure 5 holds the plot)
  x <- xs[3] + 0.12
  txt("Hazard ratio per SD", x, top - 0.85, size = 10, face = "bold")
  split_hr <- function(s) regmatches(s, regexec("^([0-9.]+) \\((.*)\\)$", s))[[1]]
  items <- list(c("lavi", v$hr_lavi), c("ee", v$hr_ee), c("trv", v$hr_trv))
  for (i in seq_along(items)) {
    k <- items[[i]][1]; hp <- split_hr(items[[i]][2]); yy <- top - 1.17 - (i - 1) * 0.44
    chip(x, yy, k, d = 0.09)
    txt(LAB[[k]], x + 0.15, yy, size = 10, just = c("left", "center"))
    txt(hp[2], value_x(xs[3]), yy, size = 12, face = "bold", just = c("right", "center"))
    txt(ci_hr(hp[3]), value_x(xs[3]), yy - 0.19, size = 10, col = C[["ink2"]], just = c("right", "center"))
  }
  grid.lines(unit(c(x, value_x(xs[3])), "in"), unit(rep(top - 2.43, 2), "in"),
             gp = gpar(col = C[["rule"]], lwd = 1))
  hp <- split_hr(v$hr_vo2)
  chip(x, top - 2.62, "comparator", d = 0.09)
  txt(sprintf("Lower peak %s", VO2), x + 0.15, top - 2.62, size = 10, just = c("left", "center"))
  txt(hp[2], value_x(xs[3]), top - 2.88, size = 12, face = "bold", just = c("right", "center"))
  txt(ci_hr(hp[3]), value_x(xs[3]), top - 3.07, size = 10, col = C[["ink2"]], just = c("right", "center"))
  txt(sprintf("Median follow-up\n%s years", v$fu_median), x, bot + 0.30, size = 10, col = C[["ink2"]])

  # Takeaway band
  grid.lines(unit(c(cx, 7.30), "in"), unit(c(0.95, 0.95), "in"), gp = gpar(col = C[["rule"]], lwd = 1))
  txt(paste0(
    "Resting diastolic indices relate to exercise capacity and ventilatory efficiency in HCM. Higher LAVi was\n",
    "associated with a steeper rise in ", VEV, " slope and with heart failure events; peak ", VO2,
    " remained the\nstrongest single predictor."),
    cx, 0.52, size = 10, just = c("left", "center"))
}

build_central_illustration <- function(results_dir, out_dir, years) {
  v <- ci_read_values(results_dir)
  width <- 7.5; height <- 5.6
  tif <- file.path(out_dir, "Central_Illustration.tif")
  pdf_path <- file.path(out_dir, "Central_Illustration.pdf")
  ragg::agg_tiff(tif, width = width, height = height, units = "in", res = 600,
                 compression = "lzw", background = "white")
  draw_central_illustration(v, years)
  grDevices::dev.off()
  grDevices::cairo_pdf(pdf_path, width = width, height = height, family = "Arial")
  draw_central_illustration(v, years)
  grDevices::dev.off()
  invisible(v)
}

# Legend: title plus a 2-3 sentence caption, then abbreviations in alphabetical
# order, as JACC requires. Numbers come from the same values as the figure.
write_central_illustration_legend <- function(v, path) {
  hr_short <- function(s) sub(" .*", "", s)
  caption <- paste0(
    "Among ", v$n_base, " patients with HCM who underwent echocardiography-aligned CPET, ",
    "higher E/e′, LAVi, and TRVmax were each associated with lower percent-predicted peak ",
    "V̇O2 and a higher V̇E/V̇CO2 slope (values per SD, adjusted for age, sex, and ",
    "body mass index). Among patients with serial CPET, higher baseline LAVi was associated with ",
    "a steeper subsequent rise in V̇E/V̇CO2 slope (the trajectory panel is schematic), and ",
    "over a median of ", v$fu_median, " years each index was associated with the composite of ",
    "heart failure hospitalization, transplantation, or death; lower peak V̇O2 ",
    "(HR ", hr_short(v$hr_vo2), " per SD) remained the strongest single predictor. ",
    "Colored squares identify each diastolic index, gray identifies peak V̇O2, and arrows ",
    "indicate the sequence of analyses."
  )
  abbrev <- paste(
    "CI = confidence interval",
    "CPET = cardiopulmonary exercise testing",
    "E/e′ = ratio of early mitral inflow velocity to early diastolic mitral annular velocity",
    "HCM = hypertrophic cardiomyopathy",
    "HR = hazard ratio",
    "LAVi = left atrial volume index",
    "q = false discovery rate-adjusted P value",
    "SD = standard deviation",
    "TRVmax = peak tricuspid regurgitation velocity",
    "V̇E/V̇CO2 = minute ventilation to carbon dioxide production relationship",
    "V̇O2 = oxygen uptake",
    sep = "; "
  )
  writeLines(c(
    "# Central Illustration legend",
    "",
    "Generated by R/08_manuscript.R from the same values as the figure.",
    "",
    "**Central Illustration. Resting Diastolic Indices, Cardiopulmonary Fitness, and Heart Failure Outcomes in Hypertrophic Cardiomyopathy**",
    "",
    caption,
    "",
    paste0(abbrev, ".")
  ), path)
  invisible(caption)
}
