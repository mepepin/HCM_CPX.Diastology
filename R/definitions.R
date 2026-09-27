# Shared analysis definitions: diastolic index sets, display labels, and ASE
# thresholds used across figures. Ported verbatim from the legacy
# `shared_defs` chunk.

diastolic_indices <- c(
  "e_e_ave",
  "e_e_lat",
  "e_e_med",
  "mv_a_dur",
  "mv_a_point",
  "mv_dec_time",
  "mv_e_a",
  "tr_max_vel",
  "la_vol_index",
  "med_peak_e_vel",
  "lat_peak_e_vel",
  "e_prime_ave"
)

label_map <- c(
  e_e_ave = "E/e' (avg)",
  e_e_lat = "E/e' (lat)",
  e_e_med = "E/e' (med)",
  mv_a_dur = "MV A dur",
  mv_a_point = "MV A point",
  mv_dec_time = "MV decel time",
  mv_e_a = "MV E/A",
  tr_max_vel = "TRVmax",
  la_vol_index = "LAVi",
  med_peak_e_vel = "Septal e'",
  lat_peak_e_vel = "Lateral e'",
  e_prime_ave = "Avg e'"
)

hcm_primary_indices <- c("e_e_ave", "la_vol_index", "tr_max_vel")

hcm_primary_labels <- c(
  e_e_ave = "E/e' (average)",
  la_vol_index = "LAVi (mL/m\u00B2)",
  tr_max_vel = "TRVmax (cm/s)"
)

hcm_primary_units <- c(
  e_e_ave = "",
  la_vol_index = " mL/m\u00B2",
  tr_max_vel = " cm/s"
)

ase_cutoffs_hcm <- tibble(
  variable = hcm_primary_indices,
  ase_cutoff = c(14, 34, 280),
  ase_direction = c(">", ">", ">")
)

ase_cutoffs_all <- tibble(variable = diastolic_indices) %>%
  left_join(ase_cutoffs_hcm, by = "variable")
