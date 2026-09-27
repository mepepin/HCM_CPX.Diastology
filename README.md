# Diastolic function and cardiopulmonary fitness in hypertrophic cardiomyopathy

Analysis code for:

> **Diastolic Dysfunction Predicts Cardiopulmonary Fitness and Clinical Outcomes in Hypertrophic Cardiomyopathy.**
> Pepin ME, O'Sullivan J, Gjermeni D, Santana EJ, Mosher B, Malunjkar S, Tso JV, Ashley E,
> Lewis E, Christle JW, Meyers J, Parikh V, Kawana M, Wheeler M, Haddad F.

This repository is the methodological supplement to the paper. It contains the complete
analysis pipeline and the aggregate results it produces. It does **not** contain patient
data.

The study relates three guideline diastolic indices — E/e′, left atrial volume index (LAVi),
and peak tricuspid regurgitation velocity (TRV<sub>max</sub>) — to percent-predicted peak
V̇O₂, V̇E/V̇CO₂ slope, their serial trajectories, and a composite of heart failure
hospitalization, transplantation, or death.

---

## Running the analysis

```bash
Rscript run_all.R          # every stage, in order (~2 minutes)
Rscript run_all.R 05 07    # only the listed stages
```

Run from the repository root. Each stage runs in a fresh R process and can only use what an
earlier stage wrote to disk, so stages are reproducible in isolation. Logs are written to
`7_Logs/<stage>.log` and the run stops at the first failing stage.

**Requirements.** R 4.6.1. Package versions are pinned in `renv.lock`; restore them with
`renv::restore()`. Figures are written through Cairo, and `pdftoppm` (poppler) is used to
rasterise a small number of PDF panels.

**Data.** The three source workbooks are patient-level and are not distributed here. Place
them in `1_Input/` as named in `2_Config/config.yml` to run the pipeline end to end. Without
them, the committed contents of `6_Results/` are the outputs of the most recent run.

## Repository layout

| Path | Contents |
|---|---|
| `run_all.R` | Entry point; runs the stages in order |
| `R/` | Pipeline stages (`01`–`08`) and the helper files they share |
| `2_Config/config.yml` | Every path and analysis parameter |
| `2_Config/manuscript_manifest.csv` | The list of submission files assembled by stage 08 |
| `6_Results/` | Aggregate results: tables, figures, and the numbers ledger |
| `1_Input/`, `5_Data/` | Patient-level; git-ignored and never committed |

## Pipeline stages

| Stage | Does |
|---|---|
| `01_import` | Reads the three source workbooks; records their SHA-256; checks the stored Excel formulas for cross-row references |
| `02_derive` | Recomputes age, BMI, BSA, lean mass, and the %-predicted V̇O₂ and heart-rate fields from same-row raw inputs |
| `03_cohort` | Cohort assembly: effort, comorbidity exclusion, echo pairing, HCM entry, baseline visit, longitudinal subcohort, outcomes |
| `04_cohort_table1` | Reference comparison by propensity overlap weighting; Figure 1; Table 1; Figure S1 |
| `05_cross_sectional` | Restricted cubic spline models; Figures 2, 3 and S2; per-SD estimates; bootstrap breakpoint search |
| `06_longitudinal` | Linear mixed models and the generalized additive refit; Figure 4; Figure S3 |
| `07_outcomes` | Kaplan-Meier and Cox models; proportional-hazards diagnostics; Figure 5; Figure S4 |
| `08_manuscript` | Central Illustration, the submission bundle, the numbers ledger, and a de-identification check over every exported file |

## How the outputs map to the paper

Stage 08 assembles the submission bundle into `6_Results/6_Manuscript/` from
`2_Config/manuscript_manifest.csv`, which is the single list of submission files.

| Paper item | File |
|---|---|
| Central Illustration | `main/CentralIllustration.pdf` |
| Figure 1 — cohort assembly | `main/Figure1_CohortFlow.pdf` |
| Figure 2 — diastolic phenotypes | `main/Figure2_ASE2025_FillingPressure.pdf` |
| Figure 3 — cross-sectional splines | `main/Figure3_Nonlinear_DiastolicIndices_RCS.pdf` |
| Figure 4 — longitudinal trajectories | `main/Figure4_Longitudinal_GAMM.pdf` |
| Figure 5 — heart failure outcomes | `main/Figure5_HeartFailure_Outcomes.pdf` |
| Table 1 — baseline characteristics | `main/Table1_Manuscript.docx` |
| Table 2 — end-point components | `main/Table2_EventComponents.csv` |
| Figure S1 — propensity overlap | `supplement/FigureS1_PropensityOverlap.pdf` |
| Figure S2 — all diastolic parameters | `supplement/FigureS2_AllDiastolic_RCS.pdf` |
| Figure S3 — LAVi threshold trajectories | `supplement/FigureS3_LAVi_Threshold_Trajectories.pdf` |
| Figure S4 — outcome sensitivity | `supplement/FigureS4_Outcomes_Sensitivity.pdf` |
| Tables S1–S5 | `supplement/TableS1…TableS5…` |

`6_Results/6_Manuscript/Numbers_Ledger.csv` lists every headline number in the paper next to
the result file it is read from, so no value in the text is transcribed by hand.

## Analysis design

- **Available-case, parameter-specific.** Each index is analysed in every patient in whom it
  was measured, drawn from a 441-patient base analytic cohort. The 261 patients with both
  E/e′ and LAVi form a prespecified consistency cohort, not the primary sample.
- **Multiplicity.** Benjamini-Hochberg within each of the cross-sectional, longitudinal, and
  time-to-event families; adjusted q values are reported alongside P values.
- **Missing covariates** are modelled rather than excluded: unrecorded status becomes its own
  category for binary covariates, and the resting LVOT gradient is paired with an indicator of
  whether it was recorded.
- **Comparators are fitted in parallel models**, never jointly with a diastolic index, so no
  estimate here is adjusted for another exposure.
- **Sensitivity analyses**: restriction to adults, the nested measurement cohorts, and
  multiple imputation. Multiple imputation is always a labelled sensitivity; every primary
  estimate is available-case.

Cohort sizes produced by the current code: 441 base analytic cohort, 261 with both E/e′ and
LAVi, 337 non-HCM reference subjects, 199 patients contributing 615 serial tests, and 436
patients with follow-up contributing 118 composite events.

## Reproducibility

The pipeline is deterministic: bootstrap and imputation seeds are fixed in
`2_Config/config.yml`, and re-running it reproduces every table and figure byte for byte.
Changing an analysis parameter means editing `config.yml` and re-running, not editing a
script — the cohort definitions have exactly one definition, in `R/cohort.R`, called from
stage 03 and nowhere else.

Stage 08 runs a de-identification check across every exported file and fails the run if a
patient identifier reaches `6_Results/`.

## License

See `LICENSE`.
