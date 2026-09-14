# HCM Diastology Manuscript: Comprehensive Revision Plan

Last reviewed: 2026-09-13

## Stage status

### Stage 1 — missingness and cohort-selection audit: complete

The audit is implemented as aggregate-only, reproducible code in
`_Scripts/R/missingness_audit.R`. It does not export identifiers and does not
edit or rerender the submission Word manuscript.

Key findings:

- The upstream flow is 440 phenotypic HCM patients, 316 with E/e' measured,
  and 261 with both E/e' and LAVi measured. The 261-patient cohort is the
  parent population for the longitudinal and outcomes subcohorts.
- TRVmax is measured in 145/261 (55.6%) and missing in 116/261 (44.4%).
- TRVmax measurement is associated with observed baseline characteristics
  (likelihood-ratio P=0.0155), although the model has weak explanatory power
  (McFadden R2=0.0438). Higher LAVi is the clearest predictor of measurement.
- The largest standardized differences for TRVmax measured versus missing are
  LAVi (0.419), V̇E/V̇CO2 slope (0.334), resting LVOT gradient (0.226), and
  female sex (0.207). Requiring TRVmax therefore selects a meaningfully
  different subset.
- Cross-sectional effect estimates are directionally similar after requiring
  TRVmax, but confidence intervals widen. The LAVi–V̇E/V̇CO2 association changes
  from P=0.0137 at N=261 to P=0.0778 at N=145.
- More concerningly, all four E/e'/LAVi-by-time interactions become nominally
  significant only in the 63-patient TRVmax-complete longitudinal subset,
  whereas none is significant in the 105-patient primary longitudinal cohort.
  This is evidence of cohort-selection sensitivity, not stronger evidence for
  requiring TRVmax.
- In outcomes models, E/e' and LAVi estimates remain directionally consistent
  in the TRVmax-complete subset, but events fall from 62 to 36.

Stage 1 decision: retain the E/e' + LAVi cohort as primary and analyze TRVmax
available-case with its denominator stated everywhere. Retain the
TRVmax-complete cohort only as a selection sensitivity analysis. Do not use a
missing-indicator term to present TRVmax as an N=261 exposure.

Audit outputs are in `2_Output/Section_06_Ancillary/`, beginning with
`Stage1_Missingness_Selection_Summary.csv` and
`Table_AnalyticCohort_RequirementFlow.csv`.

### Stage 2 — cross-sectional model presentation: complete

The primary cross-sectional hierarchy is now locked in the analysis code:

- The overall RCS P value correctly compares the age/sex/BMI-only model with
  the full spline model. The former implementation tested only the linear
  component and was relabeled and corrected before manuscript presentation.
- The nonlinearity q value is Benjamini-Hochberg adjusted across the complete
  prespecified family of six primary parameter-outcome comparisons.
- Figure 3 retains its original panel geometry, axes, typography, colors, and
  dimensions. Its boxes now report only N, overall P, and nonlinearity q.
- Piecewise overlays were removed from the primary figure because none of the
  six models shows multiplicity-adjusted evidence of nonlinearity. Breakpoint,
  Delta AIC, and curvature results remain in a supplemental diagnostics table.
- Automated checks confirm N=261 for E/e' and LAVi and N=145 for available-case
  TRVmax in both primary outcomes. Same-sample crude/adjusted checks also pass.

Stage 2 decision: describe the significant ventilatory-efficiency findings as
associations without evidence that a nonlinear curve or new threshold is
required. Do not promote a data-selected breakpoint.

### Stage 3 — comparator definition and matching: complete

The reference arm is now explicitly defined as phenotypically screened non-HCM
CPET comparators, not healthy volunteers. Any patient with an HCM diagnosis,
registry entry, or phenotype flag is excluded rather than reassigned to the
reference arm. Eligibility additionally requires no exclusionary diagnosis,
LVEF >=50%, maximum LV wall thickness <1.3 cm, and resting LVOT gradient <30
mm Hg.

The outcome-blind candidate grid evaluated propensity-score matching and
Mahalanobis matching within a propensity-score caliper, always 1:1 without
replacement and exact on sex. The prespecified rule retained the largest sample
for which every age/sex/BMI post-match absolute SMD was <0.10. The selected
specification was age/BMI Mahalanobis matching within a 0.25-SD propensity-score
caliper, yielding 177 HCM patients and 177 references. Post-match SMDs were
0.046 for age, 0 for sex, and 0.097 for BMI.

Comparator QC resolved the LVOT-gradient concern: among the 242 complete
matching-pool references, zero have an HCM flag, zero have resting LVOT
gradient >=30 mm Hg, and the maximum resting gradient is 10.67 mm Hg. Therefore
replacement with apical HCM is neither necessary nor appropriate; apical HCM
remains an HCM morphology subgroup rather than a control condition.

Stage 3 decision: freeze the balanced 177:177 comparison for Figure 1 and Table
1. Refer to this group as non-HCM references/comparators, state that gradients
are resting values, and do not describe the group as healthy controls. The
supplemental matched-cohort validation figure retains its original panel design;
its legend is now exported full width, flush left, and fully justified.

### Stage 4 — longitudinal presentation: complete

The longitudinal analysis has been rerun and audited within the fixed E/e' +
LAVi parent cohort. The primary inferential presentation is a compact linear
mixed-model table. The existing six-panel GAMM design is preserved unchanged
and copied to the supplement, and a separate supplemental observed-trajectory
figure shows the underlying within-patient measurements without a fitted
smoother.

The random-effects structure was selected consistently by CPET outcome, not
separately for each diastolic parameter. Random slopes were accepted only when
both primary E/e' and LAVi models converged without warnings, estimated a
nonzero slope variance, improved AIC by at least 2, and had a descriptive
likelihood-ratio P<0.05. This supported a patient random intercept for peak
V̇O2 and independent patient random intercepts and slopes for V̇E/V̇CO2.

Stage 4 decision: all six primary parameter-by-time interactions are
nonsignificant under the selected structures. Do not claim that a baseline LV
diastolic parameter predicts subsequent change in CPET performance.

### Stage 5 — clinical-outcomes refinement: complete

The primary outcome remains the first dated acute heart failure event, heart
transplantation, or all-cause death. The supplied outcomes sheet contains a
dated `post_acute_heart_failure` field, with complete agreement between its
binary flag and date across the full source extract, but contains no
hospitalization or admission field and no endpoint data dictionary. Therefore
the defensible term is **acute heart failure event**, not heart-failure
hospitalization. Hospitalization wording requires separate source documentation
or chart adjudication.

The V̇E/V̇CO2 proportional-hazards violation has been addressed with a single
exposure-by-log-time Cox sensitivity. The association is
time-dependent: HR 0.99 per SD at 1 year (95% CI 0.77-1.27; P=0.918), HR 1.40
at 3 years (95% CI 1.09-1.79; P=0.0079), and HR 1.64 at 5 years (95% CI
1.22-2.22; P=0.0012); exposure-by-log-time P=0.0014. The single constant HR
must not be used as the definitive V̇E/V̇CO2 prognostic estimate.

Figure 5 retains its established design: panels A and B show E/e' and LAVi,
panel C uses the V̇E/V̇CO2 cutoff of 30, panel D is the peak-V̇O2 survival
curve, and panel E contains the Cox estimates and separate surgical-myectomy
models. Event-component counts are in the embedded, flush-left, justified
legend and in a supplemental table rather than a separate panel.

Stage 5 decision: retain E/e' and LAVi as the primary prognostic results,
TRVmax as available case, V̇E/V̇CO2 as a time-varying secondary association,
and surgical myectomy as a separate secondary outcome. Do not add AF, VT/VF,
sudden death, or unvalidated ablation categories to the composite.

### Stage 6A — authoritative results ledger and denominator reconciliation: complete

The publication results ledger is available at
`2_Output/HCM_Revision_Results_Ledger.xlsx`. It records the locked cohort and
modeling decisions, primary and prespecified sensitivity results, operational
definitions, source-output provenance, and a denominator/event reconciliation.

All 17 prespecified reconciliation checks match their current analysis sources.
These include the cross-sectional parent cohort (N=261), TRVmax available-case
subset (N=145), longitudinal subcohort (105 patients; 294 peak-V̇O2 and 293
V̇E/V̇CO2 observations), outcomes subcohort (N=259; 62 events), TRVmax outcomes
subset (N=144; 36 events), 177:177 matched comparison, and the complete first-event
component sum (57 acute-HF events, 0 transplants, and 5 deaths; no same-day ties).
The workbook contains no formula errors and passes file-integrity and visual
checks.

### Stage 6B — final figure audit: complete

The full audit is documented in `2_Output/Stage6B_Figure_Audit.md`. All current
PDFs have embedded fonts. The latest main PDFs have no visible clipping,
overlapping axis text, misplaced velocity dot, or overlapping oxygen subscript.
Figure 1 and Figure 3 are ready to freeze. Figures 2, 4, and 5 require only
limited label, ordering, denominator-disclosure, or placement corrections; no
redesign is indicated. The most important content correction is to label the
Figure 5 peak-V̇O2 hazard ratio as **per SD lower peak V̇O2**, because the
modeled exposure is inverse-coded.

The supplemental audit identified legacy rendering and packaging issues:
literal `<sub>` text and `Normal/Abnormal` labels in the binary-trajectory
figures, crowded summary tables, an absent manifest-listed OMARX Figure S2,
and inconsistent supplemental numbering/embedded legends. These must be
resolved before the figure gate is marked passed.

### Stage 6B.1 — design-preserving figure corrections: complete

The limited corrections enumerated in `Stage6B_Figure_Audit.md` were applied
without changing the established panel geometry, colors, dimensions, or overall
visual design. Figure 2 prevalence rows now descend by frequency and use
age-calibrated e' terminology. Figure 4 uses LV diastolic-parameter terminology.
Figure 5 explicitly labels lower peak V̇O2, resting LVOT gradient, the VE/VCO2
cutoff of 30, and every forest-model denominator. Binary trajectory figures now
use direct threshold groups rather than Normal/Abnormal, render oxygen subscripts
correctly, and have non-overlapping summary tables. Embedded legends are flush
left and justified. The matched-comparator figure uses Reference/non-HCM
terminology and resting-gradient language.

The optional OMARX figures were removed from the submission manifest because
the expected clinical-summary figure was absent and the exploratory comparison
was not needed for the locked primary inference. Existing exploratory output was
moved to `2_Output/Optional_Exploratory_OMARX/` rather than deleted. The stale
matched-comparator duplicate was moved out of `Manuscript_Main` to
`2_Output/Archive_Prior_Figure_Placement/`. Final PDFs were rerendered, visually
audited, and confirmed to contain embedded fonts.

### Stage 6C — Paperpile-preserving manual manuscript revision: checklist ready

The section-by-section manual editing plan is available in
`MANUSCRIPT_MANUAL_EDIT_CHECKLIST.md` and is keyed to the locked results ledger.
Manual execution in Word is the next stage. Do not programmatically edit,
regenerate, or rerender the submission Word manuscript.

## Non-negotiable manuscript restriction

Do **not** programmatically edit, regenerate, or rerender the submission Word
manuscript. All manuscript edits must be made manually so Paperpile citation
fields and reference formatting are preserved. Analysis code, CSV outputs, and
figure PDFs may be regenerated independently.

## 1. Decisions now fixed

1. The primary cross-sectional cohort requires measured **E/e' and LAVi**.
2. TRVmax is **not** required for cohort entry. It is reported and modeled using
   its observed denominator.
3. Missing TRVmax is never coded as below threshold or normal.
4. E/e', LAVi, and TRVmax are modeled separately; they are not entered together
   in a primary model.
5. Primary adjustment is limited to age, sex, and BMI.
6. Structural variables, resting LVOT gradient, medications, hypertension, and
   morphology are sensitivity covariates rather than primary adjustment terms.
7. Longitudinal and clinical-outcome populations are nested subcohorts of the
   E/e' + LAVi cross-sectional cohort.
8. The primary clinical endpoint remains the first date-validated acute heart
   failure event, heart transplantation, or all-cause death.
9. Surgical myectomy remains a separate secondary endpoint and is not merged
   into the primary composite.
10. The cohort decision remains reversible through the single setting in
    `_Scripts/README_HCM_Manuscript_V2.qmd`:

    ```r
    analytic_cohort_mode <- "ee_lavi"
    ```

    Changing this to `"complete_three"` restores the E/e' + LAVi + TRVmax
    complete-case analysis.

## 2. Authoritative current populations

Use these counts consistently after the E/e' + LAVi decision:

- Phenotypic HCM entry cohort: **440 patients**.
- Cross-sectional analytic cohort: **261 patients** with E/e' and LAVi.
- TRVmax measured within the cross-sectional cohort: **145/261 patients**.
- Longitudinal subcohort after intervention censoring: **105 patients and 294
  CPETs**; 293 V̇E/V̇CO2 observations.
- TRVmax longitudinal subset: **63 patients and 176 CPETs**; 175 V̇E/V̇CO2
  observations.
- Clinical-outcomes cohort: **259 patients and 62 primary events**.
- TRVmax outcomes subset: **144 patients and 36 primary events**.
- Known-morphology cross-sectional subset: **225 patients**: 51 apical and 174
  known nonapical; morphology is unknown in 36/261.
- Known-morphology outcomes subset: **223 patients and 53 events**.
- Phenotypically screened non-HCM comparator pool: **243 patients**, of whom
  **242** have complete age, sex, and BMI for matching.
- Final balanced comparison: **177 HCM and 177 references** (total N=354).
  Post-match SMDs are 0.046 for age, 0 for sex, and 0.097 for BMI.

The active definition and N are regenerated in:

`2_Output/Section_06_Ancillary/Analytic_Cohort_Definition.csv`

## 3. Current results that guide the revision

### Cross-sectional associations

The age/sex/BMI-adjusted restricted cubic spline models use N=261 for E/e' and
LAVi and N=145 for TRVmax. The corrected overall tests compare the full spline
with the covariate-only model.

- Peak V̇O2: no clear association for E/e' (overall P=0.256), LAVi (P=0.530),
  or TRVmax (P=0.335).
- V̇E/V̇CO2 slope: associations for E/e' (overall P=0.0039), LAVi (P=0.0158),
  and TRVmax in its available-case subset (P<0.001).
- Nonlinearity q values are nonsignificant for all six relationships (range
  0.263-0.898). The principal message is therefore association, not discovery
  of a nonlinear threshold.
- Same-sample linear sensitivity models give the same practical conclusion:
  age/sex/BMI-adjusted beta per SD for V̇E/V̇CO2 slope is +1.41 for E/e'
  (P<0.001), +0.94 for LAVi (P=0.014), and +2.87 for TRVmax (P<0.001).
- Expanded adjustment reduces N to 111 for E/e'/LAVi and 60 for TRVmax. These
  models are too denominator-dependent to serve as the main analysis.

### Two-parameter phenotype description

Within the 261-patient cohort:

- Neither E/e' nor LAVi above threshold: **116**.
- E/e' alone above threshold: **27**.
- LAVi alone above threshold: **61**.
- Both above threshold: **57**.
- E/e' >14: **84/261 (32.2%)**.
- LAVi >34 mL/m2: **118/261 (45.2%)**.
- TRVmax >2.8 m/s: **17/145 (11.7%)** among those measured.

The two-parameter groups differ for V̇E/V̇CO2 slope (Kruskal-Wallis P=0.0042)
but not peak V̇O2 (P=0.57).

### Morphology sensitivity

- Adding apical-versus-known-nonapical morphology does not materially alter the
  E/e' or LAVi association with V̇E/V̇CO2 slope.
- In the known-morphology cross-sectional subset, E/e' remains associated with
  V̇E/V̇CO2 slope (P=0.002) and LAVi remains associated (P=0.002).
- Morphology should remain a sensitivity covariate. The apical subgroup is not
  large enough for stable parameter-by-phenotype interaction testing.

### Longitudinal associations

The parsimonious linear mixed models include continuous time, standardized
baseline parameter, parameter-by-time interaction, age, sex, and BMI. The
outcome-level random-effects audit selected a patient random intercept for peak
V̇O2 and independent patient random intercepts and slopes for V̇E/V̇CO2.

- E/e' and LAVi models use 105 patients/294 CPETs for peak V̇O2 and 105/293 for
  V̇E/V̇CO2 slope.
- No interaction meets P<0.05: E/e' beta=1.97 percentage points/year/SD
  (95% CI -0.14 to 4.08; P=0.067) and LAVi beta=-0.78 (95% CI -2.95 to
  1.39; P=0.479) for peak V̇O2.
- For V̇E/V̇CO2 slope, E/e' beta=0.14 units/year/SD (95% CI -0.16 to 0.45;
  P=0.359) and LAVi beta=0.24 (95% CI -0.06 to 0.55; P=0.114). The previous
  random-intercept-only LAVi P=0.0503 is not the primary result because the
  outcome-level audit supports random slopes.
- TRVmax is available case: 63 patients/176 peak-V̇O2 observations and 63/175
  V̇E/V̇CO2 observations. Its interactions are nonsignificant (P=0.326 and
  P=0.512, respectively).
- Corresponding GAMM interaction tests are also nonsignificant (P range
  0.075-0.603).
- Requiring TRVmax yields a selected 63-patient subset in which several E/e'
  and LAVi interactions become nominally significant. These are explicitly
  cohort-selection sensitivities and do not supersede the fixed primary
  105-patient analysis.

Use the linear mixed models for primary longitudinal inference. Retain GAMMs as
flexible sensitivity analyses, retain raw observed trajectories in the
supplement, and avoid claiming longitudinal worsening.

### Clinical outcomes

In the fixed E/e' + LAVi outcomes cohort of 259 patients/62 events:

- E/e': adjusted HR 1.29 per SD (95% CI 1.09-1.53; P=0.003).
- LAVi: adjusted HR 1.41 per SD (95% CI 1.14-1.73; P=0.001).
- TRVmax, available-case N=144/36 events: adjusted HR 1.25 per SD
  (95% CI 0.97-1.61; P=0.084).
- The global proportional-hazards tests are nonsignificant for E/e', LAVi, and
  TRVmax.
- The secondary V̇E/V̇CO2 Cox model violates proportional hazards (global
  P=0.003; exposure-specific P<0.001). A log-time-varying coefficient model
  resolves this: HR 0.99 at 1 year, 1.40 at 3 years, and 1.64 at 5 years;
  exposure-by-log-time P=0.0014. Present these time-specific estimates rather
  than the nonsignificant constant HR 1.11 as the prognostic interpretation.

Primary event composition:

- Acute heart failure: **57 first events**.
- Heart transplantation: **0 first events; 2 occurrences**.
- All-cause death: **5 first events; 6 occurrences**.
- Total: **62/259**.
- Date-validated surgical myectomy: **33 occurrences**, analyzed separately.
- The LVOT-gradient-adjusted surgical-myectomy models include 189 patients/19
  events for E/e' and LAVi and 108/11 for TRVmax; all three parameter estimates
  are nonsignificant.
- The source extract does not establish that acute-HF events were
  hospitalizations. Use “dated acute heart failure event” unless external
  endpoint documentation becomes available.

## 4. Remaining analysis work, in order

### A. Missingness and selection audit — completed in Stage 1

1. Produce a missingness table for every variable used in a main or sensitivity
   model, overall and by outcome/event status.
2. Compare the 261 E/e' + LAVi cohort with:
   - patients excluded for missing E/e' or LAVi;
   - the 145 with measured TRVmax;
   - the 116 without measured TRVmax.
3. Report standardized differences for age, sex, BMI, HCM morphology, wall
   thickness, resting LVOT gradient, peak V̇O2, V̇E/V̇CO2 slope, and clinical
   events where available.
4. Model `TRVmax measured` as a descriptive missingness outcome using a
   parsimonious logistic model. This is a selection diagnostic, not proof of a
   missing-at-random mechanism.
5. Retain the complete-three analysis as a prespecified sensitivity. Explicitly
   compare its effect estimates with the E/e' + LAVi primary analysis.

### B. Cross-sectional model presentation

1. Keep parameter-specific age/sex/BMI-adjusted RCS models as the primary
   cross-sectional analysis.
2. Display E/e' and LAVi at N=261; display TRVmax at N=145 in its panel and
   table row.
3. Keep the same-sample unadjusted-versus-adjusted linear table in the
   supplement to demonstrate that covariate adjustment itself does not cause
   denominator loss.
4. Move expanded structural/medication models to the supplement.
5. Simplify the Figure 3 annotation to N, overall P, and nonlinearity Q. Move
   Delta AIC, selected knots, and piecewise diagnostics to supplemental tables.

### C. Comparator matching

Completed in Stage 3.

1. The outcome-blind grid used exact sex plus propensity-score or age/BMI
   Mahalanobis matching with prespecified propensity-score calipers.
2. The largest passing specification retained 177 pairs using age/BMI
   Mahalanobis distance within a 0.25-SD propensity-score caliper.
3. All post-match absolute SMDs are <0.10: age 0.046, sex 0, BMI 0.097.
4. Of 242 complete matching-pool references, 177 were retained; of the 190
   eligible HCM candidates, 177 were retained.
5. All comparator gradients are resting measurements; zero are >=30 mm Hg and
   the maximum is 10.67 mm Hg in the eligible pool and 10.29 mm Hg in the
   matched reference arm.
6. Reproducible audit outputs are in `2_Output/Section_06_Ancillary/`:
   `Table_Matching_Candidate_Grid.csv`,
   `Table_Matching_Selected_Specification.csv`,
   `Table_Comparator_Eligibility_Flow.csv`, and
   `Table_Matched_Comparator_QC.csv`.

### D. Longitudinal analysis — completed in Stage 4

1. The linear mixed model is the main inferential model.
2. Random-effects selection was performed by CPET outcome using identical
   analysis frames. Peak V̇O2 retains a random intercept: neither primary model
   supported random slopes (E/e' delta AIC=-2.01, descriptive LRT P=0.922;
   LAVi delta AIC=-1.65, P=0.555).
3. V̇E/V̇CO2 uses an independent random intercept and slope: both primary
   models supported this structure (E/e' delta AIC=+5.76, descriptive LRT
   P=0.0053; LAVi delta AIC=+4.02, P=0.0141).
4. TRVmax follows the structure selected for its CPET outcome rather than
   receiving a parameter-specific model choice.
5. The primary table reports parameter-by-time beta, 95% CI, P value, patient
   N, observation N, analysis sample, and random-effects structure.
6. The six-panel GAMM figure is retained unchanged as a flexible supplemental
   analysis. A four-panel supplemental figure adds raw within-patient
   trajectories with no fitted model or inferential smoother.
7. Reproducible outputs are:
   - `2_Output/Table3_Longitudinal_LMM_Primary.csv` and `.docx`;
   - `2_Output/Section_04_Longitudinal/TableS_Longitudinal_RandomStructure_Audit.csv`;
   - `2_Output/Section_04_Longitudinal/TableS_Longitudinal_RandomStructure_Decision.csv`;
   - `2_Output/Section_04_Longitudinal/TableS_Longitudinal_VisitDistribution.csv`;
   - `2_Output/Section_06_Ancillary/Table_CohortSelection_Longitudinal_EeLAVi.csv`;
   - `2_Output/Manuscript_Supplemental/FigureS_Longitudinal_GAMM.pdf`; and
   - `2_Output/Manuscript_Supplemental/FigureS_Longitudinal_ObservedTrajectories.pdf`.

### E. Outcomes analysis — completed in Stage 5

1. E/e' and LAVi Cox models use the identical 259-patient/62-event cohort;
   TRVmax is available case at N=144/36 events.
2. Event components are reported as both non-mutually-exclusive occurrences
   and mutually exclusive first composite events. The 62 first events are 57
   acute heart failure events, 0 transplantations, and 5 deaths; there are no
   same-day ties. Across follow-up there are 2 transplantations and 6 deaths.
3. Source semantics were audited. The extract supports a dated acute heart
   failure event but not hospitalization wording.
4. Surgical myectomy remains separate: 33 occurred. Generic ablation was
   excluded because alcohol septal ablation cannot be distinguished reliably
   from arrhythmia ablation.
5. AF, VT/VF, sudden death, ICD placement, and pacemaker placement were not
   added to the primary composite.
6. The secondary V̇E/V̇CO2 exposure was refit with an exposure-by-log-time
   coefficient; its effect increases over follow-up and cannot be summarized by
   one time-constant HR.
7. Reproducible outputs are:
   - `2_Output/Section_05_Outcomes/Table_Endpoint_Definition_Audit.csv`;
   - `2_Output/Section_05_Outcomes/Table_Figure5_Event_Components.csv` and `.docx`;
   - `2_Output/TableS_VEVCO2_TimeVarying_Cox.csv` and `.docx`;
   - `2_Output/TableS_VEVCO2_TimeVarying_Cox_Coefficients.csv`;
   - `2_Output/Section_05_Outcomes/Table_Figure5_Myectomy_Cox.csv`; and
   - `2_Output/Figure5_HeartFailure_Outcomes.pdf`.

## 5. Missing-indicator method: exact use and limitations

### What the method does

For an incomplete adjustment covariate `C`:

1. Create `M_C = 1` when `C` is missing and 0 otherwise.
2. Replace missing `C` values with a fixed constant.
3. Include both the filled covariate and its missingness indicator in the model.

For a continuous covariate, standardize from observed values and fill missing
standardized values with 0:

```r
mu_c <- mean(dat$C, na.rm = TRUE)
sd_c <- sd(dat$C, na.rm = TRUE)

dat <- dat |>
  dplyr::mutate(
    C_missing = as.integer(is.na(C)),
    C_z = (C - mu_c) / sd_c,
    C_z_filled = dplyr::coalesce(C_z, 0)
  )

fit <- lm(outcome ~ exposure + C_z_filled + C_missing, data = dat)
```

For a binary adjustment covariate:

```r
dat <- dat |>
  dplyr::mutate(
    beta_blocker_missing = as.integer(is.na(bb_any)),
    beta_blocker_filled = dplyr::coalesce(bb_any, 0)
  )

fit <- lm(
  outcome ~ exposure + beta_blocker_filled + beta_blocker_missing,
  data = dat
)
```

For several incomplete covariates, repeat this pair for each variable:

```r
fit <- lm(
  peak_vo2 ~ e_e_z + age + Sex + BMI +
    septal_thickness_z_filled + septal_thickness_missing +
    lvot_gradient_z_filled + lvot_gradient_missing +
    beta_blocker_filled + beta_blocker_missing +
    ndhp_ccb_filled + ndhp_ccb_missing +
    hypertension_filled + hypertension_missing,
  data = cross_sectional_df
)
```

The same structure can be used syntactically in mixed and Cox models:

```r
nlme::lme(
  peak_vo2 ~ time * e_e_z + age + Sex + BMI +
    lvot_gradient_z_filled + lvot_gradient_missing,
  random = ~ 1 | ID,
  data = longitudinal_df
)

survival::coxph(
  survival::Surv(follow_up_yrs, event) ~ e_e_z + age + Sex + BMI +
    lvot_gradient_z_filled + lvot_gradient_missing,
  data = outcomes_df
)
```

### What “consistent sample size” really means

The method preserves a row only when the **outcome and primary exposure are
observed**. It can preserve N when an adjustment covariate is missing. It cannot
recover a missing outcome, missing event date, or missing physiologic exposure.

Applied to this study:

| Variable role | Missing-indicator use | Recommendation |
|---|---|---|
| E/e' or LAVi defining cohort entry | No | Require both for the primary cohort. |
| E/e' or LAVi as primary exposure | No | A missingness flag does not reconstruct the exposure. |
| TRVmax as primary exposure | Not for primary inference | Keep available-case N=145 cross-sectional and N=144 outcomes. |
| Age, sex, BMI | Not needed | Complete in the primary model samples. |
| LVOT gradient, wall thickness, medications, hypertension | Possible sensitivity only | Use paired filled-value + missing-indicator terms if a fixed-N sensitivity is desired. |
| Longitudinal CPET outcome at an unobserved visit | No | Mixed models use observed outcomes under their missing-at-random assumption; do not fill an absent outcome. |
| Event date or endpoint status | No | Resolve by source validation/adjudication, not a missing indicator. |

### Why not use it to make the TRVmax model look like N=261?

One could set missing standardized TRVmax to 0 and add `TRVmax_missing`. The
software would then report N=261, but the TRVmax slope would still be informed
only by the 145 measured values. The remaining 116 patients contribute only to
the missingness-category contrast and to other coefficients. This does not
create TRVmax information, and in an observational cohort it can bias the
exposure coefficient when measurement is related to disease severity or the
outcome. It is particularly unsuitable for a spline model because the filled
constant creates an artificial concentration at the fill value.

Therefore:

- Do not use the missing-indicator method to claim a 261-patient TRVmax effect.
- If run at all, label it a sensitivity analysis: “measured-value TRVmax slope
  with adjustment for TRVmax missingness.”
- Keep the available-case TRVmax analysis as the transparent primary TRVmax
  result.

### Statistical caution and preferable alternative

This is an observational study. Missing-indicator adjustment for incomplete
confounders can produce biased estimates, even under seemingly benign
missingness patterns. It is most defensible here as a **sensitivity analysis for
secondary adjustment covariates**, not as the primary missing-data method.

If preserving N for expanded adjustment becomes important, multiple imputation
is the preferred sensitivity:

1. Impute only variables whose missingness is plausibly missing at random given
   observed data.
2. Include the exposure, outcome, all analysis covariates, and auxiliary
   predictors of missingness in the imputation model.
3. Preserve nonlinear terms and interactions using a substantive-model-
   compatible method when needed.
4. For Cox models, include the event indicator and Nelson-Aalen cumulative
   hazard in the imputation model.
5. For longitudinal models, use a multilevel imputation procedure that respects
   repeated observations within patient.
6. Use at least as many imputations as the approximate percentage of incomplete
   cases, inspect convergence, and pool estimates with Rubin's rules.
7. Compare complete-case, missing-indicator, and multiple-imputation estimates;
   interpret disagreement as missing-data sensitivity rather than selecting the
   most favorable result.

Every missing-indicator term also consumes a model coefficient. In the Cox
models, this matters: a five-variable expanded model with five additional
missingness indicators adds ten adjustment coefficients before the exposure,
which is not compatible with only 62 events. A fixed N is not useful if achieved
by severe overfitting.

## 6. Figure revision plan

### Figure 1: cohort flow

- Keep CONSORT-style cohort assembly as the main Figure 1.
- Show the E/e' + LAVi requirement as the upstream analytic filter.
- Show N=261 as the cross-sectional parent cohort.
- Show nested longitudinal N=105/294 CPETs and outcomes N=259/62 events.
- Add a separate comparator branch: 1,015 potential references, 243 QC eligible,
  242 with complete matching variables, and the final balanced 177:177 match.
- Label the method as 1:1 exact-sex age/BMI Mahalanobis matching within a
  0.25-SD propensity-score caliper.

### Figure 2: measured parameter patterns

- Preserve the current publication design.
- Panel A: report E/e' 84/261, LAVi 118/261, and TRVmax 17/145 using observed
  denominators.
- Panel B: use the four E/e' + LAVi combinations only.
- Panels C and D: retain the same four combinations and the bottom combination
  matrix.
- Use “Neither parameter above threshold,” not “normal diastolic function.”
- Rename “Average e' (age-adjusted)” to “Average e' below age-specific reference
  limit” or “Age-calibrated average e'.”
- Keep all legends flush left and fully justified.

### Figure 3: cross-sectional associations

- Keep E/e' and LAVi panels at N=261.
- Label TRVmax panels N=145 and “available case.”
- Reduce each statistics box to N, overall P, and nonlinearity Q.
- Retain identical typography, axes, colors, spacing, and panel dimensions.
- Put detailed breakpoint diagnostics in the supplement.

### Figure 4: longitudinal analysis

- Use the compact linear mixed-model table as the preferred main-paper
  presentation; all rows report patient and observation counts.
- Retain the six-panel GAMM trajectory figure unchanged in the supplement. Its
  six interaction tests are nonsignificant.
- Include the four-panel raw observed-trajectory figure in the supplement to
  make data construction transparent. It shows lower and upper baseline halves
  with faint within-patient lines and points, but no fitted model.
- If the journal requires Figure 4 to remain a main figure, use the existing
  design without restyling and describe it as a flexible descriptive
  sensitivity, not the primary inferential result.
- Do not visually or textually imply a significant longitudinal effect.

### Figure 5: outcomes

- Keep the existing design and VE/VCO2 cutoff at 30.
- Update all panel Ns and event counts to the E/e' + LAVi outcomes cohort.
- Include clinical event components in the legend and main text, not as a
  separate figure panel.
- Keep surgical myectomy separate and explicitly labeled as such.
- Use “dated acute heart failure event”; the supplied source does not support
  “hospitalization.”
- Interpret the V̇E/V̇CO2 Cox result using the supplemental time-varying
  coefficient estimates, not the constant HR in the forest panel.

## 7. Table revision plan

### Table 1

- Reorder morphology as asymmetric septal, symmetric, and apical.
- Use “Non-obstructive HCM,” defined by measured resting LVOT gradient <30 mm Hg.
- Use “Obstructive HCM,” defined by measured resting LVOT gradient >=30 mm Hg.
- State that provocative/exercise gradients were unavailable.
- Replace the table only after matching balance is corrected.

### Main results table

Create one compact table with three clearly separated sections:

1. Cross-sectional RCS overall associations.
2. Longitudinal linear mixed-model parameter-by-time interactions.
3. Cox models for the acute HF/transplant/death composite.

Every row must report patient N; longitudinal rows must also report observation
N, and Cox rows must report event N. Mark TRVmax rows as available case.

### Supplemental tables

- Missingness and selection audit.
- Same-sample crude versus age/sex/BMI-adjusted associations.
- Expanded complete-case, missing-indicator, and—if implemented—multiple-
  imputation sensitivity models.
- Morphology sensitivity.
- Full RCS/nonlinearity diagnostics.
- Full mixed-model and GAMM statistics.
- Cox proportional-hazards diagnostics.
- Event component counts and surgical-myectomy analysis.

## 8. Manual manuscript revision plan

### Abstract

- State N=261 cross-sectional, N=105 longitudinal, and N=259 outcomes.
- Define the primary endpoint explicitly.
- State that E/e' and LAVi were required and TRVmax was available case.
- Avoid `LVDD`, “normal diastolic function,” and composite diastolic grades.
- Do not claim a significant longitudinal trajectory association.

### Introduction

End with three separate questions:

1. Are resting LV diastolic parameters associated with contemporaneous CPET
   performance?
2. Do baseline parameters predict longitudinal CPET trajectories?
3. Are baseline parameters associated with the date-validated acute
   HF/transplant/death endpoint?

### Methods

- Define the HCM entry criteria and confirm that HCM-labeled patients failing
  those criteria were excluded, not used as controls.
- Define references as phenotypically screened non-HCM CPET comparators, not
  healthy volunteers, and state that all LVOT gradients are resting. For the
  reviewer response, report that no matched reference had an HCM flag or
  resting gradient >=30 mm Hg (maximum 10.29 mm Hg; eligible-pool maximum
  10.67 mm Hg).
- Define the E/e' + LAVi cohort and available-case TRVmax approach.
- Explain why age, sex, and BMI form the minimal adjustment set.
- State that parameter models are separate to avoid collinearity.
- Explain overall P and nonlinearity Q in plain language.
- Define the linear mixed model and the outcome-level random-effects decision:
  random intercept for peak V̇O2; independent random intercept and slope for
  V̇E/V̇CO2. State that the choice required support from both primary E/e'
  and LAVi models and that likelihood-ratio tests were descriptive because
  variance components lie on a boundary.
- Define all endpoint components, censoring, and surgical myectomy. Use “dated
  acute heart failure event”; use “hospitalization” only if independently
  verified from a data dictionary or chart adjudication.
- Explain the V̇E/V̇CO2 exposure-by-log-time Cox sensitivity and report
  time-specific HRs at 1, 3, and 5 years.
- State the missing-data strategy and label missing-indicator/MI analyses as
  sensitivities if used.

### Results

- Begin with cohort flow and observed denominators.
- Report Figure 2 as measured threshold patterns, not a global diagnosis of
  diastolic dysfunction.
- Report cross-sectional E/e', LAVi, and available-case TRVmax separately.
- State that all six longitudinal interactions were not statistically
  significant under the selected random-effects structures. Report the primary
  E/e' and LAVi estimates from the compact table and label TRVmax available
  case; do not report the former LAVi P=0.0503 as the primary result.
- Report Cox HRs for E/e' and LAVi from the shared N=259/62-event cohort and
  TRVmax from N=144/36 events.
- Report every component of the clinical composite, distinguishing any
  occurrence from the first event that determined composite time.
- Report the time-varying V̇E/V̇CO2 association rather than treating its
  constant Cox HR as time invariant.
- Report morphology and missing-data analyses as sensitivities.

### Discussion organization

1. Main findings.
2. Prior literature: peak V̇O2, ventilatory efficiency, pulmonary/left-atrial
   loading, and limitations of resting Doppler parameters in HCM.
3. Original contribution: consistent association with ventilatory efficiency,
   clinical prognosis for E/e' and LAVi, and absence of confirmed longitudinal
   trajectory effects.
4. Limitations.
5. Short conclusion.

### Required limitations

- Retrospective, single-center design.
- Complete E/e' + LAVi selection may limit generalizability.
- TRVmax is missing in 116/261 and is analyzed available case.
- Missing-indicator analyses do not restore missing physiologic information and
  may be biased in observational data.
- No provoked or exercise LVOT gradients.
- Expanded adjustment can overfit and is sensitive to covariate completeness.
- Morphology is unknown in 36/261.
- The extract supports a dated acute-HF event but not hospitalization semantics;
  hospitalization wording requires external source validation.
- AF, VT/VF, sudden death, and alcohol septal ablation are not analysis-ready.
- OMARX is exploratory and does not establish a new threshold.

## 9. Decision gates before submission

1. **Comparator gate (passed):** all post-match SMDs <0.10.
2. **Missingness gate:** included-versus-excluded and TRVmax-measured-versus-
   missing comparisons complete.
3. **Endpoint gate (passed with conservative wording):** the extract supports a
   dated acute heart failure event. Hospitalization is not established and must
   not be claimed without external documentation.
4. **Model gate (longitudinal and outcomes portions passed):** outcome-level
   random-effects structures, Cox adjustment sets, and the primary/sensitivity
   hierarchy are frozen; no parameter-specific structure, selected subset, or
   invalid time-constant effect is promoted because it gives a preferred result.
5. **Figure gate (passed):** latest PDFs visually inspected for clipping,
   overlap, inconsistent padding, glyph problems, and legend justification.
6. **Manuscript gate:** every N, event count, cutoff, and endpoint definition
   reconciled manually without modifying Paperpile fields.

## 10. Final reconciliation checklist

- [x] Refit comparator matching and verify all SMDs <0.10.
- [x] Generate missingness and selection audit.
- [x] Multiple imputation is not required for the primary revision because the
      primary models use the fixed minimal adjustment set without covariate-
      driven denominator loss. Reserve it for a reviewer-requested expanded-
      adjustment sensitivity; do not use missing indicators for primary
      exposures.
- [x] Complete the V̇E/V̇CO2 proportional-hazards sensitivity.
- [x] Audit acute-HF source semantics and restrict wording to “dated acute heart
      failure event.”
- [x] Complete the longitudinal random-effects audit and freeze the
      outcome-level structures.
- [x] Generate the compact primary longitudinal table and supplemental raw
      observed-trajectory figure.
- [x] Regenerate the analysis Table 1 from the final matched population; update
      the submission manuscript manually to preserve Paperpile fields.
- [x] Reconcile 261 cross-sectional, 105 longitudinal, and 259 outcomes counts
      in the analysis outputs.
- [x] Mark every TRVmax analysis output with its observed N.
- [ ] Confirm no missing measurement is coded/described as normal.
- [ ] Replace `LVDD` with “LV diastolic parameters” where possible.
- [x] Confirm surgical myectomy is not mislabeled as all SRT.
- [x] Confirm AF/VT/VF/sudden death are not added without adjudication.
- [x] Verify generated PDF legends are flush left and fully justified; recheck
      placement after manual insertion into Word.
- [ ] Preserve all Paperpile fields by editing the Word manuscript manually.

## 11. Key methods references

- Groenwold RHH, et al. Missing covariate data in clinical research: when and
  when not to use the missing-indicator method for analysis. *CMAJ*. 2012.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC3414599/
- Greenland S, Finkle WD. A critical look at methods for handling missing
  covariates in epidemiologic regression analyses. *Am J Epidemiol*. 1995.
  https://pubmed.ncbi.nlm.nih.gov/7503045/
- Sterne JAC, et al. Multiple imputation for missing data in epidemiological and
  clinical research: potential and pitfalls. *BMJ*. 2009.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC2714692/
- White IR, Royston P. Imputing missing covariate values for the Cox model.
  *Stat Med*. 2009. https://pmc.ncbi.nlm.nih.gov/articles/PMC2998703/
- Lee KJ, et al. Framework for the treatment and reporting of missing data in
  observational studies (TARMOS). *J Clin Epidemiol*. 2021.
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8168830/

## Current source locations

- Analysis source: `_Scripts/README_HCM_Manuscript_V2.qmd`
- Figure-export helpers: `_Scripts/R/utils.R`
- Current main figures: `2_Output/Manuscript_Main/`
- Cross-sectional outputs: `2_Output/Section_03_CrossSectional/`
- Longitudinal outputs: `2_Output/Section_04_Longitudinal/`
- Outcomes outputs: `2_Output/Section_05_Outcomes/`
- QC outputs: `2_Output/Section_06_Ancillary/`
