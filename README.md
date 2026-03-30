# Modeling Left Ventricular Diastolic Dysfunction and Cardiopulmonary
Fitness in HCM
Mark E. Pepin, MD, PhD, MS, FESC
2026-03-25

# Data Import and Cohort Assembly

## Import and clean CPX data

The initial phase imports raw Cardiopulmonary Exercise Testing (CPET)
data, assigns a stable patient identifier (ID) per MRN, restricts to HCM
with adequate effort (peak RER $\geq 1.0$), and derives biometric
indices (BSA, BMI, LBMI). Baseline and cross-sectional analyses use all
eligible CPET studies; longitudinal analyses later restrict to patients
with repeated aligned testing.

    CPX tests (HCM, adequate effort): 1649 tests, 804 patients

## Comorbidity exclusions

Patients with documented CAD, COPD, or interstitial lung disease are
excluded so that exercise limitations are attributable to HCM-specific
pathology.

## Merge echocardiographic data

Each CPET is aligned with the closest resting echocardiogram within 1
week on either side of the CPX date.

    Baseline cohort: 643 patients
    Longitudinal repeated-test cohort: 1011 observations across 317 patients

## Import clinical outcomes

    Patients with outcome data: 643
    Composite HF events: 170

## Cohort flow diagram

![](README_HCM.Diastologyv2_files/figure-commonmark/flow_diagram-1.png)

**Legend.** Cohort assembly workflow for the analytic dataset.
Horizontal bars show the number of CPX studies retained at each
sequential filtering step, with labels indicating both the number of
tests and the number of unique patients remaining after each
restriction. The final bar represents the CPX studies successfully
aligned to a resting echocardiogram within the prespecified time window.

# ASE 2025 HCM Diastolic Assessment (Primary Parameter Combination Phenotypes)

To avoid collapsing HCM diastolic physiology into a binary label, we
categorized baseline studies according to the observed combination of
the 3 HCM-relevant primary parameters available in this dataset: average
E/e$'$ \>14, left atrial volume index (LAVI) \>34 mL/m$^2$, and peak
tricuspid regurgitation velocity (TR V$_{max}$) \>2.8 m/s. These
variables were selected because they correspond to the principal
HCM-related markers highlighted in the 2025 American Society of
Echocardiography update for special populations.

Combination phenotypes were assigned whenever at least 1 of the 3
primary parameters was available at baseline. Each patient was mapped to
1 of 8 mutually exclusive profiles reflecting the observed pattern of
abnormal findings: none abnormal, isolated E/e$'$ abnormality, isolated
LAVI abnormality, isolated TR V$_{max}$ abnormality, each possible
2-parameter combination, or all 3 parameters abnormal. Thus, patients
with incomplete primary-parameter data were retained in these
descriptive analyses if at least 1 observed parameter was available; the
“none abnormal” category denotes that no observed primary parameter was
abnormal, not that all 3 parameters were measured and normal.

**Primary HCM Parameter Availability:**

3 indicators available: 230 (35.8%) 2 indicators available: 243 (37.8%)
0-1 indicators available: 170 (26.4%)

**Primary Parameter Combination Phenotypes (Patients With \>=1 Primary
Parameter):**

None abnormal E/e’ only LAVI only TR Vmax only E/e’ + LAVI 327 54 96 18
60 E/e’ + TR Vmax LAVI + TR Vmax All 3 abnormal 5 14 9

**Individual Criterion Prevalence:**

E/e’ \>14: 128 / 465 (27.5%) LAVI \>34: 179 / 447 (40.0%) TR Vmax \>2.8:
46 / 374 (12.3%)

![](README_HCM.Diastologyv2_files/figure-commonmark/ase_figure-1.png)

**Figure 2. Baseline HCM Primary-Parameter Phenotypes and
Cardiopulmonary Performance.** **(A)** Frequency of each observed
combination of abnormal E/e$'$, LAVI, and TR V$_{max}$, with black dots
indicating parameters contributing to a given phenotype and gray dots
indicating parameters not present in that combination. **(B)** Peak
VO$_2$ (FRIEND 2.0 % predicted) across phenotypes. **(C)** VE/VCO$_2$
slope across the same phenotypes. Boxplots display the median and
interquartile range, and jittered points represent individual patients.
Bar and strip colors denote the number of abnormal primary parameters
within each phenotype. Abbreviations: HCM = hypertrophic cardiomyopathy;
LAVI = left atrial volume index; TR V$_{max}$ = peak tricuspid
regurgitation velocity; VE/VCO$_2$ = ventilatory efficiency slope;
VO$_2$ = oxygen consumption.

# Cross-Sectional Nonlinear Analysis

Cross-sectional nonlinear modeling is focused on the 3 HCM-relevant ASE
2025 primary parameters available in this dataset: average E/e$'$, LAVI,
and TR V$_{max}$. For each parameter, we fit restricted cubic spline
models adjusted for age, sex, and lean body mass index, then display
predicted performance across the observed range with the ASE 2025
threshold overlaid. We next repeat the same clinical framing using OMAR,
allowing data-driven thresholds to be compared directly with the ASE
reference cutoffs.

## Primary outcome: Peak VO$_2$ (FRIEND 2.0 % predicted)

![](README_HCM.Diastologyv2_files/figure-commonmark/cross_sectional_lrt_table-1.png)

## Secondary outcome: VE/VCO$_2$ slope

![](README_HCM.Diastologyv2_files/figure-commonmark/cross_sectional_vevco2-1.png)

**Figure 3. Focused Restricted Cubic Spline Partial-Effect Plots for
Predictors With the Strongest Evidence of Nonlinearity.** Panels show
adjusted restricted cubic spline partial-effect curves for the ASE
primary diastolic parameter(s) prioritized by the prespecified
nonlinearity screen, with Peak VO$_2$ (FRIEND 2.0 % predicted) in the
top row and VE/VCO$_2$ slope in the bottom row. Raw observed data are
shown as faint background points, shaded bands denote 95% confidence
intervals, and dashed vertical lines mark ASE 2025 thresholds. In-panel
boxes report the complete-case sample size, overall association p value,
nonlinearity q value, and spline-vs-linear 394AIC. The full six-panel
primary-parameter RCS display is exported separately to the supplemental
output. Abbreviations: ASE = American Society of Echocardiography; CI =
confidence interval; HCM = hypertrophic cardiomyopathy; LAVI = left
atrial volume index; TR V$_{max}$ = peak tricuspid regurgitation
velocity; VE/VCO$_2$ = ventilatory efficiency slope; VO$_2$ = oxygen
consumption.

![](README_HCM.Diastologyv2_files/figure-commonmark/figure3_summary_table-1.png)

## Sensitivity: Echo alignment window

# OMARX: Adaptive Threshold Discovery

OMARX, implemented from the newer local `omarx_v3.ipynb` notebook
engine, was used to compare a MARS-style global discovery strategy
against a confirmatory additive fit while keeping the manuscript
outcomes unchanged. The revised OMARX analysis is intended to clarify
whether Peak VO$_2$ is driven primarily by age/body-size structure and
whether VE/VCO$_2$ slope shows stronger cardiac contribution from LAVI,
TR V$_{max}$, medial septal e$'$, septal thickness, and obstruction.

Primary OMARX models use the same baseline HCM cohort but expand the
predictor set to include cardiac and body-composition variables
highlighted in the new notebook and the feedback note. Complete-case
MARS models are treated as the primary discovery analysis, additive fits
are confirmatory, and targeted QC is added for clinical-knot protection,
TRV informative missingness, BMI-vs-LBMI choice, binary obstruction
coding, predictor imputation, and bootstrap stability of the LAVI knot
in VE/VCO$_2$.

<img
src="README_HCM.Diastologyv2_files/figure-commonmark/omar_publication_figure-1.png"
style="width:100.0%" />

**Legend.** Publication-format summary of the revised OMARX analysis.
**(A)** Variable-importance profile for the BMI-adjusted OMARX MARS
model of Peak VO$_2$, showing dominant contributions from age and BMI
and no retained primary diastolic knot. **(B)** Variable-importance
profile for the corresponding VE/VCO$_2$ model, highlighting LAVI as the
dominant predictor with supplementary contributions from age, sex,
medial septal e$'$, and septal thickness. **(C)** OMARX partial-effect
plot for LAVI in the VE/VCO$_2$ model, with raw complete-case
observations in the background, the primary retained OMARX knot shown as
a solid red line, and the ASE reference value of 34 mL/m$^2$ shown as a
dashed blue line. **(D)** Robustness of the VE/VCO$_2$ LAVI threshold
across prespecified sensitivity analyses and bootstrap quality control.
Segment ranges denote the bootstrap interquartile range where
applicable.

<img
src="README_HCM.Diastologyv2_files/figure-commonmark/omar_figure-1.png"
style="width:100.0%" />

**Legend.** OMARX partial-effect plots for the retained and clinically
relevant heart predictors from the primary complete-case BMI-adjusted
models. Navy curves show the primary MARS discovery fit and gray curves
show the confirmatory additive fit. Solid red vertical lines denote
retained OMARX knots from the MARS fit, and dashed blue lines denote ASE
reference thresholds where available.

<img
src="README_HCM.Diastologyv2_files/figure-commonmark/omar_rcs_comparison-1.png"
style="width:100.0%" />

**Legend.** Direct comparison of the revised focused RCS curves and the
primary OMARX MARS partial-effect curves for the predictors that survive
the linearity screen and OMARX QC prioritization. Dashed blue vertical
lines indicate ASE cutoffs and solid red vertical lines indicate
retained OMARX knots.

## OMARX Results

For Peak VO$_2$, the primary complete-case BMI-adjusted MARS model
emphasized Age, BMI, Medial septal e’, consistent with the
interpretation that age and body-size structure dominate this outcome
more strongly than diastolic threshold structure. For VE/VCO$_2$ slope,
the primary complete-case BMI-adjusted MARS model prioritized LAVI, Age,
Female sex, Medial septal e’, supporting a stronger cardiac/diastology
signal than was seen for Peak VO$_2$. Retained MARS knots were
identified in LAVI at 33.8. In bootstrap QC, the LAVI knot was retained
in 59.5% of resamples, with median knot 34.8 mL/m$^2$ and 100.0% of
retained knots lying within 5 mL/m$^2$ of 34.

# Longitudinal GAMM Analysis

We employ generalized additive mixed models (GAMMs) via `mgcv::bam()` to
model nonlinear longitudinal trajectories of peak VO$_2$ and VE/VCO$_2$
slope over time, testing whether the trajectory diverges as a function
of baseline diastolic severity via a tensor product interaction.

    GAMM cohort: 773 observations, 293 patients

    VO2_FRIEND2_PP GAMM cohort: 773 observations, 293 patients
    VeVco2_slope GAMM cohort: 772 observations, 293 patients

![](README_HCM.Diastologyv2_files/figure-commonmark/gamm_analysis-1.png)

**Figure 4. Longitudinal Model-Derived Cardiopulmonary Trajectories
According to Baseline Diastolic Dysfunction Burden.** **(A)**
Model-derived peak VO$_2$ trajectories from a generalized additive mixed
model across quartiles of baseline diastolic dysfunction burden, defined
using the baseline LAVI-derived DD Z-score. **(B)** Model-derived
VE/VCO$_2$ slope trajectories across the same quartiles. Shaded bands
denote 95% confidence intervals. All models were adjusted for age, sex,
and lean body mass index. Abbreviations: CPET = cardiopulmonary exercise
testing; DD = diastolic dysfunction; GAMM = generalized additive mixed
model; LAVI = left atrial volume index; VE/VCO$_2$ = ventilatory
efficiency slope; VO$_2$ = oxygen consumption.

![](README_HCM.Diastologyv2_files/figure-commonmark/figure4_summary_table-1.png)

# Heart Failure Outcome Analysis

Time-to-event analysis tests whether ASE 2025 filling pressure
classification (Normal vs Elevated), continuous DD Z-scores, LVOT
gradient, and LV septal thickness predict a composite heart failure
endpoint (acute HF, chronic HF, transplant, or death).

![](README_HCM.Diastologyv2_files/figure-commonmark/hf_outcomes-1.png)

**Figure 5. Heart Failure Outcomes Across ASE Threshold Components,
Primary-Parameter Status, and Structural Obstruction Markers.** **(A)**
Kaplan-Meier estimates for the composite heart failure endpoint
stratified by abnormal baseline E/e$'$ (\>14). **(B)** Kaplan-Meier
estimates stratified by abnormal baseline LAVI (\>34 mL/m$^2$). **(C)**
Kaplan-Meier estimates stratified by abnormal baseline TR V$_{max}$
(\>2.8 m/s). **(D)** Kaplan-Meier estimates stratified by the same
binary primary-parameter grouping used in the patient-characteristics
table (`Normal` = no abnormal primary parameter, `Abnormal` = at least 1
abnormal primary parameter). **(E)** Adjusted hazard ratios from
multivariable Cox proportional hazards models for elevated filling
pressure classification, continuous diastolic dysfunction burden
expressed as the baseline LAVI-derived DD Z-score, max LVOT gradient, LV
septal thickness, age, sex, and lean body mass index. Shaded bands
denote 95% confidence intervals in Kaplan-Meier panels. Abbreviations:
CPET = cardiopulmonary exercise testing; DD = diastolic dysfunction; HCM
= hypertrophic cardiomyopathy; HR = hazard ratio; LAVI = left atrial
volume index; LVOT = left ventricular outflow tract; TR V$_{max}$ = peak
tricuspid regurgitation velocity.

# Supplemental Analyses

## All diastolic parameters: restricted cubic splines

<img
src="README_HCM.Diastologyv2_files/figure-commonmark/supplemental_rcs_all_figure-1.png"
style="width:100.0%" />

**Legend.** Supplemental restricted cubic spline models relating all
available diastolic parameters to Peak VO$_2$ and VE/VCO$_2$ slope. Each
panel shows the adjusted spline-estimated association for a single
parameter-outcome pair after adjustment for age, sex, and lean body mass
index. Shaded ribbons indicate 95% confidence intervals. Dashed blue
vertical lines indicate ASE 2025 reference cutoffs where available for
the 3 HCM primary parameters.

## All diastolic parameters: OMAR

<img
src="README_HCM.Diastologyv2_files/figure-commonmark/supplemental_omar_all_figure-1.png"
style="width:100.0%" />

**Legend.** Supplemental OMARX-derived partial dependence plots relating
all available diastolic parameters to Peak VO$_2$ and VE/VCO$_2$ slope.
Each panel shows the modeled marginal association for a single
parameter-outcome pair after adjustment for age, sex, and lean body mass
index. Solid red vertical lines denote OMARX-retained hinge locations,
and dashed blue vertical lines denote ASE 2025 reference cutoffs where
available for the 3 HCM primary parameters.

## Age-residualized diastolic indices

![](README_HCM.Diastologyv2_files/figure-commonmark/residualize-1.png)

![](README_HCM.Diastologyv2_files/figure-commonmark/residualize-3.png)

**Legend.** Distribution of age-residualized diastolic indices.
Residuals were derived from nonlinear age-adjustment models for each
index so that the resulting values reflect deviation from the expected
age-specific value. Boxplots show the median, interquartile range, and
distribution of residual variability across the included
echocardiographic parameters.

## Export missing echo audit

    Exported 164 CPX tests without echo linkage (88 patients)
