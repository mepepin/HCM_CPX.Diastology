# HCM CPX Diastology

This repository contains the Quarto analysis source for modeling left ventricular diastolic dysfunction and cardiopulmonary fitness in hypertrophic cardiomyopathy (HCM), together with a GitHub-friendly summary of the analytic workflow and major outputs.

## Analysis Summary

The analysis integrates serial cardiopulmonary exercise testing (CPX), resting echocardiography, and clinical outcomes to quantify how diastolic dysfunction relates to exercise capacity, longitudinal functional decline, and heart failure events in HCM.

### Core workflow

1. Import and clean serial CPX data, restrict to HCM with adequate effort, and derive biometric covariates.
2. Exclude major competing pulmonary and ischemic comorbidities.
3. Align each CPX study to the nearest prior echocardiogram within a prespecified window.
4. Grade diastolic dysfunction using an ASE 2025-style hierarchical framework based on available echocardiographic criteria.
5. Model cross-sectional associations between diastolic indices and peak VO2 or VE/VCO2 slope.
6. Model longitudinal peak VO2 trajectories with mixed-effects and GAMM approaches.
7. Evaluate time-to-event relationships for a composite heart failure endpoint.

### Main components

- `ASE 2025 diastolic dysfunction grading`
  Baseline studies are classified as Normal, Grade I, Grade II, Grade III, or Indeterminate using available E/e', left atrial volume index, and tricuspid regurgitation velocity data.
- `Continuous modeling`
  GAMLSS and SEAMLSS are used to derive continuous diastolic severity measures and Z-scores relative to predicted peak VO2.
- `Cross-sectional analysis`
  Natural spline models test whether individual diastolic indices improve model fit beyond age, sex, and lean body mass index.
- `Longitudinal analysis`
  Generalized additive mixed models assess whether baseline diastolic severity modifies the trajectory of peak VO2 over time.
- `Outcome analysis`
  Kaplan-Meier and Cox models test whether diastolic severity predicts a composite heart failure endpoint.

### Selected signals from the rendered analysis

- Diastolic dysfunction severity showed separation in functional capacity and ventilatory efficiency across ASE 2025 grades.
- In cross-sectional spline modeling for peak VO2, the strongest incremental signals came from mitral A duration, left atrial volume index, mitral deceleration time, and tricuspid regurgitation velocity.
- The longitudinal GAMM favored a nonlinear time-by-diastolic severity interaction over a simpler linear mixed-effects model.
- Time-to-event analyses showed separation in event-free survival by diastolic dysfunction grade, and continuous diastolic Z-scores were also associated with outcomes.

## Repository Contents

- [README_HCM.Diastology.qmd](./README_HCM.Diastology.qmd)
  Executable Quarto source for the full analysis.
- `README.md`
  This summary page for GitHub.

## Notes

- The Quarto document is preserved from the original analysis workspace and uses project-relative paths such as `../1_Input` and `../2_Output`.
- The generated figures and tables referenced in the Quarto file were produced in the source project and are not duplicated here unless added separately.
- The rendered analysis date in the Quarto source is March 25, 2026.
