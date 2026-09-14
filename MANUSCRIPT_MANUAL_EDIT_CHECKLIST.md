# Paperpile-Preserving Manual Manuscript Edit Checklist

Status: ready for manual editing after completion of analysis and figure Stage 6B.1.

This checklist is the authoritative order of operations for revising the
submission manuscript manually. Do not regenerate the Word document. Do not use
global replace across the reference list, citations, captions, or Paperpile
fields. Make prose changes section by section in a duplicate of the current
submission file and verify each numerical statement against
`2_Output/HCM_Revision_Results_Ledger.xlsx`.

## 1. Create a safe manual working copy

- [ ] In Word, use **Save As** to create a dated revision copy of the current
  manuscript.
- [ ] Confirm that in-text Paperpile citations still behave as citation fields
  in the working copy.
- [ ] Do not update all fields, unlink fields, paste the whole document into a
  new file, or rebuild the bibliography.
- [ ] Use tracked changes if the coauthor workflow requires it.
- [ ] Keep this checklist and the results ledger open beside Word.
## 2. Global terminology audit

Use Word Find to review each occurrence individually; do not use Replace All.

- [ ] Replace `LVDD` with `LV diastolic parameters`, `LV diastolic indices`, or
  the specific parameter name, depending on context.
- [ ] Replace claims of `normal diastolic function` with `neither E/e' nor LAVi
  above threshold` when referring to the Figure 2 group.
- [ ] Replace `controls` or `healthy controls` with `non-HCM references` or
  `phenotypically screened non-HCM comparators`.
- [ ] Use `non-obstructive HCM`, not `non-obstructive` alone, in formal labels.
- [ ] Qualify every analyzed LVOT gradient as `resting` unless discussing the
  absence of provoked/exercise gradients.
- [ ] Use `dated acute heart failure event`, not `heart failure
  hospitalization`, because hospitalization was not established by the source
  extract.
- [ ] Use `surgical myectomy`, not `septal reduction therapy`, for the validated
  secondary procedure endpoint. Alcohol septal ablation was not reliably
  distinguishable from arrhythmia ablation.
- [ ] Standardize physiological typography to peak V-dot O2 and V-dot E/V-dot
  CO2 slope, with the 2 as a true subscript.
- [ ] Use `age-calibrated average e'` or `average e' below the age-specific
  reference limit`; avoid the ambiguous label `age-adjusted e'`.

## 3. Numbers that must remain consistent everywhere

| Population or result | Locked value |
|---|---:|
| Phenotypic HCM entry cohort | 440 patients |
| Cross-sectional E/e' + LAVi cohort | 261 patients |
| TRVmax measured within cross-sectional cohort | 145/261 |
| Longitudinal parent subcohort | 105 patients |
| Peak-V̇O2 longitudinal observations | 294 CPETs |
| V̇E/V̇CO2 longitudinal observations | 293 CPETs |
| TRVmax longitudinal subset | 63 patients; 176 peak-V̇O2 and 175 V̇E/V̇CO2 observations |
| Clinical-outcomes cohort | 259 patients; 62 primary events |
| TRVmax outcomes subset | 144 patients; 36 events |
| Matched comparison | 177 HCM and 177 non-HCM references |
| Primary first events | 57 acute HF, 0 transplant, 5 death |
| Surgical myectomies | 33 occurrences, separate endpoint |

- [ ] Search the manuscript for each legacy cohort N and reconcile it with this
  table.
- [ ] Report TRVmax with its observed denominator every time it appears in a
  result, table, or caption.
- [ ] Never imply that missing TRVmax was below threshold or normal.

## 4. Title and short title

- [ ] Ensure the title describes LV diastolic parameters rather than a global
  diagnosis of LV diastolic dysfunction.
- [ ] Keep the title aligned with the three study questions: cross-sectional
  CPET performance, longitudinal CPET trajectories, and clinical outcomes.
- [ ] Avoid language implying that this study established a new nonlinear
  threshold.

## 5. Abstract

### Methods sentence

- [ ] State that the cross-sectional cohort required measured E/e' and LAVi.
- [ ] State that TRVmax was analyzed as available case.
- [ ] Give the three nested analytic populations: cross-sectional N=261,
  longitudinal N=105, and outcomes N=259.
- [ ] Define the primary endpoint as the first dated acute HF event, heart
  transplantation, or all-cause death.
- [ ] State that primary models adjusted for age, sex, and BMI and evaluated
  each diastolic parameter separately.

### Results sentence

- [ ] Cross-sectional: emphasize associations with V̇E/V̇CO2 slope for E/e'
  (overall P=0.0039), LAVi (P=0.0158), and available-case TRVmax (P<0.001).
- [ ] State that the corresponding peak-V̇O2 associations were not significant
  (overall P=0.256, 0.530, and 0.335).
- [ ] State that none of the six nonlinearity tests remained significant after
  multiplicity correction; do not promote a breakpoint.
- [ ] Longitudinal: state that no parameter-by-time interaction was significant
  under the selected random-effects structures.
- [ ] Outcomes: report E/e' HR 1.29 per SD (95% CI 1.09-1.53; P=0.003) and LAVi
  HR 1.41 per SD (95% CI 1.14-1.73; P=0.001).
- [ ] If TRVmax is included in the abstract, label it available case and report
  N=144/36 events, HR 1.25 (95% CI 0.97-1.61; P=0.084).
- [ ] Do not present the time-constant V̇E/V̇CO2 HR as the prognostic result;
  its proportional-hazards assumption failed.

### Conclusion sentence

- [ ] Separate the supported findings: resting LV diastolic parameters were
  associated with ventilatory efficiency, E/e' and LAVi were associated with
  the clinical composite, and longitudinal trajectory effects were not
  confirmed.

## 6. Introduction

- [ ] Preserve the literature-rich rationale, but distinguish resting Doppler
  parameters, exercise physiology, and outcomes rather than treating them as a
  single construct.
- [ ] End the Introduction with three explicit objectives:
  1. contemporaneous associations between resting LV diastolic parameters and
     CPET performance;
  2. associations of baseline parameters with longitudinal CPET trajectories;
  3. associations with the dated acute-HF/transplant/death composite.
- [ ] Avoid announcing a goal of deriving a new threshold; the corrected
  nonlinearity tests do not support that claim.

## 7. Methods

### Study population and cohort flow

- [ ] Define the phenotypic HCM entry criteria clearly.
- [ ] State explicitly that HCM-labeled patients not meeting phenotypic entry
  criteria were excluded and were not used as controls.
- [ ] Define non-HCM references as phenotypically screened CPET comparators, not
  healthy volunteers.
- [ ] State that reference-arm gradients were resting measurements; no matched
  reference had a resting LVOT gradient >=30 mm Hg. The matched maximum was
  10.29 mm Hg; the eligible-pool maximum was 10.67 mm Hg.
- [ ] Define obstructive HCM as resting LVOT gradient >=30 mm Hg and
  non-obstructive HCM as resting gradient <30 mm Hg for this analysis.
- [ ] State that provoked/exercise gradients were not uniformly available and
  therefore did not define obstruction in this registry.
- [ ] Describe the upstream E/e' + LAVi filter and the nested longitudinal and
  outcomes subcohorts exactly as shown in Figure 1.

### Echocardiographic variables

- [ ] Define E/e', LAVi, TRVmax, and age-calibrated average e' separately.
- [ ] Describe Figure 2 thresholds as parameter-specific reference thresholds,
  not a diagnosis of normal/abnormal diastolic function.
- [ ] State that TRVmax uses observed denominators and that missing values were
  not classified.

### Statistical analysis

- [ ] Explain that age, sex, and BMI form the prespecified minimal adjustment
  set.
- [ ] Explain that E/e', LAVi, and TRVmax were modeled in separate primary
  models to avoid collinearity and denominator collapse.
- [ ] Move structural measures, resting LVOT gradient, medications,
  hypertension, and morphology to explicitly labeled sensitivity analyses.
- [ ] Explain the RCS metrics in plain language:
  - the overall P value asks whether the diastolic parameter improves model fit
    beyond age, sex, and BMI;
  - the nonlinearity q value asks whether a curved relationship fits better
    than a straight line after controlling the false-discovery rate across the
    six prespecified tests.
- [ ] State that detailed knots, delta AIC, and breakpoint diagnostics are
  supplemental and that no data-selected breakpoint was promoted.
- [ ] Define the primary linear mixed models, including continuous time,
  standardized baseline parameter, parameter-by-time interaction, age, sex,
  and BMI.
- [ ] State the outcome-level random-effects decisions: patient random
  intercept for peak V̇O2; independent random intercept and slope for
  V̇E/V̇CO2. Note that variance-component likelihood-ratio tests were
  descriptive.
- [ ] Define the primary endpoint, event dating, censoring, and first-event
  convention.
- [ ] Define surgical myectomy as a separate secondary endpoint.
- [ ] Explain the V̇E/V̇CO2 exposure-by-log-time Cox analysis and why
  time-specific HRs are reported.
- [ ] State the missing-data hierarchy: E/e' + LAVi required for cohort entry;
  TRVmax available case; missing-indicator models only for secondary adjustment
  covariates, if reported; no missing-indicator substitution for an exposure,
  outcome, or event date.

## 8. Results

### Cohort and descriptive findings

- [ ] Begin with the 440-patient entry cohort, the E/e' + LAVi filter, and the
  261-patient parent cohort.
- [ ] Report the nested longitudinal and outcomes Ns immediately afterward.
- [ ] Describe the four Figure 2 groups as 116 neither threshold, 27 E/e' only,
  61 LAVi only, and 57 both.
- [ ] Report E/e' >14 as 84/261 (32.2%), LAVi >34 mL/m² as 118/261 (45.2%),
  and TRVmax >2.8 m/s as 17/145 (11.7%).
- [ ] State that the four groups differed for V̇E/V̇CO2 slope
  (Kruskal-Wallis P=0.0042) but not peak V̇O2 (P=0.57).

### Cross-sectional models

- [ ] Report E/e' and LAVi at N=261 and TRVmax at N=145 available case.
- [ ] Lead with the significant V̇E/V̇CO2 overall associations and state that
  no peak-V̇O2 overall association was significant.
- [ ] State that all six nonlinearity q values were nonsignificant
  (range 0.263-0.898), supporting association without evidence that a nonlinear
  curve or new threshold was required.
- [ ] Keep expanded complete-case models with N=111 and N=60 out of the primary
  narrative; identify them as denominator-sensitive supplemental analyses.

### Longitudinal models

- [ ] Use the compact linear mixed-model table as the primary inference.
- [ ] Report patient and observation N for every model.
- [ ] Report the primary interaction estimates:
  - E/e' with peak V̇O2: beta 1.97 percentage points/year/SD, 95% CI -0.14 to
    4.08, P=0.067;
  - LAVi with peak V̇O2: beta -0.78, 95% CI -2.95 to 1.39, P=0.479;
  - E/e' with V̇E/V̇CO2: beta 0.14 units/year/SD, 95% CI -0.16 to 0.45,
    P=0.359;
  - LAVi with V̇E/V̇CO2: beta 0.24, 95% CI -0.06 to 0.55, P=0.114.
- [ ] Label TRVmax available case and report its interaction P values as 0.326
  and 0.512.
- [ ] Do not report the former random-intercept-only LAVi P=0.0503 as the
  primary result.
- [ ] Describe GAMM trajectories as flexible sensitivity analyses; their six
  interaction P values are also nonsignificant.
- [ ] Do not claim longitudinal worsening or improvement by baseline parameter.

### Clinical outcomes

- [ ] Report E/e' and LAVi from the identical 259-patient/62-event cohort.
- [ ] Report available-case TRVmax from N=144/36 events.
- [ ] Detail the mutually exclusive first-event composition: 57 acute HF,
  0 transplant, 5 death, and no same-day multiple-component event.
- [ ] If reporting all occurrences, distinguish them from first events: two
  transplant occurrences and six death occurrences over follow-up.
- [ ] Report 33 surgical myectomies separately from the primary composite.
- [ ] State that LVOT-gradient-adjusted myectomy models were nonsignificant;
  E/e' and LAVi used N=189/19 events and TRVmax N=108/11.
- [ ] For V̇E/V̇CO2, report the time-varying sensitivity: HR 0.99 at 1 year,
  1.40 at 3 years, and 1.64 at 5 years; exposure-by-log-time P=0.0014.
- [ ] Do not interpret the forest plot's constant HR 1.11 as time invariant.

### Sensitivity analyses

- [ ] Report morphology as a sensitivity, not a primary interaction analysis.
- [ ] State that known morphology was available in 225/261 cross-sectionally
  (51 apical, 174 known nonapical) and was unknown in 36.
- [ ] State that morphology adjustment did not materially alter the E/e' or
  LAVi association with V̇E/V̇CO2.
- [ ] If missing-indicator results are included, label them as sensitivity
  analyses for incomplete adjustment covariates and do not describe them as
  recovering missing physiologic measurements.

## 9. Discussion

Use five ordered paragraphs or subsections.

1. **Main findings.** Summarize the cross-sectional ventilatory-efficiency
   associations, absence of confirmed longitudinal interactions, and outcome
   associations for E/e' and LAVi.
2. **Prior literature.** Give this the longest discussion paragraph. Integrate
   peak V̇O2, ventilatory efficiency, pulmonary and left-atrial loading, and the
   limitations of resting Doppler measurements in HCM.
3. **Original contribution.** Emphasize the common, clearly nested cohort; the
   separation of contemporaneous, longitudinal, and prognostic questions; the
   parameter-specific available-case treatment of TRVmax; and the finding that
   association did not require a new nonlinear threshold.
4. **Limitations.** Include every item below.
5. **Conclusion.** Keep short and do not overstate causality, longitudinal
   change, or threshold discovery.

### Required limitations

- [ ] Retrospective, single-center study.
- [ ] Requiring E/e' + LAVi may limit generalizability.
- [ ] TRVmax was missing in 116/261 and analyzed available case.
- [ ] Missing-indicator methods do not restore physiologic information and can
  be biased in observational data.
- [ ] Provoked and exercise LVOT gradients were unavailable.
- [ ] Expanded adjustment was constrained by covariate completeness and risk of
  overfitting.
- [ ] Morphology was unknown in 36/261.
- [ ] The source extract supports a dated acute-HF event but does not establish
  hospitalization.
- [ ] AF, VT/VF, sudden death, and alcohol septal ablation were not sufficiently
  adjudicated for endpoint inclusion.
- [ ] OMARX was exploratory and did not establish a new threshold.

## 10. Tables and figures

- [ ] Replace Table 1 with the regenerated table from the matched 177:177
  population.
- [ ] Reorder morphology rows as asymmetric septal, symmetric, and apical.
- [ ] Use `Non-obstructive HCM` and `Obstructive HCM` column labels.
- [ ] Insert the CONSORT-style cohort flow as main Figure 1.
- [ ] Confirm Figure 2 is the current descending-prevalence version with the
  E/e' + LAVi combination matrix.
- [ ] Confirm Figure 3 reports N, overall P, and nonlinearity q only.
- [ ] Decide at the manual citation/numbering pass whether the six-panel GAMM
  remains main Figure 4 or is supplemental; the compact LMM table is the
  preferred inferential presentation.
- [ ] Confirm Figure 5 panel C uses V̇E/V̇CO2 slope >30, panel D is peak V̇O2
  survival, and the legend contains the event composition and every forest
  denominator.
- [ ] Insert figures from `2_Output/Manuscript_Main/` and supplemental figures
  from `2_Output/Manuscript_Supplemental/`.
- [ ] Keep the embedded figure legends flush left and fully justified.
- [ ] Renumber supplemental figures only after deciding which binary-trajectory
  and longitudinal sensitivity figures will be submitted.

## 11. References and Paperpile safety

- [ ] Do not type over citation fields.
- [ ] When moving a sentence containing citations, move the full sentence and
  verify that each citation remains a field.
- [ ] Add or remove references through Paperpile rather than by editing the
  bibliography text.
- [ ] Do not update the bibliography until the prose and figure/table callouts
  are final.
- [ ] After the final Paperpile refresh, inspect the bibliography for duplicated
  or orphaned entries and save a new dated version.

## 12. Final numerical reconciliation

- [ ] Compare every Abstract number with the results ledger.
- [ ] Compare every Results number with the results ledger.
- [ ] Verify that Table 1 totals equal their displayed denominators.
- [ ] Verify that each longitudinal row reports both patients and observations.
- [ ] Verify that each Cox row reports patients and events.
- [ ] Confirm all threshold symbols and units: E/e' >14, LAVi >34 mL/m²,
  TRVmax >2.8 m/s, V̇E/V̇CO2 slope >30, and peak V̇O2 <80% predicted.
- [ ] Confirm that no missing value is described as normal.
- [ ] Confirm that no result calls an association nonlinear when its
  multiplicity-adjusted nonlinearity q value is nonsignificant.
- [ ] Confirm Figure 5 peak-V̇O2 HR is described as **per SD lower** fitness.
- [ ] Confirm surgical myectomy is separate from the primary composite.
- [ ] Confirm all acute-HF wording is `dated acute heart failure event` unless
  independent endpoint documentation is added.

## 13. Final submission gates

- [ ] Figure placement and numbering frozen.
- [ ] Main and supplemental table numbering frozen.
- [ ] All figure/table citations agree with the frozen numbering.
- [ ] All Ns and event counts agree with the results ledger.
- [ ] All Paperpile citations and bibliography fields remain functional.
- [ ] Clean Word comparison against the pre-revision manuscript reviewed.
- [ ] Final PDF exported manually from Word and visually inspected page by page.

## Authoritative files

- Results ledger: `2_Output/HCM_Revision_Results_Ledger.xlsx`
- Full revision rationale: `README_REVISION_PLAN.md`
- Figure audit: `2_Output/Stage6B_Figure_Audit.md`
- Main figures: `2_Output/Manuscript_Main/`
- Supplemental figures: `2_Output/Manuscript_Supplemental/`
- Main cross-sectional outputs: `2_Output/Section_03_CrossSectional/`
- Main longitudinal outputs: `2_Output/Section_04_Longitudinal/`
- Outcomes outputs: `2_Output/Section_05_Outcomes/`
- Missingness, matching, and QC outputs: `2_Output/Section_06_Ancillary/`
