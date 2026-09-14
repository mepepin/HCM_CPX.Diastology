# Stage 6B: Final Figure Audit

Audit date: 2026-09-13

## Scope and locked reference

This audit compares the current PDFs in `2_Output/Manuscript_Main/` and
`2_Output/Manuscript_Supplemental/` with the authoritative denominators,
definitions, and analysis hierarchy in `HCM_Revision_Results_Ledger.xlsx`.
It is an audit only: no figure geometry, styling, or manuscript Word file was
changed during this stage.

## Global quality-control findings

- All audited PDFs are one-page vector files and all fonts are embedded.
- The current main figures have no visible clipping, overlapping axis text,
  misplaced velocity dot, or overlapping oxygen subscript.
- Figures 2, 3, 4, and 5 have flush-left, visually justified embedded legends.
- The locked hierarchy is represented correctly: N=261 cross-sectional parent
  cohort, N=105 longitudinal subcohort, and N=259 outcomes subcohort.
- Missing TRVmax is not displayed as normal in the main figures. TRVmax uses
  N=145 cross-sectionally and N=144/36 events in its outcomes model.
- The remaining changes are content labels, denominator disclosure, and
  supplemental-asset cleanup. None requires redesigning the figures.

## Main figures

| Figure | Audit status | Findings | Required action before freeze |
|---|---|---|---|
| Figure 1: cohort flow | Pass | Shows 440 phenotypic HCM patients, the upstream E/e' + LAVi filter, N=261 parent cohort, TRVmax N=145 available case, longitudinal N=105/294 CPETs, outcomes N=259/62 events, and the 1,015 -> 243 -> 242 -> 177 comparator branch. Comparator QC and resting-gradient language are explicit. | None. |
| Figure 2: measured parameter patterns | Revise labels/order only | Primary counts are correct: E/e' 84/261, LAVi 118/261, TRVmax 17/145; four E/e' + LAVi groups are 116, 27, 61, and 57. Panels A/B are switched as requested, and C/D use the bottom combination matrix. Panel A is not ordered by decreasing prevalence and still says `Average e' (age-adjusted)`. | Sort Panel A by decreasing prevalence: age-calibrated e' 66.5%, LAVi 45.2%, E/e' 32.2%, PV S/D 27.9%, TRVmax 11.7%, MV deceleration time 10.0%, MV E/A 8.0%. Rename the e' row to `Age-calibrated average e' below age-specific reference limit` (or a compact equivalent). Preserve all other design features. |
| Figure 3: cross-sectional RCS | Pass | E/e' and LAVi panels use N=261; TRVmax panels use N=145 and are marked available case. Boxes contain only N, overall P, and nonlinearity q. Values match the ledger. The legend explains overall association and nonlinearity in plain language. | None. |
| Figure 4: longitudinal trajectories | Statistical content passes; terminology/placement unresolved | Denominators are correct: 294/105 and 293/105 for E/e' and LAVi; 176/63 and 175/63 for TRVmax. All six interaction P values are nonsignificant. Lower/upper baseline halves are shown as requested. The title still says `Diastolic Dysfunction Parameters`. The file remains in the main-figure manifest even though the revision plan designates the compact LMM table as the preferred main result and this GAMM display as supplemental. | Rename the title to `Longitudinal Cardiopulmonary Trajectories Across LV Diastolic Parameters`. Freeze the placement decision: preferred approach is the compact LMM table in the main paper and this unchanged trajectory design in the supplement. |
| Figure 5: clinical outcomes | Revise labels/denominator disclosure only | Panels C and D correctly use VE/VCO2 slope >30 and peak VO2 <80% predicted. The legend correctly defines the composite and gives 57 acute-HF, 0 transplant, and 5 death first events; myectomy is separate. Panel E labels the inverse-coded peak-VO2 exposure merely as `pVO2`, although HR 1.34 is per 1 SD **lower** peak VO2. Panel E also does not disclose all parameter-specific N/event denominators. `Max LVOT Gradient` does not state that the gradient is resting. | Rename `pVO2` to `Peak VO2 (per SD lower)` or `Peak VO2 deficit`. State the forest denominators in the legend: E/e' 259/62, LAVi 259/62, TRVmax 144/36, resting LVOT gradient 189/44, septal thickness 259/62, peak VO2 259/62, and VE/VCO2 slope 259/62. Rename `Max LVOT Gradient` to `Resting LVOT gradient`. Keep the proportional-hazards warning for VE/VCO2 and interpret the time-varying analysis. |

## Supplemental figures and packaging

| Asset | Audit status | Required action |
|---|---|---|
| Matched-cohort validation | Minor terminology revision | Replace the display abbreviation `CON` with `Non-HCM reference` or `Reference`. The caption already defines the comparator correctly and confirms resting LVOT gradient. Keep the existing design. |
| Supplemental longitudinal GAMM | Same issue as Figure 4 | Replace `Diastolic Dysfunction Parameters` with `LV Diastolic Parameters`. Preserve the six-panel design. |
| Observed longitudinal trajectories | Pass with typography cleanup | Keep as the transparent raw-data companion. Standardize axis/caption typography to V-dot O2 and V-dot E/V-dot CO2. |
| E/e', LAVi, and TRVmax binary trajectory figures | Not publication-ready | Literal HTML (`<sub>2</sub>`) appears in the rendered axes. Replace `Normal/Abnormal` with direct threshold labels (for example, `E/e' <=14` and `E/e' >14`) because one threshold does not establish globally normal or abnormal diastolic function. Widen/simplify the model-summary table so interaction labels and statistics do not collide. |
| LAVi unadjusted GAMM | Not publication-ready | The panel headings/table alignment are clipped or crowded. Apply the same typesetting corrections as the binary trajectory figures or omit this redundant asset. |
| Crude-versus-adjusted Cox figure | Direction label required | `Peak VO2` has HR >1 because the modeled exposure is per SD lower peak VO2. Relabel accordingly. Rename `Max LVOT Gradient` to `Resting LVOT gradient`. |
| All-parameter RCS figure | Usable as technical supplement | Dense but legible. Retain only if the supplement needs the full exploratory parameter family. Ensure its external legend defines overall P, nonlinearity q, and Delta AIC. |
| RCS-versus-OMARX figure | Typesetting revision required | TRVmax x-axis tick labels are crowded together. Reduce tick density or widen those facets. |
| Delta-e' coupling figure | Legend incomplete | Add a flush-left, justified explanatory legend and define the point-color scale. |
| OMARX clinical summary (expected Figure S2) | Missing | `Figure_S2_OMARX_ClinicalSummary.pdf` is listed in the figure manifest but is absent from `Manuscript_Supplemental`. Generate it or remove it from the final manifest and renumber dependent supplemental figures. |

## File-manifest cleanup

- `Figure1_MatchedControls_and_BaselineCharacteristics.*` remains in
  `Manuscript_Main` even though the current manifest also assigns the matched
  validation figure to the supplement. Remove the stale main-folder duplicate
  only after the final figure placement is frozen.
- Assign final supplemental figure numbers only after deciding whether to keep
  the OMARX and binary-trajectory assets. Several current filenames use generic
  `FigureS_` prefixes and should not be treated as submission-ready numbering.

## Stage 6B disposition

The main analytical figures are internally consistent with the locked results.
Figure 1 and Figure 3 are ready to freeze. Figures 2, 4, and 5 require limited
text/order corrections that preserve their present design. Supplemental figure
packaging and three legacy binary-trajectory PDFs require cleanup before the
figure gate can be marked passed.

## Stage 6B.1 resolution

Completed 2026-09-13. The figure gate now passes.

- Figure 2 was reordered by decreasing prevalence and now labels average e' as
  age-calibrated against an age-specific reference limit.
- Figure 1 now carries its full flush-left, justified explanatory legend in the
  submission PDF; the cohort-flow diagram itself was not restyled.
- Figure 4 now refers to LV diastolic parameters; all original six-panel
  geometry and lower/upper-half trajectories were preserved.
- Figure 5 now displays `Lower peak V̇O2`, `Resting LVOT gradient`, the
  literature cutoff of VE/VCO2 slope >30, and parameter-specific N/event counts
  in the embedded legend.
- Matched-comparator displays now use `Reference`/`non-HCM reference`, define
  obstruction from the resting LVOT gradient, and retain the original design.
- Binary E/e', LAVi, and TRVmax trajectory figures now use explicit threshold
  groups, correct oxygen typography, and a non-overlapping model-summary table.
- The crude-versus-adjusted figure now states the direction of inverse-coded
  peak V̇O2, identifies resting LVOT gradient, and uses one coherent model key.
- Figure S4 now defines the point-color scale in its flush-left, justified
  embedded legend.
- The absent/stale optional OMARX figures were removed from the submission
  manifest. Existing exploratory output was archived rather than deleted.
- The stale matched-comparator copy was removed from `Manuscript_Main` and
  archived rather than deleted.
- All changed PDFs were rerendered, inspected from high-resolution page images,
  and checked for embedded fonts. No clipping or text collisions remained.

The primary linear mixed-model table remains the preferred longitudinal
inferential presentation. Whether the six-panel GAMM display is ultimately
numbered as a main or supplemental figure should be finalized during the manual
manuscript pass so figure citations and numbering remain synchronized.
