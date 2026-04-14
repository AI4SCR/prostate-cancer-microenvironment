# RQ8: Paired High-Grade vs Low-Grade Core Analysis Within Patients

## Status: COMPLETED

## Biological Question

Within the same patient, do high-Gleason TMA cores have a systematically different TME than
low-Gleason cores — in composition, spatial diversity, and cell-cell interactions?
Which cell types show the most consistent within-patient shift from low- to high-grade?

## Rationale

99 patients have ROIs with ≥2 distinct Gleason grade groups. TMA cores are annotated as
"high" or "low" Gleason per TMA coordinate (`gleason_pattern_tma_core`). This enables a
matched within-patient comparison that controls for patient-level confounders (genetics,
treatment, co-morbidities).

This is a more powerful test than cross-patient comparisons because it removes inter-patient
variance as a confounder.

## Data

- `gleason_pattern_tma_core`: "high", "low", "only_1", "unknown"
- Filter to patients with at least one "high" and one "low" TMA core
- Use paired statistical tests (Wilcoxon signed-rank)

## Analysis Plan

### Step 1 — Identify paired patients
- Filter to patients with ≥1 "high" and ≥1 "low" TMA core
- If multiple high or low cores per patient: average the composition within each grade category

### Step 2 — Paired composition comparison
- For each cell type: Wilcoxon signed-rank test (high vs low core, matched within patient)
- Effect size: paired rank-biserial correlation
- FDR correction

### Step 3 — Paired spatial diversity comparison
- Compute per-ROI local Shannon diversity (from RQ5)
- Paired comparison: is local diversity higher in high- vs low-grade cores?

### Step 4 — Paired interaction comparison
- Use per-ROI interaction enrichment scores (from RQ6)
- Which cell-cell interactions are consistently different between high- and low-grade cores
  within the same patient?

### Step 5 — Visualize paired changes
- Lollipop plot: mean high-grade minus mean low-grade per cell type (±95% CI)
- Sankey / slope plot connecting individual patients' high vs low cores

## Completed TODOs
- [ ] Paired patient identification
- [ ] Paired composition Wilcoxon tests
- [ ] Paired spatial diversity comparison
- [ ] Paired interaction analysis
- [ ] Visualization
