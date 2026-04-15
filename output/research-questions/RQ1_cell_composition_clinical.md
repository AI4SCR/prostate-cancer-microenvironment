# RQ1: Cell Type Composition Variation Across Clinical and Pathological Groups

## Status: COMPLETED

## Biological Question

Do cell type proportions at the ROI level differ systematically across Gleason grade group,
stromogenic status, histological inflammation, cribriform architecture, and survival outcome
(PSA recurrence, clinical progression)?

This is the foundational characterization question. The draft paper describes 34 cell types and
18 niches and reports clinical associations for specific niches, but a systematic multi-group
composition analysis across all clinical variables — with effect sizes and FDR correction — is
not shown in the existing reports.

## Rationale

Cell type composition is the simplest readout from IMC data and directly connects to known biology:
- Reactive stroma / stromogenic ROIs should be enriched in CAF1 CD105+ and depleted in CAF2 AR+
- Inflamed ROIs should show elevated T/B cell proportions
- Cribriform glands should differ in luminal epithelial marker expression
- High Gleason grade should correlate with altered epithelial/stromal balance

## Data

- Cell type labels: `metadata/filtered-annotated/*.parquet` (label, main_group, meta_label)
- Clinical: `metadata/clinical.parquet`
- Filter to tumor ROIs only (is_tumor == 'yes')

## Analysis Plan

### Step 1 – Compute ROI-level cell type proportions
- Load all cell-level metadata (label column)
- For each ROI: compute fraction of cells in each label, main_group, and meta_label
- Result: DataFrame [n_rois × n_cell_types]

### Step 2 – CLR transformation
- Apply centered log-ratio (CLR) transform to proportions (add pseudocount 1e-5 before log)
- This is appropriate for compositional data and is used in the draft paper

### Step 3 – Stratified comparison across clinical groups
For each cell type and each clinical variable:
- Kruskal-Wallis test (multi-group) or Mann-Whitney U (two-group)
- Report effect size (rank-biserial correlation or eta-squared)
- FDR correction (Benjamini-Hochberg) per clinical variable

Groups to test:
- `gleason_grp`: 1–5 (ordinal Spearman correlation also tested)
- `stromogenic_smc_loss_reactive_stroma_present`: yes/no
- `inflammation`: yes/no
- `cribriform`: yes/no
- `psa_progr`: 0/1
- `clinical_progr`: 0/1

### Step 4 – Visualization
- Heatmap: mean CLR-proportions × clinical group, rows = cell types, columns = groups
- Volcano plots: effect size vs -log10(FDR) per clinical variable
- Boxplots for top significant cell types per variable

### Step 5 – Summary table
- Report FDR-significant associations (q < 0.1) with direction and effect size

## Expected Outputs

- `output/figures/RQ1/` – heatmaps, volcano plots, boxplots
- `output/tables/RQ1/` – FDR-corrected association table

## Completed TODOs

- [x] Dataset loaded successfully
- [x] Clinical variable distributions inspected
- [x] ROI-level proportion matrix computed
- [x] CLR transformation applied
- [x] Statistical tests run across all clinical groups
- [x] Figures generated

## Key Findings

- **CAF2-AR+** is the most significantly associated cell type: decreases with Gleason grade
  (Spearman rho=-0.29, FDR=4e-7), strongly depleted in cribriform ROIs (rbc=-0.38, FDR=5e-5),
  and depleted in stromogenic ROIs. This suggests loss of the AR+ smooth-muscle-like CAF state
  is a hallmark of high-grade and reactive-stroma disease.
- **Inflammation** has the richest association profile (13 meta-label associations FDR<0.1):
  B cells (rbc=+0.60) and immune-mixed (rbc=+0.57) are hugely enriched; epithelial-luminal
  and stromal compartments are depleted.
- **Gleason grade**: 13 significant associations; stromal main group decreases (rho=-0.22),
  undefined increases with Gleason.
- **PSA recurrence**: marginal signal at composition level — epithelial-luminal fraction
  decreases in PSA progressors (rbc=-0.15, FDR~0.08). Spatial analysis captures more.
- **Clinical progression**: no FDR<0.1 associations at meta-label level, but epithelial
  main group is nominally reduced in progressors.

## Limitations

- Multiple ROIs per patient → patient-level random effects not modeled; treat ROI as unit for now
- Gleason grading per core may differ from patient-level grade
- Survival follow-up varies; recurrence event rate is reasonable (psa_progr: 161/542, clinical_progr: 90/542)
