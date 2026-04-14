# RQ3: Immune Neighborhood Composition and TLS-like Structure Characterization

## Status: COMPLETED

## Biological Question

What are the local cellular neighborhoods of B cells and T cell subsets? Can we identify
TLS-like structures from the spatial data, and how do these relate to inflammation status
and clinical outcome?

## Rationale

The draft paper identifies a TLS-like B-cell/T-cell aggregate and links "composite immune
niche burden" to worse survival. However, TLS maturation markers were not in the panel.
The IMC data contains CD20 (B cells), CD3, CD4, CD8a, FoxP3 (T cell subsets), and CD68
(macrophages). We can:
1. Characterize local neighborhoods of CD20+ B cells using radius graphs
2. Test whether B-cell-rich neighborhoods are enriched for helper T cells (CD4+) vs regulatory
   T cells (FoxP3+) vs cytotoxic T cells (CD8a+) — this distinguishes mature vs suppressive TLS
3. Compute a per-ROI TLS-like score (density of B-T co-localizations)
4. Ask whether TLS-like score tracks with inflammation annotation and outcome

## Data

- Spatial coordinates: `features/spatial/filtered-annotated/*.parquet`
- Cell type labels: `metadata/filtered-annotated/*.parquet`
- Clinical: is_tumor, inflammation, psa_progr, clinical_progr

## Analysis Plan

### Step 1 – Build radius graphs per ROI
- Use ATHENA radius graph (radius=32 pixels, as in draft paper)
- For each ROI: build graph, assign cell type labels to nodes

### Step 2 – B-cell neighborhood profiles
- For each B cell (immune-B-cells(CD20+)):
  extract proportion of each cell type in its radius-32 neighborhood
- Aggregate to ROI-level mean neighborhood composition

### Step 3 – TLS-like score
- TLS-like score = fraction of B cells whose neighborhood contains ≥1 helper T cell (CD3+CD4+)
- Alternative: fraction of B cells co-localizing with any T cell subtype
- Compute per ROI

### Step 4 – Association with inflammation and outcome
- Mann-Whitney U: TLS-like score in inflamed vs non-inflamed ROIs
- Cox proportional hazards: TLS-like score vs psa_progr_time / clinical_progr_time
- Kaplan-Meier stratified by TLS-like score tertile

### Step 5 – Regulatory T cell enrichment in immune neighborhoods
- Is the ratio of FoxP3+ Tregs to total T cells higher in B-cell-rich vs B-cell-poor ROIs?
- High Treg:T ratio near B cells suggests immune suppression rather than activation

## Expected Outputs

- `output/figures/RQ3/` – neighborhood heatmaps, KM plots, scatter plots
- `output/tables/RQ3/` – TLS score per ROI, association statistics

## Completed TODOs

- [ ] Radius graph construction
- [ ] B-cell neighborhood profiles
- [ ] TLS-like score computation
- [ ] Association with inflammation and outcome
- [ ] Treg enrichment analysis

## Limitations

- Radius-32 pixel graphs depend on image resolution (~1 μm/pixel for IMC); 32 pixels ≈ 32 μm
- No TLS maturation markers (CXCL13, CCL19/21 not in panel) — spatial co-localization is proxy
- Small n for cribriform/inflammation-positive ROIs
