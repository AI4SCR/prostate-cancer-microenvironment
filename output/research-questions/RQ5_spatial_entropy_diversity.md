# RQ5: Spatial Entropy and Cell-Type Diversity as Spatial Organization Metrics

## Status: COMPLETED

## Biological Question

How heterogeneous is the local cellular microenvironment at the cell level (local Shannon diversity)?
Does spatial entropy differ between Gleason grades, stromogenic status, and survival outcomes?
Is local diversity associated with outcome beyond global composition?

## Rationale

Shannon diversity of cell types can be computed:
1. **Globally per ROI** — overall diversity of cell types within a region
2. **Locally per cell** — diversity within each cell's radius-32 neighborhood (ATHENA cell-level scores)

The local Shannon diversity captures spatial organization: a tumor with the same cell type proportions
but tightly segregated cell types will have low local diversity (cells are surrounded by similar cells),
while a well-mixed microenvironment will have high local diversity.

This distinguishes two spatial modes:
- **Segregated**: immune aggregates separate from epithelial/stromal zones → low local diversity
- **Mixed/infiltrated**: immune, stromal, and epithelial cells intermixed → high local diversity

## Data

- Cell type labels, spatial coordinates
- ATHENA: `ath.metrics.shannon(ad, attr='label', local=True, graph_key='radius32')`
- Aggregate to ROI-level: mean local Shannon diversity per ROI

## Analysis Plan

### Step 1 — Per-ROI global Shannon diversity
- H_roi = -sum(p_k * log(p_k)) over all cell types k in ROI

### Step 2 — Per-cell local Shannon diversity (radius-32 neighborhood)
- For each cell: compute Shannon diversity of cell type labels in its neighborhood
- Aggregate per ROI: mean, median, and SD of local diversity

### Step 3 — Cell-type-specific local diversity
- For each cell type: mean local diversity among cells of that type
  (e.g., are CAF1 cells in homogeneous neighborhoods or mixed ones?)

### Step 4 — Association with clinical variables
- Kruskal-Wallis / Mann-Whitney vs Gleason, stromogenic, inflammation, outcome
- Is local diversity orthogonal to global composition (partial correlation)?

### Step 5 — Patient-level max-pool and survival
- Max or mean local diversity per patient → Cox PH vs PSA recurrence, clinical progression

## Completed TODOs
- [ ] Global Shannon per ROI computed
- [ ] Local Shannon per cell computed (ATHENA)
- [ ] ROI-level aggregates
- [ ] Association tests
- [ ] Patient-level survival analysis
