# RQ4: CAF Subtype Spatial Topology Relative to Tumor Glands

## Status: COMPLETED

## Biological Question

Do CAF1 (CD105+) and CAF2 (AR+) subtypes occupy systematically different spatial positions
relative to tumor epithelium and other stromal/immune cells? Is the periglandular CAF1 niche
identifiable as a spatial topology pattern, and does it scale with clinical outcome independently
of CAF1 proportion?

## Rationale

The draft paper reports that a "periglandular CD105high CAF1 niche" has the strongest adverse
clinical association. This suggests that the spatial positioning of CAF1 cells — specifically
their proximity to luminal epithelium — matters beyond their raw abundance.

We can test this by computing:
1. Proximity index: mean distance from CAF1 cells to nearest tumor epithelial cell vs mean
   distance from CAF2 cells to nearest tumor epithelial cell
2. Interaction enrichment: observed vs expected frequency of CAF1-epithelial contacts
3. Whether CAF1–epithelial spatial proximity score (beyond proportion) stratifies outcome

## Data

- Spatial coordinates: `features/spatial/filtered-annotated/*.parquet`
- Cell type labels: `metadata/filtered-annotated/*.parquet`  
- Intensity: `features/intensity/filtered-annotated/*.parquet` (CD105, AR markers)
- Clinical: is_tumor, gleason_grp, psa_progr, clinical_progr

## Analysis Plan

### Step 1 – CAF–epithelial proximity index
- For each ROI: compute mean nearest-neighbor distance from each CAF1 (CD105+) cell to
  nearest luminal epithelial cell
- Compare with mean NN distance from CAF2 (AR+) cells to nearest luminal epithelial cell
- Per-ROI score: CAF1 proximity index = mean NN distance (lower = more periglandular)

### Step 2 – Interaction enrichment (permutation-based)
- For each ROI: count CAF1-epithelial edges in radius-32 graph
- Permute cell type labels 1000x, compute expected count
- Enrichment score = observed / expected
- Compare with CAF2-epithelial enrichment

### Step 3 – CAF1 proximity index vs outcome
- Cox model: CAF1 proximity index vs psa_progr
- Kaplan-Meier by median split
- Test independence from CAF1 proportion (partial correlation)

### Step 4 – CD105 expression gradient
- Using intensity data: is CD105 expression highest in cells adjacent to epithelium?
- Compare CD105 intensity in CAF1 cells with ≥1 epithelial neighbor vs no epithelial neighbor

## Expected Outputs

- `output/figures/RQ4/` – proximity boxplots, survival curves, scatter plots
- `output/tables/RQ4/` – proximity scores, enrichment, survival statistics

## Completed TODOs

- [ ] CAF–epithelial proximity index
- [ ] Permutation-based interaction enrichment
- [ ] Proximity index vs outcome analysis
- [ ] CD105 expression gradient analysis

## Limitations

- Radius-32 neighborhoods may not fully capture the periglandular niche geometry
- CAF subtypes are defined by discrete labels; continuous CD105/AR expression may be more
  appropriate for gradient analyses
- Multiple ROIs per patient: patient clustering not modeled for survival
