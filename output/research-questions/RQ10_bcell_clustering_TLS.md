# RQ10: B-Cell Clustering Density as Improved TLS Metric

## Status: COMPLETED

## Biological Question

Can a spatially derived B-cell cluster size or B-B contact density metric better identify
TLS-like structures than the simple "fraction of B cells with ≥1 T cell neighbor" (RQ3)?
Do patients with large B-cell clusters in their tumor microenvironment have better or worse
clinical outcomes?

## Rationale

RQ3 found that the original TLS score (fraction of B cells with ≥1 T cell neighbor) is ~0.81
across all ROIs and does not distinguish inflamed vs non-inflamed tissue well.
The more informative signal was B-B co-localization (15.5% vs 4.7% of B-cell neighbors).

A connected-component clustering approach on the spatial B-cell graph gives:
- B-cell cluster sizes per ROI
- Fraction of B cells in large clusters (≥5 or ≥10 cells)
- Max cluster size per ROI

These metrics directly quantify TLS-like spatial aggregation without relying on T cell proximity.

## Methods

- Build radius-32 graph restricted to B cells only
- Find connected components (clusters)
- Compute: max cluster size, mean cluster size, fraction in clusters ≥5, ≥10
- Test vs inflammation, Gleason, survival

## Analysis Plan

### Step 1 — B-cell spatial graph and connected components
- For each ROI: filter to B cells, build radius-32 graph, find connected components
- Extract cluster size distribution metrics

### Step 2 — Association with clinical variables
- MWU: inflamed vs non-inflamed, cribriform, stromogenic
- Spearman vs Gleason grade

### Step 3 — Patient-level survival analysis
- Max B-cell cluster size per patient → Cox PH
- Compare with RQ3's TLS score: does the cluster metric have better prognostic value?

## Completed TODOs
- [x] B-cell cluster extraction (connected components in radius-32 B-cell graph)
- [x] ROI-level associations (4 significant FDR<0.1: all vs inflammation, rbc range -0.49 to -0.63)
- [x] Patient-level Cox PH (max cluster size HR=1.002, p=0.076 for PSA recurrence — nominal only)
- [x] Comparison with RQ3 TLS score (rho=0.447 with B-B contact density; cluster metrics far more discriminative)
