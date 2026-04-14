# RQ9: Spatial Niche Reconstruction via Local Composition K-Means

## Status: COMPLETED

## Biological Question

What recurrent spatial microenvironments (niches) exist in the PCa TME?
Do niche compositions and patient-level niche fractions associate with Gleason grade,
stromogenic state, and clinical outcomes?
Does the periglandular CAF1 niche (as identified in the draft paper) predict worse outcome
when measured at the niche level rather than as a continuous proximity metric?

## Rationale

RQ4 measured the fraction of CAF1-CD105+ cells within 32px of luminal epithelium
and found a counterintuitive *protective* HR=0.72. The draft paper instead identifies
a discrete periglandular CAF1 niche cluster (one of 18 niches) and reports it as adverse.

The key difference: a continuous fraction metric averages over all CAF1 cells, while a niche
cluster identifies cells whose *entire local neighborhood* resembles the periglandular pattern.
By constructing niches via k-means on the 17-dim local composition vector of each cell,
we can replicate the draft paper's approach and test whether the niche-level fractions
show the expected clinical associations.

## Methods

For each cell, its local composition vector = fraction of each meta-label in its radius-32 neighborhood.
K-means on these vectors across all cells identifies recurrent microenvironments (niches).

Niche usage per ROI = fraction of cells assigned to each niche.
Patient-level niche usage = mean across their ROIs.

## Analysis Plan

### Step 1 — Local composition vectors
- For each cell: compute 17-dim vector of meta-label fractions within radius-32 neighborhood
- Sample up to 300k cells (memory constraint) stratified across ROIs

### Step 2 — K-means clustering
- Fit k=18 (matching draft paper) on sampled cells
- Assign all cells to niches
- Characterize niches by mean composition → identify niche identities

### Step 3 — Niche-level clinical associations
- Per-ROI niche usage matrix → Kruskal-Wallis vs Gleason, MWU vs stromogenic/inflammation
- Per-patient niche usage (mean across ROIs) → Cox PH vs PSA recurrence and clinical progression

### Step 4 — Reconcile with RQ4
- Identify the "periglandular CAF1" niche (high CAF1-CD105+, high luminal epithelium composition)
- Test: does this niche fraction predict worse outcome (as the draft paper claims)?
- Compare hazard direction with RQ4's continuous periglandular fraction

### Step 5 — Visualize
- Heatmap: k=18 niche composition matrix (niches × meta-labels)
- Bar chart: fraction of cells per niche
- Forest plot: Cox HR per niche

## Completed TODOs
- [x] Local composition vectors sampled (300k cells, stratified across 476 ROIs)
- [x] K-means k=18 fit (MiniBatchKMeans)
- [x] Niche characterization (heatmap + composition table)
- [x] Clinical associations (17 significant ROI-level, no significant Cox)
- [x] Reconciliation with RQ4 (neither continuous metric nor niche fraction shows adverse direction)
