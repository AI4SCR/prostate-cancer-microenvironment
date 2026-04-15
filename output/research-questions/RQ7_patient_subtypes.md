# RQ7: Inter-patient Heterogeneity — Patient TME Subtypes

## Status: COMPLETED

## Biological Question

Do patients cluster into distinct TME subtypes based on ROI composition?
Are these subtypes associated with Gleason grade, stromal annotations, and clinical outcome?
How much of the total variance in composition is explained by patient identity vs ROI-level variation?

## Rationale

Inter-patient heterogeneity is clinically important: if patients fall into distinct TME subtypes,
this may enable risk stratification beyond Gleason grade.

RQ2 showed that tumor epithelial phenotype has the highest ICC (0.64) — meaning patient identity
explains ~64% of variance in this cell type. A patient-level clustering approach asks whether
this patient-level signal forms discrete groups.

## Analysis Plan

### Step 1 — Patient-level composition profiles
- Aggregate ROI compositions to patient level: mean across all tumor ROIs
- Use CLR-transformed meta-label proportions (17 features)

### Step 2 — Unsupervised clustering
- Hierarchical clustering (Ward linkage) on Euclidean distance of patient CLR profiles
- Consensus clustering (bootstrap resampling) to assess stability at k=2,3,4,5 clusters
- Select optimal k by silhouette score and gap statistic

### Step 3 — Cluster characterization
- Heatmap: mean CLR composition per cluster
- Kruskal-Wallis for each cell type across clusters
- Compare with Gleason distribution, stromogenic/cribriform annotations

### Step 4 — Survival analysis per cluster
- KM curves per cluster for PSA recurrence and clinical progression
- Log-rank test; Cox PH with cluster as covariate

### Step 5 — Variance decomposition
- Quantify: what fraction of total ROI variance is between patients vs within patients?
- Extend the ICC analysis from RQ2 to understand which cell types drive inter-patient clustering

## Completed TODOs
- [ ] Patient-level CLR profiles
- [ ] Consensus clustering
- [ ] Cluster characterization
- [ ] Survival per cluster
- [ ] Variance decomposition
