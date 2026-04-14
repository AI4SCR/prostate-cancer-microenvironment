# RQ14: Patient-Level Clustering on Spatial Interaction Profiles

## Status: COMPLETED

## Biological Question
Do patients cluster more distinctly based on their spatial wiring (17×17 cell-cell interaction enrichment matrix from RQ6) than on composition alone (RQ7)?

## Key Results
- Best silhouette score: 0.050 at k=2 — lower than composition-based clustering (0.12, RQ7)
- Spatial interaction profiles are even more continuously distributed than cell-type compositions
- k=2 clusters differ modestly in inflammation status but not Gleason, stromogenic status, or survival (Cox p=0.17 for PSA recurrence)
- Spatial wiring does not add discriminative structure beyond composition alone
