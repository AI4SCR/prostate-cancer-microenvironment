# RQ18: Dense Periglandular CAF1 Niche vs Diffuse CAF1 Proximity

## Status: COMPLETED

## Biological Question
Are dense, organised CAF1-CD105+ clusters around glands (DBSCAN) structurally and prognostically distinct from diffuse periglandular proximity (RQ4)?

## Key Results
- DBSCAN (eps=40, min_samples=4) on periglandular CAF1-CD105+ cells
- Mean fraction of CAF1 in dense niche: 26.6%
- `periglandular_frac` vs `frac_caf1_in_dense_niche`: rho=0.075, p=0.10 — essentially uncorrelated
- `periglandular_frac` increases with inflammation (rbc=0.22, FDR=0.029); decreases in cribriform (rbc=-0.17, FDR=0.088)
- `max_dense_cluster_size` increases with Gleason (rho=0.12, FDR=0.041)
- Survival: all dense niche scores p>0.28 after Gleason adjustment; neither metric achieves significance
