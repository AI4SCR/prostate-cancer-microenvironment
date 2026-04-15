# RQ19: Cribriform Spatial Signature

## Status: COMPLETED

## Biological Question
Can spatial features distinguish cribriform from non-cribriform ROIs beyond cell-type composition? What is the discriminative AUC of a spatial classifier?

## Key Results
- Top univariate features (depleted in cribriform, FDR-corrected): CAF1–luminal contact enrichment (rbc=0.44, FDR=7.7e-7), CAF2-AR+ proportion (rbc=0.34, FDR=1.2e-4), CAF2–luminal contact enrichment (rbc=0.26, FDR=4.4e-3)
- Multi-feature logistic regression (5-fold CV): spatial + composition AUC = 0.725 ± 0.066
- Composition-only AUC = 0.668 ± 0.040; delta AUC = +0.057 from spatial features
- Top classifier features: CAF2-AR+ proportion, CAF1–luminal contact enrichment, luminal NN distance mean, CAF2–luminal contact enrichment
