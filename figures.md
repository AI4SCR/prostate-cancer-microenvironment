# Figure → script mapping

Which script produces each panel of each main-text figure. Blank "Script" = not implemented in this repo yet; "—" = not code-derived / no source found (see `figure_script_mapping.md` for detail).

| Figure | Panel | Description | Script |
|---|---|---|---|
| 1 | a-c | Workflow schematic, antibody panel, representative IMC images | — (BioRender + raw image crops, not code-derived) |
| 2 | a | 34-cell-type heatmap | `scripts/figures/figure2_cell_type_heatmap.R` |
| 2 | b | UMAP of all cells (compartment/cell type/patient/markers) | `scripts/figures/figure2_umap.py` (+ `scripts/data/figure2_umap_embedding.py`) |
| 3 | a | UMAP of CAF cells by subcluster | `scripts/figures/figure3_caf_umap.py` (+ `scripts/data/figure3_caf_umap_embedding.py`) |
| 3 | b | CAF subcluster marker-expression heatmap | `scripts/figures/figure3_caf_heatmap.R` |
| 3 | c-e | Representative CAF/fibrocyte ROIs | not implemented |
| 4 | a | Patient-level hierarchical clustering (P1-P6) | `scripts/figures/figure4_patient_clustering.py` |
| 4 | b | Mean metagroup distribution per patient cluster | not implemented |
| 4 | c | KM survival by patient cluster | `scripts/figures/figure4_survival.R` |
| 4 | d-e | Cox PH: overall survival / progression, per cell type | `scripts/figures/figure4_survival.R` |
| 5 | a | Niche x cell-type z-score heatmap | `scripts/figures/figure5_niche_heatmap.R` (+ `figure5_niche_clustering.py`, `figure5_niche_annotation.py`) |
| 5 | b | Niche pairwise Spearman correlation | `scripts/figures/figure5_niche_correlation.R` |
| 5 | c | Representative cores + stacked barplots | not implemented |
| 6 | a | Per-core niche-abundance heatmap | `scripts/figures/figure6_niche_abundance_heatmap.R` |
| 6 | b | KM survival by niche 6 abundance | `scripts/figures/figure6_km_niche6.R` |
| 6 | c | Stacked barplot, niches 6/16/17/18 composition | not implemented |
| 6 | d | Niche abundance vs inflammation status | `scripts/figures/figure6_inflammation_violin.R` |
| 6 | e | KM by immune-niche risk score | not implemented |
| 7 | a | Stromal niche abundance vs stromogenic status | not implemented |
| 7 | b-c | KM by niche 9 / myCAF abundance | not implemented |
| 7 | d-f | Circos plots, cell-cell interactions (niches 2/8/9) | not implemented |
| 7 | g-h | Representative cores, niche 8 / niche 9 | not implemented |
