# Data assets → status (summary)

"Reproducible" = regenerated in-repo from more primitive inputs, verified. "Staged" = read-only copy of a precomputed legacy output, found and in use. "Not found" = confirmed absent everywhere searched (checked `PCa_NHood`, `PCA_NHOODs_clean`, `PCa/`, and their newer project-specific staged siblings).

| Asset | Status | Used by |
|---|---|---|
| `metadata.parquet` | Reproducible (`export.py`) | most figure scripts (no `niche` column) |
| `clinical.parquet` | Reproducible (`export.py`) | most figure scripts |
| `intensity.parquet` | Reproducible (`export.py`) | `figure3_caf_umap_embedding.py` |
| `intensity_normalized.parquet` | Reproducible (`export.py`) | `figure2_umap.py`, `figure2_umap_embedding.py` |
| `figures/figure2/reducer_embedding.parquet` | Reproducible (`figure2_umap_embedding.py`) | `figure2_umap.py` |
| `figures/figure3/{config}/umap_embedding.parquet` | Reproducible (`figure3_caf_umap_embedding.py`) | `figure3_caf_umap.py` |
| `resources/colormaps.yaml` | Staged (ported verbatim) | most figure scripts |
| `resources/metalabels.yaml` | Staged (ported verbatim) | `figure4_metagroup_barplot.py` |
| `resources/non_marker_channels.txt` | Reproducible (`export.py`) | `figure2_umap.py`, `figure3_caf_umap.py` |
| `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated_v2.parquet` | Reproducible (verified byte-identical) | `figure5_niche_heatmap.R`, `figure6_inflammation_violin.R`, `figure6c_niche_composition_filtered.py`, `figureS4b_niche_mean_composition.py`, `figure7a_stromogenic_violin.R` |
| `prostate-cancer-microenvironment/5-niches/annotation/niche_annotations_v2.csv` | Reproducible (verified byte-identical) | `figure5_niche_heatmap.R`, `figure5_niche_correlation.R`, `figure6_niche_abundance_heatmap.R` |
| `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated.parquet` (no `_v2`) | Staged, unused | not used by any current script |
| `prostate-cancer-microenvironment/5-niches/annotation/niche_annotations.csv` (no `_v2`) | Staged, unused | not used by any current script |
| `prostate-cancer-microenvironment/5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet` | Staged | `figure4_metagroup_barplot.py`, `figure4c_patient_cluster_km.R`, `figureS3b_cluster_concordance.R`, `figureS3c_progression_km.R` |
| `prostate-cancer-microenvironment/5-niches/barplot_data/metadata_with_dendrogram_colors_label_tma_id.parquet` | Staged | `figureS3b_cluster_concordance.R` |
| `prostate-cancer-microenvironment/5-niches/frequencies/stacked_barplots/props_niche_pat_id.parquet` | Staged, unused | not used by any current script |
| `prostate-cancer-microenvironment/5-niches/frequencies/stacked_barplots/props_niche_tma_id.parquet` | Staged | `figure5_niche_correlation.R` |
| `prostate-cancer-microenvironment/5-niches/visualization/composition/niche_heatmap_data.parquet` | Staged | `figure5_niche_heatmap.R` |
| `prostate-cancer-microenvironment/5-niches/visualization/composition/mean_celltype_composition_per_niche.parquet` | Staged | `figure5_niche_heatmap.R` |
| `prostate-cancer-microenvironment/5-niches/visualization/composition/median_celltype_composition_per_niche.parquet` | Staged | `figure5_niche_heatmap.R` |
| `prostate-cancer-microenvironment/5-niches/visualization/composition/niche_abundance_stats.parquet` | Staged | `figure5_niche_heatmap.R` |
| `prostate-cancer-microenvironment/PCA_NHOODs_clean/niche_annotations_revised.xlsx` | Staged (supplied by user 2026-07-31) | `figure5_niche_annotation.py` |
| `prostate-cancer-microenvironment/PCA_NHOODs_clean/robustness/utils/` (code) | Staged | `figure5_niche_clustering.py` |
| `PCa_NHood/CellCellNeighborhoods/{cell_metadata,metadata}.parquet`, `graph_type=radius-radius=32/data.parquet` | Staged, already validated | `figure5_niche_clustering.py` |
| `prostate-cancer-microenvironment/5-niches/visualization/interactions/redo/per_niche_lfc_above_median/dataframes_v2/{niche}.parquet` | Staged (niches 2/8/9) | `figure7def_circos_plots.py` |
| `prostate-cancer-microenvironment/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas/{sample_id}.pkl` | Staged (476 files, correct location, 100% ID match) | not currently read by any script (upstream of the circos pipeline, which reads precomputed `dataframes_v2/` instead) |
| **`5-niches/frequencies/niche_frequencies_per_tma_id.parquet`** | **Not found** | `figure6_niche_abundance_heatmap.R` (Fig 6a) -- blocked |
| **`cell_annotation.parquet`** | **Not found** | `figure6_km_niche6.R`, `figure6e_immune_risk_score_km.R`, `figure7b_niche9_km.R`, `figure7c_myCAF_km.R`, `figureS6bc_niche_km.R` -- all blocked |
