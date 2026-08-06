# Figure → script mapping (summary)

| Figure | Panel | Script | Ported from | Assets |
|---|---|---|---|---|
| 1 | a-c | | | |
| 2 | a | `figure2_cell_type_heatmap.R` | `000_paper/04_heatmaps/2-cell-types-heatmap.R` | intensity_normalized.parquet, metadata.parquet |
| 2 | b | `figure2_umap.py` | `000_paper/02_umaps/0-umaps.py` | metadata.parquet, intensity_normalized.parquet, clinical.parquet, figure2 embedding cache |
| 3 | a | `figure3_caf_umap.py` | `000_paper/02_umaps/0-umaps-cafs.py` | metadata.parquet, intensity.parquet, clinical.parquet, figure3 embedding cache |
| 3 | b | `figure3_caf_heatmap.R` | `000_paper/04_heatmaps/2-cell-types-heatmap.R` (`heatmap.caf()`) | intensity_normalized.parquet, metadata.parquet |
| 3 | c-e | | | |
| 4 | a | `figure4_patient_clustering.py` | `sync_paper/05-heterogeneity/stacked-frequencies-label.py` (`pat_id` branch) | metadata.parquet, clinical.parquet |
| 4 | b | `figure4_metagroup_barplot.py` | `sync_paper/05-heterogeneity/stacked-frequencies-label.py` (tail block) | `output/figures/figure4/figure4a_patient_composition.parquet` (Fig 4a's own output), resources/metalabels.yaml |
| 4 | c | `figure4c_patient_cluster_km.R` | `sync_paper/03_survival/patient_risk_group_km.R` | `prostate-cancer-microenvironment/5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet`, clinical.parquet |
| 4 | d-e | `figure4de_cox_hazard_ratio.R` | `000_paper/03_survival/survival-proportions.r` | metadata.parquet, clinical.parquet |
| 5 | a | `figure5_niche_heatmap.R` | `sync_paper/06-spatial-niches/composition/z_score_composition_vis.R` | `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated_v2.parquet`, `.../niche_annotations_v2.csv`, `.../visualization/composition/niche_heatmap_data.parquet` |
| 5 | b | `figure5_niche_correlation.R` | `sync_paper/06-spatial-niches/abundance/correlation_frequencies.R` | `prostate-cancer-microenvironment/5-niches/frequencies/stacked_barplots/props_niche_tma_id.parquet`, `.../annotation/niche_annotations_v2.csv`, clinical.parquet |
| 5 | c | | | |
| 6 | a | `figure6_niche_abundance_heatmap.R` | `sync_paper/06-spatial-niches/abundance/heatmap_frequencies.R` | `5-niches/frequencies/niche_frequencies_per_tma_id.parquet`, clinical.parquet |
| 6 | b | `figure6_km_niche6.R` | `sync_paper/03_survival/niche_km.R` | clinical.parquet, `cell_annotation.parquet` |
| 6 | c | `figure6c_niche_composition_filtered.py` | `sync_paper/06-spatial-niches/composition/niche_composition.py` (niches 6/16/17/18) | `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated_v2.parquet` |
| 6 | d | `figure6_inflammation_violin.R` | `sync_paper/06-spatial-niches/histology/inflammation_vis.R` | `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated_v2.parquet`, clinical.parquet |
| 6 | e | `figure6e_immune_risk_score_km.R` | `sync_paper/03_survival/inflammatory_niche_risk_group.R` | clinical.parquet, `cell_annotation.parquet` |
| 7 | a | `figure7a_stromogenic_violin.R` | `sync_paper/06-spatial-niches/histology/stromogenic_vis.R` | `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated_v2.parquet`, clinical.parquet |
| 7 | b | `figure7b_niche9_km.R` | `sync_paper/03_survival/niche_km.R` (niche 9) | clinical.parquet, `cell_annotation.parquet` |
| 7 | c | `figure7c_myCAF_km.R` | `sync_paper/03_survival/celltype_km.R` (myCAF) | clinical.parquet, `cell_annotation.parquet` |
| 7 | d-f | `figure7def_circos_plots.py` | `sync_paper/06-spatial-niches/interactions/interaction_compute_circos.py` | `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated.parquet` (no `_v2`), `prostate-cancer-microenvironment/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas/{sample_id}.pkl` (niches `luminal_infiltrated`, `luminal_CAF1(CD105High)`, `tumor_CAF1(CD105High)`) |
| 7 | g-h | | | |
| S1 | a | `figureS1a_cohort_summary.R` | `000_paper/01_clinical_metadata/clinical.r` | clinical.parquet |
| S1 | b | `figureS1b_gleason_concordance.R` | `sync_paper/05-heterogeneity/patient-core-heterogeneity.R` (Gleason section) | clinical.parquet |
| S2 | a-c | `figure2_umap.py` (compartment loops) | `000_paper/02_umaps/0-umaps.py` | metadata.parquet, intensity_normalized.parquet, clinical.parquet |
| S3 | a | `figureS3a_tma_stacked_barplot.py` | `sync_paper/05-heterogeneity/stacked-frequencies-label.py` (`tma_id` branch) | metadata.parquet, clinical.parquet |
| S3 | b | `figureS3b_cluster_concordance.R` | `sync_paper/05-heterogeneity/patient-core-heterogeneity.R` (cluster section) | `prostate-cancer-microenvironment/5-niches/barplot_data/metadata_with_dendrogram_colors_label_{pat_id,tma_id}.parquet` |
| S3 | c | `figureS3c_progression_km.R` | `sync_paper/03_survival/patient_risk_group_km.R` | `prostate-cancer-microenvironment/5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet`, clinical.parquet |
| S4 | a | | | |
| S4 | b | `figureS4b_niche_mean_composition.py` | `sync_paper/06-spatial-niches/composition/niche_composition.py` | `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated_v2.parquet` |
| S5 | a | | | |
| S5 | b | `figureS5b_inflammation_km.R` | `sync_paper/03_survival/clinical_risk_group_km.R` (inflammation section) | clinical.parquet |
| S6 | a | `figureS6a_stromogenic_km.R` | `sync_paper/03_survival/clinical_risk_group_km.R` (stromogenic section) | clinical.parquet |
| S6 | b-c | `figureS6bc_niche_km.R` | `sync_paper/03_survival/niche_km.R` (niches 2/8) | clinical.parquet, `cell_annotation.parquet` |
| S7 | a | | `sync_paper/06-spatial-niches/construction/kmeans_clustering.py` (ARI robustness section) | `prostate-cancer-microenvironment/5-niches/annotation/clusters_annotated_v2.parquet` (via kmeans re-run) |
