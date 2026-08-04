# Legacy data assets (`data/legacy/`)

None of these are reproduced by any script in this repo — they're staged, read-only copies of precomputed/manually-curated outputs from the original (pre-migration) pipeline. See `figure_script_mapping.md` and `REPRODUCIBILITY.md` for context.

| Asset | Produced by | Used by |
|---|---|---|
| `5-niches/annotation/clusters_annotated_v2.parquet` | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R`, `figure6_inflammation_violin.R`, `figure6_km_niche6.R` |
| `5-niches/annotation/niche_annotations_v2.csv` | `figure5_niche_annotation.py` (from the xlsx below) | `figure5_niche_heatmap.R`, `figure5_niche_correlation.R`, `figure6_niche_abundance_heatmap.R` |
| `5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet` | unknown (precomputed, no reproducing script) | `figure4_survival.R` |
| `5-niches/frequencies/stacked_barplots/props_niche_pat_id.parquet` | unknown (precomputed, no reproducing script) | `figure6_niche_abundance_heatmap.R` |
| `5-niches/frequencies/stacked_barplots/props_niche_tma_id.parquet` | unknown (precomputed, no reproducing script) | `figure5_niche_correlation.R` |
| `5-niches/visualization/composition/niche_heatmap_data.parquet` | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/mean_celltype_composition_per_niche.parquet` | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/median_celltype_composition_per_niche.parquet` | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/niche_abundance_stats.parquet` | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `PCA_NHOODs_clean/niche_annotations_revised.xlsx` | manually curated by hand | `figure5_niche_annotation.py` |
| `PCA_NHOODs_clean/robustness/utils/` (code, not data) | legacy codebase | `figure5_niche_clustering.py` (imports `utils.clustering`, `utils.visualization`) |
| `PCa_NHood/CellCellNeighborhoods/` | unknown (undocumented upstream neighborhood-graph computation) | `figure5_niche_clustering.py` |

## Recovered legacy UMAP embeddings (`data/figures/`)

Not legacy in the "external/manually curated" sense above — extracted from the original pipeline's `reducer.pkl` files (pixi env + a pynndescent compat patch, see session history) and moved out of `data/legacy/` since they're figure-script-scoped caches, not consolidated external data.

| Asset | Produced by | Used by |
|---|---|---|
| `figures/figure2_umap/umap_embeddings.parquet` | pkl-to-parquet extraction (not a repo script) | **unused** — `figure2_umap_embedding.py` always recomputes UMAP from scratch instead |
| `figures/figure3_caf_umap/{excl_markers,caf_markers_only}/umap_embeddings.parquet` | pkl-to-parquet extraction (not a repo script) | **unused** — `figure3_caf_umap_embedding.py` always recomputes UMAP from scratch instead |
