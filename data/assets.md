# data/ asset dictionary

`DATA_DIR` (default: this directory). Every input a figure/data script
reads lives here, either reproducible (produced by a script in this repo)
or staged legacy data with no reproducing script — see REPRODUCIBILITY.md
for how the staged files were originally produced upstream of this repo.

A full, untrimmed copy of the pre-restructure `data/` (including everything
trimmed out below, e.g. the unused neighborhood-clustering parameter sweep)
is kept at `data.bk/` and is never deleted.

## clinical.parquet

Patient/core-level clinical annotations. **Reproducible** — written by
`scripts/data/export.py`.

## cells/

Per-cell tables, one row per cell, joinable on `(sample_id, object_id)`.

| file | provenance | consumers |
|---|---|---|
| `metadata.parquet` | reproducible — `export.py` | most figure scripts |
| `intensity.parquet` | reproducible — `export.py` | figure3 UMAP/heatmap scripts |
| `intensity_normalized.parquet` | reproducible — `export.py` | figure2/figure3 UMAP + heatmap scripts |
| `cell_annotation.parquet` | **staged legacy, no reproducing script** | KM/Cox survival scripts (figure6/7 + revision) |

## umap/

Ported UMAP embeddings. `UMAP.fit()` was never seeded in the legacy
pipeline, so none of these are reproducible by re-fitting — they are the
one-time ported ground truth for the published panels. See
REPRODUCIBILITY.md for the original `reducer.pkl` provenance.

| file | consumer | note |
|---|---|---|
| `all_cells.parquet` | `figure2_umap.py` | copy to `OUTPUT_FIGURES_DIR/figure2/reducer_embedding.parquet` before running |
| `caf_stromal.parquet`, `caf_markers_only.parquet`, `caf_excl_markers.parquet` | `figure3_caf_umap.py` | copy to `OUTPUT_FIGURES_DIR/figure3/<config>/umap_embedding.parquet` before running |
| `compartment_immune.parquet`, `compartment_epithelial.parquet`, `compartment_endothelial.parquet` | `figureS2_compartment_umap.py` | read directly from `DATA_DIR/umap/`, no copy step |

## niches/

All **staged legacy, no reproducing script** unless noted.

| file/dir | consumers |
|---|---|
| `clusters_annotated.parquet` | figure6/7 niche composition + violin scripts, figureS4b, figureS5a |
| `niche_annotations.csv` | figure5/6 niche heatmap/correlation/abundance (R) |
| `niche_annotations_source.xlsx` | `figure5_niche_annotation.py` (reproducing script for `clusters_annotated.parquet`, given the raw clustering output + this annotation source) |
| `niche_frequencies_per_tma_id.parquet` | `figure6_niche_abundance_heatmap.R` |
| `niche_frequencies_per_tma_id_stacked.parquet` | `figure5_niche_correlation.R` |
| `composition/*.parquet` | `figure5_niche_heatmap.R` |
| `patient_clustering/*.parquet` | figure4c, figureS3b, figureS3c, `figure4_metagroup_barplot.py` |
| `interactions/*.parquet` | `figure7def_circos_plots.py` |
| `robustness/` | legacy **code** (not data), `sys.path`-imported by `figure5_niche_clustering.py` and `figureS7_ari_robustness.py` for their k-means clustering/ARI utilities |

## neighborhoods/

Trimmed from the full neighborhood-graph sweep (~6.6G of unused
clustering-parameter configs live only in `data.bk/`) down to the 3 files
`figure5_niche_clustering.py` and `figureS7_ari_robustness.py` actually
read — the raw radius=32 neighborhood-composition graph, **staged legacy,
no reproducing script**.

| file | consumers |
|---|---|
| `cell_metadata.parquet` | figure5_niche_clustering.py, figureS7_ari_robustness.py |
| `metadata.parquet` | figure5_niche_clustering.py, figureS7_ari_robustness.py |
| `radius32_data.parquet` | figure5_niche_clustering.py, figureS7_ari_robustness.py |

## Not carried into `data/` (still in `data.bk/` only)

- `data/legacy/PCA_NHOODs_clean/` minus `niche_annotations_revised.xlsx`
  (→ `niches/niche_annotations_source.xlsx`) and `robustness/`
  (→ `niches/robustness/`) — the remaining ~65 files are old research
  code/notebooks not read by any script in this repo.
- `PCa_NHood/CellCellNeighborhoods/graph_type=radius-radius=32/{leiden-scanpy_*, scaled/, filtered/, filtered_2/}` —
  unused clustering-parameter sweep configs.
- Data required *exclusively* by `figure7_full_circos_plots.py` (raw
  per-sample anndata pickles) — that script is not supported going forward;
  see its module docstring.
- `legacy.bk/`, `figures.bk/`, `melissa-transfer/` — exact/stale duplicates
  and an unreferenced staging leftover, fully superseded by the layout above.
