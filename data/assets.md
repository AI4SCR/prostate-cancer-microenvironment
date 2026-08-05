# Data assets used by `scripts/figures/*` and `scripts/data/*`

Every external data file read by any current figure-plotting script, and
whether this repo can reproduce it. "Reproducible" = a script in this repo
regenerates it from more primitive inputs and it's been verified against
the precomputed copy; "Staged, unreproduced" = a read-only copy of a
precomputed/manually-curated legacy output with no reproducing script here;
"Live, unstaged" = data that exists on shared storage but hasn't been
copied into `data/legacy/` at all.

## `$EXPORT_DIR` assets (fully reproducible)

Produced by `scripts/data/export.py` from `$BASE_DIR` (the raw dataset).
Regenerate with `uv run python scripts/data/export.py`.

| Asset | Produced by | Used by |
|---|---|---|
| `metadata.parquet` | `scripts/data/export.py` | nearly every `scripts/figures/*` script (cell-level label/main_group/tma_id table) |
| `clinical.parquet` | `scripts/data/export.py` | nearly every `scripts/figures/*` script (patient/TMA-level clinical metadata) |
| `intensity.parquet` | `scripts/data/export.py` | `figure3_caf_umap_embedding.py` |
| `intensity_normalized.parquet` | `scripts/data/export.py` | `figure2_umap.py`, `scripts/data/figure2_umap_embedding.py` |
| `figures/figure2/reducer_embedding.parquet` | `scripts/data/figure2_umap_embedding.py` (full 2.19M-cell UMAP fit, slow) | `figure2_umap.py` |
| `figures/figure3/{config}/umap_embedding.parquet` | `scripts/data/figure3_caf_umap_embedding.py` | `figure3_caf_umap.py` |

## `resources/` assets (checked into the repo, ported from legacy)

| Asset | Produced by | Used by |
|---|---|---|
| `resources/colormaps.yaml` | ported verbatim from `000_paper/colormaps.yaml` (personal-machine path in legacy) | most `scripts/figures/*` scripts that color by `label`/`niche`/`main_group`/clinical variables |
| `resources/metalabels.yaml` | ported verbatim from `000_paper/metalabels.yaml` | `figure4_metagroup_barplot.py` |
| `resources/non_marker_channels.txt` | written by `scripts/data/export.py` from a hardcoded list (see `NON_MARKER_CHANNELS` in `prostate_cancer/utils.py`) | `figure2_umap.py`, `figure3_caf_umap.py` (via `scripts/data/*_embedding.py`) |

## `$LEGACY_DATA_DIR` assets

Corrected from the previous version of this file, which had gone stale --
`clusters_annotated_v2.parquet`/`niche_annotations_v2.csv` are now
confirmed reproducible (see below), and several newer scripts' inputs were
missing from the table entirely.

| Asset | Reproducible? | Produced by | Used by |
|---|---|---|---|
| `5-niches/annotation/clusters_annotated_v2.parquet` | **Yes** | `figure5_niche_clustering.py` + `figure5_niche_annotation.py`; verified byte-identical (all 2,051,915 cells' `niche`/`meta_niche`) -- but downstream scripts still read this staged copy directly rather than that pipeline's own output, since the dependency predates the fix (see `figure5_niche_annotation.py`'s docstring) | `figure5_niche_heatmap.R`, `figure6_inflammation_violin.R`, `figure6_km_niche6.R`, `figure6e_immune_risk_score_km.R`, `figure7a_stromogenic_violin.R`, `figure7b_niche9_km.R`, `figure7c_myCAF_km.R`, `figureS4b_niche_mean_composition.py`, `figure6c_niche_composition_filtered.py`, `figureS6bc_niche_km.R` |
| `5-niches/annotation/niche_annotations_v2.csv` | **Yes** | `figure5_niche_annotation.py` (same verification as above) | `figure5_niche_heatmap.R`, `figure5_niche_correlation.R`, `figure6_niche_abundance_heatmap.R` |
| `5-niches/annotation/clusters_annotated.parquet` (no `_v2`) | No | unknown (predates the `niche_annotations_revised.xlsx` fix) | not used by any current script. Its `niche` column disagrees with `_v2` on 73% of rows, but this is **not** a "stale, should be replaced" case -- see Discrepancies: `compute_interactions.py`/`visualize_interactions_lfc.py` (unported, not needed -- see the circos-pipeline section) genuinely computed the paper's results using this old assignment, only renaming output *files* to `_v2` niche names at the very end |
| `5-niches/annotation/niche_annotations.csv` (no `_v2`) | No | unknown, same vintage as `clusters_annotated.parquet` (no `_v2`) | not used by any current script |
| `5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet` | No | unknown (precomputed, no reproducing script) | `figure4_metagroup_barplot.py`, `figure4c_patient_cluster_km.R`, `figureS3b_cluster_concordance.R`, `figureS3c_progression_km.R` |
| `5-niches/barplot_data/metadata_with_dendrogram_colors_label_tma_id.parquet` | No | unknown (precomputed, no reproducing script) | `figureS3b_cluster_concordance.R` |
| `5-niches/frequencies/stacked_barplots/props_niche_pat_id.parquet` | No | unknown (precomputed, no reproducing script) | `figure6_niche_abundance_heatmap.R` -- **wrong-script stopgap**: this is patient-level data from `patient_heatmap.R`, the wrong legacy source for Figure 6a (paper legend is TMA/core-level); kept as-is only because the correct source's data is missing, see Missing data assets below |
| `5-niches/frequencies/stacked_barplots/props_niche_tma_id.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_correlation.R` |
| `5-niches/visualization/composition/niche_heatmap_data.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/mean_celltype_composition_per_niche.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/median_celltype_composition_per_niche.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/niche_abundance_stats.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `PCA_NHOODs_clean/niche_annotations_revised.xlsx` | No | manually curated by hand; supplied directly by the user 2026-07-31 | `figure5_niche_annotation.py` |
| `PCA_NHOODs_clean/robustness/utils/` (code, not data) | n/a | legacy codebase | `figure5_niche_clustering.py` (imports `utils.clustering`, `utils.visualization`) |
| `PCa_NHood/CellCellNeighborhoods/cell_metadata.parquet`, `metadata.parquet`, `graph_type=radius-radius=32/data.parquet` | No | unknown (undocumented upstream neighborhood-graph computation) | `figure5_niche_clustering.py` |

**Note**: `PCa_NHood/CellCellNeighborhoods/metadata.parquet` under `LEGACY_DATA_DIR` (464 rows, indexed by `tma_id`) is a *different, unrelated* file from the live, un-staged `metadata.parquet` at the same relative path under `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood/CellCellNeighborhoods/` (541 rows, indexed by short `sample_id` like `240121_003`) -- see the circos-pipeline section below. Same filename, different pipeline stage, not interchangeable.

## Circos-pipeline (Figure 7d-f) assets

`figure7def_circos_plots.py` (Fig 7d-f, niches 2/8/9) reads only the final,
already-computed per-niche interaction data -- it does NOT need the raw
anndata pickles, since the two upstream stages that would consume them
(`compute_interactions.py`, `visualize_interactions_lfc.py`) don't
themselves produce a paper panel and were skipped per direct user
instruction (their output already existed precomputed).

| Asset | Reproducible? | Produced by | Used by |
|---|---|---|---|
| `5-niches/visualization/interactions/redo/per_niche_lfc_above_median/dataframes_v2/{niche}.parquet` (niches 2/8/9 staged) | No | unknown (`compute_interactions.py` + `visualize_interactions_lfc.py`, not run in this repo -- see Discrepancies) | `figure7def_circos_plots.py` |

**Live, unstaged, not currently needed** (documented here in case
`compute_interactions.py`/`visualize_interactions_lfc.py` are ever ported):
found this session at
`/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood/CellCellNeighborhoods/`
(a different, broader live directory than the `LEGACY_DATA_DIR`-staged copy
of the same relative path, which only has a subset of its contents copied
in) -- `anndatas/{napari_sample_id}.pkl` (2,515 per-sample `AnnData`
pickles with `obsp['radius_32']`/`obsp['radius_48']` neighbor graphs) and
`metadata.parquet` (541 rows, indexed by the same short napari-style
`sample_id` as the pickle filenames -- not actually needed as a separate
input, since `clinical.parquet`'s own `napari_sample_id` column already
bridges to this repo's usual long-form `sample_id`). See Discrepancies for
the two issues that would need disclosed fixes to load/use these.

## Missing data assets (blocked panels)

Per `CLAUDE.md`'s hard constraint (added this session): when a script's
required input doesn't exist anywhere accessible, we report the gap and
leave the panel unfixed rather than writing a script to compute the
missing data ourselves. Currently one asset in this state:

| Asset | Needed by | Status |
|---|---|---|
| `niche_frequencies_per_tma_id.parquet` | `000_paper/sync_paper/06-spatial-niches/abundance/heatmap_frequencies.R` -- the correct legacy source for Figure 6a (TMA/core-level; row annotations `pat_id, os_status, disease_progr, gleason_grp, inflammation, stromogenic_smc_loss_reactive_stroma_present` match the published figure exactly, confirmed by direct user identification) | **Not found anywhere.** Checked `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood` (nothing matching `niche_frequencies*` at all) and `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/5-niches/frequencies/` (has a similarly-named `niche_frequencies_per_tma.parquet` -- no `_id` -- but its columns use an older, pre-"revised annotation" niche naming scheme, e.g. `TLS_Bcells_Tcells`, `canonical_BLepithelium`, completely different from the `_v2` niche names used throughout this repo; not the same data, not a usable substitute). Figure 6a stays on the wrong (patient-level `patient_heatmap.R`) script in the meantime -- see `open-questions.md`. |

## Recovered legacy UMAP embeddings (`data/figures/`)

Not legacy in the "external/manually curated" sense above -- extracted from
the original pipeline's `reducer.pkl` files (pixi env + a pynndescent
compat patch, see session history) and kept out of `data/legacy/` since
they're figure-script-scoped caches, not consolidated external data. Still
unused (confirmed by grep) -- re-checked this session, no change.

| Asset | Produced by | Used by |
|---|---|---|
| `figures/figure2_umap/umap_embeddings.parquet` | pkl-to-parquet extraction (not a repo script) | **unused** -- `figure2_umap_embedding.py` always recomputes UMAP from scratch instead |
| `figures/figure3_caf_umap/{excl_markers,caf_markers_only}/umap_embeddings.parquet` | pkl-to-parquet extraction (not a repo script) | **unused** -- `figure3_caf_umap_embedding.py` always recomputes UMAP from scratch instead |

## Discrepancies found while verifying this file

1. **`clusters_annotated_v2.parquet`/`niche_annotations_v2.csv` are actually
   reproducible** -- the previous version of this file marked them
   "unknown, no reproducing script," which was stale; `figure5_niche_clustering.py`
   + `figure5_niche_annotation.py` produce byte-identical output, confirmed
   in that script's own docstring. Downstream scripts still read the staged
   copy directly (a deliberate, disclosed choice, not an oversight).
2. **`clusters_annotated.parquet` (no `_v2`) initially looked stale, but
   isn't a bug to fix** -- its `niche` column disagrees with the
   verified-correct `_v2` file on 73% of rows. First instinct was that
   `compute_interactions.py` naming this exact file was an oversight and
   should be pointed at `_v2` instead. Wrong: `visualize_interactions_lfc.py`
   (read in full to check) confirms the *entire computation* -- which cells
   get grouped into which niche, all interaction-frequency/LFC math --
   genuinely ran on this old assignment to produce the paper's actual
   results; only the output *filenames* get remapped from old to new niche
   names at the very end (via a lookup on the raw k-means cluster ID,
   present in both files), before `visualize_circos_plot.py` reads them.
   Substituting `_v2` throughout the computation would have been the real
   deviation. Moot for this repo either way, since `compute_interactions.py`/
   `visualize_interactions_lfc.py` were skipped (their already-computed
   `_v2`-renamed output exists precomputed, see below).
3. **Sample-ID namespace mismatch, `anndata` pickle version incompatibility**
   -- both documented in the circos-pipeline section above. Neither
   affects any implemented script, including `figure7def_circos_plots.py`
   (reads only the final precomputed data, never touches the raw pickles or
   `clusters_annotated.parquet`) -- kept here only in case
   `compute_interactions.py`/`visualize_interactions_lfc.py` are ever
   ported for some other reason.
