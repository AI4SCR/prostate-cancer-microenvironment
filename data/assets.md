# Data assets used by `scripts/figures/*` and `scripts/data/*`

Every external data file read by any current figure-plotting script, and
where it's currently loaded from. "Reproducible" = a script in this repo
regenerates it from more primitive inputs and it's been verified against
the precomputed copy; "Staged, unreproduced" = a read-only copy of a
precomputed/manually-curated legacy output with no reproducing script here;
"Not found" = checked and confirmed absent everywhere searched.

**Search scope for this pass**: every legacy script's own hardcoded path
(even where it points to an inaccessible personal folder, e.g.
`/users/mensmeng/workspace/...`), cross-referenced against three live
directories plus their `data/legacy/`-staged copies:
- `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa_NHood/` (older, `mensmeng`-owned)
- `/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment/PCa_NHood/`
  (newer, project-specific staging -- mirrors the exact subdirectory
  structure of legacy's hardcoded `mensmeng` paths, e.g.
  `PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/`)
- `/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment/PCA_NHOODs_clean/`
  (newer, project-specific staging of the `PCA_NHOODs_clean` codebase/data)
- (also checked, since several assets live there too:
  `/work/FAC/FBM/DBC/mrapsoma/prometex/data/PCa/` and its newer sibling
  `/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment/`)

Where multiple candidate copies exist, compared by content hash (`md5sum`)
and mtime, not just filename -- results noted per asset below.

## `$EXPORT_DIR` assets (fully reproducible)

Produced by `scripts/data/export.py` from `$BASE_DIR` (the raw dataset).
Regenerate with `uv run python scripts/data/export.py`. Unaffected by this
pass -- these are already fully reproducible in-repo, not staged legacy data.

| Asset | Produced by | Used by |
|---|---|---|
| `metadata.parquet` | `scripts/data/export.py` | nearly every `scripts/figures/*` script (cell-level `label`/`main_group`/`tma_id` table -- **no `niche` column**, see the `cell_annotations` entry below) |
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

**Found, but not the right file**: `/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment/colormaps.yaml`
exists at the top level of the newer project-staged tree -- confirmed via
`diff` this is a *different, incomplete* colormaps file (missing the
`niche`/`metagroups` sections entirely, different hex values for shared
keys like `disease_progr`/`gleason_grp`). Not a newer version of
`resources/colormaps.yaml`; don't substitute it in.

## `$LEGACY_DATA_DIR` assets (`5-niches/`)

For every asset below, content was verified **byte-identical** (`md5sum`)
across all copies found: the currently-staged `data/legacy/5-niches/...`,
the older `/work/.../data/PCa/5-niches/...`, and the newer
`/work/.../data/prostate-cancer-microenvironment/5-niches/...`. So "which is
newest" doesn't matter for these -- they're the same data everywhere,
already correctly staged.

| Asset | Reproducible? | Produced by | Used by |
|---|---|---|---|
| `5-niches/annotation/clusters_annotated_v2.parquet` | **Yes** | `figure5_niche_clustering.py` + `figure5_niche_annotation.py`; verified byte-identical (all 2,051,915 cells' `niche`/`meta_niche`) -- but downstream scripts still read this staged copy directly rather than that pipeline's own output, since the dependency predates the fix (see `figure5_niche_annotation.py`'s docstring) | `figure5_niche_heatmap.R`, `figure6_inflammation_violin.R`, `figure6c_niche_composition_filtered.py`, `figureS4b_niche_mean_composition.py`, `figure7a_stromogenic_violin.R` |
| `5-niches/annotation/niche_annotations_v2.csv` | **Yes** | `figure5_niche_annotation.py` (same verification as above) | `figure5_niche_heatmap.R`, `figure5_niche_correlation.R`, `figure6_niche_abundance_heatmap.R` (blocked, see below) |
| `5-niches/annotation/clusters_annotated.parquet` (no `_v2`) | No | unknown (predates the `niche_annotations_revised.xlsx` fix) | not used by any current script. Its `niche` column disagrees with `_v2` on 73% of rows -- **not a "stale, should be replaced" case**, see Discrepancies below |
| `5-niches/annotation/niche_annotations.csv` (no `_v2`) | No | unknown, same vintage as `clusters_annotated.parquet` (no `_v2`) | not used by any current script |
| `5-niches/barplot_data/metadata_with_dendrogram_colors_label_pat_id.parquet` | No | unknown (precomputed, no reproducing script) | `figure4_metagroup_barplot.py`, `figure4c_patient_cluster_km.R`, `figureS3b_cluster_concordance.R`, `figureS3c_progression_km.R` (blocked, see below) |
| `5-niches/barplot_data/metadata_with_dendrogram_colors_label_tma_id.parquet` | No | unknown (precomputed, no reproducing script) | `figureS3b_cluster_concordance.R` |
| `5-niches/frequencies/stacked_barplots/props_niche_pat_id.parquet` | No | unknown (precomputed, no reproducing script) | not used by any current script (was `figure6_niche_abundance_heatmap.R`'s input before that script was corrected to the TMA-level `heatmap_frequencies.R` source) |
| `5-niches/frequencies/stacked_barplots/props_niche_tma_id.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_correlation.R` |
| `5-niches/visualization/composition/niche_heatmap_data.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/mean_celltype_composition_per_niche.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/median_celltype_composition_per_niche.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |
| `5-niches/visualization/composition/niche_abundance_stats.parquet` | No | unknown (precomputed, no reproducing script) | `figure5_niche_heatmap.R` |

**Also found in the newer `5-niches/` tree, not currently used by any
script** (checked, not needed): `visualization/composition/meta_niche_heatmap_data.parquet`,
`meta_niche_heatmap.png` -- a meta-niche-level (not niche-level) variant of
the z-score heatmap data. No current script reads it.

## `$LEGACY_DATA_DIR` assets (other)

| Asset | Reproducible? | Produced by | Used by |
|---|---|---|---|
| `PCA_NHOODs_clean/niche_annotations_revised.xlsx` | No | manually curated by hand; supplied directly by the user 2026-07-31 | `figure5_niche_annotation.py` |
| `PCA_NHOODs_clean/robustness/utils/` (code, not data) | n/a | legacy codebase | `figure5_niche_clustering.py` (imports `utils.clustering`, `utils.visualization`) |
| `PCa_NHood/CellCellNeighborhoods/cell_metadata.parquet`, `metadata.parquet`, `graph_type=radius-radius=32/data.parquet` | No | unknown (undocumented upstream neighborhood-graph computation) | `figure5_niche_clustering.py` -- **already validated** (output matches the published figure, 2,051,915 cells, confirmed in earlier work); not re-pointed at the newer `PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/graph_type=radius-radius=32/` copy found this pass, since that copy has a different internal substructure (different `leiden-*` param subfolders, no `filtered`/`filtered_2`) and switching an already-validated, working pipeline to an unverified alternate source isn't warranted |

**Note**: `PCa_NHood/CellCellNeighborhoods/metadata.parquet` under
`LEGACY_DATA_DIR` (464 rows, indexed by `tma_id` -- TMA-level clinical
metadata) is unrelated to the *per-cell* `cell_annotations`/`metadata.parquet`
question below (a different asset, different pipeline stage, coincidentally
similar filename).

## Circos-pipeline (Figure 7d-f) assets

`figure7def_circos_plots.py` (Fig 7d-f, niches 2/8/9) reads only the final,
already-computed per-niche interaction data -- it does not need the raw
anndata pickles, since the two upstream stages that would consume them
(`compute_interactions.py`, `visualize_interactions_lfc.py`) don't
themselves produce a paper panel and were skipped per direct user
instruction (their output already existed precomputed).

| Asset | Reproducible? | Produced by | Used by |
|---|---|---|---|
| `5-niches/visualization/interactions/redo/per_niche_lfc_above_median/dataframes_v2/{niche}.parquet` (niches 2/8/9 staged) | No | unknown (`compute_interactions.py` + `visualize_interactions_lfc.py`, not run in this repo -- see Discrepancies) | `figure7def_circos_plots.py` |

**Correction from an earlier pass -- the anndata sample-ID "mismatch" was a
wrong-location artifact, not a real discrepancy.** Previously documented:
per-sample anndata pickles found at
`/work/.../data/PCa_NHood/CellCellNeighborhoods/anndatas/` (2,515 files),
named with a short napari-style id (e.g. `240217_005`) that doesn't match
`clusters_annotated_v2.parquet`'s long-form `sample_id` -- requiring a
bridge through `clinical.parquet`'s `napari_sample_id` column, and only
matching 515/541 (95%) of files.

This pass found the actual correct location, matching legacy's own
hardcoded path exactly (`compute_interactions.py`'s
`/users/mensmeng/workspace/nhoods/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas`):
`/work/FAC/FBM/DBC/mrapsoma/prometex/data/prostate-cancer-microenvironment/PCa_NHood/final_analysis/evaluation/proportion/CellCellNeighborhoods/anndatas/`
-- 476 files, named with the long-form `sample_id` directly (e.g.
`231204_ibl_x2y4_29_11.pkl`). Verified: **100% overlap** with
`clusters_annotated_v2.parquet`'s 476 unique `sample_id` values, zero
missing either direction. No ID-bridging needed at all at this location.
Still requires the same `anndata`-version compatibility monkeypatch to load
(see Discrepancies) -- that part of the earlier finding stands.

The older, napari-ID-named location is a different (larger, 2,515-file)
dataset from an earlier/different processing run; not what legacy's own
path structure points to, and not needed now that the correct location is
identified.

## Missing data assets (blocked panels)

Per `CLAUDE.md`'s hard constraint: when a script's required input doesn't
exist anywhere accessible, report the gap and leave the panel unfixed
rather than computing the missing data ourselves. Two assets confirmed
genuinely absent, searched across all locations listed at the top of this
file (not just the two originally checked):

| Asset | Needed by | Status |
|---|---|---|
| `5-niches/frequencies/niche_frequencies_per_tma_id.parquet` | `figure6_niche_abundance_heatmap.R`, ported from `000_paper/sync_paper/06-spatial-niches/abundance/heatmap_frequencies.R` (the correct Figure 6a source, row annotations confirmed matching the published figure exactly) | **Not found anywhere**, including the newer `5-niches/frequencies/` tree checked this pass (which has `props_niche_tma.parquet`/`props_niche_pat_id.parquet`/etc., but no `_tma_id` variant, and all of them use an older, pre-"revised annotation" niche naming scheme, e.g. `TLS_Bcells_Tcells` -- confirmed not usable, same issue as the old `PCa/5-niches/frequencies/` copy). `figure6_niche_abundance_heatmap.R` is a strict verbatim port and fails at exactly this `read_parquet()` call. |
| `cell_annotations` / `cell_annotation.parquet` (sync_paper's raw per-cell input, distinct from `clusters_annotated_v2.parquet` -- described in `sync_paper/config/paths.yaml` as "Cell annotations" vs. clusters_annotated's "Annotated cell clusters") | `figure6e_immune_risk_score_km.R`, `figure6_km_niche6.R`, `figure7b_niche9_km.R`, `figureS6bc_niche_km.R`, `figure7c_myCAF_km.R` (all read `paths$cell_annotations` in their `sync_paper` sources) | **Not found anywhere**, including both `PCa_NHood` trees and `PCA_NHOODs_clean`. This repo's `EXPORT_DIR/metadata.parquet` is used as the closest analog (same conceptual role: raw per-cell `label`/`main_group` table) -- confirmed it lacks the `niche` column these scripts need, so all five scripts fail inside `compute_label_frequency()`, matching what the actual `cell_annotation.parquet` (if it lacks `niche` too, as its "raw annotations" vs. "annotated clusters" naming suggests) would also do. |

## Recovered legacy UMAP embeddings (`data/figures/`)

Not legacy in the "external/manually curated" sense above -- extracted from
the original pipeline's `reducer.pkl` files (pixi env + a pynndescent
compat patch, see session history) and kept out of `data/legacy/` since
they're figure-script-scoped caches, not consolidated external data. Still
unused (confirmed by grep).

| Asset | Produced by | Used by |
|---|---|---|
| `figures/figure2_umap/umap_embeddings.parquet` | pkl-to-parquet extraction (not a repo script) | **unused** -- `figure2_umap_embedding.py` always recomputes UMAP from scratch instead |
| `figures/figure3_caf_umap/{excl_markers,caf_markers_only}/umap_embeddings.parquet` | pkl-to-parquet extraction (not a repo script) | **unused** -- `figure3_caf_umap_embedding.py` always recomputes UMAP from scratch instead |

## Discrepancies found while verifying this file

1. **Anndata sample-ID mismatch was a wrong-location artifact, not a real
   discrepancy** -- see the Circos-pipeline section above for the full
   correction. The actual path matching legacy's own hardcoded location has
   100% ID overlap with `clusters_annotated_v2.parquet`, no bridging needed.
2. **`anndata` pickle version incompatibility still applies** at the
   correct location too -- pickled with an older `anndata` whose internal
   `AnnDataFileManager` used the state key `_adata`; the installed
   `anndata` 0.12.10 expects `_adata_ref` and raises `KeyError` otherwise.
   Confirmed these objects aren't actually file-backed (`_filename`/`_file`
   are `None`), so it's a pure key-rename compatibility issue. Verified the
   correct-location pickles load fine with the same monkeypatch used
   before.
3. **`clusters_annotated.parquet` (no `_v2`) is not stale, just used
   differently** -- its `niche` column disagrees with the verified-correct
   `_v2` file on 73% of rows. `visualize_interactions_lfc.py` (read in full)
   confirms the entire interaction computation genuinely ran on this old
   assignment to produce the paper's actual results; only the output
   *filenames* get remapped to `_v2` niche names at the very end, before
   `visualize_circos_plot.py` reads them. Moot for this repo either way,
   since `compute_interactions.py`/`visualize_interactions_lfc.py` were
   skipped (their output already exists precomputed as `dataframes_v2/`).
4. **A look-alike `colormaps.yaml` exists in the live data tree but isn't
   the right file** -- see the `resources/` section above.
